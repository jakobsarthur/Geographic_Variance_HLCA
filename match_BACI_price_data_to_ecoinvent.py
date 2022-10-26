import numpy as np
import pandas as pd
import os
import ray
import json
import warnings
import pickle
import time
import datatable as dt
from functools import reduce
from operator import or_
import pdb


def Match_BACI_data_to_ecoinvent(
        PRO,
        path_to_BACI_data,
        outputDir=None,
        USD_EURO_exr=None,
        Nsamples=3000,
        mapping_data_dir='../mapping_data/',
        region_mapping_dict_name='ecoinvent35-baci_region_mapping.json',
        commodity_mapping_dict_name='ecoinvent35-HS12_mapping.json',
        exio_base_country_mapping_file='eco_baci_country_mapping.xlsx',
        n_cores=None
        ):
    """
    This function matches BACI price distribution data to ecoinvent activties
    based on their reference product and geography. The distribution is the
    volume (in tons) weighted price distribution of all trade flows exported from
    the activity's region (geography). For every actitvity with a match it outputs
    a pickle file in the outputDir. It relies on Ray multiprocessing.
    
    Input:

    PRO                     Pandas DataFrame with ecoinvent process meta data.
                            output of ecospold2matrix:
                            https://github.com/majeau-bettez/ecospold2matrix
    path_to_BACI_data       Path to BACI trade data csv file
    outputDir               Directory to save the output pickles
    USD_EURO_exr            Exchange rate USD to EURO to use. the USD price will
                            divided by this number to obtain the euro price.
                            Default is 1.
    Nsamples                Number of samples to draw from the price distribution
    mapping_data_dir        Directory with the mapping dictionaries. Default:
                            '../mappind_data'
    region_mapping_dict_name
                            Name of region mapping dict in the above directory.
                            Default: 'ecoinvent35-baci_region_mapping.json',
    commodity_mapping_dict_name
                            Name of the commodity mapping dictionary in the
                            mapping directory. Default: 'ecoinvent35-HS12_mapping.json',
    exio_base_country_mapping_file
                            Name of the xlsx file with a mapping between EXIOBASE
                            regions, BACI codes and country names
                            Default: 'eco_baci_country_mapping.xlsx'
    n_cores                 Number of cores to use for the multiprocessing.
                            Default: N_max-1, where N_max = # available cpu's
   
    """
    #Save time for overall walltime
    t0 = time.time()

    # Check if necessary inputs have been given
    if outputDir==None:
        raise Exception("Please provide a path to save the pickle files")
    if USD_EURO_exr==None:
        warnings.warn("No exchange given. Default rate of 1 will be used")
        USD_EURO_exr = 1 
    
    # Check if mapping dictionaries are available
    mapping_data = os.path.realpath(mapping_data_dir)
    if region_mapping_dict_name==None or commodity_mapping_dict_name==None:
        if not ((region_mapping_dict_name in os.listdir(mapping_data)) and
                (commodity_mapping_dict_name in os.listdir(mapping_data))):
            raise Exception("Please provide directory with mapping dictionaries\
                    between activities and regions and HS12 commodities. As well\
                    as the correct file names if different from dedault")
    # Check if BACI file exists:
    if not os.path.isfile(path_to_BACI_data):
        raise Exception("Could not find '{}'. Please provide a valid path to\
                          BACI data file".format(path_to_BACI_data))
    
    # Load dictionaries
    eco_HS12_mapping, eco_baci_region_mapping_dic = read_mapping_dicts(
            mapping_data, region_mapping_dict_name, commodity_mapping_dict_name)

    #Load in country_mapping to exio regions:
    exio_country_mapping = read_exio_country_mapping(mapping_data_dir,
                                                 exio_base_country_mapping_file)
 
    # Load BACI data:
    baci_data = readDataBACI_dt(path_to_BACI_data, USD_EURO_exr)

    # Check if outDir exists and if not make it.
    outputDir = os.path.realpath(outputDir)
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
        print('Created the output directory {}'.format(outputDir))

    # Prepare dictionary for fast row lookup:
    print('prepare lookup dictionary for faster row lookup...')
    baci_i_row_dic = {}
    for reg in np.unique(baci_data[:,'i']):
        baci_i_row_dic[reg] = np.where(baci_data[:,'i'].to_numpy().flatten() == reg)[0].tolist()

    # If not specified, use 1 less than the available number of cores.
    if n_cores==None:
        n_cores = os.cpu_count()-1
    
    print("Initializing Ray multiprocessing with {} cores".format(n_cores))
    ray.shutdown()
    ray.init(num_cpus=n_cores)
   
    eco_baci_region_mapping_dic_id = ray.put(eco_baci_region_mapping_dic)
    eco_HS12_mapping_id = ray.put(eco_HS12_mapping)
    baci_data_id = ray.put(baci_data)
    baci_i_row_dic_id = ray.put(baci_i_row_dic)
    exio_country_mapping_id = ray.put(exio_country_mapping)
    


    results_ids = [get_baci_price_data.remote(act.Index, act,
                                   eco_baci_region_mapping_dic_id,
                                   eco_HS12_mapping_id, baci_data_id,
                                   baci_i_row_dic_id, exio_country_mapping_id,
                                   Nsamples, outputDir)
            for act in PRO.loc[(PRO['unitName']=='kg') &
            (PRO['cpc'].str.split(':').str.get(0).str.len() >= 4)].itertuples()]
    _ = ray.get(results_ids)
    #for i,act in enumerate(PRO.loc[(PRO['unitName']=='kg') &
    #    (PRO['cpc'].str.split(':').str.get(0).str.len() >= 4)].itertuples()):
    #    # only iterate over processes with the right units and with a cpc 
    #    # code that has at least 4 digits

    #    if i%100 ==0:
    #        print('{} of total {}'.format(i+1, len(PRO.loc[(PRO['unitName']=='kg') &
    #                     (PRO['cpc'].str.split(':').str.get(0).str.len() >=4)])))
    #    
    #    get_baci_price_data.remote(act.Index, act,
    #                               eco_baci_region_mapping_dic,
    #                               eco_HS12_mapping, baci_data,
    #                               baci_i_row_dic, exio_country_mapping,
    #                               Nsamples, outputDir)
    
    


#    for i,act in enumerate(PRO[100:110].append(PRO.loc[t]).itertuples()):
#        # only iterate over processes with the right units and with a cpc 
#        # code that has at least 4 digits
#
#        get_baci_price_data(act.Index, act,
#                                   eco_baci_region_mapping_dic,
#                                   eco_HS12_mapping, baci_data,
#                                   baci_i_row_dic, exio_country_mapping,
#                                   Nsamples, outputDir)
    print('Done matching activties to BACI price data')
    # Close ray remote
    print('Shutting down Ray multiprocessing')
#    ray.shutdown()
    dt = time.time()-t0
    print("Elapsed wall time: {} hours, {}, minutes, and {:.0f} seconds.".format(
            dt//3600, (dt%3600)//60, dt%60))
    

def readDataBACI(path_to_BACI_data, USD_EURO_exr=1):
    """
    Input:
    path_to_BACI_data   path to BACI data file for year to use

    USD_EURO_exr        USD per EURO exchange rate default is 1. See EUROSTAT:
                        https://ec.europa.eu/eurostat/databrowser/bookmark/bc107428-d077-4a9a-86a5-a4f045cf63e9?lang=en

    """
    print('reading in BACI data...')
    baci_data = pd.read_csv(path_to_BACI_data, sep=',',
                            dtype={'t':int, 'i':int,
                                   'j':int, 'k':str,
                                   'v':float, 'q':float})
    baci_data['p'] = baci_data['v']/baci_data['q']
    baci_data['p_euro'] = baci_data['p']/USD_EURO_exr
    return baci_data

def readDataBACI_dt(path_to_BACI_data, USD_EURO_exr=1):
    """
    Input:
    path_to_BACI_data   path to BACI data file for year to use

    USD_EURO_exr        USD per EURO exchange rate default is 1. See EUROSTAT:
                        https://ec.europa.eu/eurostat/databrowser/bookmark/bc107428-d077-4a9a-86a5-a4f045cf63e9?lang=en

    """
    print('reading in BACI data...')
    baci_data = dt.fread(path_to_BACI_data, sep=',')
    baci_data[:,'k']=baci_data[:,dt.as_type(dt.f.k, str)]
    baci_data['p'] = baci_data[:, dt.f.v/dt.f.q]
    baci_data['p_euro'] = baci_data[:,dt.f.p/USD_EURO_exr]
    return baci_data


def read_mapping_dicts(mapping_dir, region_mapping_dict_name,
                commodity_mapping_dict_name):

    with open(os.path.join(mapping_dir, commodity_mapping_dict_name), 'r') as fh:
        eco_HS12_mapping = json.load(fh)
    with open(os.path.join(mapping_dir, region_mapping_dict_name), 'r') as fh:
        eco_baci_region_mapping_dic = json.load(fh)

    return eco_HS12_mapping, eco_baci_region_mapping_dic


def read_exio_country_mapping(mapping_dir, mapping_file):
    country_mapping = pd.read_excel(
                            os.path.join(mapping_dir, mapping_file),
                            sheet_name='Ecoinvent-BACI96', header=0,
                            index_col=0, keep_default_na=False, na_values=['']
                            )
    return country_mapping





@ray.remote
def get_baci_price_data(proc_index, act, eco_baci_region_mapping_dic,
                        eco_HS12_mapping, baci_data, baci_i_row_dic,
                        exio_country_mapping, draw_nsamples=3000, outDir=None):
    """Maps baci prices to ecoinvent processes. Uses volume weighted price
    price distribution from which it draws 'draw_nsamples=3000' samples.
    Not all ecoinvent activity reference products have 3000 flows which means
    multiple prices will be sampled more often. This however will not affect the
    distribution in a significant way.
    Input:
    proc_index          the UUID of the ecoinvent activity
    act                 all meta data from the PRO dataframe for this activity
    eco_baci_region_mapping_dic
                        dict mapping ecoinvent regions to BACI region codes
    eco_HS12_mapping    dict mapping ecoinvent act to HS12 commodities
    baci_data           Dataframe with BACI data
    baci_i_row_dic      Dictionary mapping rows for each exporting country
    exio_country_mapping
                        Dataframe with mapping between EXIOBASE regions and BACI
                        country codes
    draw_nsamples       Numberof samples to draw from the price distribution
    outDir              Directory to save the pickle files
    """
    try:
        eco_HS12_mapping[proc_index]
    except KeyError:
        pdb.set_trace()

    # Check if process has a BACI product mapping
    if eco_HS12_mapping[proc_index] != None:
        hs12_codes = eco_HS12_mapping[proc_index]
        # get matching baci regions from dictionary
        if act.geography != 'RoW':
            baci_regs = eco_baci_region_mapping_dic[act.geography]
        else:
            baci_regs = eco_baci_region_mapping_dic[proc_index]

        # get matching exiobase regions from exiobase-baci mapping
        exioregs_mapping = exio_country_mapping.set_index('CODE_BACI').loc[
                          baci_regs, 'DESIRE code']
        exioregs_mapping = exioregs_mapping.reset_index()
        #some baci regions are duplicate with NAN DESIRE code so drop these
        exioregs_mapping.dropna(inplace=True)
        exioregs_mapping.set_index('DESIRE code', inplace=True)

        # Make empty dictonary for process meta data and price data
        price_dic = {}
        # store process meta data
        price_dic['activityName'] = act.activityName
        price_dic['geography'] = act.geography
        price_dic['productName'] = act.productName
        price_dic['cpc'] = act.cpc
        price_dic['unitName'] = act.unitName
        price_dic['pricesByCountry'] = {}


        all_unique_regs = []
        # Now loop over unique exriobase regions to get a associated price distribution
        for exioreg in exioregs_mapping.index.unique():
            try:
                regs_baci = exioregs_mapping.loc[exioreg, 'CODE_BACI'].unique()
            # If only a single value is returned the above returns a attribute error
            except AttributeError:
                regs_baci = [exioregs_mapping.loc[exioreg, 'CODE_BACI']]
            if len(hs12_codes) <= 4:
                indices_regs = [row for reg in regs_baci if reg in baci_i_row_dic.keys() for row in baci_i_row_dic[reg]]
                indices_prods = [row for prod in hs12_codes for row in np.where(baci_data[indices_regs,'k'].to_numpy().flatten() ==prod)[0]]
                total_mask = list(set(indices_regs).intersection(indices_prods))
                total_mask.sort()

                prices_euro = baci_data[indices_regs, 'p_euro'][indices_prods,:].to_numpy().flatten()
                weights = baci_data[indices_regs,'q'][indices_prods,:].to_numpy().flatten() 
            else:
                prices_euro = baci_data[isin('i', regs_baci) &
                        isin('k', hs12_codes),'p_euro'].to_numpy().flatten()
                weights = baci_data[isin('i', regs_baci) &
                        isin('k', hs12_codes),'q'].to_numpy().flatten()
            
            try:
                avg, std = weighted_avg_and_std(prices_euro, weights)
                sample_price = np.random.choice(prices_euro, size=draw_nsamples,
                                                p=weights/weights.sum())
                price_dic['pricesByCountry'][exioreg] = {}

                price_dic['pricesByCountry'][exioreg]['prices_euro'] = prices_euro
                price_dic['pricesByCountry'][exioreg]['weights'] = weights
                price_dic['pricesByCountry'][exioreg]['price_sample'] = sample_price
                price_dic['pricesByCountry'][exioreg]['price_baci_mean'] = avg
                price_dic['pricesByCountry'][exioreg]['price_baci_std'] = std
                price_dic['pricesByCountry'][exioreg]['price_percentiles'] = np.percentile(sample_price, [2.5,16,50,84,97.5])
                price_dic['pricesByCountry'][exioreg]['price_baci_min'] = np.min(prices_euro)
                price_dic['pricesByCountry'][exioreg]['price_baci_max'] = np.max(prices_euro)
                price_dic['pricesByCountry'][exioreg]['nr_baci_flows'] = len(prices_euro)

            except ZeroDivisionError:
                # print('Zero Divission Error, zero volume flows')
                pass

            except ValueError:
                print(proc_index)
                print(price_dic)
                pdb.set_trace()

        # Check id the prices by country dictionary is not empty and save
        if price_dic['pricesByCountry']:
            with open(os.path.join(outDir, '{}.pickle'.format(proc_index)), 'wb') as fh:
                pickle.dump(price_dic, fh)
        # Otherwise delete dictionary and continue
        else:
            del price_dic
    return 0


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))

def isin(column, iterable):
    content = [dt.f[column] == entry
               for entry in iterable]
    return reduce(or_, content)
