# import pandas as pd
import numpy as np
import sys
sys.path.append('/home/jakobs/Documents/IndEcol/OASES/pylcaio/src/')
import pylcaio

import hybrid_lcaio

import os
import gzip
import pickle

import scipy.sparse
from pypardiso import spsolve, factorized
import time
import ray
import h5py
from filelock import FileLock

##number of runs
#N_runs= 10000
#constant_price = False
#use_ecoinvent_price = False # this overrides the constant price flag
#constant_region_var_price = False # this overrides the above two flags
#
#N_cores = 2
#
#lcaio_object_path = '/home/jakobs/Documents/IndEcol/OASES/EuropeanConsumptionFootprint/Databases/plcaio_object_STAM_saved_both_strategies_True_ecoinvent_3_8_EXIOBASE_3_8_2.pickle'
#directory_baci_price_data = '/home/jakobs/Documents/IndEcol/OASES/EuropeanConsumptionFootprint/Databases/price_data_per_region/2015/'
## C_io_path = '/home/jakobs/Documents/IndEcol/Data/EXIOBASE/exiobase3_6/Characterization_EB36.xlsx'
#price_data_feather_file = '/home/jakobs/Documents/IndEcol/OASES/EuropeanConsumptionFootprint/Databases/price_dataframe_ecoinvent_3_8_2015_regionally_aggregated.ftr'
#results_file_name = '/home/jakobs/Documents/IndEcol/OASES/EuropeanConsumptionFootprint/Databases/Results_geography_MC_BACI_price_{}_runs_3FPs_ProdVolWeights_constantprice_mean_{}_ecoinventprice_{}_constant_region_{}_only_Cu_impacts.hdf5'.format(N_runs, constant_price, use_ecoinvent_price, constant_region_var_price)
#hdb_object_file_name = '/home/jakobs/Documents/IndEcol/OASES/EuropeanConsumptionFootprint/Databases/HDB_STAM_ecoinvent_3_8_exiobase_3_8_2.pickle'
#double_counting_method = 'binary'
#exiobase_version = '3_8_2'
#ecoinvent_version = '3_8'


def initiate_hdb(
        ecoinvent_version=None,
        exiobase_version=None,
        hdb_object_file_name=None,
        lcaio_object_path=None,
        Database_dir=None,
        ):

    #Check if the necessary arguments have been given:
    arguments_not_passed = []
    if ecoinvent_version is None: arguments_not_passed.append('ecoinvent_version')
    if exiobase_version  is None: arguments_not_passed.append('exiobase_version')
    if Database_dir      is None: arguments_not_passed.append('Database_dir')
    if arguments_not_passed:
        raise NameError('please provide the following arguments: \n {}'.format(
                        '\n'.join(arguments_not_passed)))

    if not hdb_object_file_name is None and os.path.exists(hdb_object_file_name):
        print('loading Hybrid Database Object...')
        with open(hdb_object_file_name, 'rb') as fh:
            hdb = pickle.load(fh)
        calc_lca_impacts_flag = False

    else:
        #check if lcaio_object_path exists:
        if lcaio_object_path is None or not os.path.exists(lcaio_object_path):
            raise FileNotFoundError('Please provide a valid lcaio_object_path'+
                                    ' or a valid hdb_object_file_name')
        # check if hdb_object_file_name was passed
        if hdb_object_file_name is None:
            hdb_object_file_name = os.path.join(Database_dir, 
                            'HDB_ecoinvent_{}_exiobase_{}.pickle'.format(
                            ecoinvent_version, exiobase_version)
                            )
            print('No hdb_object_file_name was passed. hdb_object_file_name set to {}'.format(hdb_object_file_name))

        print("Reading in lcaio object, this takes a while...")
        with gzip.open(lcaio_object_path, 'rb') as fh:
            lcaio_object = pickle.load(fh)

        print("Creating hybrid database object...")
        hdb = hybrid_lcaio.Hybrid_LCAIO(lcaio_object)
        hdb.impact_names_lca = lcaio_object.extended_impact_names_IW_eco

        # because first time load, need to calculate lca impacts too
        hdb.calc_lca_impacts_flag = True

#        print('Reading in IO Characterisation factors...')
#        C_io = pd.read_excel(C_io_path, index_col=0, header=0)
#        C_io = C_io.T
#        Impact_names_io = C_io.index.str.rstrip().values
#        C_io = scipy.sparse.csr_matrix(C_io)
        print('Calculate IO multiplyers...')
        hdb.calc_M_io(lcaio_object.C_io, lcaio_object.extended_impact_names_IW_exio)

        # No delete lcaio_object to clear up some memory
        del lcaio_object

        print('Calculate L_lca, this takes a while...')
        hdb.calc_L_lca()

        # save ecoinvent and exiobase versions
        hdb.ecoinvent_version = ecoinvent_version
        hdb.exiobase_version = exiobase_version
        hdb.Database_dir = Database_dir

        print('saving hybrid db as {}...'.format(hdb_object_file_name))
        with open(hdb_object_file_name, 'wb') as fh:
            pickle.dump(hdb, fh)

    return hdb, hdb_object_file_name


def prepare_calculations(
        hdb,
        MC_par_dic,
        double_counting_method=None,
        results_file_path=None,
        directory_baci_price_data=None,
        price_data_feather_file=None
        ):

    #Check if the necessary arguments have been given:
    arguments_not_passed = []
    if double_counting_method       is None: arguments_not_passed.append('double_counting_method')
    if directory_baci_price_data    is None: arguments_not_passed.append('directory_baci_price_data')
    if price_data_feather_file      is None: arguments_not_passed.append('price_data_feather_file')
    if arguments_not_passed:
        raise NameError('please provide the following arguments: \n {}'.format(
                        '\n'.join(arguments_not_passed)))

        # If no results file name is passed, create a standard name and path based on parameters
    if results_file_path is None:
        results_file_name = 'Impacts_Hybrid_MC_ecoinvent_{}_exiobase_{}_double_counting_correction_{}.hdf5'.format(
                hdb.ecoinvent_version, hdb.exiobase_version,double_counting_method)
        results_file_path = os.path.join(hdb.Database_dir,results_file_name)
        print('Results will be saved in {}'.format(results_file_path))

    # Checking if valid results file name
    if os.path.exists(results_file_path):
        print('You are about to overwrite the file {}'.format(results_file_path))
        valid_input = False
        while not valid_input:
            user_input = input("Do you want to do this? Y/N  ")
            if user_input == 'Y' or user_input == 'y':
                valid_input = True
                print('Will overwrite file {}'.format(results_file_path))
                os.remove(results_file_path)
                continue
            elif user_input == 'N' or user_input == 'n':
                print('Please provide a different file name and run again!')
                return 0


    # Set double counting correction method before starting calculations:
    hdb.double_counting_method = double_counting_method


    # in case that the regions and prices vary build a dictionary with the price
    # distributions per region. Else build a fixed geography cut-off matrix
    if not MC_par_dic['constant_region_var_price'] and not isinstance(hdb.price_dict, dict):
        print('Building price dictionary...')
        hdb.build_price_dictionary(path_baci_price_data=directory_baci_price_data,
                               price_data_df_file=price_data_feather_file)
    elif not isinstance(hdb.priceless_Cu, scipy.sparse.csr.csr_matrix):
        print('Loading price Data Frame from feather file...')
        hdb.load_price_array(price_data_df_file=price_data_feather_file)
        print('Calculating priceless Cu...')
        hdb.build_priceless_Cu()

    print('Initialising the results hdf5 file...')
    f = h5py.File(results_file_path, 'w')
    f.attrs['date created'] = time.ctime()
    f.attrs['Name'] = 'Hybrid (MRIO) impacts for ecoinvent {} and EXIOBASE {}'.format(hdb.ecoinvent_version, hdb.exiobase_version)
    f.attrs['double counting method'] = hdb.double_counting_method
    if MC_par_dic['constant_region_var_price']:
        f.attrs['Price Data'] = 'As modelled in DOI: 10.3389/frsus.2021.666209, fixed geographies'
    elif MC_par_dic['use_ecoinvent_price']:
        f.attrs['Price data'] = 'Fixed Ecoinvent Price'
    elif MC_par_dic['constant_price']:
        f.attrs['Price data'] = 'Median prices per region taken from BACI where possible,\
                or median price from distributions moddeled as described in DOI: 10.3389/frsus.2021.666209'
    else:
        f.attrs['Price data'] = 'BACI, where possible, rest modelled as described in DOI: 10.3389/frsus.2021.666209'


    impacts = f.create_group('io_impacts')
    for index in MC_par_dic['indices']:
        impacts.create_dataset(hdb.impact_names_io[index],
                               shape=(MC_par_dic['N_runs'],len(hdb.PRO)),
                               dtype='f8')
    f.close()

def run_monte_carlo(hdb, MC_par_dic, results_file_path, N_cores, hdb_object_file_name):
    if MC_par_dic['constant_region_var_price']:
        #don't use ray multiplrocessing
        t01 = time.time()
        _ = [MC_footprint(hdb, MC_par_dic['indices'], results_file_path,
                          constant_price=MC_par_dic['constant_price'],
                          use_ecoinvent_price=MC_par_dic['use_ecoinvent_price'],
                          constant_region_var_price=MC_par_dic['constant_region_var_price'],
                          FU=MC_par_dic['FU'], run_nr=i, total_runs=MC_par_dic['N_runs']) for i in range(MC_par_dic['N_runs'])]
        print('Done calculating results')
        dt1 = time.time()-t01
        print("Elapsed wall time for Calculating {} runs: {} hours, {}, minutes, and {:.0f} seconds.".format(
                N_runs, dt1//3600, (dt1%3600)//60, dt1%60))
    else:
        # initiate ray
        print('initialising Ray...')
        ray.init(num_cpus=N_cores, _memory=5000*1024*1024,
                object_store_memory= 10000 * 1024 * 1024)
        print('Serialising hybrid database...')
        hdb_id = ray.put(hdb)

        t01 = time.time()
        print('Calculating results...')
        dummy_ids = [MC_footprint_ray.remote(hdb_id, MC_par_dic['indices'], results_file_path,
            constant_price=MC_par_dic['constant_price'], use_ecoinvent_price=MC_par_dic['use_ecoinvent_price'],
            constant_region_var_price=MC_par_dic['constant_region_var_price'], FU=MC_par_dic['FU'], run_nr=i, total_runs=MC_par_dic['N_runs']) for i in range(MC_par_dic['N_runs'])]
        dummy_resulst = ray.get(dummy_ids)
        del dummy_resulst
        print('Done calculating results')
        dt1 = time.time()-t01
        print("Elapsed wall time for Calculating {} runs: {} hours, {}, minutes, and {:.0f} seconds.".format(
                N_runs, dt1//3600, (dt1%3600)//60, dt1%60))
        print('Shut down Ray...')
        ray.shutdown()

    if hdb.calc_lca_impacts_flag:
        print('Calculating LCA impacts...')
        hdb.calc_lca_impacts(FU=MC_par_dic['FU'])
        del hdb.calc_lca_impacts_flag

        print('Saving lca impacts to file')
        f = h5py.File(results_file_path, 'a')
        lca_impacts = f.create_dataset('lca_impacts', data=hdb.lca_impacts)
        f.close()

    # save the results file path in the hdb object
    # (specify dc method to allow for multiple results files using same hdb obj)
    file_name_string = 'MC_results_file_{}'.format(hdb.double_counting_method)
    setattr(hdb, file_name_string, results_file_path)
    time_stamp_string = 'time_stamp_{}'.format(hdb.double_counting_method)
    setattr(hdb, time_stamp_string, time.ctime())
    print('saving hybrid db as {}...'.format(hdb_object_file_name))
    with open(hdb_object_file_name, 'wb') as fh:
        pickle.dump(hdb, fh)

    f = h5py.File(results_file_name, 'a')
    f.attrs['hdb_file_name'] = hdb_object_file_name
    f.close()

    # Write Meta data to text file with the same name as results file
    meta_data_file_path = results_file_path.replace('.hdf5','.txt')
    print('writing meta data to {}'.format(meta_data_file_path))
    with open(meta_data_file_path, 'w') as fh:
        fh.write('Results file generated at : {}\n'.format(time.ctime()))
        fh.write('HDB Object file name : {}\n'.format(hdb_object_file_name))
        fh.write('Monte Carlo parameters :\n')
        for par,value in MC_par_dic.items():
            fh.write('{} : {}\n'.format(par,value))
        fh.write('double counting method : {}\n'.format(hdb.double_counting_method))
        fh.write('ecoinvent version : {}'.format(hdb.ecoinvent_version))
        fh.write('exiobase version : {}'.format(hdb.exiobase_version))
        fh.write('Calculated impacts:\n')
        for index in MC_par_dic['indices']:
            fh.write('Index {} : {}\n'.format(index,hdb.impact_names_io[index]))
 #   '''

#    result_ids = [MC_footprint(hdb, indices,
#                 constant_price=constant_price, FU=FU) for i in range(10)]
    print('Done! {}'.format(time.ctime()))
    return 0

@ray.remote
def MC_footprint_ray(obj, indices, file_name,
                 constant_price=False,
                 use_ecoinvent_price=False,
                 constant_region_var_price=False,
                 FU=None, run_nr=None, total_runs=None):
    if run_nr%100==0:
        print('Run {} of {}. {}'.format(run_nr,total_runs, time.ctime()))
    Cu = obj.create_Cu(constant_price=constant_price,
                       use_ecoinvent_price=use_ecoinvent_price,
                       constant_region_var_price=constant_region_var_price)
    if FU is not None:
        impacts = (obj.M_io[indices,:].dot(Cu).dot(FU)).A
    else:
        impacts = (obj.M_io[indices,:].dot(Cu)).A
    with FileLock(file_name+'.lock'):
            f = h5py.File(file_name, 'a')
            for i,index in enumerate(indices):
                f['io_impacts/{}'.format(obj.impact_names_io[index]
                                         )][run_nr,:] = impacts[i,:]
            f.close()
    return 0

def MC_footprint(obj, indices, file_name,
                 constant_price=False,
                 use_ecoinvent_price=False,
                 constant_region_var_price=False,
                 FU=None, run_nr=None, total_runs=None):
    if run_nr%100==0:
        print('Run {} of {}. {}'.format(run_nr,total_runs, time.ctime()))
    Cu = obj.create_Cu(constant_price=constant_price,
                       use_ecoinvent_price=use_ecoinvent_price,
                       constant_region_var_price=constant_region_var_price)
    if FU is not None:
        impacts = (obj.M_io[indices,:].dot(Cu).dot(FU)).A
    else:
        impacts = (obj.M_io[indices,:].dot(Cu)).A
    with FileLock(file_name+'.lock'):
            f = h5py.File(file_name, 'a')
            for i,index in enumerate(indices):
                f['io_impacts/{}'.format(obj.impact_names_io[index]
                                         )][run_nr,:] = impacts[i,:]
            f.close()
    return 0



# if __name__=="__main__":
#     print('Starting at: {}'.format(time.ctime()))
#     #Save time for overall walltime
#     t0 = time.time()
#     indices=[5,22,31,34,35]
#     main(N_runs=N_runs, indices=indices, constant_price=constant_price,
#             use_ecoinvent_price=use_ecoinvent_price,
#             constant_region_var_price=constant_region_var_price, FU=None,
#             double_counting_method)
#     dt = time.time()-t0
#     print("Total elapsed wall time: {} hours, {}, minutes, and {:.0f} seconds.".format(
#             dt//3600, (dt%3600)//60, dt%60))
#     print('Current time {}'.format(time.ctime()))
