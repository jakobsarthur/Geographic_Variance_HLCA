import pandas as pd
import numpy as np

import sys
sys.path.append('/home/jakobs/Documents/IndEcol/OASES/ecospold2matrix/')
sys.path.append('/home/jakobs/Documents/IndEcol/OASES/pymrio/')
sys.path.append('/home/jakobs/Documents/IndEcol/OASES/pylcaio/src/')
import ecospold2matrix as e2m
import pymrio
import pylcaio

sys.path.append('/home/jakobs/Documents/IndEcol/OASES/Geographic_uncertainty/scripts')
from aggregate import build_MultiIndex_Aggregation_Matrix, Aggregation_dic_from_file

import hybrid_lcaio

import os
import gzip
import pickle
import pandas as pd

from importlib import reload
import pdb
import scipy.sparse
from pypardiso import spsolve, factorized
import matplotlib.pyplot as plt
import scipy.stats as stats
import time
import json
import ray
import random
import h5py
from filelock import FileLock

#number of runs
N_runs= 10000
constant_price = True
use_ecoinvent_price = False # this overrides the constant price flag
constant_region_var_price = False # this overrides the above two flags

N_cores = 2

lcaio_object_path = '/home/jakobs/Documents/IndEcol/OASES/Geographic_uncertainty/Databases/plcaio_object_STAM.pickle'
directory_baci_price_data = '/home/jakobs/Documents/IndEcol/OASES/Geographic_uncertainty/Databases/price_data_per_region/'
C_io_path = '/home/jakobs/Documents/IndEcol/Data/EXIOBASE/exiobase3_6/Characterization_EB36.xlsx'
price_data_feather_file = '/home/jakobs/Documents/IndEcol/OASES/pylcaio/src/Databases/ecoinvent3.5_exiobase3_2012/price_dataframe.ftr'
results_file_name = '/home/jakobs/Documents/IndEcol/OASES/Geographic_uncertainty/Databases/Results_geography_MC_BACI_price_{}_runs_3FPs_ProdVolWeights_constantprice_mean_{}_ecoinventprice_{}_constant_region_{}_only_Cu_impacts.hdf5'.format(N_runs, constant_price, use_ecoinvent_price, constant_region_var_price)
hdb_object_file_name = '/home/jakobs/Documents/IndEcol/OASES/Geographic_uncertainty/Databases/HDB_STAM.pickle'


def main(N_runs, indices, constant_price=False, use_ecoinvent_price=False,
         constant_region_var_price=False, FU=None):
    
    # Checking if valid results file name
    if os.path.exists(results_file_name):
        print('You are about to overwrite the file {}'.format(results_file_name))
        valid_input = False
        while not valid_input:
            user_input = input("Do you want to do this? Y/N  ")
            if user_input == 'Y' or user_input == 'y':
                valid_input = True
                print('Will overwrite file {}'.format(results_file_name))
                os.remove(results_file_name)
                continue
            elif user_input == 'N' or user_input == 'n':
                print('Please provide a different file name and run again!')
                return 0
    if os.path.exists(hdb_object_file_name):
        print('loading Hybrid Database Object...')
        with open(hdb_object_file_name, 'rb') as fh:
            hdb = pickle.load(fh)
        calc_lca_impacts_flag = False

    else:
        calc_lca_impacts_flag = True
        print("Reading in lcaio object, this takes a while...")
        with gzip.open(lcaio_object_path, 'rb') as fh:
            lcaio_object = pickle.load(fh)
    
        print("Creating hybrid database object...")
        hdb = hybrid_lcaio.Hybrid_LCAIO(lcaio_object)
        del lcaio_object
        
        print('Calculate L_lca, this takes a while...')
        hdb.calc_L_lca()

        print('Reading in IO Characterisation factors...')
        C_io = pd.read_excel(C_io_path, index_col=0, header=0)
        C_io = C_io.T
        Impact_names_io = C_io.index.str.rstrip().values
        C_io = scipy.sparse.csr_matrix(C_io)
        print('Calculate IO multiplyers...')
        hdb.calc_M_io(C_io, Impact_names_io)

    if not constant_region_var_price and not isinstance(hdb.price_dict, dict):
        print('Building price dictionary...')
        hdb.build_price_dictionary(path_baci_price_data=directory_baci_price_data,
                               price_data_df_file=price_data_feather_file)
    elif not isinstance(hdb.priceless_Cu, scipy.sparse.csr.csr_matrix):
        print('Loading price Data Frame from feather file...')
        hdb.load_price_array(price_data_df_file=price_data_feather_file)
        print('Calculating priceless Cu...')
        hdb.build_priceless_Cu()

    print('Initialising the results hdf5 file...')
    f = h5py.File(results_file_name, 'w')
    f.attrs['date created'] = time.ctime()
    f.attrs['Name'] = 'Hybrid (MRIO) impacts for ecoinvent 3.5 and EXIOBASE 3.6'
    if constant_region_var_price:
        f.attrs['Price Data'] = 'As modelled in DOI: 10.3389/frsus.2021.666209, fixed geographies'
    elif use_ecoinvent_price:
        f.attrs['Price data'] = 'Fixed Ecoinvent Price'
    elif constant_price:
        f.attrs['Price data'] = 'Median prices per region taken from BACI where possible,\
                or median price from distributions moddeled as described in DOI: 10.3389/frsus.2021.666209'
    else:
        f.attrs['Price data'] = 'BACI, where possible, rest modelled as described in DOI: 10.3389/frsus.2021.666209'

    # f.attrs['']
    impacts = f.create_group('io_impacts')
    for index in indices:
        impacts.create_dataset(hdb.impact_names_io[index],
                               shape=(N_runs,len(hdb.PRO)),
                               dtype='f8')
    f.close()

#    '''
    if constant_region_var_price:
        #don't use ray multiplrocessing
        t01 = time.time()
        _ = [MC_footprint(hdb, indices, results_file_name,
                          constant_price=constant_price,
                          use_ecoinvent_price=use_ecoinvent_price,
                          constant_region_var_price=constant_region_var_price,
                          FU=FU, run_nr=i) for i in range(N_runs)]
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
        dummy_ids = [MC_footprint_ray.remote(hdb_id, indices, results_file_name,
            constant_price=constant_price, use_ecoinvent_price=use_ecoinvent_price,
            constant_region_var_price=constant_region_var_price, FU=FU, run_nr=i) for i in range(N_runs)]
        dummy_resulst = ray.get(dummy_ids)
        del dummy_resulst
        print('Done calculating results')
        dt1 = time.time()-t01
        print("Elapsed wall time for Calculating {} runs: {} hours, {}, minutes, and {:.0f} seconds.".format(
                N_runs, dt1//3600, (dt1%3600)//60, dt1%60))
        print('Shut down Ray...')
        ray.shutdown()

    if calc_lca_impacts_flag==True:
        print('Calculating LCA impacts...')
        hdb.calc_lca_impacts(FU=FU)

    print('Saving lca impacts to file')
    f = h5py.File(results_file_name, 'a')
    lca_impacts = f.create_dataset('lca_impacts', data=hdb.lca_impacts)
    f.close()
    hdb.MC_results_file = results_file_name
    hdb.time_stamp = time.ctime()
    print('saving hybrid db as {}...'.format(hdb_object_file_name))
    with open(hdb_object_file_name, 'wb') as fh:
        pickle.dump(hdb, fh)

    f = h5py.File(results_file_name, 'a')
    f.attrs['hdb_file_name'] = hdb_object_file_name
    f.close()

 #   '''

#    result_ids = [MC_footprint(hdb, indices,
#                 constant_price=constant_price, FU=FU) for i in range(10)]
    print('Done!')
    return 0

@ray.remote
def MC_footprint_ray(obj, indices, file_name,
                 constant_price=False,
                 use_ecoinvent_price=False,
                 constant_region_var_price=False,
                 FU=None, run_nr=None):
    if run_nr%100==0:
        print('Run {} of {}. {}'.format(run_nr,N_runs, time.ctime()))
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
                 FU=None, run_nr=None):
    if run_nr%100==0:
        print('Run {} of {}. {}'.format(run_nr,N_runs, time.ctime()))
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


if __name__=="__main__":
    print('Starting at: {}'.format(time.ctime()))
    #Save time for overall walltime
    t0 = time.time()
    indices=[4,8,31]
    main(N_runs=N_runs, indices=indices, constant_price=constant_price,
            use_ecoinvent_price=use_ecoinvent_price,
            constant_region_var_price=constant_region_var_price, FU=None)
    dt = time.time()-t0
    print("Total elapsed wall time: {} hours, {}, minutes, and {:.0f} seconds.".format(
            dt//3600, (dt%3600)//60, dt%60))
    print('Current time {}'.format(time.ctime()))
