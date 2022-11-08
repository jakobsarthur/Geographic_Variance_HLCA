import pandas as pd
import numpy as np

import sys
# sys.path.append('/home/jakobs/Documents/IndEcol/OASES/ecospold2matrix/')
# sys.path.append('/home/jakobs/Documents/IndEcol/OASES/pymrio/')
# sys.path.append('/home/jakobs/Documents/IndEcol/OASES/pylcaio/src/')
# import ecospold2matrix as e2m
# import pymrio
# import pylcaio


import os
import gzip
import pickle
import pandas as pd

import pdb
import scipy.sparse
from pypardiso import spsolve, factorized
import matplotlib.pyplot as plt
import scipy.stats as stats
import time
import numba
import json
import ray
import random



class Hybrid_LCAIO(object):


    def __init__(self, lcaio_object, double_counting_method):
        self.A_lca = lcaio_object.A_ff
        self.F_lca = lcaio_object.F_f
        self.C_lca = lcaio_object.C_f
        self.A_io = lcaio_object.A_io
        self.F_io = lcaio_object.F_io
        self.C_io = None


        self.PRO = lcaio_object.PRO_f
        self.dictRoW = lcaio_object.dictRoW
        self.countries_per_region = lcaio_object.countries_per_regions
        self.regions_of_IO = lcaio_object.regions_of_IO
        self.hybridized_processes = lcaio_object.hybridized_processes
        self.list_not_to_hyb = lcaio_object.list_not_to_hyb
        self.hybridized_processes_boolean_mask = self.PRO.index.isin(
                                                      self.hybridized_processes)

        self.region_weights_df = (lcaio_object.H*lcaio_object.Geo).sum(
                                                                 level='region')
        self.H = lcaio_object.H
        #Set concordance to zero for processes not to hybridize
        self.H.loc[:,self.list_not_to_hyb] = 0
        self.H = scipy.sparse.csr_matrix(self.H)
        self.Fixed_geo_concordance = self.H.multiply(lcaio_object.Geo).tocsr() 

        self.double_counting = lcaio_object.double_counting_method
        if self.double_counting == 'STAM':
            self.correction_matrix_STAM = lcaio_object.correct_STAM
        elif self.double_counting == 'binary':
            self.correction_matrix_binary = lcaio_object.correct_binary

        self.sectors_of_IO = lcaio_object.sectors_of_IO
        self.regions_of_IO = lcaio_object.regions_of_IO

        self.Nsec = len(self.sectors_of_IO)
        self.Nreg = len(self.regions_of_IO)
        self.NregionSecotrs_io = self.Nsec*self.Nreg
        self.make_io_region_index_dict()
        self.make_country_industry_IO_mapping()

        self.price_dict = None
        self.price_data_array = None
        self.priceless_Cu = None

        del lcaio_object

    def make_io_region_index_dict(self):
        self.io_region_index_mapping = {}
        for i,region in enumerate(self.regions_of_IO):
            indices = np.arange(i*self.Nsec,i*self.Nsec+self.Nsec)
            self.io_region_index_mapping[region] = indices

    def make_country_industry_IO_mapping(self):
        self.country_industry_IO_mapping = {}
        for i,j in enumerate(pd.MultiIndex.from_product([self.regions_of_IO,
                                        self.sectors_of_IO]).to_flat_index()):
            self.country_industry_IO_mapping[j] = i


    def generate_country_vector(self):
        '''
        This method is obsolete as it does not include the production volumes
        of the IO industries/countries'''
        # copy IO region vec (pd series)
        new_region_vec = self.PRO.io_geography.copy()
        if not hasattr(self, 'hybrid_and_region_mask'):
            self.hybridized_and_region_mask = np.logical_and(
                    ~new_region_vec.isin(self.regions_of_IO),
                    new_region_vec.index.isin(self.hybridized_processes))
        empty_RoW = []
        for index in new_region_vec[self.hybridized_and_region_mask].index:
            if new_region_vec.loc[index] in self.countries_per_region.keys():
                new_region_vec.loc[index] = random.choice(
                                            self.countries_per_region[
                                            new_region_vec.loc[index]])
            else:
                try:
                    new_region_vec.loc[index] = random.choice(
                                                self.dictRoW[
                                                new_region_vec.loc[index]])
                except ValueError:
                    print('Undefined RoW region for: ',
                            index, new_region_vec.loc[index])
                    empty_RoW.append(index)
                    continue
        return new_region_vec

    def build_price_dictionary(self, path_baci_price_data=None,
                                     price_data_df_file=None):
        """
        Method builds a dictionary containing all price data for all
        hybridised processes, ordered by country (for aggregate region processes
        ) and using BACI price data for those that have it available and
        currently still a fixed price distribution for those without BACI prices
        . These fixed price distributions are taken from the paper
        10.3389/frsus.2021.666209 .

        path_baci_price_data       directory contining the pickle files with
                                   baci price pickle files
        price_data_df_file         feather file with price distributions for
                                   those processes without a BACI price.
        """
        self.price_data_df = pd.read_feather(price_data_df_file).set_index('index')
        self.price_dict = {}
        self.processes_with_baci_price = []
        for row in self.PRO.loc[self.hybridized_processes].itertuples():
            index = row.Index
            try:
                with open(os.path.join(path_baci_price_data,
                                       '{}.pickle'.format(index)), 'rb') as fh:
                    price_dic = pickle.load(fh)
                    self.price_dict[index] = price_dic
                region_weights = self.region_weights_df[index]
                self.price_dict[index]['regions'] = self.region_weights_df.index.to_numpy()
                available_regions = set(self.price_dict[index]['pricesByCountry'].keys())
                regions_to_set_to_zero = set(self.price_dict[index]['regions']) - available_regions
                region_weights.loc[regions_to_set_to_zero] = 0
                self.price_dict[index]['region_weights'] = region_weights.to_numpy()/region_weights.sum()
                self.processes_with_baci_price.append(index)

            except FileNotFoundError:
                self.price_dict[index] = {}
                self.price_dict[index]['pricesByCountry'] = {}
                region = self.PRO.loc[index,'io_geography']
                if not region in self.regions_of_IO:
                    if region in self.countries_per_region.keys():
                        countries = self.countries_per_region[region]
                    else:
                        countries = self.dictRoW[region]
                else:
                    countries = [region]

                for country in countries:
                    self.price_dict[index]['pricesByCountry'][country] = {}
                    self.price_dict[index]['pricesByCountry'][country]['price_sample'] = self.price_data_df.loc[index].to_numpy()
                self.price_dict[index]['region_weights'] = self.region_weights_df[index].to_numpy()
                self.price_dict[index]['regions'] = self.region_weights_df.index.to_numpy()

    def get_region_and_price(self, constant_price=False, use_ecoinvent_price=False):
        if self.price_dict==None:
            raise AttributeError("Please initiate the price dict by\
                    by running the method\
                    build_price_dictionary(path_baci_price_data=None,\
                    price_data_df_file=None. Else run create_Cu with\
                    use_ecoinvent_price=True)")
        new_region_vec = self.PRO.io_geography.copy()
        new_price_vec = self.PRO.price.copy()
        for row in self.PRO.loc[self.hybridized_processes].itertuples():
            index = row.Index
            region = row.io_geography
            if region in self.regions_of_IO:
                new_region = region
            else:
                try:
                    new_region = np.random.choice(
                            self.price_dict[index]['regions'],
                            p=self.price_dict[index]['region_weights'])
                except ValueError:
                    # There are a few processes with zero production volumes
                    # (these are empty industries) (2 cases)
                    # just pick region from the available price data
                    new_region = np.random.choice(
                            [*self.price_dict[index]['pricesByCountry']])
            new_region_vec.loc[index] = new_region
            if not use_ecoinvent_price:
                try:
                    if not constant_price:
                        # use a random price from the price sample
                        # note that this already includes the price distribution
                        new_price = np.random.choice(self.price_dict[index][
                                    'pricesByCountry'][new_region]['price_sample'])
                    else:
                        # use the mean price of the price sample (but not
                        # region specific!)
                        # (that means BACI price  if available and ecoinvent
                        # price otherwise)
                        new_price = np.mean(self.price_data_df.loc[index])
                    new_price_vec.loc[index] = new_price
                except KeyError:
                    print(index, region)
                    pdb.set_trace()
        if use_ecoinvent_price:
            # just use the ecoinvent prices
            new_price_vec = self.PRO.price.copy()
        return new_price_vec, new_region_vec

    def get_variable_price(self):
        if not isinstance(self.price_data_array, np.ndarray):
            raise AttributeError("Please load the price DataFrame by running the\
                    method load_price_df(price_data_df_file).")
        shape = self.price_data_array.shape
        return self.price_data_array[np.arange(shape[0]),
                                     np.random.randint(shape[1],size=shape[0])]

    def load_price_array(self, price_data_df_file):
        price_data_df = pd.read_feather(price_data_df_file).set_index('index')
        self.price_data_array = price_data_df.to_numpy()

    def get_column_indices(self, region_vec):
        column_indices = [self.country_industry_IO_mapping[x] for x in zip(
            region_vec.loc[self.hybridized_processes_boolean_mask],
            self.PRO.loc[self.hybridized_processes_boolean_mask,
                                                            'ProductTypeName'])]
        return column_indices

    def build_uncorrected_priceless_Cu(self, region_vec):
        column_indices = self.get_column_indices(region_vec)
        Cu = np.zeros((self.NregionSecotrs_io, len(self.PRO)))
        Cu[:,self.hybridized_processes_boolean_mask] = self.A_io[:,
                                                       column_indices].toarray()
        return Cu

    def build_priceless_Cu(self):
        uncorrected_Cu = self.A_io.dot(self.Fixed_geo_concordance)
        self.priceless_Cu = self.correction_matrix_STAM.multiply(uncorrected_Cu).tocsr()


    def create_Cu(self, constant_price=False, use_ecoinvent_price=False,
                  constant_region_var_price=False):
        """Note that the use_ecoinvent_price flag overrules the
        constant_price flag.

        constant_price flag False means use price samples to draw from

        constant_price flag True takes the mean of the price distributions

        use_ecoinvent_price flag overrules both options and uses the fixed
        prices from ecoinvent

        constant_region_var_price flag overrules all the above. It keeps the
        regions fix as per pylcaio (weighted average of production volumes) and
        varies the price according to a predefined dataframe also used in
        build_price_dictionary. For this option the method load_price_array must be
        run first.
        """
        
        if constant_region_var_price:
            if self.priceless_Cu==None:
                print("Building priceless Cu...")
                self.build_priceless_Cu()
            price_vec = self.get_variable_price()
            Cu = self.priceless_Cu.multiply(price_vec).tocsr()
        else:
            price_vec, region_vec = self.get_region_and_price(
                                        constant_price=constant_price,
                                        use_ecoinvent_price=use_ecoinvent_price)
            Cu = self.build_uncorrected_priceless_Cu(region_vec)
            Cu *= price_vec
            Cu = self.correction_matrix_STAM.multiply(Cu).tocsr()
        return Cu

    def calc_L_lca(self):
        A_lca = scipy.sparse.eye(self.A_lca.shape[0])-self.A_lca
        self.L_lca = scipy.sparse.csr_matrix(scipy.sparse.linalg.spsolve(A_lca,
                                                np.eye(A_lca.shape[0])))

    def calc_lca_impacts(self, FU=None):
        if FU is not  None:
            self.lca_impacts = (self.C_lca.dot(self.F_lca).dot(FU)).A
        else:
            self.lca_impacts = (self.C_lca.dot(self.F_lca).dot(self.L_lca)).A

    def calc_M_io(self, C_io, impact_names_io):
        self.C_io = C_io
        self.impact_names_io = impact_names_io
        L_io = scipy.sparse.csr_matrix(
                        scipy.linalg.solve(np.eye(self.A_io.shape[0])
                                -self.A_io.todense(),
                                np.eye(self.A_io.shape[0])))
        self.M_io = C_io.dot(self.F_io).dot(L_io)







