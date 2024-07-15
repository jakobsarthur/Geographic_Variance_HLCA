# Geographic_Variance_HLCA

[![DOI](https://zenodo.org/badge/455208079.svg)](https://zenodo.org/doi/10.5281/zenodo.12731487)



Repo to build hybrid lca database with ecoinvent en exiobase, including geographic and price variance. 
**Relevant Publication:** ***Accepted at Journal of Industrial Ecology***

The repository builds on a fork of the package [pylcaio[(https://github.com/MaximeAgez/pylcaio/tree/master). The relevant fork is found [here](https://github.com/OASES-project/pylcaio).
Otherwise it requires:

- Scipy
- Numpy
- [Ray distributed processing](https://github.com/ray-project/ray). 
- pypardiso
- [pylcaio](https://github.com/OASES-project/pylcaio)
- [ecospold2matrix](https://github.com/majeau-bettez/ecospold2matrix)
- pymrio
- feather

## File description:

### Script files

- **hybrid_lcaio.py:** The core of this package, a class that deals with the building of the hybrid model and (re-) calculating of the the Cut-off matrix. No need to use directly.
- **run_hybrid_MC.py:** Script to run a MC simulation to estimate the gepographic and price variance. Can be run in stand alone mode, but is called from the notebook ***run_hybrid_montecarlo.ipynb***
- **match_BACI_price_data_to_ecoinvent.py:** Script to match the BACI price daata to econivent processes. (price data for ecoinvent 3.5 and BACI for the year 2012 can be found here [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12744397.svg)](https://doi.org/10.5281/zenodo.12744397)).



### [Notebooks](https://github.com/jakobsarthur/Geographic_Variance_HLCA/tree/master/notebooks)

- **RunPylcio.ipynb:** Notebook to generate a pylcaio object (used for all the concordances).
- **run_hybrid_montecarlo.ipynb:** Notebook to run a MC simulation to estimate the gepographic and price variance
- **Analyse_MC)results.ipynb:** A working notebook that can be used to analyse the MC simulation data and produce plot. **WARNING:** not a cleaned notebook, **use with care**!!!

### [mapping_data](https://github.com/jakobsarthur/Geographic_Variance_HLCA/tree/master/mapping_data)
Mapping files between ecoinvent processes (v3.5 and v3.8) and BACI for products (HS12 classification) and regional mappings.
- eco_baci_country_mapping.xlsx : excel file mapping the regions in ecoinvent to country codes
- ecoinvent35-baci_region_mapping.json : json file mapping the activity_product uuids to a list of country codes
- ecoinvent35-HS12_mapping.json : json file mapping the activity_product uuids to a list of HS12 commodity codes
- ecoinvent38-baci_region_mapping.json : json file mapping the activity_product uuids to a list of country codes
- ecoinvent38-HS12_mapping.json : json file mapping the activity_product uuids to a list of HS12 commodity codes


