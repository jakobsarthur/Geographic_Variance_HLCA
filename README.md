# Geographic_Variance_HLCA


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
- **match_BACI_price_data_to_ecoinvent.py:** Script to match the BACI price daata to econivent processes.

### [Notebooks](https://github.com/jakobsarthur/Geographic_Variance_HLCA/tree/master/notebooks)

- **RunPylcio.ipynb:** Notebook to generate a pylcaio object (used for all the concordances).
- - **run_hybrid_montecarlo.ipynb:** Notebook to run a MC simulation to estimate the gepographic and price variance
