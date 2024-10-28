## Comparison of results across simulators

This is a collection of scripts serving to compare results from different spiking neural network simulators.

The scripts have been used, in particular, to compare the outcome of synapse plasticity simulations and neural network simulations run with:

* Arbor (see [this](https://github.com/jlubo/arbor_network_consolidation) repo),
* Brian 2 (see [this](https://github.com/jlubo/brian_synaptic_plasticity_stc) repo and [this](https://github.com/jlubo/brian_network_plasticity) repo),
* Stand-alone simulator for synaptic memory consolidation (see [this](https://github.com/jlubo/memory-consolidation-stc) repo).


### Overview

* __plotSimResultsComparisonMeanSEM.py__ - enables to compare the time course of specific quantities across simulators, in particular, showing the mean and standard error of the mean (SEM) of the quantities from batches of multiple trials
* __plotSpikeRasterComp.py__ - enables to compare spike raster data across simulators
* __spikeRasterReformat.py__ - reformats spike raster data for easier comparison
* __isiAnalysisSpikeRaster.py__ - extracts the distribution of inter-spike intervals from spike raster data
* __utilityFunctions.py__ - diverse utility functions needed by other modules (cf. [here](https://github.com/jlubo/memory-consolidation-stc/blob/main/analysis/utilityFunctions.py)) 

