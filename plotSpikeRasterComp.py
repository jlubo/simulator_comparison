#################################################################
### Script to plot a comparison of two spike raster datasets. ###
#################################################################

### Copyright 2022-2024 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
### Contact: mail[at]jlubo.net

import numpy as np
import matplotlib.pyplot as plt
import json

#####################################
# plotRaster
# Plots a spike raster, highlighting different neuronal subpopulations by different colors
# config: configuration parameters in JSON format
# spikes_stacked_standalone: two-dimensional array containing the neuron and time of each spike (from stand-alone simulation)
# spikes_stacked_arbor: two-dimensional array containing the neuron and time of each spike (from Arbor simulation)
# store_filepath [optional]: path to the graphics file to be stored
# figure_fmt [optional]: format of resulting graphics file
def plotRaster(config, spikes_stacked_standalone, spikes_stacked_arbor, store_filepath = "merged_spike_raster.png"):
	N_CA = int(config["populations"]["N_CA"])
	N_exc = int(config["populations"]["N_exc"]) 
	N_inh = int(config["populations"]["N_inh"])
	N = N_exc + N_inh
	
	#plt.axhspan(0, N_CA, color='k', alpha=0.5) # background color
	#plt.axhspan(N_exc, N, color='r', alpha=0.5) # background color
	
	mask_CA_standalone = (spikes_stacked_standalone[1] < N_CA)
	mask_exc_standalone = np.logical_and(spikes_stacked_standalone[1] >= N_CA, spikes_stacked_standalone[1] < N_exc)
	mask_inh_standalone = (spikes_stacked_standalone[1] >= N_exc)

	mask_CA_arbor = (spikes_stacked_arbor[1] < N_CA)
	mask_exc_arbor = np.logical_and(spikes_stacked_arbor[1] >= N_CA, spikes_stacked_arbor[1] < N_exc)
	mask_inh_arbor = (spikes_stacked_arbor[1] >= N_exc)
	
	marker_type = '.' # ',' 
	marker_size = 1
	plt.plot(spikes_stacked_standalone[0][mask_CA_standalone], spikes_stacked_standalone[1][mask_CA_standalone], marker_type, label="Stand-alone", color='blue', alpha=0.4, markersize=marker_size)
	plt.plot(spikes_stacked_standalone[0][mask_exc_standalone], spikes_stacked_standalone[1][mask_exc_standalone], marker_type, color='blue', alpha=0.4, markersize=marker_size)
	plt.plot(spikes_stacked_standalone[0][mask_inh_standalone], spikes_stacked_standalone[1][mask_inh_standalone], marker_type, color='blue', alpha=0.4, markersize=marker_size)
	
	plt.plot(spikes_stacked_arbor[0][mask_CA_arbor], spikes_stacked_arbor[1][mask_CA_arbor], marker_type, label="Arbor", color='red', alpha=0.4, markersize=marker_size)
	plt.plot(spikes_stacked_arbor[0][mask_exc_arbor], spikes_stacked_arbor[1][mask_exc_arbor], marker_type, color='red', alpha=0.4, markersize=marker_size)
	plt.plot(spikes_stacked_arbor[0][mask_inh_arbor], spikes_stacked_arbor[1][mask_inh_arbor], marker_type, color='red', alpha=0.4, markersize=marker_size)
	
	print("Spikes in the first 500 ms (stand-alone sim.):", len(spikes_stacked_standalone[0][spikes_stacked_standalone[0] < 500]))
	print("Spikes in the first 500 ms (Arbor sim.):", len(spikes_stacked_arbor[0][spikes_stacked_arbor[0] < 500]))

	for ts in spikes_stacked_standalone[0][spikes_stacked_standalone[1] == 3]: # draw vertical lines for spikes of neuron #4 in stand-alone simulation
		plt.axvline(ts, ls='dotted', c="purple", lw=0.7, alpha=0.4)
	
	plt.xlim(0, 0.2)
	plt.ylim(-0.5, N-0.5)
	plt.xlabel('Time (s)', fontsize=12)
	plt.ylabel('Neuron index', fontsize=12)
	#ylabels=np.linspace(0, N, num=5, endpoint=False)
	#plt.yticks(ylabels, ylabels, rotation='horizontal')
	plt.yticks([0,1,2,3])
	plt.legend()
	plt.tight_layout()
	plt.savefig(store_filepath, dpi=800)
	plt.close()
	
#####################################
if __name__ == '__main__':

	timestamp_standalone = "22-07-05_14-16-22"
	timestamp_arbor = "22-07-04_16-35-12"

	spikes_stacked_standalone = np.loadtxt(timestamp_standalone + "_reformatted_spikes.txt").transpose()
	spikes_stacked_arbor = np.loadtxt(timestamp_arbor + "_reformatted_spikes.txt").transpose()
	config = json.load(open(timestamp_arbor + "_config.json", "r"))

	if spikes_stacked_standalone.size != 0 and spikes_stacked_arbor.size != 0:
		plotRaster(config, spikes_stacked_standalone, spikes_stacked_arbor, store_filepath = "merged_spike_raster.png")
	else:
		print("Missing spike data for plotting...")
