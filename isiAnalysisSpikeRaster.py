################################################################################################
### Script to analyze the inter-spike intervals from spike raster data. Requires SciPy 1.16. ###
################################################################################################

### Copyright 2022-2024 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
### Contact: mail[at]jlubo.net

import numpy as np
import pandas as pd
from scipy import pi
import matplotlib.pyplot as plt
from pylab import figure
import os
from pathlib import Path
import json

# extractIntitalPeriodSpikes
# Extracts the distribution (histogram) of interspike intervals from a given period
# of spike raster data
# simulator: indicates the simulator from which the data stems ("Stand-alone", "Arbor", or "Brian")
# dt: duration of one timestep (in s)
# t_max: initial period to be considered (in s)
# return: two arrays, the first containing the (mean) ISI value for each bin, 
#         and the second containing the relative frequency of this value
def extractIsiDist(simulator, dt, t_max):


	if simulator == "Stand-alone":
		file_suffix = "_spike_raster.txt" # suffix of the filename for spike raster data files
		col_sep = "\t\t" # separator for data columns
	elif simulator == "Arbor" or simulator == "Brian":
		file_suffix = "_spikes.txt"
		col_sep = " "
	else:
		raise ValueError("Unknown simulator: '" + simulator + "'.") 

	# search in current directory for spike raster data files and use them all
	rawpaths = Path(".")
	df = None
	for x in sorted(rawpaths.iterdir()):

		full_path = str(x)
		tpath = os.path.split(full_path)[1] # take tail
		if file_suffix in tpath:
			print("Reading", tpath)
			df_new = pd.read_table(tpath, header=None, sep=col_sep, engine='python') #, nrows=59
			if df is None:
				df = df_new
			else:
				df = df.append(df_new)

	if df is None:
		print("No data for", simulator, "found. Exiting...")
		exit()

	# unit conversion
	if simulator == "Arbor" or simulator == "Brian":
		df[df.columns[0]] /= 1000 # convert ms to s

	# extract spikes of neuron 0 in the intital period (defined by t_max)
	df_initial_period = df[np.logical_and(df[df.columns[0]] < t_max, df[df.columns[1]] == 0)]
	#df = df_initial_period
	print("Spikes of neuron 0 in the first", t_max, "s:")
	print("  " + simulator + ":", len(df_initial_period))

	spike_times = df_initial_period[df_initial_period.columns[0]].to_numpy()
	isis = spike_times[1:] - spike_times[0:-1]
	print(isis)

	dec = -int(np.around(np.log10(dt)))
	print("Relevant decimal places:", dec)
	#bins=np.linspace(np.around(np.min(isis), dec)-5*dt,np.around(np.max(isis), dec)+5*dt, 10)
	isis_bin_edges = np.arange(np.around(np.min(isis)-2.5*dt, dec),np.around(np.max(isis)+2.5*dt, dec), dt)

	#isis_hist, isis_bin_edges = np.histogram(isis, bins=10, range=None, normed=None)
	isis_hist, isis_bin_edges = np.histogram(isis, bins=isis_bin_edges, normed=None)

	isis_bin_size = isis_bin_edges[1] - isis_bin_edges[0]
	isis_bins = isis_bin_edges[0:-1] + isis_bin_size/2
	
	print("Edges:", isis_bin_edges)
	print("Mean of bins:", isis_bins)
	return isis_bins, isis_bin_edges, isis_hist


# Find and read config file
timestamp_arbor, timestamp_standalone = None, None
rawpaths = Path(".")
config = None
for path in sorted(rawpaths.iterdir()):
	tpath = os.path.split(str(path))[1] # take tail
	if not path.is_dir():
		if "_config.json" in tpath:
			config = json.load(open(tpath, "r"))
if config is None:
	raise FileNotFoundError("No config file found.")

# Set constants
dt = config["simulation"]["dt"] # duration of one time bin (in ms)
t_max = 0.1 #config["simulation"]["runtime"] # initial period to be considered (in s)
t_ref = config["neuron"]["t_ref"] # duration of the refractory period (in ms)

# Extract ISI data
isis_bins_arbor, isis_bin_edges_arbor, isis_hist_arbor = extractIsiDist("Arbor", dt/1000, t_max)
isis_bins_standalone, isis_bin_edges_standalone, isis_hist_standalone = extractIsiDist("Stand-alone", dt/1000, t_max)

# Create histogram plot for Arbor data
figure(figsize=(8, 6), dpi=150)
ax1 = plt.subplot(211)
plt.title('Arbor')
plt.xlabel('Inter-spike interval (ms)')
plt.ylabel('Relative frequency')
plt.xticks(isis_bins_arbor*1000)
#plt.hist(isis_hist_arbor, bins=isis_bin_edges_arbor*1000)
#plt.step(isis_bins_arbor*1000, isis_hist_arbor, "-")
plt.stairs(isis_hist_arbor, isis_bin_edges_arbor*1000)
#inv = ax1.transData.inverted()
#x, y1 = inv.transform((1040, 5300))
#_, y2 = inv.transform((1040, 4300))
#_, y3 = inv.transform((1040, 3300))
#ax1.text(x, y1, "t_max = " + "{:.0f}".format(t_max*1000) + " ms")
#ax1.text(x, y2, "t_ref = " + "{:.3f}".format(t_ref) + " ms")
#ax1.text(x, y3, "dt = " + "{:.3f}".format(dt) + " ms")
ax1.annotate("t_max = " + "{:.0f}".format(t_max*1000) + " ms", (1260, 1050), xycoords='figure pixels')
ax1.annotate("t_ref = " + "{:.3f}".format(t_ref) + " ms", (1260, 1000), xycoords='figure pixels')
ax1.annotate("dt = " + "{:.3f}".format(dt) + " ms", (1260, 950), xycoords='figure pixels')

# Create histogram plot for Stand-alone data
dist_overlap = max(np.amin(isis_bins_arbor), np.amin(isis_bins_standalone)) < min(np.amax(isis_bins_arbor), np.amax(isis_bins_standalone)) # tests if distributions overlap
if dist_overlap:
	ax2 = plt.subplot(212, sharex=ax1)
	plt.xticks(np.concatenate((isis_bins_arbor, isis_bins_standalone))*1000)
	print("Distributions overlap. Using shared x-axis.")
else: # separate axes else
	ax2 = plt.subplot(212)
	plt.xticks(isis_bins_standalone*1000)
	print("Distributions do not overlap. Using separate x-axes.")

plt.title('Stand-alone')
plt.xlabel('Inter-spike interval (ms)')
plt.ylabel('Relative frequency')
#plt.hist(isis_hist_standalone, bins=isis_bin_edges_standalone*1000)
#plt.step(isis_bins_standalone*1000, isis_hist_standalone, "-")
plt.stairs(isis_hist_standalone, isis_bin_edges_standalone*1000)

#plt.subplot_tool()
plt.tight_layout()
plt.savefig("isi_comparison_neuron_0.png", bbox_inches='tight', dpi=300)
#plt.show()
plt.close()
