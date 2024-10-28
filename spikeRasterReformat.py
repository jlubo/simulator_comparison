#################################################################################################################
### Script to reformat spike raster data from different simulators such that they can be compared more easily ###
#################################################################################################################

### Copyright (c) Jannik Luboeinski 2022-2024
### License: Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
### Contact: mail[at]jlubo.net

import numpy as np
import os
import time
import sys
from pathlib import Path
from shutil import copy2
from utilityFunctions import *

new_plots = True # defines if new spike raster plots shall be created using gnuplot
exc_pop_size = 4 # number of neurons in the excitatory population
inh_pop_size = 0 # number of neurons in the inhibitory population
core_size = 1 # number of neurons in the cell assembly
tot_size = exc_pop_size + inh_pop_size # number of neurons in the whole network
core = np.arange(core_size) # indices of the neurons in the cell assembly

####################################
# spikeRasterReformat
# Reformats spike raster data from different simulators such that they can be compared more easily. 
# To this end, shifts by the time of the first spike, and determines the spiking of neurons in defined binning
# periods. Furthermore, computes the number of spikes in each binning period computes and average firing rates
# of defined populations.
# timestamp: the timestamp of the simulation data
# spike_raster_file: the name of the spike raster file
# period_duration [optional]: binning period in units of seconds, default is 2e-6 s
# col_sep [optional]: specifies the column separator string, default is one space symbol
# time_unit [optional]: specifies the time unit of the given data, default is "ms"
# output_dir [optional]: relative path to the output directory, default is .
def spikeRasterReformat(timestamp, spike_raster_file, period_duration=2e-6, col_sep=' ', time_unit='ms',
                        output_dir='.'):

	t0 = time.time()
	epsilon = 0.1*period_duration # very small number used to avoid error with the floor function

	# read the first and the last line
	with open(spike_raster_file, 'rb') as f:
		first_line = f.readline().decode()
		f.seek(-2, os.SEEK_END)
		while f.read(1) != b'\n': # seek last line
			f.seek(-2, os.SEEK_CUR)
		last_line = f.readline().decode()
	if time_unit == "ms":
		time_factor = 1/1000
	elif time_unit == "s":
		time_factor = 1
	else:
		raise ValueError(f"Unsupported time unit '{time_unit}'.")
	# shift if the first spike time is less than or greater than 0
	t_first = float(first_line.split(col_sep)[0])
	if abs(t_first) > epsilon:
		t_shift = -t_first
		print("Shifting by t_first =", t_first, time_unit)
	else:
		t_shift = 0.0
	# compute number of periods to consider
	num_periods_tot = int(np.ceil((float(last_line.split(col_sep)[0]) + t_shift) * time_factor / period_duration)) + 1
	
	# count lines
	with open(spike_raster_file) as f:
		num_rows = sum(1 for _ in f)
	print("num_rows =", num_rows)

	# time and spikes of the last period
	last_period = -1
	last_spikes = []

	# counters per period for the different subpopulations
	counterCA = np.zeros(num_periods_tot, dtype=int)
	counterCtrl = np.zeros(num_periods_tot, dtype=int)
	counterInh = np.zeros(num_periods_tot, dtype=int)
	counterAll = np.zeros(num_periods_tot, dtype=int)
	
	# input of spike data
	f = open(spike_raster_file)

	# output of reformatted spike data
	fout = open(os.path.join(output_dir, timestamp + '_reformatted_spikes.txt'), 'w')

	# read, reformat, and write all data
	for line in f:
		
		# read out current line, shift and scale time
		row = line.split(col_sep)
		t = (float(row[0]) + t_shift) * time_factor

		n = int(row[1])

		current_period = int(np.floor((t + epsilon) / period_duration))

		# update counters
		counterAll[current_period] += 1
		if n in core:
			counterCA[current_period] += 1
		elif n < exc_pop_size:
			counterCtrl[current_period] += 1
		else:
			counterInh[current_period] += 1

		# new period
		if current_period > last_period:

			# write spikes of last period to file in a sorted manner
			for spiking_neuron in sorted(last_spikes):
				t_reformatted = round(last_period*period_duration, 4)
				fout.write(str(t_reformatted) + "\t\t" + str(spiking_neuron) + "\n")

			last_period = current_period
			last_spikes = []

		# add number of spiking neuron to list
		last_spikes.append(n)
		
	f.close()
	fout.close()
	

	# output of spikes per time bin
	fout_binned = open(os.path.join(output_dir, timestamp + '_spike_number.txt'), 'w')
	for i in range(num_periods_tot):
		fout_binned.write(str(round(i*period_duration, 4)) + "\t\t" + \
		                  str(counterCA[i]) + "\t\t" + str(counterCtrl[i]) + "\t\t" + str(counterInh[i]) + "\t\t" + str(counterAll[i]) + "\n")
	fout_binned.close()

	# determine and print the elapsed time
	time_el = round(time.time()-t0) # elapsed time in seconds
	time_el_str = "Elapsed time: "
	if time_el < 60:
		time_el_str += str(time_el) + " s"
	else:
		time_el_str += str(time_el // 60) + " m " + str(time_el % 60) + " s"
	print(time_el_str)

	# write population data (size, number of spikes, firing rates) to file
	fout = open(os.path.join(output_dir, timestamp + '_firing_rates.txt'), 'w')
	if core_size > 0:
		fout.write("nu(core) = " + str(np.sum(counterCA) / (num_periods_tot*period_duration) / core_size) + \
		           ", n_spikes(core) = " + str(np.sum(counterCA)) + "\n")
	else:
		fout.write("nu(core) = n/a, n_spikes(core) = n/a\n")
	if exc_pop_size - core_size > 0:
		fout.write("nu(control) = " + str(np.sum(counterCtrl) / (num_periods_tot*period_duration) / (exc_pop_size-core_size)) + \
		           ", n_spikes(control) = " + str(np.sum(counterCtrl)) + "\n")
	else:
		fout.write("nu(control) = n/a, n_spikes(control) = n/a\n")
	if inh_pop_size > 0:
		fout.write("nu(inh.) = " + str(np.sum(counterInh) / (num_periods_tot*period_duration) / (inh_pop_size)) + \
		           ", n_spikes(inh.) = " + str(np.sum(counterInh)) + "\n")
	else:
		fout.write("nu(inh.) = n/a, n_spikes(inh.) = n/a\n")
	fout.write("nu(all) = " + str(np.sum(counterAll) / (num_periods_tot*period_duration) / (tot_size)) + \
	           ", n_spikes(all) = " + str(np.sum(counterAll)) + "\n")
	fout.write("core_size = " + str(core_size) + "\n")
	fout.write("exc_pop_size = " + str(exc_pop_size) + "\n")
	fout.write("inh_pop_size = " + str(inh_pop_size) + "\n")
	fout.write("tot_size = " + str(tot_size) + "\n")
	fout.close()


######################################
# dirRecursion
# Walks recursively through a directory looking for spike raster data;
# if data are found, calls `spikeRasterReformat()`
# directory: the directory to consider
# output_dir: relative path to the output directory
# new_plots [optional]: specifies if new plots shall be created, defualt is True
def dirRecursion(directory, output_dir, new_plots = True):

	def copyParamsFile(hpath, timestamp, output_dir):
		params_file = os.path.join(hpath, timestamp + "_PARAMS.txt")
		json_file = os.path.join(hpath, timestamp + "_config.json")
		if os.path.exists(params_file):
			copy2(params_file, output_dir)
		elif os.path.exists(json_file):
			copy2(json_file, output_dir)
		else:
			print("Warning: " + hpath + ": no parameter file found.")

	rawpaths = Path(directory)

	print("Reading directory " + directory)
	rawpaths = Path(directory)

	for x in sorted(rawpaths.iterdir()):

		dest_file = ""

		full_path = str(x)
		hpath = os.path.split(full_path)[0] # take head
		tpath = os.path.split(full_path)[1] # take tail

		if not x.is_dir():

			# data from stand-alone implementation (https://github.com/jlubo/memory-consolidation-stc)
			if hasTimestamp(tpath) and "_spike_raster.txt" in tpath:
				timestamp = tpath.split("_spike_raster.txt")[0]
				spikeRasterReformat(timestamp, full_path, col_sep="\t\t", time_unit="s", output_dir=output_dir)
				copyParamsFile(hpath, timestamp, output_dir)

			# data from other implementation (e.g., https://github.com/jlubo/arbor_network_consolidation)
			elif hasTimestamp(tpath) and "_spikes.txt" in tpath:
				timestamp = tpath.split("_spikes.txt")[0]
				spikeRasterReformat(timestamp, full_path, col_sep=" ", time_unit="ms", output_dir=output_dir)
				copyParamsFile(hpath, timestamp, output_dir)

		else:
			if hasTimestamp(tpath):
				dirRecursion(directory + os.sep + tpath, output_dir)


###############################################
# main:

# set and create output directory
output_dir = "./reformatted"
#output_dir = "../" 

if not os.path.exists(output_dir):
	os.mkdir(output_dir)

print("Output directory:", output_dir)

# walk through directories and analyze data
dirRecursion('.', output_dir) 
