#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil
import os.path

##########################################################################################
# RETRIEVE USER INPUTS
##########################################################################################

#=========================================================================================
# create parser
#=========================================================================================
version_nb = "0.0.1"
parser = argparse.ArgumentParser(prog = 'botlzmann_inversion', usage='', add_help = False, formatter_class = argparse.RawDescriptionHelpFormatter, description =\
'''
**************************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@bioch.ox.ac.uk)
git: https://github.com/jhelie/boltzmann_inversion
**************************************************

This script the opposite of the log of densities. The densities in each column of the input
file are normalised with respect to themselves.

NB: for quick plotting with xmgrace first change NaN values in output file to 0.

[ USAGE ]

Option	      Default  	Description                    
-----------------------------------------------------
-f			: xvg file(s)
-o			: name of outptut file
-col		1	: ignore the first X column(s) in the input file
--comments	@,#	: lines starting with these characters will be considered as comment

Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
 
''')

#options
parser.add_argument('-f', nargs=1, dest='xvgfilename', help=argparse.SUPPRESS, required=True)
parser.add_argument('-o', nargs=1, dest='output_file', default=["no"], help=argparse.SUPPRESS)
parser.add_argument('--col', nargs=1, dest='col_start', default=[1], type=int, help=argparse.SUPPRESS)
parser.add_argument('--comments', nargs=1, dest='comments', default=['@,#'], help=argparse.SUPPRESS)

#other options
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

args = parser.parse_args()
args.xvgfilename = args.xvgfilename[0]
args.output_file = args.output_file[0]
args.col_start = args.col_start[0]
args.comments = args.comments[0].split(',')

#=========================================================================================
# import modules (doing it now otherwise might crash before we can display the help menu!)
#=========================================================================================

#generic science modules
try:
	import numpy as np
except:
	print "Error: you need to install the numpy module."
	sys.exit(1)
try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.colors as mcolors
	mcolorconv = mcolors.ColorConverter()
	import matplotlib.cm as cm				#colours library
	import matplotlib.ticker
	from matplotlib.ticker import MaxNLocator
	from matplotlib.font_manager import FontProperties
	fontP=FontProperties()
except:
	print "Error: you need to install the matplotlib module."
	sys.exit(1)
try:
	import pylab as plt
except:
	print "Error: you need to install the pylab module."
	sys.exit(1)

#=======================================================================
# sanity check
#=======================================================================
	
if not os.path.isfile(args.xvgfilename):
	print "Error: file " + str(args.xvgfilename) + " not found."
	sys.exit(1)
if args.output_file == "no":
	args.output_file = args.xvgfilename[:-4].split('/')[-1] + '_bi.xvg'
elif args.output_file[-3:] != 'xvg':
	args.output_file += '.xvg'

##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

def load_xvg():
	
	global x_axis
	global s_species
	global nb_species
	global data_energy
			
	#get file content
	filename = args.xvgfilename
	with open(filename) as f:
		lines = f.readlines()
	
	#determine legends and nb of lines to skip
	s_species = []
	tmp_nb_rows_to_skip = 0
	for l_index in range(0,len(lines)):
		line = lines[l_index]
		if line[0] in args.comments:
			tmp_nb_rows_to_skip += 1
			if len(line.split("legend \"")) > 1:
				s_species.append(line.split("legend \"")[-1][:-2])
	
	#get data
	data_density = np.loadtxt(filename, skiprows = tmp_nb_rows_to_skip)
	x_axis = np.copy(data_density[:,0])
	data_density = data_density[:,args.col_start:]
	s_species = s_species[args.col_start-1:]
	nb_species = np.shape(data_density)[1]
	data_density /= np.sum(data_density, axis = 0)
	data_energy = -np.log(data_density)
	data_energy -= np.min(data_energy, axis = 0)
	data_energy[np.isinf(data_energy)] = np.nan
	
	return
def write_xvg():
	#open files
	filename_xvg = os.getcwd() + '/' + str(args.output_file)
	output_xvg = open(filename_xvg, 'w')
	
	#general header
	output_xvg.write("# [free energy estimation - written by botlzmann_inversion v" + str(version_nb) + "]\n")
	output_xvg.write("# input file: " + str(args.xvgfilename) + "\n")	
	#xvg metadata
	output_xvg.write("@ title \"Relative free energies\"\n")
	output_xvg.write("@ xaxis label \"z distance to bilayer center (Angstrom)\"\n")
	output_xvg.write("@ yaxis label \"relative free energies (kT)\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length " + str(nb_species) + "\n")
	for s_index in range(0, nb_species):
		#debug
		print str(s_species[s_index])
		output_xvg.write("@ s" + str(s_index) + " legend \"" + str(s_species[s_index]) + "\"\n")
	#data
	for x_index in range(0, np.shape(data_energy)[0]):
		results = str(x_axis[x_index])
		for s_index in range(0, nb_species):
			results += "	" + str(data_energy[x_index, s_index])
		output_xvg.write(results + "\n")		
	output_xvg.close()	

	return

##########################################################################################
# MAIN
##########################################################################################

load_xvg()
write_xvg()

#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check result in file '" + str(args.output_file)
print ""
sys.exit(0)
