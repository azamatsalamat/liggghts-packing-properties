# LIGGGHTS packing density, mean and RMS diamter calculator

import numpy as np
import pandas as pd
from random import uniform

# Enter the input filenames in this array separating with comma
fname = 'test'
filenames = [fname+'.particles']

# initializing the output table
output = pd.DataFrame(columns=['Total_Density', 'Region', 'Lower_boundary', 'Upper_boundary', 'Packing_Density',
								'Mean_Diameter', 'RMS_Diameter'])

row = pd.DataFrame([[0, 0, 0, 0, 0, 0, 0]], columns=['Total_Density', 'Region', 'Lower_boundary', 'Upper_boundary',
								'Packing_Density','Mean_Diameter', 'RMS_Diameter'])

# assigning geometry (cylinder)
cylinder_h = 0		#m (cylinder height, to be calculated later based on the packing height)
cylinder_r = 0.25	#m (cylinder radius)
t = 0				#m (thickness of the cylinder)
origin_x = 0.25		#m (x-coordinate of the cylinder origin)
origin_y = 0.25		#m (y-coordinate of the cylinder origin)

def find_rdf(data):
	# cylinder height is found based on the average height of the packing
	global cylinder_h
	data = data[data['z'] < 2]
	max_heights = pd.DataFrame(columns=['local_zmax'])

	# divide the cylinder cross-section into 10 sectors in both directions creating multiple square sectors
	n_grid = 10

	# find the maximum z-coordinate of a particle in each square sector and store them to max_heights
	for i in np.around(np.arange(-cylinder_r, cylinder_r, 2*cylinder_r/n_grid), 3):
		for j in np.around(np.arange(-cylinder_r, cylinder_r, 2*cylinder_r/n_grid), 3):
			local_data = data[(data['x'] > i) & (data['y'] > j) & (data['x'] < i + 2*cylinder_r/n_grid) & (data['y'] < j + 2*cylinder_r/n_grid)]
			if len(local_data.index) >= 10:
				zmax = local_data['z'].max() - t
				max_heights = max_heights.append({'local_zmax': zmax}, ignore_index=True)

	# find the average maximum height among all square sectors and use it as cylinder height
	cylinder_h = max_heights['local_zmax'].mean()
	print("The max z-position is", cylinder_h)

	print("Calculating RDF")

	# create output table for RDF data
	RDF_output = pd.DataFrame(columns=['r', 'RDF'])

	r = 0
	# average number density of particles
	rho_0 = len(data)/(np.pi*cylinder_r**2*cylinder_h)

	# 10 regions to be analyzed
	r_delta = 0.1*cylinder_r

	while r <= cylinder_r:
		x_coord = 0
		y_coord = 0

		# number of particles in a particular region
		N = len(data[(np.sqrt((data['x']-x_coord)**2+(data['y']-y_coord)**2) < r_delta + r) & (np.sqrt((data['x']-x_coord)**2+(data['y']-y_coord)**2) > r)])

		RDF_row = pd.DataFrame([[0, 0]], columns=['r', 'RDF'])
		RDF_row.loc[0, 'r'] = r
		RDF_row.loc[0, 'RDF'] = N/(4*np.pi*r**2*r_delta*rho_0) # formula to calculate RDF
		RDF_output = RDF_output.append(RDF_row)
		print("Finished for r=", r)
		r = r + r_delta
	
	return RDF_output

# function to truncate values
def trunc(values, decimals = 0):
    return np.trunc(values*10**decimals)/(10**decimals)

# calling functions
for i in range(len(filenames)):
	print("Analyzing", filenames[i])

	# read table from the input file
	data = pd.read_table(filenames[i], delim_whitespace=True, names=['id', 'type', 'x', 'y',
 					'z', 'vx', 'vy', 'vz', 'fx', 'fy', 'fz', 'radius', 'mass'], header=None, skiprows=list(range(9)))
	data = data[['id', 'x', 'y', 'z', 'radius', 'mass']]

	output = find_rdf(data)

outfilename = fname + '_RDF.csv'
output.to_csv(path_or_buf=outfilename, sep=' ')