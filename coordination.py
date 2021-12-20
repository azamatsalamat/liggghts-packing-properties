# LIGGGHTS packing density, mean and RMS diamter calculator

import numpy as np
import pandas as pd
from random import uniform

# Enter the input filenames in the fname variable
fname = 'test'
filenames = [fname+'.particles']

# assigning geometry (cylinder)
cylinder_h = 0		#m (cylinder height, to be calculated later based on the packing height)
cylinder_r = 0.25	#m (cylinder radius)
t = 0				#m (thickness of the cylinder)
origin_x = 0.25		#m (x-coordinate of the cylinder origin)
origin_y = 0.25		#m (y-coordinate of the cylinder origin)

def find_coordination_number(data):
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

	print("Calculating CN")

	# create a new table for CN data
	CN_output = pd.DataFrame(columns=['CN'])

	# calculate the number of neighbors (coordinatio number) for each particle
	for i in range(len(data)):
		CN_row = pd.DataFrame([[0]], columns=['CN'])

		# neighbors are all particles located within 0.01 of the particle radius (subject to change)
		distance = data.iloc[i, data.columns.get_loc('radius')]*1.005
		x_coord = data.iloc[i, data.columns.get_loc('x')]
		y_coord = data.iloc[i, data.columns.get_loc('y')]
		z_coord = data.iloc[i, data.columns.get_loc('z')]
		neighbors = data[(np.sqrt((data['x']-x_coord)**2+(data['y']-y_coord)**2+(data['z']-z_coord)**2) < distance+data['radius'])]

		CN_row.loc[0, 'CN'] = len(neighbors) - 1
		CN_output = CN_output.append(CN_row, ignore_index=True)
		print("Calculated for a particle", i)

	CN_numbers = CN_output.CN.unique()
	CN_dist = pd.DataFrame(columns=['CN', 'Occurrences'])
	
	# calculate the frequency of each unique coordination number and store in a table
	for i in range(len(CN_numbers)-1):
		CN_distrow = pd.DataFrame([[0, 0]], columns=['CN', 'Occurrences'])
		CN_distrow.loc[0, 'CN'] = CN_numbers[i]
		CN_distrow.loc[0, 'Occurrences'] = len(CN_output[(CN_output['CN'] == CN_numbers[i])])
		CN_dist = CN_dist.append(CN_distrow)

	return CN_dist

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

	output = find_coordination_number(data)

outfilename = fname + '_CN.csv'
output.to_csv(path_or_buf=outfilename, sep=' ')