# LIGGGHTS packing density, mean and RMS diamter calculator

import numpy as np
import pandas as pd
from random import uniform

# Enter the input filenames in this array separating with comma
filenames = ['test.particles']

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

# function to calculate total density
def find_total_density(data):
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

	#Monte-Carlo integration to find packing density (one minus void fraction)
	ntry = 200000 # number of trials (increasing this will improve accuracy but also increase computation time)
	hit = 0

	for i in range (0, ntry):
		# pick a random location in cylindrical coordinates
		rtry = uniform(0, 0.8*cylinder_r) # radius (calculate only up to 80% of the cylinder radius to eliminate possible wall effects)
		atry = uniform(0, 2*np.pi) # angle
		ztry = uniform(0.2*cylinder_h, 0.8*cylinder_h) # height (limit between 20% and 80% of the total height to eliminate possible wall effects)

    	# converting to cartesian coordinates
		xtry = rtry*np.cos(atry) + origin_x
		ytry = rtry*np.sin(atry) + origin_y

		# compare the distance between all particles and the random location to detect if it is located inside a particle
		distance = (xtry-data['x'])**2+(ytry-data['y'])**2+(ztry-data['z'])**2
		if (distance < data['radius']**2).any():
			hit += 1

	total_rho = hit/ntry # total packing density
	print('The total packing density is', total_rho)
	# adding the result to the output table
	row.loc[0, 'Total_Density'] = total_rho

# function to find density in a particular region (between two fractions of cylinder height)
def find_local_density(data, lower_fraction, upper_fraction):
	region = data[(data['z'] < t + upper_fraction*cylinder_h) & (data['z'] > t + lower_fraction*cylinder_h)]

	ntry = 100000 # increase this to improve accuracy; reduce to save computational costs
	hit = 0

	for i in range (0, ntry):
		# random lcoation in cylindrical coordinates
		rtry = uniform(0, cylinder_r*0.8) # radius
		atry = uniform(0, 2*np.pi) # angle
		ztry = uniform(lower_fraction*cylinder_h, upper_fraction*cylinder_h) # height

    	# converting to cartesian coordinates
		xtry = rtry*np.cos(atry) + origin_x
		ytry = rtry*np.sin(atry) + origin_y

		distance = (xtry-region['x'])**2+(ytry-region['y'])**2+(ztry-region['z'])**2
		if (distance < region['radius']**2).any():
			hit += 1

	rho = hit/ntry # local packing density
	print("The density for this region is", rho)

	row.loc[0, 'Packing_Density'] = rho

# function to find average diameter in a particular region (between two fractions of cylinder height)
def find_average_diameter(data, lower_fraction, upper_fraction):
	region = data[(data['z'] < t + upper_fraction*cylinder_h) & (data['z'] > t + lower_fraction*cylinder_h)
					& (np.sqrt(data['x']**2 + data['y']**2) < 0.8*cylinder_r)]

	average_d = (region['radius']*2).mean() # mean diameter in the region

	print("The average diameter for this region is", average_d, "m")

	row.loc[0, 'Mean_Diameter'] = average_d

# function to find rms diameter in a prticular region (between two fractions of cylinder height)
def find_rms_diameter(data, lower_fraction, upper_fraction):
	region = data[(data['z'] < t + upper_fraction*cylinder_h) & (data['z'] > t + lower_fraction*cylinder_h)
					& (np.sqrt(data['x']**2 + data['y']**2) < 0.8*cylinder_r)]

	rms_d = np.sqrt(((region['radius']*2)**2).mean()) # RMS diameter in the region

	print("The RMS diameter for this region is", rms_d, "m")

	row.loc[0, 'RMS_Diameter'] = rms_d

# function to truncate values
def trunc(values, decimals = 0):
    return np.trunc(values*10**decimals)/(10**decimals)

# fractions of the cylinder height (put them in increasing order)
fractions = [0.001, 0.2, 0.4, 0.6, 0.8]

# calling functions
for i in range(len(filenames)):
	print("Analyzing", filenames[i])

	# read table from the input file
	data = pd.read_table(filenames[i], delim_whitespace=True, names=['id', 'type', 'x', 'y',
 					'z', 'vx', 'vy', 'vz', 'fx', 'fy', 'fz', 'radius', 'mass'], header=None, skiprows=list(range(9)))
	data = data[['id', 'x', 'y', 'z', 'radius', 'mass']]

	region_id = len(fractions)-1

	find_total_density(data)

	# analyze for each region
	for i in range(len(fractions)-1):
		row.loc[0, 'Region'] = region_id
		region_id = region_id - 1
		print("\n")
		print("Analyzing region between", fractions[i], "and", fractions[i+1], "of the total height:")

		find_local_density(data, fractions[i], fractions[i+1])
		find_average_diameter(data, fractions[i], fractions[i+1])
		find_rms_diameter(data, fractions[i], fractions[i+1])

		row.loc[0, 'Lower_boundary'] = fractions[i]
		row.loc[0, 'Upper_boundary'] = fractions[i+1]
		output = output.append(row, ignore_index=True) # add the results to the final output table

# printing and saving the final output data to a file
print("\n", output)
output.to_csv(path_or_buf='test.csv', sep=' ')
