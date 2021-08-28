import numpy as np
import pandas as pd
from random import uniform

#No need in these lines below
#parameter_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0, 1.5, 2.0]
#vibration_param = 'amp'

#filenames = [None]*len(parameter_values)
#for i in range(len(parameter_values)):
#	filenames[i] = '{}{}.particles'.format(vibration_param, parameter_values[i])

#Enter the filenames in this array
filenames = ['test.particles']

#initializing the output data table
output = pd.DataFrame(columns=['Total_Density', 'Region', 'Lower_boundary', 'Upper_boundary', 'Packing_Density',
								'Mean Diameter', 'RMS Diameter'])

row = pd.DataFrame([[0, 0, 0, 0, 0, 0, 0]], columns=['Total_Density', 'Region', 'Lower_boundary', 'Upper_boundary',
								'Packing_Density','Mean Diameter', 'RMS Diameter'])

#assigning values according to the geometry
cylinder_h = 0		#m
cylinder_r = 0.25	#m
t = 0				#m
origin_x = 0.25		#m
origin_y = 0.25		#m

#function to calculate total density
def find_total_density(data):
	#finding the cylinder height
	global cylinder_h
	data = data[data['z'] < 2]
	max_heights = pd.DataFrame(columns=['local_zmax'])

	n_grid = 10

	for i in np.around(np.arange(-cylinder_r, cylinder_r, 2*cylinder_r/n_grid), 3):
		for j in np.around(np.arange(-cylinder_r, cylinder_r, 2*cylinder_r/n_grid), 3):
			local_data = data[(data['x'] > i) & (data['y'] > j) & (data['x'] < i + 2*cylinder_r/n_grid) & (data['y'] < j + 2*cylinder_r/n_grid)]
			if len(local_data.index) >= 10:
				zmax = local_data['z'].max() - t
				max_heights = max_heights.append({'local_zmax': zmax}, ignore_index=True)

	cylinder_h = max_heights['local_zmax'].mean()
	print("The max z-position is", cylinder_h)

	#Monte-Carlo integration to find 1-void fraction
	ntry = 200000
	hit = 0

	for i in range (0, ntry):
		#cylindrical coordinates
		rtry = uniform(0, 0.8*cylinder_r)
		atry = uniform(0, 2*np.pi)#angle
		ztry = uniform(0.2*cylinder_h, 0.8*cylinder_h)#height

    	#converting to cartesian coordinates
		xtry = rtry*np.cos(atry) + origin_x
		ytry = rtry*np.sin(atry) + origin_y

		distance = (xtry-data['x'])**2+(ytry-data['y'])**2+(ztry-data['z'])**2
		if (distance < data['radius']**2).any():
			hit += 1

	total_rho = hit/ntry
	print('The total packing density is', total_rho)
	row.loc[0, 'Total_Density'] = total_rho

#function to find density in a particular region
def find_local_density(data, lower_fraction, upper_fraction):
	region = data[(data['z'] < t + upper_fraction*cylinder_h) & (data['z'] > t + lower_fraction*cylinder_h)]

	ntry = 100000
	hit = 0

	for i in range (0, ntry):
		#cylindrical coordinates
		rtry = uniform(0, cylinder_r*0.8)
		atry = uniform(0, 2*np.pi)#angle
		ztry = uniform(lower_fraction*cylinder_h, upper_fraction*cylinder_h)#height

    	#converting to cartesian coordinates
		xtry = rtry*np.cos(atry) + origin_x
		ytry = rtry*np.sin(atry) + origin_y

		distance = (xtry-region['x'])**2+(ytry-region['y'])**2+(ztry-region['z'])**2
		if (distance < region['radius']**2).any():
			hit += 1

	rho = hit/ntry
	print("The density for this region is", rho)

	row.loc[0, 'Packing_Density'] = rho

#function to find average diameter in a particular region
def find_average_diameter(data, lower_fraction, upper_fraction):
	region = data[(data['z'] < t + upper_fraction*cylinder_h) & (data['z'] > t + lower_fraction*cylinder_h)
					& (np.sqrt(data['x']**2 + data['y']**2) < 0.2)]

	average_d = (region['radius']*2000000/50/15).mean()

	print("The average diameter for this region is", average_d, "micron")

	row.loc[0, 'Mean Diameter'] = average_d

#function to find rms diameter in a prticular region
def find_rms_diameter(data, lower_fraction, upper_fraction):
	region = data[(data['z'] < t + upper_fraction*cylinder_h) & (data['z'] > t + lower_fraction*cylinder_h)
					& (np.sqrt(data['x']**2 + data['y']**2) < 0.2)]

	rms_d = np.sqrt(((region['radius']*2000000/50/15)**2).mean())

	print("The rms diameter for this region is", rms_d, "micron")

	row.loc[0, 'RMS Diameter'] = rms_d

#function to call all previous functions
def full_analyze(data, lower_fraction, upper_fraction):
	print("\n")
	print("Analyzing region between", lower_fraction, "and", upper_fraction, "of total height: \n")

	find_local_density(data, lower_fraction, upper_fraction)
	find_average_diameter(data, lower_fraction, upper_fraction)
	find_rms_diameter(data, lower_fraction, upper_fraction)

	row.loc[0, 'Lower_boundary'] = lower_fraction
	row.loc[0, 'Upper_boundary'] = upper_fraction

	global output
	output = output.append(row, ignore_index=True)

#function to truncate values
def trunc(values, decimals = 0):
    return np.trunc(values*10**decimals)/(10**decimals)

#calling functions
#fractions = [0.001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
fractions = [0.001, 0.2, 0.4, 0.6, 0.8]

for i in range(len(filenames)):
	print("Analyzing", filenames[i])

	data = pd.read_table(filenames[i], delim_whitespace=True, names=['id', 'type', 'x', 'y',
 					'z', 'vx', 'vy', 'vz', 'fx', 'fy', 'fz', 'radius', 'mass'], header=None, skiprows=list(range(9)))
	data = data[['id', 'x', 'y', 'z', 'radius', 'mass']]

	#row.loc[0, 'Parameter'] = trunc(parameter_values[i], 2)

	region_id = len(fractions)-1

	find_total_density(data)

	for i in range(len(fractions)-1):
		row.loc[0, 'Region'] = region_id
		region_id = region_id - 1
		full_analyze(data, fractions[i], fractions[i+1])

#printing and saving the final output data to a file
print("\n", output)
output.to_csv(path_or_buf='test.csv', sep=' ')
