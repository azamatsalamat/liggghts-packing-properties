# LIGGGHTS: packing properties calculator

Calculates the packing density, mean, RMS particle diameters, coordination number and RDF distributions for a packing simulated in the LIGGGHTS package (in cylindrical containers)

## Files:

### density.py
The main script file. Reads data from input file (test.particles). Finds the packing height by dividing the cross-section into rectangular cross sections and averaging the maximum z-coordinates for all segments. Contains functions: find_total_density, find_local_density, find_average_diameter, find_rms_diameter. Density is calculated using Monte-Carlo integration method to take into account particles overlapping. For calculating the diameters and local density, the cylindrical packing is divided into regions across its height (z-axis). 

### coordination.py
Reads data from input file (test.particles). Main function: find_coordination_number. The coordination number (CN) distribution is calculated by counting neighbors within 0.01 of the particle diameter. Then the frequency of each unique coordination number is found.

### rdf.py
Reads data from input file (test.particles). Main functionL find_rdf. The radial distribution function (RDF) is a function of probability of finding a particle within a certain region. 

### test.particles
Example input file. Dump file from LIGGGHTS package. Contains data for all particles on their location, velocities, forces, radius and mass. 

### test.csv
Example output file. Stores the calculated values for total density, local densities and diameters for each region. 

### test_CN.csv
Example output file. Stores the values for coordination numbers and their frequency. 

### test_RDF.csv
Example output file. Stores the values for radial distance and RDF. 

## How to use:

Prerequisites:
- Python
- Python libraries (numpy, pandas)

Put the input files into the same folder with the script file. Enter their names into the 'filenames' array (line 9). Enter the region boudnaries (fractions of the packing height) of interest into the 'fractions' array (line 123). Run the script and obtain the output ('test.csv' file, line 155)
