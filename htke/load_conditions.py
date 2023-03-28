import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class Conditions():

	"""
	Load excel file with user defined experimental conditions and smooth data if required.
	
	Excel file must include:
	
	File must be titled "Conditions.xlsx"
	Experiment: Experiment title, eg "Exp1 - R1"
	A: Name of component A
	B: Name of component B
	C: Name of component C - This is the catalyst
	[A]0: Initial concentration of A
	[B]0: Iniital concentration of B
	{C]0: Initial concentration of C
	Initial Conv: SPKA conversion of initial point (ie 10, 20, 30)
	SPKA: Number of SPKA reactions per profile (including t0)
	Interval Size: SPKA interval size (ie 10%, 20% steps)
	tR (min): Residence time (minutes)
	"""


	def read(processed_ir_data):
	
		"""
		Reads and processes "Conditions.xlsx"
		
		Parameter
		---------
		processed_ir_data: Dataframe output by Peaks()
		
		Return
		------
		exerimental_data: Dataframe containing exmperimental data sepcificed in "Conditions.xlsx" and processed_ir_data
		Will also check that numbr of IR datapoints matches the number of specified conditions.
		"""
	
		# Read excel sheet
		conditions = pd.read_excel("Conditions.xlsx")  

		# Expand SPKA conditions for each reaction
		experiments = conditions['Experiment'].unique()

		data=[]
		for var in experiments:
			   
			# For single reaction
			tmp = conditions.loc[conditions['Experiment'] == var]
			
			# Experiments will always start at 0% conversion
			spka_start = [0]
			
			# Create and append df for each SPKA datapoint
			# The next point will start at the specified 'Initial Conversion' and go up by 'Interval Size' * number of SPKA points - 2 (-1 for t0 and -1  as must be one less than in excel)
			spka_points = [tmp['Initial Conv'].iloc[0] + tmp['Interval Size'].iloc[0] * x for x in range(0, tmp['SPKA'].iloc[0] - 2)]
			
			# Make sure it starts at spka_start
			spka = pd.DataFrame(spka_start + spka_points, columns = ['SPKA'])
			
			# Add this to the dataframe
			tmp = tmp.drop(['SPKA'],axis=1) # drop manually filled column
			tmp = tmp.append(spka, ignore_index=True).ffill()
			
			# Take first row as t0 conditions and duplicate
			t0_cond = pd.DataFrame(tmp.iloc[0,:]).T
			pd.concat([t0_cond,tmp]).reset_index()
			
			data.append(tmp)

		processed_conditions=pd.concat(data).reset_index(drop=True).fillna(value='0')

		# Check that you have the same number of peaks as experimental conditions
		if len(processed_ir_data) != len(processed_conditions):
			print('You have a problem: IR datapoints = ',len(processed_ir_data),', Number of conditions = ',len(processed_conditions))
		else:
			print('Inputs seem good: IR Datapoints = ',len(processed_ir_data),', Number of conditions = ',len(processed_conditions))

		# Merge IR and Conditions dataframes
		experimental_data = pd.concat([processed_conditions,processed_ir_data],axis=1)

		return experimental_data
		
		
	def read_no_init_cond(processed_ir_data):
	
		"""
		Reads and processes "Conditions.xlsx"
		
		NOTE: Does not require an 'Initial Conditions' column in Conditions.xlsx, to be used with old data. The SPKA intervals start from 0 and go up in the stated intervals.
		
		Parameter
		---------
		processed_ir_data: Dataframe output by Peaks()
		
		Return
		------
		exerimental_data: Dataframe containing exmperimental data sepcificed in "Conditions.xlsx" and processed_ir_data
		Will also check that numbr of IR datapoints matches the number of specified conditions.
		"""
	
		# Read excel sheet
		conditions = pd.read_excel("Conditions.xlsx")  

		# Expand SPKA conditions for each reaction
		experiments = conditions['Experiment'].unique()

		data=[]
		for var in experiments:
			   
			# For single reaction
			tmp = conditions.loc[conditions['Experiment'] == var]
			
			# Find interval size
			interval_size = tmp['Interval Size'].iloc[0]
			
			# Create and append df for each SPKA datapoint
			spka_points = pd.DataFrame(range(0,tmp['SPKA'].iloc[0] - 1), columns=['SPKA']) * interval_size # Must be one less than in excel sheet
			tmp = tmp.drop(['SPKA'],axis=1) # drop manually filled column
			tmp = tmp.append(spka_points, ignore_index=True).ffill()
			
			# Take first row as t0 conditions and duplicate
			t0_cond = pd.DataFrame(tmp.iloc[0,:]).T
			pd.concat([t0_cond,tmp]).reset_index()
			
			data.append(tmp)

		processed_conditions=pd.concat(data).reset_index(drop=True).fillna(value='0')

		# Check that you have the same number of peaks as experimental conditions
		if len(processed_ir_data) != len(processed_conditions):
			print('You have a problem: IR datapoints = ',len(processed_ir_data),', Number of conditions = ',len(processed_conditions))
		else:
			print('Inputs seem good: IR Datapoints = ',len(processed_ir_data),', Number of conditions = ',len(processed_conditions))

		# Merge IR and Conditions dataframes
		experimental_data = pd.concat([processed_conditions,processed_ir_data],axis=1)

		return experimental_data


	def linear_correction(experimental_data):
		"""
		Performs a linear correction by plotting Peak Property data vs SPKA Conversion. 
		
		Parameters
		----------
		experimental_data =  dataframe output by conditions.read()
		
		Return
		------
		experimental_data = dataframe with 'Peak Property' reaplaced with smoothed line. Original data now in 'Raw Peak Property'.		
		"""
		
		# Find unique reaction numbers
		reaction_list = experimental_data['Experiment'].unique()
		
		df = []
		for var in reaction_list:
			
			# Find a single experiment
			tmp = experimental_data[experimental_data['Experiment'] == var].reset_index(drop = True)

			# Define x and y - without the first point which is the t0
			x = tmp['SPKA'].astype(float)[1:]
			y = tmp['Peak Property'][1:]

			# Find linear fit
			a,b = np.polyfit(x, y, 1)

			# Create the smoothed y data
			best_fit_line = []
			best_fit_line = a * x + b

			# Find t0 point
			t0 = pd.Series(experimental_data[experimental_data['Experiment'] == var]['Peak Property'].iloc[0])

			# Add the t0 and first point back in
			best_fit_line = pd.concat([t0, best_fit_line])

			# Add to dataframe
			tmp['Raw Peak Property'] = tmp['Peak Property']
			tmp['Peak Property'] = best_fit_line

			# Append to list
			df.append(tmp)
					
		experimental_data = pd.concat(df).reset_index(drop = True)
		
		return experimental_data
			
			
	def t0_correction(experimental_data, no_reactions, points_per_reaction):
		
		"""
		Finds all t0 values and averages them, discards outliers.
		
		Parameter
		---------
		experimental_data: dataframe containing experimental data from read().
		no_reaction: number of reactions
		points_per_reaction: points per reaction, including t0
		
		Returns
		--------
		experimental_data: Dataframe with original t0 data replaced
		"""
		
		# Find unique reaction numbers
		reaction_list = experimental_data['Experiment'].unique()

		# Create list of t0 points
		df = [experimental_data[experimental_data['Experiment'] == var]['Peak Property'].iloc[0] for var in reaction_list]

		# Remove any outliers: based on one std dev away from mean
		df_no_outliers = [x for x in df if (x > np.mean(df) - np.std(df))]
		df_no_outliers = [x for x in df_no_outliers if (x < np.mean(df) + np.std(df))]
		average_t0 = np.mean(df_no_outliers)

		# Find the index for 'Peak Property'
		peak_property_index = [x for x, s in enumerate(experimental_data.columns) if 'Peak Property' in s]

		# Replace original t0s with the average, cannot chain so must use a single iloc 
		for var in range(0, no_reactions):
			experimental_data.iloc[var * points_per_reaction, peak_property_index] = average_t0

		return experimental_data


	def plot_corrections(experimental_data):
		"""
		Plots the raw vs corrected IR data points.
		"""
		
		# Find unique reaction numbers
		reaction_list = experimental_data['Experiment'].unique()
		
		for var in reaction_list:
			
			# Find a single experiment
			tmp = experimental_data[experimental_data['Experiment'] == var]
	
			plt.figure(figsize=(10,5))
			plt.scatter(tmp['SPKA'][1:], tmp['Raw Peak Property'][1:], label = 'Uncorrected Raw Data')
			plt.scatter(tmp['SPKA'][1:], tmp['Peak Property'][1:], label = 'Corrected Data')
			plt.title(label = str(var))
			plt.legend()
			plt.show()