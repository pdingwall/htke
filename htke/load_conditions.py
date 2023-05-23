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
		For any number of points per reaction.
		
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

			# Add in the y-intercept
			x[len(x) + 1] = 100
			y[len(y) + 1] = 0
			
			# Define weights - Create an array equal to the length of tmp (ie points per reaction)
			w = np.full(len(tmp), 1.0)

			# Change the final element (for point (100, 0) to be very small)
			w[-1] = 0.0000001

			# Find linear fit - weighted so that theoretical point (100, 0) stays at (100, 0)
			a,b = np.polyfit(x, y, 1, w = 1/w)

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
		

	def linear_correction_wlr(experimental_data):
		"""
		Performs a weighted linear regression correction by plotting Peak Property data vs SPKA Conversion.
		For five points per reaction only.
		Errors are defined from experiments GL-06-50 and GL-06-57.
		
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

			# Add in the y-intercept
			x[len(x) + 1] = 100
			y[len(y) + 1] = 0
			
			# Define weights - Std Dev of each point determined from GL-06-57, fifth point must be (100, 0) so make std dev artifically tiny
			w = np.array([0.081114, 0.074629, 0.037825, 0.032234, 0.0000001])

			# Find linear fit - weights = 1/stddev
			a,b = np.polyfit(x, y, 1, w = 1/w)

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
		


	def linear_correction_wlr_test(experimental_data):
		"""
		Performs a weighted linear regression correction by plotting Peak Property data vs SPKA Conversion.
		Estiamtes the error at each point.
		
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


			# Weighted least squares
			# generate the augmented feature matrix (bias + feature)
			X = np.c_[np.ones(x.shape[0]),x]

			# solution of linear regression
			w_lr = np.linalg.inv(X.T @ X) @ X.T @ y

			# calculate residuals
			res = y - X @ w_lr

			# estimate the covariance matrix
			C = np.diag(res**2)

			# solution of weighted linear regression
			w_wlr = np.linalg.inv(X.T @ np.linalg.inv(C) @ X) @ (X.T @ np.linalg.inv(C) @ y)


			# Create the smoothed y data
			best_fit_line = []
			best_fit_line = w_wlr[1] * x + w_wlr[0]

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


	def plot_corrections(experimental_data, points_per_reaction, reactions_per_system = 3):
		
		"""
		Plots the raw vs corrected IR data points.
		"""
		
		# Find the number of systems
		number_of_systems = len(experimental_data) / (points_per_reaction * reactions_per_system)
		
		fig, ax = plt.subplots(int(number_of_systems), int(reactions_per_system), figsize = (14,8))
		fig.tight_layout(w_pad = 5, h_pad = 5) # Makes spacing better

		for var_row in range(0, int(number_of_systems), 1):

			# Create a list of the rows in which new systems start, to be use in iloc based slicing
			starting_point = [var for var in range(0, len(experimental_data), (points_per_reaction) * reactions_per_system)]

			# Which system are we finding? Determined by function parameter.
			var = starting_point[var_row]

			# Create dataframe of individual system, all relative to the starting point
			tmp = experimental_data.iloc[var : var + points_per_reaction * reactions_per_system, :].reset_index(drop = True)

			# Plot first column - still includes t0 so slice must start at +1 for each row
			ax[var_row, 0].scatter(tmp['SPKA'][1 : points_per_reaction], 
								   tmp['Raw Peak Property'][1 : points_per_reaction], label = 'Uncorrected Raw Data')
			ax[var_row, 0].scatter(tmp['SPKA'][1 : points_per_reaction], 
								   tmp['Peak Property'][1 : points_per_reaction], label = 'Corrected Data')

			# Plot second column
			ax[var_row, 1].scatter(tmp['SPKA'][points_per_reaction + 1 : 2 * points_per_reaction], 
								   tmp['Raw Peak Property'][points_per_reaction + 1 : 2 * points_per_reaction], label = 'Uncorrected Raw Data')
			ax[var_row, 1].scatter(tmp['SPKA'][points_per_reaction + 1 : 2 * points_per_reaction], 
								   tmp['Peak Property'][points_per_reaction + 1 : 2 * points_per_reaction], label = 'Corrected Data')

			# Plot third column
			ax[var_row, 2].scatter(tmp['SPKA'][2 * points_per_reaction + 1 : 3 * points_per_reaction], 
								   tmp['Raw Peak Property'][2 * points_per_reaction + 1 : 3  * points_per_reaction], label = 'Uncorrected Raw Data')
			ax[var_row, 2].scatter(tmp['SPKA'][points_per_reaction + 1 : 2 * points_per_reaction], 
								   tmp['Peak Property'][2 * points_per_reaction + 1 : 3 * points_per_reaction], label = 'Corrected Data')

			# Set Styles
			for var2 in range(0, reactions_per_system):
				# Set Titles
				ax[var_row, var2].set_title(str(tmp['Experiment'][var2 * points_per_reaction]))
					
				# Set labels
				ax[var_row, var2].set_ylabel('Peak Property')
				ax[var_row, var2].set_xlabel('SPKA Conversion')
				
				# Add legends
				ax[var_row, var2].legend()				

				# Set ylim
				ax[var_row, var2].set_ylim([0, float(tmp['Peak Property'].max()) * 1.05])


	def plot_corrections_old(experimental_data):
		"""
		Plots the raw vs corrected IR data points in a single column.
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
			plt.xlabel('SPKA Conversion')
			plt.ylabel('Peak Property')
			plt.legend()
			plt.show()