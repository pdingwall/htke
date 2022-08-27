import pandas as pd
import numpy as np

class Conditions():

	"""Load excel file with user defined experimental conditions. Must include:
	
	File must be titled "Conditions.xlsx"
	Experiment: Experiment title, currently formatted "Exp1 - R1"
	A: Name of component A
	B: Name of component B
	C: Name of component C - This is the catalyst
	[A]0: Initial concentration of A
	[B]0: Iniital concentration of B
	{C]0: Initial concentration of C
	SPKA: Number of SPKA reactions per profile (including t0)
	Interval Size: SPKA interval size (ie 10%, 20% steps)
	tR (min): Residence time (minutes)
		
	"""

	def read(processed_ir_data):
	
		"""Reads and processes "Conditions.xlsx"
		
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
		
		
		