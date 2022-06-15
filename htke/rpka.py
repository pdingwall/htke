import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from scipy.optimize import minimize

class RPKA():

	def __init__(self, spka_data):
		self.spka_data = spka_data
		
		# Find unique reaction numbers
		reaction_list = self.spka_data['Experiment'].unique()
		self.reaction_list = reaction_list
		
		# Create new dataframe from rpka_data
		rpka_data = self.spka_data.drop(columns=['Interval Size',
                                   'tR (min)',
                                   'SPKA',
                                   'Relative Time',
                                   'Peak Property',
                                   'Method',
                                   'SPKA Conversion',
                                   'SPKA Ideal t0 Concentration',
                                   'Normalised IR Concentration'])
		
		# Determine parameters needed for RPKA
		rpka_data['[Excess]'] = rpka_data['[B]0'] - rpka_data['[A]0']
		rpka_data['[B]'] = rpka_data['[A]'] + rpka_data['[Excess]']
		rpka_data['[C]'] = rpka_data['[C]0']
		
		self.rpka_data = rpka_data

	# For diff_excess - first experiment
	#def func_exp1(self, x):
		#return self.rate_exp_a/(self.yconc_exp_a**x) - self.x_exp_a

	# For diff_excess - second experiment
	#def func_exp2(self, x):
		#return self.rate_exp_b/(self.yconc_exp_b**x) - self.x_exp_b

	# Residuals for RPKA
	def residuals(self, x):
		return ((self.rate_exp_a/(self.yconc_exp_a**x)) - (self.rate_exp_b/(self.yconc_exp_b**x)))

	def sum_residuals(self, x):
		return sum(self.residuals(x)**2)


	def diff_excess(self):
		"""
		Perform a different excess analysis.
		
		Parameters
		----------
		spka_data must be passed to RPKA()
		
		Returns
		-------
		rpka_data: Dataframe with component orders.		
		"""
		
		# Determine which experiments look at a change in which reagent
		rpka_experiments = {}

		for var in range(1, len(self.reaction_list)):
			
			# Select the first experiment and the nth experiment
			tmp = self.rpka_data.loc[self.rpka_data['Experiment'].isin([self.reaction_list[0],self.reaction_list[var]])]

			# Drop all columns which have only a single unique value, the second column is now the reagent we want - THIS MIGHT CAUSE PROBLEMS IN THE FUTURE
			nunique = tmp.nunique()
			cols_to_drop = nunique[nunique == 1].index
			tmp = tmp.drop(cols_to_drop, axis=1)
			
			# The second columnd will contain the compound which is changing
			rpka_reagent = tmp.columns[1]
			rpka_reagent = rpka_reagent[rpka_reagent.find("[")+1:rpka_reagent.find("]")]
			
			rpka_experiments[self.reaction_list[var]] = rpka_reagent

		# RPKA different excess
		for var in range(1, len(self.reaction_list)):
			
			# exp_a will always be R1
			exp_a = self.reaction_list[0]
			exp_b = self.reaction_list[var]
			
			# Define substrate to process
			rpka_reagent = rpka_experiments[self.reaction_list[var]]
			
			# Check whether the reagent is the limiting reagnet (should be A, the first entry) or not as axes will be different
			if rpka_reagent == rpka_experiments[self.reaction_list[1]]:

				reagent_on_y_axis = '[' + rpka_experiments[self.reaction_list[1]] + ']'
				reagent_on_x_axis = '[' + rpka_experiments[self.reaction_list[var+1]] + ']'

			else:

				reagent_on_y_axis = '[' + rpka_experiments[self.reaction_list[var]] + ']'
				reagent_on_x_axis = '[' + rpka_experiments[self.reaction_list[1]]  + ']'   

			# Define exp_a
			x_exp_a = self.rpka_data[reagent_on_x_axis].loc[self.rpka_data['Experiment'] == exp_a].to_numpy()
			rate_exp_a = self.rpka_data['Rate'].loc[self.rpka_data['Experiment'] == exp_a].to_numpy()
			yconc_exp_a = self.rpka_data[reagent_on_y_axis].loc[self.rpka_data['Experiment'] == exp_a].to_numpy()
			# These must be passed back up to the residuals functions
			self.x_exp_a = x_exp_a
			self.rate_exp_a = rate_exp_a
			self.yconc_exp_a = yconc_exp_a
			
			# Define exp_b
			x_exp_b = self.rpka_data[reagent_on_x_axis].loc[self.rpka_data['Experiment'] == exp_b].to_numpy()
			rate_exp_b = self.rpka_data['Rate'].loc[self.rpka_data['Experiment'] == exp_b].to_numpy()
			yconc_exp_b = self.rpka_data[reagent_on_y_axis].loc[self.rpka_data['Experiment'] == exp_b].to_numpy()
			# These must be passed back up to the residuals functions
			self.x_exp_b = x_exp_b
			self.rate_exp_b = rate_exp_b
			self.yconc_exp_b = yconc_exp_b

			# Calculate order
			reagent_order = np.round(minimize(self.sum_residuals, 1).x,2)
			self.rpka_data['Order in ' + rpka_reagent] = pd.Series(reagent_order)
			
			# Plot
			fig = plt.figure(figsize=(10,3))
			plt.title('Order in ' + rpka_reagent + ' = ' + self.rpka_data['Order in ' + rpka_reagent].iloc[0].astype(str))
			plt.xlabel(reagent_on_x_axis)
			plt.ylabel('Rate/' + reagent_on_y_axis + '^x')
			plt.scatter(x_exp_a, rate_exp_a/yconc_exp_a**self.rpka_data['Order in ' + rpka_reagent].iloc[0], label = exp_a)
			plt.scatter(x_exp_b, rate_exp_b/yconc_exp_b**self.rpka_data['Order in ' + rpka_reagent].iloc[0], label = exp_b)   
			plt.legend()  

		return self.rpka_data.fillna(method="ffill")
	
	