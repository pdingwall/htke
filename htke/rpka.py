import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.optimize import minimize

class RPKA():
	
	"""Class to perform RPKA"""

	def __init__(self, spka_data, points_per_reaction, reactions_per_system = 3):	
	
		"""
		Create RPKA profile using the SPKA data
		
		Parameters
		----------
		spka_data: output dataframe from SPKA function
		points_per_reaction: points per reaction
		reactions_per_system: Number of different reactions per system, default is 3
		"""
		
		self.spka_data = spka_data

		
		# Create new dataframe from rpka_data
		rpka_data = self.spka_data.drop(columns=['Interval Size',
                                   'tR (min)',
                                   'SPKA',
                                   'Relative Time',
                                   'Peak Property',
                                   'Method',
                                   'SPKA Conversion',
                                   'SPKA Ideal t0 Concentration',
                                   'Normalised IR Concentration',
                                   'Initial Conv'])
		
		# If no experimental smoothing then will be no 'Raw Peak Property' column, this will skip error
		try:
			rpka_data = rpka_dat.drop(columns=['Raw Peak Property'])
		except:
			# Do nothing and move on
			pass
		
		# Determine parameters needed for RPKA
		rpka_data['[Excess]'] = rpka_data['[B]0'] - rpka_data['[A]0']
		rpka_data['[B]'] = rpka_data['[A]'] + rpka_data['[Excess]']
		rpka_data['[C]'] = rpka_data['[C]0']
		
		# Add the reagent names, if possible
		try:
			# Read reagents.xlsx to get reagent names
			reagent_names = pd.read_excel('reagents.xlsx', header = None, names = ['Abbreviation','Reagent Name'])

			# Need this as a dictionary
			reagent_names = dict(zip(reagent_names['Abbreviation'], reagent_names['Reagent Name']))
			
			# Map the keys and values
			rpka_data["A"] = rpka_data["A"].map(reagent_names)
			rpka_data["B"] = rpka_data["B"].map(reagent_names)
			rpka_data["C"] = rpka_data["C"].map(reagent_names)    
			
		except:
			print("'reagents.xlsx' could not be found, reagents will remain unnamed")
		
		self.rpka_data = rpka_data
				
		# Currently, we run three reactions per system
		self.reactions_per_system = reactions_per_system

		# t0 data is no longer included
		self.points_per_system = points_per_reaction -1

		# Determine the number of systems in this run
		self.number_of_systems = len(rpka_data) / (self.points_per_system * self.reactions_per_system)
		
	def choose_system(self, system_number):
	
		"""
		For internal use.
		
		Returns a dataframe of a single system. Starts at 0!
		Works on the assumption that there are three reactions per system studied.
		
		Parameters
		----------
		system_number: The system to return, in numerical order, from rpka_data
		
		Returns
		-------
		system: Dataframe of a single system
		"""
		
		# Create a list of the rows in which new systems start, to be use in iloc based slicing
		starting_point = [var for var in range(0,len(self.rpka_data),(self.points_per_system) * self.reactions_per_system)]
		
		# Which system are we finding? Determined by function parameter.
		var = starting_point[system_number]
		
		# Create dataframe of individual system, all relative to the starting point
		system = self.rpka_data.iloc[var : var + self.points_per_system * self.reactions_per_system, :].reset_index(drop = True)
			
		return system

	
	def sum_residuals_b(self, x):
		
		"""
		Returns the sum of the residuals to find the order in B (b), Rate/[B]b must be on y-axis. To be used by minimiser.
		"""

		# Residuals for finding b (order in B), y-axis is Rate/[B]b
		residuals = ((self.standard_rate.to_numpy() / (self.standard_conc_B.to_numpy() ** x)) 
					 - (self.diff_B_rate.to_numpy() / (self.diff_B_conc_B.to_numpy() ** x)))

		# Sum the residuals ** 2
		sum_residuals = np.sum(residuals ** 2)
		
		return sum_residuals


	def straight_line_a(self, x):
		
		"""
		Returns the summed r2 of two kinetic profiles to find the order in A (a), [A]a must be on x-axis. To be used by minimiser.
		"""
		
		# This might need to change
		b = self.order_in_B
		
		# Find the r2 for the standard reaction
		rxn_standard_r2 = np.corrcoef((self.standard_conc_A ** x).astype(float),
									  (self.standard_rate / self.standard_conc_B ** b).astype(float), 1)[0,1]

		# Find the r2 for the change in B reaction
		rxn_diff_B_r2 = np.corrcoef((self.diff_B_conc_A ** x).astype(float),
									  (self.diff_B_rate / self.diff_B_conc_B ** b).astype(float), 1)[0,1]
		
		r2_to_min = -(rxn_standard_r2 + rxn_diff_B_r2)
		
		return r2_to_min

	
	def sum_residuals_c(self, x):
		
		"""
		Returns the sum of the residauls to find the order in C (c), Rate[C]c must be on y-axis. To be used by the minimiser.
		"""

		# Residuals for finding c (order in C), y-axis is Rate/[C]c
		residuals = ((self.standard_rate.to_numpy() / (self.standard_conc_C.to_numpy() ** x)) 
					 - (self.diff_C_rate.to_numpy() / (self.diff_C_conc_C.to_numpy() ** x)))

		# Sum the residuals ** 2
		sum_residuals = np.sum(residuals ** 2)
		
		return sum_residuals
	
		
	def diff_excess(self, min_meth = 'trust-constr', lower_bound = 0, upper_bound = 2, initial_guess = 2):
    
		"""
		Perform a different excess analysis. Reactions MUST be performed in this order:
		Standard
		Different [B]
		Different [C] - where C is catalyst

		Parameters
		----------
		method: Minimiser method: trust-constr, TNC (use boundaries), BFGS, CG, (do not use boundaries)
		lower_bound: lower bound for minimiser
		upper_bound: upper bound for minimiser
		initial_guess: initial guess for minimiser
		
		Returns
		-------
		rpka_data: original dataframe with columnns for each reagent order
		"""
		
		# Set the bounds for the minimiser
		bnds = ((lower_bound, upper_bound),)		
		
		df_a = []
		df_b = []
		df_c = []
		
		# Need a for loop to work through each system
		for var in range(0, self.number_of_systems.astype(int), 1):
			
			# Choose the system
			tmp = self.choose_system(var)

			# Create some variables to pass up to funtions for minimiser
			# Standard Reaction is run first
			self.standard_rate = tmp['Rate'].iloc[0 : self.points_per_system]
			self.standard_conc_A = tmp['[A]'].iloc[0 : self.points_per_system]
			self.standard_conc_B = tmp['[B]'].iloc[0 : self.points_per_system]
			self.standard_conc_C = tmp['[C]'].iloc[0 : self.points_per_system]

			# Reaction with different [B] will always be run second
			self.diff_B_rate = tmp['Rate'].iloc[self.points_per_system : 2 * self.points_per_system]
			self.diff_B_conc_A = tmp['[A]'].iloc[self.points_per_system : 2 * self.points_per_system]
			self.diff_B_conc_B = tmp['[B]'].iloc[self.points_per_system : 2 * self.points_per_system]

			# Reaction with different [C] will always be run third
			self.diff_C_rate = tmp['Rate'].iloc[2 * self.points_per_system : 3 * self.points_per_system]
			self.diff_C_conc_A = tmp['[A]'].iloc[self.points_per_system : 2 * self.points_per_system]
			self.diff_C_conc_C = tmp['[C]'].iloc[2 * self.points_per_system : 3 * self.points_per_system]

			# Determine order in B. Must be first as order_in_B is required to find order in A
			self.order_in_B = minimize(self.sum_residuals_b, initial_guess, method = min_meth, bounds=bnds).x

			# Determine order in A
			order_in_A = minimize(self.straight_line_a, initial_guess, method = min_meth, bounds=bnds).x

			# Determing order in C - lower bound set to 0.5 as catalyst unlikely to be less than this
			order_in_C = minimize(self.sum_residuals_c, initial_guess, method = min_meth, bounds=((0.5, upper_bound),)).x

			# Create a series of each order of length tmp to append in loop
			a = pd.Series(np.round(order_in_A, 2)).repeat(len(tmp))
			b = pd.Series(np.round(self.order_in_B, 2)).repeat(len(tmp))
			c = pd.Series(np.round(order_in_C, 2)).repeat(len(tmp))

			df_a.append(a)
			df_b.append(b)
			df_c.append(c)
		
		# Is this needed? Think it keeps it internal to this function, likely a good idea
		rpka_data = self.rpka_data
		
		# Orders
		rpka_data['Order in A'] = pd.concat(df_a, ignore_index = True)
		rpka_data['Order in B'] = pd.concat(df_b, ignore_index = True)
		rpka_data['Order in C'] = pd.concat(df_c, ignore_index = True)
		
		# Data for plotting in python or excel
		rpka_data['[A]a'] = rpka_data['[A]'].astype(float) ** rpka_data['Order in A']
		rpka_data['[B]b'] = rpka_data['[B]'].astype(float) ** rpka_data['Order in B']
		rpka_data['Rate/[B]b'] = rpka_data['Rate'].astype(float) / rpka_data['[B]b']
		rpka_data['Rate/[C]c'] = rpka_data['Rate'].astype(float) / rpka_data['[C]'] ** rpka_data['Order in C'].astype(float)
		
		return rpka_data
		
		
	def check_results(self, final_rpka_data):

		"""
		Shows the RPKA profiles for the auto-found orders as a visual check.

		Parameters
		----------
		final_rpka_data: Needs the output dataframe from rpka.diff_excess()
		
		Returns
		-------
		Plots of Rate/[B]b vs [A]a and Rate/[C]c vs [A] for each system

		"""  
		
		# Create subplots
		fig, ax = plt.subplots(self.number_of_systems.astype(int), 2, figsize = (14,14))
		fig.tight_layout(w_pad = 5, h_pad = 5) # Makes spacing better

		# Set the row, will be up to max number of systems
		for var_row in range(0, self.number_of_systems.astype(int), 1):

			# Cannot use choose_system as needs rpka_data dataframe as output from diff_excess()
			# Create a list of the rows in which new systems start, to be use in iloc based slicing
			starting_point = [var for var in range(0, len(final_rpka_data), (self.points_per_system) * self.reactions_per_system)]

			# Which system are we finding? Determined by function parameter.
			var = starting_point[var_row]

			# Create dataframe of individual system, all relative to the starting point
			tmp = final_rpka_data.iloc[var : var + self.points_per_system * self.reactions_per_system, :].reset_index(drop = True)

			# Define what to plot
			# Standard reaction - the first of the three systems
			rxn_standard_Aa = tmp['[A]a'].iloc[0 : self.points_per_system]
			rxn_standard_RateBb = tmp['Rate/[B]b'].iloc[0 : self.points_per_system]

			rxn_standard_A = tmp['[A]'].iloc[0 : self.points_per_system]
			rxn_standard_RateCc = tmp['Rate/[C]c'].iloc[0 : self.points_per_system]

			# Change in [B] - the second of the three systems
			rxn_diff_B_Aa = tmp['[A]a'].iloc[self.points_per_system : 2 * self.points_per_system]
			rxn_diff_B_RateBb = tmp['Rate/[B]b'].iloc[self.points_per_system : 2 * self.points_per_system]

			# Change in [C] - the third of the three systems
			rxn_diff_C_A = tmp['[A]'].iloc[2 * self.points_per_system : 3 * self.points_per_system]
			rxn_diff_C_RateCc = tmp['Rate/[C]c'].iloc[2 * self.points_per_system : 3 * self.points_per_system]

			# Find the reaction names
			rxn_standard_name = str('Standard: ' + tmp['Experiment'].iloc[0])
			rxn_diff_B_name = str('Diff B: ' + tmp['Experiment'].iloc[self.points_per_system])
			rxn_diff_C_name = str('Diff C: ' + tmp['Experiment'].iloc[2 * self.points_per_system])

			# Plot order in A and B - first column
			ax[var_row, 0].scatter(rxn_standard_Aa, rxn_standard_RateBb, label = rxn_standard_name)
			ax[var_row, 0].scatter(rxn_diff_B_Aa, rxn_diff_B_RateBb, label = rxn_diff_B_name)

			# Plot order in C - second column
			ax[var_row, 1].scatter(rxn_standard_A, rxn_standard_RateCc, label = rxn_standard_name)
			ax[var_row, 1].scatter(rxn_diff_C_A, rxn_diff_C_RateCc, label = rxn_diff_C_name)

			# Set titles
			ax[var_row, 0].set_title(str(tmp['A'][0]) + ' ^ ' + str(tmp['Order in A'][0]) + ' (linearity) : ' +
									 str(tmp['B'][0]) + ' ^ ' + str(tmp['Order in B'][0]) + ' (overlay)')
			ax[var_row, 1].set_title(str(tmp['C'][0]) + ' ^ ' + str(tmp['Order in C'][0]) + ' (overlay)')

			# Set labels
			ax[var_row, 0].set_ylabel('Rate/[B]b')
			ax[var_row, 0].set_xlabel('[A]a')

			ax[var_row, 1].set_ylabel('Rate/[C]c')
			ax[var_row, 1].set_xlabel('[A]')
			
	        # Add legends
			ax[var_row, 0].legend()
			ax[var_row, 1].legend()


	def check_results_one_system(self, final_rpka_data):

		"""
		Shows the RPKA profiles for the auto-found orders as a visual check, for a singel system only.

		Parameters
		----------
		final_rpka_data: Needs the output dataframe from rpka.diff_excess()
		
		Returns
		-------
		Plots of Rate/[B]b vs [A]a and Rate/[C]c vs [A] for each system

		"""  
		
		# Create subplots
		fig, ax = plt.subplots(1, 2, figsize = (14,8))

		# Set the row, will be up to max number of systems
		for var_row in range(0, self.number_of_systems.astype(int), 1):

			# Cannot use choose_system as needs rpka_data dataframe as output from diff_excess()
			# Create a list of the rows in which new systems start, to be use in iloc based slicing
			starting_point = [var for var in range(0,len(final_rpka_data),(self.points_per_system) * self.reactions_per_system)]

			# Which system are we finding? Determined by function parameter.
			var = starting_point[var_row]

			# Create dataframe of individual system, all relative to the starting point
			tmp = final_rpka_data.iloc[var : var + self.points_per_system * self.reactions_per_system, :].reset_index(drop = True)

			# Define what to plot
			# Standard reaction - the first of the three systems
			rxn_standard_Aa = tmp['[A]a'].iloc[0 : self.points_per_system]
			rxn_standard_RateBb = tmp['Rate/[B]b'].iloc[0 : self.points_per_system]

			rxn_standard_A = tmp['[A]'].iloc[0 : self.points_per_system]
			rxn_standard_RateCc = tmp['Rate/[C]c'].iloc[0 : self.points_per_system]

			# Change in [B] - the second of the three systems
			rxn_diff_B_Aa = tmp['[A]a'].iloc[self.points_per_system : 2 * self.points_per_system]
			rxn_diff_B_RateBb = tmp['Rate/[B]b'].iloc[self.points_per_system : 2 * self.points_per_system]

			# Change in [C] - the third of the three systems
			rxn_diff_C_A = tmp['[A]'].iloc[2 * self.points_per_system : 3 * self.points_per_system]
			rxn_diff_C_RateCc = tmp['Rate/[C]c'].iloc[2 * self.points_per_system : 3 * self.points_per_system]

			# Find the reaction names
			rxn_standard_name = str('Standard: ' + tmp['Experiment'].iloc[0])
			rxn_diff_B_name = str('Diff B: ' + tmp['Experiment'].iloc[self.points_per_system])
			rxn_diff_C_name = str('Diff C: ' + tmp['Experiment'].iloc[2 * self.points_per_system])

			# Plot order in A and B - first column
			ax[0].scatter(rxn_standard_Aa, rxn_standard_RateBb, label = rxn_standard_name)
			ax[0].scatter(rxn_diff_B_Aa, rxn_diff_B_RateBb, label = rxn_diff_B_name)

			# Plot order in C - second column
			ax[1].scatter(rxn_standard_A, rxn_standard_RateCc, label = rxn_standard_name)
			ax[1].scatter(rxn_diff_C_A, rxn_diff_C_RateCc, label = rxn_diff_C_name)

			# Set titles
			ax[0].set_title(str(tmp['A'][0]) + ' ^ ' + str(tmp['Order in A'][0]) + ' (linearity) : ' +
									 str(tmp['B'][0]) + ' ^ ' + str(tmp['Order in B'][0]) + ' (overlay)')
			ax[1].set_title(str(tmp['C'][0]) + ' ^ ' + str(tmp['Order in C'][0]) + ' (overlay)')

			# Set labels
			ax[0].set_ylabel('Rate/[B]b')
			ax[0].set_xlabel('[A]a')

			ax[1].set_ylabel('Rate/[C]c')
			ax[1].set_xlabel('[A]')
			
	        # Add legends
			ax[0].legend()
			ax[1].legend()
		
			
	def manual(self, system_number, a, b, c):
		
		"""
		Perform a manual different excess analysis for a single system.

		Parameters
		----------
		system_number: Which system, starting at 0
		a: Order in A - Looking for linearity
		b: Order in B - Looking for overlay
		c: Order in C - Looking for overlay
		
		Returns
		-------
		Plots of Rate/[B]b vs [A]a and Rate/[C]c vs [A]

		"""  
		tmp = self.choose_system(system_number)
		
		# To plot
		rxn_standard_Aa = tmp['[A]'].iloc[0 : self.points_per_system] ** a
		rxn_standard_RateBb = tmp['Rate'].iloc[0 : self.points_per_system] / tmp['[B]'].iloc[0 : self.points_per_system] ** b

		rxn_standard_A = tmp['[A]'].iloc[0 : self.points_per_system]
		rxn_standard_RateCc = tmp['Rate'].iloc[0 : self.points_per_system] / tmp['[C]'].iloc[0 : self.points_per_system] ** c

		# Change in [B] - the second of the three systems
		rxn_diff_B_Aa = tmp['[A]'].iloc[self.points_per_system : 2 * self.points_per_system] ** a
		rxn_diff_B_RateBb = tmp['Rate'].iloc[self.points_per_system : 2 * self.points_per_system] / tmp['[B]'].iloc[self.points_per_system : 2 * self.points_per_system] ** b

		# Change in [C] - the third of the three systems
		rxn_diff_C_A = tmp['[A]'].iloc[2 * self.points_per_system : 3 * self.points_per_system]
		rxn_diff_C_RateCc = tmp['Rate'].iloc[2 * self.points_per_system : 3 * self.points_per_system] / tmp['[C]'].iloc[2 * self.points_per_system : 3 * self.points_per_system] ** c

		# Find the reaction names
		rxn_standard_name = str('Standard: ' + tmp['Experiment'].iloc[0])
		rxn_diff_B_name = str('Diff B: ' + tmp['Experiment'].iloc[self.points_per_system])
		rxn_diff_C_name = str('Diff C: ' + tmp['Experiment'].iloc[2 * self.points_per_system])

		# Create subplots
		fig, ax = plt.subplots(1, 2, figsize=(15,5))

		# Plot order in A and B - first column
		ax[0].scatter(rxn_standard_Aa, rxn_standard_RateBb, label = rxn_standard_name)
		ax[0].scatter(rxn_diff_B_Aa, rxn_diff_B_RateBb, label = rxn_diff_B_name)

		# Plot order in C - second column
		ax[1].scatter(rxn_standard_A, rxn_standard_RateCc, label = rxn_standard_name)
		ax[1].scatter(rxn_diff_C_A, rxn_diff_C_RateCc, label = rxn_diff_C_name)

		# Set titles
		ax[0].set_title(str(tmp['A'][0]) + ' ^ ' + str(a) + ' (linearity) : ' + str(tmp['B'][0]) + ' ^ ' + str(b) + ' (overlay)')
		ax[1].set_title(str(tmp['C'][0]) + ' ^ ' + str(c) + ' (overlay)')

		# Set labels
		ax[0].set_ylabel('Rate/[B]b')
		ax[0].set_xlabel('[A]a')

		ax[1].set_ylabel('Rate/[C]c')
		ax[1].set_xlabel('[A]')
		
		# Add legends
		ax[0].legend()
		ax[1].legend()