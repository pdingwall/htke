import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from scipy.optimize import minimize

class SPKA():

	"""Class to create SPKA profiles."""
	
	def __init__(self, experimental_data):
		self.experimental_data = experimental_data
		
		# Find unique reaction numbers
		reaction_list = experimental_data['Experiment'].unique()
		self.reaction_list = reaction_list
		
	def spka(self, sm_monitored):
	
		"""
		Create SPKA profiles. Must pass only a single peak from Peaks().
		
		Parameter
		---------
		experimental_data: Dataframe output by Conditions().
		sm_monitored: Must be a limiting reagent. Formatted the same as in "Conditions.xlsx" (ie '[A]0')
		
		Return
		------
		spka_data: Dataframe containing calculated SPKA profiles of rate and concentration
		"""
		
		# Stop that damn slice error
		pd.options.mode.chained_assignment = None  # default='warn'

		data = []

		for reaction in self.reaction_list:
				
			# Just single reaction
			tmp = self.experimental_data.loc[self.experimental_data['Experiment'] == reaction]

			# Normalise the IR concentration against the t0 value, set to ideal concentration
			norm_ir_conc = tmp['Peak Property'] / tmp['Peak Property'].iloc[0] * tmp[sm_monitored].iloc[0]

			# Drop the first row as this is t0 data and should not be part of the SPKA profile
			tmp.drop(index=tmp.index[0], axis=0, inplace = True)

			# Calculate the SPKA Conversion
			spka_conv = 1 - (tmp['SPKA'])/100 

			# Calculate the ideal starting concentration for each datapoint
			# Change this if have the actual t0s for each datapoint
			spka_ideal_conc = spka_conv * tmp[sm_monitored]

			# Add these to dataframe
			tmp['SPKA Conversion'] = spka_conv
			tmp['SPKA Ideal t0 Concentration'] = spka_ideal_conc

			# Find the residence time
			residence_time = tmp['tR (min)'].iloc[0]

			# Calculate the SPKA rate data
			spka_rate = (spka_ideal_conc - norm_ir_conc) / residence_time

			# Calculate the SPKA concentration data
			spka_conc = (spka_ideal_conc - norm_ir_conc) / 2 + norm_ir_conc

			# Add these to dataframe
			tmp['Normalised IR Concentration'] = norm_ir_conc
			tmp['Rate'] = spka_rate
			tmp[sm_monitored[:3]] = spka_conc
					   
			data.append(tmp)
			
		spka_data = pd.concat(data).reset_index(drop=True)

		return spka_data


	def plot(self, spka_data, points_per_reaction, reactions_per_system = 3):
	
		"""
		Plot spka_data
		
		parameters
		----------
		spka_data: Dataframe output by spka()
		points_per_reaction: points per reaction as determined earlier in sheet, include t0
		reaction_per_system: default of three
		
		Return
		------
		An indiviudal plot for each reaction.
		"""

		# spka_data has the t0 removed so need to update the points_per_reaciton
		updated_points_per_reaction = points_per_reaction.astype(int) - 1

		# Find the number of systems
		number_of_systems = len(spka_data) / (updated_points_per_reaction * reactions_per_system)

		fig, ax = plt.subplots(int(number_of_systems), int(reactions_per_system), figsize = (14,14))
		fig.tight_layout(w_pad = 5, h_pad = 5) # Makes spacing better

		for var_row in range(0, int(number_of_systems), 1):

			# Create a list of the rows in which new systems start, to be use in iloc based slicing
			starting_point = [var for var in range(0, len(spka_data), (updated_points_per_reaction) * reactions_per_system)]

			# Which system are we finding? Determined by function parameter.
			var = starting_point[var_row]

			# Create dataframe of individual system, all relative to the starting point
			tmp = spka_data.iloc[var : var + updated_points_per_reaction * reactions_per_system, :].reset_index(drop = True)
			
			# Plot first column
			ax[var_row, 0].scatter(tmp['[A]'][0 : updated_points_per_reaction], 
								   tmp['Rate'][0 : updated_points_per_reaction])

			# Plot second column
			ax[var_row, 1].scatter(tmp['[A]'][updated_points_per_reaction : 2 * updated_points_per_reaction], 
								   tmp['Rate'][updated_points_per_reaction : 2 * updated_points_per_reaction])

			# Plot third column
			ax[var_row, 2].scatter(tmp['[A]'][2 * updated_points_per_reaction : 3 * updated_points_per_reaction], 
								   tmp['Rate'][2 * updated_points_per_reaction : 3  * updated_points_per_reaction])

			# Set Styles
			for var2 in range(0, reactions_per_system):
				# Set Titles
				ax[var_row, var2].set_title(str(tmp['Experiment'][var2 * updated_points_per_reaction]))
					
				# Set labels
				ax[var_row, var2].set_ylabel('Rate')
				ax[var_row, var2].set_xlabel('[A]')

				# Set lims
				ax[var_row, var2].set_ylim([0, float(tmp['Rate'].max()) * 1.1])
				ax[var_row, var2].set_xlim([0, 0.1])

	
	def plot_old(self, spka_data):
	
		"""
		Plot spka_data
		
		parameters
		----------
		spka_data: Dataframe output by spka()
		
		Return
		------
		An indiviudal plot for each reaction
		"""
		
		for reaction in self.reaction_list:
			tmp = spka_data.loc[spka_data['Experiment'] == reaction]

			fig = plt.figure(figsize=(10,3))
			ax = fig.subplots()
			ax.title.set_text(reaction)
			ax.scatter(tmp['[A]'],tmp['Rate'], color='r')

			#  Find line of best fit => np.polyfit(x, y, 1)
			a, b = np.polyfit(tmp['[A]'].astype(float), tmp['Rate'].astype(float), 1)
			plt.plot(tmp['[A]'], a * tmp['[A]'] + b)

			
	def compare(self, sm_monitored):
	
		"""
		Compare SPKA profiles constructed using different peak properties.
		Note: This requires Peaks.compare() and that the dataframe 'compare' is used in Condition.read()
				
		Parameters
		----------
		sm_monitored: Must be a limiting reagent. Formatted the same as in "Conditions.xlsx" (ie '[A]0')
				
		Return
		r2: Dataframe of r2 values for each reaciton and peak property
		------
		
		"""
		# Define the peak_properties to compare
		peak_properties = ['Prominence', 'Experimental Area', 'Fitted Area']
		
		# Stop that damn slice error
		pd.options.mode.chained_assignment = None  # default='warn'
		data = []

		# Create SPKA profiles for everypeak property
		for reaction in self.reaction_list:

			# Just single reaction
			tmp = self.experimental_data.loc[self.experimental_data['Experiment'] == reaction]
			
			for var in peak_properties:

				# Normalise the IR concentration against the t0 value, set to ideal concentration
				norm_ir_conc = tmp[var] / tmp[var].iloc[0] * tmp[sm_monitored].iloc[0]

				# Calculate the SPKA Conversion - using [1:] ignores the first row so don't have to drop from dataframe
				spka_conv = 1 - (tmp['SPKA'][1:])/100 

				#Calculate the ideal starting concentration for each datapoint
				#Change this if have the actual t0s for each datapoint
				spka_ideal_conc = spka_conv * tmp[sm_monitored]

				# Add these to dataframe
				tmp['SPKA Conversion'] = spka_conv
				tmp['SPKA Ideal t0 Concentration'] = spka_ideal_conc

				# Find the residence time
				residence_time = tmp['tR (min)'].iloc[0]

				# Calculate the SPKA rate data
				spka_rate = (spka_ideal_conc - norm_ir_conc) / residence_time

				# Calculate the SPKA concentration data
				spka_conc = (spka_ideal_conc - norm_ir_conc) / 2 + norm_ir_conc

				# Add these to dataframe
				#tmp['Normalised IR Concentration'] = norm_ir_conc
				tmp['Rate - ' + var] = spka_rate
				tmp[sm_monitored[:3] + ' - ' + var] = spka_conc

			data.append(tmp)

		# Define spka_data
		spka_data = pd.concat(data).dropna().reset_index(drop=True)
		
		# Plot the data
		# One plot per reaction
		for reaction in self.reaction_list:
			tmp = spka_data.loc[spka_data['Experiment'] == reaction]
			fig = plt.figure(figsize=(10,3))
			ax = fig.subplots()
			
			# Each peak property, with a best fit line, in each plot
			for var in peak_properties:
				ax.title.set_text(reaction)
				ax.scatter(tmp[sm_monitored[:3] + ' - ' + var], tmp['Rate - ' + var], label = var)
				ax.legend()

				#  Find line of best fit => np.polyfit(x, y, 1)
				a, b = np.polyfit(tmp[sm_monitored[:3] + ' - ' + var].astype(float), tmp['Rate - ' + var].astype(float), 1)
				ax.plot(tmp[sm_monitored[:3] + ' - ' + var], a * tmp[sm_monitored[:3] + ' - ' + var] + b)
		
		# Find r2 values
		df3 = []
		for reaction in self.reaction_list:
				
			# Just one reaction
			tmp = spka_data.loc[spka_data['Experiment'] == reaction]
			
			# For loop along each peak property
			df =[]
			for var in peak_properties:
				r2 = np.corrcoef(tmp[sm_monitored[:3] + ' - ' + var].astype(float), tmp['Rate - ' + var].astype(float))[0,1]
				df.append(r2)
				df2 = pd.DataFrame(df)
			df3.append(df2)
				
		# Must transpose df to be the same orientation as previous dfs
		final = pd.concat(df3, axis=1).T.reset_index(drop=True)
		final.columns = list(peak_properties)
		final.loc['Sum'] = final.sum()

		# Append Reaction name to DF
		reaction_list2 = list(self.reaction_list)
		reaction_list2.append('Sum')
		final['Reaction'] = reaction_list2
		
		return final		