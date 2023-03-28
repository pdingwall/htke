import pandas as pd
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_prominences
from scipy.optimize import curve_fit
from scipy import integrate as intg
from scipy.optimize import leastsq
from scipy.optimize import minimize
from scipy.optimize import minimize_scalar

class Peaks():

	"""Process raw IR data:
	
	prominence_all: Find the prominence of peaks of every wavelength in raw ir_data
	prominence: Find the prominence of peaks of a single wavelength in raw ir_data
	height: Find the height of peaks of a single wavelength in raw ir_data
	exp_area: Find the areas of the peaks in the experimental data
	fitted_area_sp: Examine a single peak and fit a Gaussian curve for a single wavelength
	fitted_area: Find the areas of the fitted idealised Gaussian peaks
	plot: Plot a single picked peaks against the raw data
	"""

	
	def __init__(self, ir_data):
		self.ir_data = ir_data
		
		# Create list of all the IR peaks monitored in dataframe
		peak_list = [x for x in list(ir_data) if 'Peak' in x]
		self.peak_list = peak_list
		

	def prominence_finder(self, var):
	# Used by other function to automatically find the prominence threshold
		
		# Read "Conditions.xlsx"
		quick_cond = pd.read_excel("Conditions.xlsx")

		# Calculate the number of reactions this specifies
		no_reactions = len(quick_cond['Experiment']) * quick_cond['SPKA'].iloc[0]		

		# Find the prominence
		x = self.ir_data[self.peak_of_interest]
		peaks = find_peaks(x, prominence = var)
		
		# Number of peaks found for given prominence
		no_peaks = len(peaks[0])
		
		# The absolute difference between the number of peaks found and the number of reactions specified
		abs_diff_prom = abs(no_peaks - no_reactions)
		
		return abs_diff_prom	


	def auto_prom(self, peak_of_interest, lower_bound = 0, upper_bound = 1):
		
		"""Find the prominence threshold. Requires "Conditions.xlsx" to exist.
		
		Parameter
		---------
		peak_of_interest: Single wavelength to integrate (ie 'Peak at 1704 cm-1') 
		lower_bound = 0
		upper_bound = 1
		"""
		
		self.peak_of_interest = peak_of_interest
		
		# Read "Conditions.xlsx"
		quick_cond = pd.read_excel("Conditions.xlsx")

		# Calculate the number of reactions this specifies
		no_reactions = len(quick_cond['Experiment'])
		
		# Find the number of points_per_reaction - should be the same for each reaction
		points_per_reaction = quick_cond['SPKA'][0]	
		
		# Golden > Bounded > Brent ??? Could also look into tolerance
		mins = minimize_scalar(self.prominence_finder, method='Golden', bounds=(lower_bound, upper_bound))
		
		return mins.x, no_reactions, points_per_reaction


	def height_finder(self, var):
	# Used by other function to automatically find the height threshold
		
		# Read "Conditions.xlsx"
		quick_cond = pd.read_excel("Conditions.xlsx")

		# Calculate the number of reactions this specifies
		no_reactions = len(quick_cond['Experiment']) * quick_cond['SPKA'].iloc[0]		
		
		# Find the height
		x = self.ir_data[self.peak_of_interest]
		peaks = find_peaks(x, height = var)
		
		# Number of peaks found for given height
		no_peaks = len(peaks[0])
		
		# The absolute difference between the number of peaks found and the number of reactions specified
		abs_diff_prom = abs(no_peaks - no_reactions)
		
		return abs_diff_prom	


	def auto_height(self, peak_of_interest, lower_bound = 0, upper_bound = 1):
		
		"""Find the height threshold. Requires "Conditions.xlsx" to exist.
		
		Parameter
		---------
		peak_of_interest: Single wavelength to integrate (ie 'Peak at 1704 cm-1') 
		lower_bound = 0
		upper_bound = 1
		"""
		
		self.peak_of_interest = peak_of_interest
		
		mins = minimize_scalar(self.height_finder, method='bounded', bounds=(lower_bound, upper_bound))
		
		return mins.x
	
		
	def prominence_all(self, peak_threshold):
	
		"""Find the prominence of peaks of every wavelength in raw ir_data
		
		Parameter
		---------
		peak_threshold : Peak threshold
		
		Returns
		-------
		processed_ir_data : Dataframe of relative time and prominence for each wavelength		
		"""
		
		# For loop through each defined peak
		processed_ir_data = []
		list_of_dfs = []
				
		for var in self.peak_list:
			
			# Find single peak
			x = self.ir_data[var]

			# Find peak height and position
			peaks = find_peaks(x, prominence = peak_threshold)
			peak_prominence = peaks[1]
			peak_pos = self.ir_data['Relative Time'][peaks[0]]

			# Add to dataframe
			tmp = pd.DataFrame(peak_pos).reset_index(drop=True)
			tmp2 = pd.DataFrame(peak_prominence).rename(columns = {'prominences':var}) # CHANGE THIS IF CHANGING METHOD
			tmp = pd.concat([tmp,tmp2],axis=1).rename(columns={'Relative Time':'Relative Time ' + var})

			list_of_dfs.append(tmp)

		processed_ir_data = pd.concat(list_of_dfs, axis=1)
				
		# Plot all picked peaks
		for var in self.peak_list:
			fig = plt.figure(figsize=(20,3))
			ax = fig.subplots()
			ax.title.set_text(var)
			ax.plot(self.ir_data['Relative Time'], self.ir_data[var])
			ax.scatter(processed_ir_data['Relative Time ' + var], processed_ir_data[var], color='r')		
		
		return processed_ir_data


	def prominence(self, peak_threshold, peak_of_interest):
	
		"""Find the prominence of peaks of a single wavelength in raw ir_data
		
		Parameter
		---------
		peak_threshold : Peak threshold
		peak_of_interest: Single wavelength to integrate (ie 'Peak at 1704 cm-1') 
		
		Returns
		-------
		processed_ir_data : Dataframe of relative time and prominence for a chosen wavelength		
		"""
			
		# Find single peak
		x = self.ir_data[peak_of_interest]

		# Find peak height and position
		peaks = find_peaks(x, prominence = peak_threshold)
		peak_prominence = peaks[1]
		peak_pos = self.ir_data['Relative Time'][peaks[0]]

		# Add to dataframe
		tmp = pd.DataFrame(peak_pos).reset_index(drop=True)
		tmp2 = pd.DataFrame(peak_prominence).rename(columns = {'prominences':'Peak Property'})
		processed_ir_data = pd.concat([tmp,tmp2],axis=1)
		
		# Drop the 'base' columns
		processed_ir_data.drop([col for col in processed_ir_data.columns if 'base' in col], axis=1, inplace=True)
		
		# Add the method
		processed_ir_data['Method'] = 'prominence'

		return processed_ir_data
		
		
	def height(self, peak_threshold, peak_of_interest):
	
		"""Find the height of peaks of a single wavelength in raw ir_data
		
		Parameter
		---------
		peak_threshold : Peak threshold
		peak_of_interest: Single wavelength to integrate (ie 'Peak at 1704 cm-1') 
		
		Returns
		-------
		processed_ir_data : Dataframe of relative time and heights of peaks for a chosen wavelength		
		"""
			
		# Find single peak
		x = self.ir_data[peak_of_interest]

		# Find peak height and position
		peaks = find_peaks(x, height = peak_threshold)
		peak_height = peaks[1]
		peak_pos = self.ir_data['Relative Time'][peaks[0]]

		# Add to dataframe
		tmp = pd.DataFrame(peak_pos).reset_index(drop=True)
		tmp2 = pd.DataFrame(peak_height).rename(columns = {'peak_heights':'Peak Property'})
		processed_ir_data = pd.concat([tmp,tmp2],axis=1)

		# Add the method
		processed_ir_data['Method'] = 'height'
		
		return processed_ir_data
	
	
	def exp_area(self, peak_threshold, residence_time, peak_of_interest, time_adjust_before = 0, time_adjust_after = 0):
		
		"""Find the areas of the peaks in the experimental data
		
		Parameter
		---------
		peak_threshold : Peak threshold
		residence_time: Residence time of a single peak
		peak_of_interest: Single wavelength to integrate (ie 'Peak at 1704 cm-1')
		time_adjust_before: Time to remove from front of peak
		time_adjust_after: Time to add after peak
		
		Returns
		-------
		experimental_area		
		"""
		
		# Requires prominence data (don't use prominence function as don't want returned processed_ir_data)
		# Find single peak
		x = self.ir_data[peak_of_interest]

		# Find peak height and position
		peaks = find_peaks(x, prominence = peak_threshold)
		peak_prominence = peaks[1]
		peak_pos = self.ir_data['Relative Time'][peaks[0]]
		
		# Make it a list
		list_of_peaks = list(peak_pos)
		
		df = []
				
		for var in list_of_peaks:
	
			# Find a single peak
			single_peak = self.ir_data.loc[(self.ir_data['Relative Time'] >= var - (residence_time/2 - time_adjust_before)) 
										 & (self.ir_data['Relative Time'] <= var + (residence_time/2 + time_adjust_after))]
			
			# Find peak area of the experimental data
			exp_peak_area = np.trapz(single_peak[peak_of_interest])
			
			df.append(pd.Series(exp_peak_area))
		
		experimental_area = pd.concat(df)
			
		# Finalise dataframe
		tmp = pd.DataFrame(peak_pos).reset_index(drop=True)
		tmp2 = pd.DataFrame(experimental_area, columns=['Peak Property']).reset_index(drop=True)
		processed_ir_data = pd.concat([tmp,tmp2],axis=1)
		
		# Add the method
		processed_ir_data['Method'] = 'experimental area'
		
		return processed_ir_data


	def gaussian(self, x, a, b, c):
		"""
		For use internally, defines a gaussian function to fit to experimental data.
		"""
		
		return a*np.exp(-np.power(x - b, 2)/(2*np.power(c, 2)))


	def fitted_area(self, peak_threshold, residence_time, peak_of_interest, time_adjust_before = 0, time_adjust_after = 0, p0 = [5, -1.5, 2]):
		
		"""Find the areas of the fitted idealised Gaussian peaks
		
		Parameter
		---------
		peak_threshold : Peak threshold
		residence_time: Residence time of a single peak
		peak_of_interest: Single wavelength to integrate (ie 'Peak at 1704 cm-1')
		time_adjust_before: Time to remove from front of peak
		time_adjust_after: Time to add after peak
		
		Returns
		-------
		experimental_area
		"""
		
		# Requires prominence data (don't use prominence function as don't want returned processed_ir_data)
		# Find single peak
		x = self.ir_data[peak_of_interest]

		# Find peak height and position
		peaks = find_peaks(x, prominence = peak_threshold)
		peak_prominence = peaks[1]
		peak_pos = self.ir_data['Relative Time'][peaks[0]]
		
		# Make it a list
		list_of_peaks = list(peak_pos)
		
		df = []
				
		for var in list_of_peaks:
	
			# Find a single peak
			single_peak = self.ir_data.loc[(self.ir_data['Relative Time'] >= var - (residence_time/2 - time_adjust_before)) 
										 & (self.ir_data['Relative Time'] <= var + (residence_time/2 + time_adjust_after))]
			   
			# Fit Gaussian and find parameters
			X = np.arange(0, len(single_peak))
			pars, cov = curve_fit(f=self.gaussian, xdata=X, ydata=single_peak[peak_of_interest], p0=[5, -1, 2], bounds=(-np.inf, np.inf))

			# Find area (function, starting point, end point, args=arguments to use in function)
			peak_area = intg.quad(self.gaussian, 0, len(single_peak)-1, args=(pars[0],pars[1],pars[2]))
					
			df.append(pd.Series(peak_area[0]))
		
		fitted_area = pd.concat(df)
			
		# Finalise dataframe
		tmp = pd.DataFrame(peak_pos).reset_index(drop=True)
		tmp2 = pd.DataFrame(fitted_area, columns=['Peak Property']).reset_index(drop=True)
		processed_ir_data = pd.concat([tmp,tmp2],axis=1)
		
		# Add the method
		processed_ir_data['Method'] = 'fitted area'
		
		return processed_ir_data


	def fitted_area_sp(self, peak_threshold, residence_time, peak_of_interest, time_adjust_before = 0, time_adjust_after = 0, picked_peak = 0, p0 = [5, -1.5, 2]):
		
		"""Examine a single peak and fit a Gaussian curve 
		
		Parameter
		---------
		peak_threshold : Peak threshold
		residence_time: Residence time of a single peak
		peak_of_interest: Single wavelength to integrate (ie 'Peak at 1704 cm-1')
		time_adjust_before: Time to remove from front of peak
		time_adjust_after: Time to add after peak
		picked_peak: index of peak to examine
		p0: Initial guesses for gaussian function, if blank is [5, -1.5, 2]
		
		Returns
		-------
		Plot of fitted area	for a single peak
		"""
		
		# Requires prominence data (don't use prominence function as don't want returned processed_ir_data)
		# Find single peak
		x = self.ir_data[peak_of_interest]

		# Find peak height and position
		peaks = find_peaks(x, prominence = peak_threshold)
		peak_prominence = peaks[1]
		peak_pos = self.ir_data['Relative Time'][peaks[0]]
		
		# Make it a list
		list_of_peaks = list(peak_pos)
			
		# Find a single peak - Adjusted slightly to account for tailing
		single_peak = self.ir_data.loc[(self.ir_data['Relative Time'] >= list_of_peaks[picked_peak] - (residence_time/2 - time_adjust_before)) 
									 & (self.ir_data['Relative Time'] <= list_of_peaks[picked_peak] + (residence_time/2 + time_adjust_after))]

		# Fit Gaussian 
		X = np.arange(0, len(single_peak))
		pars, cov = curve_fit(f=self.gaussian, xdata=X, ydata=single_peak[peak_of_interest], p0 = p0, bounds=(-np.inf, np.inf))
		
		# Plot curves
		plt.scatter(X, single_peak[peak_of_interest], s=20, label='Data')
		x = np.linspace(0, len(single_peak), 60)
		plt.plot(x, self.gaussian(x, *pars), linestyle='--', linewidth=2, color='black')
		plt.fill_between(x, self.gaussian(x, *pars), color='blue', alpha = 0.3) # fill the integrated area
		plt.xlabel("Index (see above chart for corresponding time)")
		plt.ylabel("Peak Area")
		plt.show()		


	def exp_area_sp(self, peak_threshold, residence_time, peak_of_interest, time_adjust_before = 0, time_adjust_after = 0, picked_peak = 0):

		"""Examine a single peak and integrate it 
		
		Parameter
		---------
		peak_threshold : Peak threshold
		residence_time: Residence time of a single peak
		peak_of_interest: Single wavelength to integrate (ie 'Peak at 1704 cm-1')
		time_adjust_before: Time to remove from front of peak
		time_adjust_after: Time to add after peak
		picked_peak: index of peak to examine
		
		Returns
		-------
		Plot of experimental area for a single peak	
		"""
		
		x = self.ir_data[peak_of_interest]

		# Find peak height and position
		peaks = find_peaks(x, prominence = peak_threshold)
		peak_prominence = peaks[1]
		peak_pos = self.ir_data['Relative Time'][peaks[0]]

		# Make it a list
		list_of_peaks = list(peak_pos)

		# Find a single peak
		single_peak = self.ir_data.loc[(self.ir_data['Relative Time'] >= list_of_peaks[picked_peak] - (residence_time/2 - time_adjust_before)) 
								& (self.ir_data['Relative Time'] <= list_of_peaks[picked_peak] + (residence_time/2 + time_adjust_after))]

		# Find peak area of the experimental data
		exp_peak_area = np.trapz(single_peak[peak_of_interest])

		# Plot data
		plt.scatter(single_peak['Relative Time'], single_peak[peak_of_interest], label='Data')

		# Plot start and end points of integration
		int_thresh_x = (single_peak['Relative Time'].iloc[0], single_peak['Relative Time'].iloc[-1])
		int_thresh_y = (single_peak[peak_of_interest].iloc[0], single_peak[peak_of_interest].iloc[-1])

		plt.scatter(int_thresh_x, int_thresh_y, marker = 'x', s=100, color = 'black')

		# Plot Prominence point
		plt.scatter(peak_pos.iloc[picked_peak], single_peak[single_peak['Relative Time'] == peak_pos.iloc[picked_peak]][peak_of_interest], color = 'red', s = 50)

		# Fill integrated area
		plt.fill_between(single_peak['Relative Time'], single_peak[peak_of_interest],    
			#where = ((single_peak['Relative Time'] >= int_thresh_x[0]) & (single_peak['Relative Time'] <= int_thresh_x[1])),
			color='blue', alpha = 0.3)

		plt.xlabel("Time (min)")
		plt.ylabel("Peak Area")
		plt.show()


	def area_finder(self, residence_time):
				
		# Calculate area
		area = Peaks.exp_area(self, self.prom_thresh, residence_time, self.peak_of_interest, time_adjust_before = 0, time_adjust_after = 0)
		
		df = []

		for var in range(0, self.no_reactions * self.points_per_reaction, self.points_per_reaction):
			each_r2 = np.corrcoef(area.iloc[1 + var : 10 + var, 0],area.iloc[1 + var : 10 + var, 1])[0,1]
			df.append(each_r2)

		# Take absolute values of r2
		r2 = [abs(ele) for ele in df]

		# Sum all the r2 (keep this negative for minimisation)
		sum_r2 = -sum(r2)

		return sum_r2


	def auto_area(self, peak_of_interest, prom_thresh, no_reactions, points_per_reaction, lower_bound = 1, upper_bound = 10):
		
		"""Find the experimental area threshold.
		
		Parameter
		---------
		peak_of_interest: Single wavelength to integrate (ie 'Peak at 1704 cm-1') 
		lower_bound = 1
		upper_bound = 10
		"""
		
		self.peak_of_interest = peak_of_interest
		self.prom_thresh = prom_thresh
		self.no_reactions = no_reactions
		self.points_per_reaction = points_per_reaction
		
		mins = minimize_scalar(self.area_finder, method='bounded', bounds=(lower_bound, upper_bound), options={'xatol' : 1e-5})
		
		return mins.x

	
	def compare(self, peak_height, peak_threshold, residence_time, peak_of_interest, no_reactions, points_per_reaction, time_adjust_before = 0, time_adjust_after = 0, p0 = [5, -1.5, 2]):
		
		"""
		Compare prominence, height, experimental area, and fitted area
		
		Parameter
		---------
		peak_height: peak height threshold
		peak_threshold: peak prominence threshold
		residence_time: Residence time of a singel peak
		peak_of_interest: Single wavelength to integrate (ie 'Peak at 1704 cm-1')
		no_reaction: Number of reactions performed
		points_per_reaction: Number of SPKA peaks in reaction (including t0)
		time_adjust_before: Time to remove from front of peak
		time_adjust_after: Time to add after peak
		p0: Initial guesses for gaussian function, if blank is [5, -1.5, 2]

		Returns
		--------
		Plot of all normalised peak properties with best fit lines
		final: Dataframe of r2 values for best fit lines
		compare: Dataframe containing the peak properties for prominence, height, experimental area, and fitted area
		"""
	
		# Run each function. Use the relative time 
		prominence = Peaks.prominence(self, peak_threshold, peak_of_interest)
		height = Peaks.height(self, peak_height, peak_of_interest).iloc[:,1]
		exp_area = Peaks.exp_area(self, peak_threshold, residence_time, peak_of_interest, time_adjust_before, time_adjust_after).iloc[:,1]
		fitted_area = Peaks.fitted_area(self, peak_threshold, residence_time, peak_of_interest, time_adjust_before, time_adjust_after).iloc[:,1]
		tmp = [prominence,height,exp_area,fitted_area]
		compare = pd.concat(tmp, axis=1)
		
		# Drop method columns
		compare.drop([col for col in compare.columns if 'Method' in col], axis=1, inplace=True)
		
		# Rename columnswith their methods, define manually
		col_headings = ('Relative Time', 'Prominence', 'Height', 'Experimental Area', 'Fitted Area')
		compare.columns = col_headings
		
		# Normalise the data against the largest value in each column
		df =[]
		for var in range(1,len(compare.columns)):
			tmp = compare.iloc[:,var]/max(compare.iloc[:,var])
			df.append(tmp)
		normalised = pd.concat(df, axis=1)

		# Calculate the r2 for each reaction (not including the t0 peak)		
		# For loop along each reaction
		df3 = []
		for var in range(0, no_reactions * points_per_reaction, points_per_reaction):
			df = []
			
			# For loop along along each peak property
			for var2 in range(0,len(normalised.columns)):
				r2 = np.corrcoef(compare.iloc[1 + var : 10 + var, 0],normalised.iloc[1 + var : 10 + var, var2])[0,1]
				df.append(r2)
				df2 = pd.DataFrame(df)
			df3.append(df2)

		# Must transpose df to be the same orientation as previous dfs
		final = pd.concat(df3, axis=1).T.reset_index(drop=True)
		final.columns = list(normalised.columns)
		final.loc['Sum'] = final.sum()
		final
		
		# Plot all picked peaks
		fig = plt.figure(figsize=(20,6))
		for var in range(0,4):
			plt.scatter(compare['Relative Time'], normalised.iloc[:,var], label = normalised.columns[var])
			plt.legend()
		
		# Keep the line colours the same as the points
		line_colours = ('b','orange','g','r')

		# Create best fit line through each experiment (minus t0)
		for var in range(0, no_reactions * points_per_reaction, points_per_reaction):
			a, b = np.polyfit(compare.iloc[1+var:10+var,0],normalised.iloc[1+var:10+var,:], 1)
			
			# Plot lines
			for var2 in range(0, 4): # Change 4 to number of methods examined
				plt.plot(compare.iloc[1+var:10+var,0], a[var2] * compare.iloc[1+var:10+var,0] + b[var2], color=line_colours[var2])			

		return final, compare


	def compare_no_height(self, peak_threshold, residence_time, peak_of_interest, no_reactions, points_per_reaction, time_adjust_before = 0, time_adjust_after = 0, p0 = [5, -1.5, 2]):
		
		"""
		Compare prominence, experimental area, and fitted area
		
		Parameter
		---------
		peak_threshold: peak prominence threshold
		residence_time: Residence time of a singel peak
		peak_of_interest: Single wavelength to integrate (ie 'Peak at 1704 cm-1')
		no_reaction: Number of reactions performed
		points_per_reaction: Number of SPKA peaks in reaction (including t0)
		time_adjust_before: Time to remove from front of peak
		time_adjust_after: Time to add after peak
		p0: Initial guesses for gaussian function, if blank is [5, -1.5, 2]

		Returns
		--------
		Plot of all normalised peak properties with best fit lines
		final: Dataframe of r2 values for best fit lines
		compare: Dataframe containing the peak properties for prominence, experimental area, and fitted area
		"""
	
		# Run each function. Use the relative time 
		prominence = Peaks.prominence(self, peak_threshold,peak_of_interest)
		exp_area = Peaks.exp_area(self, peak_threshold, residence_time, peak_of_interest, time_adjust_before, time_adjust_after).iloc[:,1]
		fitted_area = Peaks.fitted_area(self, peak_threshold, residence_time, peak_of_interest, time_adjust_before, time_adjust_after).iloc[:,1]
		tmp = [prominence,exp_area,fitted_area]
		compare = pd.concat(tmp, axis=1)
		
		# Drop method columns
		compare.drop([col for col in compare.columns if 'Method' in col], axis=1, inplace=True)
		
		# Rename columns with their methods, define manually
		col_headings = ('Relative Time', 'Prominence', 'Experimental Area', 'Fitted Area')
		compare.columns = col_headings
		
		# Normalise the data against the largest value in each column
		df =[]
		for var in range(1,len(compare.columns)):
			tmp = compare.iloc[:,var]/max(compare.iloc[:,var])
			df.append(tmp)
		normalised = pd.concat(df, axis=1)

		# Calculate the r2 for each reaction (not including the t0 peak)		
		# For loop along each reaction
		df3 = []
		for var in range(0, no_reactions * points_per_reaction, points_per_reaction):
			df = []
			
			# For loop along along each peak property
			for var2 in range(0,len(normalised.columns)):
				r2 = np.corrcoef(compare.iloc[1 + var : 10 + var, 0],normalised.iloc[1 + var : 10 + var, var2])[0,1]
				df.append(r2)
				df2 = pd.DataFrame(df)
			df3.append(df2)

		# Must transpose df to be the same orientation as previous dfs
		final = pd.concat(df3, axis=1).T.reset_index(drop=True)
		final.columns = list(normalised.columns)
		final.loc['Sum'] = final.sum()
		final
		
		# Plot all picked peaks
		fig = plt.figure(figsize=(20,6))
		for var in range(0,3):
			plt.scatter(compare['Relative Time'], normalised.iloc[:,var], label = normalised.columns[var])
			plt.legend()
		
		# Keep the line colours the same as the points
		line_colours = ('b','g','r')

		# Create best fit line through each experiment (minus t0)
		for var in range(0, no_reactions * points_per_reaction, points_per_reaction):
			a, b = np.polyfit(compare.iloc[1+var : 10+var, 0],normalised.iloc[1+var : 10+var, :], 1)
			
			# Plot lines
			for var2 in range(0, 3): # Change 3 to number of methods examined
				plt.plot(compare.iloc[1+var : 10+var,0], a[var2] * compare.iloc[1+var : 10+var,0] + b[var2], color=line_colours[var2])			

		return final, compare
		

	def plot(self, processed_ir_data, peak_of_interest):
	
		"""Plot a single picked peaks against the raw data"""
		
		# What is the name of the peak
		peak_list = [x for x in list(processed_ir_data) if 'Peak' in x]
		
		# Plot all picked peaks
		fig = plt.figure(figsize=(20,3))
		ax = fig.subplots()
		ax.title.set_text(peak_of_interest)
		ax.plot(self.ir_data['Relative Time'], self.ir_data[peak_of_interest])
		ax.scatter(processed_ir_data.iloc[:,0],processed_ir_data.iloc[:,1], color='r')
		
		#ax.plot(self.ir_data['Relative Time'], self.ir_data[peak_list[1]])
		#ax.scatter(processed_ir_data[peak_list[0]],processed_ir_data[peak_list[1]], color='r')

