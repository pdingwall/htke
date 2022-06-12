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

	"""Process raw IR data"""
	
	def __init__(self, ir_data):
		self.ir_data = ir_data
		
		# Create list of all the IR peaks monitored
		peak_list = [x for x in list(ir_data) if 'Peak' in x]
		self.peak_list = peak_list
		
	def prominence(self, peak_threshold):
	
		"""Find the prominence of peaks of every wavelength in raw ir_data
		
		Parameter:
		----------
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

		processed_ir_data = pd.concat(list_of_dfs,axis=1)

		return processed_ir_data

		
	def gaussian(x, a, b, c):
		return a*np.exp(-np.power(x - b, 2)/(2*np.power(c, 2)))
	
	def exp_area(self, peak_threshold, residence_time, peak_of_interest):
		
		"""Find the areas of the experimental peak data"""
		
		# Requires prominence data (don't use prominence function as don't want returned processed_ir_data)
		for var in self.peak_list:
			
			# Find single peak
			arse = self.ir_data[var]

			# Find peak height and position
			peaks = find_peaks(arse, prominence = peak_threshold)
			peak_prominence = peaks[1]
			peak_pos = self.ir_data['Relative Time'][peaks[0]]	
		
		# Make it a list
		list_of_peaks = list(peak_pos)
		
		experimental_area = []
		
		for var in list_of_peaks:
    
			# Find a single peak
			single_peak = self.ir_data.loc[(self.ir_data['Relative Time'] >= var - (residence_time/2-0.7)) 
										 & (self.ir_data['Relative Time'] <= var + (residence_time/2+0.7))]
			   
			# Find peak area of the experimental data
			exp_peak_area = np.trapz(single_peak[peak_of_interest])
			
			experimental_area.append(exp_peak_area)
	
		#processed_ir_data = pd.concat( Need to finish this
		
		return experimental_area



	# Currently doesn't work
	def fitted_area_sp(self, peak_threshold, residence_time, peak_of_interest, picked_peak=0):
		
		"""Examine a single peak and fit a Gaussian curve for a single wavelength"""
		
		# Requires prominence data (don't use prominence function as don't want returned processed_ir_data)
		for var in self.peak_list:
			
			# Find single peak
			arse = self.ir_data[var]

			# Find peak height and position
			peaks = find_peaks(arse, prominence = peak_threshold)
			peak_prominence = peaks[1]
			peak_pos = self.ir_data['Relative Time'][peaks[0]]
		
		
		# Find a single peak - Adjusted slightly to account for tailing
		single_peak = self.ir_data.loc[(self.ir_data['Relative Time'] >= list(peak_pos)[0] - (residence_time/2-0.7)) 
									& (self.ir_data['Relative Time'] <= list(peak_pos)[0] + (residence_time/2+0.7))]

		# Fit Gaussian 
		X = np.arange(0, len(single_peak))
		#x = self.ir_data[peak_of_interest]
		pars, cov = curve_fit(f=self.gaussian, xdata=X, ydata=single_peak[peak_of_interest], p0=[5, -1.5, 2], bounds=(-np.inf, np.inf))
		
		# Plot curves
		plt.scatter(X, single_peak[peak_of_interest], s=20, label='Data')
		x = np.linspace(0, len(single_peak), 60)
		plt.plot(x, gaussian(x, *pars), linestyle='--', linewidth=2, color='black')
		plt.xlabel("Index (see above chart for corresponding time)")
		plt.ylabel("Peak Area")
		plt.show()
		
		


	def plot(self,processed_ir_data):
	
		"""Plot the picked peaks against the raw data"""
				
		# Plot all picked peaks
		for var in self.peak_list:
			fig = plt.figure(figsize=(20,3))
			ax = fig.subplots()
			ax.title.set_text(var)
			ax.plot(self.ir_data['Relative Time'],self.ir_data[var])
			ax.scatter(processed_ir_data['Relative Time ' + var],processed_ir_data[var], color='r')