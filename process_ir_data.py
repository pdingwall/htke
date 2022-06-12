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
		tmp2 = pd.DataFrame(peak_prominence)#.rename(columns = {'prominences':peak_of_interest})
		processed_ir_data = pd.concat([tmp,tmp2],axis=1)
		
		# Drop the 'base' columns
		processed_ir_data.drop([col for col in processed_ir_data.columns if 'base' in col], axis=1, inplace=True)

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
		tmp2 = pd.DataFrame(peak_height)#.rename(columns = {'heights':peak_of_interest})
		processed_ir_data = pd.concat([tmp,tmp2],axis=1)

		return processed_ir_data
	
	
	def exp_area(self, peak_threshold, residence_time, peak_of_interest):
		
		"""Find the areas of the peaks in the experimental data
		
		Parameter
		---------
		peak_threshold : Peak threshold
		residence_time: Residence time of a single peak
		peak_of_interest: Single wavelength to integrate (ie 'Peak at 1704 cm-1')
		
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
			single_peak = self.ir_data.loc[(self.ir_data['Relative Time'] >= var - (residence_time/2-0.7)) 
										 & (self.ir_data['Relative Time'] <= var + (residence_time/2+0.7))]
			   
			# Find peak area of the experimental data
			exp_peak_area = np.trapz(single_peak[peak_of_interest])
			
			df.append(pd.Series(exp_peak_area))
		
		experimental_area = pd.concat(df)
			
		# Finalise dataframe
		tmp = pd.DataFrame(peak_pos).reset_index(drop=True)
		tmp2 = pd.DataFrame(experimental_area, columns=['Experimental Area']).reset_index(drop=True)
		processed_ir_data = pd.concat([tmp,tmp2],axis=1)
		
		return processed_ir_data


	def gaussian(self, x, a, b, c):
		return a*np.exp(-np.power(x - b, 2)/(2*np.power(c, 2)))


	def fitted_area(self, peak_threshold, residence_time, peak_of_interest, p0 = [5, -1.5, 2]):
		
		"""Find the areas of the fitted idealised Gaussian peaks
		
		Parameter
		---------
		peak_threshold : Peak threshold
		residence_time: Residence time of a single peak
		peak_of_interest: Single wavelength to integrate (ie 'Peak at 1704 cm-1')
		
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
			single_peak = self.ir_data.loc[(self.ir_data['Relative Time'] >= var - (residence_time/2-0.7)) 
										 & (self.ir_data['Relative Time'] <= var + (residence_time/2+0.7))]
			   
			# Fit Gaussian and find parameters
			X = np.arange(0, len(single_peak))
			pars, cov = curve_fit(f=self.gaussian, xdata=X, ydata=single_peak[peak_of_interest], p0=[5, -1, 2], bounds=(-np.inf, np.inf))

			# Find area (function, starting point, end point, args=arguments to use in function)
			peak_area = intg.quad(self.gaussian, 0, len(single_peak)-1, args=(pars[0],pars[1],pars[2]))
					
			df.append(pd.Series(peak_area[0]))
		
		fitted_area = pd.concat(df)
			
		# Finalise dataframe
		tmp = pd.DataFrame(peak_pos).reset_index(drop=True)
		tmp2 = pd.DataFrame(fitted_area, columns=['Fitted Area']).reset_index(drop=True)
		processed_ir_data = pd.concat([tmp,tmp2],axis=1)
		
		return processed_ir_data


	def fitted_area_sp(self, peak_threshold, residence_time, peak_of_interest, picked_peak = 0, p0 = [5, -1.5, 2]):
		
		"""Examine a single peak and fit a Gaussian curve 
		
		Parameter
		---------
		peak_threshold : Peak threshold
		residence_time: Residence time of a single peak
		peak_of_interest: Single wavelength to integrate (ie 'Peak at 1704 cm-1')
		picked_peak: index of peak to examine
		p0: Initial guesses for gaussian function, if blank is [5, -1.5, 2]
		
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
			
		# Find a single peak - Adjusted slightly to account for tailing
		single_peak = self.ir_data.loc[(self.ir_data['Relative Time'] >= list_of_peaks[picked_peak] - (residence_time/2-0.7)) 
									 & (self.ir_data['Relative Time'] <= list_of_peaks[picked_peak] + (residence_time/2+0.7))]

		# Fit Gaussian 
		X = np.arange(0, len(single_peak))
		pars, cov = curve_fit(f=self.gaussian, xdata=X, ydata=single_peak[peak_of_interest], p0 = p0, bounds=(-np.inf, np.inf))
		
		# Plot curves
		plt.scatter(X, single_peak[peak_of_interest], s=20, label='Data')
		x = np.linspace(0, len(single_peak), 60)
		plt.plot(x, self.gaussian(x, *pars), linestyle='--', linewidth=2, color='black')
		plt.xlabel("Index (see above chart for corresponding time)")
		plt.ylabel("Peak Area")
		plt.show()		
	
	#def compare(self, peak_threshold, residence_time, peak_of_interest, p0 = [5, -1.5, 2]):
	
		
	def plot(self,processed_ir_data, peak_of_interest):
	
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