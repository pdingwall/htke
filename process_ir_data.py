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
	
	
	def exp_area(self, peak_threshold, residence_time, peak_of_interest, time_adjust_before = 0.7, time_adjust_after = 0.7):
		
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
		tmp2 = pd.DataFrame(experimental_area, columns=['Experimental Area']).reset_index(drop=True)
		processed_ir_data = pd.concat([tmp,tmp2],axis=1)
		
		return processed_ir_data


	def gaussian(self, x, a, b, c):
		return a*np.exp(-np.power(x - b, 2)/(2*np.power(c, 2)))


	def fitted_area(self, peak_threshold, residence_time, peak_of_interest, time_adjust_before = 0.7, time_adjust_after = 0.7, p0 = [5, -1.5, 2]):
		
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
		tmp2 = pd.DataFrame(fitted_area, columns=['Fitted Area']).reset_index(drop=True)
		processed_ir_data = pd.concat([tmp,tmp2],axis=1)
		
		return processed_ir_data


	def fitted_area_sp(self, peak_threshold, residence_time, peak_of_interest, time_adjust_before = 0.7, time_adjust_after = 0.7, picked_peak = 0, p0 = [5, -1.5, 2]):
		
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
		single_peak = self.ir_data.loc[(self.ir_data['Relative Time'] >= list_of_peaks[picked_peak] - (residence_time/2 - time_adjust_before)) 
									 & (self.ir_data['Relative Time'] <= list_of_peaks[picked_peak] + (residence_time/2 + time_adjust_after))]

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
	
	
	def compare(self, peak_height, peak_threshold, residence_time, peak_of_interest, no_reactions, points_per_reaction, time_adjust_before = 0.7, time_adjust_after = 0.7, p0 = [5, -1.5, 2]):
		
		"""
		Compare prominence, height, experimental area, and fitted area
		
		Parameter
		---------
		peak_height: peak height threshold
		peak_threshold: peak prominence threshold
		residence_time: Residence time of a singel peak
		peak_of_interest: Single wavelength to integrate (ie 'Peak at 1704 cm-1')
		time_adjust_before: Time to remove from front of peak
		time_adjust_after: Time to add after peak
		no_reaction: Number of reactions performed
		points_per_reaction: Number of SPKA peaks in reaction (including t0)
		p0: Initial guesses for gaussian function, if blank is [5, -1.5, 2]

		Retturns
		--------
		Plot of all normalised peak properties with best fit lines
		Dataframe of r2 values for best fit lines
		"""
	
		# Run each function. Use the relative time 
		prominence = Peaks.prominence(self, peak_threshold,peak_of_interest)
		height = Peaks.height(self, peak_height,peak_of_interest).iloc[:,1]
		exp_area = Peaks.exp_area(self, peak_threshold, residence_time, peak_of_interest, time_adjust_before, time_adjust_after).iloc[:,1]
		fitted_area = Peaks.fitted_area(self, peak_threshold, residence_time, peak_of_interest, time_adjust_before, time_adjust_after).iloc[:,1]
		tmp = [prominence,height,exp_area,fitted_area]
		compare = pd.concat(tmp, axis=1)
		
		# Normalise the data against the largest value in each column
		df =[]
		for var in range(1,len(compare.columns)):
			tmp = compare.iloc[:,var]/max(compare.iloc[:,var])
			df.append(tmp)
		normalised = pd.concat(df, axis=1)

		# Calculate the r2 for each reaction (not including the t0 peak)
		df3 = []
		# For loop along each reaction
		for var in range(0, no_reactions * points_per_reaction, points_per_reaction):
			df = []
			
			# For loop along along each peak property
			for var2 in range(0,len(normalised.columns)):
				r2 = np.corrcoef(compare.iloc[1+var:10+var,0],normalised.iloc[1+var:10+var,var2])[0,1]
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
			for var2 in range(0, no_reactions):
				plt.plot(compare.iloc[1+var:10+var,0], a[var2] * compare.iloc[1+var:10+var,0] + b[var2], color=line_colours[var2])			

		return final				
	

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