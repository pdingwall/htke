import pandas as pd
import numpy as np
import os
import glob
from scipy import sparse
from scipy.sparse.linalg import spsolve

class Data:
    
    """Class to read IR data"""
    
    def read():
        
        """
		Will glob all txt files and return them in a dataframe.
		This is specifically for IR data"""
        
        globbed_files = glob.glob("*.txt")

        frame=[]
        
        for txt in globbed_files:

            tmp = pd.read_csv(txt, sep="\t")

            # Add filename as column (make sure that filename.txt is the experiment number)
            exp_no=os.path.splitext(txt)[0]
            tmp['Exp No']=exp_no

            # Change hh:mm:ss format into seconds
            tmp['Relative Time'] = tmp['Relative Time'].str.split(':').apply(lambda x :int(x[0]) * 3600 + int(x[1]) * 60 + int(x[2]))

            # Then back into more useful minutes
            tmp['Relative Time'] = tmp['Relative Time'].astype(int)/60

            frame.append(tmp)

        ir_data = pd.concat(frame, ignore_index=True) 
        
        return ir_data
        
    def plot(ir_data):
        
        """Visaulise the raw IR data"""
        
        fig = ir_data.plot(x='Relative Time',figsize=(20,6))
                
        return fig
		
		
    def baseline_correction(data, peak_of_interest, lam = 200, p = 0.01, niter = 10):
        
        """
        data: supply ir_data[peak_of_interest]
		peak_of_interest: peak of interest
		lam: Affects resolution, lower is more detailed, higher is more general 
		p: Affects resolution, higher is more detailed
		niter: Affects depth, higher takes baseline into the peak (keep at 10 or below)
        From: https://stackoverflow.com/questions/57350711/baseline-correction-for-spectroscopic-data
        """

        #this is the code for the fitting procedure        
        L = len(data)
        w = np.ones(L)
        D = sparse.diags([1,-2,1], [0,-1,-2], shape=(L,L-2))

        for jj in range(int(niter)):
            W = sparse.spdiags(w, 0, L, L)
            Z = W + lam * D.dot(D.transpose())
            z = spsolve(Z, w*data.astype(np.float64))
            w = p * (data > z) + (1-p) * (data < z)

        return pd.Series(z)