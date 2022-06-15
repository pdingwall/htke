import pandas as pd
import numpy as np
import os
import glob

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