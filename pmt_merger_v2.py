# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 10:31:40 2018

@author: turkowyd
"""

import pandas as pd
import json
import os
import numpy as np

os.chdir("E:\\Software\\python_script_depository")
import mistral as ms
    
    
  
# Giving the file path    
file_folder = input('Paste the path to the folder: ')
os.chdir(file_folder)

# Searching for pmt files in folder
list_of_cells = []

for f in os.listdir(file_folder):
    if "roi" in f and f.endswith(".pmt"):
        list_of_cells.append(f)

# Loading pmt files into pandas DataFrames
list_of_pmts = list()

for f in list_of_cells:
    temp = ms.load(f)
    list_of_pmts.append(temp)

# Changing DAG and Track_id values to avoid duplications
for i in range(1, len(list_of_pmts)):
    list_of_pmts[i]['Particle']['DAG'] = list_of_pmts[i]['Particle']['DAG'] + max(list_of_pmts[i-1]['Particle']['DAG']) +1
    list_of_pmts[i]['Particle']['Track_id'] = list_of_pmts[i]['Particle']['Track_id'] + max(list_of_pmts[i-1]['Particle']['Track_id']) +1

for i in range(len(list_of_pmts)):
    cell_class = np.zeros((len(list_of_pmts[i]['Particle']['DAG']),1))
    cell_class = cell_class + i
#    cell_class = pd.DataFrame(cell_class, columns = ['cellid'])
#    cell_dict = {'Particle': cell_class}

    list_of_pmts[i]['Particle']['cellid'] = cell_class
    

# Appending multiple DataFrames into a single DataFrame
final_pmt = list_of_pmts[0].copy()

# Changing the indexing (otherwise Mistral won't recognize tracks)
for i in range(1, len(list_of_pmts)):
    final_pmt['Particle'] = final_pmt['Particle'].append(list_of_pmts[i]['Particle'])
    final_pmt['Track'] = final_pmt['Track'].append(list_of_pmts[i]['Track'])
    
    
particle_index_list = list(range(len(final_pmt['Particle'])))
track_index_list = list(range(len(final_pmt['Track'])))    
   

final_pmt['Particle'].index = particle_index_list
final_pmt['Track'].index = track_index_list


# Export DataFrame into pmt file
ms.dump(final_pmt, 'merged.pmt')

print('LOOK FOR MERGED.pmt FILE')
print('LOOK FOR MERGED.pmt FILE')
print('LOOK FOR MERGED.pmt FILE')
