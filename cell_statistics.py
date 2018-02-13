# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 17:32:01 2018

@author: turkowyd
"""

import numpy as np
import pandas as pd
import mistral as ms
import os
import math


#Number of frames in the movie
frames = 10000

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
    
    
# Magic part

size = list()
D3_average = list()
MD_average = list()
MSD1_average = list()
D3_std = list()
MD_std = list()
MSD1_std = list()
D3_median = list()
MD_median = list()
MSD1_median = list()

for pmt in list_of_pmts:
    no_of_particles = len(pmt['Particle'])
    size.append(no_of_particles)
    
    D3 = list()
    MD = list()
    MSD1 = list()
    for i in range(len(pmt['Track'])):
#        print(pmt['Track']['size'])
        if int(pmt['Track']['size'][i]) > 6:
            D3.append(pmt['Track']['D3'][i])
            MD.append(pmt['Track']['MD'][i])
            MSD1.append(pmt['Track']['MSD1'][i])
    D3_average.append(np.mean(D3))
    MD_average.append(np.mean(MD))
    MSD1_average.append(np.mean(MSD1))
    D3_std.append(np.std(D3))
    MD_std.append(np.std(MD))
    MSD1_std.append(np.std(MSD1)) 
    D3_median.append(np.median(D3))
    MD_median.append(np.median(MD))
    MSD1_median.append(np.median(MSD1))

size = np.array(size)

with open('asdf.txt') as f:
    for i in range(len(size)):
        for s, mean, median, std in zip(size, MD_average, MD_median, MD_std):
            f.write(str([s, mean, median, std]))