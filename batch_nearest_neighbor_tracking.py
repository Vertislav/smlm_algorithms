# -*- coding: utf-8 -*-
"""
Created on Sun Feb  4 18:46:17 2018

@author: VerTislav
"""

import numpy as np
import pandas as pd
import os
import math
import mistral as ms


#Segmentizer in nm
segmentizer = 1300

# Destination folder
file_folder = input('Please, describe the folder with files:')
os.chdir(file_folder)



#Creates a class for particles
class Particle:
    def __init__(self, x, y, t, i, ident):
        self.x = x
        self.y = y
        self.t = t
        self.i = i
        self.ident = ident
        self.track_id = -1
    
        
    def future(self, particle, segmentizer):
        self.nbrs = []
        self.nn = []
        self.dist = []
        self.index = []        
        for p in particle:
            if p.t == self.t + 1:
                self.nbrs.append(p)
        if not self.nbrs:
            self.nbrs.append(imag_particle)

        for n in self.nbrs:
            self.dist.append(np.hypot(self.x-n.x, self.y-n.y))
            self.index = np.where(self.dist == min(self.dist))
        if min(self.dist) < segmentizer:
            self.nn = self.nbrs[self.index[0][0]]
        else:
            self.nn = ghost_particle
 
           
    def past(self, particle, segmentizer):
        self.past_nbrs = []
        self.past_nn = []
        self.past_dist = []
        self.past_index = []
        for p in particle:
            if p.t == self.t - 1:
                self.past_nbrs.append(p)
        if not self.past_nbrs:
            self.past_nbrs.append(imag_particle)
            
        for n in self.past_nbrs:
            self.past_dist.append(np.hypot(self.x-n.x, self.y-n.y))
            self.past_index = np.where(self.past_dist == min(self.past_dist))
        if min(self.past_dist) < segmentizer:
            self.past_nn = self.past_nbrs[self.past_index[0][0]]
        else:
            self.past_nn = ghost_particle
            
            
# Creates a class for track
class Track:
    def __init__(self, seed):
        self.seed = seed
        self.ident = []
        self.track =[]
        self.track.append(seed)
        
        
    def create_tracks(self):
        if self.track[-1].nn.ident is None:
            pass
        else:
            self.track.append(self.track[-1].nn)
            self.create_tracks()
 
           
    def length(self):
        self.size = len(self.track)

        
    def plotter(self):
        self.x_list = []
        self.y_list = []
        for t in self.track:
            self.x_list.append(t.x)
            self.y_list.append(t.y)

            
    def jd(self):
        self.steps = []
        self.MD = []
        for i in range(len(self.x_list)-1):
            self.steps.append(np.hypot(self.x_list[i+1] - self.x_list[i], self.y_list[i+1]- self.y_list[i]))
        self.MD = np.mean(self.steps)



# Create an imaginary particle (to avoid empty neighbors lists)
imag_particle = Particle(-100000, -100000, 0, 0, 0)

# Ghost particle ;) (particles without a nearest neighbor have a ghost particle partner, avoiding errors again)
ghost_particle = Particle(None, None, None, None, None)



# Searching for localizationfiles in folder
list_of_cells = []

for f in os.listdir(file_folder):
    if "roi" in f and f.endswith(".txt"):
        list_of_cells.append(f)

#List of lists for all localization files        
list_of_rois = list()

for f in list_of_cells:
    list_of_rois.append(np.loadtxt(f))


roi = 1    
for f in list_of_rois:
    # Creating Particle class objects from localization file
    particle =[Particle(f[i, 0], f[i, 1], f[i, 2], f[i, 3], i) for i in range(len(f))]
    
    # Executing a method which will assign to each particle list of particles from the next frame (if exist). Each particle will get a nearest neighbor in next frame (if possible)    
    for p in particle:
        p.future(particle, segmentizer)
        
    # Executing a method which will assign to each particle list of particles from the previous frame (if exist). Each particle will get a nearest neighbor in previous frame (if possible)    
    for p in particle:
        p.past(particle, segmentizer)
        
    # Create Track objects with seeds
    track = [Track(p) for p in particle if p.past_nn.ident is None]
    
    #Gives to track objects unique ident number
    for i in range(len(track)):
        track[i].ident = i
        
    #Create tracks
    for t in track:
        t.create_tracks()
        
    #Calculates size of the track    
    for t in track:
        t.length()
    
    #Creates a lists of x and y coordinates. Useful for plotting and MSD/JD computations
    for t in track:
        t.plotter()
     
     # Computes the jump distance  (JD)  
    for t in track:
        t.jd()
        if math.isnan(t.MD):
            t.MD = 0
            
    for t in track:
        for p in t.track:
            p.track_id = t.ident
            
    export = {'Particle' : pd.DataFrame(columns=['x', 'y', 't', 'Intensity', 'id', 'Track_id']), 'Track': pd.DataFrame(columns=['id', 'MD', 'length', 't_start']) }

    particle_x = []
    particle_y = []
    particle_t = []
    particle_i = []
    particle_ident = []
    particle_track_id = []
    
    for p in particle:
        particle_x.append(p.x)
        particle_y.append(p.y)
        particle_t.append(p.t)
        particle_i.append(p.i)
        particle_ident.append(p.ident)
        particle_track_id.append(p.track_id)
    
    
    export['Particle']['x'] = particle_x
    export['Particle']['y'] = particle_y
    export['Particle']['t'] = particle_t
    export['Particle']['Intensity'] = particle_i
    export['Particle']['id'] = particle_ident
    export['Particle']['Track_id'] = particle_track_id
    
    track_id = []
    track_md = []
    track_size = []
    track_start = []
    
    for t in track:
        track_id.append(t.ident)
        track_md.append(t.MD)
        track_size.append(t.size)
        track_start.append(t.seed.t)
        
    export['Track']['id'] = track_id
    export['Track']['MD'] = track_md
    export['Track']['length'] = track_size
    export['Track']['t_start'] = track_start
    
    # Export DataFrame into pmt file
    ms.dump(export, 'nn_roi' + str(roi) + '.pmt')
    roi += 1
    
    
    
# Merging pmts into a single file

# Searching for pmt files in folder
list_of_cells = []

for f in os.listdir(file_folder):
    if "nn_roi" in f and f.endswith(".pmt"):
        list_of_cells.append(f)
        
# Loading pmt files into pandas DataFrames
list_of_pmts = list()

for f in list_of_cells:
    temp = ms.load(f)
    list_of_pmts.append(temp)

# Changing DAG and Track_id values to avoid duplications
for i in range(1, len(list_of_pmts)):
    list_of_pmts[i]['Particle']['Track_id'] = list_of_pmts[i]['Particle']['Track_id'] + max(list_of_pmts[i-1]['Particle']['Track_id']) +1

for i in range(len(list_of_pmts)):
    cell_class = np.zeros((len(list_of_pmts[i]['Particle']['Track_id']),1))
    cell_class = cell_class + i

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
ms.dump(final_pmt, 'NN_merged.pmt')

print('LOOK FOR NN_MERGED.pmt FILE')
print('LOOK FOR NN_MERGED.pmt FILE')
print('LOOK FOR NN_MERGED.pmt FILE')