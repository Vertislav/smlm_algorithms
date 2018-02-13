# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 21:55:10 2018

@author: VerTislav
"""


import pandas as pd
import os
import numpy as np
import mistral as ms
import matplotlib.pyplot as plt
import math
import glob

class Particle:
    def __init__(self, x, y, t, Intensity, DAG, p_app, t_mod, Track_id):
        self.x = x
        self.y = y
        self.t = t
        self.Intensity = Intensity
        self.DAG = DAG
        self.p_app = p_app
        self.t_mod = t_mod
        self.Track_id = Track_id



class Track:
    def __init__(self, D3, MD, MSD1, R2_3, t_start, t_start_mod, size):
        self.D3 = D3
        self.MD = MD
        self.MSD1 = MSD1
        self.R2_3 = R2_3
        self.t_start = t_start
        self.t_start_mod = t_start_mod
        self.size = size


        
class Cell:
    def __init__(self):
        self.origin = []
        self.particles = []
        self.x_pos = []
        self.y_pos = []
        self.new_x = []
        self.new_y = []

        
    def rotate(self, origin, px, py, angle):
        ox, oy = origin      
        qx = ox + math.cos(-np.deg2rad(angle)) * (px - ox) - math.sin(-np.deg2rad(angle)) * (py - oy)
        qy = oy + math.sin(-np.deg2rad(angle)) * (px - ox) + math.cos(-np.deg2rad(angle)) * (py - oy)
        return qx, qy

    
    def load(self, particle):
        for p in particle:
            self.particles.append(p)
            self.x_pos.append(p.x)
            self.y_pos.append(p.y)

            
    def find_centroid(self):
        self.origin = [np.mean(self.x_pos), np.mean(self.y_pos)]

            
    def rotate_grid_search(self, start = 0, stop = 180, jump = 1):
        self.rotation_matrix_x = np.zeros((len(self.x_pos), int(stop/jump)))
        self.rotation_matrix_y = np.zeros((len(self.y_pos), int(stop/jump)))
        self.std_y = []
        for angle in range(start, stop, jump):
            self.rotation_matrix_x[:, int(angle/jump)], self.rotation_matrix_y[:, int(angle/jump)] = self.rotate(self.origin, self.x_pos[:], self.y_pos[:], angle)
        for i in range(int(stop/jump)):
            self.std_y.append(np.std(self.rotation_matrix_y[:,i]))
        self.rotate_angle = np.where(self.std_y == min(self.std_y))
        for i in range(len(self.rotation_matrix_x)):
            self.new_x.append(self.rotation_matrix_x[i, self.rotate_angle[0][0]])
            self.new_y.append(self.rotation_matrix_y[i, self.rotate_angle[0][0]])

            
    def norm_positions(self):
        self.new_origin = [np.mean([np.min(self.new_x), np.max(self.new_x)]), np.mean([np.min(self.new_y), np.max(self.new_y)])]
        self.norm_x = (self.new_x - self.new_origin[0]) / (np.max(self.new_x - self.new_origin[0]))
        self.norm_y = (self.new_y - self.new_origin[1]) / (np.max(self.new_y - self.new_origin[1]))
        self.norm_origin = [0,0]

        
    def reflect(self):
        # Vertical rotation
        self.vert_x = self.norm_x
        self.vert_y = -(self.norm_y)
        # Horizontal rotation
        self.horiz_x = -(self.norm_x)
        self.horiz_y = self.norm_y
        # Vertical and horizontal rotation
        self.verthor_x = -(self.norm_x)
        self.verthor_y = -(self.norm_y)


            
class SuperCell:
    def __init__(self, particle, track, cell, particle_array, track_array):
        self.particle = particle
        self.track = track
        self.cell = cell
        self.particle_array = particle_array
        self.track_array = track_array

        
    def noname(self):
        self.norm_particle =[Particle(self.cell.norm_x[i], self.cell.norm_y[i], self.particle_array[i, 0], self.particle_array[i, 1], self.particle_array[i, 2], self.particle_array[i, 3], self.particle_array[i, 4], self.particle_array[i, 5]) for i in range(len(self.cell.norm_x))]
        self.norm_particle_vert =[Particle(self.cell.vert_x[i], self.cell.vert_y[i], self.particle_array[i, 0], self.particle_array[i, 1], self.particle_array[i, 2], self.particle_array[i, 3], self.particle_array[i, 4], self.particle_array[i, 5]) for i in range(len(self.cell.vert_x))]
        self.norm_particle_horiz =[Particle(self.cell.horiz_x[i], self.cell.horiz_y[i], self.particle_array[i, 0], self.particle_array[i, 1], self.particle_array[i, 2], self.particle_array[i, 3], self.particle_array[i, 4], self.particle_array[i, 5]) for i in range(len(self.cell.horiz_x))]
        self.norm_particle_verthor =[Particle(self.cell.verthor_x[i], self.cell.verthor_y[i], self.particle_array[i, 0], self.particle_array[i, 1], self.particle_array[i, 2], self.particle_array[i, 3], self.particle_array[i, 4], self.particle_array[i, 5]) for i in range(len(self.cell.verthor_x))]

        
    def export(self):
        # Creating an export dictionaries
        self.export_norm = {'Particle' : pd.DataFrame(columns = ['x', 'y', 't', 'Intensity', 'DAG', 'p_app', 't_mod', 'Track_id']),
                           'Track' : pd.DataFrame(columns = ['D3', 'MD', 'MSD1', 'R2_3', 't_start', 't_start_mod', 'size'])}
            
        self.export_vert = {'Particle' : pd.DataFrame(columns = ['x', 'y', 't', 'Intensity', 'DAG', 'p_app', 't_mod', 'Track_id']),
                           'Track' : pd.DataFrame(columns = ['D3', 'MD', 'MSD1', 'R2_3', 't_start', 't_start_mod', 'size'])}
            
        self.export_horiz = {'Particle' : pd.DataFrame(columns = ['x', 'y', 't', 'Intensity', 'DAG', 'p_app', 't_mod', 'Track_id']),
                           'Track' : pd.DataFrame(columns = ['D3', 'MD', 'MSD1', 'R2_3', 't_start', 't_start_mod', 'size'])}
            
        self.export_vert_hor = {'Particle' : pd.DataFrame(columns = ['x', 'y', 't', 'Intensity', 'DAG', 'p_app', 't_mod', 'Track_id']),
                           'Track' : pd.DataFrame(columns = ['D3', 'MD', 'MSD1', 'R2_3', 't_start', 't_start_mod', 'size'])}
            
        self.list_reflections = [self.norm_particle, self.norm_particle_vert, self.norm_particle_horiz, self.norm_particle_verthor]
        self.list_exports = [self.export_norm, self.export_vert, self.export_horiz, self.export_vert_hor]
            
        for cell, cell_dict in zip(self.list_reflections, self.list_exports):
            cell_dict['Particle']['x'] = [p.x for p in cell]
            cell_dict['Particle']['y'] = [p.y for p in cell]
            cell_dict['Particle']['t'] = [p.t for p in cell]
            cell_dict['Particle']['Intensity'] = [p.Intensity for p in cell]
            cell_dict['Particle']['DAG'] = [p.DAG for p in cell]
            cell_dict['Particle']['p_app'] = [p.p_app for p in cell]
            cell_dict['Particle']['t_mod'] = [p.t_mod for p in cell]
            cell_dict['Particle']['Track_id'] = [p.Track_id for p in cell]
                
        for t in self.track:
            self.d3 = t.D3
            self.md = t.MD
            self.msd1 = t.MSD1
            self.r2_3 = t.R2_3
            self.t_start = t.t_start
            self.t_start_mod = t.t_start_mod
            self.size = t.size
                
        for cell_dict in self.list_exports:
            cell_dict['Track']['D3'] = self.track_array[:,0]
            cell_dict['Track']['MD'] = self.track_array[:,1]
            cell_dict['Track']['MSD1'] = self.track_array[:,2]
            cell_dict['Track']['R2_3'] = self.track_array[:,3]
            cell_dict['Track']['t_start'] = self.track_array[:,4]
            cell_dict['Track']['t_start_mod'] = self.track_array[:,5]
            cell_dict['Track']['size'] = self.track_array[:,6]
            
        # Changing DAG and Track_id values to avoid duplications
        for i in range(1, len(self.list_exports)):
            self.list_exports[i]['Particle']['DAG'] = self.list_exports[i]['Particle']['DAG'] + max(self.list_exports[i-1]['Particle']['DAG']) +1
            self.list_exports[i]['Particle']['Track_id'] = self.list_exports[i]['Particle']['Track_id'] + max(self.list_exports[i-1]['Particle']['Track_id']) +1
        
        # Adding a cellid attribute to distinguish rotations
        for i in range(len(self.list_exports)):
            self.cell_class = np.zeros((len(self.list_exports[i]['Particle']['DAG']),1))
            self.cell_class = self.cell_class + i
            self.list_exports[i]['Particle']['cellid'] = self.cell_class
        
        # Appending multiple DataFrames into a single DataFrame
        self.final_pmt = self.list_exports[0].copy()
        
        # Changing the indexing (otherwise Mistral won't recognize tracks)
        for i in range(1, len(self.list_exports)):
            self.final_pmt['Particle'] = self.final_pmt['Particle'].append(self.list_exports[i]['Particle'])
            self.final_pmt['Track'] = self.final_pmt['Track'].append(self.list_exports[i]['Track'])
            
        self.particle_index_list = list(range(len(self.final_pmt['Particle'])))
        self.track_index_list = list(range(len(self.final_pmt['Track'])))  
        
        self.final_pmt['Particle'].index = self.particle_index_list
        self.final_pmt['Track'].index = self.track_index_list

                      
    def vibrio(self):
        self.left_orientation = {'Particle' : pd.DataFrame(columns = ['x', 'y', 't', 'Intensity', 'DAG', 'p_app', 't_mod', 'Track_id']),
                           'Track' : pd.DataFrame(columns = ['D3', 'MD', 'MSD1', 'R2_3', 't_start', 't_start_mod', 'size'])}
        self.right_orientation = {'Particle' : pd.DataFrame(columns = ['x', 'y', 't', 'Intensity', 'DAG', 'p_app', 't_mod', 'Track_id']),
                           'Track' : pd.DataFrame(columns = ['D3', 'MD', 'MSD1', 'R2_3', 't_start', 't_start_mod', 'size'])}
        
        # Creating an export dictionaries
        self.export_norm = {'Particle' : pd.DataFrame(columns = ['x', 'y', 't', 'Intensity', 'DAG', 'p_app', 't_mod', 'Track_id']),
                           'Track' : pd.DataFrame(columns = ['D3', 'MD', 'MSD1', 'R2_3', 't_start', 't_start_mod', 'size'])}
            
        self.export_vert = {'Particle' : pd.DataFrame(columns = ['x', 'y', 't', 'Intensity', 'DAG', 'p_app', 't_mod', 'Track_id']),
                           'Track' : pd.DataFrame(columns = ['D3', 'MD', 'MSD1', 'R2_3', 't_start', 't_start_mod', 'size'])}
            
        self.export_horiz = {'Particle' : pd.DataFrame(columns = ['x', 'y', 't', 'Intensity', 'DAG', 'p_app', 't_mod', 'Track_id']),
                           'Track' : pd.DataFrame(columns = ['D3', 'MD', 'MSD1', 'R2_3', 't_start', 't_start_mod', 'size'])}
            
        self.export_vert_hor = {'Particle' : pd.DataFrame(columns = ['x', 'y', 't', 'Intensity', 'DAG', 'p_app', 't_mod', 'Track_id']),
                           'Track' : pd.DataFrame(columns = ['D3', 'MD', 'MSD1', 'R2_3', 't_start', 't_start_mod', 'size'])}
            
        self.list_reflections = [self.norm_particle, self.norm_particle_vert, self.norm_particle_horiz, self.norm_particle_verthor]
        self.list_exports = [self.export_norm, self.export_vert, self.export_horiz, self.export_vert_hor]
            
        for cell, cell_dict in zip(self.list_reflections, self.list_exports):
            cell_dict['Particle']['x'] = [p.x for p in cell]
            cell_dict['Particle']['y'] = [p.y for p in cell]
            cell_dict['Particle']['t'] = [p.t for p in cell]
            cell_dict['Particle']['Intensity'] = [p.Intensity for p in cell]
            cell_dict['Particle']['DAG'] = [p.DAG for p in cell]
            cell_dict['Particle']['p_app'] = [p.p_app for p in cell]
            cell_dict['Particle']['t_mod'] = [p.t_mod for p in cell]
            cell_dict['Particle']['Track_id'] = [p.Track_id for p in cell]
                
        for t in self.track:
            self.d3 = t.D3
            self.md = t.MD
            self.msd1 = t.MSD1
            self.r2_3 = t.R2_3
            self.t_start = t.t_start
            self.t_start_mod = t.t_start_mod
            self.size = t.size
                
        for cell_dict in self.list_exports:
            cell_dict['Track']['D3'] = self.track_array[:,0]
            cell_dict['Track']['MD'] = self.track_array[:,1]
            cell_dict['Track']['MSD1'] = self.track_array[:,2]
            cell_dict['Track']['R2_3'] = self.track_array[:,3]
            cell_dict['Track']['t_start'] = self.track_array[:,4]
            cell_dict['Track']['t_start_mod'] = self.track_array[:,5]
            cell_dict['Track']['size'] = self.track_array[:,6]
            
        self.left_orientation = self.export_norm
        self.left_orientation['Particle'].append(self.export_vert['Particle'])
        self.left_orientation['Track'].append(self.export_vert['Track'])
        
        self.right_orientation = self.export_horiz
        self.right_orientation['Particle'].append(self.export_vert_hor['Particle'])
        self.right_orientation['Track'].append(self.export_vert_hor['Track'])


        
class DensPlot:
    def __init__(self, final_pmt, track_length = 6, md_thres = 300):
        self.final_pmt = final_pmt
        self.track_length = track_length
        self.track_index = np.array(self.final_pmt['Track'].index)
        self.track_MD = np.array(self.final_pmt['Track']['MD'])
        self.track_D3 = np.array(self.final_pmt['Track']['D3'])
        self.track_size = np.array(self.final_pmt['Track']['size'])
        self.particle_trackid = np.array(self.final_pmt['Particle']['Track_id'])
        self.particle_x = np.array(self.final_pmt['Particle']['x'])
        self.particle_y = np.array(self.final_pmt['Particle']['y'])
        self.particle_D3 = []
        self.particle_MD = []
        self.particle_size = []
        self.md_thres = md_thres
        
        for i in range(len(self.particle_x)):
            self.particle_D3.append(self.track_D3[int(self.particle_trackid[i])])
            self.particle_MD.append(self.track_MD[int(self.particle_trackid[i])])
            self.particle_size.append(self.track_size[int(self.particle_trackid[i])])
            
        self.particle_x = np.reshape(self.particle_x, (len(self.particle_x), 1))
        self.particle_y = np.reshape(self.particle_y, (len(self.particle_y), 1))
        self.particle_D3 = np.reshape(self.particle_D3, (len(self.particle_D3), 1))
        self.particle_MD = np.reshape(self.particle_MD, (len(self.particle_MD), 1))
        self.particle_size = np.reshape(self.particle_size, (len(self.particle_size), 1))       
        self.particle_array = np.concatenate((self.particle_x, self.particle_y, self.particle_D3, self.particle_MD, self.particle_size), axis = 1)
            
    def long_tracks(self):
        self.tracks_long = []
        for p in self.particle_array:
            if p[4] > self.track_length:
                self.tracks_long.append(p)
        if len(self.tracks_long) == 0:
            self.tracks_long.append(self.particle_array[-1])
                   
        self.tracks_long = np.array(self.tracks_long)
        
    def slow_track(self):
        self.tracks_slow = []
        for t in self.tracks_long:
            if t[3] < self.md_thres:
                self.tracks_slow.append(t)
        if len(self.tracks_slow) == 0:
            self.tracks_slow.append(self.tracks_long[-1])
            
        self.tracks_slow = np.array(self.tracks_slow)
        
    def heatmap(self, bin_size = [40,100]): 
        self.bin_size = bin_size
        self.heatmap_array = np.zeros((self.bin_size[0], self.bin_size[1]))
        self.hist2d_weights = np.array(np.histogram2d(self.tracks_long[:,0], self.tracks_long[:,1], bins=self.bin_size, weights=self.tracks_long[:,3], range=[[-1,1],[-1,1]]))
        self.hist2d_counts = np.array(np.histogram2d(self.tracks_long[:,0], self.tracks_long[:,1], bins=self.bin_size, range=[[-1,1],[-1,1]]))
        self.hist2d_MD_average = self.hist2d_weights[0]/ self.hist2d_counts[0]
#        self.dataframe_MD_average = pd.DataFrame(self.hist2d_MD_average, columns = self.hist2d_counts[2][1:], index = self.hist2d_counts[1][1:])

    def slow_heatmap(self, bin_size = [40,100]):
        self.bin_size = bin_size
        self.slow_heatmap_array = np.zeros((self.bin_size[0], self.bin_size[1]))
        self.slow_hist2d_weights = np.array(np.histogram2d(self.tracks_slow[:,0], self.tracks_slow[:,1], bins=self.bin_size, weights=self.tracks_slow[:,3], range=[[-1,1],[-1,1]]))
        self.slow_hist2d_counts = np.array(np.histogram2d(self.tracks_slow[:,0], self.tracks_slow[:,1], bins=self.bin_size, range=[[-1,1],[-1,1]]))
        self.slow_hist2d_MD_average = self.slow_hist2d_weights[0]/ self.slow_hist2d_counts[0]
        
        
        
        
# Giving the file path    
file_folder = input('Paste the path to the folder: ')
os.chdir(file_folder)

# Searching for pmt files in folder
all_pmts = list()
for filename in glob.glob(str(file_folder) + "\**\*.pmt", recursive = True):
    all_pmts.append(filename)
    
list_of_cells = list()
for pmt in all_pmts:
    if "roi" in pmt and pmt.endswith(".pmt"):
        list_of_cells.append(pmt)
        
# Loading pmt files into pandas DataFrames
list_of_pmts = list()

for f in list_of_cells:
    list_of_pmts.append(ms.load(f))    
    

# Transferring pandas DataFrames into Particle class, loading all data.
supercells = list()
    
for f in list_of_pmts:
    x = f['Particle']['x']
    y = f['Particle']['y']
    t = f['Particle']['t']
    intensity = f['Particle']['Intensity']
    dag = f['Particle']['DAG']
    p_app = f['Particle']['p_app']
    t_mod = f['Particle']['t_mod']
    track_id = f['Particle']['Track_id']
    for i in range(len(track_id)):
        if math.isnan(track_id[i]) == True:
            track_id[i] = 0
    particle = [Particle(x[i], y[i], t[i], intensity[i], dag[i], p_app[i], t_mod[i], track_id[i]) for i in range(len(x))]
    
    t = np.reshape(t, ((len(t), 1))) 
    intensity = np.reshape(intensity, ((len(intensity), 1)))
    dag = np.reshape(dag, ((len(dag), 1)))
    p_app = np.reshape(p_app, ((len(p_app), 1)))
    t_mod = np.reshape(t_mod, ((len(t_mod), 1))) 
    track_id = np.reshape(track_id, ((len(track_id), 1))) 
    particle_array = np.concatenate([t, intensity, dag, p_app, t_mod, track_id], axis = 1)
    
    
    d3 = f['Track']['D3']
    md = f['Track']['MD']
    msd1 = f['Track']['MSD1']
    r2_3 = f['Track']['R2_3']
    t_start = f['Track']['t_start']
    t_start_mod = f['Track']['t_start_mod']
    size = f['Track']['size']
    track = [Track(d3[i], md[i], msd1[i], r2_3[i], t_start[i], t_start_mod[i], size[i]) for i in range(len(d3))]
    
    d3 = np.reshape(d3, ((len(d3), 1)))
    md = np.reshape(md, ((len(md), 1)))
    msd1 = np.reshape(msd1, ((len(msd1), 1)))
    r2_3 = np.reshape(r2_3, ((len(r2_3), 1)))
    t_start = np.reshape(t_start, ((len(t_start), 1)))
    t_start_mod = np.reshape(t_start_mod, ((len(t_start_mod), 1)))
    size = np.reshape(size, ((len(size), 1)))
    track_array = np.concatenate([d3, md, msd1, r2_3, t_start, t_start_mod, size], axis = 1)
    

    
    # Create a Cell class for particle, to perform all computations, rotations, normalizations, etc.
    cell = Cell()
    cell.load(particle)
    cell.find_centroid()
    cell.rotate_grid_search()
    cell.norm_positions()
    cell.reflect()

    
    supercell = SuperCell(particle, track, cell, particle_array, track_array)
    supercells.append(supercell)
 
for s in supercells:
    s.noname()
    s.export()
    
# Bin size
bin_size = [40, 100]
    
denses = list()

for sc in supercells:
    dense = DensPlot(sc.final_pmt)
    dense.long_tracks()
    dense.slow_track()
    dense.heatmap(bin_size=bin_size)
    dense.slow_heatmap(bin_size = bin_size)
    denses.append(dense)
    
counts = np.zeros((bin_size[0],bin_size[1]))
weights = np.zeros((bin_size[0],bin_size[1]))

slow_counts = np.zeros((bin_size[0],bin_size[1]))
slow_weights = np.zeros((bin_size[0],bin_size[1]))


counts_lists = list()
weigh_lists = list()

for d in denses:
    counts += d.hist2d_counts[0]
    weights += d.hist2d_weights[0]
    slow_counts += d.slow_hist2d_counts[0]
    slow_weights += d.slow_hist2d_weights[0]
    counts_lists.append(d.hist2d_counts [0])
    weigh_lists.append(d.hist2d_weights[0])
    


md_averages = weights / counts
slow_md_averages = slow_weights / slow_counts

pd_column = np.array(range(-int(bin_size[1]/2), int(bin_size[1]/2)))/10
pd_row = np.array(range(-int(bin_size[0]/2), int(bin_size[0]/2)))/50

md_averages_dataframe = pd.DataFrame(md_averages, columns = pd_column, index = pd_row)
slow_md_averages_dataframe = pd.DataFrame(slow_md_averages, columns = pd_column, index = pd_row)  

md_averages_dataframe.to_csv("averaged_heatmap.csv")
slow_md_averages_dataframe.to_csv("averaged_slow_heatmap.csv")
             
plt.imshow(md_averages, cmap="hot", vmin= 0, vmax= 700)
plt.colorbar(label = "Jump Distance [nm]")
plt.xlabel("Long Axis")
plt.ylabel("Short Axis")
plt.savefig('heatmap_100_40bins.png')
plt.close()

plt.imshow(counts, cmap="coolwarm")
plt.colorbar(label = "Number of licalizations")
plt.xlabel("Long Axis")
plt.ylabel("Short Axis")
plt.savefig('heatmap_counts.png')
plt.close()

plt.imshow(slow_md_averages, cmap="hot", vmin= 0, vmax= 300)
plt.colorbar(label = "Jump Distance [nm]")
plt.xlabel("Long Axis")
plt.ylabel("Short Axis")
plt.savefig('slow_heatmap_100_40bins_md300.png')
plt.close()
    
