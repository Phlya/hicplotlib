# -*- coding: utf-8 -*-
"""
Created on Tue Apr 29 17:53:23 2014

@author: ilya
"""
import hicplotlib
import matplotlib.colors as mcolors

a = hicplotlib.HiCPlot()
a.set_chromosomes(['2L', '2R', '3L', '3R', '4', 'X'])
a.set_chromosomes_boundaries([(0, 1151), (1151, 2209), (2209, 3437), 
                              (3437, 4833), (4833, 4901), (4901, 6023)])
a.set_resolution(20000)
a.read_and_combine_two_matrices('/media/ilya/SSD/HiC/Kc167/IC-heatmap-20K.mtx',
                                '/media/ilya/SSD/HiC/S2/IC-heatmap-20K.mtx',
                                ('Kc167', 'S2'))
c = mcolors.ColorConverter().to_rgb
a.make_cmap([c('white'), c('yellow'),
     0.4, c('yellow'), c('red'),
     0.8, c('red'), c('black')])
a.apply_cmap()
a.plot_whole_genome_heatmap()