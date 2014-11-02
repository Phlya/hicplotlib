# -*- coding: utf-8 -*-
"""
Created on Tue Apr 29 17:53:23 2014

@author: ilya
"""
#import matplotlib

#matplotlib.use('SVG')

import hicplotlib
import numpy as np
import matplotlib.pyplot as plt


a = hicplotlib.HiCPlot()
a.set_chromosomes(['2L', '2R', '3L', '3R', '4', 'X'])
a.set_chromosomes_boundaries([(0, 1151), (1151, 2209), (2209, 3437), 
                              (3437, 4833), (4833, 4901), (4901, 6023)])
a.set_resolution(20000)
#a.read_two_matrices(('/media/ilya/SSD/HiC/S2/IC-heatmap-20K.mtx', 
#                     '/media/ilya/SSD/HiC/S2LD/IC-heatmap-20K.mtx'), 
#                    ('S2', 'S2LD'))
a.read_matrix('/home/ilya/Документы/biology/Drosophila_cells_Hi-C/BG3/IC-heatmap-20K.mtx', 'BG3')
a.plot_whole_genome_heatmap(colormap='jet', 
                            diagonal_markers={10**5:'grey', 10**6:'black'})
#a.hide_bads()
#c = mcolors.ColorConverter().to_rgb
#a.make_cmap(
#    [c('white'), c('yellow'),
#     0.4, c('yellow'), c('red'),
#     0.8, c('red'), c('black')])
#a.apply_cmap()
#a.data = a.data/np.sum(a.data)
#a.data2 = a.data2/np.sum(a.data2)
#data = (a.data-a.data2)/(a.data+a.data2)
#absmax = np.max(np.abs(data))
#print data[0:100, 0:100]
#a.data = a._normalize_array(a.data)
#a.data2 = a._normalize_array(a.data2)
#data = np.tril(a.data)+np.triu(a.data2)
#data = (a.data-a.data2)/(a.data+a.data2)
#a.plot_by_chromosome_heatmaps(compare='Triangles', savepath='/home/ilya/Документы/biology/Drosophila_cells_Hi-C/HiC_results/heatmaps/BG3vsOSC/', format='svg', vmin=0)#, vmax=absmax)
#x1, y1 = a.scale_plot(a.data, color='g', plot=True, label='S2')#, stat=np.median)
#x2, y2 = a.scale_plot(a.data2, color='b', plot=True, label='S2LD')#, stat=np.median)

#s2x, s2y, s2opt, s2cov = a.fit(x1, y1, minx=5*10**5, maxx=3*10**6)
#print 'a*e^(-b*x)+c'
#print 'S2 fit:', s2opt
##print zip(s2xfit, s2yfit)
#s2ldx, s2ldy, s2ldopt, s2ldcov = a.fit(x2, y2, minx=5*10**5, maxx=3*10**6)
#print 'S2LD fit:', s2ldopt
##print zip(s2ldxfit, s2ldyfit)
#
#plt.plot(s2x, s2y, 'g--', label='S2 exponential fit')
#plt.plot(s2ldx, s2ldy, 'b--', label='S2LD exponential fit')
#def exp_func(x, a, b, c):
#    return a * np.exp(-b * x) + c
    
#yfit = exp_func(x1, 1, 1, 1)
#plt.plot(x1, yfit, 'g--', label=)


#plt.legend()
#plt.xlim(4*10**4, 10**7)
#plt.ylim(10**-3, 20)
#plt.title('Scale plot')
#plt.grid()
#plt.xlabel('Distance (bp)')
#plt.ylabel('Contact probability')
#a.show_figure()
#a._set_axlabels()
#a.ax.locator_params(nbins=15)
#a.show_figure()
a.save_figure(savepath='/home/ilya/Документы/biology/Drosophila_cells_Hi-C/HiC_results/heatmaps/BG3/all-all.svg', format='svg')