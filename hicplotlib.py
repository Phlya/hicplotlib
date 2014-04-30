# -*- coding: utf-8 -*-
"""
Created on Tue Apr 29 17:34:20 2014

@author: ilya
"""
from os import path
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import host_subplot, make_axes_locatable
from math import ceil
#from smooth import smooth
import numpy as np
from scipy import ndimage

import track

class HiCPlot:
    def __init__(self):
        self.name = ''
        self.name2 = ''
        self.cmap='jet'
        self.data=None
        self.resolution=1
        self.chromosomes = None
        self.boundaries = None
    def read_matrix(self, matrix_file, name=''):
        self.data = np.loadtxt(matrix_file)
        if name is not None:
            self.name = name
    def read_and_combine_two_matrices(self, matrixfile1, matrixfile2, 
                                      names=None):
        self.data = np.triu(np.loadtxt(matrixfile1), 1) + \
                    np.tril(np.loadtxt(matrixfile2), 1)
        if names is not None:
            self.name, self.name2 = names
    def get_chromosomes(self, chromosomes):
        self.chromosomes = chromosomes
    def get_chromosomes_boundaries(self, boundaries):
        self.boundaries = boundaries
    def make_cmap(self, seq):
        """
        Return a LinearSegmentedColormap
        seq: a sequence of floats and RGB-tuples. The floats should be 
        increasing and in the interval (0,1).
        Thank you to unutbu at stackoverflow, see
        http://stackoverflow.com/a/16836182/1304161
        """
        seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
        cdict = {'red': [], 'green': [], 'blue': []}
        for i, item in enumerate(seq):
            if isinstance(item, float):
                r1, g1, b1 = seq[i - 1]
                r2, g2, b2 = seq[i + 1]
                cdict['red'].append([item, r1, r2])
                cdict['green'].append([item, g1, g2])
                cdict['blue'].append([item, b1, b2])
        self.cmap = mcolors.LinearSegmentedColormap('CustomMap', cdict)
    def apply_cmap(self, cmap=None):
        if cmap is None:
            cmap=self.cmap
        plt.register_cmap(name=cmap.name, cmap=cmap)
        plt.set_cmap(self.cmap)
        
    def _nice_ticks(self, tick_val, tick_pos):
        value = float(tick_val)/10**6
        return str("%.2f" % value) +'M'

    def _get_title(self):
        if self.name and self.name2:
            return self.name+'/'+self.name2
        elif self.name:
            return self.name
        else:
            return ''
            
    def _make_axlabels(self):
        if self.name and self.name2:
            self.ax.set_xlabel(self.name, size=20)
            self.ax.set_ylabel(self.name2, size=20)
        else:
            return
    
    def plot_whole_genome_heatmap(self, data=None, savepath=False, 
                                  picture_type='svg'):
        if data is None:
            data=self.data
        self.ax = host_subplot(111)
        self.locations = [np.mean(i)*self.resolution for i in self.boundaries]
        length, height = data.shape
        plt.imshow(data, interpolation='none', origin='lower',
               extent=[0, length*self.resolution, 0, length*self.resolution])
        bar = plt.colorbar()
        bar.set_label(label=u'$log_2(число\ чтений)$', size=20)
        self.ax.xaxis.set_tick_params(length=5, direction='out')
        self.ax.yaxis.set_tick_params(length=5, direction='out')
        
        self.ax.xaxis.set_tick_params(direction='out', length=5)
        self.axlabels = self.ax.get_xticklabels()
        for label in self.axlabels:
            label.set_rotation(30)
        self.ax.yaxis.set_tick_params(direction='out', length=5)
        
        self.ax.xaxis.set_major_formatter(mticker.FuncFormatter(self._nice_ticks))
        self.ax.yaxis.set_major_formatter(mticker.FuncFormatter(self._nice_ticks))
        
        self.ax.tick_params(axis='x', which='minor', bottom='on')
        self.ax.tick_params(axis='y', which='minor', left='on')
        self.ax.tick_params(axis='x', which='major', bottom='on', top='on', 
                            labelbottom='on', labelsize=15)
        self.ax.tick_params(axis='y', which='major', left='on', right='on', 
                            labelleft='on', labelsize=15)
        self.ax2 = self.ax.twin()
        self.ax2.set_xticks(self.locations)
        self.ax2.set_xticklabels(self.chromosomes, fontsize=20)
        self.ax2.set_yticks(self.locations)
        self.ax2.set_yticklabels(self.chromosomes, fontsize=20)
        self.ax2.xaxis.set_tick_params(length=0)
        self.ax2.yaxis.set_tick_params(length=0)
        plt.suptitle(self._get_title(), fontsize=20)
        self._make_axlabels()
        if not savepath:
            plt.show()
        else:
            plt.savefig(savepath, picture_type)
            plt.clf()
            
    def plot_chromosome_pair_heatmap(self, name, name2=None, data=None,
                                    savepath=False, picture_type='svg'):
        if data is None:
            data=self.data
        n = self.chromosomes.index(name)
        start, end = self.boundaries[n]
        if name2 != None and name2 != name:
            n2 = self.chromosomes.index(name2)
            start2, end2 = self.boundaries[n2]
        else:
            name2, start2, end2 = name, start, end
        chrdata = data[start:end, start2:end2]
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.ax.xaxis.set_major_formatter(mticker.FuncFormatter(self._nice_ticks))
        self.ax.yaxis.set_major_formatter(mticker.FuncFormatter(self._nice_ticks))
        self.ax.xaxis.set_tick_params(top='off', direction='out', length=5)
        self.ax.yaxis.set_tick_params(right='off', direction='out', length=5)
        plt.imshow(chrdata, origin='lower', interpolation='none',
                   extent=[0, (end-start)*self.resolution,
                           0, (end2-start2)*self.resolution])
        self.bar = plt.colorbar()
        self.bar.set_label(label=u'$log_2(N\ of\ reads)$', size=20)
        plt.title(name+'-'+name2)
        if not savepath:
            plt.show()
        else:
            plt.savefig(savepath, format=picture_type)
            plt.clf()
        
    def plot_by_chromosome_heatmaps(self, data=None, only_intrachromosome=True,
                                    exclude_names = ('M'),
                                    savepath=False, picture_type='svg'):
        if data is None:
            data=self.data
        for chromosome in self.chromosomes:
            if only_intrachromosome:
                if savepath:
                    filepath = path.join(path.abspath(savepath), 
                                         chromosome+'-'+chromosome+'.svg')
                self.plot_chromosome_pair_heatmap(chromosome, savepath=filepath)
            else:
                for chromosome1 in self.chromosomes:
                    if savepath:
                        filepath = path.join(path.abspath(savepath), 
                                             chromosome+'-'+chromosome1+'.svg')
                    self.plot_chromosome_pair_heatmap(chromosome, chromosome1,
                                                      savepath=filepath)