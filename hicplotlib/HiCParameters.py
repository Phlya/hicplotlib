# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
from math import ceil

class HiCParameters(object):
    '''
    A class to set up Hi-C experiment properties for use in other classes:
    resolution, chromosome names and lengths. Has some functions to convert
    between nt-based and bin-based coordinates, so both can be obtained from
    the attributes.
    '''
    def __init__(self, resolution=0, chromosomes=[], lengths=[],
                 boundaries=[]):
        self.resolution = 0
        self.chromosomes = []
        self.lengths_bp = []
        self.starts_bp = []
        self.ends_bp = []
        self.lengths = []
        self.starts = []
        self.ends = []
        self.boundaries_bp = []
        self.boundaries_bp_dict = dict()
        self.boundaries = []
        if resolution:
            self.set_resolution(resolution)
        if chromosomes:
            self.set_chromosomes(chromosomes)
        if lengths:
            self.set_chromosomes_lengths(lengths)
        if boundaries:
            self.set_chromosomes_boundaries(boundaries)

    def set_resolution(self, resolution):
        """
        Specify resolution of your data, e.g. 20000 = 20 Kb
        """
        self.resolution = resolution

    def set_chromosomes(self, chromosomes):
        """
        Set names of chromosomes from a list in the order they appear in your
        data.
        """
        self.chromosomes = chromosomes

    def set_chromosomes_boundaries(self, boundaries):
        """
        Set coordinates of chromosomes from a list in the order they appear in
        your data. Coordinates should be supplied as number of bins at used
        resolution!
        """
        self.boundaries = boundaries

    def set_chromosomes_lengths(self, lengths):
        '''
        Set lengths of chromosomes in bp. Accepts a list (or a tuple) of
        numbers, ints or floats. At the same time calculates starts and ends
        (and boundaries, which it returns) in bp in a Hi-C matrix.
        If resolutions was set before, calls calculate_boundaries() and
        returns boundaries in bins. Sets self.lengths_bp, self.starts_bp,
        self.ends_bp and self.boundaries_bp (if resolution was set, also all
        of those in bins).
        '''
        self.lengths_bp = lengths
        self.ends_bp = list(np.cumsum(self.lengths_bp))
        self.starts_bp = list(np.cumsum([0]+self.lengths_bp[:-1]))
        self.boundaries_bp = zip(self.starts_bp, self.ends_bp)
        self.__make_boundaries_bp_dict()
        if self.resolution:
            return self.calculate_boundaries()
        return self.boundaries_bp
    
    def set_chromosomes_from_chrfile(self, chrfile):
        '''
        Set names and lengths of chromosomes in bp taking them from a chrfile. A
        chrfile looks like
        chr1\t100
        chr2\t50
        '''
        lengths_bp = []
        chromosomes = []
        with open(chrfile) as f:
            for line in f.readlines():
                name, length_bp = line.split()
                chromosomes.append(name)
                lengths_bp.append(int(length_bp))
        self.set_chromosomes_lengths(lengths_bp)
        self.set_chromosomes(chromosomes)

    def __make_boundaries_bp_dict(self):
        '''
        Make a dict of chrname:(0, chrlength)
        '''
        for name, length in zip(self.chromosomes, self.lengths_bp):
            self.boundaries_bp_dict[name] = (0, length)

    def calculate_boundaries(self):
        '''
        After setting lengths_bp, call this to calculate the boundaries in a Hi-C
        matrix. Set resolution before.
        '''
        self.lengths = [int(ceil(i/self.resolution)) for i in self.lengths_bp]
        self.starts = list(np.cumsum([0]+self.lengths[:-1]))
        self.ends = list(np.cumsum(self.lengths))
        self.boundaries = zip(self.starts, self.ends)
        return self.boundaries

    def calculate_approx_boundaries_bp(self, boundaries=None, resolution=None):
        '''
        Calculate boundaries, starts and ends in bp based on boundaries in bins
        and specified resolution. Doesn't set those values to attributes, only
        returns them!
        '''
        if boundaries is None:
            boundaries = self.boundaries
        if resolution is None:
            resolution = self.resolution
        boundaries_bp = ([(start * resolution, end * resolution)
                                  for start, end in boundaries])
        starts_bp = [i[0] for i in boundaries_bp]
        ends_bp = [i[1] for i in boundaries_bp]
        return boundaries_bp, starts_bp, ends_bp

    def extract_settings(self, obj):
        '''
        Updates the object's settings according to self.
        '''
        obj.resolution = self.resolution
        obj.chromosomes = self.chromosomes
        if self.boundaries is None:
            self.calculate_boundaries()
        obj.boundaries = self.boundaries
        if self.lengths_bp is not None:
            obj.lengths_bp = self.lengths_bp
            obj.boundaries_bp = self.boundaries_bp
    
    def rearrange_chromosomes(self, data, newchrorder):
        '''
        Rearrange a data array. Assumes chromosomes are ordered as in 
        self.chromosomes. Rearranges based on newchrorder. Returns new data.
        Also can be used to get only a subset of chromosomes in specified 
        order. Based on http://stackoverflow.com/a/23455019/1304161
        '''
        neworder = [self.chromosomes.index(i) for i in newchrorder]
        def rearrange(l):
            return [l[i] for i in neworder]
        boundaries = np.concatenate([np.arange(l, u) for l, u in rearrange(self.boundaries)])
        return data[np.ix_(boundaries, boundaries)]
        
    def get_chromosome_boundaries(self, name):
        i = self.chromosomes.index(name)
        return self.boundaries[i]
    
    def get_chromosome_pair_boundaries(self, name, name2=None):
        if name is None:
            raise ValueError('Specify at least one chromosome name')
        if name2 is None:
            name2 = name
        start, end = self.get_chromosome_boundaries(name)
        start2, end2 = self.get_chromosome_boundaries(name2)
        return start, end, start2, end2