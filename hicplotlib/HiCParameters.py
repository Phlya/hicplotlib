# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 21:30:09 2015

@author: ilya
"""
import numpy as np
from math import ceil

class HiCParameters(object):
    def __init__(self, resolution=None, chromosomes=None, lengths=None,
                 boundaries=None):
        self.resolution = None
        self.chromosomes = None
        self.lengths_bp = None
        self.starts_bp = None
        self.ends_bp = None
        self.lengths = None
        self.starts = None
        self.ends = None
        self.boundaries_bp = None
        self.boundaries = None
        if resolution is not None:
            self.set_resolution(resolution)
        if chromosomes is not None:
            self.set_chromosomes(chromosomes)
        if lengths is not None:
            self.set_chromosomes_lengths(lengths)
        if boundaries is not None:
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
        if self.resolution is not None:
            return self.calculate_boundaries()
        return self.boundaries_bp

    def __make_boundaries_bp_dict(self):
        self.boundaries_bp_dict = dict()
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
        boundaries_bp = ([(start * resolution, end * resolution) for start, end
                                                                 in boundaries])
        starts_bp = [i[0] for i in boundaries_bp]
        ends_bp = [i[1] for i in boundaries_bp]
        return boundaries_bp, starts_bp, ends_bp

    def extract_settings(self, obj):
        '''
        Gets all the needed information from a HiCParameters object. Uses
        obj.settings if no arguments are supplied, if settings are provided
        (as a HiCParameters object), sets it to obj.settings and applies it.
        '''
        obj.resolution = obj.settings.resolution
        obj.chromosomes = obj.settings.chromosomes
        if obj.settings.boundaries is None:
            obj.settings.calculate_boundaries()
        obj.boundaries = obj.settings.boundaries
        if obj.settings.lengths_bp is not None:
            obj.lengths_bp = obj.settings.lengths_bp
            obj.boundaries_bp = obj.settings.boundaries_bp