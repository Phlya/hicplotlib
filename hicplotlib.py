# -*- coding: utf-8 -*-
"""
This is a module to help with plotting Hi-C data. Everything is inside
hicplotlib.HiCPlot class.
"""
from __future__ import division
from os import path
#from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import matplotlib.cm as cmap
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import host_subplot#, make_axes_locatable
#from matplotlib.ticker import MultipleLocator
from math import ceil
#from smooth import smooth
import numpy as np
from scipy import ndimage

#import track

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
        if self.resolution is not None:
            return self.calculate_boundaries()
        return self.boundaries_bp

    def calculate_boundaries(self):
        '''
        After setting lengths, call this to calculate the boundaries in a Hi-C
        matrix. Set resolution before.
        '''
        self.lengths = [int(ceil(i/self.resolution)) for i in self.lengths_bp]
        self.starts = list(np.cumsum(self.lengths))
        self.ends = list(np.cumsum([0]+self.lengths[:-1]))
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

class HiCPlot(object):
    def __init__(self, settings=None):
        if settings is not None:
            self.settings = settings
        self.extract_settings()
        self.name = ''
        self.name2 = ''
        self.cmap = cmap.get_cmap()
        self.data = None
        self.data_dict = {}

    def extract_settings(self, settings=None):
        '''
        Gets all the needed information from a HiCParameters object. Uses
        self.settings if no arguments are supplied, if settings are provided
        (as a HiCParameters object), sets it to self.settings and applies it.
        '''
        if settings is not None:
            self.settings = settings
        self.resolution = self.settings.resolution
        self.chromosomes = self.settings.chromosomes
        if self.settings.boundaries is None:
            self.settings.calculate_boundaries()
        self.boundaries = self.settings.boundaries
        if self.settings.lengths_bp is not None:
            self.lengths_bp = self.settings.lengths_bp
            self.boundaries_bp = self.settings.boundaries_bp

    def read_matrix(self, matrix_file, name='', *args, **kwargs):
        """
        Load a matrix from a file. Uses np.loadtxt(), so all additional
        arguments will be passed to it. It is also recommended to provide a
        **name** argument, as it will then appear on graphs and can be used in
        **get_data**.
        """
        self.data = np.loadtxt(matrix_file, *args, **kwargs)
        if name is not None:
            self.name = name
            self.data_dict[name] = self.data

    def read_two_matrices(self, matrix_files, names=None, *args, **kwargs):
        """
        A shortcut to load 2 matrices. Uses **read_matrix**, so all additional
        arguments are passed to it. One can specify names as a list or a tuple.
        """
        self.data, self.data2 = (np.loadtxt(matrix_files[0], *args, **kwargs),
                                 np.loadtxt(matrix_files[1], *args, **kwargs))
        if names is not None:
            self.name, self.name2 = names
            self.data_dict[self.name] = self.data
            self.data_dict[self.name2] = self.data2

    def read_intervals(self, interval_file, chrom=False, names=False, *args,
                       **kwargs):
        '''
        Read an interval file, bed-style, but not necessarily with chromosome
        data (e.g. in genomic coordinates if you concatenate chromosomes). If
        chrom=True, then first column is read as chromosomes. Returns a pandas
        DataFrame object with 'Start' and 'End' columns (and optionally
        'Chromosome' as the first one).
        '''
        import pandas as pd
        if names:
            columns = names
        elif chrom:
            columns = 'Chromosome', 'Start', 'End'
        else:
            columns = 'Start', 'End'
        intervals = pd.read_csv(interval_file, sep='\t', names=columns)
        return intervals

    def get_data(self, name):
        """
        Retrieve data by it's name. Calls self.data_dict[name]
        """
        return self.data_dict[name]

    def _precalculate_TADs_in_array(self, array):
        '''
        Calculates greendale statistics for TADs calculation. These may be
        reused with multiple gammas very fast
        '''
        from greendale import segment
        k = array.sum(axis=0)
        pass_mask = k != 0
        Wcomm = segment.normalized_weights_by_segment(array)
        Wnull = segment.normalized_weights_by_segment(np.outer(k,k))
        return Wcomm, Wnull, pass_mask, len(array)

    def _calculate_TADs(self, parameters, gamma):
        '''
        Calculate TADs based on greendale statistics (parameters) and a
        specified gamma value. Returns a pandas DataFrame with columns 'Start'
        and 'End' with coordinates in bins.
        '''
        import pandas as pd
        from greendale import segment
        Wcomm, Wnull, pass_mask, length = parameters
        starts, scores = segment.potts_segmentation(Wcomm, Wnull, gamma,
                                                    pass_mask=pass_mask)
        pos = np.r_[starts, length]
        domains = zip(pos[:-1], pos[1:])
        domains = np.array([(i[0], i[1]) for i in domains if i[1]-i[0]>2])
        domains = pd.DataFrame(domains, columns=('Start', 'End'))
        return domains

    def find_TADs(self, data=None, gammalist=range(10, 110, 10)):
        '''
        Finds TADs in data with a list of gammas. Returns a pandas DataFrame
        with columns 'Start', 'End' and 'Gamma'. Use genome_intervals_to_chr on
        the returned object to get coordinates in bed-style format and not in
        coordinates in bins of concatenated genome.
        '''
        import pandas as pd
        if data is None:
            data=self.data
        parameters  = self._precalculate_TADs_in_array(data)
        domains = pd.DataFrame()
        for g in gammalist:
            domains_g = self._calculate_TADs(parameters, g)
            domains_g['Gamma'] = g
            domains_g['Start'] *= self.resolution
            domains_g['End'] -= 1
            domains_g['End'] *= self.resolution
            domains = domains.append(domains_g, ignore_index=True)
        domains[['Start', 'End']] = domains[['Start', 'End']].astype(int)
        domains = domains[['Start', 'End', 'Gamma']]
        return domains

    def find_TADs_by_chromosomes(self, data=None, gammadict={}):
        '''
        Åpply TAD finding to each chromosome separately. As required gamma
        varies very much with size of supplied matrix for calculation, you
        should supply different gammas for each chromosome in a
        gammadict{chromosome_name:[gamma1, gamma2, ...],
                  ...}
        Returns a pandas DataFrame with columns 'Chromosome', 'Start', 'End' and
        'Gamma'
        '''
        import pandas as pd
        if data is None:
            data=self.data
        domains = pd.DataFrame()
        for i, chrname in enumerate(self.chromosomes):
            start, end = self.boundaries[i]
            chrdata = data[start:end, start:end]
            parameters  = self._precalculate_TADs_in_array(chrdata)
            domains_chr = pd.DataFrame()
            for g in gammadict[chrname]:
                domains_g = self._calculate_TADs(parameters, g)
                domains_g['Chromosome'] = chrname
                domains_g['Gamma'] = g
                domains_g['Start'] *= self.resolution
                domains_g['End'] -= 1
                domains_g['End'] *= self.resolution
                domains_chr = domains_chr.append(domains_g, ignore_index=True)
            domains = domains.append(domains_chr, ignore_index=True)
        domains[['Start', 'End']] = domains[['Start', 'End']].astype(int)
        domains = domains[['Chromosome', 'Start', 'End', 'Gamma']]
        return domains.reset_index(drop=True)

    def genome_coordinate_to_chr(self, coordinate, coordinates_from_bins=True):
        '''
        Recalculate a coordinate from genome-based to chromosome-based.
        '''
        if coordinates_from_bins:
            boundaries_bp, starts_bp, ends_bp = (
                                self.settings.calculate_approx_boundaries_bp())
        else:
            boundaries_bp = self.boundaries_bp
            starts_bp, ends_bp = self.starts_bp, self.ends_bp
        for i, boundaries in enumerate(boundaries_bp):
            start, end = boundaries
            if start<=coordinate<end:
                chrname = self.chromosomes[i]
                chrcoordinate = coordinate - starts_bp[i]
                break
            elif (i == len(boundaries_bp)-1 and
                            0 < coordinate - ends_bp[i] < self.resolution):
                if self.settings.boundaries is None:
                    self.settings.calculate_boundaries()
                chrname = self.chromosomes[i]
                chrcoordinate = self.boundaries[i][1]
                break
            else:
                continue
        if not chrname:
            raise ValueError('Coordinate out of range')
        return chrname, chrcoordinate

    def _remove_interchr_intervals(self, intervals):
        '''
        Accepts a pandas DataFrame with at least 2 columns 'Start_chromosome'
        and 'End_chromosome', checks all rows to contain the same value in these
        columns and removes the ones, where is is not true. Returns the object,
        identical to supplied intervals, but with only one column 'Chromosome'
        instead  of two chromosome columns.
        '''
        unmatched = intervals['Start_chromosome']!=intervals['End_chromosome']
        if any(unmatched):
            intervals = intervals[np.logical_not(unmatched)]
        intervals.drop('End_chromosome', inplace=True, axis=1)
        intervals.rename(columns={'Start_chromosome':'Chromosome'},
                             inplace=True)
        return intervals

    def genome_intervals_to_chr(self, intervals, remove_crossborder=True):
        '''
        Recalculate coordinates of intervals from genome-based to
        chromosome-based. Intervals is a pandas DataFrame with at least 2
        columns, 'Start' and 'End', with coordinates in a concatenated genome.
        '''
        import pandas as pd
        new_intervals = pd.DataFrame(columns=['Start_chromosome', 'Start',
                                              'End_chromosome', 'End'],
                                            index=intervals.index)
        for i in intervals.index:
            startchr, start = self.genome_coordinate_to_chr(intervals.ix[i,
                                                                    'Start'])
            endchr, end = self.genome_coordinate_to_chr(intervals.ix[i, 'End'])
            new_intervals.iloc[i]=startchr, start, endchr, end
        cols = [col for col in intervals.columns if col not in ['Start', 'End']]
        new_intervals[cols] = intervals[cols]
        if remove_crossborder:
            new_intervals = self._remove_interchr_intervals(new_intervals)
        return new_intervals

    def write_TADs(self, domains, basename):
        for g in set(domains['Gamma']):
            domains_g = domains[domains['Gamma']==g]
            domains_g.to_csv(basename+'_g'+str(g)+'.bed', header=False,
                           index=False, sep='\t', cols=['Chromosome', 'Start',
                                                        'End'])
    def _list_files(self, path):
        '''
        List all files in a directory
        '''
        import os
        files = []
        for name in os.listdir(path):
            if os.path.isfile(os.path.join(path, name)):
                files.append(os.path.join(path, name))
        return files

    def read_TADs(self, path, basename='TADs',
                  names=['Chromosome', 'Start', 'End']):
        from os import path as p
        import pandas as pd
        files = self._list_files(path)
        files = [f for f in files if p.basename(f).startswith(basename)]
        print 'Reading TADs from these files:'
        print files
        domains = pd.DataFrame(columns=names)
        for f in files:
            g = float(f.split('.')[0].split('_g')[1])
            domains_g = self.read_intervals(f, names=names)
            domains_g['Gamma'] = g
            domains = domains.append(domains_g, ignore_index=True)
        return domains

    def make_inter_intervals(self, intervals, shorten_by_resolution=True):
        '''
        Accepts a bed-style pandas DataFrame with columns 'Chromosome', 'Start',
        'End'. Returns the file with the same format, but containing coordinates
        of intervals between specified (doesn't include telomeric regions not
        covered with intervals!). If shorten_by_resolution, subtracts
        self.resolution from ends (useful for intervals, acquired from Hi-C
        data, such as TADs).
        '''
        start = intervals[:-1][['End', 'Chromosome']]
        start = start.reset_index(drop=True).rename(columns={
                                              'Chromosome':'Start_chromosome',
                                              'End':'Start'})

        end = intervals[1:][['Start', 'Chromosome']]
        if shorten_by_resolution:
            end -= self.resolution
        end = end.reset_index(drop=True).rename(columns={
                                                'Chromosome':'End_chromosome',
                                                'Start':'End'})
        inter_intervals = start.merge(end, left_index=True, right_index=True)
        inter_intervals = self._remove_interchr_intervals(inter_intervals)
        return inter_intervals[['Chromosome', 'Start','End']]

    def make_interTADs(self, domains):
        '''
        Makes inter-TADs from TADs DataFrame (columns 'Chromosome', 'Start',
        'End', 'Gamma').
        '''
        import pandas as pd
        interTADs = pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Gamma'])
        for g in set(domains['Gamma']):
            interTADs_g = self.make_inter_intervals(domains[domains['Gamma']==g])
            interTADs_g['Gamma'] = g
            interTADs = interTADs.append(interTADs_g, ignore_index=True)
        return interTADs

    def describe_TADs(self, domains, functions=None):
        '''
        Group TADs gy 'Gamma' and returns statistics by group (applies it to
        'Length' column, which is calculated is not supplied). Includes count,
        np.median, np.mean, np.min, np.max and coverage by default. Coverage
        calculates genome coverage of TADs based on true chromosome lengths. If
        supplied with functions argument, you can add any other functions to
        that statistics.
        Uses DataFrame.groupby().aggregate() methods.
        '''
        if functions is None:
            functions = []
        elif callable(functions):
            functions = [functions]
        elif type(functions) is not list:
            functions = list(functions)
        else:
            pass
        def count(x):
            return len(x)
        def coverage(x):
            return np.sum(x)/self.boundaries_bp[-1][1]
        if 'Length' not in domains.columns:
            domains['Length'] = domains['End']-domains['Start']
            domains['Length'] = domains['Length'].astype(int)
        stats = domains.groupby('Gamma')['Length'].agg([count, np.median, np.mean, np.min,
                                       np.max, coverage] + list(functions))
        return stats

    def plot_TADs_length_distribution(self, domains, *args, **kwargs):
        '''
        Plots size distribution of TADs using DataFrame.hist() method by 'Gamma'
        . Creates a column for each possible size in the data. All arguments are
        passed to the .hist() method. Returns a list of lists of axes.
        '''
        if 'Length' not in domains.columns:
            domains['Length'] = domains['End']-domains['Start']
            domains['Length'] = domains['Length'].astype(int)
        bins = range(-self.resolution/2, max(domains['Length'])+
                                           3*self.resolution/2, self.resolution)
        axes = domains['Length'].hist(by=domains['Gamma'], bins=bins,
                                                                *args, **kwargs)
        for ax in axes.flatten():
            ax.set_xlim([-self.resolution/2,
                         max(domains['Length'])+self.resolution/2])
        return axes

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
        return self.cmap

    def apply_cmap(self, cmap=None):
        """
        Make matplotlib use specified colourmap. By default will use
        **self.cmap**
        """
        if cmap is None:
            cmap = self.cmap
        else:
            plt.set_cmap(cmap)
            return
        plt.register_cmap(name=cmap.name, cmap=cmap)
        plt.set_cmap(cmap)

    def _nice_ticks(self, tick_val, tick_pos):
        '''
        Makes use of labelling of coordinates with M for Megabase and K for
        Kilobase, looks nice.
        Inspired by Nikolay Vyahhi's approach at StackOverflow, see
        http://stackoverflow.com/a/7039989/1304161
        '''
        value = float(tick_val)
        return '%1.0fM' % (value * 1e-6) if value >= 10**6 and \
                                            str(int(value))[-6:-1] == '00000' \
        else '%1.1fM' % (value * 1e-6) if value >= 10**6 and \
                                            str(int(value))[-5:-1] == '0000' \
        else '%1.2fM' % (value * 1e-6) if value >= 10**6 and \
                                            str(int(value))[-4:-1] == '000' \
        else '%1.3fM' % (value * 1e-6) if value >= 10**6 and \
                                            str(int(value))[-3:-1] == '00' \
        else '%1.0fK' % (value * 1e-3) if value < 10**6 and \
                                            str(int(value))[-3:-1] == '00' \
        else '%1.0f' % value

    def _make_title(self):
        if self.name and self.name2:
            return self.name+'/'+self.name2
        elif self.name:
            return self.name
        else:
            return ''

    def _set_axlabels(self):
        if self.name and self.name2:
            self.ax.set_ylabel(self.name, size=15)
            self.ax.set_xlabel(self.name2, size=15)
        else:
            return

    def set_zoom(self, x1, x2, y1, y2):
        """
        Specify visible region of **self.ax**
        """
        self.ax.set_xlim(x1, x2)
        self.ax.set_ylim(y1, y2)
        self.update_figure()

    def show_figure(self, *args, **kwargs):
        """
        Just runs plt.show(). Passes all arguments to it.
        """
        plt.show(*args, **kwargs)

    def update_figure(self, *args, **kwargs):
        """
        Just runs plt.draw(). Passes all arguments to it.
        """
        plt.draw(*args, **kwargs)

    def save_figure(self, savepath, *args, **kwargs):
        """
        Saves figure to **savepath**. Just runs plt.savefig() for that. Passes
        all arguments to it.
        """
        plt.savefig(savepath, *args, **kwargs)

    def clear_figure(self):
        """
        Just runs plt.clf()
        """
        plt.clf()

    def hide_bads(self, data=None, mask_zeroes=False, cmap=None, col='w',
                  alpha=0.0):
        """
        Masks all values that are not finite (NaN, inf, neginf). Can also mask
        zero values. At the same time uses cmap.set_bad to hide masked values.
        """
        if data is None:
            data = self.data
        if cmap is None:
            cmap = self.cmap
        data = np.ma.masked_where(~np.isfinite(data), data)
#        data = np.ma.masked_where(np.isinf(data), data)
#        if substitute == 'min':
#            substitute = np.min(data)
#        data[np.isinf(data)] = substitute
#        data[np.isneginf(data)] = substitute

        if mask_zeroes:
            data = np.ma.masked_equal(data, np.min(data))
#        else:
#            data[data == 0] = substitute
        data = np.ma.array(data, mask=np.isnan(data))
        cmap.set_bad(col, alpha)
#        apply_cmap()
        return data

    def _plot_combined_heatmap(self, data, data2, *args, **kwargs):
        if 'origin' in kwargs:
            if kwargs['origin'] == 'lower':
                data, data2 = data2, data
#        data = self._normalize_array(data)
#        data2 = self._normalize_array(data2)
        data = np.triu(self.hide_bads(data, col=(0,0,0.5), alpha=1))
        data2 = np.tril(self.hide_bads(data2, col=(0,0,0.5), alpha=1))
        plt.imshow(data+data2, *args, **kwargs)
        self.colorbar = plt.colorbar()
        self.colorbar.set_label(u'$log_2(N\ of\ reads)$')
#        data[np.tril_indices_from(data)] = np.ma.masked
#        plt.imshow(data, *args, **kwargs)
#        plt.colorbar(orientation='vertical')

#        data2[np.triu_indices_from(data2)] = np.ma.masked
#        plt.imshow(data2, *args, **kwargs)
#        plt.colorbar(orientation='horizontal')

    def _normalize_array(self, array):
        '''
        Normalizes an array so that sum of each row or column equals 1
        '''
        k = np.mean(np.sum(array, axis=1))
        return array/k

    def _jaccard(self, data, data2):
        return (data-data2)/(data+data2)

    def _chisq(self, data, data2):
        return (data-data2)**2/(data+data2)

    def normalize(self, data=None):
        '''
        Normalizes data so that sum of each row or column equals 1
        '''
        if data is None:
            try:
                self.data = self._normalize_array(self.data)
            except:
                print 'No data! Aborting normalization'
                return
        else:
            return self._normalize_array(data)
        try:
            self.data2 = self._normalize_array(self.data2)
        except:
            print 'No data2, continuing'

    def make_triangle(self, data=None):
        """
        Returns one half of a heatmap rotated by 45 degrees.
        """
        if data is None:
            data = self.data
        print 'Making triangle'
        data = np.tril(data)
        print 'Rotating'
        data = ndimage.rotate(data, 45, order=0, reshape=True, prefilter=False,
                              cval=0)
        data = self.hide_bads(data, mask_zeroes=True)
        data = data[~np.all(data == 0, axis=1)]
        return data

    def set_region_value(self, start, end, value=0, data=None):
        """
        Sets a specified **value** to all interactions of a region
        (**start**, **end**).
        Can be useful to hide regions with huge duplications or other
        abnormalities.
        """
        if data is None:
            data = self.data
        data[start:end, 0:data.shape[1]] = value
        data[0:data.shape[0], start:end] = value
        return data

    def remove_chromosomes(self, chrnames, data=None, chromosomes=None,
                           boundaries=None):
        """
        Removes all chromosomes with specifed **chrnames**. Now can only be used
        with the rightmost chromosomes, as it doesn't recalculate boundaries.
        You can do it yourself, though, it should work.
        """
        #TODO recalculate boundaries!
        if data is None:
            data = self.data
        if chromosomes is None:
            chromosomes = self.chromosomes
        if boundaries is None:
            boundaries = self.boundaries
        for chrname in chrnames:
            try:
                i = chromosomes.index(chrname)
            except ValueError:
                print 'Chromosome', chrname, 'not found!'
                continue
            start, end = boundaries[i]
            data = np.delete(data, range(start, end), 0)
            data = np.delete(data, range(start, end), 1)
            del boundaries[i]
            del chromosomes[i]
        return data

    def plot_whole_genome_heatmap(self, data=None, data2=None, triangle=False,
                                  diagonal_markers=False, compare=False,
                                  normalize=False, savepath=False,
                                  format='svg', colormap=False, log=True,
                                  colorbar=True, figsize=False,
                                  *args, **kwargs):
        if data is None:
            data = np.log2(self.data+1)
        else:
            if log:
                data = np.log2(data+1)
        if data2 is None:
            try:
                data2 = np.log2(self.data2+1)
            except AttributeError:
                pass
        if not colormap:
            colormap = self.cmap
        else:
            colormap = cmap.get_cmap(colormap)
        if not figsize:
            figsize = (15, 10)
        length, height = data.shape
        if triangle:
            data = self.make_triangle(data)
            try:
                data2 = self.make_triangle(data2)
            except:
                pass
        if normalize:
            data = self._normalize_array(data)
            if data2 is not None:
                data2 = self._normalize_array(data2)

        def determine_aspect(shape, extent):
            dx = (extent[1] - extent[0]) / float(shape[1])
            dy = (extent[3] - extent[2]) / float(shape[0])
            return dx / dy
        extent = [0, length*self.resolution,
                  0, length*self.resolution]
        aspect = determine_aspect(data.shape, extent)
        fig = plt.figure(figsize=figsize)
        fig.set_dpi(72)
        self.ax = host_subplot(111)
        self.ax.xaxis.set_tick_params(length=5, direction='out')
        self.ax.yaxis.set_tick_params(length=5, direction='out')

        self.ax.xaxis.set_tick_params(direction='out', length=5)
        self.ax.yaxis.set_tick_params(direction='out', length=5)

        self.ax.xaxis.set_major_formatter(mticker.FuncFormatter(
                                          self._nice_ticks))
        self.ax.yaxis.set_major_formatter(mticker.FuncFormatter(
                                          self._nice_ticks))
#        self.ax.minorticks_on()
#        self.ax.xaxis.set_minor_locator(MultipleLocator(200000))
#        self.ax.tick_params(axis='x', which='minor', bottom='on', top='off',
#                            direction='out', length=2, width=0.01)
#        self.ax.tick_params(axis='y', which='minor', left='on', right='off',
#                            direction='out', length=2, width=0.01)
        self.ax.tick_params(axis='x', which='major', bottom='on', top='on',
                            labelbottom='on', labelsize=12)
        self.ax.tick_params(axis='y', which='major', left='on', right='on',
                            labelleft='on', labelsize=12)
        self.ax.set(adjustable='box-forced')
        self.axChrLabels = self.ax.twin()
        self.locations = [sum(i)/2*self.resolution for i in self.boundaries]
        self.axChrLabels.set_xticks(self.locations)
        if triangle:
            self.axChrLabels.yaxis.set_visible(False)
            self.ax.yaxis.set_visible(True)
            self.axChrLabels.xaxis.tick_bottom()
        else:
            self.axChrLabels.yaxis.tick_right()
            self.axChrLabels.xaxis.tick_top()
        self.axChrLabels.set_xticklabels(self.chromosomes, fontsize=15)
        self.axChrLabels.set_yticks(self.locations)
        self.axChrLabels.set_yticklabels(self.chromosomes, fontsize=15)
        self.axChrLabels.xaxis.set_tick_params(length=0)
        self.axChrLabels.yaxis.set_tick_params(length=0)
        if diagonal_markers:
#            im = plt.imshow(data, visible=False, cmap=colormap, interpolation='none')
#            self.colorbar = plt.colorbar(im)
            norm = plt.Normalize(data[np.isfinite(data)].min(),
                                 data[np.isfinite(data)].max())
            data = colormap(norm(data))
            for multiple in diagonal_markers.keys():
                for start, end in self.boundaries:
                    for i in range(start, end, multiple/self.resolution):
                        data[i, i] = mcolors.ColorConverter().to_rgba(
                                                    diagonal_markers[multiple])
            im = ScalarMappable(norm, colormap)
            im.set_array(data)
            if colorbar:
                self.colorbar = plt.colorbar(im)
        plt.suptitle(self._make_title(), fontsize=15)
        if not compare:
            data = self.hide_bads(data, col='blue', alpha=1.0)
            self.image = plt.imshow(data, interpolation='none', origin='lower',
                         extent=extent, aspect=aspect, cmap=colormap,
                         *args, **kwargs)
            self.barlabel = u'$log_2(N\ of\ reads)$'
        elif compare == 'Jaccard':
            self.image = plt.imshow(self._jaccard(data, data2),
                                    interpolation='none', origin='lower',
                                    extent=extent, cmap='RdBu', *args, **kwargs)
            self.barlabel = u'Jaccard'
        elif compare == 'Chi-Square':
            self.image = plt.imshow(self._chisq(data, data2),
                                    interpolation='none', origin='lower',
                                    extent=extent, cmap='RdBu', *args, **kwargs)
            self.barlabel = u'Chi-squared'
        elif compare == 'Triangles':
            self._plot_combined_heatmap(data, data2, extent=extent,
                                        origin='lower', interpolation='none',
                                        *args, **kwargs)
            self._set_axlabels()
#        if compare != 'Triangles':
#            divider = make_axes_locatable(self.ax)
#            self.cax = divider.append_axes("right", size="5%", pad=0.05)
#            self.colorbar = plt.colorbar()
#            self.colorbar.set_label(self.barlabel, size=15)
        if not diagonal_markers:
            if colorbar:
                self.colorbar = plt.colorbar()
        self.colorbar.set_label(self.barlabel, size=15)

    def plot_chromosome_pair_heatmap(self, name, name2=None, data=None,
                                     compare=False, log=True, *args, **kwargs):
        n = self.chromosomes.index(name)
        start, end = self.boundaries[n]
        if name2 != None and name2 != name:
            n2 = self.chromosomes.index(name2)
            start2, end2 = self.boundaries[n2]
        else:
            name2, start2, end2 = name, start, end

        if data is None:
            if log:
                data = np.log2(self.data)
            else:
                data = self.data
            chrdata = data[start:end, start2:end2]
        elif type(data) is list or type(data) is tuple and compare:
            if len(data) == 2:
                data, data2 = data
        else:
            print 'Something wrong with data!'
            return ValueError
        if compare:
            try:
                data2 = self.data2
                if log:
                    data2 = np.log2(data2+1)
            except:
                print 'Can not compare, no data2 present!'
        self.ax = host_subplot(111)
        self.ax.xaxis.set_major_formatter(mticker.FuncFormatter(
                                          self._nice_ticks))
        self.ax.yaxis.set_major_formatter(mticker.FuncFormatter(
                                          self._nice_ticks))
        self.ax.xaxis.set_tick_params(top='off', direction='out', length=5)
        self.ax.yaxis.set_tick_params(right='off', direction='out', length=5)
        if not compare:
            plt.imshow(chrdata, origin='lower', interpolation='none',
                       extent=[0, (end-start)*self.resolution,
                               0, (end2-start2)*self.resolution],
                       *args, **kwargs)
            self.bar = plt.colorbar()
            self.bar.set_label(label=u'$log_2(N\ of\ reads)$', size=15)
        elif compare == 'Triangles':
            chrdata = data[start:end, start2:end2]
            chrdata2 = data2[start:end, start2:end2]
            self._plot_combined_heatmap(chrdata, chrdata2, origin='lower',
                                        interpolation='none',
                                        extent=[0, (end-start)*self.resolution,
                                            0, (end2-start2)*self.resolution],
                                        *args, **kwargs)
            self._set_axlabels()
        else:
            print 'This type of comparison is not implemented for chromosome '+\
                  'pair heatmaps'
            return
        plt.title(name+'-'+name2)

    def plot_by_chromosome_heatmaps(self, data=None, only_intrachromosome=True,
                                    exclude_names=('M'), compare=False,
                                    log=True, savepath=False, format='svg',
                                    *args, **kwargs):
        if data is None:
            data = self.data
        for chromosome in self.chromosomes:
            if chromosome not in exclude_names:
                if only_intrachromosome:
                    self.plot_chromosome_pair_heatmap(chromosome,
                                                      compare=compare,
                                                      log=log, *args, **kwargs)
                    if savepath:
                        filepath = path.join(path.abspath(savepath),
                                             chromosome+'-'+chromosome+'.'+
                                             format)
                        self.save_figure(filepath)
                        self.clear_figure()
                    else:
                        self.show_figure()
                else:
                    for chromosome1 in self.chromosomes:
                        if chromosome1 not in exclude_names:
                            self.plot_chromosome_pair_heatmap(chromosome,
                                                              chromosome1,
                                                              compare=compare,
                                                              log=log,
                                                              *args, **kwargs)
                            if savepath:
                                filepath = path.join(path.abspath(savepath),
                                                     chromosome+'-'+
                                                     chromosome1+'.'+
                                                     format)
                                self.save_figure(filepath, format)
                                self.clear_figure()
                            else:
                                self.show_figure()

    def _nonzero_stat_by_diagonal(self, data, func=np.mean):
        """
        Calculate some statistic (by default, mean) for each diagonal of a
        dataset. Return x-values and y-values.
        """
        x, y = [], []
        for d in xrange(len(data)):
            diagonal = np.diagonal(data, d)
            diagonal = diagonal[np.nonzero(diagonal)]
            x.append(self.resolution*d)
            y.append(func(diagonal))
        return np.array(x), np.array(y)

    def zero_interchromosomal_interactions(self, data=None):
        """
        Set all trans-interactions to 0
        """
        if data is not None:
            data = self.data
        for i in self.boundaries:
            data[i[0]:i[1], 0:i[0]] = 0
            data[i[0]:i[1], i[1]:-1] = 0
        return data

    def scale_plot(self, data=None, plot=True, stat=np.mean, *args, **kwargs):
        """
        Plot a scale plot. Doesn't seem to be working properly if compared to
        hiclib's plotScaling function, do not use this one!
        """
        if data is None:
            data = self.data
#        data[np.isneginf(data)]=0
        data = self.zero_interchromosomal_interactions(data)
#        self.ax=plt.subplot(111)
#        self.ax.imshow(data)
#        plt.show()
        x, y = self._nonzero_stat_by_diagonal(data, func=stat)
        if plot:
            try:
                self.ax.plot(x, y, *args, **kwargs)
            except AttributeError:
                self.ax = plt.subplot(111)
                self.ax.plot(x, y, *args, **kwargs)
            plt.xscale('log')
            plt.yscale('log')
        return x, y

    def fit(self, x, y, func='exp', minx='min', maxx='max', guess=None):
        from scipy.optimize import curve_fit
        def exp_func(x, a, b, c):
            return a * np.exp(-b * x) + c

        if func == 'exp':
            func = exp_func
            if guess is None:
                guess = (1, 1, 1)
        else:
            import inspect
            n = len(inspect.getargspec(func)['args'])-1
            guess = np.ones(n)

        if minx != 'min':
            x = x[x > minx]
            y = y[x > minx]
        if maxx != 'max':
            x = x[x < maxx]
            y = y[x < maxx]

        popt, pcov = curve_fit(func, x, y, guess)
        return x, exp_func(x, *popt), popt, pcov