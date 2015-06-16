# -*- coding: utf-8 -*-
"""
This is a class for working with genomic intervals (bed- and bedgraph files) and
searching for TADs. For some operations uses pybedtools. Uses greendale for TADs
mapping.
"""
from __future__ import division
import pandas as pd
import numpy as np
import pybedtools as pbt

def _list_files(path):
    '''
    List all files in a directory
    '''
    import os
    files = []
    for name in os.listdir(path):
        if os.path.isfile(os.path.join(path, name)):
            files.append(os.path.join(path, name))
    return files

class GenomicIntervals(object):
    def __init__(self, settings=None):
        if settings is not None:
            self.settings = settings
            self.settings.extract_settings(self)
        
    def read_intervals(self, interval_file, bedtool=False,
                       chrom=True, names=False, header=True):
        '''
        Read an interval file, bed-style, but not necessarily with chromosome
        data (e.g. in genomic coordinates if you concatenate chromosomes). If
        chrom=True, then first column is read as chromosomes. Returns a pandas
        DataFrame object with 'Start' and 'End' columns (and optionally
        'Chromosome' as the first one). If bedtool=True, then return a BedTool.
        '''
        if bedtool:
            intervals = pbt.BedTool(interval_file)
        else:
            if names:
                columns = names
            elif chrom:
                columns = 'Chromosome', 'Start', 'End'
            else:
                columns = 'Start', 'End'
            intervals = pd.read_csv(interval_file, sep='\t', names=columns,
                                    header=header)
        return intervals
    
    def read_bedgraph(self, bedgraph_file, bedtool=True, skiprows=1):
        '''
        Read a bedgraph file, returns a BedTool or an interval-style pandas
        DataFrame with an additional column 'Score'. For UCSC tracks use
        skiprows=1 to omit the header.
        '''
        if not skiprows:
            if bedtool:
                return pbt.BedTool(bedgraph_file)
            else:
                return pd.read_csv(bedgraph_file,
                                  names=['Chromosome', 'Start', 'End', 'Score'],
                                  skiprows=skiprows, delim_whitespace=True)
        else:
            data = pd.read_csv(bedgraph_file,
                               names=['Chromosome', 'Start', 'End', 'Score'],
                               skiprows=skiprows, delim_whitespace=True)
            if not bedtool:
                return data
            else:
                return self.bedtool_from_df(data)
    
    def make_bins(self, binsize=1000, bedtool=True):
        '''
        Split all chromosomes into bins by binsize. Returns a BedTool or an
        interval-style pandas DataFrame.
        '''
        if 'lengths_bp' not in dir(self) or 'chromosomes' not in dir(self):
            raise AttributeError('Please, set chromosomes with lengths to make '
                                 'bins')
        a = pbt.BedTool()
        a = a.window_maker(w=binsize, g=self.settings.boundaries_bp_dict)
        if not bedtool:
            return pd.read_table(a.fn, names=['Chromosome', 'Start', 'End'])
        else:
            return a
#        bins = pd.DataFrame(columns=['Chromosome', 'Start', 'End'])
#        for name, length in zip(self.chromosomes, self.lengths_bp):
#            bin_boundaries = range(0, length, binsize)
#            if length > bin_boundaries[-1]:
#                bin_boundaries.append(length)
#            bins_chr = pd.DataFrame({'Chromosome':name,
#                                     'Start':bin_boundaries[:-1],
#                                     'End':bin_boundaries[1:]})
#            bins = bins.append(bins_chr, ignore_index=True)
#        return bins[['Chromosome', 'Start', 'End']]
    
    def bedtool_from_df(self, df):
        '''
        Make a BedTool from a pandas DataFrame.
        '''
        return pbt.BedTool(df.to_string(header=False, index=False),
                           from_string=True).sort()
    
    def binarize_bedgraph(self, bedgraph_bedtool, binsize, function='mean',
                          bedtool=True, *args, **kwargs):
        '''
        Map bedgraph to bins of specified size. Apply any function available for
        `bedtools map` function. Return a BedTool or interval-style pandas
        DataFrame with an additional column 'Score'. All arguments are passed to
        the `BedTool.map` function.
        NOTICE: if an interval crosses the bin border, it will be counted in
        both bins by default; bedtools doesn't support splitting score in
        mapping AFAIK.
        '''
        bins_bedtool = self.make_bins(binsize)
        mapped = bins_bedtool.map(bedgraph_bedtool, c=4, o=function, null=0,
                                  *args, **kwargs)
        if not bedtool:
            return pd.read_table(mapped.fn, names=['Chromosome', 'Start', 'End',
                                                   'Score'])
        else:
            return mapped

    def genome_coordinate_to_chr(self, coordinate, from_bins=True):
        '''
        Recalculate a coordinate from genome-based to chromosome-based.
        '''
        if from_bins:
            boundaries_bp, starts_bp, ends_bp = (
                                self.settings.calculate_approx_boundaries_bp())
        else:
            boundaries_bp = self.boundaries_bp
            starts_bp, ends_bp = self.settings.starts_bp, self.settings.ends_bp
        chrcoordinate = None
        for i, boundaries in enumerate(boundaries_bp):
            start, end = boundaries
            chrname = self.chromosomes[i]
            if start<=coordinate<end:
                chrcoordinate = coordinate - starts_bp[i]
                break
            elif (i == len(boundaries_bp)-1 and
                            0 < coordinate - ends_bp[i] < self.resolution):
                if self.settings.boundaries is None:
                    self.boundaries = self.settings.calculate_boundaries()
                chrcoordinate = self.boundaries[i][1]
                break
            else:
                continue
        if chrcoordinate is None:
            raise ValueError('Coordinate '+ str(coordinate) +' out of range')
        return chrname, chrcoordinate

    def chr_interval_to_genome(self, chromosome, start, end,
                                 from_bins=True):
        '''
        Recalculate a coordinate from chromosome-based to genome-based.
        '''
        if from_bins:
            boundaries_bp, starts_bp, ends_bp = (
                                self.settings.calculate_approx_boundaries_bp())
        else:
            starts_bp, ends_bp = self.settings.starts_bp, self.settings.ends_bp
        chrn = self.chromosomes.index(chromosome)
        diff = starts_bp[chrn]
        return start+diff, end+diff
    
    def _chr_interval_to_genome(self, coordinate, from_bins=True):
        start, end = self.chr_interval_to_genome(coordinate['Chromosome'],
                                                 coordinate['Start'],
                                                 coordinate['End'],
                                                 from_bins)
        return pd.Series({'Start_gen':start, 'End_gen':end})
      
    def chr_intervals_to_genome(self, intervals, from_bins=True): 
        from functools import partial
        f = partial(self._chr_interval_to_genome, from_bins=from_bins)
        return intervals.apply(f, axis=1)[['Start_gen', 'End_gen']]\
                                                              / self.resolution

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
            intervals = intervals[~unmatched]
        intervals = intervals.drop('End_chromosome', axis=1)
        intervals = intervals.rename(columns={'Start_chromosome':'Chromosome'})
        return intervals

    def genome_intervals_to_chr(self, intervals, remove_crossborder=True):
        '''
        Recalculate coordinates of intervals from genome-based to
        chromosome-based. Intervals is a pandas DataFrame with at least 2
        columns, 'Start' and 'End', with coordinates in a concatenated genome.
        '''
        new_intervals = pd.DataFrame(columns=['Start_chromosome', 'Start',
                                              'End_chromosome', 'End'],
                                            index=intervals.index)
        for i in intervals.index:
            startchr, start = self.genome_coordinate_to_chr(intervals.ix[i,
                                                                    'Start'])
            try:
                endchr, end = self.genome_coordinate_to_chr(intervals.ix[i, 'End'])
            except ValueError:
                if i == intervals.index[-1]:
                    endchr = self.chromosomes[-1]
                    end = self.boundaries_bp[-1]
                else:
                    raise ValueError('Not last end is out of boundaries')
            new_intervals.iloc[i]=startchr, start, endchr, end
        cols = [col for col in intervals.columns if col not in ['Start', 'End']]
        new_intervals[cols] = intervals[cols]
        if remove_crossborder:
            new_intervals = self._remove_interchr_intervals(new_intervals)
        return new_intervals
        
    def make_inter_intervals(self, intervals, shorten_by_resolution=False):
        '''
        Accepts a bed-style pandas DataFrame with columns 'Chromosome', 'Start',
        'End'. Returns the file with the same format, but containing coordinates
        of intervals between provided (doesn't include telomeric regions not
        covered with intervals!). If shorten_by_resolution, subtracts
        self.resolution from ends (useful for intervals, acquired from Hi-C
        data, such as TADs).
        '''
        intervals = intervals.sort(columns=['Chromosome', 'Start'])
        start = intervals[:-1][['End', 'Chromosome']]
        start = start.reset_index(drop=True).rename(columns={
                                              'Chromosome':'Start_chromosome',
                                              'End':'Start'})

        end = intervals[1:][['Start', 'Chromosome']]
        if shorten_by_resolution:
            end['Start'] -= self.resolution
        end = end.reset_index(drop=True).rename(columns={
                                                'Chromosome':'End_chromosome',
                                                'Start':'End'})
        inter_intervals = start.merge(end, left_index=True, right_index=True)
        inter_intervals = self._remove_interchr_intervals(inter_intervals)
        return inter_intervals[['Chromosome', 'Start','End']]

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
        from greendale import segment
        Wcomm, Wnull, pass_mask, length = parameters
        starts, scores = segment.potts_segmentation(Wcomm, Wnull, gamma,
                                                    pass_mask=pass_mask)
        pos = np.r_[starts, length]
        domains = zip(pos[:-1], pos[1:])
        domains = np.array([(i[0], i[1]) for i in domains if i[1]-i[0]>2])
        domains = pd.DataFrame(domains, columns=('Start', 'End'))
        return domains

    def find_TADs(self, data, gammalist=range(10, 110, 10)):
        '''
        Finds TADs in data with a list of gammas. Returns a pandas DataFrame
        with columns 'Start', 'End' and 'Gamma'. Use genome_intervals_to_chr on
        the returned object to get coordinates in bed-style format and not in
        coordinates of concatenated genome.
        '''
        parameters  = self._precalculate_TADs_in_array(data)
        domains = pd.DataFrame(columns=('Start', 'End', 'Gamma'))
        for g in gammalist:
            domains_g = self._calculate_TADs(parameters, g)
            domains_g['Gamma'] = g
            domains_g['Start'] *= self.resolution
#            domains_g['End'] -= 1
            domains_g['End'] *= self.resolution
            domains = domains.append(domains_g, ignore_index=True)
        domains[['Start', 'End']] = domains[['Start', 'End']].astype(int)
        domains = domains[['Start', 'End', 'Gamma']]
        return domains

    def find_TADs_by_chromosomes(self, data, gammadict={}):
        '''
        Åpply TAD finding to each chromosome separately. As required gamma
        varies very much with size of supplied matrix for calculation, you
        should supply different gammas for each chromosome in a
        gammadict{chromosome_name:[gamma1, gamma2, ...],
                  ...}
        Returns a pandas DataFrame with columns 'Chromosome', 'Start', 'End'
        and 'Gamma'
        '''
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
#                domains_g['End'] -= 1
                domains_g['End'] *= self.resolution
                domains_chr = domains_chr.append(domains_g, ignore_index=True)
            domains = domains.append(domains_chr, ignore_index=True)
        domains[['Start', 'End']] = domains[['Start', 'End']].astype(int)
        domains = domains[['Chromosome', 'Start', 'End', 'Gamma']]
        return domains.reset_index(drop=True)

    def write_TADs(self, domains, basename):
        for g in set(domains['Gamma']):
            domains_g = domains[domains['Gamma']==g]
            domains_g.to_csv(basename+'_g'+str(g)+'.bed', header=False,
                           index=False, sep='\t', cols=['Chromosome', 'Start',
                                                        'End'])

    def read_TADs(self, path, basename='TADs',
                  names=['Chromosome', 'Start', 'End'], listfiles=False):
        from os import path as p
        files = _list_files(path)
        files = [f for f in files if p.basename(f).startswith(basename)]
        if listfiles:
            print 'Reading TADs from these files:'
            print files
        domains = pd.DataFrame(columns=names)
        for f in files:
            g = float(f.split('.')[0].split('_g')[1])
            domains_g = self.read_intervals(f, names=names)
            domains_g['Gamma'] = g
            domains = domains.append(domains_g, ignore_index=True)
        return domains[names+['Gamma']]

    def make_interTADs(self, domains):
        '''
        Makes inter-TADs from TADs DataFrame (columns 'Chromosome', 'Start',
        'End', 'Gamma').
        '''
        interTADs = pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Gamma'])
        for g in set(domains['Gamma']):
            interTADs_g = self.make_inter_intervals(domains[domains['Gamma']==g])
            interTADs_g['Gamma'] = g
            interTADs = interTADs.append(interTADs_g, ignore_index=True)
        return interTADs

    def describe_TADs(self, domains, functions=None):
        '''
        Group TADs gy 'Gamma' and returns length statistics by group. Includes
        count, np.median, np.mean, np.min, np.max and coverage by default.
        Coverage calculates genome coverage of TADs based on true chromosome
        lengths. If supplied with functions argument, you can add any other
        functions to that statistics. Functions should take 1 argument and they
        are applied to a pdDataFrame column with lengths of TADs.
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
        lengths = pd.DataFrame(columns=['Gamma', 'Length'])
        lengths['Gamma']=domains['Gamma']
        lengths['Length'] = domains['End']-domains['Start']
        lengths['Length'] = lengths['Length'].astype(int)
        stats = lengths.groupby('Gamma')['Length'].agg([count, np.median,
                                                        np.mean, np.min,
                                                        np.max, coverage] +
                                                       list(functions))
        return stats

    def plot_TADs_length_distribution(self, domains, show=True, *args, **kwargs):
        '''
        Plots size distribution of TADs using DataFrame.hist() method by 'Gamma'
        . Creates a column for each possible size in the data. All arguments are
        passed to the .hist() method. Returns a list of lists of axes.
        '''
        import matplotlib.pyplot as plt
        if 'Length' not in domains.columns:
            domains['Length'] = domains['End']-domains['Start']
            domains['Length'] = domains['Length'].astype(int)
        bins = range(int(-self.resolution/2),
                     int(max(domains['Length'])+3*self.resolution/2),
                     int(self.resolution))
        axes = domains['Length'].hist(by=domains['Gamma'], bins=bins,
                                                                *args, **kwargs)
        for ax in axes.flatten():
            ax.set_xlim([-self.resolution/2,
                         max(domains['Length'])+self.resolution/2])
        if show:
            plt.show()
        return axes
    
    def compare_intervals(self, intervals1, intervals2, precision_both=None,
                          precision_each=None):
        '''
        Accepts two pandas datasets with intervals (chromosome-based) and
        returns identical, present in only the first and only the second
        dataset. You can specify precision_both or precision_each for
        approximate comparison. The first relates to the sum of differences in
        both coordinates (e.g. precision_both=40000 will see all of these
        examples as identical, as well as perfectly matching intervals:
        (('chr1', 0, 200000), ('chr1', 20000, 200000)) -> sum(differences) = 20000
        (('chr1', 0, 200000), ('chr1', 0, 240000)) -> sum(differences) = 40000
        (('chr1', 20000, 200000), ('chr1', 0, 180000)) -> sum(differences) = 40000
        (('chr1', 20000, 200000), ('chr1', 0, 180000)) -> sum(differences) = 40000
        etc.
        The second relates to each difference, so each of them mustn't exceed
        the precision.
        Default for both - 0 (exact comparison). If one is specified, the other
        one is not used; both can't be used simultaneously.
        '''
        if precision_both is not None and precision_each is not None:
            raise ValueError, 'Please specify only 1 type of precision'
        elif precision_both is not None:
            mode = 'both'
        elif precision_each is not None:
            mode = 'each'
        elif precision_both is None and precision_each is None:
            mode = 'exact'
        else:
            raise ValueError, 'Something wrong with precision specification'
        
        def remove_identical(ds1_list, ds2_list, left=True, right=True):
            '''
            Idea comes from here:
            http://stackoverflow.com/a/29464365/1304161
            '''
            ds1 = set(ds1_list)
            ds2 = set(ds2_list)
            if left:
                intervals1_unique = pd.DataFrame(list(ds1.difference(ds2)),
                                                 columns=intervals1.columns)
                intervals1_unique.sort(columns = ['Chromosome', 'Start', 'End'],
                                   inplace=True)
                if not right:
                    return intervals1_unique
            if right:
                intervals2_unique = pd.DataFrame(list(ds2.difference(ds1)),
                                                 columns=intervals2.columns)
                intervals2_unique.sort(columns = ['Chromosome', 'Start', 'End'],
                                   inplace=True)
                if not left:
                    return intervals2_unique
            return intervals1_unique, intervals2_unique
            
        ds1 = list(tuple(line) for line in intervals1.values)
        ds2 = list(tuple(line) for line in intervals2.values)
        
        if mode == 'exact':
            shared = pd.merge(intervals1, intervals2,
                              on=['Chromosome', 'Start', 'End'], how='inner')
            shared.sort(columns = ['Chromosome', 'Start', 'End'], inplace=True)
            intervals1_unique, intervals2_unique = remove_identical(ds1, ds2)
            return shared, intervals1_unique, intervals2_unique
            
        elif mode == 'both':
            
            def is_similar(coord1, coord2):
                chrom1, start1, end1 = coord1
                chrom2, start2, end2 = coord2
                return chrom1==chrom2 and abs(start1-start2) + abs(end1-end2)\
                                                              <= precision_both
                            
        elif mode == 'each':
            
            def is_similar(coord1, coord2):
                chrom1, start1, end1 = coord1
                chrom2, start2, end2 = coord2
                return chrom1==chrom2 and all(np.abs(np.array(
                                 (start1-start2, end1-end2))) <= precision_each)
                
        shared1 = []
        shared2 = []
        for coord1 in ds1:
            for coord2 in ds2:
                if is_similar(coord1, coord2):
                    shared1.append(coord1)
                    shared2.append(coord2)
                    break
                
        shared1 = pd.DataFrame(shared1, columns=intervals1.columns)
        shared2 = pd.DataFrame(shared2, columns=intervals2.columns)
        shared1_list = [tuple(line) for line in shared1.values]
        shared2_list = [tuple(line) for line in shared2.values]
        intervals1_unique = remove_identical(ds1, shared1_list, right=False)
        intervals2_unique = remove_identical(ds2, shared2_list, right=False)
        return shared1, shared2, intervals1_unique, intervals2_unique
    
    def get_density(self, interval, data):
        '''
        Get sum of all interactions (i.e. density) in a square from interval.
        '''
        start, end = tuple(interval)
        if any(np.isnan([start, end])):
               return np.NaN
        start, end = int(start), int(end)
        return 1.0*np.sum(data[start:end, start:end]) 

    def get_density_frac(self, interval, data):
        start, end = tuple(interval)
        if any(np.isnan([start, end])):
               return np.NaN
        start, end = int(start), int(end)
        return 1.0*np.sum(data[start:end, start:end])/np.sum(data[start:end])/2        
            
    def get_intervals_density(self, intervals, data, triangle=True, norm=False,
                              frac=False):
        '''
        Calculates "density" of all intervals from heatmap data, as measured by
        number of ligation junctions inside the intervals. By default only
        uses half of the square matrix. Optionally normalizes as reads per
        million (*norm*). Optionally normalizes sum in an interval by the sum
        of all columns of the interval (*frac*).
        '''
        from functools import partial
        
        if triangle:
            data = np.triu(data)
        if norm:
            data /= np.sum(data)
            data *= 10**6
        ints = self.chr_intervals_to_genome(intervals)
        if frac:
            f = partial(self.get_density_frac, data=data)
        else:
            f = partial(self.get_density, data=data)
        ints = ints.apply(f, axis=1)
        return ints

    def get_border_strength(self, coordinates, data):
        start1, end1, start2, end2 = np.array(coordinates).astype(int)
        sum1 = np.sum(data[start1:end1, start1:end1])
        sum2 = np.sum(data[start2:end2, start2:end2])
        sum_inter = np.sum(data[start1:end1, start2:end2])
        return (sum1 + sum2)/sum_inter
    
    def get_borders_strength(self, intervals, data):
        '''
        Calculate strength of domain borders. Divides sum of interactions
        within two neighbouring domains by the sum of interactions between
        those domains.
        '''
        from functools import partial
        ints = self.chr_intervals_to_genome(intervals)
        ints1 = ints[1:].reset_index(drop=True)
        ints1.columns = ['Start1', 'End1']
        ints2 = ints[:-1].reset_index(drop=True)
        ints2.columns = ['Start2', 'End2']
        ints = pd.merge(ints1, ints2, left_index=True, right_index=True)
        ends = self.settings.ends
        crosschr = []
        for i in ints.index:
            if np.searchsorted(ends, ints['Start1'][i]) !=\
               np.searchsorted(ends, ints['Start2'][i]):
                crosschr.append(i)
        ints.drop(crosschr, axis=0, inplace=True)
        f = partial(self.get_border_strength, data=data)
        borders = self.make_inter_intervals(intervals)
        borders['Strength'] = ints.apply(f, axis=1)
        return pd.merge(borders, ints, left_index=True, right_index=True)
        