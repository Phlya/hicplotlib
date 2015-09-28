# -*- coding: utf-8 -*-

from __future__ import division
import pandas as pd
import numpy as np
import pybedtools as pbt
from functools import partial
try:
    from greendale import segment
except ImportError:
    print 'Please install greendale (https://bitbucket.org/nvictus/greendale).'
    print 'It is used for TAD calling.'

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

def _precalculate_TADs_in_array(array):
    '''
    Calculates greendale statistics for calling TADs. These may be
    reused with multiple gammas very fast.
    '''
    k = array.sum(axis=0)
    pass_mask = k != 0
    Wcomm = segment.normalized_weights_by_segment(array)
    Wnull = segment.normalized_weights_by_segment(np.outer(k, k))
    return Wcomm, Wnull, pass_mask, len(array)

def _calculate_TADs(Wcomm, Wnull, pass_mask, length, gamma,
                    segmentation='potts', write_g=True):
    '''
    Calculate TADs based on greendale statistics (*parameters*) and a
    specified gamma value. Returns a pandas DataFrame with columns 'Start'
    and 'End' with coordinates in bins.
    '''
    if segmentation == 'potts':
        starts, scores = segment.potts_segmentation(Wcomm, Wnull, gamma,
                                                    pass_mask=pass_mask)
    elif segmentation == 'armatus':
        starts, scores = segment.armatus_segmentation(Wcomm, gamma,
                                                      pass_mask=pass_mask)
    else:
        raise ValueError, 'Unsupported segmentation, use potts or armatus'
    pos = np.r_[starts, length]
    domains = zip(pos[:-1], pos[1:], scores)
    domains = pd.DataFrame(domains, columns=('Start', 'End', 'Score'))
    if write_g:
        domains['Gamma'] = gamma
    return domains.reset_index(drop=True)

class GenomicIntervals(object):
    """
    This is a class for working with genomic intervals (bed- and bedgraph files)
    and searching for TADs. For some operations uses pybedtools. Uses greendale for TADs
    mapping.
    """
    def __init__(self, settings=None):
        if settings is not None:
            self.settings = settings
            self.settings.extract_settings(self)

    def read_intervals(self, interval_file, bedtool=False,
                       chrom=True, names=False, header=True, sort=True):
        '''
        Read an interval file, bed-style, but not necessarily with chromosome
        data (e.g. in genomic coordinates if you concatenate chromosomes). If
        chrom=True, then first column is read as chromosomes. Returns a pandas
        DataFrame object with 'Start' and 'End' columns (and optionally
        'Chromosome' as the first one). If bedtool=True, then return a BedTool. By default
        sorts intervals by all specified columns (or, if not specified, by ['Chromosome',
        'Start', 'End']) (or just ['Start', 'End'] if chrom=False).
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
            if sort:
                intervals = intervals.sort(list(columns))
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
            raise AttributeError('Please, set chromosomes with lengths to make'
                                 ' bins')
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
        Map bedgraph to bins of specified size. Apply any function available
        for *bedtools map* function. Return a BedTool or interval-style pandas
        DataFrame with an additional column 'Score'. All arguments are passed
        to the *BedTool.map* function.
        NOTICE: if an interval crosses the bin border, it will be counted in
        both bins by default; bedtools doesn't support splitting score in
        mapping AFAIK.
        '''
        bins_bedtool = self.make_bins(binsize)
        mapped = bins_bedtool.map(bedgraph_bedtool, c=4, o=function, null=0,
                                  *args, **kwargs)
        if not bedtool:
            return pd.read_table(mapped.fn, names=['Chromosome', 'Start',
                                                   'End', 'Score'])
        else:
            return mapped

    def genome_coordinate_to_chr(self, coordinate, from_bins=True):
        '''
        Recalculate a coordinate from genome-based to chromosome-based.
        '''
        if from_bins:
            boundaries_bp, starts_bp, ends_bp = self.settings.calculate_approx_boundaries_bp()
        else:
            boundaries_bp = self.boundaries_bp
            starts_bp, ends_bp = self.settings.starts_bp, self.settings.ends_bp
        chrcoordinate = None
        for i, boundaries in enumerate(boundaries_bp):
            start, end = boundaries
            chrname = self.chromosomes[i]
            if start <= coordinate < end:
                chrcoordinate = coordinate - starts_bp[i]
                break
            elif i == len(boundaries_bp)-1 and 0 < coordinate - ends_bp[i] < self.resolution:
                if self.settings.boundaries is None:
                    self.boundaries = self.settings.calculate_boundaries()
                chrcoordinate = self.boundaries[i][1]
                break
            else:
                continue
        if chrcoordinate is None:
            raise ValueError('Coordinate '+ str(coordinate) +' out of range')
        return chrname, chrcoordinate

    def chr_interval_to_genome(self, chromosome, start, end, from_bins=True):
        '''
        Recalculate a coordinate from chromosome-based to genome-based.
        '''
        if from_bins:
            boundaries_bp, starts_bp, _ = self.settings.calculate_approx_boundaries_bp()
        else:
            starts_bp = self.settings.starts_bp
        chrn = self.chromosomes.index(chromosome)
        diff = starts_bp[chrn]
        return start+diff, end+diff

    def _chr_interval_to_genome(self, coordinate, from_bins=True,
                                in_cols=['Chromosome', 'Start', 'End'],
                                out_cols=['Start', 'End']):
        start, end = self.chr_interval_to_genome(coordinate[in_cols[0]],
                                                 coordinate[in_cols[1]],
                                                 coordinate[in_cols[2]],
                                                 from_bins)
        return pd.Series({out_cols[0]:start, out_cols[1]:end})

    def chr_intervals_to_genome(self, intervals, from_bins=True,
                                in_cols=['Chromosome', 'Start', 'End'],
                                out_cols=['Start', 'End']):
        f = partial(self._chr_interval_to_genome, from_bins=from_bins, in_cols=in_cols,
                    out_cols=out_cols)
        return intervals.apply(f, axis=1)[out_cols] / self.resolution

    def _remove_interchr_intervals(self, intervals,
                                   in_chr_names=['Start_chromosome', 'End_chromosome'],
                                   out_chr_name='Chromosome'):
        '''
        Accepts a pandas DataFrame with at least 2 columns 'Start_chromosome' and 
        'End_chromosome' (or other, depending on **in_chr_names**), checks all rows to
        contain the same value in these columns and removes the ones, where is is not
        true. Returns the object, identical to supplied intervals, but with only one
        column 'Chromosome' (or another, depending on **out_chr_name** instead  of two
        chromosome columns.
        '''
        unmatched = intervals[in_chr_names[0]] != intervals[in_chr_names[1]]
        if any(unmatched):
            intervals = intervals[~unmatched]
        intervals = intervals.drop(in_chr_names[1], axis=1)
        intervals = intervals.rename(columns={in_chr_names[0]:out_chr_name})
        return intervals

    def genome_intervals_to_chr(self, intervals, remove_crosschr=True,
                                in_cols=['Start', 'End'],
                                out_cols=['Start_chromosome', 'Start',
                                          'End_chromosome', 'End'],
                                single_chr_col='Chromosome'):
        '''
        Recalculate coordinates of intervals from genome-based to
        chromosome-based. Intervals is a pandas DataFrame with at least 2
        columns, 'Start' and 'End', with coordinates in a concatenated genome.
        '''
        new_intervals = pd.DataFrame(columns=out_cols,
                                     index=intervals.index)
        for i in intervals.index:
            startchr, start = self.genome_coordinate_to_chr(intervals.ix[i, in_cols[0]])
            try:
                endchr, end = self.genome_coordinate_to_chr(intervals.ix[i, in_cols[1]])
            except ValueError:
                if i == intervals.index[-1]:
                    endchr = self.chromosomes[-1]
                    end = self.boundaries_bp[-1][-1]
                else:
                    raise ValueError('Not last end is out of boundaries')
            new_intervals.iloc[i] = [startchr, start, endchr, end]
        new_intervals = new_intervals.reset_index(drop=True)
        cols = [c for c in intervals.columns if c not in in_cols]
        new_intervals[cols] = intervals[cols].reset_index(drop=True)
        if remove_crosschr:
            new_intervals = self._remove_interchr_intervals(new_intervals,
                                                            out_chr_name=single_chr_col)
        return new_intervals

    def make_inter_intervals(self, intervals, shorten_by_resolution=False):
        '''
        Accepts a bed-style pandas DataFrame with columns 'Chromosome',
        'Start', 'End'. Returns the file with the same format, but containing
        coordinates of intervals between provided (doesn't include telomeric
        regions not covered with intervals!). If shorten_by_resolution,
        subtracts self.resolution from ends (useful for intervals, acquired
        from Hi-C data, such as TADs).
        '''
        intervals = intervals.sort(columns=['Chromosome', 'Start'])
        start = intervals[:-1][['End', 'Chromosome']]
        start = start.reset_index(drop=True).rename(columns={'Chromosome':'Start_chromosome',
                                                             'End':'Start'})

        end = intervals[1:][['Start', 'Chromosome']]
        if shorten_by_resolution:
            end['Start'] -= self.resolution
        end = end.reset_index(drop=True).rename(columns={'Chromosome':'End_chromosome',
                                                         'Start':'End'})
        inter_intervals = start.merge(end, left_index=True, right_index=True)
        inter_intervals = self._remove_interchr_intervals(inter_intervals)
        return inter_intervals[['Chromosome', 'Start', 'End']]

    def find_TADs(self, data, gammalist=range(10, 110, 10), segmentation='potts',
                  minlen=3, drop_gamma=False, n_jobs='auto'):
        '''
        Finds TADs in data with a list of gammas. Returns a pandas DataFrame
        with columns 'Start', 'End' and 'Gamma'. Use genome_intervals_to_chr on
        the returned object to get coordinates in bed-style format and not in
        coordinates of concatenated genome.
        If *drop_gamma*, drops the 'Gamma' column (useful when using 1 gamma)
        '''
        if n_jobs is 'auto': #Empirical values on my computer; with >8 Gb memory try increasing n_jobs
            if segmentation == 'potts':
                n_jobs = 3
            elif segmentation == 'armatus':
                n_jobs = 6
        if ~np.isfinite(data).any():
            print 'Non-finite values in data, substituting them with zeroes'
            data[~np.isfinite(data)] = 0
        Wcomm, Wnull, pass_mask, length = _precalculate_TADs_in_array(data)
        f = _calculate_TADs
        if n_jobs >= 1:
            from joblib import Parallel, delayed
            domains = Parallel(n_jobs=n_jobs, max_nbytes=1e6)(
                              delayed(f)(Wcomm, Wnull, pass_mask, length, g, segmentation)
                                                                       for g in gammalist)
        elif n_jobs is None or n_jobs == False or n_jobs == 0:
            domains = []
            for g in gammalist:
                domains_g = f(Wcomm, Wnull, pass_mask, length, g, segmentation)
                domains.append(domains_g)
        domains = pd.concat(domains, ignore_index=True)
        domains = domains.query('End-Start>='+str(minlen)).copy()
        domains = domains.sort(columns=['Gamma', 'Start', 'End'])
        domains.reset_index(drop=True, inplace=True)
        domains[['Start', 'End']] = domains[['Start', 'End']].astype(int)
        domains[['Start', 'End']] *= self.resolution
        domains = domains[['Start', 'End', 'Score', 'Gamma']]
        if drop_gamma:
            domains.drop('Gamma', axis=1, inplace=True)
        domains = self.genome_intervals_to_chr(domains).reset_index(drop=True)
        return domains

    def find_TADs_by_chromosomes(self, data, gammadict={}, minlen=3):
        '''
        Åpply TAD finding to each chromosome separately. As required gamma
        varies very much with size of supplied matrix for calculation, you
        should supply different gammas for each chromosome in a
        *gammadict={chromosome_name:[gamma1, gamma2, ...], ...}*.
        Returns a pandas DataFrame with columns 'Chromosome', 'Start', 'End'
        and 'Gamma'.
        '''
        raise DeprecationWarning, 'Will be deprecated unless a suitable use-case is found'
        domains = []
        for i, chrname in enumerate(self.chromosomes):
            start, end = self.boundaries[i]
            chrdata = data[start:end, start:end]
            parameters = _precalculate_TADs_in_array(chrdata)
            domains_chr = []
            for g in gammadict[chrname]:
                domains_g = _calculate_TADs(*parameters, gamma=g)
                domains_g = domains_g.query('End-Start>='+str(minlen)).copy()
                domains_g.reset_index(drop=True, inplace=True)
                domains_g['Chromosome'] = chrname
                domains_g['Gamma'] = g
                domains_g['Start'] *= self.resolution
                domains_g['End'] *= self.resolution
            domains = pd.concat(domains, ignore_index=True)
        domains = pd.concat(domains_chr, ignore_index=True)
        domains[['Start', 'End']] = domains[['Start', 'End']].astype(int)
        domains = domains[['Chromosome', 'Start', 'End', 'Score', 'Gamma']]
        return domains

    def write_TADs(self, domains, path, *args, **kwargs):
        '''
        Save TADs in files, one gamma in one file. *Path* specifies the
        location and, if required, beginning of the filenames. All other args
        and kwargs are passed on to DF.to_csv method.
        '''
        for g in set(domains['Gamma']):
            domains_g = domains[domains['Gamma'] == g]
            domains_g.to_csv(path+'_g'+str(g)+'.bed', *args, **kwargs)

    def read_TADs(self, path, basename='TADs',
                  names=['Chromosome', 'Start', 'End'], listfiles=False):
        '''
        Reads TADs from files written in the style of write_TADs function.
        '''
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
        interTADs = []
        for g, data in domains.groupby('Gamma'):
            interTADs_g = self.make_inter_intervals(data.copy())
            interTADs_g['Gamma'] = g
        interTADs = pd.concat(interTADs, ignore_index=True)
        return interTADs

    def describe_TADs(self, domains, functions=None, feature='Length'):
        '''
        Groups TADs gy 'Gamma' and returns statistics by group. Includes
        count, np.median, np.mean, np.min, np.max and coverage by default.
        Coverage calculates genome coverage of TADs based on true chromosome
        lengths. If supplied with functions argument, you can add any other
        functions to that statistics. Functions should take 1 argument and they
        are applied to a pd.DataFrame column with the feature of TADs.
        If *feature* supplied, doesn't use length but rather that specified
        column of the dataset (and thus doesn't calculate coverage as it
        probably doesn't make sense).
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

        data = pd.DataFrame(columns = ['Gamma', feature])
        data['Gamma'] = domains['Gamma']
        funcs = [count, np.median, np.mean, np.min, np.max]
        if feature == 'Length' and 'Length' not in domains.columns:
            data[feature] = domains['End']-domains['Start']
            funcs.append(coverage)
        else:
            data[feature] = domains[feature]

        data[feature] = data[feature].astype(float)
        stats = data.groupby('Gamma')[feature].agg(funcs + functions)
        return stats

    def plot_TADs_distribution(self, domains, feature='Length', group='Gamma',
                               kind='hist', bins='autolength', autoxlim=True,
                               *args, **kwargs):
        '''
        By default plots size distribution of TADs using DataFrame.plot()
        method by 'Gamma'. Creates a column for each possible size in the
        data. All arguments are passed to the .plot() method.
        Alternatevely, any other column can be used as grouping, and any other
        column can be used for calculation. If *kind* is other than 'hist',
        sns.factorplot with that *kind* is used.
        Returns an sns.FacetGrid if.
        '''
        from matplotlib import pyplot as plt
        import seaborn as sns
        domains = domains.copy()
        if feature == 'Length' and 'Length' not in domains.columns:
            domains['Length'] = domains['End']-domains['Start']
            domains['Length'] = domains['Length'].astype(int)
        if kind == 'hist':
            if bins == 'autolength':
                bins = range(int(-self.resolution/2),
                             int(max(domains[feature])+3*self.resolution/2),
                             int(self.resolution))
            else:
                bins = 10
            g = sns.FacetGrid(data=domains, col=group, col_wrap=3)
            g.map(plt.hist, feature, bins=bins)

            if autoxlim:
                def set_xlim_ax(*args, **kwargs):
                    plt.xlim([-self.resolution/2,
                              max(domains[feature])+self.resolution*3/2])
                g.map(set_xlim_ax)
            return g
        elif kind in ['point', 'bar', 'box', 'violin', 'strip']:
            return sns.factorplot(x=group, y=feature, data=domains,
                                  kind=kind, *args, **kwargs)
        elif kind in ['count']:
            return sns.factorplot(x=group, data=domains,
                                  kind=kind, *args, **kwargs)
        else:
            raise ValueError, 'Unsupported kind'

    def compare_intervals(self, intervals1, intervals2, spec_funcs=[],
                          precision_both=None, precision_each=None,
                          precision_center=None, precision_length=None,
                          N=False):
        '''
        Accepts two pandas datasets with intervals (chromosome-based) and
        returns identical, present in only the first and only the second
        dataset. You can specify precision_both, precision_each,
        precision_length and/or precision_center for approximate comparison.
        The first relates to the sum of differences in both coordinates - e.g.
        precision_both=40000 will see all of these examples as identical, as
        well as perfectly matching intervals:
        (('chr1', 0, 200000), ('chr1', 20000, 200000)) -> sum(differences) = 20000
        (('chr1', 0, 200000), ('chr1', 0, 240000)) -> sum(differences) = 40000
        (('chr1', 20000, 200000), ('chr1', 0, 180000)) -> sum(differences) = 40000
        (('chr1', 20000, 200000), ('chr1', 0, 180000)) -> sum(differences) = 40000
        etc.
        The second relates to each difference, so each of them mustn't exceed
        the precision.
        Center relates to distance between centers of the intervals.
        Length relates to lengths of the intervals.
        Other functions can be specified in a list of *spec_funcs*. They need
        to accept two arguments: coordinates1 and coordinates2, tuples of
        coordinates (chr, start, end). They should return a boolean value.
        If any functions are specified, they are used together with default
        functions as identified by various precision arguments. If none of
        presicion arguments are specified, only special functions will be used.
        You don't need to compare chromosomes inside the functions, it is done
        anyway.
        Default - exact comparison.
        If compares exactly, returns 3 DataFrames: shared, unique1 and unique2.
        If not exactly, returns 4 DataFrames: shared1, shared2, unique1 and
        unique2. If *N*, returns only the numbers of intervals in each
        DataFrame.
        '''
        truefunc = lambda x, y: True
        if precision_both:
            def both_func(coord1, coord2):
                chrom1, start1, end1 = coord1
                chrom2, start2, end2 = coord2
                return abs(start1-start2) + abs(end1-end2) <= precision_both
        else:
            both_func = truefunc

        if precision_each:
            def each_func(coord1, coord2):
                chrom1, start1, end1 = coord1
                chrom2, start2, end2 = coord2
                return all(np.abs(np.array((start1-start2, end1-end2))) <= precision_each)
        else:
            each_func = truefunc

        if precision_center:
            def center_func(coord1, coord2):
                chrom1, start1, end1 = coord1
                chrom2, start2, end2 = coord2
                center1 = (start1+end1)/2
                center2 = (start2+end2)/2
                return abs(center2 - center1) <= precision_center
        else:
            center_func = truefunc

        if precision_length:
            def length_func(coord1, coord2):
                chrom1, start1, end1 = coord1
                chrom2, start2, end2 = coord2
                len1 = end1-start2
                len2 = end2-start2
                return abs(len2-len1) <= precision_length
        else:
            length_func = truefunc

        funcs = [both_func, each_func, center_func, length_func] + spec_funcs

        def remove_identical(ds1_list, ds2_list):
            '''
            Idea comes from here:
            http://stackoverflow.com/a/29464365/1304161
            '''
            ds1 = set(ds1_list)
            ds2 = set(ds2_list)
            intervals1_unique = pd.DataFrame(list(ds1.difference(ds2)),
                                             columns=intervals1.columns)
            intervals1_unique.sort(columns=['Chromosome', 'Start', 'End'],
                                   inplace=True)
            intervals2_unique = pd.DataFrame(list(ds2.difference(ds1)),
                                             columns=intervals2.columns)
            intervals2_unique.sort(columns=['Chromosome', 'Start', 'End'],
                                   inplace=True)
            return intervals1_unique, intervals2_unique

        ds1 = list(tuple(line) for line in intervals1.values)
        ds2 = list(tuple(line) for line in intervals2.values)

        if not any([precision_each, precision_both, precision_center,
                    precision_length]) and not spec_funcs:
            shared = pd.merge(intervals1, intervals2,
                              on=['Chromosome', 'Start', 'End'], how='inner')
            shared.sort(columns=['Chromosome', 'Start', 'End'], inplace=True)
            intervals1_unique, intervals2_unique = remove_identical(ds1, ds2)
            if not N:
                return shared, intervals1_unique, intervals2_unique
            else:
                return len(shared), len(intervals1_unique),\
                                                         len(intervals2_unique)

        shared1 = []
        shared2 = []
        for coord1 in ds1:
            for coord2 in ds2:
                if coord1[0] != coord2[0]:
                    continue
                if all([f(coord1, coord2) for f in funcs]):
                    shared1.append(coord1)
                    shared2.append(coord2)
                    break

        shared1 = pd.DataFrame(shared1, columns=intervals1.columns)
        shared2 = pd.DataFrame(shared2, columns=intervals2.columns)

        shared1_list = [tuple(line) for line in shared1.values]
        shared2_list = [tuple(line) for line in shared2.values]

        intervals1_unique = remove_identical(ds1, shared1_list)[0]
        intervals2_unique = remove_identical(ds2, shared2_list)[0]

        if not N:
            return shared1, shared2, intervals1_unique, intervals2_unique
        else:
            return len(shared1), len(shared2), len(intervals1_unique),\
                                                         len(intervals2_unique)
    
    def get_sum_in_rect(self, intervals, data):
        start1, end1, start2, end2 = [int(i) for i in intervals]
        return np.sum(data[start1:end1, start2:end2])
    
    def get_density(self, interval, data, frac=False, norm='square'):
        '''
        Get sum of all interactions (i.e. density) in a square from *interval*.
        If *frac*, divides the sum of interactions by all interactions of the
        region.
        '''
        start, end = tuple(interval)
        if any(np.isnan([start, end])):
            return np.NaN
        start, end = int(start), int(end)
        if not frac:
            s = np.sum(data[start:end, start:end])/2
        else:
            s = np.sum(data[start:end, start:end]) /2 / \
                       np.sum(data[start:end])
        if norm == 'length':
            return s/(end-start)
        elif norm == 'square':
            return s*2/(end-start)**2
        elif not norm:
            return s

    def get_intervals_density(self, intervals, data, data_norm=True,
                              intervals_norm=False, frac=False,
                              in_cols=['Chromosome', 'Start', 'End']):
        '''
        Calculates "density" of all intervals from heatmap data, as measured by
        number of ligation junctions inside the intervals. By default only
        uses half of the square matrix. Optionally normalizes as reads per
        million (*data_norm*). Optionally normalizes sum in an interval by the
        sum of all columns of the interval (*frac*) and/or by the length or
        squared length of each interval (*intervals_norm={'length'|'square'}*).
        '''
        if data_norm:
            data /= np.sum(data)            
        ints = self.chr_intervals_to_genome(intervals, in_cols)
        if frac:
            f = partial(self.get_density, data=data, frac=True,
                        norm=intervals_norm)
        else:
            f = partial(self.get_density, data=data, norm=intervals_norm)
        return ints.apply(f, axis=1)

    def get_border_strength(self, coordinates, data, data_norm=True):
        '''
        Calculate insulation strength between two intervals (e.g. TADs).
        Calculated by dividing the sum of interactions inside two regions by
        the sum of interactions between the regions. Omits the spacer between
        them.
        *Coordinates* - (start1, end1, start2, end2) in bin numbers.
        (For the sake of convenience to use in DF.apply() method)
        *data* - Hi-C interactions matrix
        '''
        start1, end1, start2, end2 = coordinates
        if data_norm:
            data /= np.sum(data)
        sum1 = self.get_density((start1, end1), data)
        sum2 = self.get_density((start2, end2), data)
        sum_inter = np.sum(data[start1:end1, start2:end2])/2
        return pd.Series({'Start1':start1, 'End1':end1, 'Start2':start2, 'End2':end2,
                          'Density1':sum1, 'Density2':sum2, 'InterDensity':sum_inter,
                          'Strength':(sum1 + sum2)/sum_inter})

    def get_borders_strength(self, intervals, data):
        '''
        Calculate strength of domain borders. By divides sum of
        interactions within two neighbouring domains by the sum of
        interactions between those domains (**self.get_border_strength()**).
        '''
        ints = intervals.copy()
        ints['Density'] = self.get_intervals_density(ints, data, data_norm=False)
        ints.rename(columns={'Start':'Start1', 'End':'End1', 'Density':'Density1'},
                    inplace=True)
        grouped = ints.groupby('Chromosome')
        ints[['Start2', 'End2', 'Density2']] = grouped[['Start1', 'End1',
                                                        'Density1']].shift(-1)
        ints.dropna(inplace=True)
        ints[['Start', 'End']] = ints[['End1', 'Start2']]
        f = partial(self.get_sum_in_rect, data=data)
        ints_temp = pd.DataFrame()
        ints_temp[['Start1', 'End1']] = self.chr_intervals_to_genome(
                                        ints[['Chromosome', 'Start1', 'End1']],
                                        in_cols=['Chromosome', 'Start1', 'End1'],
                                        out_cols=['Start1', 'End1'])
        ints_temp[['Start2', 'End2']] = self.chr_intervals_to_genome(
                                        ints[['Chromosome', 'Start2', 'End2']],
                                        in_cols=['Chromosome', 'Start2', 'End2'],
                                        out_cols=['Start2', 'End2'])
        ints['interDensity'] = ints_temp[['Start1', 'End1',
                                               'Start2', 'End2']].apply(f, axis=1)
        ints['Strength'] = (ints['Density1'] + ints['Density2']) / ints['interDensity']
        return ints
        
    def get_insulation_scores(self, data, arm_length=3, gap=0):
        '''
        Calculate insulation scores. The value is equivalent to border strength, but
        calculated for the dataset **data**. Uses pairs of intervals with length of
        **arm_length** and a gap between them with length **gap**. Both are in number of
        bins, returned DataFrame is in coordinates.
        '''
        max_coord = len(data)+1
        starts1 = range(max_coord-arm_length-gap)
        ends1 = [i+arm_length for i in starts1]
        ints = pd.DataFrame({'Start':starts1, 'End':ends1})
        ints *= self.resolution
        ints = self.genome_intervals_to_chr(ints)
        ints['Density1'] = self.get_intervals_density(ints, data, data_norm=False)
        ints = ints.rename(columns={'Start':'Start1', 'End':'End1'})
        grouped = ints.groupby('Chromosome')
        ints[['Start2', 'End2', 'Density2']] = grouped[['Start1', 'End1',
                                                        'Density1']].shift(-arm_length-gap)
        ints_temp = pd.DataFrame()
        ints_temp[['Start1', 'End1']] = self.chr_intervals_to_genome(
                                        ints[['Chromosome', 'Start1', 'End1']],
                                        in_cols=['Chromosome', 'Start1', 'End1'],
                                        out_cols=['Start1', 'End1'])
        ints_temp[['Start2', 'End2']] = self.chr_intervals_to_genome(
                                        ints[['Chromosome', 'Start2', 'End2']],
                                        in_cols=['Chromosome', 'Start2', 'End2'],
                                        out_cols=['Start2', 'End2'])
        ints_temp = ints_temp.dropna()
        f = partial(self.get_sum_in_rect, data=data)
        ints = ints.reset_index(drop=True)
        ints['interDensity'] = ints_temp.apply(f, axis=1).reset_index(drop=True)
        ints['Insulation'] = (ints['Density1']+ints['Density2'])/ints['interDensity']
        return ints