# -*- coding: utf-8 -*-
"""
This is a module to help with plotting Hi-C data.
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
import numpy as np
from scipy import ndimage

class HiCPlot(object):
    def __init__(self, settings=None):
        if settings is not None:
            self.settings = settings
        self.settings.extract_settings(self, self.settings)
        self.name = ''
        self.name2 = ''
        self.cmap = cmap.get_cmap()
        self.data = None
        self.data_dict = {}

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
        
    def get_data(self, name):
        """
        Retrieve data by it's name. Calls self.data_dict[name]
        """
        return self.data_dict[name]
    
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
        if savepath:
            plt.ioff()
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
            norm = plt.Normalize(data[np.isfinite(data)].min(),
                                 data[np.isfinite(data)].max())
            data = colormap(norm(data))
            for multiple in diagonal_markers.keys():
                for start, end in self.boundaries:
                    for i in range(start, end, multiple//self.resolution):
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
        if savepath:
            plt.savefig(savepath)

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