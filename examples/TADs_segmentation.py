
# coding: utf-8

# In[1]:

import numpy as np
import hicplotlib as hpl
import matplotlib.pyplot as plt
import seaborn as sns #Recommend upgrading to new seaborn 0.6
import pandas as pd
from __future__ import division
sns.set_color_codes()
get_ipython().magic(u'pylab inline')
#I suggest you run it without inline mode as this way figures are quite cluttered


# In[2]:

settings = hpl.HiCParameters()
settings.resolution=20000
settings.set_chromosomes(['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX'])
settings.set_chromosomes_lengths([23011544, 21146708, 24543557, 27905053, 1351857, 22422827])
#Can be also done by parsing a file with chromosome names and lengths with settings.set_chromosomes_from_chrfile
gi = hpl.GenomicIntervals(settings)


# Loading a Hi-C matrix, do it any way you prefer to get your data.

# In[3]:

s2 = np.loadtxt('../../../Documents/biology/Drosophila_cells_Hi-C/S2/IC-heatmap-20K.mtx')


# This is how it looks, just a symmetrical matrix. The diagonal has been removed here, so zeroes in the middle.

# In[4]:

s2


# Now let's call the TADs. It is done via greendale (https://bitbucket.org/nvictus/greendale/) which needs to be installed separately. Let's first use Potts segmentation

# In[5]:

s2_tads = gi.genome_intervals_to_chr(gi.find_TADs(s2, gammalist=list(np.linspace(0, 0.9, 10))+range(1, 10)+range(10, 155, 5), segmentation='potts'))
s2_tads['Segmentation'] = 'Potts'


# And now let's try Armatus segmentation (implementation of the Armatus algorithm from (Filippova et al, 2013). This takes much longer.

# In[6]:

s2_tads_arm = gi.genome_intervals_to_chr(gi.find_TADs(s2, gammalist=list(np.linspace(0, 0.9, 10))+range(1, 10)+range(10, 155, 5), segmentation='armatus'))
s2_tads_arm['Segmentation'] = 'Armatus'


# Let's now have a look at those TADs.

# In[17]:

s2_tads.reset_index(drop=True)


# In[18]:

s2_tads_arm.reset_index(drop=True)


# Seems like it worked, but we only saw the head and tail of the tables... Let's just now see some statistics across all gamma values we used. Coverage shows percentage of the genome covered by TADs.

# In[19]:

gi.describe_TADs(s2_tads)


# In[20]:

gi.describe_TADs(s2_tads_arm)


# Nice! Now let's visualise some of these things nicely.

# In[21]:

s2_tads_all = pd.concat([s2_tads, s2_tads_arm])
s2_tads_all.sort(['Segmentation', 'Gamma', 'Chromosome', 'Start', 'End'], inplace=True)


# In[22]:

s2_tads_all['Length'] = s2_tads_all['End'].astype(int)-s2_tads_all['Start'].astype(int)


# Let's see how many TADs we get depending on gamma and segmentation method.

# In[23]:

sns.countplot(x='Gamma', hue='Segmentation', data=s2_tads_all)
plt.xticks(rotation=90)
plt.show()


# Let's compare the lengthes of TADs depending on gamma with these two segmentation methods.

# In[24]:

sns.pointplot(x='Gamma', y='Length', hue='Segmentation', data=s2_tads_all, zorder=15)
sns.stripplot(x='Gamma', y='Length', hue='Segmentation', data=s2_tads_all, jitter=True, zorder=1, alpha=0.5)
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[:2], labels[:2], loc='upper right')
plt.xticks(rotation=90)
plt.yscale('log')
plt.show()


# Now calculating TAD borders strength. Let's use pandas.groupby.
# 
# Strength by default is calculated as ratio of sum of intra-TAD interactions of two neighbouring TADs to inter-TAD interactions of those TADs with each other.

# In[25]:

s2_l = []
grouped = s2_tads_all.groupby(['Gamma', 'Segmentation'])
for (g, s), d in grouped:
    if len(d.index)>1:
        a = gi.get_borders_strength(d.copy(), s2)
    a['Gamma'] = g
    a['Segmentation'] = s
    s2_l.append(a)
s2_strengths_all = pd.concat(s2_l)


# Now let's compare the strengths we are getting.

# In[26]:

sns.pointplot(x='Gamma', y='Strength', hue='Segmentation', data=s2_strengths_all, zorder=15)
sns.stripplot(x='Gamma', y='Strength', hue='Segmentation', data=s2_strengths_all, jitter=True, zorder=1, alpha=0.5)
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[:2], labels[:2], loc='upper left')
plt.xticks(rotation=90)
plt.yscale('log')
plt.show()


# Make your own conclusions...
