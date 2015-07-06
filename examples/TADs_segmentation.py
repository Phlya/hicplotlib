
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

gammas = list(np.linspace(0, 0.9, 10))+range(1, 10)+range(10, 155, 5)
s2_tads = gi.genome_intervals_to_chr(gi.find_TADs(s2, gammalist=gammas, segmentation='potts'))
s2_tads['Segmentation'] = 'Potts'


# And now let's try Armatus segmentation (implementation of the Armatus algorithm from (Filippova et al, 2013)). This takes much longer.

# In[6]:

s2_tads_arm = gi.genome_intervals_to_chr(gi.find_TADs(s2, gammalist=gammas, segmentation='armatus'))
s2_tads_arm['Segmentation'] = 'Armatus'


# Let's combine the data from both algorithms in a long-style DataFrame so it's easier to analyse it.

# In[7]:

s2_tads_all = pd.concat([s2_tads, s2_tads_arm])
s2_tads_all.sort(['Segmentation', 'Gamma', 'Chromosome', 'Start', 'End'], inplace=True)


# Sometimes the algorithms produce TADs with no interactions inside from empty regions of the map. Let's calculate TADs densities (as a sum of interactions inside the TAD) and only keep TADs with Density>0. We'll use the Density values themselves later.

# In[8]:

s2_tads_all['Density'] = gi.get_intervals_density(s2_tads_all[['Chromosome', 'Start', 'End']], s2)
s2_tads_all = s2_tads_all.query('Density>0').copy() #Copying to shut up the warning later


# Let's now have a look at those TADs by segmentation method separately.

# In[9]:

s2_tads_all.query('Segmentation=="Potts"')


# In[10]:

s2_tads_all.query('Segmentation=="Armatus"')


# Seems like it worked, but we only saw the head and tail of the tables... Let's just now see some statistics across all gamma values we used. First let's have a look at lengths of called TADs. Coverage shows percentage of the genome covered by TADs.

# In[11]:

gi.describe_TADs(s2_tads_all.query('Segmentation=="Potts"'))


# In[12]:

gi.describe_TADs(s2_tads_all.query('Segmentation=="Armatus"'))


# OK, now let's use another column for statistics - Score.

# In[13]:

gi.describe_TADs(s2_tads_all.query('Segmentation=="Potts"'), feature='Score')


# In[14]:

gi.describe_TADs(s2_tads_all.query('Segmentation=="Armatus"'), feature='Score')


# Nice! Now let's visualise some of these things easily in a pretty way with seaborn.

# Let's see how many TADs we get depending on gamma and segmentation method.

# In[15]:

sns.countplot(x='Gamma', hue='Segmentation', data=s2_tads_all)
plt.xticks(rotation=90)
plt.show()


# Let's compare the lengths of TADs depending on gamma with these two segmentation methods.

# In[16]:

s2_tads_all['Length'] = s2_tads_all['End']-s2_tads_all['Start']
s2_tads_all['Length'] = s2_tads_all['Length'].astype(int)


# In[17]:

sns.pointplot(x='Gamma', y='Length', hue='Segmentation', data=s2_tads_all, zorder=15)
sns.stripplot(x='Gamma', y='Length', hue='Segmentation', data=s2_tads_all, jitter=True, zorder=1, alpha=0.5)
handles, labels = plt.gca().get_legend_handles_labels() #Have to do this because of a weird legend behavior otherwise...
plt.legend(handles[:2], labels[:2], loc='upper right')
plt.xticks(rotation=90)
plt.yscale('log')
plt.show()


# And again, let's look at the Score in the same way.

# In[18]:

sns.pointplot(x='Gamma', y='Score', hue='Segmentation', data=s2_tads_all, zorder=15)
sns.stripplot(x='Gamma', y='Score', hue='Segmentation', data=s2_tads_all, jitter=True, zorder=1, alpha=0.5)
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[:2], labels[:2], loc='lower left')
plt.xticks(rotation=90)
plt.show()


# Not really visible for Armatus, let's plot it separately now.

# In[19]:

sns.pointplot(x='Gamma', y='Score', data=s2_tads_all.query('Segmentation=="Armatus"'), zorder=15)
sns.stripplot(x='Gamma', y='Score', data=s2_tads_all.query('Segmentation=="Armatus"'), color='b', jitter=True, zorder=1, alpha=0.5)
handles, labels = plt.gca().get_legend_handles_labels()
plt.xticks(rotation=90)
plt.ylim(ymin=0)
plt.title('Armatus')
plt.show()


# Now let's compare TAD densities, i.e. sum of interactions inside each TAD, across gamma values. As you might remember, we have already calculated it before to remove empty TADs.

# In[20]:

sns.pointplot(x='Gamma', y='Density', hue='Segmentation', data=s2_tads_all, zorder=15)
sns.stripplot(x='Gamma', y='Density', hue='Segmentation', data=s2_tads_all, jitter=True, zorder=1, alpha=0.5)
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[:2], labels[:2], loc='upper right')
plt.xticks(rotation=90)
plt.ylim(ymin=0)
plt.show()


# Looks a little similar to an exaggerated plot with lengths we got earlier. Well, obviously, the way we defined density is strongly dependent on TAD length! Let's try to fix it: we shall divide the Density by Length^2 and hopefully get rid of this bias at least on average. You can actually do it from the very beginning by passing norm='square' to gi.get_intervals_density() method.

# In[21]:

s2_tads_all['Density'] /= s2_tads_all['Length']**2


# Both algorithms are somehow based on "density" calculation, so we see a clear relationship between Gamma and Density, especially for the more stable Potts segmentation. High density values seem to be a good clue for finding an optimal Gamma value.

# In[22]:

sns.pointplot(x='Gamma', y='Density', hue='Segmentation', data=s2_tads_all, zorder=15)
sns.stripplot(x='Gamma', y='Density', hue='Segmentation', data=s2_tads_all, jitter=True, zorder=1, alpha=0.5)
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[:2], labels[:2], loc='upper right')
plt.xticks(rotation=90)
plt.ylim(ymin=0)
plt.show()


# Now calculating TAD borders strength. Let's use pandas.groupby to apply it to all Gamma values and both segmentation methods and keep information about them at the same time.
# 
# Strength by default is calculated as ratio of sum of intra-TAD interactions of two neighbouring TADs to inter-TAD interactions of those TADs with each other.

# In[23]:

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

# In[24]:

sns.pointplot(x='Gamma', y='Strength', hue='Segmentation', data=s2_strengths_all, zorder=15)
sns.stripplot(x='Gamma', y='Strength', hue='Segmentation', data=s2_strengths_all, jitter=True, zorder=1, alpha=0.5)
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[:2], labels[:2], loc='upper left')
plt.xticks(rotation=90)
plt.yscale('log')
plt.ylim(ymin=0.4)
plt.show()


# Interestingly, the dependency is not quite the same as for Density.
# 
# Overall, obviously, this dataset requires very different gamma values in these two algorithms: <=1 for Armatus and 70-100 for Potts. Remember to take into account number of TADs you are getting as well as their properties, such as length. Unfortunately there is no single rule and gamma has to be defined empirically. The best judge here is our eye when comparing TAD calls to the heatmap. This will be covered in a different example.
