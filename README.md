hicplotlib
==========
A draft for a library to analyze and plot Hi-C interactions data. For now works well with whole-genome and by-chromosomal plotting and produces reasonably nice pictures. The good thing is that everything is based on matplotlib and all plotting options can be modified "externally", and all objects of a figure can be easily accessed. No documentation for now, but names of methods and comments should be pretty straightforward.

Provides some analysis tools, such as TADs finding using greendale, their descriptive statistics and comparison, basic observed/expected calculation. Can calculate TADs density, i.e. sum of Hi-C interactions, and estimate TAD borders strength.

Installation: pip install git+git://github.com/Phlya/hicplotlib.git
