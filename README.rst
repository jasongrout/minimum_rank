Minimum Rank Sage Library
=========================

One goal of this library is to be able to be used both by loading the necessary files in the Sage notebook via load statements, as well as installation as a python module in Sage.  Thus, it is important that the names in the separate submodules do not trample on each other (since the load command imports everything at the top level)

See http://sage.cs.drake.edu/home/pub/69/ for an example of how to load and use this library in Sage.  In particular, the following code in a Sage notebook cell will load this library::

  URL='http://github.com/jasongrout/minimum_rank/raw/minimum_rank_1_1_3/'
  files=['Zq_c.pyx','Zq.py','zero_forcing_64.pyx','zero_forcing_wavefront.pyx','minrank.py', 'inertia.py']
  for f in files:
      load(URL+f)
