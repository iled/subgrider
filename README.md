# subgrider
Manipulations and statistics of grids in the GSLIB format

---

**subgrider** is a tool to perform some simple manipulations, data extraction, and computation of some statistics, of data structured in the [Geo-EAS format](http://www.gslib.com/gslib_help/format.html) (a.k.a. GSLIB format). This is a regular grid in 3-D.

It was developed during the work for my MS thesis, entitled *Applying Spatial Bootstrap and Bayesian Update in uncertainty assessment at oil reservoir appraisal stages* (original title is in Portuguese), defended in [Instituto Superior Técnico](http://www.ist.eu), December 2010, Lisbon, Portugal. The corresponding extended abstract (in English) is available in arXiv: [arXiv:1702.04450](https://arxiv.org/abs/1702.04450). It had some further development in the following year, while I was working with integration of coarse and fine data (BGeost - [Liu and Journel, 2009](http://www.sciencedirect.com/science/article/pii/S0098300408001386)).

After all this time, I finally make it available here, mainly for educational purposes (maybe just for myself). The code is obviously poorly mainted, it only includes a simple text interface, which is in Portuguese, as well as the few comments in the code (eventually I may translate it to English). I cannot be held liable for any of its use or consequences of its use. It is released under GPLv3. I was not using git at the time, but I tried to to recreate some of the development history, including the [binaries generated](https://github.com/iled/subgrider/releases) at that time.

I am happy to receive and reply to any [issue](https://github.com/iled/subgrider/issues) on it.

# Features

- Extract data
  - Vertical line, given coordinates (x, y)
  - Multiple vertical lines, given file with coordinates (x, y)
  - Subgrid, given a cuboid volume within the initial grid (*subgrider*!)

- Data manipulation
  - Create a mask (0 or 1)
  - Bend and unbend, given a map with the shifts
  - Upscale (create coarse data (blocks) in BGeost format)
  - Convert BGeost format to GSLIB point set

- Sampling
  - Random sampling of vertical line, given file with coordinates (x, y)

- Statistics
  - Mean and variance of simulated maps
    - Input - several grids with the same dimensions
    - Output - one grid with the same dimension which values are the mean or variance

- Bayesian update
  - Percentile
  - Quantile
  - Likelihood
