# uranography

## Introduction

`uranography` is a collection of classes that supply tools for making
`bokeh` plots of data on the celestial sphere. It is primarily intended for 
use in interactive environments such as jupyter notebooks or dashboards.

The `uranography` python module provides a collection of classes, each of which supports a different map projection.
Some of these map projections include controls that let them simulate traditional tools used in astronomy:
the Lambert azimutal equal area projection simulates a [planisphere](https://en.wikipedia.org/wiki/Planisphere), and the othrogrophic projection simulates an [armillary sphere](https://en.wikipedia.org/wiki/Armillary_sphere).
Each of these classes includes:

- a member (`plot`) that is a perfectly normal instance of [`bokeh.plotting.Figure`](https://docs.bokeh.org/en/latest/docs/reference/plotting/figure.html#bokeh.plotting.figure). Through this instance, users can create and modify the full high level API provided by `bokeh`, and the full set of `bokeh` interactive tools can be applied.
- `bokeh` [`transforms`](https://docs.bokeh.org/en/latest/docs/reference/transform.html) that apply the map projection in the client. By using `bokeh` client-side transforms to handle map projections, `uranography` supports the application and interactive adjustment of map projection parameters without commutication with a server (or `python` process of any sort). This means that updates to projections (e.g. with a control slider in a web page) do not incur any overhead for communication with a server, thereby allowing smooth motion and interaction even with slow network connections. Furthermore, plots can be saved as strings and either loaded into browsers direcly from disk, or embedded into other web pages, and the plots will remain active and interactive even if the `python` process that produced the plots no longer exists.
- Methods for plotting [`healpix`](https://healpix.jpl.nasa.gov/) and [healsparse](https://github.com/LSSTDESC/healsparse) arrays.
- Methods for adding a variety of features commonly used in maps of the sky, including:
  
  - graticules in equatorial coordinates (R.A. and declination),
  - graticules in horizon coordinates (altitude and azimuth),
  - the ecliptic plane,
  - the horizon (or circles of any other altitude),
  - circles of a sphere (great or small circes) with arbitrary centers and radii, or arcs of such circles beginning and ending at arbitrary bearings from the center.
  - the galactic plane, and
  - stars from the [Yale Bright Star Catalog](http://tdc-www.harvard.edu/catalogs/bsc5.html).

These methods add standard named glyphs, data sources, and renderer models to the instance of `bokeh.plotting.Figure`, so they can be selected by name from the `plot` member and adjusted and refined using the standard `bokeh` API.

## Installation

First, clone the `uranography` repository:

```
git clone git@github.com:lsst/uranography.git
cd uranography
```

Create a conda environment for it:

```
conda create --channel conda-forge --name uranography --file requirements.txt python=3.11
```

If you want to run tests, install the test requirements as well:

```
conda activate uranography
conda install -c conda-forge --file=test-requirements.txt
```

Install the `uranography` project into this environment:

```
cd <wherever you checked out uranography>/uranography
pip install -e .
```

Install the kernel from the new jupyter environment:

```
python -m ipykernel install --user --name=uranography
```

To make the documentation:

```
cd docs
make html
```

The root of the refenece documentation will then be `docs/_build/html/source/index.html`.

For tutorial documentation, follow the tutorial in `uranography/notebooks/uranography.ipynb` in `jupyter`.