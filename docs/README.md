# uranography

## Description

Documentation available at https://uranography.lsst.io

`uranography` is a collection of tools for making plots on the celestial sphere with `bokeh`.
It is primarily intended for  use in interactive environments such as jupyter notebooks or dashboards.

It supports a handful of map projections, and interactive tools that adjust them in ways that simulate traditional tools
used in astronomy:
the Lambert azimutal equal area projection simulates a [planisphere](https://en.wikipedia.org/wiki/Planisphere),
and the othrogrophic projection simulates an [armillary sphere](https://en.wikipedia.org/wiki/Armillary_sphere).

The `uranography` module includes tools to produce:
- `bokeh` [`transforms`](https://docs.bokeh.org/en/latest/docs/reference/transform.html) that apply the map projection in the client.
- projections of [`healpix`](https://healpix.jpl.nasa.gov/) and [healsparse](https://github.com/LSSTDESC/healsparse) arrays.  
- graticules in equatorial coordinates (R.A. and declination),
- graticules in horizon coordinates (altitude and azimuth),
- the ecliptic plane,
- the horizon (or circles of any other altitude),
- circles of a sphere (great or small circes) with arbitrary centers and radii, or arcs of such circles beginning and ending at arbitrary bearings from the center.
- the galactic plane, and
- stars from the [Yale Bright Star Catalog](http://tdc-www.harvard.edu/catalogs/bsc5.html).

These methods add `bokeh` glyphs, data sources, and renderer models with standard or manually assigned assigned names to an instance of `bokeh.plotting.Figure`, so they can be selected by name and adjusted and refined using the standard `bokeh` API.

More documentation is available at [uranography.lsst.io](https://uranography.lsst.io).