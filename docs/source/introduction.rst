Introduction
============

``uranography`` is a collection of classes that supply tools for making
``bokeh`` plots of data on the celestial sphere. It is primarily intended for 
use in interactive environments such as jupyter notebooks or dashboards.

The ``uranography`` python module provides a collection of classes, each of which supports a different map projection.
Some of these map projections include controls that let them simulate traditional tools used in astronomy:
the Lambert azimutal equal area projection simulates a `planisphere <https://en.wikipedia.org/wiki/Planisphere>`_, and the othrogrophic projection simulates an `armillary sphere <https://en.wikipedia.org/wiki/Armillary_sphere>`_.
Each of these classes includes:

- a member (``plot``) that is a perfectly normal instance of `bokeh.plotting.Figure <https://docs.bokeh.org/en/latest/docs/reference/plotting/figure.html#bokeh.plotting.figure>`_. Through this instance, users can create and modify the full high level API provided by ``bokeh``, and the full set of ``bokeh`` interactive tools can be applied.
- ``bokeh`` `transforms <https://docs.bokeh.org/en/latest/docs/reference/transform.html>`_ that apply the map projection in the client. By using ``bokeh`` client-side transforms to handle map projections, ``sphremap`` supports the application and interactive adjustment of map projection parameters without commutication with a server (or ``python`` process of any sort). This means that updates to projections (e.g. with a control slider in a web page) do not incur any overhead for communication with a server, thereby allowing smooth motion and interaction even with slow network connections. Furthermore, plots can be saved as strings and either loaded into browsers direcly from disk, or embedded into other web pages, and the plots will remain active and interactive even if ``python`` process that produced the plots no longer exists.
- Methods for plotting `healpix <https://healpix.jpl.nasa.gov/>`_ and `healsparse <https://github.com/LSSTDESC/healsparse>`_ arrays.
- Methods for adding a variety of features commonly used in maps of the sky, including:
  
  - graticules in equatorial coordinates (R.A. and declination),
  - graticules in horizon coordinates (altitude and azimuth),
  - the ecliptic plane,
  - the horizon (or circles of any other altitude),
  - circles of a sphere (great or small circes) with arbitrary centers and radii, or arcs of such circles beginning and ending at arbitrary bearings from the center.
  - the galactic plane, and
  - stars from the `Yale Bright Star Catalog <http://tdc-www.harvard.edu/catalogs/bsc5.html>`_.

These methods add standard named glyphs, data sources, and renderer models to the instance of ``bokeh.plotting.Figure``, so they can be selected by name from the ``plot`` member and adjusted and refined using the standard ``bokeh`` API.

The ability to simulate a planisphere or armillary sphere make plots created by ``uranography`` particularly useful for planning observing.
Where map projections commonly used in cosmology (e.g. Mollweide) heavily distort features and periodicities of primary importance to observation planning (e.g. the horizon),
such features are more naturally represented in a simulated planisphere.
When used with fully interactive sliders, an orthographic projection helps a human build a mental three dimensional picture of the sky difficult to achieve with other map projections, or with non-interactive orthographic projections.

For example, compare distortion of important circles on the sky (e.g. the horizon and ecliptic) in the commonly used Mollweide projection:

.. raw:: html
    :file: mollweide.html

with their appearance and motion on a smoothly interactive image of a `planisphere <https://en.wikipedia.org/wiki/Planisphere>`_ :

.. raw:: html
    :file: planisphere.html

or with a mental picture created with the aide of a virtual `armillary sphere <https://en.wikipedia.org/wiki/Armillary_sphere>`_, after playing around with the different sliders:

.. raw:: html
    :file: asphere.html

All three of these figures were created using ``uranography``.