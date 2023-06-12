"""Base classes for sky maps made using bokeh.
"""

from collections import OrderedDict, namedtuple
from collections.abc import Iterable
from copy import deepcopy
from warnings import warn
import uuid
import time

import numpy as np
import pandas as pd
import healpy as hp
import bokeh
import bokeh.plotting
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
import astropy.visualization
import panel as pn
from IPython.display import display

try:
    import healsparse
except ModuleNotFoundError:
    pass

from .readjs import read_javascript


ProjSliders = namedtuple("ProjSliders", ["alt", "az", "mjd"])

# Angles of excactly 90 degrees result in edge conditions, e.g.
# horizon circles with gaps depending on how trig is rounded.
# Define an "almost 90" to get consistent behaviour.
ALMOST_90 = np.degrees(
    np.arccos(2 * np.finfo(float).resolution)  # pylint: disable=E1101
)


class SphereMap:
    """Base for maps of the sphere.

    Parameters
    ----------
    plot : `bokeh.plotting._figure.figure`, optional
        Figure to which to add the map, by default None
    mjd : `float`, optional
        The Modified Julian Date
    location: `EarthLocation` or `str`, optional
        The location of the observatory, defaults to lsst.
    """

    alt_limit = 0
    transform_js_fnames = ("coord_utils.js", "hp.js")
    transform_js_call = "return hpTransform()"
    update_js_fnames = ("coord_utils.js", "hp.js")
    update_js_command = "updateHpData()"
    max_star_glyph_size = 15
    proj_slider_keys = ["mjd"]
    default_title = ""
    x_col = "x_hp"
    y_col = "y_hp"
    default_graticule_line_kwargs = {"color": "darkgray", "name": "graticule_glyph"}
    default_ecliptic_line_kwargs = {"color": "green", "name": "ecliptic_glyph"}
    default_galactic_plane_line_kwargs = {
        "color": "blue",
        "name": "galactic_plane_glyph",
    }
    default_horizon_line_kwargs = {
        "color": "black",
        "line_width": 6,
        "name": "horizon_glyph",
    }

    def __init__(self, plot=None, mjd=None, location="Cerro Pachon"):
        self.mjd = Time.now().mjd if mjd is None else mjd

        if isinstance(location, EarthLocation):
            self.location = location
        else:
            self.location = EarthLocation.of_site(location)

        if plot is None:
            self.plot = bokeh.plotting.figure(
                frame_width=512,
                frame_height=512,
                match_aspect=True,
                title=self.default_title,
            )
        else:
            self.plot = plot

        self._viewable = None
        self._figure = None
        self._hpix_cmap = None
        self.notebook_handle = None
        self.star_data = None
        self.plot.axis.visible = False
        self.plot.grid.visible = False

        self.visible_slider_names = []
        self.update_functions = []
        self.force_update_time = bokeh.models.Div(
            text=str(time.time_ns()), name="update_time", visible=False
        )
        self.add_sliders()

    @property
    def healpix_data(self):
        """The healpix data source (deprecated!)"""
        warn('Use self.plot.select(name="hpix_ds") instead', DeprecationWarning)
        data_source = self.plot.select(name="hpix_ds")
        return data_source

    @property
    def healpix_renderer(self):
        """The healpix renderers (deprecated!)"""
        warn('Use self.plot.select(name="hpix_renderer") instead', DeprecationWarning)
        renderer = self.plot.select(name="hpix_renderer")
        return renderer

    @property
    def healpix_glyph(self):
        """The healpix data source (deprecated!)"""
        warn(
            'Use self.plot.select(name="hpix_renderer").glyph instead',
            DeprecationWarning,
        )
        renderers = self.plot.select(name="hpix_renderer")
        try:
            glyph = renderers.glyph
        except AttributeError:
            glyph = [r.glyph for r in renderers]
        return glyph

    @property
    def star_data_source(self):
        """Return the star data source (deprecated)."""
        warn(
            'Use self.plot.select(name="bright_star_ds") instead',
            DeprecationWarning,
        )
        star_data = self.plot.select(name="bright_star_ds")
        return star_data

    @property
    def healpix_cmap(self):
        """Return the bokeh color map dict used for bokeh."""
        cmap_transform = self.plot.select(name="hpix_color_transform")
        if len(cmap_transform) == 0 and self._hpix_cmap is not None:
            cmap = self._hpix_cmap
        else:
            cmap = {"field": "value", "transform": cmap_transform}
        return cmap

    @property
    def figure(self):
        """Return the bokeh figure for the map.

        Returns
        -------
        figure : `bokeh.models.layouts.LayoutDOM`
            A viewable bokeh figure.
        """
        if self._figure is None:
            shown_sliders = [self.sliders[sn] for sn in self.visible_slider_names]
            self._figure = bokeh.layouts.column(
                self.force_update_time, self.plot, *shown_sliders
            )

        return self._figure

    @property
    def viewable(self):
        """A panel viewable that shows the figure."""
        if self._viewable is None:
            self._viewable = pn.Row(self.figure)

        return self._viewable

    def update(self, viewable=None):
        """Trigger the javascript update callback.

        Parameters
        ==========
        viewable : `panel.viewable.Viewable`
            If not None, push the notebook to the frontend set by the viewable.
            Defaults to self.viewable
        """
        self.force_update_time.text = str(time.time_ns())
        if viewable is None:
            viewable = self._viewable

        if viewable is not None:
            pn.io.push_notebook(viewable)

    @property
    def lst(self):
        """Return the Local Sidereal Time."""
        lst = (
            Time(self.mjd, format="mjd", location=self.location)
            .sidereal_time("mean")
            .deg
        )
        return lst

    @lst.setter
    def lst(self, value):
        """Modify the MJD to match the LST, keeping the same (UT) day."""
        mjd_start = np.floor(self.mjd)
        lst_start = (
            Time(mjd_start, format="mjd", location=self.location)
            .sidereal_time("mean")
            .deg
        )
        self.mjd = mjd_start + ((value - lst_start) % 360) / 360.9856405809225

    @property
    def update_js(self):
        """Return javascript code to update the plots.

        Returns
        -------
        js_code : `str`
            Javascript code to update the bokeh model.
        """
        js_code = "\n".join(read_javascript(fn) for fn in self.update_js_fnames)
        js_code = "\n".join(
            [js_code, self.update_js_command, "data_source.change.emit()"]
        )
        return js_code

    def proj_transform(self, proj_coord, data_source=None, column_name=None):
        """Return a bokeh projection transformer.

        Parameters
        ----------
        proj_coord : `str`
            'x' or 'y', the projection coodinate to compute
        data_source : `bokeh.models.ColumnDataSource`
            The data source to project. Must have either 'ra' and 'decl'
            columns, or 'alt' and 'az'. All in degrees. Defaults to None,
            in which case the transform takes a column of ra, decl pairs.
        column_name : `str`
            The name of the column with ra, decl pairs. If None, columns
            with 'ra' and 'decl' or 'alt' and 'az' names instead. Defaults
            to None.

        Returns
        -------
        transform : `dict`
            With keys `field` and `transform`, suitable for passing to
            bokeh plotting functions.
        """
        try:
            mjd_slider = self.sliders["mjd"]
        except KeyError:
            mjd_slider = {"value": self.mjd}

        try:
            center_alt_slider = self.sliders["alt"]
        except KeyError:
            center_alt_slider = {"value": 90}

        try:
            center_az_slider = self.sliders["az"]
        except KeyError:
            center_az_slider = {"value": 0}

        js_code = "\n".join(read_javascript(fn) for fn in self.transform_js_fnames)
        js_code = "\n".join([js_code, self.transform_js_call])
        js_transform = bokeh.models.CustomJSTransform(
            args=dict(
                data_source=data_source,
                center_alt_slider=center_alt_slider,
                center_az_slider=center_az_slider,
                mjd_slider=mjd_slider,
                lat=self.location.lat.deg,
                lon=self.location.lon.deg,
                proj_coord=proj_coord,
            ),
            v_func=js_code,
        )

        if column_name is None:
            for column_name in ("ra", "alt"):
                if column_name in data_source.data:
                    break

        coord_transform = bokeh.transform.transform(column_name, js_transform)
        return coord_transform

    def x_transform(self, column_name):
        """Return a bokeh projection transformer for x

        Parameters
        ----------
        column_name : `str`
            The name of the column with ra, decl pairs

        Returns
        -------
        transform : `dict`
            With keys `field` and `transform`, suitable for passing to
            bokeh plotting functions.
        """
        coord_transform = self.proj_transform("x", column_name=column_name)
        return coord_transform

    def y_transform(self, column_name):
        """Return a bokeh projection transformer for y

        Parameters
        ----------
        column_name : `str`
            The name of the column with ra, decl pairs

        Returns
        -------
        transform : `dict`
            With keys `field` and `transform`, suitable for passing to
            bokeh plotting functions.
        """
        coord_transform = self.proj_transform("y", column_name=column_name)
        return coord_transform

    def eq_to_horizon(self, ra, decl, degrees=True, cart=True):  # pylint: disable=C0103
        """Convert equatorial to horizon coordinates

        Parameters
        ----------
        ra : `numpy.ndarray`
            Values for Right Ascension
        decl : `numpy.ndarray`
            Values for declination
        degrees : bool, optional
            Values are in degrees (if False, values are in radians),
            by default True
        cart : bool, optional
            Return cartesion coordinates rather than alt, az, by default True

        Returns
        -------
        coords : `list` [`np.ndarray`]
            Either alt, az (if cart=False) with az measured east of north,
            or x, y with +x pointing west and +y pointing north

        Azimuth is east of north
        """

        if isinstance(ra, pd.Series):
            ra = ra.values

        if isinstance(decl, pd.Series):
            decl = decl.values

        angle_unit = u.deg if degrees else u.rad
        eq_coords = SkyCoord(ra=ra * angle_unit, dec=decl * angle_unit, frame="icrs")
        obs_time = Time(self.mjd, format="mjd", scale="utc")
        horizon_frame = astropy.coordinates.AltAz(
            obstime=obs_time, location=self.location
        )
        horizon_coords = eq_coords.transform_to(horizon_frame)

        if cart:
            zd = np.pi / 2 - horizon_coords.alt.rad  # pylint: disable=C0103
            x = -zd * np.sin(horizon_coords.az.rad)  # pylint: disable=C0103
            y = zd * np.cos(horizon_coords.az.rad)  # pylint: disable=C0103

            invisible = horizon_coords.alt.deg < self.alt_limit
            x[invisible] = np.nan
            y[invisible] = np.nan

            return x, y

        if degrees:
            alt = horizon_coords.alt.deg
            az = horizon_coords.az.deg  # pylint: disable=C0103
            return alt, az
        else:
            alt = horizon_coords.alt.rad
            az = horizon_coords.az.rad  # pylint: disable=C0103
            return alt, az

    def make_healpix_data_source(
        self, hpvalues, nside=32, bound_step=1, nest=False, name="hpix_ds"
    ):
        """Make a data source of healpix values, corners, and projected coords.

        Parameters
        ----------
        hpvalues : `numpy.ndarray` or `healsparse.HealSparseMap`
            Healpixel values (RING pixel ordering unless
            a HealSparseMap is provided.)
        nside : int, optional
            healpixel nside for display, by default 32
        bound_step : int, optional
            number of boundary points for each side of each healpixel,
            by default 1
        nest : `bool`, optionol
            Is the healpix array provided in NEST ordering? Defaults to False.
            (if hpvalues is a HealSparseMap, nest is always True)
        name : `str`, option
            The bokeh name for the data source model, defaults to `hpix_ds`

        Returns
        -------
        hpix_datasource : `bokeh.models.ColumnDataSource`
            Data source for healpixel values and bounds.
        """
        try:
            # If this doesn't throw an exception, we
            # were passed an instance of healsparse.
            if nside < hpvalues.nside_sparse:
                values = hpvalues.degrade(nside)
            else:
                values = hpvalues
            hpids = values.valid_pixels
            finite_values = values[hpids]
            nside = values.nside_sparse
            nest = True
        except AttributeError:
            order = "NEST" if nest else "RING"
            values = np.copy(hpvalues)
            try:
                values[np.isnan(values)] = hp.UNSEEN
            except TypeError:
                # If isnan fails, value must have been of a type without
                # nans, so we don't need to do anything.
                pass

            values = hp.ud_grade(hpvalues, nside, order_in=order, order_out=order)
            values[values == hp.UNSEEN] = np.nan
            hpids = np.isfinite(values).nonzero()[0]
            finite_values = values[hpids]

        npix = len(hpids)
        npts = npix * 4 * bound_step

        hpix_bounds_vec = hp.boundaries(nside, hpids, bound_step, nest=nest)
        # Rearrange the axes to match what is used by hp.vec2ang
        hpix_bounds_vec_long = np.moveaxis(hpix_bounds_vec, 1, 2).reshape((npts, 3))
        ra, decl = hp.vec2ang(  # pylint: disable=C0103
            hpix_bounds_vec_long, lonlat=True
        )
        center_ra, center_decl = hp.pix2ang(nside, hpids, nest=nest, lonlat=True)

        # in hpix_bounds, each row corresponds to a healpixels, and columns
        # contain lists where elements of the lists correspond to corners.

        hpix = bokeh.models.ColumnDataSource(
            {
                "hpid": hpids,
                "value": finite_values,
                "center_ra": center_ra.tolist(),
                "center_decl": center_decl.tolist(),
                "x_hp": hpix_bounds_vec[:, 0, :].tolist(),
                "y_hp": hpix_bounds_vec[:, 1, :].tolist(),
                "z_hp": hpix_bounds_vec[:, 2, :].tolist(),
                "ra": ra.reshape(npix, 4 * bound_step).tolist(),
                "decl": decl.reshape(npix, 4 * bound_step).tolist(),
            },
            name=name,
        )

        # adding the projection column must come before
        # adding the callback (done by _finish_data_source),
        # or the callback to update the healpix corners
        # will not be created.
        self._add_projection_columns(hpix, nside)

        self.connect_controls(hpix)

        return hpix

    def _add_zero_projection_columns(self, hpix):
        """Add columns of sets of dummy coords to hold projection values."""
        num_healpix = len(hpix.data["ra"])
        points_per_healpix = len(hpix.data["ra"][0])
        zero_coord = [[0.0] * points_per_healpix] * num_healpix
        for column in self.x_col, self.y_col:
            if column not in hpix.data:
                hpix.add(deepcopy(zero_coord), name=column)
        hpix.add(deepcopy(zero_coord), name="z_orth")
        return hpix

    def _add_projection_columns(
        self, hpix, nside, projector=None  # pylint: disable=W0613
    ):
        """Add columns with projection values."""
        if projector is None:
            return self._add_zero_projection_columns(hpix)

        coord_cols = ["x_hp", "y_hp", "z_hp"]
        n_hpix = len(hpix.data["ra"])
        corners_per_hpix = len(hpix.data["ra"][0])
        corners = (
            hpix.to_df()
            .loc[:, ["hpid"] + coord_cols]
            .explode(column=coord_cols, ignore_index=True)
        )
        assert len(corners) == n_hpix * corners_per_hpix

        x, y = projector.vec2xy(  # pylint: disable=C0103
            corners.x_hp.astype(float),
            corners.y_hp.astype(float),
            corners.z_hp.astype(float),
        )

        hpix.add(x.reshape(n_hpix, corners_per_hpix).tolist(), name=self.x_col)
        hpix.add(y.reshape(n_hpix, corners_per_hpix).tolist(), name=self.y_col)

        return hpix

    def _make_graticule_points(
        self,
        min_decl=-80,
        max_decl=80,
        decl_space=20,
        min_ra=0,
        max_ra=360,
        ra_space=30,
        step=1,
        name="eq_graticule_ds",
    ):
        """Create points that define graticules

        Parameters
        ----------
        min_decl : int, optional
            Decl. of minimum R.A. graticule
            and ends of declination graticules (deg),
            by default -80
        max_decl : int, optional
            Decl. of maximum R.A. graticulas
            and ends of declination graticules (deg),
            by default 80
        decl_space : int, optional
            Spacing of decl. graticules (deg), by default 20
        min_ra : int, optional
            R.A. of first R.A. graticule (deg), by default 0
        max_ra : int, optional
            R.A. of last R.A. graticule (deg), by default 360
        ra_space : int, optional
            Spacing of R.A. graticules (deg), by default 30
        step : int, optional
            Spacing of points along graticules, by default 1
        name : `str`
            Name for the bokeh model, defaults to eq_graticule_ds.

        Returns
        -------
        graticule_points : `bokeh.models.ColumnDataSource`
            Bokeh data sources defining points in graticules.
        """

        # Bohek puts gaps in lines when there are NaNs in the
        # data frame. We will be platting many graticules, and we do not
        # want them connected, so define a "stop" element to separate
        # different graticules.
        stop_df = pd.DataFrame(
            {
                "decl": [np.nan],
                "ra": [np.nan],
                "grat": None,
            }
        )
        graticule_list = []

        # Define points in each declination graticule
        for decl in np.arange(min_decl, max_decl + decl_space, decl_space):
            ra_steps = np.arange(0, 360 + step)
            this_graticule = pd.DataFrame(
                {
                    "grat": f"decl{decl}",
                    "decl": decl,
                    "ra": ra_steps,
                }
            )

            graticule_list.append(this_graticule)
            graticule_list.append(stop_df)

        # Define points in each R.A. graticule
        for ra in np.arange(min_ra, max_ra + step, ra_space):  # pylint: disable=C0103
            decl_steps = np.arange(min_decl, max_decl + step, step)
            this_graticule = pd.DataFrame(
                {
                    "grat": f"ra{ra}",
                    "decl": decl_steps,
                    "ra": ra,
                }
            )

            graticule_list.append(this_graticule)
            graticule_list.append(stop_df)

        graticule_points = bokeh.models.ColumnDataSource(
            pd.concat(graticule_list), name=name
        )

        self.connect_controls(graticule_points)
        return graticule_points

    def _make_horizon_graticule_points(
        self,
        min_alt=0,
        max_alt=80,
        alt_space=20,
        min_az=0,
        max_az=360,
        az_space=30,
        step=1,
        name="hz_gratidule_ds",
    ):
        """Create points that define graticules

        Parameters
        ----------
        min_alt : int, optional
            Alt. of minimum az graticule
            and ends of alt graticules (deg),
            by default 0
        max_alt : int, optional
            Alt. of maximum az graticulas
            and ends of alt graticules (deg),
            by default 80
        alt_space : int, optional
            Spacing of alt. graticules (deg), by default 20
        min_az : int, optional
            Az of first azimuth graticule (deg), by default 0
        max_ra : int, optional
            Az of last azimuth graticule (deg), by default 360
        az_space : int, optional
            Spacing of azimuth graticules (deg), by default 30
        step : int, optional
            Spacing of points along graticules, by default 1

        Returns
        -------
        graticule_points : `bokeh.models.ColumnDataSource`
            Bokeh data sources defining points in graticules.
        """

        # Bohek puts gaps in lines when there are NaNs in the
        # data frame. We will be platting many graticules, and we do not
        # want them connected, so define a "stop" element to separate
        # different graticules.
        stop_df = pd.DataFrame(
            {
                "decl": [np.nan],
                "ra": [np.nan],
                "alt": [np.nan],
                "az": [np.nan],
                "grat": None,
            }
        )
        graticule_list = []

        # Define points in each alt graticule
        for alt in np.arange(min_alt, max_alt + alt_space, alt_space):
            radius = 90 - alt
            start_bear = 0
            end_bear = 360 + step
            this_graticule = pd.DataFrame(
                self._make_horizon_circle_points(
                    90, 0, radius, start_bear, end_bear, step
                ).data
            )
            this_graticule["grat"] = f"Alt{alt}"

            graticule_list.append(this_graticule)
            graticule_list.append(stop_df)

        for az in np.arange(min_az, max_az + step, az_space):  # pylint: disable=C0103
            radius = 90
            this_graticule = pd.DataFrame(
                self._make_horizon_circle_points(
                    0, az + 90, radius, 0, 360 + step, step
                ).data
            )
            this_graticule.query(
                f"(alt > {min_alt}) and (alt <= {max_alt}) and (abs(az-{az}) < 1)",
                inplace=True,
            )
            this_graticule.sort_values(by="alt", inplace=True)
            this_graticule["grat"] = f"Az{az}"

            graticule_list.append(this_graticule)
            graticule_list.append(stop_df)

        graticule_points = bokeh.models.ColumnDataSource(
            pd.concat(graticule_list), name=name
        )
        self.connect_controls(graticule_points)
        return graticule_points

    def _make_circle_points(
        self,
        center_ra,
        center_decl,
        radius=90.0,
        start_bear=0,
        end_bear=360,
        step=1,
        name="circle",
    ):
        """Create points along a circle or arc on a sphere

        Parameters
        ----------
        center_ra : `float`
            R.A. of the center of the circle (deg.).
        center_decl : `float`
            Decl. of the center of the circle (deg.).
        radius : float, optional
            Radius of the circle (deg.), by default 90.0
        start_bear : int, optional
            Bearing (E. of N.) of the start of the circle (deg.), by default 0
        end_bear : int, optional
            Bearing (E. of N.) of the end of the circle (deg.), by default 360
        step : int, optional
            Spacing of the points along the circle (deg.), by default 1
        name : `str`, optional
            Name for the bokeh model

        Returns
        -------
        circle : `bokeh.models.ColumnDataSource`
            Bokeh data source for points in the circle.
        """
        ras = []
        decls = []

        bearing_angles = np.arange(start_bear, end_bear + step, step) * u.deg
        center_coords = SkyCoord(center_ra * u.deg, center_decl * u.deg)
        circle_coords = center_coords.directional_offset_by(
            bearing_angles, radius * u.deg
        )
        bearings = bearing_angles.value.tolist()
        ras = circle_coords.ra.deg.tolist()
        decls = circle_coords.dec.deg.tolist()

        circle = bokeh.models.ColumnDataSource(
            data={
                "bearing": bearings,
                "ra": ras,
                "decl": decls,
            },
            name=name,
        )

        self.connect_controls(circle)
        return circle

    def _make_horizon_circle_points(
        self,
        alt=ALMOST_90,
        az=0,  # pylint: disable=C0103
        radius=90.0,
        start_bear=0,
        end_bear=360,
        step=1,
        name="horizon_ds",
    ):
        """Define points in a circle with the center defined in horizon coords.

        Parameters
        ----------
        alt : `float`
            Altitude of circle center, by default 90.
        az : `float`
            Azimuth of circle center, by default 0.
        radius : `float`, optional
            Radius of the circle (deg.), by default 90.0
        start_bear : int, optional
            Bearing of the start of the circle (deg.), by default 0
        end_bear : int, optional
            Bearing of the end of the circle (deg.), by default 360
        step : int, optional
            Spacing of points along the circle., by default 1
        name : str, optional
            Name for the bokeh model.

        Returns
        -------
        circle : `bokeh.models.ColumnDataSource`
            Bokeh data source with points along the circle.
        """

        center = SkyCoord(
            alt=alt * u.deg,
            az=az * u.deg,
            obstime=astropy.time.Time(self.mjd, format="mjd", scale="utc"),
            location=self.location,
            frame="altaz",
        ).transform_to(astropy.coordinates.ICRS)

        center_ra = center.ra.deg
        center_decl = center.dec.deg

        eq_circle_points = self._make_circle_points(
            center_ra, center_decl, radius, start_bear, end_bear, step
        )
        ra = np.array(eq_circle_points.data["ra"])  # pylint: disable=C0103
        decl = np.array(eq_circle_points.data["decl"])
        alt, az = self.eq_to_horizon(ra, decl, degrees=True, cart=False)

        circle_data = dict(eq_circle_points.data)
        circle_data["alt"] = alt.tolist()
        circle_data["az"] = az.tolist()

        circle = bokeh.models.ColumnDataSource(data=circle_data, name=name)

        return circle

    def _make_points(self, points_data, ds_name="points_ds"):
        """Create a bokeh data source with locations of points on a sphere.

        Parameters
        ----------
        points_data : `Iterable` , `dict` , or `pandas.DataFrame`
            A source of data (to be passed to `pandas.DataFrame`)
            Must contain the following columns or keys:

            ``"ra"``
                The Right Ascension in degrees.
            ``"decl"``
                The declination in degrees.
        ds_name : `str`, optional
            Name for the bokeh model, defaults to points_ds

        Returns
        -------
        point : `bokeh.models.ColumnDataSource`
            A data source with point locations, including projected coords.
        """

        points_df = pd.DataFrame(points_data)

        try:
            name = points_df["name"]
        except KeyError:
            name = [str(uuid.uuid4())] * len(points_df)

        try:
            glyph_size = points_df["glyph_size"]
        except KeyError:
            glyph_size = [10] * len(points_df)

        data = {
            "name": name,
            "ra": points_df.ra.tolist(),
            "decl": points_df.decl.tolist(),
            "glyph_size": glyph_size,
        }

        # Add any additional data provided
        for column_name in points_df.columns:
            if column_name not in data:
                data[column_name] = points_df[column_name].to_list()

        points = bokeh.models.ColumnDataSource(data=data, name=ds_name)

        self.connect_controls(points)
        return points

    def _make_marker_data_source(
        self,
        ra=None,  # pylint: disable=C0103
        decl=None,
        name="anonymous",
        glyph_size=5,
        min_mjd=None,
        max_mjd=None,
        ds_name="marker_ds",
    ):
        """Add one or more circular marker(s) to the map.

        Parameters
        ----------
        ra : `float` or `Iterable`, optional
            R.A. of the marker (deg.), by default None
        decl : `float` or `Iterable`, optional
            Declination of the marker (deg.), by default None
        name : `str` or `Iterable` , optional
            Name for the thing marked, by default "anonymous"
        glyph_size : `int` or `Iterable`, optional
            Size of the marker, by default 5
        min_mjd : `float`
            Earlist time for which to show the marker
        max_mjd : `float`
            Latest time for which to show the marker
        ds_name : `str`, optional
            Name for the bokeh data source, defaults to marker_ds.

        Returns
        -------
        data_source : `bokeh.models.ColumnDataSource`
            A data source with marker locations, including projected coords.
        """
        ras = ra if isinstance(ra, Iterable) else [ra]
        decls = decl if isinstance(decl, Iterable) else [decl]
        if len(ras) > 0:
            glyph_sizes = (
                glyph_size
                if isinstance(glyph_size, Iterable)
                else [glyph_size] * len(ras)
            )
            names = [name] * len(ras) if isinstance(name, str) else name
        else:
            glyph_sizes = np.array([])
            names = np.array([])

        data = {
            "ra": ras,
            "decl": decls,
            "name": names,
            "glyph_size": glyph_sizes,
        }

        if (min_mjd is not None) or (max_mjd is not None):
            if len(ras) == 0:
                data["in_mjd_window"] = np.array([])
            else:
                data["in_mjd_window"] = [1] * len(ras)

        if min_mjd is not None:
            if not isinstance(min_mjd, Iterable):
                min_mjd = [min_mjd]
            if len(ras) < 1:
                min_mjd = np.array([])

            data["min_mjd"] = min_mjd

            for marker_index, this_min_mjd in enumerate(min_mjd):
                if self.mjd < this_min_mjd:
                    data["in_mjd_window"][marker_index] = 0

        if max_mjd is not None:
            if not isinstance(max_mjd, Iterable):
                max_mjd = [max_mjd]
            if len(ras) < 1:
                max_mjd = np.array([])

            data["max_mjd"] = max_mjd

            for marker_index, this_max_mjd in enumerate(max_mjd):
                if self.mjd > this_max_mjd:
                    data["in_mjd_window"][marker_index] = 0

        data_source = self._make_points(data, ds_name=ds_name)

        self.connect_controls(data_source)
        return data_source

    def add_sliders(self):
        """Add (already defined) sliders to the map."""

        # In this, the base class, there are no sliders
        # to add, so just initialize the dictionary.
        self.sliders = OrderedDict()

    def add_mjd_slider(self):
        """Add a slider to control the MJD."""
        if "mjd" not in self.sliders:
            self.sliders["mjd"] = bokeh.models.Slider(
                start=self.mjd - 1,
                end=self.mjd + 1,
                value=self.mjd,
                step=0.25 / (24 * 12),
                title="MJD",
                name="mjd_slider",
            )

            self.visible_slider_names.append("mjd")

    def _create_js_update_func(self, data_source):
        """Set the javascript update functions for each slider

        Parameters
        ----------
        data_source : `bokeh.models.ColumnDataSource`
            The bokeh data source to update.
        """
        update_func = bokeh.models.CustomJS(
            args=dict(
                data_source=data_source,
                center_alt_slider={"value": 90},
                center_az_slider={"value": 0},
                mjd_slider=self.sliders["mjd"],
                lat=self.location.lat.deg,
                lon=self.location.lon.deg,
            ),
            code=self.update_js,
        )
        self.update_functions.append(update_func)
        return update_func

    def set_js_update_func(self, data_source):
        """Set the javascript update functions for each slider

        Parameters
        ----------
        data_source : `bokeh.models.ColumnDataSource`
            The bokeh data source to update.
        """
        update_func = self._create_js_update_func(data_source)

        self.force_update_time.js_on_change("text", update_func)

        for proj_slider_key in self.proj_slider_keys:
            try:
                self.sliders[proj_slider_key].js_on_change("value", update_func)
            except KeyError:
                pass

    def set_emit_update_func(self, data_source):
        """Set a javascript update functions for each slider to do nothing

        Parameters
        ----------
        data_source : `bokeh.models.ColumnDataSource`
            The bokeh data source to update.
        """
        update_func = bokeh.models.CustomJS(
            args={"data_source": data_source},
            code="data_source.change.emit()",
        )

        self.force_update_time.js_on_change("text", update_func)

        for proj_slider_key in self.proj_slider_keys:
            try:
                self.sliders[proj_slider_key].js_on_change("value", update_func)
            except KeyError:
                pass

    def show(self):
        """Show the map using bokeh only."""
        self.notebook_handle = bokeh.io.show(self.figure, notebook_handle=True)
        self.update()
        bokeh.io.push_notebook()

    def notebook_display(self):
        """Use panel to show the figure in within a notebook."""
        display(self.viewable)
        self.update()

    def add_healpix(
        self,
        data,
        cmap=None,
        nside=16,
        bound_step=1,
        ds_name="hpix_ds",
        name="hpix_renderer",
    ):
        """Add healpix values to the map

        Parameters
        ----------
        data : `numpy.ndarray`
            Healpixel values (RING pixel ordering)
        cmap : `dict`, optional
            With the following keys:

            ``"field"``
                The field to map to color (`str`)
            ``"transform"``
                The transform from value to color
                (`bokeh.models.mappers.ColorMapper`)

            by default None
        nside : `int`, optional
            Healpix nside to use for display, by default 16
        bound_step : `int`, optional
            number of boundary points for each side of each healpixel,
            by default 1
        ds_name : `str`, optional
            The bokeh name for the data source, defaults to hpix_ds.
        name : `str`, optional
            The bokeh name for the bokeh glyph, defaults ho hpix_glyph.

        Returns
        -------
        data_sounce : `bokeh.models.ColumnDataSource`
            The data source with the healpix values and bounds.
        cmap : `dict`
            With the following keys:

            ``"field"``
                The field to map to color (`str`)
            ``"transform"``
                The transform from value to color
                (`bokeh.models.mappers.ColorMapper`)

        hp_glyph : `bokeh.models.glyphs.Patches`
            The bokeh glyphs for the plotted patches.
        """
        if isinstance(data, bokeh.models.DataSource):
            data_source = data
        else:
            data_source = self.make_healpix_data_source(
                data, nside, bound_step, name=ds_name
            )

        self._add_projection_columns(data_source, nside)

        if cmap is None:
            cmap = make_zscale_linear_cmap(
                data_source.data["value"], cmap_name="hpix_color_transform"
            )
            self._hpix_cmap = cmap

        hpgr = self.plot.patches(
            xs=self.x_col,
            ys=self.y_col,
            fill_color=cmap,
            line_color=cmap,
            source=data_source,
            name=name,
        )

        hp_glyph = hpgr.glyph

        return data_source, cmap, hp_glyph

    def add_graticules(self, graticule_kwargs=None, line_kwargs=None):
        """Add graticules to the map

        Parameters
        ----------
        graticule_kwargs : dict, optional
            Keywords to be passed to ``SphereMap.make_graticule_points``,
            by default {}
        line_kwargs : dict, optional
            Keywords to be passed to ``bokeh.plotting.figure.line``,
            by default {}

        Returns
        -------
        graticules : ` `bokeh.models.ColumnDataSource`
            The bokeh data source with points defining the graticules.
        """
        if graticule_kwargs is None:
            graticule_kwargs = {}

        if line_kwargs is None:
            line_kwargs = {}

        graticule_points = self._make_graticule_points(**graticule_kwargs)
        kwargs = deepcopy(self.default_graticule_line_kwargs)
        kwargs.update(line_kwargs)
        self.plot.line(
            x=self.proj_transform("x", graticule_points),
            y=self.proj_transform("y", graticule_points),
            source=graticule_points,
            **kwargs,
        )
        return graticule_points

    def add_horizon_graticules(self, graticule_kwargs=None, line_kwargs=None):
        """Add graticules to the map

        Parameters
        ----------
        graticule_kwargs : dict, optional
            Keywords to be passed to ``SphereMap.make_graticule_points``,
            by default {}
        line_kwargs : dict, optional
            Keywords to be passed to ``bokeh.plotting.figure.line``,
            by default {}

        Returns
        -------
        graticules : ` `bokeh.models.ColumnDataSource`
            The bokeh data source with points defining the graticules.
        """
        if graticule_kwargs is None:
            graticule_kwargs = {}

        if line_kwargs is None:
            line_kwargs = {}

        graticule_points = self._make_horizon_graticule_points(**graticule_kwargs)
        kwargs = deepcopy(self.default_graticule_line_kwargs)
        kwargs.update(line_kwargs)
        self.plot.line(
            x=self.proj_transform("x", graticule_points),
            y=self.proj_transform("y", graticule_points),
            source=graticule_points,
            **kwargs,
        )
        return graticule_points

    def add_circle(self, center_ra, center_decl, circle_kwargs=None, line_kwargs=None):
        """Draw a circle on the map.

        Parameters
        ----------
        center_ra : `float`
            R.A. of the center of the circle (deg.)
        center_decl : `float`
            Decl. of the center of the circle (deg.)
        circle_kwargs : dict, optional
            Keywords to be passed to ``SphereMap.make_circle_points``,
            by default {}
        line_kwargs : dict, optional
            Keywords to be passed to ``bokeh.plotting.figure.line``,
            by default {}

        Returns
        -------
        circle_points : `bokeh.models.ColumnDataSource`
            The bokeh data source with points defining the circle.
        """
        if circle_kwargs is None:
            circle_kwargs = {}

        if line_kwargs is None:
            line_kwargs = {}

        if "name" not in line_kwargs:
            line_kwargs["name"] = "sphere_circle_glyph"

        circle_points = self._make_circle_points(
            center_ra, center_decl, **circle_kwargs
        )
        self.plot.line(
            x=self.proj_transform("x", circle_points),
            y=self.proj_transform("y", circle_points),
            source=circle_points,
            **line_kwargs,
        )
        return circle_points

    def add_horizon(
        self,
        zd=ALMOST_90,  # pylint: disable=C0103
        data_source=None,
        circle_kwargs=None,
        line_kwargs=None,
    ):
        """Add a circle parallel to the horizon.

        Parameters
        ----------
        zd : int, optional
            Zenith distance of the circle (deg), by default (almost) 90
        data_source : `bokeh.models.ColumnDataSource`, optional
            Bokeh data source for points on the circle,
            None if the should be generated.
            By default, None
        circle_kwargs : dict, optional
            Keywords to be passed to ``SphereMap.make_circle_points``,
            by default {}
        line_kwargs : dict, optional
            Keywords to be passed to ``bokeh.plotting.figure.line``,
            by default {}

        Returns
        -------
        circle_points : `bokeh.models.ColumnDataSource`
            The bokeh data source with points defining the circle.
        """
        if circle_kwargs is None:
            circle_kwargs = {}

        if line_kwargs is None:
            line_kwargs = {}

        if data_source is None:
            circle_points = self._make_horizon_circle_points(
                90, 0, radius=zd, **circle_kwargs
            )
            if "mjd" in self.sliders:
                self.set_emit_update_func(circle_points)
        else:
            circle_points = data_source

        kwargs = deepcopy(self.default_horizon_line_kwargs)
        kwargs.update(line_kwargs)
        self.plot.line(
            x=self.proj_transform("x", circle_points),
            y=self.proj_transform("y", circle_points),
            source=circle_points,
            **kwargs,
        )
        return circle_points

    def add_marker(
        self,
        ra=None,  # pylint: disable=C0103
        decl=None,
        name="anonymous",
        glyph_size=5,
        min_mjd=None,
        max_mjd=None,
        data_source=None,
        circle_kwargs=None,
    ):
        """Add one or more circular marker(s) to the map.

        Parameters
        ----------
        ra : `float` or `Iterable`, optional
            R.A. of the marker (deg.), by default None
        decl : `float` or `Iterable`, optional
            Declination of the marker (deg.), by default None
        name : `str` or `Iterable` , optional
            Name for the thing marked, by default "anonymous"
        glyph_size : `int` or `Iterable`, optional
            Size of the marker, by default 5
        min_mjd : `float` or `Iterable`, optional
            Earliest time for which to show the marker.
        max_mjd : `float` or `Iterable`, optional
            Latest time for which to show the marker.
        data_source : `bokeh.models.ColumnDataSource`, optional
            Data source for the marker, None if a new one is to be generated.
            By default, None
        circle_kwargs : dict, optional
            Keywords to be passed to ``bokeh.plotting.figure.circle``,
            by default {}

        Returns
        -------
        data_source : `bokeh.models.ColumnDataSource`
            A data source with marker locations, including projected coords.
        """
        if circle_kwargs is None:
            circle_kwargs = {}

        if "name" not in circle_kwargs:
            circle_kwargs["name"] = "circle_glyph"

        if data_source is None:
            data_source = self._make_marker_data_source(
                ra, decl, name, glyph_size, min_mjd, max_mjd
            )

        self.plot.circle(
            x=self.proj_transform("x", data_source),
            y=self.proj_transform("y", data_source),
            size="glyph_size",
            source=data_source,
            **circle_kwargs,
        )

        return data_source

    def _make_patches_data_source(self, patches_data, name=None):
        """Create a bokeh data source with locations of points on a sphere.

        All patches must have the same number of vertices.

        Parameters
        ----------
        patches_data : `pandas.DataFrame` or `dict`
            Must contain the following columns or keys:

            ``"ra"``
                The Right Ascension in degrees.
            ``"decl"``
                The declination in degrees.
        name : `str`, optional
            The name for the bokeh model.

        Returns
        -------
        data_source : `bokeh.models.ColumnDataSource`
            A data source with patch coordinates.
        """
        warn(
            "Create the data source directly with bokeh.models.ColumnDataSource",
            DeprecationWarning,
        )

        patches_df = pd.DataFrame(patches_data)
        ra = np.stack(patches_df.ra.values)  # pylint: disable=C0103
        decl = np.stack(patches_df.decl.values)

        wide_shape = ra.shape

        ra = ra.flatten()  # pylint: disable=C0103
        decl = decl.flatten()

        data = {
            "ra": ra.reshape(*wide_shape).tolist(),
            "decl": decl.reshape(*wide_shape).tolist(),
        }

        # Add any additional data provided
        for column_name in patches_df.columns:
            if column_name not in data:
                data[column_name] = patches_df[column_name].to_list()

        points = bokeh.models.ColumnDataSource(data=data, name=name)

        self.connect_controls(points)
        return points

    def add_patches(
        self,
        patches_data=None,
        data_source=None,
        patches_kwargs=None,
    ):
        """Add one or more patches to the map.

        Parameters
        ----------
        patches_data : `pandas.DataSource`
            Source of data.
        glyph_size : `int` or `Iterable`, optional
            Size of the marker, by default 5
        min_mjd : `float` or `Iterable`, optional
            Earliest time for which to show the marker.
        max_mjd : `float` or `Iterable`, optional
            Latest time for which to show the marker.
        data_source : `bokeh.models.ColumnDataSource`, optional
            Data source for the marker, None if a new one is to be generated.
            By default, None
        patches_kwargs : dict, optional
            Keywords to be passed to ``bokeh.plotting.figure.circle``,
            by default {}

        Returns
        -------
        data_source : `bokeh.models.ColumnDataSource`
            A data source with marker locations, including projected coords.
        """
        warn(
            "Use self.plot.patches directly instead",
            DeprecationWarning,
        )
        if patches_kwargs is None:
            patches_kwargs = {}

        if data_source is None:
            data_source = self._make_patches_data_source(patches_data)

        self.plot.patches(
            xs=self.proj_transform("x", data_source),
            ys=self.proj_transform("y", data_source),
            source=data_source,
            **patches_kwargs,
        )

        return data_source

    def add_stars(
        self,
        points_data=None,
        data_source=None,
        mag_limit_slider=False,
        star_kwargs=None,
    ):
        """Add stars to the map

        Parameters
        ----------
        points_data : `Iterable` , `dict` , or `pandas.DataFrame`
            A source of data (anything that can be passed to
            `pandas.DataFrame`)
            Must contain the following columns or keys:

            ``"ra"``
                The Right Ascension in degrees.
            ``"decl"``
                The declination in degrees.
        data_source : `bokeh.models.ColumnDataSource`, optional
            The bokeh data source to use (None to generate a new one).
            By default, None.
        mag_limit_slider : `bool` , optional
            Generate a slider limiting the magnitude of stars to plot,
            by default False
        star_kwargs : `dict` , optional
            _description_, by default {}

        Returns
        -------
        data_source : `bokeh.models.ColumnDataSource`
            The bokeh data source with points defining star locations.
        """
        if star_kwargs is None:
            star_kwargs = {}

        if "name" not in star_kwargs:
            star_kwargs["name"] = "stars_glyph"

        self.star_data = pd.DataFrame(points_data)
        if data_source is None:
            star_data_source = self._make_points(
                self.star_data, ds_name="bright_star_ds"
            )
        else:
            star_data_source = data_source

        self.plot.star(
            x=self.proj_transform("x", star_data_source),
            y=self.proj_transform("y", star_data_source),
            size="glyph_size",
            source=star_data_source,
            **star_kwargs,
        )

        if mag_limit_slider:
            mag_slider = bokeh.models.Slider(
                start=0,
                end=6.5,
                value=3,
                step=0.5,
                title="Magnitude limit for bright stars",
                name="mag_limit_slider",
            )
            mag_slider.on_change("value", self.limit_stars)

            self.sliders["mag_limit"] = mag_slider
            self.visible_slider_names.append("mag_limit")

        return star_data_source

    def limit_stars(self, attr, old_limit, mag_limit):  # pylint: disable=W0613
        """Apply a magnitude limit to mapped stars

        Parameters
        ----------
        attr : `str`
            Attribute of the slider to use (ignored)
        old_limit : `float`
            Old value for the magnitude limit (ignored)
        mag_limit : `float`
            Now value for the magnitude limit

        Note
        ----
        This method is intended to be called as a callback by bokeh.
        """
        star_data = self.star_data.query(f"Vmag < {mag_limit}").copy()
        star_data.loc[:, "glyph_size"] = (
            self.max_star_glyph_size
            - (self.max_star_glyph_size / mag_limit) * star_data["Vmag"]
        )
        stars = self._make_points(star_data)
        star_data_source = self.plot.select(name="bright_star_ds")
        star_data_source.data = dict(stars.data)

    def add_ecliptic(self, **kwargs):
        """Map the ecliptic.

        Returns
        -------
        points : `bokeh.models.ColumnDataSource`
            The bokeh data source with points on the ecliptic.
        """
        ecliptic_pole = SkyCoord(
            lon=0 * u.degree, lat=90 * u.degree, frame="geocentricmeanecliptic"
        ).icrs
        line_kwargs = deepcopy(self.default_ecliptic_line_kwargs)
        line_kwargs.update(kwargs)
        points = self.add_circle(
            ecliptic_pole.ra.deg, ecliptic_pole.dec.deg, line_kwargs=line_kwargs
        )
        return points

    def add_galactic_plane(self, **kwargs):
        """Map the galactic plane

        Returns
        -------
        points : `bokeh.models.ColumnDataSource`
            The bokeh data source with points on the galactic plane.
        """
        galactic_pole = SkyCoord(l=0 * u.degree, b=90 * u.degree, frame="galactic").icrs
        line_kwargs = deepcopy(self.default_galactic_plane_line_kwargs)
        line_kwargs.update(kwargs)
        points = self.add_circle(
            galactic_pole.ra.deg, galactic_pole.dec.deg, line_kwargs=line_kwargs
        )
        return points

    def decorate(self):
        """Add graticules, the ecliptic, and galactic plane to the map."""
        self.add_graticules()
        self.add_ecliptic()
        self.add_galactic_plane()

    def connect_controls(self, data_source):
        """Connect map controls to a data source so it updates.

        Parameters
        ----------
        data_source : `Bokeh.models.DataSource`
            A data source that needs to update automatically in response
            to changes to the map controls.
        """


class MovingSphereMap(SphereMap):
    """Parent class for maps that adjust based controls."""

    def connect_controls(self, data_source):
        """Connect map controls to a data source so it updates.

        Parameters
        ----------
        data_source : `Bokeh.models.DataSource`
            A data source that needs to update automatically in response
            to changes to the map controls.
        """
        if self.x_col in data_source.data:
            self.set_js_update_func(data_source)
        else:
            self.set_emit_update_func(data_source)


def make_zscale_linear_cmap(
    values, field_name="value", palette="Inferno256", cmap_name=None, **kwargs
):
    """Automatic linear scaling for images.

    Parameters
    ----------
    values : `numpy.array`
        Values for which to createa scale.
    field_name : `str`
        The bokeh column name for which to create a color scale
    palette : `str`
        The name of the bokeh palette to use.

    Returns
    -------
    cmap : `dict`
        The bokeh `DataSpec` dict to apply the color mapper.
    """
    zscale_interval = astropy.visualization.ZScaleInterval(**kwargs)
    if np.any(np.isfinite(values)):
        if len(np.unique(values)) == 1:
            value = np.nanmax(values)
            if value == 0.0:
                scale_limits = [-1.0, 0.1]
            else:
                scale_limits = [0.0, value * 1.1]
        if len(np.unique(values)) == 2:
            padding = (np.nanmax(values) - np.nanmin(values)) * 0.1
            scale_limits = [np.nanmin(values) - padding, np.nanmax(values) + padding]
        else:
            scale_limits = zscale_interval.get_limits(values)
    else:
        scale_limits = [0, 1]
    cmap = bokeh.transform.linear_cmap(
        field_name, palette, scale_limits[0], scale_limits[1]
    )
    cmap["transform"].name = cmap_name

    return cmap


def split_healpix_by_resolution(hp_map, nside_low=8, nside_high=None):
    """Split a healpix map into two at different resolutions.

    Parameters
    ----------

    hp_map : `numpy.ndarray`
        The original healpix map values.
    nside_low : `int`
        The healpix nside for the low-resolution map
    nside_high : `int`
        THe healpix nside for the high-resolution map

    Returns
    -------
    hp_map_high : `numpy.ndarray` or `healsparse.HealSparseMap`
        The high resolution healpixels. If an `ndarray`, other values
        have a value of `healpy.UNSEEN`.
    hp_map_low : `numpy.ndarray` or `healsparse.HealSparseMap`
        The low resolution healpixels. If an `ndarray`, other values
        have a value of `healpy.UNSEEN`.
    """
    if nside_high is None:
        nside_high = hp.npix2nside(hp_map.shape[0])
    else:
        hp_map = hp.ud_grade(hp_map, nside_high)

    hp_map_low = hp.ud_grade(hp_map, nside_low)
    hp_map_rev = hp.ud_grade(hp_map_low, nside_high)

    # 1 flags healpixels where the high resolution hpix value can be
    # completely reconstructed from the low resolution one. In this case,
    # that means all high resolution pixels within the low resolution
    # one have the same value.
    rev_matches_high = np.where(hp_map == hp_map_rev, 1, 0)
    use_low_low = np.where(hp.ud_grade(rev_matches_high, nside_low) == 1, 1, 0)
    use_low_high = np.where(hp.ud_grade(use_low_low, nside_high) == 1, 1, 0)

    try:
        hp_map_high = healsparse.HealSparseMap(
            nside_coverage=nside_high,
            healpix_map=hp.reorder(
                np.where(use_low_high == 0, hp_map, hp.UNSEEN), r2n=True
            ),
        )
        hp_map_low = healsparse.HealSparseMap(
            nside_coverage=nside_low,
            healpix_map=hp.reorder(
                np.where(use_low_low == 1, hp_map_low, hp.UNSEEN), r2n=True
            ),
        )
    except NameError:
        hp_map_high = np.where(use_low_high == 0, hp_map, hp.UNSEEN)
        hp_map_low = np.where(use_low_low == 1, hp_map_low, hp.UNSEEN)

    return hp_map_high, hp_map_low
