"""Interactive sky map that works like an armillary sphere."""

from typing import Literal

import astropy.units as u
import bokeh
import numpy as np
import panel as pn
from astropy.coordinates import SkyCoord
from astropy.time import Time
from IPython.display import display

from .readjs import read_javascript
from .sphere import rotate_cart
from .spheremap import MovingSphereMap

CoordinateSystem = Literal["horizon", "equatorial", "both"]


class ArmillarySphere(MovingSphereMap):
    """Orthographic projection of the sky, presented like an armillary sphere.

    Parameters
    ----------
    plot : `bokeh.plotting.figure`, optional
        Figure to which to add the map, by default None
    mjd : `float`, optional
        The Modified Julian Date
    location: `EarthLocation` or `str`, optional
        The location of the observatory, defaults to lsst.
    """

    x_col = "x_orth"
    y_col = "y_orth"
    update_js_fnames = ("coord_utils.js", "orthographic.js")
    update_js_command = "updateOrthoData()"
    transform_js_fnames = ("coord_utils.js", "orthographic.js")
    transform_js_call = "return orthoTransform()"
    proj_slider_keys = ["alt", "az", "mjd", "up"]
    default_title = "Orthographic projection"
    default_coordinates: CoordinateSystem = "horizon"

    def to_orth_zenith(self, hpx, hpy, hpz):
        """Convert healpy vector coordinates to orthographic coordinates

        Parameters
        ----------
        hpx : `numpy.ndarray`
            Healpy vector x coordinates
            x=1, y=0, z=0 corresponds to R.A.=0 deg, Decl=0 deg.
            x=-1, y=0, z=0 corresponds to R.A.=180 deg, Decl=0 deg.
        hpy : `numpy.ndarray`
            Healpy vector y coordinates
            x=0, y=1, z=0 corresponds to R.A.=90 deg, Decl=0 deg.
            x=0, y=-1, z=0 corresponds to R.A.=270 deg, Decl=0 deg.
        hpz : `numpy.ndarray`
            Healpy vector z coordinates
            x=0, y=0, z=1 corresponds to Decl=90 deg.
            x=0, y=0, z=-1 corresponds to Decl=-90 deg.

        Returns
        -------
        x : `numpy.ndarray`
            Orthographic x coordinate (positive to the right)
        y : `numpy.ndarray`
            Orthographic y coordinate (positive up)
        z : `numpy.ndarray`
            Orthographic z coordinate (positive toward the viewer)
        """
        x1, y1, z1 = rotate_cart(0, 0, 1, -90, hpx, hpy, hpz)  # pylint: disable=C0103
        x2, y2, z2 = rotate_cart(1, 0, 0, self.location.lat.deg + 90, x1, y1, z1)  # pylint: disable=C0103

        npole_x1, npole_y1, npole_z1 = rotate_cart(0, 0, 1, -90, 0, 0, 1)
        npole_x2, npole_y2, npole_z2 = rotate_cart(
            1, 0, 0, self.location.lat.deg + 90, npole_x1, npole_y1, npole_z1
        )
        x3, y3, z3 = rotate_cart(npole_x2, npole_y2, npole_z2, -self.lst, x2, y2, z2)  # pylint: disable=C0103

        # x3, y3, z3 have the center right, now rotate it so that north is "up"
        # the 2-3 transform is about the axis through the n pole, so
        # the n pole is the same in 3 an 2.

        # Find the direction of the north pole, angle form +x axis toward
        # +y axis
        npole_x3, npole_y3 = npole_x2, npole_y2
        orient = np.degrees(np.arctan2(npole_y3, npole_x3))

        # To the n pole on the y axis, we must rotate it the rest of the 90 deg
        x4, y4, z4 = rotate_cart(0, 0, 1, 90 - orient, x3, y3, z3)  # pylint: disable=C0103

        # In astronomy, we are looking out of the sphere from the center to the
        # back (which naturally results in west to the right).
        # Positive z is out of the screen behind us, and we are at the center,
        # so to visible part is when z is negative (coords[2]<=0).
        # So, set the points with positive z to NaN so they are
        # not shown, because they are behind the observer.

        # Use np.finfo(z3[0]).resolution instead of exactly 0, because the
        # assorted trig operations result in values slightly above or below
        # 0 when the horizon is in principle exactly 0, and this gives an
        # irregularly dotted/dashed appearance to the horizon if
        # a cutoff of exactly 0 is used.

        try:
            orth_invisible = z4 > np.finfo(z4.dtype).resolution  # pylint: disable=E1101
        except ValueError:
            # If z4 is somehow an integer
            orth_invisible = z4 > 0

        x4[orth_invisible] = np.nan
        y4[orth_invisible] = np.nan
        z4[orth_invisible] = np.nan

        return x4, y4, z4

    def _add_projection_columns(self, hpix, nside=None, projector=None):
        """Adds pre-calculated projection columns for this projection."""
        hp_vec_cols = ["x_hp", "y_hp", "z_hp"]
        corners = hpix.to_df().loc[:, ["hpid"] + hp_vec_cols].explode(column=hp_vec_cols, ignore_index=True)
        corners["x_orth"], corners["y_orth"], corners["z_orth"] = self.to_orth_zenith(
            corners["x_hp"], corners["y_hp"], corners["z_hp"]
        )

        hpix_data = corners.groupby("hpid").agg(lambda x: x.tolist())

        for column in ["x_orth", "y_orth", "z_orth"]:
            hpix.add(hpix_data[column].tolist(), name=column)

        return hpix

    def set_js_update_func(self, data_source):
        """Set the javascript update functions for each slider

        Parameters
        ----------
        data_source : `bokeh.models.ColumnDataSource`
            The bokeh data source to update.
        """

        update_func = bokeh.models.CustomJS(
            args=dict(
                data_source=data_source,
                center_alt_slider=self.sliders["alt"],
                center_az_slider=self.sliders["az"],
                mjd_slider=self.sliders["mjd"],
                up_selector=self.sliders["up"],
                lat=self.location.lat.deg,
                lon=self.location.lon.deg,
            ),
            code=self.update_js,
        )

        self.update_functions.append(update_func)
        self.force_update_time.js_on_change("text", update_func)

        for proj_slider_key in self.proj_slider_keys:
            try:
                self.sliders[proj_slider_key].js_on_change("value", update_func)
            except KeyError:
                pass

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

        try:
            up_selector = self.sliders["up"]
        except KeyError:
            up_selector = {"value": "zenith is up"}

        js_code = "\n".join(read_javascript(fn) for fn in self.transform_js_fnames)
        js_code = "\n".join([js_code, self.transform_js_call])
        js_transform = bokeh.models.CustomJSTransform(
            args=dict(
                data_source=data_source,
                center_alt_slider=center_alt_slider,
                center_az_slider=center_az_slider,
                mjd_slider=mjd_slider,
                up_selector=up_selector,
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

    def add_sliders(self, center_alt=90, center_az=180):
        """Add (already defined) sliders to the map."""
        super().add_sliders()
        self.sliders["alt"] = bokeh.models.Slider(
            start=-90,
            end=90,
            value=center_alt,
            step=1,
            width=None,
            title="Altitude of map center at date & time from site",
        )
        self.sliders["az"] = bokeh.models.Slider(
            start=-360,
            end=360,
            value=center_az,
            step=1,
            width=None,
            title="Azimuth of map center at date & time from site",
        )
        self.sliders["mjd"] = bokeh.models.Slider(
            start=self.mjd - 1,
            end=self.mjd + 1,
            value=self.mjd,
            step=1.0 / (24 * 60),
            width=None,
            title="MJD",
        )
        self.sliders["up"] = bokeh.models.Select(
            value="zenith is up", options=["zenith is up", "north is up"], width=None
        )
        self.sliders["up"].visible = False

        self.visible_slider_names.append("alt")
        self.visible_slider_names.append("az")
        self.visible_slider_names.append("up")
        self.visible_slider_names.append("mjd")

        if self.default_coordinates == "equatorial":
            self.sliders["az"].visible = False
            self.sliders["alt"].visible = False
            self.sliders["up"].value = "north is up"

        if self.default_coordinates in ("equatorial", "both"):
            self.add_eq_sliders()

    def add_eq_sliders(self):
        """Add sliders to control the central RA and Decl of the map."""

        # Get RA and Decl positions corresponding to the current
        # alt/az/mjd
        eq_coords = SkyCoord(
            alt=self.sliders["alt"].value * u.deg,
            az=self.sliders["az"].value * u.deg,
            obstime=Time(self.sliders["mjd"].value, format="mjd"),
            location=self.location,
            frame="altaz",
        ).transform_to("icrs")
        initial_ra = np.round(eq_coords.ra.deg)
        initial_decl = np.round(eq_coords.dec.deg)

        self.sliders["ra"] = bokeh.models.Slider(
            start=0, end=360, value=initial_ra, step=1, title="R.A. of map center", width=None
        )

        self.sliders["decl"] = bokeh.models.Slider(
            start=-90, end=90, value=initial_decl, step=1, title="Declination of map center", width=None
        )

        # If sliders for both coords are visible,
        # default to showing a choice for up.
        if self.sliders["alt"].visible:
            self.sliders["up"].visible = True
        else:
            self.sliders["up"].value = "north is up"

        # bokeh js object to track whether we are in the middle
        # of an update, used to prevent multiple redundant callbacks.
        update_guard = bokeh.models.Div(text="false", visible=False)

        coord_update_js = read_javascript("match_eq_horizon_sliders.js")

        match_horizon_to_eq_callback = bokeh.models.CustomJS(
            args=dict(
                alt_slider=self.sliders["alt"],
                az_slider=self.sliders["az"],
                ra_slider=self.sliders["ra"],
                decl_slider=self.sliders["decl"],
                mjd_slider=self.sliders["mjd"],
                lat_deg=self.location.lat.deg,
                lon_deg=self.location.lon.deg,
                update_guard=update_guard,
                coord_to_update="horizon",
            ),
            code=coord_update_js,
        )
        for coord_name in ("ra", "decl"):
            self.sliders[coord_name].js_on_change("value", match_horizon_to_eq_callback)

        match_eq_to_horizon_callback = bokeh.models.CustomJS(
            args=dict(
                alt_slider=self.sliders["alt"],
                az_slider=self.sliders["az"],
                ra_slider=self.sliders["ra"],
                decl_slider=self.sliders["decl"],
                mjd_slider=self.sliders["mjd"],
                lat_deg=self.location.lat.deg,
                lon_deg=self.location.lon.deg,
                update_guard=update_guard,
                coord_to_update="eq",
            ),
            code=coord_update_js,
        )
        for coord_name in ("alt", "az", "mjd"):
            self.sliders[coord_name].js_on_change("value", match_eq_to_horizon_callback)

        self.visible_slider_names.append("ra")
        self.visible_slider_names.append("decl")

    def notebook_display(self):
        """Use panel to show the figure in within a notebook."""
        template = pn.Template(
            """
                {% extends base %}
                {% block postamble %}
                <script src="https://unpkg.com/gpu.js@latest/dist/gpu-browser.min.js"></script>
                {% endblock %}
            """
        )
        template.add_panel("main", self.viewable)
        display(template)
        self.update()
