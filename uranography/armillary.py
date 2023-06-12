"""Interactive sky map that works like an armillary sphere.
"""
import numpy as np
import panel as pn
import bokeh
from IPython.display import display

from .spheremap import MovingSphereMap
from .sphere import rotate_cart


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
    proj_slider_keys = ["alt", "az", "mjd"]
    default_title = "Armillary Sphere"

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
        x2, y2, z2 = rotate_cart(  # pylint: disable=C0103
            1, 0, 0, self.location.lat.deg + 90, x1, y1, z1
        )

        npole_x1, npole_y1, npole_z1 = rotate_cart(0, 0, 1, -90, 0, 0, 1)
        npole_x2, npole_y2, npole_z2 = rotate_cart(
            1, 0, 0, self.location.lat.deg + 90, npole_x1, npole_y1, npole_z1
        )
        x3, y3, z3 = rotate_cart(  # pylint: disable=C0103
            npole_x2, npole_y2, npole_z2, -self.lst, x2, y2, z2
        )

        # x3, y3, z3 have the center right, now rotate it so that north is "up"
        # the 2-3 transform is about the axis through the n pole, so
        # the n pole is the same in 3 an 2.

        # Find the direction of the north pole, angle form +x axis toward
        # +y axis
        npole_x3, npole_y3 = npole_x2, npole_y2
        orient = np.degrees(np.arctan2(npole_y3, npole_x3))

        # To the n pole on the y axis, we must rotate it the rest of the 90 deg
        x4, y4, z4 = rotate_cart(  # pylint: disable=C0103
            0, 0, 1, 90 - orient, x3, y3, z3
        )

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
        corners = (
            hpix.to_df()
            .loc[:, ["hpid"] + hp_vec_cols]
            .explode(column=hp_vec_cols, ignore_index=True)
        )
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

    def add_sliders(self, center_alt=90, center_az=180):
        """Add (already defined) sliders to the map."""
        super().add_sliders()
        self.sliders["alt"] = bokeh.models.Slider(
            start=-90,
            end=90,
            value=center_alt,
            step=np.pi / 180,
            title="center alt",
        )
        self.sliders["az"] = bokeh.models.Slider(
            start=-90, end=360, value=center_az, step=np.pi / 180, title="center Az"
        )
        self.sliders["mjd"] = bokeh.models.Slider(
            start=self.mjd - 1,
            end=self.mjd + 1,
            value=self.mjd,
            step=1.0 / (24 * 60),
            title="MJD",
        )

        self.visible_slider_names.append("alt")
        self.visible_slider_names.append("az")
        self.visible_slider_names.append("mjd")

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
