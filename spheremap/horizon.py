"""Interactive sky map in horizon (alt/az) coordinates."""

import bokeh
from .spheremap import MovingSphereMap


class HorizonMap(MovingSphereMap):
    """Horizon map of the sky, with zenith at center. Not equal area.

    Parameters
    ----------
    plot : `bokeh.plotting.figure`, optional
        Figure to which to add the map, by default None
    mjd : `float`, optional
        The Modified Julian Date
    location: `EarthLocation` or `str`, optional
        The location of the observatory, defaults to lsst.
    """

    x_col = "x_hz"
    y_col = "y_hz"
    update_js_fnames = ("coord_utils.js", "horizon.js")
    update_js_command = "updateHorizonData()"
    transform_js_fnames = ("coord_utils.js", "horizon.js")
    transform_js_call = "return horizonTransform()"
    proj_slider_keys = ["mjd"]
    default_title = "Horizon"

    def _add_projection_columns(self, hpix, nside=None, projector=None):
        """Adds pre-calculated projection columns for this projection."""
        coord_cols = ["ra", "decl"]
        n_hpix = len(hpix.data["ra"])
        corners_per_hpix = len(hpix.data["ra"][0])
        corners = (
            hpix.to_df()
            .loc[:, ["hpid"] + coord_cols]
            .explode(column=coord_cols, ignore_index=True)
        )
        assert len(corners) == n_hpix * corners_per_hpix

        x_hz, y_hz = self.eq_to_horizon(
            corners["ra"].astype(float), corners["decl"].astype(float)
        )

        hpix.add(x_hz.reshape(n_hpix, corners_per_hpix).tolist(), name="x_hz")
        hpix.add(y_hz.reshape(n_hpix, corners_per_hpix).tolist(), name="y_hz")

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
                center_alt_slider=90,
                center_az_slider=0,
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

    def add_sliders(self):
        """Add sliders to the map."""
        super().add_sliders()
        self.sliders["mjd"] = bokeh.models.Slider(
            start=self.mjd - 1,
            end=self.mjd + 1,
            value=self.mjd,
            step=1.0 / (24 * 60),
            title="MJD",
        )

        self.visible_slider_names.append("mjd")
