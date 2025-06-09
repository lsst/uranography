"""Interactive planisphere (Lambert Azimuthal Equal Area projection)."""

import bokeh
import healpy as hp
import numpy as np

from .spheremap import SphereMap


class Planisphere(SphereMap):
    """Lambert azimuthal equal area projection of the sky,
    presented like a planisphere.

    Parameters
    ----------
    plot : `bokeh.plotting.figure`, optional
        Figure to which to add the map, by default None
    mjd : `float`, optional
        The Modified Julian Date
    location: `EarthLocation` or `str`, optional
        The location of the observatory, defaults to lsst.
    """

    x_col = "x_laea"
    y_col = "y_laea"
    update_js_fnames = ("coord_utils.js", "laea.js")
    update_js_command = "updateLAEAData()"
    transform_js_fnames = ("coord_utils.js", "laea.js")
    transform_js_call = "return laeaTransform()"
    default_title = "Planisphere"

    def __init__(self, plot=None, mjd=None, location="Cerro Pachon", laea_limit_mag=88.0):
        super().__init__(plot, mjd, location)
        self.laea_limit_mag = laea_limit_mag

    @property
    def laea_rot(self):
        """Return the `rot` tuple to be used in the Lambert EA projection

        Returns
        -------
        rot : `tuple` [`float`]
            The `rot` tuple to be passed to `healpy.projector.AzimuthalProj`.
        """
        rot = (0, -90, 0) if self.location.lat.deg < 0 else (0, 90, 180)
        return rot

    @property
    def laea_limit(self):
        """Return the lat. furthest from the center for the LAEA projection.

        Returns
        -------
        `limit` : `float`
            The maximum (or minimum) value for the latitude shown in the
            Lambert Azimuthal Equal Area plot.
        """
        limit = self.laea_limit_mag if self.location.lat.deg < 0 else -1 * self.laea_limit_mag
        return limit

    def _add_projection_columns(self, hpix, nside, projector=None):
        """Adds pre-calculated projection columns for this projection."""
        proj = hp.projector.AzimuthalProj(rot=self.laea_rot, lamb=True)
        proj.set_flip("astro")
        hpix = super()._add_projection_columns(hpix, nside, proj)

        # Healpixels at the opposite pole from the center behave badly, so
        # hide them.
        for hp_idx, decl in enumerate(hpix.data["center_decl"]):
            if (self.location.lat.deg < 0 and decl > self.laea_limit) or (
                self.location.lat.deg > 0 and decl < self.laea_limit
            ):
                for corner_idx in range(len(hpix.data["x_laea"][hp_idx])):
                    hpix.data["x_laea"][hp_idx][corner_idx] = np.nan
                    hpix.data["y_laea"][hp_idx][corner_idx] = np.nan

        return hpix

    def label_ra_graticules(self, graticule_points, **text_kwargs):
        """Add R.A. graticule labels.

        Parameters
        ----------
        graticule_points : `bokeh.models.ColumnDataSource`
            The data source for the graticules themselves
        **text_kwargs
            Additional parameters passed to `bokeh.figure.text`
        """

        graticules = graticule_points.to_df()
        ras = np.sort(graticules.loc[graticules["grat"].str.startswith("ra").astype("bool"), "ra"].unique())[
            :-1
        ]
        eps = np.finfo(np.float32).eps
        for ra in ras:
            text = f"{ra}\u00b0\n{int(ra/15)} hr"
            if self.laea_rot[1] < 0:
                sin_ang, cos_ang = np.sin(np.radians(ra)), np.cos(np.radians(ra))
                decl = 90 - eps
            else:
                sin_ang, cos_ang = np.sin(np.radians(-ra)), np.cos(np.radians(-ra))
                decl = -90 + eps

            if sin_ang > eps:
                horizontal_anchor = "right"
            elif sin_ang < -eps:
                horizontal_anchor = "left"
            else:
                horizontal_anchor = "center"

            if cos_ang < -eps:
                vertical_anchor = "top"
            elif cos_ang > eps:
                vertical_anchor = "bottom"
            else:
                vertical_anchor = "center"

            label_ds = bokeh.models.ColumnDataSource({"coords": [(ra, decl)], "text": [text]})
            self.plot.text(
                x=self.x_transform("coords"),
                y=self.y_transform("coords"),
                text="text",
                source=label_ds,
                anchor=f"{vertical_anchor}_{horizontal_anchor}",
                **text_kwargs,
            )

    def label_decl_graticules(self, graticule_points, **text_kwargs):
        """Add declination graticule labels.

        Parameters
        ----------
        graticule_points : `bokeh.models.ColumnDataSource`
            The data source for the graticules themselves
        **text_kwargs
            Additional parameters passed to `bokeh.figure.text`
        """

        graticules = graticule_points.to_df()
        ra = 135 if self.laea_rot[1] > 0 else 225
        decls = np.sort(
            graticules.loc[graticules["grat"].str.startswith("decl").astype("bool"), "decl"].unique()
        )
        label_ds = bokeh.models.ColumnDataSource(
            {"coords": [(ra, d) for d in decls], "text": [f"{d}\u00b0" for d in decls]}
        )
        self.plot.text(
            x=self.x_transform("coords"),
            y=self.y_transform("coords"),
            text="text",
            source=label_ds,
            anchor="center_center",
            **text_kwargs,
        )
