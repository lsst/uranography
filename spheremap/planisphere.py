"""Interactive planisphere (Lambert Azimuthal Equal Area projection)."""


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

    def __init__(
        self, plot=None, mjd=None, location="Cerro Pachon", laea_limit_mag=88.0
    ):
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
        limit = (
            self.laea_limit_mag
            if self.location.lat.deg < 0
            else -1 * self.laea_limit_mag
        )
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
                    hpix.data["x_laea"][hp_idx][corner_idx] = np.NaN
                    hpix.data["y_laea"][hp_idx][corner_idx] = np.NaN

        return hpix
