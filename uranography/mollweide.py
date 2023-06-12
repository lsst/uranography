"""Interactive sky map in a Mollweide projection."""

import numpy as np
import healpy as hp
from .spheremap import SphereMap


class MollweideMap(SphereMap):
    """Mollweide projection of the sky.

    Parameters
    ----------
    plot : `bokeh.plotting.figure`, optional
        Figure to which to add the map, by default None
    mjd : `float`, optional
        The Modified Julian Date
    location: `EarthLocation` or `str`, optional
        The location of the observatory, defaults to lsst.
    """

    x_col = "x_moll"
    y_col = "y_moll"
    update_js_fnames = ("coord_utils.js", "mollweide.js")
    update_js_command = "updateMollweideData()"
    transform_js_fnames = ("coord_utils.js", "mollweide.js")
    transform_js_call = "return mollweideTransform()"
    default_title = "Mollweide"

    def _add_projection_columns(self, hpix, nside, projector=None):
        """Adds pre-calculated projection columns for this projection."""
        proj = hp.projector.MollweideProj()
        proj.set_flip("astro")
        hpix = super()._add_projection_columns(hpix, nside, proj)

        resol = np.degrees(hp.nside2resol(nside))
        num_pix = len(hpix.data["ra"])
        for i in range(num_pix):
            center_ra = hpix.data["center_ra"][i]
            center_decl = hpix.data["center_decl"][i]

            # Skip any healpixes not close to the discontinuity
            ra_resol = resol / np.cos(np.radians(center_decl))
            if np.abs((center_ra + 180) % 360 - 180) < 2 * ra_resol:
                continue

            # If there are any x points with a different sign than the center,
            # set them to nan
            for j, x in enumerate(hpix.data["x_moll"][i]):  # pylint: disable=C0103
                center_ra_sign = 1 if center_ra % 360 < 180 else -1
                # Remember, we are looking out from the Earth, so
                # positive RA is left, negative right
                if np.sign(x) == center_ra_sign:
                    hpix.data["x_moll"][i][j] = np.NaN
                    hpix.data["y_moll"][i][j] = np.NaN

        return hpix
