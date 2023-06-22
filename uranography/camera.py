"""Tool for computing an outline of a camera footprint."""

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u


class CameraFootprintPerimeter(object):
    """Compute outline of an exposure with a camera."""

    def __init__(self, file_name):
        """Read the vertices of the polygon defining the camera footprint.

        Parameters
        ----------
        file_name : `str`
            The name of the file with the camera footprint vertices

        Returns
        -------
        footprint : `CameraFootprintPerimeter`
            The camera footprint instance
        """

        self.vertices = pd.read_csv(
            file_name, delim_whitespace=True, header=0, index_col=False
        )

        if "angle" not in self.vertices.columns:
            self.vertices["angle"] = np.degrees(
                np.arctan2(self.vertices.y, self.vertices.x)
            )

        if "r" not in self.vertices.columns:
            self.vertices["r"] = np.hypot(self.vertices.y, self.vertices.x)

    def single_eq_vertices(self, ra, decl, rotation=0):  # pylint: disable=C0103
        """Compute vertices for a single pair of equatorial coordinates

        Parameters
        ----------
        ra : `float`
            The R.A. (in degrees)
        decl : `float`
            The declination (in degrees)
        rotation : `float`
            The camera rotation (in degrees)

        Returns
        -------
        ra : `numpy.ndarray`
            An array of the R.A. of the vertices of the polygon surrounding
            the camera footprint (degrees).
        decl : `numpy.ndarray`
            An array of the declinations of the vertices of the polygon
            surrounding the camera footprint (degrees).
        """
        center = SkyCoord(ra, decl, unit="deg")

        # rotation matches the sense used by
        # rubin_sim.utils.camera_footprint.LsstCameraFootprint
        eq_vertices = center.directional_offset_by(
            (self.vertices.angle.values + rotation) * u.deg,
            self.vertices.r.values * u.deg,
        )
        ra = eq_vertices.ra.deg  # pylint: disable=C0103
        decl = eq_vertices.dec.deg
        return ra, decl

    def __call__(self, ra, decl, rotation=0):  # pylint: disable=C0103
        """Compute vertices for a single pair of equatorial coordinates

        Parameters
        ----------
        ra : `np.ndarray`
            The R.A. of pointings (in degrees)
        decl : `np.ndarray`
            The declination of pointings (in degrees)
        rotation : `float` or `np.ndarray`
            The camera rotation(s) (in degrees)

        Returns
        -------
        ra : `numpy.ndarray`
            An array of the R.A. of the vertices of the polygon surrounding
            the camera footprints (degrees).
        decl : `numpy.ndarray`
            An array of the declinations of the vertices of the polygon
            surrounding the camera footprints (degrees).
        """
        if np.isscalar(ra):
            return self.single_eq_vertices(ra, decl, rotation)

        if np.isscalar(rotation):
            rotation = np.full_like(ra, rotation)

        vertex_ras, vertex_decls = [], []
        for this_ra, this_decl, this_rotation in zip(ra, decl, rotation):
            this_vertex_ra, this_vertex_decl = self.single_eq_vertices(
                this_ra, this_decl, this_rotation
            )
            vertex_ras.append(this_vertex_ra)
            vertex_decls.append(this_vertex_decl)

        return vertex_ras, vertex_decls
