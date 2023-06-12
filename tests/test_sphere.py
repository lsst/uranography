import unittest

import numpy as np
import healpy as hp
from astropy.coordinates import SkyCoord

from numpy.random import default_rng

from uranography import sphere


def _random_point_on_sphere(rng):
    # Not slick or efficient, but it works
    mag = np.inf
    while mag > 1:
        vec = rng.uniform(-1, 1, 3)
        mag = np.linalg.norm(vec)

    vec = vec / mag
    return vec


class test_sphere(unittest.TestCase):
    def test_offset_sep_bear(self):
        rng = default_rng(6563)
        num_test_values = 100

        for i in range(num_test_values):
            vec1 = _random_point_on_sphere(rng)
            vec2 = _random_point_on_sphere(rng)
            ra, decl = hp.vec2ang(np.array([vec1, vec2]), lonlat=True)

            # Make sure at least one point is near each pole
            if i == 0:
                decl[0] = np.degrees(np.arccos(1.23e-13))
            elif i == 1:
                decl[0] = np.degrees(-1 * np.arccos(2.45e-14))

            coords = SkyCoord(ra, decl, unit="deg")
            separation = coords[0].separation(coords[1]).deg
            bearing = coords[0].position_angle(coords[1]).deg

            test_ra, test_decl = sphere.offset_sep_bear(
                ra[0], decl[0], separation, bearing, degrees=True
            )

            self.assertAlmostEqual(test_ra % 360, ra[1] % 360)
            self.assertAlmostEqual(test_decl % 360, decl[1] % 360)

            test_ra, test_decl = sphere.offset_sep_bear(
                np.radians(ra[0]),
                np.radians(decl[0]),
                np.radians(separation),
                np.radians(bearing),
                degrees=False,
            )

            self.assertAlmostEqual(np.degrees(test_ra) % 360, ra[1] % 360)
            self.assertAlmostEqual(np.degrees(test_decl) % 360, decl[1] % 360)

    def test_rotate_cart(self):
        rng = default_rng(4861)
        num_test_values = 100

        for i in range(num_test_values):
            vec0 = _random_point_on_sphere(rng)
            axis = _random_point_on_sphere(rng)
            angle = rng.uniform(-360, 360)
            vec1 = np.array(
                sphere.rotate_cart(
                    axis[0], axis[1], axis[2], angle, vec0[0], vec0[1], vec0[2]
                )
            )
            self.assertAlmostEqual(np.linalg.norm(vec1), 1.0)
            self.assertAlmostEqual(axis.dot(vec0), axis.dot(vec1))

        np.testing.assert_almost_equal(
            np.array(sphere.rotate_cart(0, 0, 1, 90, 1, 0, 0)), np.array((0, 1, 0))
        )

        np.testing.assert_almost_equal(
            np.array(sphere.rotate_cart(0, 0, 1, 180, 1, 0, 0)), np.array((-1, 0, 0))
        )

        np.testing.assert_almost_equal(
            np.array(sphere.rotate_cart(0, 0, 1, 90, 1, 0, 0)), np.array((0, 1, 0))
        )

        np.testing.assert_almost_equal(
            np.array(sphere.rotate_cart(1, 0, 0, 90, 0, 1, 0)), np.array((0, 0, 1))
        )


if __name__ == "__main__":
    unittest.main()
