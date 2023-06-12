import unittest
import healpy as hp
import numpy as np
import bokeh
from astropy.coordinates import SkyCoord, AltAz
from astropy.time import Time
import astropy.units as u
from uranography.api import ArmillarySphere
from .helpers import exercise_map_class, RNG, TEST_MJD, make_simple_map


class TestArmillarySphere(unittest.TestCase):
    def test_init(self):
        ArmillarySphere(mjd=TEST_MJD, location="Cerro Pachon")
        ArmillarySphere(mjd=TEST_MJD, location="APO")

    def test_to_orth_zenith(self):
        # make test data
        nside = 8
        npts = hp.nside2npix(nside)
        hpids = np.arange(npts)
        hp_x, hp_y, hp_z = hp.pix2vec(nside, hpids)

        test_map = ArmillarySphere(mjd=TEST_MJD, location="Cerro Pachon")
        x, y, z = test_map.to_orth_zenith(hp_x, hp_y, hp_z)

        ra, decl = hp.vec2ang(np.stack([hp_x, hp_y, hp_z]).T, lonlat=True)
        eq_coords = SkyCoord(ra * u.deg, decl * u.deg, frame="icrs")
        horizon_frame = AltAz(
            obstime=Time(test_map.mjd, format="mjd", scale="utc"),
            location=test_map.location,
            pressure=0.0,
        )
        hz_coords = eq_coords.transform_to(horizon_frame)
        ref_z = -np.sin(hz_coords.alt.rad)
        ref_x = np.sin(-1 * hz_coords.az.rad) * np.cos(hz_coords.alt.rad)
        ref_y = np.cos(hz_coords.az.rad) * np.cos(hz_coords.alt.rad)

        invisible = ref_z > 0
        ref_x[invisible] = np.nan
        ref_y[invisible] = np.nan
        ref_z[invisible] = np.nan

        # The ArmillarySphere class ignores distortions that
        # astropy does not, se we don't expect them to be
        # the same to high precision.
        np.testing.assert_allclose(x, ref_x, atol=0.01)
        np.testing.assert_allclose(y, ref_y, atol=0.01)
        np.testing.assert_allclose(z, ref_z, atol=0.01)

    def test__add_projection_columns(self):
        test_map = ArmillarySphere(mjd=TEST_MJD, location="Cerro Pachon")
        nside = 8
        hpvalues = RNG.random(hp.nside2npix(nside))

        # make_healpix_data_source calls _add_projection_columns
        data_source = test_map.make_healpix_data_source(hpvalues, nside=nside)
        self.assertIsInstance(data_source, bokeh.models.ColumnDataSource)

    def test_set_js_update_func(self):
        test_points, test_map = make_simple_map(ArmillarySphere)
        test_map.add_mjd_slider()
        test_map.set_js_update_func(test_points)

    # This test fails, possibly because it can't load gpujs in the test case.
    @unittest.skip("Known failure")
    def test_many(self):
        exercise_map_class(ArmillarySphere)
