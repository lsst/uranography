import unittest
import healpy as hp
import bokeh
from uranography.api import HorizonMap
from .helpers import exercise_map_class, RNG, TEST_MJD, make_simple_map


class TestHorizonMap(unittest.TestCase):
    def test_init(self):
        HorizonMap(mjd=TEST_MJD, location="Cerro Pachon")
        HorizonMap(mjd=TEST_MJD, location="APO")

    def test__add_projection_columns(self):
        test_map = HorizonMap(mjd=TEST_MJD, location="Cerro Pachon")
        nside = 8
        hpvalues = RNG.random(hp.nside2npix(nside))

        # make_healpix_data_source calls add_projection_columns
        data_source = test_map.make_healpix_data_source(hpvalues, nside=nside)
        self.assertIsInstance(data_source, bokeh.models.ColumnDataSource)

    def test_set_js_update_func(self):
        test_points, test_map = make_simple_map(HorizonMap)
        test_map.add_mjd_slider()
        test_map.set_js_update_func(test_points)

    def test_add_sliders(self):
        test_map = HorizonMap(mjd=TEST_MJD)
        test_map.add_sliders()

    def test_many(self):
        exercise_map_class(HorizonMap)
