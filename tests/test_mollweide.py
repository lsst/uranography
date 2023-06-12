import unittest
import healpy as hp
import bokeh
from uranography.api import MollweideMap
from .helpers import exercise_map_class, RNG, TEST_MJD


class TestMollweideMap(unittest.TestCase):
    def test_init(self):
        MollweideMap(mjd=TEST_MJD, location="Cerro Pachon")
        MollweideMap(mjd=TEST_MJD, location="APO")

    def test__add_projection_columns(self):
        test_map = MollweideMap(mjd=TEST_MJD, location="Cerro Pachon")
        nside = 8
        hpvalues = RNG.random(hp.nside2npix(nside))

        # make_healpix_data_source calls _add_projection_columns
        data_source = test_map.make_healpix_data_source(hpvalues, nside=nside)
        self.assertIsInstance(data_source, bokeh.models.ColumnDataSource)

    def test_many(self):
        exercise_map_class(MollweideMap)
