import unittest
from uranography.api import Planisphere
from .helpers import exercise_map_class, TEST_MJD


class TestPlanisphere(unittest.TestCase):
    def test_init(self):
        Planisphere(mjd=TEST_MJD, location="Cerro Pachon")
        Planisphere(mjd=TEST_MJD, location="APO")

    def test_laea_rot(self):
        south_map = Planisphere(mjd=TEST_MJD, location="Cerro Pachon")
        self.assertEqual(south_map.laea_rot, (0, -90, 0))

        north_map = Planisphere(mjd=TEST_MJD, location="APO")
        self.assertEqual(north_map.laea_rot, (0, 90, 180))

    def test_laea_limit(self):
        south_map = Planisphere(
            mjd=TEST_MJD, location="Cerro Pachon", laea_limit_mag=88.0
        )
        self.assertEqual(south_map.laea_limit, 88)

        north_map = Planisphere(mjd=TEST_MJD, location="APO", laea_limit_mag=88.0)
        self.assertEqual(north_map.laea_limit, -88)

    def test_many(self):
        exercise_map_class(Planisphere)
