import unittest
from uranography.stars import load_bright_stars


class test_stars(unittest.TestCase):
    def test_load_stars(self):
        stars = load_bright_stars()
        self.assertGreater(len(stars), 9000)
        for column_name in ("name", "ra", "decl", "Vmag"):
            self.assertIn(column_name, stars.columns)
        self.assertGreaterEqual(stars.ra.min(), 0)
        self.assertLessEqual(stars.ra.min(), 1)
        self.assertGreaterEqual(stars.ra.max(), 359)
        self.assertLessEqual(stars.ra.max(), 360)
        self.assertGreaterEqual(stars.decl.min(), -90)
        self.assertLessEqual(stars.decl.min(), -88)
        self.assertGreaterEqual(stars.decl.max(), 89)
        self.assertLessEqual(stars.decl.max(), 90)
        self.assertGreaterEqual(stars.Vmag.min(), -2)
        self.assertLessEqual(stars.Vmag.max(), 10)


if __name__ == "__main__":
    unittest.main()
