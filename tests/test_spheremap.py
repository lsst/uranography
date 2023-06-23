import unittest
from copy import deepcopy

import numpy as np
import pandas as pd
import healpy as hp
import bokeh.plotting
import panel as pn
import astropy.coordinates
from astropy.coordinates import EarthLocation, SkyCoord
import astropy.units as u
from astropy.time import Time

from uranography.api import SphereMap
from .helpers import TEST_MJD, TEST_STARS, RNG, exercise_map_class, make_simple_map

NUM_POINTS = 10


class TestSphereMap(unittest.TestCase):
    def test_init(self):
        SphereMap(plot=bokeh.plotting.figure())
        SphereMap(mjd=TEST_MJD)
        SphereMap(location="APO")
        SphereMap(location=EarthLocation.of_site("CTIO"))

    def test_figure(self):
        test_map = SphereMap(mjd=TEST_MJD)
        self.assertIsInstance(test_map.figure, bokeh.models.layouts.LayoutDOM)

    def test_viewable(self):
        test_map = SphereMap(mjd=TEST_MJD)
        self.assertIsInstance(test_map.viewable, pn.viewable.Viewable)

    @unittest.skip("Not implemented")
    def test_update(self):
        raise NotImplementedError()

    def test_lst(self):
        test_map = SphereMap(mjd=TEST_MJD)
        lst = test_map.lst
        orig_mjd = test_map.mjd
        test_map.lst = (lst + 123.4) % 360
        self.assertEqual(np.floor(orig_mjd), np.floor(test_map.mjd))
        test_map.mjd = orig_mjd
        self.assertAlmostEqual(lst, test_map.lst)

    def test_update_js(self):
        test_map = SphereMap(mjd=TEST_MJD)
        js = test_map.update_js
        self.assertIsInstance(js, str)

    def test_proj_transform(self):
        test_map = SphereMap(mjd=TEST_MJD)
        test_points = bokeh.models.ColumnDataSource(data=TEST_STARS)
        for coord in ("x", "y"):
            transform = test_map.proj_transform(coord, test_points)
            self.assertIsInstance(
                transform["transform"], bokeh.models.transforms.Transform
            )

    def test_eq_to_horizon(self):
        test_map = SphereMap(mjd=TEST_MJD, location="Cerro Pachon")
        test_points = bokeh.models.ColumnDataSource(data=TEST_STARS)

        ra, decl = test_points.data["ra"], test_points.data["decl"]

        test_alt_rad, test_az_rad = test_map.eq_to_horizon(
            np.radians(ra), np.radians(decl), degrees=False, cart=False
        )
        ref_coord_altaz = SkyCoord(
            alt=test_alt_rad * u.rad,
            az=test_az_rad * u.rad,
            obstime=Time(TEST_MJD, format="mjd", scale="utc"),
            location=EarthLocation.of_site("Cerro Pachon"),
            frame="altaz",
        )
        ref_coord_eq = ref_coord_altaz.transform_to(astropy.coordinates.ICRS)
        visible = test_alt_rad > 0

        np.testing.assert_allclose(ra, ref_coord_eq.ra.deg)
        np.testing.assert_allclose(decl, ref_coord_eq.dec.deg)

        test_alt_deg, test_az_deg = test_map.eq_to_horizon(
            ra, decl, degrees=True, cart=False
        )
        np.testing.assert_allclose(test_alt_deg, np.degrees(test_alt_rad))
        np.testing.assert_allclose(test_az_deg, np.degrees(test_az_rad))

        test_x, test_y = test_map.eq_to_horizon(ra, decl)
        self.assertTrue(np.all(np.isnan(test_x[~visible])))
        self.assertTrue(np.all(np.isnan(test_y[~visible])))
        np.testing.assert_allclose(
            np.hypot(test_x, test_y)[visible], np.pi / 2 - test_alt_rad[visible]
        )
        np.testing.assert_allclose(
            np.arctan2(-test_x, test_y)[visible] % (2 * np.pi), test_az_rad[visible]
        )

    def test_make_healpix_data_source(self):
        # Test with a plain healpy array
        nside = 8
        bound_step = 1
        npix = hp.nside2npix(nside)
        hpvalues = RNG.random(npix)

        # Make one value bad to verify that it does not get included
        hpvalues[3] = np.nan

        test_map = SphereMap(mjd=TEST_MJD, location="Cerro Pachon")

        data_source = test_map.make_healpix_data_source(hpvalues, nside=nside)
        self.assertIsInstance(data_source, bokeh.models.ColumnDataSource)

        self.assertEqual(len(data_source.data["value"]), npix - 1)
        self.assertEqual(len(data_source.data["ra"][5]), 4 * bound_step)

    def test_make_graticule_points(self):
        test_map = SphereMap(mjd=TEST_MJD)
        graticule_points = test_map._make_graticule_points()
        self.assertIsInstance(graticule_points, bokeh.models.ColumnDataSource)

    def test_make_horizon_graticule_points(self):
        test_map = SphereMap(mjd=TEST_MJD)
        graticule_points = test_map._make_horizon_graticule_points()
        self.assertIsInstance(graticule_points, bokeh.models.ColumnDataSource)

    def test_make_circle_points(self):
        test_map = SphereMap(mjd=TEST_MJD)
        test_ra, test_decl = 12.3, -2.2
        test_radius = 5
        circle_points = test_map._make_circle_points(test_ra, test_decl, test_radius)
        self.assertIsInstance(circle_points, bokeh.models.ColumnDataSource)

        center = SkyCoord(test_ra * u.deg, test_decl * u.deg, frame="icrs")
        for _, row in circle_points.to_df().iterrows():
            point = SkyCoord(row.ra * u.deg, row.decl * u.deg, frame="icrs")
            self.assertAlmostEqual(center.separation(point).deg, test_radius)
            self.assertAlmostEqual(
                center.position_angle(point).deg % 360, row.bearing % 360
            )

    def test_make_horizon_circle_points(self):
        test_map = SphereMap(mjd=TEST_MJD)
        test_alt, test_az = 88.0, -2.2
        test_radius = 65.0
        circle_points = test_map._make_horizon_circle_points(
            test_alt, test_az, test_radius
        )
        self.assertIsInstance(circle_points, bokeh.models.ColumnDataSource)

    def test_make_points(self):
        test_map = SphereMap(mjd=TEST_MJD)
        stars = deepcopy(TEST_STARS)

        test_map._make_points(pd.DataFrame(stars))

        points = test_map._make_points(stars)
        self.assertTrue("glyph_size" in points.data)

        del stars["name"]
        points = test_map._make_points(stars)
        self.assertTrue("name" in points.data)

        stars["foo"] = np.arange(len(stars["ra"]))
        stars["bar"] = np.arange(len(stars["ra"]))
        points = test_map._make_points(stars)
        self.assertTrue("foo" in points.data)
        self.assertTrue("bar" in points.data)

    def test_make_marker_data_source(self):
        test_map = SphereMap(mjd=TEST_MJD)
        ds = test_map._make_marker_data_source(
            ra=TEST_STARS["ra"], decl=TEST_STARS["decl"], name=TEST_STARS["name"]
        )
        self.assertIsInstance(ds, bokeh.models.ColumnDataSource)

    def test_add_sliders(self):
        test_map = SphereMap(mjd=TEST_MJD)
        test_map.add_sliders()

    def test_add_mjd_slider(self):
        test_map = SphereMap(mjd=TEST_MJD)
        test_map.add_mjd_slider()
        slider = test_map.sliders["mjd"]
        self.assertTrue(slider.start < TEST_MJD < slider.end)

    def test_set_js_update_func(self):
        test_points, test_map = make_simple_map(SphereMap)
        test_map.add_mjd_slider()
        test_map.set_js_update_func(test_points)

    def test_set_emit_update_func(self):
        test_points, test_map = make_simple_map(SphereMap)
        test_map.add_mjd_slider()
        test_map.set_emit_update_func(test_points)

    def test_show(self):
        test_map = make_simple_map(SphereMap)[1]
        test_map.show()

    def test_notebook_display(self):
        test_map = make_simple_map(SphereMap)[1]
        test_map.notebook_display()

    def test_add_healpix(self):
        # Test with a plain healpy array
        nside = 8
        npix = hp.nside2npix(nside)
        hpvalues = RNG.random(npix)

        # Make one value bad to verify that it does not get included
        hpvalues[3] = np.nan

        test_map = SphereMap(mjd=TEST_MJD, location="Cerro Pachon")
        data_source, cmap, hp_glyph = test_map.add_healpix(hpvalues, nside=nside)
        self.assertIsInstance(data_source, bokeh.models.ColumnDataSource)
        self.assertIsInstance(cmap["transform"], bokeh.models.mappers.ColorMapper)
        self.assertIsInstance(hp_glyph, bokeh.models.glyphs.Patches)

        test_map_2 = SphereMap(mjd=TEST_MJD, location="Cerro Pachon")

        test_map_2.add_healpix(hpvalues, cmap, nside=nside, bound_step=2)

    def test_add_graticules(self):
        test_map = SphereMap(mjd=TEST_MJD)
        test_map.add_graticules()

    def test_add_horizon_graticules(self):
        test_map = SphereMap(mjd=TEST_MJD)
        test_map.add_horizon_graticules()

    def test_add_circle(self):
        test_map = SphereMap(mjd=TEST_MJD)
        test_map.add_circle(TEST_STARS["ra"][0], TEST_STARS["decl"][1])

    def test_add_horizon(self):
        test_map = SphereMap(mjd=TEST_MJD)
        test_map.add_horizon()

    def test_add_marker(self):
        test_map = SphereMap(mjd=TEST_MJD)
        test_map.add_marker(
            ra=TEST_STARS["ra"], decl=TEST_STARS["decl"], name=TEST_STARS["name"]
        )

    def test_make_patches_data_source(self):
        patch1_ra = [30, 30, 50, 50]
        patch1_decl = [20, 50, 50, 20]
        patch2_ra = [130, 130, 150, 150]
        patch2_decl = [-20, -50, -50, -20]

        patch_data = {
            "ra": [patch1_ra, patch2_ra],
            "decl": [patch1_decl, patch2_decl],
            "foo": [1, 10],
        }

        test_map = SphereMap(mjd=TEST_MJD)
        patch_ds = test_map._make_patches_data_source(patch_data)

        self.assertTrue("foo" in patch_ds.data)

    def test_add_patches(self):
        patch1_ra = [30, 30, 50, 50]
        patch1_decl = [20, 50, 50, 20]
        patch2_ra = [130, 130, 150, 150]
        patch2_decl = [-20, -50, -50, -20]

        patch_data = {
            "ra": [patch1_ra, patch2_ra],
            "decl": [patch1_decl, patch2_decl],
            "foo": [1, 10],
        }

        test_map = SphereMap(mjd=TEST_MJD)
        test_map.add_patches(patch_data)

    def test_add_stars(self):
        test_map = SphereMap(mjd=TEST_MJD)
        test_map.add_stars(TEST_STARS, mag_limit_slider=False)

        test_map = SphereMap(mjd=TEST_MJD)
        test_map.add_stars(TEST_STARS, mag_limit_slider=True)

    def test_limit_stars(self):
        test_map = SphereMap(mjd=TEST_MJD)
        test_map.add_stars(TEST_STARS, mag_limit_slider=True)
        mag_limit = -0.5
        test_map.limit_stars(1, 99, mag_limit)

        star_data_source = test_map.plot.select(name="bright_star_ds")
        for mag in star_data_source.data["Vmag"]:
            self.assertTrue(mag < mag_limit)

    def test_add_ecliptic(self):
        test_map = SphereMap(mjd=TEST_MJD)
        test_map.add_ecliptic()

    def test_add_galactic_plane(self):
        test_map = SphereMap(mjd=TEST_MJD)
        test_map.add_galactic_plane()

    def test_decorate(self):
        test_map = SphereMap(mjd=TEST_MJD)
        test_map.decorate()

    def test_many(self):
        exercise_map_class(SphereMap)
