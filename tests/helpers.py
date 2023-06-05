from tempfile import TemporaryDirectory
from pathlib import Path
import numpy as np
import healpy as hp
import pandas as pd
import bokeh.io

TEST_MJD = 60300
RNG = np.random.default_rng(seed=6563)
TEST_STARS = {
    "name": [
        "Sirius",
        "Canopus",
        "Arcturus",
        "Alpha Centauri",
        "Vega",
        "Capella",
        "Rigel",
        "Procyon",
        "Achemar",
        "Betelgeuse",
    ],
    "ra": [101.3, 96.0, 213.9, 219.9, 279.2, 79.2, 78.6, 114.8, 24.4, 88.8],
    "decl": [-16.7, -52.7, 19.2, -60.8, 38.8, 46.0, -8.2, 5.2, -57.2, 7.4],
    "Vmag": [-1.46, -0.72, -0.04, -0.01, 0.03, 0.08, 0.12, 0.38, 0.46, 0.5],
}


def exercise_map_class(MapClass):
    test_map = MapClass(mjd=TEST_MJD, location="Cerro Pachon")
    test_map.add_mjd_slider()

    nside = 8
    hpdata = RNG.random(hp.nside2npix(nside))
    test_map.add_healpix(hpdata, nside=nside)
    test_map.decorate()
    test_map.add_horizon()
    test_map.add_marker(
        30, 10, name="Something", glyph_size=15, circle_kwargs={"color": "orange"}
    )

    stars = pd.DataFrame(TEST_STARS)
    test_map.add_stars(stars)

    with TemporaryDirectory() as test_dir:
        out_path = Path(test_dir)
        png_fname = out_path.joinpath("test_plot.png")
        bokeh.io.export_png(test_map.figure, filename=png_fname)


def make_simple_map(cls):
    test_map = cls(mjd=TEST_MJD)
    test_points = bokeh.models.ColumnDataSource(data=TEST_STARS)
    test_map.plot.asterisk(
        test_map.proj_transform("x", test_points),
        test_map.proj_transform("y", test_points),
        source=test_points,
    )
    return test_points, test_map
