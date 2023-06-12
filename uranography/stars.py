"""Load the Yale Bright Star catalog into a pandas.DataFrame."""

from collections import OrderedDict
import os
import pandas as pd

BSC5_URL = "http://tdc-www.harvard.edu/catalogs/bsc5.dat.gz"


def load_bright_stars(fname=None):
    """Read the Yale Bright Star Catalog into a pandas.DataFrame.

    Parameters
    ----------
    fname : `str`, optional
        Name of file from which to load the catalog, by default None

    Returns
    -------
    bright_stars : `pandas.DataFrame`
        The catalog of bright stars.
    """
    if fname is None:
        try:
            fname = os.environ["BSC5_FNAME"]
        except KeyError:
            fname = BSC5_URL

    ybs_columns = OrderedDict(
        (
            ("HR", (0, 4)),
            ("name", (4, 14)),
            ("RA_hour", (75, 77)),
            ("RA_min", (77, 79)),
            ("RA_sec", (79, 83)),
            ("decl_sign", (83, 84)),
            ("decl_deg", (84, 86)),
            ("decl_min", (86, 88)),
            ("decl_sec", (88, 90)),
            ("Vmag", (102, 107)),
        )
    )

    compression = "gzip" if fname.endswith(".gz") else "infer"

    bright_stars = pd.read_fwf(
        fname,
        colspecs=[ybs_columns[k] for k in ybs_columns],
        names=[k for k in ybs_columns],
        compression=compression,
    )
    bright_stars["ra"] = (360 / 24) * (
        bright_stars.RA_hour + (bright_stars.RA_min + bright_stars.RA_sec / 60.0) / 60.0
    )
    bright_stars["decl"] = (
        bright_stars.decl_deg
        + (bright_stars.decl_min + bright_stars.decl_sec / 60.0) / 60.0
    )
    southern_stars = bright_stars.decl_sign == "-"
    bright_stars.loc[southern_stars, "decl"] = (
        -1 * bright_stars.loc[southern_stars, "decl"]
    )
    return bright_stars
