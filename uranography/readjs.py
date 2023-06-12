"""Load javascript code for use in callbacks.

The main work of this submodule is locating the relevant file within the
uranography package.
"""

import importlib.resources


def read_javascript(fname):
    """Read javascript source code from the current package.

    Parameters
    ----------
    fname : `str`
        The name of the file from which to load js source code.

    Return
    ------
    js_code : `str`
        The loaded source code.
    """
    root_package = __package__.split(".", maxsplit=1)[0]

    try:
        js_path = importlib.resources.files(root_package).joinpath("js").joinpath(fname)
        with importlib.resources.as_file(js_path) as js_file_path:
            with open(js_file_path, "r", encoding="UTF-8") as js_io:
                js_code = js_io.read()
    except AttributeError as e:  # pylint: disable=C0103
        # If we are using an older version of importlib, we need to do
        # this instead:
        if e.args[0] != "module 'importlib.resources' has no attribute 'files'":
            raise e

        with importlib.resources.path(root_package, ".") as root_path:
            full_name = root_path.joinpath("js").joinpath(fname)
            with open(full_name, "r", encoding="UTF-8") as js_io:
                js_code = js_io.read()

    return js_code
