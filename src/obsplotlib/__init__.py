from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("obsplotlib")
except PackageNotFoundError:
    # If the package is not installed, don't add __version__
    pass


