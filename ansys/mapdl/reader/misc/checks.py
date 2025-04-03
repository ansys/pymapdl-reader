ERROR_GRAPHICS_REQUIRED = (
    "Graphics are required for this method. Please install the ``graphics`` target "
    " to use this method. You can install it by running `pip install ansys-mapdl-reader[graphics]`"
    " or `pip install ansys-mapdl-reader[all]`."
)
"""Message to display when graphics are required for a method."""

ERROR_SCIPY_REQUIRED = (
    "SciPy is required for this method. Please install the ``all`` target "
    " to use this method. You can install it by running `pip install ansys-mapdl-reader[all]`."
)
"""Message to display when SciPy is required for a method."""

__GRAPHICS_AVAILABLE = None
"""Global variable to store the result of the graphics imports."""
__SCIPY_AVAILABLE = None
"""Global variable to store the result of the SciPy imports."""


def run_if_graphics_required():
    """Check if graphics are available."""
    global __GRAPHICS_AVAILABLE
    if __GRAPHICS_AVAILABLE is None:
        try:
            # Attempt to perform the imports
            import pyvista  # noqa F401
            import trame  # noqa F401
            import vtk  # noqa F401

            __GRAPHICS_AVAILABLE = True
        except (ModuleNotFoundError, ImportError):
            __GRAPHICS_AVAILABLE = False

    if __GRAPHICS_AVAILABLE is False:
        raise ImportError(ERROR_GRAPHICS_REQUIRED)


def run_if_scipy_required():
    """Check if graphics are available."""
    global __SCIPY_AVAILABLE
    if __SCIPY_AVAILABLE is None:
        try:
            # Attempt to perform the imports
            from scipy import sparse  # noqa F401
            from scipy.sparse import linalg  # noqa F401

            __SCIPY_AVAILABLE = True
        except (ModuleNotFoundError, ImportError):
            __SCIPY_AVAILABLE = False

    if __SCIPY_AVAILABLE is False:
        raise ImportError(ERROR_GRAPHICS_REQUIRED)


def graphics_required(method):
    """Decorate a method as requiring graphics.
    Parameters
    ----------
    method : callable
        Method to decorate.
    Returns
    -------
    callable
        Decorated method.
    """

    def wrapper(*args, **kwargs):
        run_if_graphics_required()
        return method(*args, **kwargs)

    return wrapper


def scipy_required(method):
    """Decorate a method as requiring graphics.
    Parameters
    ----------
    method : callable
        Method to decorate.
    Returns
    -------
    callable
        Decorated method.
    """

    def wrapper(*args, **kwargs):
        run_if_scipy_required()
        return method(*args, **kwargs)

    return wrapper


def are_graphics_available():
    """Check if graphics are available."""
    try:
        run_if_scipy_required()
        return True
    except ImportError:
        return False
