"""Errors specific to ansys-mapdl-reader"""

class VersionError(ValueError):
    """Raised when MAPDL is the wrong version"""

    def __init__(self, msg='Invalid MAPDL version'):
        ValueError.__init__(self, msg)


class NoDistributedFiles(FileNotFoundError):
    """Unable to find any distributed result files """

    def __init__(self, msg='Unable to find any distributed result files.'):
        FileNotFoundError.__init__(self, msg)
