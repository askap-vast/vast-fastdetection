class Error(Exception):
    """Base class of other error classes."""
    def __init__(self, message):
        Exception.__init__(self)
        self.message = message

    def __str__(self):
        return f"{self.message}"


class NoInputError(Error):
    def __init__(self, message="Input file(s) not available"):
        super().__init__(message)

class NoStdMapError(Error):
    def __init__(self, message="std map not available"):
        super().__init__(message)

class NoMapError(Error):
    def __init__(self, message="None of chisquare, peak or gaussian map is available"):
        super().__init__(message)
