class Error(Exception):
    """Base class of other error classes."""
    def __init__(self, message):
        Exception.__init__(self)
        self.message = message

    def __str__(self):
        return f"{self.message}"


class NoInputError(Error):
    def __init__(self, message="No input file(s) available"):
        super().__init__(message)

