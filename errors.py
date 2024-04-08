class HexExistsError(Exception):
    def __init__(self, message):
        super(HexExistsError, self).__init__(message)
        
class IncorrectCoordinatesError(Exception):
    def __init__(self, message):
        super(IncorrectCoordinatesError, self).__init__(message)

class MistmatchError(Exception):
    def __init__(self, message):
        super(MistmatchError, self).__init__(message)
    