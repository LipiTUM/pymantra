class StatsTestException(ValueError):
    def __init__(self, message):
        self.message = message


class FoldChangeException(ValueError):
    def __init__(self, message):
        self.message = message


class CorrelationException(ValueError):
    def __init__(self, message, is_partial: bool = False):
        self.message = message
        self.is_partial = is_partial


class NodeTypeError(ValueError):
    def __init__(self, message):
        self.message = message


class GraphComponentError(ValueError):
    def __init__(self, message):
        self.message = message


class R2Warning(RuntimeWarning):
    pass
