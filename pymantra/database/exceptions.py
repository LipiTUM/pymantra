class IncorrectNodeType(ValueError):
    """
    Exception for incorrect node types. This can either mean that an invalid
    node type is encountered or that an invalid number of node types for a
    single node was found
    """
    def __init__(self, message: str):
        self.message = message


class IncorrectEdgeType(ValueError):
    """
    Exception for incorrect edge types. This can either mean that an invalid
    edge type is encountered or that an invalid number of edge types for a
    single edge was found
    """
    def __init__(self, message: str):
        self.message = message
