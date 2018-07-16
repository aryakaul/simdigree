class Dummy:

    """
    Class for use when given a pedigree to model.
    I create a 'dummy' for each individual with key information about that
    individual
    """

    def __init__(self, parents, founder, children=[]):
        self.parents = parents
        self.founder = founder
        self.children = children

    def __str__(self):
        return("Parents: %s\nFounder?: %s\n Children: %s\n" % (self.get_parents(), self.is_founder(), self.get_children()))

    def is_founder(self):
        return self.founder

    def get_parents(self):
        return self.parents

    def get_children(self):
        return self.children

    def add_child(self, child):
        copy = list(self.children)
        copy.append(child)
        self.children = copy

class DummyPair:

    """
    Class for use when given a pedigree to model.
    I create a 'DummyPair' for each married pair found in the original
    pedigree
    """

    def __init__(self, parents, children=[]):
        self.pair = parents
        self.children = children

    def __eq__(self, other):
        if self.pair == other.pair:
            return True
        return False

    def __str__(self):
        return("Pair: %s\tNo.Children: %s" % (self.get_pair(), self.get_num_children()))

    def get_num_children(self):
        return len(self.get_children())

    def get_children(self):
        return self.children

    def add_child(self, child):
        copy = list(self.children)
        copy.append(child)
        self.children = copy

    def get_pair(self):
        return self.pair
