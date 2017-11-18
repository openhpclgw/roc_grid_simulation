import itertools as it

class Coord(object):
    def __init__(self, i, j):
        self.i = i
        self.j = j
        self.__coord = (i, j)

    def __getitem__(self, key):
        return self.__coord[key]

    def __str__(self):
        return str(self.__coord)

    def distance(self, other):
        return abs(self.i-other.i) + abs(self.j-other.j)

    def is_neighbor(self, other):
        return self.distance(other) == 1

    def is_h(self, other):
        return abs(self.j-other.j) == 1

    # this should be enough to run the logic for now
    def __gt__(self, o):
        return self.i > o.i or (self.i == o.i and self.j > o.j)


class BoundingBox(object):
    # bbox : (left, top, width, height)
    def __init__(self, left, top, width, height):
        # just in case for backward compat
        self.tup = (left, top, width, height)

        self.left = left
        self.top = top
        self.width = width
        self.height = height

        self.coord = Coord(left, top)

    def __iter__(self):
        for i, j in it.product(range(self.top, self.top+self.height),
                               range(self.left, self.left+self.width)):
            yield i, j

    def __str__(self):
        return str(self.tup)

    def __repr__(self):
        return str(self)

    # intersection operator
    def __and__(self, other):
        small_j = self if self.left <= other.left else other
        big_j = self if small_j==other else other

        j_overlap = min(big_j.width,
                        small_j.left+small_j.width-big_j.left) # +:overlap

        small_i = self if self.top <= other.top else other
        big_i = self if small_i==other else other

        i_overlap = min(big_i.height,
                        small_i.top+small_i.height-big_i.top)
        print(self)
        print(other)
        print(i_overlap)
        print(j_overlap)

        if i_overlap>0 and j_overlap>0:
            return BoundingBox(big_j.left, big_i.top,
                               j_overlap, i_overlap)
        else:
            return None
    
    def lo(self, o):
        return BoundingBox(self.left+o, self.top,
                           self.width, self.height)
    def to(self, o):
        return BoundingBox(self.left, self.top+o,
                           self.width, self.height)

    def __getitem__(self, key):
        return self.tup[key]
