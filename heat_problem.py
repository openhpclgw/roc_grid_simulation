import numpy as np
import itertools as it

class HeatProblem(object):
    # bbox : (left, top, width, height)
    def __init__(self, N, source, sink, conductance, src_val=1., sink_val=0.):
        self.N = N
        self.source = source
        self.sink = sink
        self.conductance = conductance
        self.src_val = src_val
        self.sink_val = sink_val

    def __iter_bbox(self, bbox):
        for i,j in it.product(range(bbox[1], bbox[1]+bbox[3]),
                              range(bbox[0], bbox[0]+bbox[2])):
            yield (i,j)

    def source_iter(self):
        for idx in self.__iter_bbox(self.source):
            yield idx

    def sink_iter(self):
        for idx in self.__iter_bbox(self.sink):
            yield idx

    def __is_in_bbox(self, bbox, point):
        i, j = point
        left, top, width, height = bbox

        if j < left or j >= left+width:
            return False
        if i < top or i >= top+height:
            return False

        return True

    def is_source(self, p):
        return self.__is_in_bbox(self.source, p)

    def is_sink(self, p):
        return self.__is_in_bbox(self.sink, p)

    def gen_matrix(self):
        mat = np.zeros((self.N, self.N))
        for (i,j) in self.__iter_bbox(self.source):
            mat[i][j] = self.src_val
        for (i,j) in self.__iter_bbox(self.sink):
            mat[i][j] = self.sink_val

        return mat
        

