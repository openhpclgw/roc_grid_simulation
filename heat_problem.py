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
        
    def numerical_solve(self, num_steps):
        N = self.N
        # grid = self.gen_matrix()
        # grid2 = self.gen_matrix()
        grids = (self.gen_matrix(), self.gen_matrix())
        c = self.conductance

        for step in range(num_steps):
            for (i,j) in it.product(range(N), range(N)):
                if self.is_source((i,j)):
                    continue
                if self.is_sink((i,j)):
                    continue
                ing=step%2
                outg=1-ing
                tmp_sum = 0.
                if i-1 >= 0:
                    tmp_sum += grids[ing][i-1][j]
                if i+1 < N:
                    tmp_sum += grids[ing][i+1][j]
                if j-1 >= 0:
                    tmp_sum += grids[ing][i][j-1]
                if j+1 < N:
                    tmp_sum += grids[ing][i][j+1]
                grids[outg][i][j] = (tmp_sum)/4.

            abs_delta = 0.
            for (i,j) in it.product(range(N), range(N)):
                abs_delta += abs(grids[0][i][j]-grids[1][i][j])
            print(abs_delta)
            if abs_delta == 0:
                break

        return grids[num_steps%2]
