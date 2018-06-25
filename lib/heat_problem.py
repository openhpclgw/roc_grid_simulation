import numpy as np
import itertools as it
from common import *

# This is to be used for creating intermediate heat problems for
# virtualization
class PointHeatSource(object):
    def __init__(self, coord, val):
        self.coord = coord
        self.val = val


class HeatProblem(object):
    # bbox : (left, top, width, height)
    def __init__(self, N, sources, sinks, resistance,
                 src_val=1., sink_val=0.):

        def scalar_iter(obj):
            if isinstance(obj, list):
                for i in obj:
                    yield i
            else:
                yield obj

        self.N = N

        self.sources = [BoundingBox(*s) for s in scalar_iter(sources)]
        self.source_idxs = set()
        for s in self.sources:
            self.source_idxs |= {idx for idx in s}

        self.sinks = [BoundingBox(*s) for s in scalar_iter(sinks)]
        self.sink_idxs = set()
        for s in self.sinks:
            self.sink_idxs |= {idx for idx in s}

        self.resistance = resistance
        self.src_val = src_val
        self.sink_val = sink_val

    def get_r(self, n1, n2):
        return self.resistance(self.N, n1, n2)

    # def source_iter(self):
        # if isinstance(self.sources, list):
            # for s in self.sources:
                # yield s
        # else:
            # yield self.sources

    # def sink_iter(self):
        # if isinstance(self.sinks, list):
            # for s in self.sinks:
                # yield s
        # else:
            # yield self.sinks

    def __is_in_bbox(self, bbox, point):
        i, j = point
        left, top, width, height = bbox

        if j < left or j >= left+width:
            return False
        if i < top or i >= top+height:
            return False

        return True

    def is_source(self, p):
        return p in self.source_idxs

    def is_sink(self, p):
        return p in self.sink_idxs

    def gen_matrix(self, initial_values=None):
        if initial_values is not None:
            assert initial_values.shape == (self.N, self.N)
            return np.copy(initial_values)
        else:
            mat = np.zeros((self.N, self.N))
            for i, j in self.source_idxs:
                mat[i][j] = self.src_val
            for i, j in self.sink_idxs:
                mat[i][j] = self.sink_val
            return mat

    # this method is very naive and can use a lot of optimizations
    def numerical_solve(self, num_steps, epsilon=1e-9,
                        initial_values=None, report_progress=False):
        N = self.N
        grids = (self.gen_matrix(initial_values),
                 self.gen_matrix(initial_values))

        for step in range(num_steps):
            for i, j in it.product(range(N), range(N)):
                if self.is_source((i, j)):
                    continue
                if self.is_sink((i, j)):
                    continue
                ing = step % 2
                outg = 1-ing
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
            for i, j in it.product(range(N), range(N)):
                abs_delta += abs(grids[0][i][j]-grids[1][i][j])
            if report_progress:
                print('Step ', step, ' Delta ', abs_delta)
            if abs_delta < epsilon:
                break

        return grids[num_steps % 2]
