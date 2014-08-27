# Copyright (c) 2011, Mathias Laurin
#
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.


import numpy as np
import math


def vector_to_xyz(vector):
    return vector[0], vector[1], vector[2]


def normalize(vector):
    nvector = np.linalg.norm(vector)
    return vector / nvector


def distance(position1, position2):
    return np.sqrt(np.sum((position2 - position1) * (position2 - position1)))


def rad2deg(a):
    return ((a * 180.0 / math.pi) + 360.0) % 360.0


def deg2rad(a):
    return ((a * math.pi / 180.0))


def angle2d(vector1, vector2):
    vector1 = normalize(vector1)
    vector2 = normalize(vector2)
    return (math.atan2(vector2[1], vector2[0]) -
            math.atan2(vector1[1], vector1[0]))


def dangle2d(vector1, vector2):
    return rad2deg(angle2d(vector1, vector2))


def angle_OBSOLETE(vector1, vector2):
    vector1 = normalize(vector1)
    vector2 = normalize(vector2)
    dotproduct = np.dot(vector1, vector2)
    try:
        angle = math.acos(dotproduct)
    except ValueError:
        # Rounding error
        if dotproduct > 1.0:
            angle = 0.0
        elif dotproduct < 1.0:
            angle = math.pi
        else:
            raise
    return angle


def angle(vector1, vector2):
    return math.atan2(np.linalg.norm(np.cross(vector1, vector2)),
                      np.dot(vector1, vector2))


def dangle(vector1, vector2):
    ''' return angle between vector1 and vector2 in degrees '''
    return rad2deg(angle(vector1, vector2))


def make_homogeneous(matrix, vector):
    matrix = np.column_stack((matrix, vector))
    z = np.zeros(matrix.shape[1])
    z[z.shape[0] - 1] = 1.0
    matrix = np.row_stack((matrix, z))
    return matrix


class BasisChange:
    def __init__(self, u, v, w, origin=np.array([0.0, 0.0, 0.0])):
        self._basis = np.matrix([u, v, w]).T
        self.set_origin(origin)

    def __make_homogeneous(self, origin):
        return make_homogeneous(self._basis.I, -origin)

    def set_origin(self, origin):
        ''' create homo changing matrix from _basis and origin '''
        # reset origin
        self._changing_matrix = \
            self.__make_homogeneous(np.array([0.0, 0.0, 0.0]))
        # obtain position of origin in new basis
        new_origin = self.transform(origin)
        # set new origin
        self._changing_matrix = self.__make_homogeneous(new_origin)

    def reset_origin(self):
        self.set_origin(np.array([0.0, 0.0, 0.0]))

    def transform(self, vector):
        vector = np.array(vector)
        vector.resize(4)
        vector += np.array([0.0, 0.0, 0.0, 1.0])
        new_vector = self._changing_matrix * np.matrix(vector).T
        return np.resize(new_vector, 3)
