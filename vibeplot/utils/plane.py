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
from .vecalgebra import normalize


class Plane:
    def __init__(self, coord1, coord2, coord3):
        self.normal, self.distance_origin =\
                self._calculate_normal(coord1, coord2, coord3)

    def _calculate_normal(self, coord1, coord2, coord3):
        """ return Hessian normal and distance from origin """
        n = normalize(np.cross(coord2 - coord1, coord3 - coord1))
        d = -np.sum(n * coord1)
        return n, d

    def _eq_zero(self, number):
        threshold = 1e-10
        return abs(number) < threshold

    @property
    def intercepts(self):
        return -self.distance_origin / self.normal

    def intersection_normal(self, coord1):
        ''' intersection of line passing by coord1 and normal to this plane
        return coord1 if coord1 is in this plane'''
        return self.intersection(coord1, coord1 + self.normal)

    def intersection(self, coord1, coord2):
        ''' intersection of line passing by coord1, coord2 and Plane
        return coord1 if line in this Plane or None if line parallel '''
        denominator = np.dot(self.normal, coord1 - coord2)
        if self._eq_zero(denominator):  # either line in plane or parallel
            return coord1 if self.in_plane(coord1) else None
        numerator = (np.dot(self.normal, coord1) +
                     self.distance_origin)
        u = numerator / denominator
        return coord1 + u * (coord2 - coord1)

    def in_plane(self, coord_other):
        ''' return True if point at coord_other is in Plane) '''
        return self._eq_zero(np.dot(self.normal,
                                    np.array(coord_other)) +
                             self.distance_origin)

    def distance_from(self, coord_other):
        ''' return the distance between point at coord_other and Plane) '''
        return (np.dot(self.normal,
                       np.array(coord_other)) +
                self.distance_origin)

    def dihedral_angle(self, plane_other):
        ''' return dihedral angle between this and plane_other, in radians '''
        dot = np.dot(self.normal, plane_other.normal)
        try:
            dihed = math.acos(dot)
        except ValueError:
            # Rounding error
            if dot > 1.0:
                dihed = 0.0
            elif dot < 1.0:
                dihed = math.pi
            else:
                raise
        return dihed
