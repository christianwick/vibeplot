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

"""
graph utilities, see http://www.python.org/doc/essays/graphs/
"""


def walk_vertices(graph, start, order, path=None):
    """depth-first search: walk order nodes from start"""
    if path is None:
        path = []
    path = path + [start]
    if len(path) == order:
        return [path]
    paths = []
    for node in graph[start]:
        if node not in path:
            newpaths = walk_vertices(graph, node, order, path)
            for newpath in newpaths:
                paths.append(newpath)
    return paths


def filter_doubles(paths):
    """path A to B is the same as B to A, remove double"""
    newpaths = []
    for idx, path in enumerate(paths):
        path.reverse()
        if path not in paths[idx + 1:]:
            newpaths.append(path)
    return newpaths


def make_paths(graph, order):
    paths = []
    for node in graph.iterkeys():
        res = walk_vertices(graph, node, order)
        for path in res:
            paths.append(path)
    return filter_doubles(paths)
