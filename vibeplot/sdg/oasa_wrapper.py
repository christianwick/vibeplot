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


import oasa


def to_oasa(vbmol):
    """convert vibeplot's molecular graph to oasa.molecule"""
    mol = oasa.molecule()
    for atom in vbmol.atoms:
        vertex = mol.create_vertex()
        # ignore coords: they need to be re-calculated!
        #vertex.coords = (atom.x, atom.y, atom.z)
        vertex.symbol = atom.symbol
        vertex.idx = atom.index  # allows maping back onto vibeplot.graph
        mol.add_vertex(v=vertex)
    edges = set()
    for atom, v in vbmol.graph.iteritems():
        for other in v:
            e = [other.index, atom.index]
            e.sort()
            e = tuple(e)
            if e not in edges:
                edges.add(e)
                mol.add_edge(e[0], e[1])
    return mol


def denovo(vbmol):
    """de novo generation with oasa"""
    mol = to_oasa(vbmol)
    oasa.coords_generator.calculate_coords(mol, bond_length=1)
    atoms = vbmol.atoms
    for atom in mol.atoms:
        x, y, z = atom.coords
        atoms[atom.idx].x_paper = x
        atoms[atom.idx].y_paper = y
        assert(atom.symbol == atoms[atom.idx].symbol)
        assert(z == 0)

__all__ = ['denovo']
