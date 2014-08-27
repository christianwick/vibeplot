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

from rdkit import Chem
from rdkit.Chem import AllChem
from vibeplot.datafile.molfile import Molfile


def sdg(vbmol, algo='conversion'):
    mol = Chem.MolFromMolBlock(Molfile.mol_to_molfile(vbmol), removeHs=False)
    if algo == 'conversion':
        AllChem.GenerateDepictionMatching3DStructure(mol, mol)
    elif algo == 'de novo':
        AllChem.Compute2DCoords(mol)
    else:
        raise NameError
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        x, y, z = mol.GetConformer().GetAtomPosition(idx)
        vbmol.atoms[idx].x_paper = x
        vbmol.atoms[idx].y_paper = y
        assert(atom.GetSymbol() == vbmol.atoms[idx].symbol)
        assert(atom.GetIdx() == vbmol.atoms[idx].index)
        assert(z == 0)

def conversion(vbmol):
    """3D -> 2D conversion with RDKit"""
    sdg(vbmol, algo='conversion')

def denovo(vbmol):
    """de novo generation with RDKit"""
    sdg(vbmol, algo='de novo')

__all__ = ['denovo', 'conversion']
