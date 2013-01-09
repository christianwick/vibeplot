import numpy as np
import vibeplot.chemistry as chem


class def_dict(dict):
    def __missing__(self, key):
        return 0


def header():
    return "\n".join(["", "vibeplot", ""])

def counts_line(mol):
    return ''.join(
        ["%(atom_count)3i",
         "%(bond_count)3i",
         "%(atom_lists)3i",
         "%(zero)3i",     # obsolete
         "%(chiral)3i",   # chiral
         "%(stext)3i",    # stext
         "%(reac)3i",     # reaction components + 1
         "%(reac_ct)3i",  # reactants count
         "%(prod_ct)3i",  # product count
         "%(inter_ct)3i", # intermediates count
         "%(prop)3i",     # nb lines of additional properties
         "%(ctab_vers)6s",
        ]) % def_dict(atom_count=len(mol.atoms),
                      bond_count=mol.count_bonds(),
                      ctab_vers="V2000")

def atom_block(mol):
    return "\n".join([
        ''.join(
            ["%(x)10.4f",
             "%(y)10.4f",
             "%(z)10.4f",
             " ",
             "%(symbol)-3s",    # symbol
             "%(zero)2i",       # mass difference
             "%(zero)3i",       # charge
             "%(zero)3i",       # stereo parity
             "%(hydrogens)3i",  # hydrogen count + 1
             "%(zero)3i",       # stereo care
             "%(valence)3i",    # valence
             "%(zero)3i",       # H0 designator
             "%(zero)3i",       # reaction cmpt type
             "%(zero)3i",       # reaction cmpt nb
             "%(mapping)3i",    # atom-atom mapping number
             "%(zero)3i",       # inversion/retension flag
             "%(zero)3i",       # exact change
            ]) % def_dict(x=atom.x_paper,
                          y=atom.y_paper,
                          z=0,
                          symbol=atom.symbol)
        for atom in mol.atoms])

def bond_block(mol):
    def generate_bond_list():
        bond_list = []
        for parent in mol.graph:
            for child in mol.graph[parent]:
                if [child, parent] not in bond_list:
                    bond_list.append([parent, child])
        return bond_list

    bond_block = ["".join(
        ["%(first_atom)3i",
         "%(second_atom)3i",
         "%(bond_type)3i",
         "%(zero)3i",  # bond stereo
         "%(zero)3i",  # not used
         "%(zero)3i",  # bond topology
         "%(zero)3i",  # reacting center status
        ]) % def_dict(first_atom=atom1.index + 1,
                      second_atom=atom2.index + 1,
                      bond_type=1)
        for atom1, atom2 in generate_bond_list()]
    bond_block.sort()
    return "\n".join(bond_block)

def properties_block():
    return "M  END"

def eof():
    return "$$$$"


def write_molfile(mol):
    """Create string with the molecule in the molfile format.

    Parameters
    ----------
    mol: 'vibeplot.chemistry.Molecule'

    Returns
    -------
    molfile: string

    Notes
    -----
    molfile format: DOI: 10.1021/ci00007a012

    """
    return "\n".join([header(),
                      counts_line(mol),
                      atom_block(mol),
                      bond_block(mol),
                      properties_block(),
                      eof()])


def parse_molfile(filename):
    """Parse file with the molfile format.

    Parameters
    ----------
    filename: str

    Returns
    -------
    molecule: 'vibeplot.chemistry.Molecule'
    spectrum: dict
        None

    """
    atom_list = []
    graph = {}
    with open(filename, "r") as f:
        for skip in range(3):
            f.readline()
        header = f.readline().split()
        for enum_atom in range(int(header[0])):
            atom_line = f.readline().split()
            atom = chem.Atom(symbol=atom_line[3],
                             index=enum_atom,
                             xyz=np.array(atom_line[0:3], dtype=float))
            atom_list.append(atom)
        for enum_bond in range(int(header[1])):
            bond_line = f.readline().split()
            child_idx = int(bond_line[0]) - 1
            parent_idx = int(bond_line[1]) - 1
            parent = atom_list[parent_idx]
            child = atom_list[child_idx]
            graph.setdefault(parent, []).append(child)
            graph.setdefault(child, []).append(parent)
        return chem.Molecule(graph), None


if __name__ == "__main__":
    print (parse_molfile("tmp/test.sdl"))

