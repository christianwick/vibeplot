class def_dict(dict):
    def __missing__(self, key):
        return 0


class Molfile(object):
    """
    create string with the molecule in the molfile format

    molfile format: DOI: 10.1021/ci00007a012
    """
    
    counts_line_fmt = ''.join(
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
        ])
    atom_block_fmt = ''.join(
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
        ])
    bond_block_fmt = ''.join(
        ["%(first atom)3i",
         "%(second atom)3i",
         "%(bond type)3i",
         "%(zero)3i",  # bond stereo
         "%(zero)3i",  # not used
         "%(zero)3i",  # bond topology
         "%(zero)3i",  # reacting center status
        ])

    def __init__(self, mol):
        self.mol = mol

    def generate_bond_list(self):
        bond_list = []
        for parent in self.mol.graph:
            for child in self.mol.graph[parent]:
                if [child, parent] not in bond_list:
                    bond_list.append([parent, child])
        return bond_list

    def header(self):
        return "\n".join(["", "vibeplot", ""])

    def counts_line(self):
        return Molfile.counts_line_fmt % def_dict(
            {"atom_count": len(self.mol.atoms),
             "bond_count": self.mol.count_bonds(),
             "ctab_vers": "V2000"})

    def atom_block(self):
        return "\n".join([
            Molfile.atom_block_fmt % def_dict(
                {"x": atom.x_paper,
                 "y": atom.y_paper,
                 "z": 0,
                 "symbol": atom.symbol,})
            for atom in self.mol.atoms])

    def bond_block(self):
        bond_block = [Molfile.bond_block_fmt %
                      def_dict({"first atom": atom1.index + 1,
                                "second atom": atom2.index + 1,
                                "bond type": 1,})
                      for atom1, atom2 in self.generate_bond_list()]
        bond_block.sort()
        return "\n".join(bond_block)

    def properties_block(self):
        return "M  END"

    def eof(self):
        return "$$$$"

    @classmethod
    def mol_to_molfile(self, mol):
        molfile = Molfile(mol)
        return "\n".join([molfile.header(), molfile.counts_line(),
                          molfile.atom_block(),
                          molfile.bond_block(),
                          molfile.properties_block(),
                          molfile.eof()])
