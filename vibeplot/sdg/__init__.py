import oasa_wrapper as oasa
import rdkit_wrapper as rdkit


import vibeplot.utils.vecalgebra as va


def sdg(molecule, sdg_algo=oasa.denovo, vdw_tolerance=1.1):
    graph = {}
    for idx, atom in enumerate(molecule.atoms):
        if atom not in graph:
            # make sure every atom is in graph
            graph[atom] = []
        for other_atom in molecule.atoms[idx+1:]:
            if (va.distance(atom.coords, other_atom.coords) <
                    (atom.radius + other_atom.radius) * vdw_tolerance):
                graph.setdefault(atom, []).append(other_atom)
                graph.setdefault(other_atom, []).append(atom)
    molecule.graph = graph
    sdg_algo(molecule)
