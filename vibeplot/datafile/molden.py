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
import vibeplot.chemistry as chem
from base import AbstractDataFile


def pm_to_au(x):
    return float(x) * 0.01890

def au_to_pm(x):
    return float(x) * 52.918


class MoldenFile(AbstractDataFile):
    ''' Parse file in molden format, for info, see
    http://www.cmbi.ru.nl/molden/molden_format.html '''

    headers = {'[Atoms]': 'atoms',
               '[FR-NORM-COORD]': 'displacements',
               '[FREQ]': 'frequency',
               '[GEOMETRIES]': 'geometries'}

    def __init__(self, filename='molden.input'):
        self.freq_list = []
        super(MoldenFile, self).__init__(filename)

    def _parse_atom_section(self, symbol, index, mass, x, y, z):
        symbol = symbol.capitalize()
        atom = chem.Atom(
            **{'symbol': symbol,
               'index': int(index) - 1,
               'mass': int(mass),
              })
        atom.x = au_to_pm(x)
        atom.y = au_to_pm(y)
        atom.z = au_to_pm(z)
        self.atoms.append(atom)

    def _parse_frequency_section(self, wavenumber):
        self.freq_list.append(float(wavenumber))

    def _parse_displacements_header_section(self, v, vibration_idx):
        self.vibration_idx = int(vibration_idx)
        self.atom_idx = 0
        vibration = chem.Vibration()
        vibration.frequency = self.freq_list[self.vibration_idx - 1]
        self.vibrations.append(vibration)

    def _parse_displacements_section(self, x, y, z):
        xyz = np.array(
            (au_to_pm(x),
             au_to_pm(y),
             au_to_pm(z)))
        xyz += self.atoms[self.atom_idx].coords
        self.atoms[self.atom_idx].normal_coords[
            self.freq_list[self.vibration_idx - 1]] = xyz
        self.atom_idx += 1

    def parse(self):
        with open(self.filename, 'r') as f:
            for line in f:
                # Skip empty lines
                if not line.strip():
                    continue
                firstword = line.split()[0]
                # Detect section
                if line.startswith('['):
                    header = firstword
                    section = MoldenFile.headers.get(header)
                    continue
                # Parse line
                if section is 'atoms':
                    callback = self._parse_atom_section
                elif section is 'frequency':
                    callback = self._parse_frequency_section
                elif section is 'displacements':
                    if firstword == 'vibration':
                        callback = self._parse_displacements_header_section
                    else:
                        callback = self._parse_displacements_section
                else:
                    callback = None
                if callback:
                    callback(*line.split())
        vibrations = {}
        for vib in self.vibrations:
            if vib.frequency > chem.Molecule.THRESHOLD:
                vibrations[vib.key] = vib
        return chem.Molecule(self.atoms, vibrations)
