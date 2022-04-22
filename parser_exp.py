import json
import nmrglue
from rdkit import Chem
from bs4 import BeautifulSoup as bs

class nmrMLParser:
    def __init__(self):
        self.contents=None
        self.path=None

    def _load_nmrml(self, path: str):
        self.path = path
        with open(self.path, 'r') as file:
            content = file.readlines()
            content = "".join(content)
            self.contents = bs(content, "lxml")

    def load(self, path: str):

        self._load_nmrml(path)
        self._set_exp_protocol()
        self._set_atom_type()
        self._reformat_peaks_table()

    def _set_atom_type(self):
        spectrum = []
        if 'h_nmr' in self.contents:
            spectrum.append('h_nmr')
        if 'c_nmr' in self.contents:
            spectrum.append('c_nmr')

        self.atom_types = spectrum

    def _set_exp_protocol(self):
        st = None
        freq = None
        solv = None

        try:
            for standard in self.contents.findAll('chemicalshiftstandard'):
                st = standard.attrs['name']
        except:
            pass
        try:
            for frequency in self.contents.findAll('effectiveexcitationfield'):
                freq = frequency.attrs['value']
        except:
            pass

        try:
            for solvent in self.contents.findAll('fieldfrequencylock'):
                solv = solvent.attrs['fieldfrequencylock']
        except:
            pass

        self.exp_protocol = {'Standard': st, 'Frequency': freq, 'Solvent': solv}    

    def _reformat_peaks_table(self):
        assign = []
        idx = []

        for at in self.contents.findAll('atoms'):
            at = at.attrs['atomrefs'].split(' ')
            if len(at) == 1:
                idx.append(int(at[0]))
            else:
                idx.append(list(map(int, at)))

        for peak in self.contents.findAll('multiplet'):
            assign.append(float(peak.attrs['center']))

        self.peak_dict = {'RDKit Index': idx,
                          'Shift': assign}


        
class JSONParser:
    def __init__(self):
        self.contents = None
        self.path = None

    def _load_json(self, path: str, index):
        self.path = path
        with open(path) as f:
            js = json.load(f)
        self.contents = js[index]

    def load(self, path: str, index):

        self._load_json(path, index)
        self._set_atom_type()
        self._set_exp_protocol()
        self._reformat_peaks_table()

    def _set_atom_type(self):
        spectrum = []
        if 'h_nmr' in self.contents:
            spectrum.append('h_nmr')
        if 'c_nmr' in self.contents:
            spectrum.append('c_nmr')

        self.atom_types = spectrum

    def _set_exp_protocol(self):

        st = None
        freq = None
        solv = None

        for atom_type in self.atom_types:
            try:
                st = self.contents[atom_type]['reference']
            except:
                pass
            try:
                freq = self.contents[atom_type]['frequency']
            except:
                pass
            try:
                solv = self.contents[atom_type]['solvent']
            except:
                pass
 
        self.exp_protocol = {'Standard': st, 'Frequency': freq, 'Solvent': solv}    


    def _reformat_peaks_table(self):
        idx = []
        assign = []

        for atom_type in self.atom_types:

            peaks = self.contents[atom_type]['spectrum']

            for peak in peaks:
                idx.append(peak['rdkit_index'])
                assign.append(peak['shift'])


        self.peak_dict = {'RDKit Index': idx,
                          'Shift': assign}

    def _make_combined_peak_matrix(self):

        smi = self.contents['smiles']
        mol = Chem.MolFromSmiles(smi)
        mol = Chem.AddHs(mol)

        matrix = Chem.GetAdjacencyMatrix(mol)

        numheavyatoms = mol.GetNumHeavyAtoms()

        combined_matrix = []
        

        for atom in range(numheavyatoms):
            group = [self.peak_dict['Shift'][atom], 'H', 'H', 'H']
            for i in range(numheavyatoms, len(matrix), 1):
                g_spot = 1
                if matrix[atom][i] == 1:
                    group[g_spot] = self.peak_dict['Shift'][i]
                    g_spot += 1
            combined_matrix.append(group)



class JDXParser:

    def __init__(self):
        self.contents = None
        self.path = None

    def _load_jdx(self, path: str):
        self.contents = nmrglue.jcampdx.read(path)[0]
        self.mol = Chem.MolFromMolBlock(self.contents['_datatype_LINK'][0]['$MOLFILE'][0])

    def load(self, path: str):

        self._load_jdx(path)
        self._set_atom_type()
        self._reformat_peaks_table()
        self._set_exp_protocol()

    def _set_exp_protocol(self):

        st = None
        freq = None
        solv = None

        for atom_type in self.atom_types:

            try:
                st = self.contents['.SHIFTREFERENCE'][0].split(',')[1]
            except:
                pass
            try:
                freq = self.contents['.OBSERVEFREQUENCY'][0]
            except:
                pass
            try:
                solv = self.contents[atom_type]['solvent']
            except:
                pass

        self.exp_protocol = {'Standard': st, 'Frequency': freq, 'Solvent': solv}    


    def _set_atom_type(self):
        self.atom_types = self.contents['_datatype_NMRPEAKTABLE'][0]['.OBSERVENUCLEUS']
    
    def _reformat_peaks_table(self):
        peaks = self.contents['_datatype_NMRPEAKASSIGNMENTS'][0]['PEAKASSIGNMENTS'][0].split('\n')

        idx = []
        assign = []
        for peak in range(1,len(peaks)):

            new_format = peaks[peak].strip('()').split(',')
            idx.append(int(new_format[2].strip(' <').strip('>')))
            assign.append(float(new_format[0]))

        self.peak_dict = {'RDKit Index': idx,
                          'Shift': assign}





