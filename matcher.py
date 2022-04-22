import json
import glob,os
import pandas as pd
from sneasy import parser_exp
from rdkit import Chem
from munkres import Munkres
import numpy as np

class Matcher:
    def __init__(self):
        self.id = None
        self.experimental_contents = None
        self.predicted_contents = None
        self.experimental_path = None
        self.predicted_path = None
        self.path = None
        self.smiles = None
        self.database_reference = None


    def _load_separate_exp(self, path: str):

        extension = os.path.splitext(path)[-1].lower()

        if 'nmrml' in extension:
            parser = parser_exp.nmrMLParser()
            parser.load(path)
            experimental_contents = parser.peak_dict
            experimental_contents['protocol'] = parser.exp_protocol

        if 'json' in extension:
            parser = parser_exp.JSONParser()
            parser.load(path, self.index)
            experimental_contents = parser.peak_dict
            experimental_contents['protocol'] = parser.exp_protocol


        if 'jdx' in extension:
            parser = parser_exp.JDXParser()
            parser.load(path)
            experimental_contents = parser.peak_dict
            experimental_contents['protocol'] = parser.exp_protocol

        return experimental_contents

    def _update_calc_load(self):

        self.atom_types = ['h', 'c']

        for atom_type in self.atom_types:
            predict = self.predicted_contents[atom_type+'_nmr_assignments']['spectrum']
            for atom in predict:
                atom['experimental_rdkit_shift'] = None
                atom['experimental_munkres_shift'] = None

            self.predicted_contents[atom_type+'_nmr_assignments']['spectrum'] = predict

    def load_exp(self):

        exp_dict = {}

        for atom_type in self.atom_types:
            path = self.experimental_path[atom_type]
            exp_dict[atom_type] = self._load_separate_exp(path)

        self.experimental_contents = exp_dict

    def load_calc(self, path:str):

        self.predicted_path = path
        #self.id = os.path.split(path)[-1].split('.')[0]
        with open(path, 'r+') as file:
            js = json.load(file)
        self.predicted_contents = js

        if js['np_id'] is not None:
          self.id = js['np_id']
        if self.id is None:
          self.id = js['alternative_id']

        #self._update_calc_load()

    def _find_jeol_exp(self):

        # Load in csv of exp location
        df2 = pd.read_csv('resources/jeol_entries.csv')

        if self.id.startswith('NP'):

            h_list = glob.glob('exp/jeol/nmrMLs_for_protons/1H_for_JSviewer_NPMRDID_'+self.id+'*')
            c_list = glob.glob('exp/jeol/nmrMLs_for_carbons/13C_for_JSviewer_NPMRDID_'+self.id+'*')

            experimental_path = [h_list[0], c_list[0]]

        if self.id.startswith('JEOL'):
            #use JEOL to find file
            entry = df2.loc[df2['JEOL ID'] == int(self.id.split('L')[-1])]

            experimental_path = [list(entry['H FILENAME'])[0],list(entry['C FILENAME'])[0]]

        h_path = experimental_path[0].strip()
        c_path = experimental_path[1].strip()

        self.experimental_path = {'h': h_path, 'c': c_path}

        self.atom_types = ['h','c']

    def _find_json_exp(self):

        df = pd.read_csv('resources/json_curator_entries.csv')

        if '-' in self.predicted_path.split('/')[-1]:
            inchikey = self.predicted_path.split('/')[-1].split('.')[0]
        else:
            self.smiles = self.predicted_contents['smiles']
            inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(self.smiles))

        if inchikey is not None:
            entry = df.loc[df['INCHIKey'] == inchikey]
            path = 'exp/'+list(entry['FILENAME'])[0]

        self.experimental_path = {'h': path, 'c': path}
        self.index = int(list(entry['Index'])[0])
        self.atom_types = ['h', 'c']

    def _find_jcamp_exp(self):
        df1 = pd.read_csv('resources/jcamp_entries.csv')
        df2 = pd.read_csv('resources/jeol_np_key.csv')

        # use JEOL to find NP
        if self.id.startswith('JEOL'):
            nme = df2.loc[df2['JEOL ID'] == int(jeol.split('L')[-1])]
            np_id = list(nme['NP_MRD_ID'])[0]

            entry = df1.loc[df1['NPMRD ID'] == np_id]
            paths = list(entry['FILENAME'])

        if self.id.startswith('NP'):
            entry = df1.loc[df1['NPMRD ID'] == np_id]
            paths = list(entry['FILENAME'])

        h_path = None
        c_path = None
        for path in paths:

            if '1HNMR' in path:
                h_path = path
            if '13CNMR' in path:
                c_path = path

        self.experimental_path = {'h': h_path, 'c': c_path}
        self.atom_types = ['h', 'c']

    def _find_npmrd_from_smiles(self):

        df1 = pd.read_csv('resources/npmrd.csv')

        nme = df1.loc[df1['Can Smiles'] == self.predicted_contents['smiles']]
        if nme is not None:
            self.id = nme
        

    def find_exp(self):

        try:
            self._find_jeol_exp()
            print(self.experimental_path)
        except:
            pass
        try:
            self._find_jcamp_exp()
            print(self.experimental_path)
        except:
            pass
        try:
            self._find_json_exp()
            print(self.experimental_path)
        except:
            pass
        try:
            self.load_exp()
        except:
            pass


    def save(self, fmt='json'):
        path = path.strip()
        extension = os.path.splitext(path)[-1].lower()

        if fmt == 'nmrML':
            return save_nmrml(path)

        if fmt == 'json':
            return save_json(path)

        if fmt == 'jdx':
            return save_jdx(path)
        
        raise IOError('Extension {} not recognized.'.format(extension))

    def _connect_predicted_shifts(self):

        js = self.predicted_contents

        all_shift = []

        for i in js['c_nmr_assignments']['spectrum']:
            all_shift.append([i['nwchem_atom_index']])

        h_id = []

        for i in js['h_nmr_assignments']['spectrum']:
            h_id.append([i['predicted_h_shift'], i['assoc_heavy_atom_nw_index'], i['nwchem_atom_index']])

        all_indices = []
        for i in all_shift:
            # Fill out H's with shifts
            idx_list = [i[0]]
            for j in h_id:
                if i[0] == j[1]:
                    idx_list.append(j[2])
                    i.append(j[0])
            for shift in js['c_nmr_assignments']['spectrum']:
                if i[0] == shift['nwchem_atom_index']:
                    i[0] = shift['predicted_c_shift']

            all_indices.append(idx_list)
        for i in all_shift:
            if len(i) == 1:
                i.append(10000.0)

        return {'Shifts': all_shift, 'Indices': all_indices}

    def _connect_experimental_shifts(self):

        exp = self.experimental_contents

        smi = self.smiles
        mol = Chem.MolFromSmiles(smi)
        mol = Chem.AddHs(mol)

        matrix = Chem.GetAdjacencyMatrix(mol)
        numheavyatoms = mol.GetNumHeavyAtoms()

        collapsed_matrix = []
        for i in range(len(matrix)):
            li = [g for g, x in enumerate(matrix[i]) if x == 1]
            nu_li = []
            for nu in li:
                if nu >= numheavyatoms:
                    nu_li.append(nu)
            collapsed_matrix.append(nu_li)

        all_shift = []

        for i in range(len(exp['c']['RDKit Index'])):
            all_shift.append([exp['c']['RDKit Index'][i]])      

        h_id = []

        for i in range(len(exp['h']['RDKit Index'])):
            if isinstance(exp['h']['RDKit Index'][i], int):
                h_id.append([exp['h']['Shift'][i], exp['h']['RDKit Index'][i]])
            else:
                for j in exp['h']['RDKit Index'][i]:
                    h_id.append([exp['h']['Shift'][i], j])


        all_indices = []
        for i in all_shift:
            idx_list = [i[0]]
            for idx in collapsed_matrix[i[0]-1]:
                for group in h_id:
                    if group[1]-1 == idx:
                        i.append(group[0])
                        idx_list.append(group[1])
            if len(idx_list) == 1:
                    i.append(10000.0)
            all_indices.append(idx_list)

            for g in range(len(exp['c']['RDKit Index'])):
                if exp['c']['RDKit Index'][g] == i[0]:
                    i[0] = exp['c']['Shift'][g]

        return {'Shifts': all_shift, 'Indices': all_indices}

    def _build_cost_matrix(self, atom_type):


        exp_shift = self.experimental_contents[atom_type]['Shift']

        calc_shift = []

        for i in self.predicted_contents[atom_type+'_nmr_assignments']['spectrum']:
            calc_shift.append(i['predicted_'+atom_type+'_shift'])

        matrix = []
        
        for exp in range(len(exp_shift)):
            line = []
            for calc in range(len(calc_shift)):

                dif = np.abs(calc_shift[calc] - exp_shift[exp])
                line.append(dif)

            matrix.append(line)
        
        self.cost_matrix = matrix
        self.exp_shift = exp_shift
        self.calc_shift = calc_shift

    def _connect_shifts(self):
        try:
            result = self._connect_predicted_shifts()
            self.connected_pred = result
        except:
            pass

        try:
            result = self._connect_experimental_shifts()
            self.connected_exp = result
        except:
            pass

    def _match_random(self):

        # read in rdkit matches

        def swap_random(seq):
            idx = range(len(seq))
            i1, i2 = random.sample(idx, 2)
            seq[i1], seq[i2] = seq[i2], seq[i1]

        for atom_type in self.atom_types:
            predict = self.predicted_contents[atom_type+'_nmr_assignments']['spectrum']
            indexes = []
            # 1. Randomly choose row in matrix, randomly choose element in row.

            while len(indexes) < n_exp:

                row = random.choice(range(len(matrix)))
                element = random.choice(range(len(matrix[row])))
                monte = random.random()
                if len(indexes) == n_calc-1:
                    element = np.where(matrix[row] != 0)[0][0]
                    monte = 1
                if matrix[row][element] != 0:

                    if matrix[row][element]/(np.sum(matrix[row])) >= monte:
                        indexes.append([row, element])
                        for i in matrix:
                            i[element] = 0
                        for i in matrix[row]:
                            i = 0
            indexes = list(sorted(indexes))

            reorder_calc_shield = []
            k = 0
            while k <= len(indexes)-1:
                reorder_calc_shield.append(calc_shield[indexes[k][1]])
                k += 1

            return reorder_calc_shield

    def _match_rdkit(self):

        for atom_type in self.atom_types:
            predict = self.predicted_contents[atom_type+'_nmr_assignments']['spectrum']
            exp_peak_dict = self.experimental_contents[atom_type]

            for atom in predict:
                exp_shift = None
                rdkit_index = int(atom['rdkit_atom_index'])
                try:
                    for i in range(len(exp_peak_dict['RDKit Index'])):
                        if isinstance(exp_peak_dict['RDKit Index'][i],int):
                            if rdkit_index == exp_peak_dict['RDKit Index'][i]:
                                exp_shift = exp_peak_dict['Shift'][i]
                                continue
                        else:
                            if rdkit_index in exp_peak_dict['RDKit Index'][i]:
                                exp_shift = exp_peak_dict['Shift'][i]
                                continue
                except:
                    continue 
                atom['experimental_rdkit_shift'] = exp_shift

            self.predicted_contents[atom_type+'_nmr_assignments']['spectrum'] = predict
            self.predicted_contents[atom_type+'_nmr_assignments']['exp protocol'] = self.experimental_contents[atom_type]['protocol']

    def _match_munkres(self):

        for atom_type in self.atom_types:
            # Build cost matrix
            self._build_cost_matrix(atom_type)
            munk = Munkres()
            # Obtain indices after munkres matching
            indexes = munk.compute(self.cost_matrix)

            # Fill arrays with shifts in order according to indexes
            exp_shifts = []
            calc_shifts = []
            for index in indexes:
                exp_shifts.append(self.exp_shift[index[0]])
                calc_shifts.append(self.calc_shift[index[1]])

            # Zipping the exp shifts and calc shifts?
            pair_list = []
            for shift in range(len(exp_shifts)):
                pair = [exp_shifts[shift], calc_shifts[shift]]
                pair_list.append(pair)
            
            predict = self.predicted_contents[atom_type+'_nmr_assignments']['spectrum']
            for atom in predict:
                exp_shift = None
                for i in range(len(pair_list)):
                    if pair_list[i][1] == float(atom['predicted_'+atom_type+'_shift']):
                        exp_shift = pair_list[i][0]
                        break
                atom['experimental_munkres_shift'] = exp_shift
        
        self.predicted_contents[atom_type+'_nmr_assignments']['spectrum'] = predict
        self.predicted_contents[atom_type+'_nmr_assignments']['exp protocol'] = self.experimental_contents[atom_type]['protocol']

    def match(self):
        try:
            self._match_rdkit()
        except:
            pass
        try:
            self._match_munkres()
        except:
            pass

    def save_json(self):
        '''
        Write generated NWChem configuration to file.

        '''
        with open(self.predicted_path, 'w') as outfile:
            json.dump(self.predicted_contents, outfile, indent=4)

        outfile.close()

    def score(self):

        self.atom_types = ['h', 'c']

        for atom_type in self.atom_types:
            predict = self.predicted_contents[atom_type+'_nmr_assignments']['spectrum']
            rd_error = 0
            munkres_error = 0
            rd_count = 0
            munkres_count = 0
            for atom in predict:
                if atom['experimental_rdkit_shift'] is not None:
                    exp_shift = atom['experimental_rdkit_shift']
                    calc_shift = atom['predicted_'+atom_type+'_shift']
                    rd_error += abs(exp_shift - calc_shift)
                    rd_count += 1
                if atom['experimental_munkres_shift'] is not None:
                    exp_shift = atom['experimental_rdkit_shift']
                    calc_shift = atom['predicted_'+atom_type+'_shift']
                    munkres_error += abs(exp_shift - calc_shift)
                    munkres_count += 1
            if rd_error != 0:
                rd_mae = rd_error/rd_count
            else:
                rd_mae = None
            if munkres_error != 0:
                munkres_mae = munkres_error/munkres_count
            else:
                munkres_mae = None
            self.predicted_contents[atom_type+'_nmr_assignments']['munkres_mean_average_error'] = munkres_mae
            self.predicted_contents[atom_type+'_nmr_assignments']['rdkit_mean_average_error'] = rd_mae

    def set_p(self, ref_path):
        self._set_db_reference(ref_path)
        self._set_smiles()

def update_experimental_values(predict_path, ref_path):

    match = Matcher()

    # load predicted shifts
    match.load(predict_path)
    match.set_p(ref_path = ref_path)

    # find experimental data
    match.find_experimental_data()
    match.set_atom_type()

    #use rdkit indices to match experimental shifts
    match.match_rdkit_indices()

    #save file
    match.save_json()


