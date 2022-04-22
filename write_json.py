from ast import Index
from distutils.archive_util import make_archive
from unittest import removeHandler

from sklearn.datasets import make_multilabel_classification
import sneasy.parser_calc as parse
import shutil, os, glob
import json
import pandas as pd
from rdkit import Chem
import time
import requests
from openbabel import pybel


class JSONWriter:

    def __init__(self):
        pass

    def _xyz2mol(self):  
        mols = list(pybel.readfile("xyz",self.results['geometry'][0]))
        mo = mols[0]
        mo_block = mo.write("mol") 
        return Chem.MolFromMolBlock(mo_block, removeHs=False)

    def _find_id(self):
        '''Find NP-MRD ID, To do: find JEOL or HMDB too'''
        self.ID = None
        self.alt_ID = None
        id_list = pd.read_csv('resources/npmrd.csv')
        alt_id_list = pd.read_csv('resources/CH-NMR-NP_reformat.csv',sep='\t')

        SMILES = Chem.MolToSmiles(Chem.MolFromSmiles(self.results['smiles']))
        # Find NPMRD ID
        try:
            self.ID = list(id_list.loc[id_list['Can Smiles'] == SMILES]['NP_MRD_ID'])[0]
        except IndexError:
            pass
        # FIND JEOL ID
        try:
            self.alt_ID = list(alt_id_list.loc[alt_id_list['SMILES'] == SMILES]['JEOL ID'])[0]
        except IndexError:
            pass

    def _get_canonical_indices(self):

        mol = self._xyz2mol()
        nwchem_indices = self.results['shielding']['index']

        # define substructure (canonicalized version of mol)

        can_smi = Chem.MolToSmiles(Chem.MolFromSmiles(self.results['smiles']))
        can_mol = Chem.MolFromSmiles(can_smi)
        mol_H = Chem.AddHs(can_mol)
    
        lst = mol.GetSubstructMatch(mol_H)

        select_list = []
        # find the matching rdkit index 
        for idx in nwchem_indices:
            for g in lst:
                if g == idx - 1:
                    rdkit_i = lst.index(g)
            select_list.append(rdkit_i+1)

        self.results['rdkit_idx'] = select_list


    def _sort_nwchem(self):

        connectivity = self.results['connectivity']
        heavy_atom_index = []
        for pair in range(len(connectivity)):
            if connectivity[pair][0] == 'C':
                heavy_atom_index.append([connectivity[pair][2]])
            if connectivity[pair][0] == 'H':
                if connectivity[pair][1] == 'C':
                    heavy_atom_index.append([connectivity[pair][2], connectivity[pair][3]])
                else:
                    heavy_atom_index.append([connectivity[pair][2], connectivity[pair][1]])

        return heavy_atom_index

    def _dump_json(self, outfilename):

        with open(outfilename, 'w') as outfile:
            json.dump(self.compiled_json, outfile, indent=4)

        outfile.close()

    def _h_nmr_factory(self):

        INPUT = self.json_input
        data=[]
        for idx in range(len(INPUT)):
            if "h_shielding" in INPUT[idx]:
                data.append(INPUT[idx])

        return data

    def _c_nmr_factory(self):
        # Collects cnmr data together into a dictionary for each atom
        INPUT = self.json_input
        data=[]
        for idx in range(len(INPUT)):
            if "c_shielding" in INPUT[idx]:
                data.append(INPUT[idx])
        return data

    def _json_structuring(self):

        comp = {
                'np_id': self.ID,
                'alternative_id': self.alt_ID,
                'inchi_key': Chem.MolToInchiKey(Chem.MolFromSmiles(self.results['smiles'])),
                'smiles': self.results['smiles'],
                'dft_protocol': self.results['protocol'],
                'conformers': 'ffxtb, 1 conformers',
                'conversion_version': 1.0,
                'c_nmr_assignments':{
                    'solvent': self.results['protocol']['solvation'][0],
                    'spectrum': self._c_nmr_factory()
                },
                'h_nmr_assignments':{
                    'solvent': self.results['protocol']['solvation'][0],
                    'spectrum': self._h_nmr_factory()
                }
                }
        self.compiled_json = comp
    
    def _find_heavy_atom(self, idx):
        heavy_idx = None
        results = self.results
        for item in range(len(results['sorted_idx'])):
            if idx == results['sorted_idx'][item][0]:
                heavy_idx = results['sorted_idx'][item][1]
                break
        return heavy_idx     

    def _make_listdict(self):

        results = self.results

        m_c = -1.0175
        b_c = 185.65
        m_h = -1.093
        b_h = 32.031

        compiled_list = []

        for idx in range(len(results['shielding']['index'])):

            atom = results['shielding']['atom'][idx]
            if 'C' in atom:
                INPUT = {
                    'nwchem_atom_index': results['shielding']['index'][idx],
                    'rdkit_atom_index': results['rdkit_idx'][idx],
                    'c_shielding': results['shielding']['shielding'][idx],
                    'predicted_c_shift': (float(results['shielding']['shielding'][idx])-b_c)/m_c,
                    'experimental_rdkit_shift': None,
                    'experimental_munkres_shift': None
                    }
            if atom == 'H':
                INPUT = {
                    'nwchem_atom_index': results['shielding']['index'][idx],
                    'rdkit_atom_index': results['rdkit_idx'][idx],
                    'assoc_heavy_atom_nw_index': self._find_heavy_atom(results['shielding']['index'][idx]),
                    'h_shielding': results['shielding']['shielding'][idx],
                    'predicted_h_shift': (float(results['shielding']['shielding'][idx])-b_h)/m_h,
                    'experimental_rdkit_shift': None,
                    'experimental_munkres_shift': None
                    }

            compiled_list.append(INPUT)
    
        return compiled_list

    def _parse_calc(self, path: str):
        self.path = path
        nwcparse = parse.NWChemParser()
        nwcparse.load(path)
        try:
            self.results = nwcparse.parse()
        except RuntimeError:
            pass

        split_path = os.path.split(self.path)
        self.basename = split_path[-1].split('.')[0]
        smi_path = split_path[0].split('output')[0] + 'input/' + self.basename + '.smi'
        with open(smi_path, 'r') as f:
            self.results['smiles'] = f.readlines()[0]

    def _build_json(self):

        #xyz_list = glob.glob((os.path.split(self.path))[0]+'/*xyz')

        #for xyz in xyz_list:
        #    if 'cosmo' in xyz:
        #        xyz_list.remove(xyz)

        # makes mol from xyz using pybel and rdkit
        self.results['mol'] = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(self.results['smiles'])))
        # get canonical rdkit indices from mol
        self._get_canonical_indices()

        # Find NP-MRD ID
        self._find_id()

        self.results['sorted_idx'] = self._sort_nwchem()

        self.json_input = self._make_listdict()

    def _write_json(self):

        if self.ID is not None:
            NAME = self.ID
        elif self.alt_ID is not None:
            NAME = self.alt_ID
        else:
            NAME = Chem.MolToInchiKey(Chem.MolFromSmiles(self.results['smiles']))

        self.directory = 'npmrd_calc/'+NAME

        if os.path.isdir(self.directory):
            #write json into path
            pass
        else:
            os.mkdir(self.directory)

        self._json_structuring()

        self._dump_json(self.directory+'/'+NAME+'.json')   

    def _find_conformers(self):

        path = os.path.split(self.path)

        conf_path = path[0].split('shielding')[0]+'crest/'+self.basename+'/crest_conformers.xyz'

        shutil.copy2(conf_path, self.directory)

    def consolidate_data(self, path:str):
        # grabbing crest conformers and 

        self._parse_calc(path)

        self._build_json()

        self._write_json()

        self._find_conformers()


