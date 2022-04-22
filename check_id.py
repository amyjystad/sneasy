import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import json
from sklearn.linear_model import LinearRegression
import pandas as pd
from rdkit import Chem
import math

files = glob.glob('calc/json_curator/*/*/*json')

df = pd.read_csv('id_smiles_key/npmrd.csv')
df2 = pd.read_csv('jeol_np_key.csv')

def clear_shift_entries(js):
  atom_type = ['c', 'h']
  for at in atom_type:
    js[at+'_nmr']['munkres_mean_average_error'] = None
    js[at+'_nmr']['rdkit_mean_average_error'] = None
    for atom in js[at+'_nmr']['spectrum']:
      atom['exp_shift_munkres_match'] = None
      atom['exp_shift_rdkit_match'] = None
      if 'experimental_shift' in atom:
        del atom['experimental_shift']
  if 'alt_id' in js:
    del js['alt_id']


def save_json(path, js):
   with open(path, 'w') as outfile:
      json.dump(js, outfile, indent=4)

   outfile.close()

for fi in files:
  print(fi)
  with open(fi) as f:
    js = json.load(f)

  smi = js['smiles']

  ID = js['np_id']

  ALT_ID = js['alternative_id']

  smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi))

  row = df.loc[df['Can Smiles'] == smi]

  print(ID, ALT_ID)
  js['np_id'] = None
  js['alternative_id'] = None

  if row.empty:
    row2 = df2.loc[df2['Canonical SMILES'] == smi]
    if row2.empty:
      pass
    else:
      print(row2)
      js['np_id'] = list(row2['NPMRD ID'])[0]
      js['alternative_id'] = list(row2['JEOL ID'])[0]
  else:
    print(row)
    js['np_id'] = list(row['NP_MRD_ID'])[0]
    if math.isnan(list(row['JEOL_ID'])[0]):
      js['alternative_id'] = None
    else:
      js['alternative_id'] = 'JEOL'+ str(int(list(row['JEOL_ID'])[0]))

  clear_shift_entries(js)
  save_json(fi, js)


