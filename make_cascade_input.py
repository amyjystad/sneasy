import matcher


# Make pickle file with data for cascade input. contains mol_id, Mol(rdkit mol object), n_atoms, atom_index, Shift 

m = matcher.Matcher()

mol_id = []
Mol = []
n_atoms = []
atom_index = []
Shift = []

#iterate 


dict = {'mol_id': mol_id,
        'Mol': Mol,
        'n_atoms': n_atoms,
        'Shift': Shift}

df = pd.DataFrame(dict)

df.to_pickle('jeol_input.p')

# load model

import sys
sys.path.append('../../')

from keras.models import Model, load_model

from nfp.layers import (MessageLayer, GRUStep, Squeeze, EdgeNetwork,
                               ReduceAtomToMol, ReduceBondToAtom,
                               GatherAtomToBond, ReduceAtomToPro)
from nfp.models import GraphModel


model = load_model(filepath, custom_objects={'GraphModel': GraphModel,
                                             'Squeeze': Squeeze,
                                             'GatherAtomToBond': GatherAtomToBond,
                                             'ReduceBondToAtom': ReduceBondToAtom,
                                             'ReduceAtomToPro': ReduceAtomToPro})

# process data from pickle file

def rbf_expansion(distances, mu=0, delta=0.1, kmax=256):
    k = np.arange(0, kmax)
    logits = -(np.atleast_2d(distances).T - (-mu + delta * k))**2 / delta
    return np.exp(logits)

def atomic_number_tokenizer(atom):
    return atom.GetNumRadicalElectrons()

def _compute_stacked_offsets(sizes, repeats):
    return np.repeat(np.cumsum(np.hstack([0, sizes[:-1]])), repeats)

def process_data(batch_data):
        batch_data['distance_rbf'] = rbf_expansion(batch_data['distance'])

        offset = _compute_stacked_offsets(
            batch_data['n_pro'], batch_data['n_atom'])

        offset = np.where(batch_data['atom_index']>=0, offset, 0)
        batch_data['atom_index'] += offset

        del batch_data['n_atom']
        del batch_data['n_bond']
        del batch_data['distance']

        return batch_data

# we have our model.input and our input for each molecule
# need to mimic model.input but fill it in by matching model.input[x].name
# with item in input directory

[ -26.093033],
       [  42.375546],
       [  63.423428],
       [ -14.78653 ],
       [ -39.473724],
       [-115.10062 ],
       [ -38.87964 ],
       [ -63.362305],
       [  52.70403 ],
       [  73.39171 ],
       [ -75.365265],
       [-117.84874 ]



