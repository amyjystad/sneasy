import argparse
from scipy.spatial import distance
import pandas as pd
import numpy as np
from munkres import Munkres, print_matrix
import os
import matplotlib.pyplot as plt  # To visualize
from sklearn.linear_model import LinearRegression
import glob


def convert_to_shifts(atom_list, array):
    C = 187.614272
    H = 31.694003

    shift_c = []
    shift_h = []
    shield_c = []
    shield_h = []
    g = 0
    while g <= len(array)-1:
        if atom_list[g] == 'C':
            shift_c.append(C-array[g])
            shield_c.append(array[g])
        if atom_list[g] == 'H':
            shift_h.append(H-array[g])
            shield_h.append(array[g])
        g += 1

    return shift_c, shift_h, shield_c, shield_h

def load_file(filename):

    df = pd.read_csv(filename, sep='\t')

    df_H = df.loc[df['Atom'] == 'H']
    df_C = df.loc[df['Atom'] == 'C']

    return df_H, df_C

def separate_exp(atom_list, array, atom='C'):
    exp = []

    g = 0
    while g <= len(array)-1:
        if atom_list[g] == atom and array[g] >= 0:
            exp.append(array[g])
        g += 1

    return exp


def read_shifts(filename):

    df = pd.read_csv(filename, sep='\t')
    atom = list(df.Atom)
    exp = df['Exp Shift'].replace('Nan', 0)
    method = df.columns[2:]

    return df, atom, exp, method


def idx(filename, exp, calc_shift, calc_shield, name, atm='C'):

    matrix = np.zeros((len(exp), len(calc_shift)))

    i = 0
    while i <= len(exp)-1:
        j = 0
        while j <= len(calc_shift)-1:
            matrix[i][j] += np.abs(exp[i]-calc_shift[j])
            j += 1
        i += 1

    m = Munkres()
    indexes = m.compute(matrix)

    print(indexes)
    new_calc = []
    k = 0
    while k <= len(indexes)-1:
        print(k, calc_shield[indexes[k][1]])
        new_calc.append(calc_shield[indexes[k][1]])
        k += 1
    df = pd.DataFrame(list(zip(exp, new_calc)),
                      columns=['Exp Shift', 'Calc Shield'])

    method = str(name).split(' ')[0]
    opt = str(name).split(' ')[1]
    basename = 'new_' + \
        os.path.splitext(os.path.basename(filename))[0] + '_' + method+'_'+opt+'_'+atm+'.csv'

    df.to_csv(basename, sep='\t')


def linear(dir, atom):
    name = str(dir)+'/*'+str(atom) + '.csv'
    all_files = glob.glob(name)

    li = []

    for filename in all_files:
        df = pd.read_csv(filename, sep='\t', index_col=None, header=0)
        df = df[df['Exp Shift'] != 0]
        li.append(df)

    frame = pd.concat(li, axis=0, ignore_index=True)

    #expeirmental
    X = frame.iloc[:, 1].values.reshape(-1, 1)  # values converts it into a numpy array
    #calculated
    Y = frame.iloc[:, 2].values.reshape(-1, 1)  # -1 means that calculate the dimension of rows, but have 1 column
    linear_regressor = LinearRegression()  # create object for the class
    reg = linear_regressor.fit(X, Y)  # perform linear regression
    Y_pred = linear_regressor.predict(X)  # make predictions

    m = reg.coef_
    b = reg.intercept_
    calc_shift = (Y - b)/m

    dif = np.abs(np.subtract(X, calc_shift))
    mae = pd.DataFrame(dif, columns=['MAE'])
    mae = mae['MAE']

    mae = mae[mae.between(mae.quantile(0), mae.quantile(.99))]
    print(dir,'\t',atom,'\t',np.average(mae), '\t', np.std(mae))

#    return np.average(mae)


if __name__ == '__main__':
#    parser = argparse.ArgumentParser(description='')
#    parser.add_argument('infile', type=argparse.FileType('r'), nargs='+', help='path to input .tsv file')

#    args = parser.parse_args()

    exp_files = glob.glob('exp/*Shift')
    calc_files = glob.glob('amber/opt/*Opt')

#    exp_files = args.infile_exp
#    calc_files = args.infile_calc

    # Make nested dictionary of exp, calc shield and calc shifts for all ten molecules for all methods ()

    d_h = {}
    d_c = {}
    for g in exp_files:
        for h in calc_files:
            if ((g).split('/')[-1]).split('_')[0] == ((h).split('/')[-1]).split('_')[0]:
                IDS = ((h).split('/')[-1]).split('_')[0] + "_" + ((h).split('/')[-1]).split('_')[1]
                exp = load_file(g)
                calc_shield = load_file(h)
                calc_shift = convert_to_shifts_CONS(calc_shield)

                d[IDS] = {"exp": exp, "calc_shield": calc_shield, "work": calc_shift}

    df = pd.DataFrame(data=d)
    for g in files:
        IDS = ((g).split('/')[-1]).split('.')[0]
        df_h, df_c = load_file(g)
        for i in list(df_h):
                METHOD = i
                exp = df[i]
                calc_shift = convert_to_shifts_CONS(calc_shield)

        d_h[IDS][METHOD] = {"exp": exp, "calc_shield": calc_shield, "work": calc_shift}
        d_c[IDS][METHOD] = 
    #FIRST GUESS

    for g in args.infile:
        df, atom_list, exp, methods = read_shifts(g)
        exp_c, exp_h = separate_exp(atom_list, exp, atom='C'), separate_exp(atom_list, exp, atom='H')
        for i in methods:
            meth = str(i)
            shift_c, shift_h, shield_c, shield_h  = convert_to_shifts(atom_list, df[meth])
            idx(g.name, exp_c, shift_c, shield_c, i, atm='C')
            idx(g.name, exp_h, shift_h, shield_h, i, atm='H')



for i in subdirs:
    mae = 0
    converge = 0.01
    if i.endswith('no') or i.endswith('opt'):
        linear(i,'H')
        linear(i,'C')




