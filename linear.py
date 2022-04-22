import numpy as np
import matplotlib.pyplot as plt  # To visualize
import pandas as pd  # To read data
from sklearn.linear_model import LinearRegression
import os

import pandas as pd
import glob

import cost_matrix

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

    r_2 = reg.score(X,Y)
    m = reg.coef_
    b = reg.intercept_


#    plt.scatter(X, Y)
#    plt.plot(X, Y_pred, color='red')
#    nme = dir.split('/')[1] + '_' + atom
#    plt.title(nme)
#    plt.show()

    calc_shift = (Y - b)/m

    dif = np.abs(np.subtract(X, calc_shift))
    mae = pd.DataFrame(dif, columns=['MAE'])
    mae = mae['MAE']



#    print(len(mae))
#    mae = mae[mae.between(mae.quantile(0), mae.quantile(.99))]
#    print(len(mae))
#    print(dir,atom)
    print(dir,'\t',atom,'\t',np.average(mae), '\t', np.std(mae))
#    print(np.std(mae))
#    print(np.max(mae))
#    print(np.min(mae))

#    print(m[0][0])
#    print(b[0])
#    print(r_2)

#subdirs = [x[0] for x in os.walk('.')]

d = '.'
subdirs = [os.path.join(d, o) for o in os.listdir(d) if os.path.isdir(os.path.join(d,o))]
print(subdirs)



for i in subdirs:
    mae = 0
    converge = 0.01
    if i.endswith('no') or i.endswith('opt'):
        linear(i,'H')
        linear(i,'C')

