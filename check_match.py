from sklearn.cluster import SpectralClustering as CLUSTER
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import json
from sklearn.linear_model import LinearRegression
import pandas as pd
import math

files = glob.glob('calc/json/*/*/*json')

at = 'c'
err = 10
method = 'munkres'

c_upper_bound = 67.13807584744352

exp_good = []
calc_good = []
exp_bad = []
calc_bad = []
file_good = []
file_bad = []

def distance(x1, y1, a, b, c):

    d = abs((a * x1 + b*y1 + c)) / (math.sqrt(a * a + b*b))
    return d


def linear(X, Y):

    X = np.array(X).reshape((-1,1))
    Y = np.array(Y)
    linear_regressor = LinearRegression()  # create object for the class
    reg = linear_regressor.fit(X, Y)  # perform linear regression

    m = reg.coef_
    b = reg.intercept_
 
    Y_pred = linear_regressor.predict(X)
    return m,b,Y_pred

def collect_dots(FILE):

  with open(FILE) as f:
    js = json.load(f)

  if js[at+'_nmr'][method+'_mean_average_error'] <= c_upper_bound:
    for atom in js[at+'_nmr']['spectrum']:
          if atom['exp_shift_'+method+'_match'] is not None:
              shift.append(atom['exp_shift_'+method+'_match'])
              shield.append(atom[at+'_shielding'])
              pred_shift.append(atom['predicted_'+at+'_shift'])
              X.append([atom['exp_shift_'+method+'_match'],atom['predicted_'+at+'_shift']])
  else:
    pass


def sort(FILE, m,y, b):
  with open(FILE) as f:
    js = json.load(f)
  for atom in js[at+'_nmr']['spectrum']:
      if atom['exp_shift_'+method+'_match'] is not None:
        if distance(atom['exp_shift_'+method+'_match'],atom[at+'_shielding'], m,y, b) < err:
          exp_good.append(atom['exp_shift_'+method+'_match'])
          calc_good.append(atom['predicted_'+at+'_shift'])
          file_good.append(FILE)
        else:
          exp_bad.append(atom['exp_shift_'+method+'_match'])
          calc_bad.append(atom['predicted_'+at+'_shift'])
          file_bad.append(FILE)

def plot(shield,shift):
  plt.plot(shield, shift, 'o')
  plt.show()

shift = []
shield =[]
X=[]
pred_shift = []
big_error = []
reg_error = []

for fi in files:
  try:
    collect_dots(fi)
  except:
    continue

#m1, b1, pred_shift = linear(shift, shield)

#m1 = -1.0931
#b1 = 32.031
m1 = -1.0175
b1 = 185.65

for fi in files:
  try:
    sort(fi, -1*m1, 1, -1*b1)
  except:
    continue

dict = {'experimental': shift, 'calculated':pred_shift}
df = pd.DataFrame(dict)
df.dropna()

plt.scatter(exp_bad, calc_bad, color='red')
plt.scatter(exp_good, calc_good, color='green')
#sns.jointplot(data=df,x="experimental", y='calculated', kind="hex")
plt.show()

error = []

for i in range(len(calc_good)):
  if shift[i] != None:
     delta = np.abs(calc_good[i] - exp_good[i])
     error.append(delta)
     
mae = pd.DataFrame(error, columns=['MAE'])
mae = mae['MAE']

fil = []
for i in file_good:
 if i not in fil:
  if i not in file_bad:
   fil.append(i)

print(len(mae))
#mae = mae[mae.between(mae.quantile(0), mae.quantile(.99))]
print(np.average(mae), '\t', np.std(mae), '\t', len(error), '\t')

#er = {'big errors': big_error}
#df = pid.DataFrame(er)
#df.to_csv('c_errors.csv')

#dict = {'file': file_bad, 'exp': exp_bad, 'calc': calc_bad}
#df2 = pd.DataFrame(dict)
#df2.to_csv('switches_'+method+'.csv')

