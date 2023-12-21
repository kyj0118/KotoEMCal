# for data loading
import uproot
import numpy as np
from scipy.sparse import coo_matrix

# for training
import xgboost
from sklearn.multioutput import MultiOutputRegressor
from sklearn.model_selection import train_test_split

# for memory & cpu usage check 
import psutil 
import os 

# etc
import sys
import time

def memchk(message: str = 'debug'):
    pid = os.getpid()
    py  = psutil.Process(pid)
    rss = py.memory_info().rss / 2. **20
    print(f"[{message}] memory usage: {rss: 10.5f} MB")
    return rss


time_start = time.time()

model_PATH = 'model'
model_name_x = '{}/model0.json'.format(model_PATH)
model_name_y = '{}/model1.json'.format(model_PATH)

# Define the branches you want to read
branches = ["nScintHit", "ScintHit.e", "ScintHit.ModuleID","nCsIHit", "CsIHit.e", "CsIHit.CellID", "PrimaryParticle.px", "PrimaryParticle.py", "PrimaryParticle.pz"]

# List of ROOT files
#root_files = [f"../Geant4FiberW/root/electron1000MeV_{i:04d}.root:tree" for i in range(1, 21)]
root_files = [f"../Geant4FiberW/root/electron1000MeV_step_{i:04d}.root:tree" for i in range(1, 21)]

memchk('Before data loading')

# Use uproot to concatenate the data from all files
data = uproot.concatenate(root_files, branches, library="np")

nMaxScint=24*16
nMaxCsI=25

# Now data is a dictionary where keys are branch names and values are concatenated arrays
nScintHit = data["nScintHit"]
ScintHit_e = data["ScintHit.e"]
ScintHit_ModuleID = data["ScintHit.ModuleID"]

nCsIHit = data["nCsIHit"]
CsI_e = data["CsIHit.e"]
CsI_CellID = data["CsIHit.CellID"]

px = data["PrimaryParticle.px"]
py = data["PrimaryParticle.py"]
pz = data["PrimaryParticle.pz"]

# Initialize lists to store data for COO matrix (features) and targets
data_values = []
targets = []

total_channels = nMaxScint + nMaxCsI
feature_matrix = np.zeros((len(nScintHit), total_channels))

# Process the concatenated data
for i in range(len(nScintHit)):
    # Scintillator features
    for j in range(nScintHit[i]):
        channel_index = ScintHit_ModuleID[i][j]
        feature_matrix[i, channel_index] = ScintHit_e[i][j]

    # CsI features
    for j in range(nCsIHit[i]):
        channel_index = CsI_CellID[i][j] + nMaxScint
        feature_matrix[i, channel_index] = CsI_e[i][j]


    # Calculate px/pz and py/pz for the current event
    px_pz = px[i] / pz[i]
    py_pz = py[i] / pz[i]

    # Append the targets
    targets.append([px_pz, py_pz])

# Convert targets list to numpy array
targets = np.array(targets)

# Hyperparameters
ne = int(300)
md = int(100)
ss = float(0.8)
lr = float(0.02)
gamma = float(0.0)

modelx = xgboost.XGBRegressor(n_estimators=ne, learning_rate=lr, gamma=gamma, subsample=ss,colsample_bytree=1, max_depth=md)
modely = xgboost.XGBRegressor(n_estimators=ne, learning_rate=lr, gamma=gamma, subsample=ss,colsample_bytree=1, max_depth=md)

modelx.load_model(model_name_x)
modely.load_model(model_name_y)


predictionx = modelx.predict( feature_matrix )
predictiony = modely.predict( feature_matrix )

foutname='result.txt'
fout = open(foutname,'w')
for i in range(0, predictionx.size):
    fout.write('{} {} {} {}\n'.format(predictionx[i],predictiony[i], targets[i,0], targets[i,1]))
fout.close()

