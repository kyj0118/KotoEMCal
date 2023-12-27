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


if (len(sys.argv)) is not 5:
    print ('usage: python train.py [1. targetIndex(int)] [2. energy] [3. xy position] [4. particle]')
    exit()
targetIndex = int(sys.argv[1])
particleEnergy = int(sys.argv[2])
genPosition = str(sys.argv[3])
particle = str(sys.argv[4])

model_PATH = 'model'
memlog_PATH = 'memory_log'
timelog_PATH = 'time_log'

os.system('mkdir -p {}'.format(model_PATH))
os.system('mkdir -p {}'.format(memlog_PATH))
os.system('mkdir -p {}'.format(timelog_PATH))

model_name = f"{model_PATH}/{particle}{particleEnergy}MeV_position{genPosition}_{targetIndex}.json"
model_name_tmp = f"{model_PATH}/{particle}{particleEnergy}MeV_position{genPosition}_{targetIndex}_tmp.json"

flog_time=f'{timelog_PATH}/time_{particle}{particleEnergy}MeV_position{genPosition}_{targetIndex}.txt'
flog_memory=f'{memlog_PATH}/memory_{particle}{particleEnergy}MeV_position{genPosition}_{targetIndex}.txt'


# Define the branches you want to read
branches = ["nScintHit", "ScintHit.e", "ScintHit.ModuleID","nCsIHit", "CsIHit.e", "CsIHit.CellID", "PrimaryParticle.px", "PrimaryParticle.py", "PrimaryParticle.pz"]

# List of ROOT files
#root_files = [f"../Geant4FiberW/root/electron1000MeV_{i:04d}.root:tree" for i in range(1, 21)]
fname = f"../Geant4FiberW/root/{particle}{particleEnergy}MeV_position{genPosition}"
root_files = [f"{fname}_{i:04d}.root:tree" for i in range(1, 21)]

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
row_indices = []
col_indices = []
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

# save model every 20 tree boosting
nestimators_loop = int(20)

print("file read time :", time.time() - time_start )
time_before_fit = time.time()
rss_dataloading = memchk('After data loading')

xgb = xgboost.XGBRegressor(n_estimators=nestimators_loop, learning_rate=lr, gamma=gamma, subsample=ss,colsample_bytree=1, max_depth=md,verbosity=3,n_jobs=1)

ntrees=0
if (os.path.isfile(model_name) == True):
    xgb.load_model(model_name)
    ntrees=xgb.get_booster().num_boosted_rounds()

ntrainiter=int((ne-ntrees)/nestimators_loop)
if ntrainiter==0 :
    print ('ntrees : {} , n iteration : {}'.format(ntrees,ntrainiter))
    exit()
    
print ('ntrees : {} , n iteration : {}'.format(ntrees,ntrainiter))

for _ in range(ntrainiter):
    if ntrees == 0 :
        xgb.fit(feature_matrix, targets[:,targetIndex])
    else :
        xgb.fit(feature_matrix, targets[:,targetIndex],xgb_model=model_name)
        
    xgb.save_model(model_name_tmp)
    os.system('mv -f {} {}'.format(model_name_tmp,model_name))
    
    ntrees=xgb.get_booster().num_boosted_rounds()
    print ('ntrees : ', ntrees)
    rss_training = memchk('After training')
    fout_memory = open(flog_memory,'a')
    fout_memory.write('{0} {1:.1f}\n'.format(ntrees,rss_training))
    fout_memory.close()
    
    training_time = time.time() -  time_before_fit
    fout_time = open(flog_time,'a')
    fout_time.write('{0} {1:.1f}\n'.format(ntrees,training_time))
    fout_time.close()
