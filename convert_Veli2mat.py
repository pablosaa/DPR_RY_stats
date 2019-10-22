import scipy.io
import numpy as np

dprns = np.load('../data/sr_pns.npy')
dprns_phase = np.load('../data/sr_DPRns_phase.npy')
dprans_phase = np.load('../data/sr_DPRans_phase.npy')

dprns_type  = np.load('../data/sr_type.npy')
dprans = np.load('../data/sr_pes.npy')
dprhip = np.load('../data/sr_hip.npy')
ry=np.load('../data/gr_pp_ipoli.npy')

scipy.io.savemat('../data/DPR_RADOLAN.mat',mdict={'DPRns':dprns,'DPRns_ph':dprns_phase,'DPRans_ph':dprans_phase,'DPRns_ty':dprns_type,'DPRans':dprans,'DPR_hip':dprhip,'RY':ry})
#scipy.io.savemat('./mamatdata.mat',mdict={'arr':arr,'rra':arr+1})

