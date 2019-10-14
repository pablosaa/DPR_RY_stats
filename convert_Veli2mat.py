import scipy.io
import numpy as np

dprns = np.load('sr_pns.npy')
dprns_phase = np.load('sr_phase.npy')
dprns_type  = np.load('sr_type.npy')
dprans = np.load('sr_pes.npy')
dprhip = np.load('sr_hip.npy')
ry=np.load('gr_pp_ipoli.npy')

scipy.io.savemat('./DPR_RADOLAN.mat',mdict={'DPRns':dprns,'DPRns_ph':dprns_phase,'DPRns_ty':dprns_type,'DPRans':dprans,'DPR_hip':dprhip,'RY':ry})
#scipy.io.savemat('./mamatdata.mat',mdict={'arr':arr,'rra':arr+1})

