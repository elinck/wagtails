import moments, sys, os, matplotlib, numpy as np
from moments import Misc,Spectrum,Numerics,Manips,Integration,Demographics1D,Demographics2D
from matplotlib import pyplot as plt
import pandas as pd
os.chdir("/media/burke/bigMac/ethan/moments_revision/")

# set print options
np.set_printoptions(precision=3)

# build dd
fs=moments.Spectrum.from_file("allopatric_parental_sfs.txt")
fs=fs.fold()
fs.pop_ids("m")
ns=fs.sample_sizes
projections = [16,16]
fs.mask[1,:]=True #mask singletons in all populations
fs.mask[:,1]=True
fs=fs.sample()
ns=fs.sample_sizes

def IM(params, ns):
    s,nPer,nAlb,tSplit,m12,m21 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    nu1_func = lambda t: s*(nPer/s)**(t/tSplit)
    nu2_func = lambda t: (1-s)*(nAlb/(1-s))**(t/tSplit)
    nu_func = lambda t: [nu1_func(t), nu2_func(t)]
    fs.integrate(nu_func, tSplit, dt_fac=0.01,
                m=np.array([[0 , m12 ] , [ m21 , 0]]))
    return  fs

upper_bound = [1,10,10,10,10,10]
lower_bound = [1e-4,1e-4,1e-4,1e-4,1e-4,1e-4]

top_params=pd.read_csv("/home/gese8011/elinck_bin/demo_inference/IM_modelparams.txt",sep="\t",header=None)
top_params=top_params.sort_values(6,ascending=False)
poptg=np.array(top_params)[0,:]
poptg=poptg[range(8)]
poptg=moments.Inference.optimize(poptg, fs, IM, lower_bound=lower_bound, upper_bound=upper_bound, verbose=True, maxiter=1)

model=IM(poptg, ns)
ll_model=moments.Inference.ll_multinom(model,fs)
theta = moments.Inference.optimal_sfs_scaling(model, fs)
L=9.89e8*(14908552/37526336)
Nref=theta/(4*2.3e-9*L)
nPer=poptg[1]*Nref
nAlb=poptg[2]*Nref
tSplit=poptg[3]*2*Nref
m12=(poptg[4]/(2*Nref)*nPer)
m21=(poptg[5]/(2*Nref)*nAlb)
out1=[nPer,nAlb,tSplit,m12,m21,ll_model,theta]
out1="\t".join(map(str,out1))+"\n"
f=open("/home/gese8011/elinck_bin/demo_inference/IM_realparams_boots.txt",'a')
f.write(out1)
f.close()
out2="\t".join(map(str,poptg))+"\t"+str(ll_model)+"\t"+str(theta)+"\n"
f=open("/home/gese8011/elinck_bin/demo_inference/IM_modelparams_boots.txt",'a')
f.write(out2)
f.close()

#can run this script in parallel with eg
#parallel -N0 -j 25 python IM_bootstraps.py ::: {1..25}
