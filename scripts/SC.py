import allel
import moments, sys, os, matplotlib, numpy as np
from moments import Misc,Spectrum,Numerics,Manips,Integration,Demographics1D,Demographics2D
from matplotlib import pyplot as plt
os.chdir("/home/gese8011/elinck_bin/demo_inference")

# load sfs
fs=moments.Spectrum.from_file("allopatric_parental_unlinked_sfs.txt")
fs=fs.fold()
fs.pop_ids("m")
ns=fs.sample_sizes
projections = [16,16]
fs.mask[1,:]=True #mask singletons in all populations
fs.mask[:,1]=True

def SC(params, ns):
    s,nPer,nAlb,t1,t2,m12,m21 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    nu1_func = lambda t: s*(nPer/s)**(t/t2)
    nu2_func = lambda t: (1-s)*(nAlb/(1-s))**(t/t2)
    nu_func = lambda t: [nu1_func(t), nu2_func(t)]
    fs.integrate([nPer,nAlb], t1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs.integrate(nu_func, t2, dt_fac=0.01, m=np.array([[0, m12] , [m21, 0]]))
    return  fs

upper_bound = [1,10,10,10,10,10,10]
lower_bound = [1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4]

for i in range(5):
        poptg=[np.random.uniform(lower_bound[x],upper_bound[x]) for x in range(7)]
        poptg=moments.Inference.optimize_log(poptg, fs, SC,lower_bound=lower_bound,
                                             upper_bound=upper_bound,verbose=True,
                                             maxiter=5)
        model=SC(poptg, ns)
        ll_model=moments.Inference.ll_multinom(model,fs)
        theta = moments.Inference.optimal_sfs_scaling(model, fs)
        L=950000000*(549501/37526336)
        Nref=theta/(4*2.3e-9*L)
        nPer=poptg[1]*Nref
        nAlb=poptg[2]*Nref
        t1=poptg[3]*2*Nref
        t2=poptg[4]*2*Nref
        m12=poptg[5]/(2*Nref)
        m21=poptg[6]/(2*Nref)
        out1=[nPer,nAlb,t1,t2,m12,m21,ll_model,theta]
        out1="\t".join(map(str,out1))+"\n"
        f=open("/home/gese8011/elinck_bin/demo_inference/SC_realparams.txt",'a')
        f.write(out1)
        f.close()
        out2="\t".join(map(str,poptg))+"\t"+str(ll_model)+"\t"+str(theta)+"\n"
        f=open("/home/gese8011/elinck_bin/demo_inference/SC_modelparams.txt",'a')
        f.write(out2)
        f.close()
