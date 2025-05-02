#module load python3.12-anaconda/2024.06-1;conda activate moments_jcbn; python

# import packages that'll be used
import moments
from moments import Numerics
from moments import Integration
from moments import Spectrum
import dadi
from dadi import Misc
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
from datetime import datetime

###
i=int(sys.argv[1])
print("now processing", i)
###

indat = pd.read_csv("VA_KY_pairs.txt",header=0,delimiter="\t")
pop1=indat.pop1[i]
pop2=indat.pop2[i]
root_dadi="dadi_objects_duos"
namsam = ["probs",pop1, pop2,"delim"]
paths = [root_dadi,".".join(namsam)]
"/".join(paths)
fs_file = "/".join(paths)
###
Pair_name = "_".join([pop1, pop2])
###
### collect the metadata
root_meta="./L_meta_objects_duos"
namsam2 = ["probs",pop1, pop2,"meta"]
paths2 = [root_meta,".".join(namsam2)]
"/".join(paths2)
meta_file = "/".join(paths2)
inmet = pd.read_csv(meta_file,header=0,delimiter="\t")
pop1_n1=inmet.an1_ne[0]
pop2_n2=inmet.an2_ne[0]
L = inmet.L[0]

print("fs file is", fs_file )
print("Pair_name is", Pair_name )
print("names are", pop1,pop2 )

#constants
mu = 2.8e-9 #from Keightley et al. 2014
g = 0.06666667 #equals 15 gen/year --- See Pool 2015
iterations=100
#####

####
print('loading data')
dd = dadi.Misc.make_data_dict(fs_file) #reads in genomalicious SNP file
data = pd.read_csv(fs_file, sep="\t", nrows=1)

#####
#setting pop id's and projections from if/else
pop_id=["A_pop1",  "B_pop2"]
pools=[pop1_n1, pop2_n2 ] 

print('folding sfs')
fs_folded = Spectrum.from_data_dict(dd, pop_ids=pop_id, projections=pools, polarized=False) #takes data dict and folds
ns = fs_folded.sample_sizes #gets sample sizes of dataset
S = fs_folded.S()

####
#opening output file to give column names
PMmod=open('split_migasym/%s_output.migasymm.txt' % Pair_name,'w')
PMmod.write(
            str("Pair_name")+'\t'+ #print pair name
            str("fs_name")+'\t'+ #double checking fs_lines[y] is working as I want it to
            str("L")+'\t'+ #double checking L is working as I want it to
            str("pop1_size")+'\t'+ #nu1
            str("pop2_size")+'\t'+ #nu2
            str("divergence_time")+'\t'+ #divergence T
            str("mig_pop1")+'\t'+ #Migration ij
            str("mig_pop2")+'\t'+ #Migration ji
            str("Mij")+'\t'+
            str("Mji")+'\t'+
            str("theta")+'\t'+
            str("nu1")+'\t'+
            str("nu2")+'\t'+
            str("Ts")+'\t'+
            str("m12")+'\t'+
            str("m21")+'\t'+
            str("fs_sanitycheck")+'\t'+
            str("-2LL_model")+'\t'+
            str("AIC")+'\n')
PMmod.close()

####
#opening output file to give column names
#model with symmetric migration from Moments bitbucket
def split_migasym_moments(params, ns, pop_ids=None):
    if pop_ids is not None and len(pop_ids) != 2:
        raise ValueError("pop_ids must be a list of two population IDs")
    nu1, nu2, T, m12, m21 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T, m=np.array([[0, m12], [m21, 0]]))
    fs.pop_ids = pop_ids
    return fs

func_moments = split_migasym_moments

#boundaries
#bounded to same standards as theta param script
upper_bound = [10, 10, 5, 50, 50]
lower_bound = [1e-5, 1e-5, 1e-5, 0, 0]

# run X optimizations from uniform-sampled starting params
for i in range(int(iterations)): #iterations is imported from sys. argument #1
	print("starting optimization "+str(i))
	#Parameters are set above " nu1, nu2, T, m12, m21"
	#This number is 5. i.e., count parameters. 
	params = len(["nu1", "nu2", "Ts", "m12", "m21"]) #for use in AIC calculation
	#Start the run by picking random parameters from a uniform distribution.
	prior=[np.random.uniform(lower_bound[x],upper_bound[x]) for x in range(params)]
	#This is the optimization step for moments.
	#popt is the prior.
	#fs folded is a tranform SFS by folding it. The original SFS loaded in sys.arg #1 is quasi-folded. i.e. polarized to reference genome.
	#Folding it is a must, because reference is not 100% ancestral
	popt=moments.Inference.optimize_log(prior, fs_folded, func_moments,
	lower_bound=lower_bound, upper_bound=upper_bound,
	verbose=False, maxiter=100,
	)
	#This is the moments function.
	model = func_moments(popt, ns)
	#Calculate log likelihood of the model fit
	ll_model=moments.Inference.ll_multinom(model, fs_folded)
	#Now calculate AIC of model fit
	aic = 2*params - 2*ll_model
	print('Maximum log composite likelihood: {0}'.format(ll_model))
	#Now estimate theta from model fit
	theta = moments.Inference.optimal_sfs_scaling(model, fs_folded)
	#reducing complexity of calculations to follow by adding variables in lieu of expressions/esoteric df calls
	Nref= theta/(4*mu*L)
	nu1=popt[0]
	nu2=popt[1]
	Ts=popt[2]
	m12=popt[3]
	m21=popt[4]
	#Now calculate Ts from Model fit
	divergence_time = 2*Nref*Ts*g #calculates divergence time in years
	#Now calculate Migration rate (fraction of migrants that move between pops)
	Mij = m12/(2*Nref) #actual migration rate for populations from 2 into 1
	Mji = m21/(2*Nref) #actual migration rate for populations from 1 into 2
	#Now we are estimated the nominal migrant count based on Mij
	mig_pop1 = Mij*(nu1*Nref) #number of individuals going i to j: migrants = Mij*nu1*Nref
	mig_pop2 = Mji*(nu2*Nref) #number of individuals going j to i: migrants = Mij*nu2*Nref
	#Now estimate population size
	pop1_size = nu1*Nref #pop1 size
	pop2_size = nu2*Nref #pop2 size
	#Open the output file
	PMmod=open('split_migasym/%s_output.migasymm.txt' % Pair_name,'a')
	#Dumping output ot outfile
	PMmod.write(
	str(Pair_name)+'\t'+ #print pair name
	str(fs_file)+'\t'+ #double checking fs is the right one
	str(L)+'\t'+ #double checking L is working as desired
	str(pop1_size)+'\t'+ #nu1
	str(pop2_size)+'\t'+ #nu2
	str(divergence_time)+'\t'+ #divergence T
	str(mig_pop1)+'\t'+ #Migrants ij
	str(mig_pop2)+'\t'+ #Migrants ji
	str(Mij)+'\t'+ #migration rate ij
	str(Mji)+'\t'+ #migration rate ji
	str(theta)+'\t'+
	str(nu1)+'\t'+
	str(nu2)+'\t'+
	str(Ts)+'\t'+
	str(m12)+'\t'+
	str(m21)+'\t'+
	str(S)+'\t'+ #sanity check... should give number of segregating sites in SFS
	str(ll_model)+'\t'+
	str(aic)+'\n')
	PMmod.close()

print("Moments finished running")