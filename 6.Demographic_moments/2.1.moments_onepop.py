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
import numpy as numpy
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
#opening output file to give column names
PMmod=open('null/%s_output.null.txt' % Pair_name,'w')
PMmod.write(
            str("Pair_name")+'\t'+ #print pair name
            str("fs_name")+'\t'+ #double checking fs_lines[y] is working as I want it to
            str("L")+'\t'+ #double checking L is working as I want it to
            str("theta")+'\t'+ #double checking L is working as I want it to
            str("nuB")+'\t'+ #double checking L is working as I want it to
            str("nuF")+'\t'+ #double checking L is working as I want it to
            str("T")+'\t'+ #double checking L is working as I want it to
            str("fs_sanitycheck")+'\t'+
            str("-2LL_model")+'\t'+
            str("AIC")+'\n')
PMmod.close()


print('defining functions')
####
def bottlegrowth_split_mig(params, ns, pop_ids=None):
	 """
	 params = (nuB, nuF, m, T, Ts)
	 ns = [n1, n2]

	 Instantanous size change followed by exponential growth then split with
	 migration.

	 - nuB: Ratio of population size after instantanous change to ancient
		population size
	 - nuF: Ratio of contempoary to ancient population size
	 - m: Migration rate between the two populations (2*Na*m).
	 - T: Time in the past at which instantaneous change happened and growth began
		(in units of 2*Na generations)
	 - Ts: Time in the past at which the two populations split.
	 - n1, n2: Sample sizes of resulting Spectrum.

	 :param ns: List of population sizes in first and second populations.
	 :param pop_ids: List of population IDs.
	 """
	 if pop_ids is not None and len(pop_ids) != 2:
		  raise ValueError("pop_ids must be a list of two population IDs")
	 nuB, nuF, m, T, Ts = params
	 nu_func = lambda t: [nuB * numpy.exp(numpy.log(nuF / nuB) * t / T)]
	 sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	 fs = moments.Spectrum(sts)
	 fs.integrate(nu_func, T - Ts, dt_fac=0.01)
	 # we split the population
	 fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	 nu0 = nu_func(T - Ts)[0]
	 nu_func = lambda t: 2 * [nu0 * numpy.exp(numpy.log(nuF / nu0) * t / Ts)]
	 fs.integrate(nu_func, Ts, m=numpy.array([[0, m], [m, 0]]))
	 fs.pop_ids = pop_ids
	 return fs


def bottlegrowth(params, ns, pop_ids=None):
	 """
	 params = (nuB, nuF, T)

	 ns = [n1, n2]

	 Instantanous size change followed by exponential growth with no population
	 split.

	 - nuB: Ratio of population size after instantanous change to ancient
		population size
	 - nuF: Ratio of contempoary to ancient population size
	 - T: Time in the past at which instantaneous change happened and growth began
		(in units of 2*Na generations)
	 - n1, n2: Sample sizes of resulting Spectrum.

	 :param params: List of parameters, (nuB, nuF, T).
	 :param ns: List of population sizes in first and second populations.
	 :param pop_ids: List of population IDs.
	 """
	 if pop_ids is not None and len(pop_ids) != 2:
		  raise ValueError("pop_ids must be a list of two population IDs")
	 nuB, nuF, T = params
	 return bottlegrowth_split_mig((nuB, nuF, 0, T, 0), ns, pop_ids=pop_ids)

func_moments=bottlegrowth

upper_bound = [10, 10, 5]
lower_bound = [1e-5, 1e-5, 1e-5, 0]


for i in range(int(iterations)): #iterations is imported from sys. argument #1
	print("starting optimization "+str(i))
	#Parameters are set above " nu1, nu2, T, m12, m21"
	#This number is 5. i.e., count parameters. 
	params = len(["nuB", "nuF", "T"]) #for use in AIC calculation
	#Start the run by picking random parameters from a uniform distribution.
	prior=[numpy.random.uniform(lower_bound[x],upper_bound[x]) for x in range(params)]
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
	Nref= theta/(4*mu*L)
	nuB=popt[0]
	nuF=popt[1]
	T=popt[2]
	#Now calculate Ts from Model fit	#Open the output file
	PMmod=open('null/%s_output.null.txt' % Pair_name,'a')
	#Dumping output ot outfile
	PMmod.write(
	str(Pair_name)+'\t'+ #print pair name
	str(fs_file)+'\t'+ #double checking fs is the right one
	str(L)+'\t'+ #double checking L is working as desired
	str(theta)+'\t'+
	str(nuB)+'\t'+
	str(nuF)+'\t'+
	str(T)+'\t'+
	str(S)+'\t'+ #sanity check... should give number of segregating sites in SFS
	str(ll_model)+'\t'+
	str(aic)+'\n')
	PMmod.close()

print("Moments finished running")