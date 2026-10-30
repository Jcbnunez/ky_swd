#module load python3.12-anaconda/2024.06-1;conda activate moments_jcbn; python

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
idx = int(sys.argv[1])
print("now processing", idx)

###
indat = pd.read_csv("VA_KY_pairs.txt", header=0, delimiter="\t")

pop1 = indat.pop1[idx]
pop2 = indat.pop2[idx]

root_dadi = "dadi_objects_duos"

namsam = ["probs", pop1, pop2, "delim"]
paths = [root_dadi, ".".join(namsam)]
fs_file = "/".join(paths)

Pair_name = "_".join([pop1, pop2])

###
### collect metadata
root_meta = "./L_meta_objects_duos"

namsam2 = ["probs", pop1, pop2, "meta"]
paths2 = [root_meta, ".".join(namsam2)]
meta_file = "/".join(paths2)

inmet = pd.read_csv(meta_file, header=0, delimiter="\t")

pop1_n1 = inmet.an1_ne[0]
pop2_n2 = inmet.an2_ne[0]
L = inmet.L[0]

print("fs file is", fs_file)
print("Pair_name is", Pair_name)
print("names are", pop1, pop2)

# constants
mu = 2.8e-9      # Keightley et al. 2014
g = 0.06666667   # 15 generations/year
iterations = 100

print("loading data")

dd = dadi.Misc.make_data_dict(fs_file)
data = pd.read_csv(fs_file, sep="\t", nrows=1)

pop_id = ["A_pop1", "B_pop2"]
pools = [pop1_n1, pop2_n2]

print("folding sfs")

fs_folded = Spectrum.from_data_dict(
    dd,
    pop_ids=pop_id,
    projections=pools,
    polarized=False
)

ns = fs_folded.sample_sizes
S = fs_folded.S()

##########################################
# Output file
##########################################

PMmod = open(
    'split_nomig/%s_output.nomig.txt' % Pair_name,
    'w'
)

PMmod.write(
    "Pair_name\t"
    "fs_name\t"
    "L\t"
    "pop1_size\t"
    "pop2_size\t"
    "divergence_time\t"
    "theta\t"
    "nu1\t"
    "nu2\t"
    "Ts\t"
    "fs_sanitycheck\t"
    "-2LL_model\t"
    "AIC\n"
)

PMmod.close()

print("defining functions")

##########################################
# Strict isolation model
##########################################
def split_nomig_moments(params, ns, pop_ids=None):
    if pop_ids is not None and len(pop_ids) != 2:
        raise ValueError(
            "pop_ids must be a list of two population IDs"
        )
    nu1, nu2, T = params
    sts = moments.LinearSystem_1D.steady_state_1D(
        ns[0] + ns[1]
    )
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(
        fs,
        ns[0],
        ns[1]
    )
    # no migration after split
    fs.integrate(
        [nu1, nu2],
        T
    )
    fs.pop_ids = pop_ids
    return fs
func_moments = split_nomig_moments


##########################################
# parameter bounds
##########################################

upper_bound = [10, 10, 5]
lower_bound = [1e-5, 1e-5, 1e-5]

print("optimization loop")

for rep in range(iterations):
    prior = [
        np.random.uniform(lower_bound[x], upper_bound[x])
        for x in range(3)
    ]
    popt = moments.Inference.optimize_log(
        prior,
        fs_folded,
        func_moments,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        verbose=False,
        maxiter=100,
    )
    params = 3
    model = func_moments(popt, ns)
    ll_model = moments.Inference.ll_multinom(
        model,
        fs_folded
    )
    aic = 2 * params - 2 * ll_model
    print(
        "Maximum log composite likelihood: {0}".format(
            ll_model
        )
    )
    theta = moments.Inference.optimal_sfs_scaling(
        model,
        fs_folded
    )
    Nref = theta / (4 * mu * L)
    nu1 = popt[0]
    nu2 = popt[1]
    Ts = popt[2]
    divergence_time = 2 * Nref * Ts * g
    pop1_size = nu1 * Nref
    pop2_size = nu2 * Nref
    PMmod = open(
        'split_nomig/%s_output.nomig.txt' % Pair_name,
        'a'
    )
    PMmod.write(
        str(Pair_name) + '\t' +
        str(fs_file) + '\t' +
        str(L) + '\t' +
        str(pop1_size) + '\t' +
        str(pop2_size) + '\t' +
        str(divergence_time) + '\t' +
        str(theta) + '\t' +
        str(nu1) + '\t' +
        str(nu2) + '\t' +
        str(Ts) + '\t' +
        str(S) + '\t' +
        str(ll_model) + '\t' +
        str(aic) + '\n'
    )
    PMmod.close()

print("Moments finished running")
