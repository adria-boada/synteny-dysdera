#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# moments_Personal.py
#
# 25 de juny 2024  <adria@molevol-OptiPlex-9020>

help_msg = """
This is meant to be a general use script to run 'moments' to fit any model on an
AFS/joint-SFS with one to five populations. The user will have to edit information
about their allele frequency spectrum and provide a custom model.

The optimization routine runs a user-defined number of rounds, each with a
user-defined or predefined number of replicates. The starting parameters are
initially random, but after each round is complete the parameters of the best
scoring replicate from that round are used to generate perturbed starting
parameters for the replicates of the subsequent round. The arguments controlling
steps of the optimization algorithm (maxiter) and perturbation of starting
parameters (fold) can be supplied by the user for more control across rounds.
The user can also supply their own set of initial parameters, or set custom
bounds on the parameters (upper_bound and lower_bound) to meet specific model
needs. This flexibility should allow these scripts to be generally useful for
model-fitting with any data set.

This script must be in the same working directory as Optimize_Functions.py, which
contains all the functions necessary.
"""

import sys
#>import os
import numpy as np
import moments
#>import pylab
#>from datetime import datetime
import Optimize_Functions
import Models_2D

# AVAILABLE TESTING/FITTING MODELS
# --------------------------------
available_models = {
    "no_mig": {
        "func": Models_2D.no_mig,
        "param": 3,
        "labels": "nu1, nu2, T",
    },
    "sym_mig": {
        "func": Models_2D.sym_mig,
        "param": 4,
        "labels": "nu1, nu2, m, T"
    },
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if __name__ == '__main__':

    # Parse command line arguments.
    import argparse
    parser = argparse.ArgumentParser(
        description=help_msg,
        # make sure the 'help_msg' is not automatically
        # wrapped at 80 characters (manually assign newlines).
        formatter_class=argparse.RawTextHelpFormatter,
    )
    # REQUIRED PARAMETERS
    # -------------------
    # Fitxer SNP partint de VCF segons el format tipificat als exemples de
    # 'moments_Pipeline'.
    parser.add_argument("--snps_file", required=True,
                        type=argparse.FileType("r"),
                        help="The SNPs file.")
    # 'pop_ids' is a list which should match the populations
    # headers of your SNPs file columns
    parser.add_argument("--pop_ids", required=True, nargs="+",
                        metavar="ID", type=str,
                        help="Population IDs (matching those within the "+
                        "SNPs file).")
    # Projection sizes, not in individuals but in ALLELES.
    parser.add_argument("--projections", metavar="SIZES",
                        required=True, nargs="+", type=int,
                        help="Projection sizes in ALLELES instead "+
                        "of individuals.")
    parser.add_argument("--rounds", required=True, type=int,
                        help="Amount of rounds of optimization.")
    parser.add_argument("--replicates", required=True,
                        nargs="+", type=int,
                        help="Replicates for each round; must be a list "+
                        "of ints of length 'ROUNDS'. Separate values "+
                        "with blank space.")
    parser.add_argument("--fit_models", required=True, nargs="+",
                        choices=["no_mig", "sym_mig"],
                        help="The name of the model we want to fit.")

    # OPTIONAL PARAMETERS
    # -------------------
    parser.add_argument("--polarized", action="store_const", const=True,
                        default=False,
                        help="Flag for polarized SFS input (default: "+
                        "not polarized with outgroup).")
    parser.add_argument("--independent_runs", required=False,
                        default=1, type=int,
                        help="The number of independent runs of optimization "+
                        "performed.")
#>    parser.add_argument("--upper_bounds", required=False, default=None,
#>                        nargs="+", type=int,
#>                        help="A list with the upper bounds for the model's "+
#>                        "parameters; must be a list of ints of length "+
#>                        "'MODEL PARAMETERS'. Separate values with blank space "+
#>                        "(default: parameters unbound).")
#>    parser.add_argument("--lower_bounds", required=False, default=None,
#>                        nargs="+", type=int,
#>                        help="A list with the lower bounds for the model's "+
#>                        "parameters; must be a list of ints of length "+
#>                        "'MODEL PARAMETERS'. Separate values with blank space "+
#>                        "(default: parameters unbound).")
    parser.add_argument("--max_iters", required=False, type=int, nargs="+",
                        default=None,
                        help="Steps of the optimization algorithm "+
                        "(default: The reverse of a geometric series "+
                        "starting with '50').")
    parser.add_argument("--fold", required=False, type=int, nargs="+",
                        default=None,
                        help="The perturbation of the starting parameters; "+
                        "must be a list of ints of length 'ROUNDS'. "+
                        "Separate values with a blank space "+
                        "(default: from 'ROUNDS' to '1').")

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.

    # CHECK COMMAND LINE PARAMETERS
    # -----------------------------
    if len(args.replicates) != int(args.rounds):
        raise Exception(
            "The length of the list 'REPLICATES' " +
            "({}) ".format(len(args.replicates)) +
            "should be the number of rounds " +
            "({}).".format(int(args.rounds)))
    if args.max_iters:
        max_iters = args.max_iters
        if len(max_iters) != int(args.rounds):
            raise Exception(
                "The length of the list 'MAX_ITERS' " +
                "({}) ".format(len(max_iters)) +
                "should be the number of rounds " +
                "({}).".format(int(args.rounds)))
    else:
        # Setup default max_iters:
        max_iters = [50]
        while len(max_iters) < int(args.rounds):
            max_iters.append(np.ceil(
                # Geometric series starting from 50.
                max_iters[0] * ((1/2) ** len(max_iters))
            ))
        max_iters = max_iters[::-1]
    if args.fold:
        algor_fold = args.fold
        if len(algor_fold) != int(args.rounds):
            raise Exception(
                "The length of the list 'FOLD' " +
                "({}) ".format(len(algor_fold)) +
                "should be the number of rounds " +
                "({}).".format(int(args.rounds)))
    else:
        # Setup default fold:
        algor_fold = list(range(int(args.rounds), 0, -1))

    # READ SFS FILE AND PRINT SUMMARY
    # -------------------------------
    # Avoid having the file open in memory through 'argparse'.
    args.snps_file.close()
    # Start by importing the data required for the joint-site freq. spectrum.
    # Create a Python dictionary from the SNPs file.
    snps_data_dict = moments.Misc.make_data_dict(args.snps_file.name)

    # Convert this dictionary into an AFS object.
    # [polarized = False] creates folded spectrum object.
    fs = moments.Spectrum.from_data_dict(
        snps_data_dict,
        pop_ids=args.pop_ids,
        projections=args.projections,
        polarized=args.polarized)

    # Print useful information about the AFS or joint-site FS.
    print("= = = = = = = = = = = = = = = = =")
    print("Overview of the SFS file... '{}'".format(args.snps_file.name))
    print("Projection: {}".format(args.projections))
    print("Sample sizes: {}".format(fs.sample_sizes))
    # Sum of the sites rounded to two decimal places.
    print("Sum of sites in the SFS: {}".format(np.around(fs.S(), 2)))
    print("Steps of the optimization algorithm, per round: " +
          str(max_iters))
    print("Perturbation of the starting parameters, per round: " +
          str(algor_fold))
    print("= = = = = = = = = = = = = = = = =")

    # START OPTIMIZATION ALGORITHM
    # ----------------------------
    for independent_run in range(args.independent_runs):
        prefix = "Exec{}".format(independent_run)
        for mod in args.fit_models:

            Optimize_Functions.Optimize_Routine(
                fs=fs,
                outfile=str(prefix),
                model_name="{0}_{1}".format(mod, prefix),
                func=available_models[mod]["func"],
                rounds=int(args.rounds),
                param_number=int(available_models[mod]["param"]),
                fs_folded=not bool(args.polarized),
                reps=args.replicates,
                maxiters=max_iters,
                folds=algor_fold,
                param_labels=available_models[mod]["labels"],
            )

