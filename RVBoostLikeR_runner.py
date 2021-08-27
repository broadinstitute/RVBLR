#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os, sys, re

if sys.version_info[0] != 3:
    print("This script requires Python 3")
    exit(1)

import datetime
import subprocess
import argparse
import logging
FORMAT = "%(asctime)-15s: %(levelname)s %(module)s.%(name)s.%(funcName)s %(message)s"
logger = logging.getLogger('ctat_mutations')
logging.basicConfig(stream=sys.stderr, format=FORMAT, level=logging.INFO)


sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../../PyLib"]))
from Pipeliner import Pipeliner, Command, run_cmd, ParallelCommandList

RVB_UTILDIR = os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "util"])
CTAT_UTILDIR = os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../../src"])


def main():

    parser = argparse.ArgumentParser(description="wrapper for running rvboost-like-R", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--input_vcf", type=str, required=True, help="input vcf file")

    parser.add_argument("--work_dir", type=str, required=True, help="working directory name for intermediates")

    parser.add_argument("--attributes", type=str, required=False, help="vcf info attributes to use for scoring purposes",
                        default="DJ,ReadPosRankSum,QD,FS,ED,PctExtPos,RS")
    

    parser.add_argument("--score_threshold", type=float, required=False, default=0.05, help="score threshold for filtering rvboost results")


    parser.add_argument("--output_filename", required=True, help="name of output file containing final list of variants")
    
    args = parser.parse_args()
    
    
    ## prep for run
    output_dir = args.work_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    checkpts_dir = os.path.join(output_dir, "__chckpts")

    pipeliner = Pipeliner(checkpts_dir)
    
    ## build pipeline

    # extract feature matrix
    matrix_file = os.path.join(output_dir, "features.matrix")

    # ensure we capture the common variant annots
    if "RS" not in args.attributes.split(","):
        args.attributes += ",RS"
    
    
    cmd = " ".join([ os.path.join(CTAT_UTILDIR, "annotated_vcf_to_feature_matrix.py"),
                    "--vcf", args.input_vcf,
                    "--features", args.attributes,
                    "--output", matrix_file])
    pipeliner.add_commands([Command(cmd, "feature_extraction_to_matrix.ok")])
    
    
    # run rvboost-like-R
    boost_scores_file = matrix_file + ".RVBLR_var_scores"
    cmd = " ".join([ os.path.join(RVB_UTILDIR, "RVBoostLike.R"),
                    "--feature_matrix", matrix_file,
                    "--attributes", args.attributes,
                    "--boosting_score_threshold", str(args.score_threshold),
                    "--output", boost_scores_file])
        
    pipeliner.add_commands([Command(cmd, "rvboost_core.ok")])

        
    # generate filtered vcf file:
    
    cmd = " ".join([ os.path.join(CTAT_UTILDIR, "extract_boosted_vcf.py"),
                    "--vcf_in", args.input_vcf,
                    "--boosted_variants_matrix", boost_scores_file,
                    "--vcf_out", args.output_filename])

    pipeliner.add_commands([Command(cmd, "filt_qscore_{}.ok".format(args.score_threshold))])
    
    
    
    pipeliner.run()

    print("-done. See RVBLR vcf output file: {}".format(args.output_filename))
    
    sys.exit(0)
    


if __name__ == '__main__':
    main()
    
