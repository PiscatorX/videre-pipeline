#!/usr/bin/env python

import glob
import argparse

class  MetaHitPairs(object):

    def  __init__(self):
        
        parser = argparse.ArgumentParser("wrapper for megahit grouping forward and left reads")
        parser.add_argument('-c','--cmd',  required=False, help = "full command for metahit") 
        parser.add_argument('-f','--fwd_glob', required=False, help = "wildcard for forward reads")
        parser.add_argument('-r','--rev_glob', required=False, help = "wildcard for reverse reads")
        parser.add_argument('-c','--cmd', help = "commands for megahit")
        args = parser.parse_args()
        self.args = args

    def megahit(self):


        print(args)




megahit =  MetaHitPairs()
megahit.megahit()
