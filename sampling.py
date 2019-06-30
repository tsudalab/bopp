from __future__ import print_function
from subprocess import Popen, PIPE
from math import *
import random,os
import numpy as np
import random as pr
from copy import deepcopy
#from types import IntType, ListType, TupleType, StringTypes
import itertools
import time
import math
import argparse
import subprocess
from keras.preprocessing import sequence
from keras.models import model_from_json
#from mpi4py import MPI
import sys,os 
from multiprocessing import Pool
from HMUtil.HMpara import * 
from MDutil.embedding import *
from HPCUtil.HPCtool import *
from GRURNN.execRNN import *
from GRURNN.execRNN import *
import time

from parameter import *

###############init the service#################
os.system(AMBERtoolconfig)
os.system(GMXconfig)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
import tensorflow as tf
tf.logging.set_verbosity(tf.logging.ERROR)
##################################################
#predefine amino acid list. B and space is for token and padding.
#aalist=["B","A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","X"," "]
#val=aalist
#aalen=56 #56



# evaluation step using homology modeling and MD simulation
def evaluate(seqsel,wd,seqfn,HPC,HPCtype,HMpath,HMlib,groupid,qstatcmd):
    #Call Homology modeling from here
    jobid=execHM(seqsel,wd,seqfn,HPC,HPCtype,HMpath,HMlib,groupid)
    if HPC==True:
        waitcheck(jobid,qstatcmd)
    #Call MD validating from here
    if HPC==True:
        execMD()
    return valpospep


# updating the neural network
def RNNupdate(model,seq):

    return


#main process come from here
def main():
    global aalen,aalist,val,wd,geninter,genepoch,peplength,numcore,cutoffrate
    #call generator and classifier
    pool = Pool(processes=numcore)
    pepgenerate=pool.map(actcrit,range(genepoch))
    tmppepgen=[]
    qstatcmd="qstat -a "
    #flattening the output list out of the ranks
    for x in range(0,len(pepgenerate)):
        for y in range(0,len(pepgenerate[x])): 
            print(pepgenerate[x][y])
            tmppepgen.append(pepgenerate[x][y])
    # removing the duplicates between ranks
    tmppepgen=list(set(tmppepgen))
    #save all to file
    fn=open(wd+"/genpep.txt","w")
    for x in range(0,len(tmppepgen)):
        fn.write(tmppepgen[x]+" \n")
    fn.close()
    print(tmppepgen)
    print(len(tmppepgen))

main()
