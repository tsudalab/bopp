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
#checkparam()
os.system(AMBERtoolconfig)
os.system(GMXconfig)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
import tensorflow as tf
tf.logging.set_verbosity(tf.logging.ERROR)
##################################################




# evaluation step using homology modeling and MD simulation
def evaluateMD(seqsel,wd,seqfn,HPC,HPCtype,HMpath,HMlib,groupid,qstatcmd):
    #Call Homology modeling from here
    jobid=execHM(seqsel,wd,seqfn,HPC,HPCtype,HMpath,HMlib,groupid)
    if HPC==True:
        waitcheck(jobid,qstatcmd)
    #Call MD validating from here
    jobid=execMD(len(seqsel),wd,acpype,GMXpath,ntomp,mpicall)
    if HPC==True:
        waitcheck(jobid,qstatcmd)
    
    return valpospep,valnegpep

#user-defined evaluation:
def evaluate(seqsel,wd,seqfn,HPC,HPCtype,HMpath,HMlib,groupid,qstatcmd):
    #user defined evaluation from here    
    ##########Explanation for the input of evaluation procedure#########
    #seqsel is the list of the selected of the sequence candidates by Neural network
    #wd is current working directory
    #seqfn is the file containing the sequence list
    #HPC is boolean variable to control the activation of using HPC 
    #HPCtype is the variable to check which HPC is using
    #HMlib is the variable indicating the library to I-TASSER
    #groupid is the groupid on HPC system
    #qstatcmd is the checking status of job on HPC
    return valpospep,valnegpep


#main process come from here
def main():
    global aalen,aalist,val,wd,geninter,genepoch,peplength,numcore,cutoffrate
    genmodroot=genmod
    clasmodroot=clasmod
    #sampling the actor-critic model
    #call generator and classifier
    pool = Pool(processes=numcore)
    pepgenerate=pool.map(actcrit,range(genepoch))
    tmppepgen=[]
    #flattening the output list out of the ranks
    for x in range(0,len(pepgenerate)):
        for y in range(0,len(pepgenerate[x])): 
            print(pepgenerate[x][y])
            mppepgen.append(pepgenerate[x][y])
    # removing the duplicates between ranks
    tmppepgen=list(set(tmppepgen))
    save all to file
    fn=open(wd+"/genpep.txt","w")
    for x in range(0,len(tmppepgen)):
        fn.write(tmppepgen[x]+" \n")
    fn.close()
    print(tmppepgen)
    print(len(tmppepgen))
    #call evaluation from here:
    if not(usereval):
        valpospep,valnegpep=evaluateMD(seqsel,wd,seqfn,HPC,HPCtype,HMpath,HMlib,groupid,qstatcmd)
    else:
        valpospep,valnegpep=evaluate(seqsel,wd,seqfn,HPC,HPCtype,HMpath,HMlib,groupid,qstatcmd)


main()
