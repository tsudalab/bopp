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
from mpi4py import MPI
import sys,os 
from multiprocessing import Pool
from HMUtil.HMpara import * 
from MDutil.embedding import *
import time


#predefine amino acid list. B and space is for token and padding.
aalist=["B","A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","X"," "]
val=aalist
aalen=56 #56

#global generated peptide sequences
#pepgen=[]

# sequence generating subroutine
def actor(actormod,pepgennum,peplength,cntint):
    end1="X"
    end2="X"
    seq=[]
    intseq=[]
    print("Playing with sequence generator")
    random.seed(cntint)
    for i in range(0,pepgennum):
        #randomization to create input for generator NN
        get_int=[[int(np.random.rand()*20.0+1.0)] for t in range(peplength)] 
        x=np.reshape(get_int,(1,len(get_int)))
        x_pad= sequence.pad_sequences(x, maxlen=aalen, dtype='int32',padding='post', truncating='pre', value=22.0)
        predictions=actormod.predict(x_pad)
        #identify the generated sequence by generator NN
        tmpintseq=[np.argmax(predictions[0][cnt]) for cnt in range(0,len(predictions[0]))]
        seq.append(''.join([aalist[tmpintseq[cnt]] for cnt in range(0,peplength)]))
        intseq.append(tmpintseq[:peplength])
    return seq,intseq

# classification subroutine
def critic(criticmod,intseq,cutoffrate,ep,wd):
    cri_seq=[]
    #open file for save
    os.system("mkdir "+wd+"/ep"+str(ep))
    fnall=open(wd+"/ep"+str(ep)+"/genseqprob.txt","w")
    fncut=open(wd+"/ep"+str(ep)+"/cutseqprob.txt","w")
    #iterating through all the generated sequence
    for i in range(0,len(intseq)):
        x=np.reshape(intseq[i],(1,len(intseq[i])))
        #print(x)
        x_pad= sequence.pad_sequences(x, maxlen=aalen, dtype='int32',padding='post', truncating='pre', value=22.0)
        predictions=criticmod.predict(x_pad)
        # save all the generated sequence with probability
        fnall.write(''.join([aalist[intseq[i][cnt]] for cnt in range(0,len(intseq[i]))])+" "+str(predictions[0][0][0])+" \n")
        # check the cutoffrate and save only the new ones
        if predictions[0][0][0]>cutoffrate:
            tmpseq=''.join([aalist[intseq[i][cnt]] for cnt in range(0,len(intseq[i]))])
            print("Sequence "+tmpseq+" with probability at "+str(predictions[0][0][0]))
            cri_seq.append(tmpseq)
            if tmpseq in dbpall:
                print("Unfortunately, sequence is already found in databases!")
            else:
                fncut.write(''.join([aalist[intseq[i][cnt]] for cnt in range(0,len(intseq[i]))])+" "+str(predictions[0][0][0])+" \n")
    fnall.close()
    fncut.close()
    return cri_seq

# parallelizing generator and classifier
#def actcrit(pepgen,actormod,criticmod,geninter,peplength,cutoffrate):
def actcrit(cntint):
    global wd
    print("Dealing with epoch "+str(cntint))
    actormod=loadRNN("GRURNN","model-1Feb2019-GRU256-64")
    criticmod=loadRNN("GRURNN","AMPcls-GRU256-64") 
    gen_seq,intseq=actor(actormod,geninter,peplength,cntint)
    cri_seq=critic(criticmod,intseq,cutoffrate,cntint,wd)
    #pepgen.extend(cri_seq)
    #print(pepgen)
    return cri_seq 

# neural network model loading
def loadRNN(path,filename):
    json_file = open(path+"/"+filename+".json","r")
    RNNjson = json_file.read()
    json_file.close()
    loadRNN = model_from_json(RNNjson)
    loadRNN.load_weights(path+"/"+filename+".h5")
    return loadRNN

#check for waiting job
def waitcheck(superkonjobid,qstatcmd):
    leftjob=1
    while leftjob>0:
        os.system(qstatcmd+" > currjoblist.txt ")
        f=open("currjoblist.txt","r")
        cjl=f.readlines()
        f.close()
        leftjob=0
        for x in cjl:
            for y in superkonjobid:
                if y in x:
                    leftjob=leftjob+1
        if leftjob>0:
            time.sleep(300)
    return

# evaluation step using homology modeling and MD simulation
def evaluate(seqsel,wd,seqfn,HPC,HPCtype,HMpath,HMlib,qstatcmd):
    #Call Homology modeling from here
    jobid=execHM(seqsel,wd,seqfn,HPC,HPCtype,HMpath,HMlib)
    if HPC==True:
        waitcheck(jobid,qstatcmd)
    #Call MD validating from here
    
    return valpospep


# updating the neural network
def RNNupdate(model,seq):

    return


#main process come from here
if __name__ == "__main__":
    #pepgen=[]
    aalist=["B","A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","X"," "]
    cutoffrate=0.999
    numcore=8
    peplength=20
    geninter=1000
    genepoch=1000
    wd="test3-0_999"
    HPC=True
    HPCtype=1
    HMpath="/work/gk73/k73003/software/I-TASSER5.1/I-TASSERmod/runI-TASSER.pl"
    HMlib="/work/gk73/k73003/software/I-TASSER5.1/libdir"
    GMXpath="/home/k0055/k005503/software/gromacs/bin/gmx_mpi"
    GMXconfig="source /home/k0055/k005503/software/gromacs/bin/GMXRC.bash"
    AMBERtoolconfig="source ~/software/amber16/amber.sh"
    acpype="../acpype/scripts/acpype.py"
    os.system(AMBERtoolconfig)
    os.system(GMXconfig)
    ntomp=6 #needed for system run setting with openmp larger than 6 cores
    mpicall="mpijob -np 12 " #calling the mpi process
    
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


