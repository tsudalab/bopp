from __future__ import print_function
import numpy as np
import os
import sys
import multiprocessing as mp
from keras.models import Sequential
from keras.layers import Dense, Activation,TimeDistributed,MaxPooling1D
from keras.layers import LSTM,GRU
from keras.layers.embeddings import Embedding
from keras.optimizers import RMSprop, Adam
from keras.utils.data_utils import get_file
from keras.layers import Dropout
import numpy as np
import random
import sys
from keras.utils.np_utils import to_categorical
from keras.preprocessing import sequence
from keras.models import model_from_json
from random import sample
from parameter import *
#predefine amino acid list. B and space is for token and padding.
aalist=["B","A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","X"," "]


def updatemodel(seqinp,maxlnpep):
    #remove the tailing nextline character and adding token and padding seq
    #print(aalist)
    #tmp=seqinp.strip()
    tmp=seqinp.strip()+"X"
    while len(tmp)<=maxlnpep:
        tmp=tmp+" "
    #create the coding array for peptide sequence
    coding=[]
    seqid=[]
    for x in range(0,maxlnpep+1):
        #print(tmp[x])
        tmpctgr=to_categorical(aalist.index(tmp[x]), num_classes=len(aalist),dtype="int32")
        #print(tmpctgr)
        coding.append(tmpctgr) 
        seqid.append(aalist.index(tmp[x]))
    print("length of coding is "+str(len(coding)))
    #print(seqid)
    return seqid,coding

def loaddata(csvpath,csvpathneg,maxlnpep): 
    #load and read file
    f=open(csvpath,'r')
    ln=f.readlines()[1:]
    lenln=len(ln)
    clnpep=[]
    clncoding=[]
    f.close()  
    fn=open(csvpathneg,'r')
    lnn=fn.readlines()[1:]
    lenlnn=len(lnn)
    fn.close()
    #maxlnpep define the maximum length of the peptide
    #clnpep is the format sequence array
    #clean the line and add begin token and space padding to the data array
    datacutoff=0
    f=open("../AMP-data/RNN-dropoutdata-3Mar2019-GRU256-64.csv","w")
    #seqlist=sample(range(0,lenln),lenln-1000)
    seqlist=sample(range(0,lenln),lenln)  
    seqlistneg=sample(range(0,lenlnn),lenlnn)
    for i in range(0,lenln):
        print("process sequence "+str(i)+" over "+str(lenln))
        if (len(ln[i])<=maxlnpep)&(i in seqlist):
            frmseq,frmcod=seqfrmat(ln[i],maxlnpep)
            #frmcod=to_categorical(1, num_classes=2,dtype="int32")
            frmcod=[[1]]
            clnpep.append(frmseq)
            clncoding.append(frmcod)
        else:
            f.write(ln[i].strip()+"X"+"\n")
    for i in range(0,lenlnn):
        print("process negative sequence "+str(i)+" over "+str(lenlnn))
        if (len(ln[i])<=maxlnpep)&(i in seqlistneg):
            frmseq,frmcod=seqfrmat(lnn[i],maxlnpep)
            #frmcod=to_categorical(0, num_classes=2,dtype="int32")
            frmcod=[[0]]
            clnpep.append(frmseq)
            clncoding.append(frmcod)
        else:
            f.write(lnn[i].strip()+"X"+"\n")
    f.close()
    return clnpep,clncoding

def save_model(model):
    # serialize model to JSON
    model_json = model.to_json()
    with open("AMPcls-GRU256-64.json", "w") as json_file:
        json_file.write(model_json)
    # serialize weights to HDF5
    model.save_weights("AMPcls-GRU256-64.h5")
    print("Saved model to disk")
    return

if __name__ == "__main__":
    maxlnpep=55
    nproc=4
    X_data,Y_data=loaddata("../AMP-data/AMP-data-clean.csv","../AMP-data/nAMP_natrl.csv",maxlnpep)
    X=np.array((X_data))
    Y= np.array((Y_data))
    #initialize NN model
    model = Sequential()
    aalstln=len(aalist)
    print(aalstln)
    dataln=X.shape[1]
    print(dataln)
    print(X.shape)
    print(Y.shape)
    model.add(Embedding(input_dim=aalstln, output_dim=len(aalist), input_length=dataln,mask_zero=False))
    model.add(GRU(output_dim=256, activation='tanh',return_sequences=True))
    model.add(Dropout(0.2))
    model.add(Dense(1, activation='sigmoid')) 
    model.add(MaxPooling1D(pool_size=52))
    optimizer=Adam(lr=0.00001) # try much smaller one 0.001 0.00001
    print(model.summary())
    model.compile(loss='binary_crossentropy', optimizer=optimizer, metrics=['accuracy'])
    model.fit(X,Y,epochs=3000, batch_size=512,validation_split=0.1)
    save_model(model) 



