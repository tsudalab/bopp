import os
from HPCUtil.HPCtool import *
from parameter import *
#from mpi4py import MPI  

def writefasta(seq,indname,wd): #write selected sequence in FASTA format
    f=open(wd+"/"+indname+"/seq.fasta","w") 
    f.write("> seq "+indname+" \n")
    f.write(seq+"\n")
    f.close()
    return 

def writeHM2HPC(indname,wd,HPCtype,HMpath,HMlib): #append the I-TASSER command to submit script 
    f=open(wd+"/"+indname+"/"+indname+".jsub","a") 
    f.write(HMpath+" -libdir "+HMlib+" -seqname "+indname+" -datadir "+wd+"/indname"+" -runstyle gnuparallel \n")
    f.close()
    return 

#if __name__ == "__main__":
def execHMHPC(seqsel,wd,seqfn,HPC,HPCtype,HMpath,HMlib,groupid):  #executing the Homology Modeling stage 
    f=open(wd+"/"+seqfn,"r")
    ln=f.readlines()
    superkonjobid=[]
    for x in range(0,len(ln)):
        os.system("mkdir "+wd+"/"+"p3-"+str(x))
        writefasta(ln[x].strip(),"p3-"+str(x),wd)
        if HPC==True:
            currwd=os.getcwd()
            os.system("cd "+wd+"/"+indname+"/")
            writeHM2HPC("p3-"+str(x),wd,HPCtype)
            HPCsubmit(indname,wd,HPCtype,groupid)
            os.system("cd "+currwd)
            print("The job "+idname+" has been dispatched to HPC.")
            superkonjobid.append(indname)
        else:
            print("No HPC is selected. Proceed the HM on this machine.")
            currwd=os.getcwd()
            os.system("cd "+wd+"/"+indname+"/")
            os.system(HMpath+" -libdir "+HMlib+" -seqname "+indname+" -datadir "+wd+"/indname"+" -runstyle gnuparallel")
            os.system("cd "+currwd)
    return superkonjobid

def execHMserial(seqsel,wd,seqfn,HPC,HPCtype,HMpath,HMlib,groupid):  #executing the Homology Modeling stage 
    f=open(wd+"/"+seqfn,"r")
    ln=f.readlines()
    superkonjobid=[]
    for x in range(0,len(ln)):
        os.system("mkdir "+wd+"/"+"p3-"+str(x))
        writefasta(ln[x].strip(),"p3-"+str(x),wd)
        if HPC==True:
            currwd=os.getcwd()
            os.system("cd "+wd+"/"+indname+"/")
            writeHM2HPC("p3-"+str(x),wd,HPCtype)
            HPCsubmit(indname,wd,HPCtype,groupid)
            os.system("cd "+currwd)
            print("The job "+idname+" has been dispatched to HPC.")
            superkonjobid.append(indname)
        else:
            print("No HPC is selected. Proceed the HM on this machine.")
            currwd=os.getcwd()
            os.system("cd "+wd+"/"+indname+"/")
            os.system(HMpath+" -libdir "+HMlib+" -seqname "+indname+" -datadir "+wd+"/indname"+" -runstyle gnuparallel")
            os.system("cd "+currwd)
    return superkonjobid

    