import os
#from mpi4py import MPI  

def writefasta(seq,indname,wd):
    f=open(wd+"/"+indname+"/seq.fasta","w") 
    f.write("> seq "+indname+" \n")
    f.write(seq+"\n")
    f.close()
    return 

def writesubmitoakforest(indname,wd,HPCtype,HMpath,HMlib):
    if HPCtype==1:
        print("Going to write OakForest style queue script")
        f=open(wd+"/"+indname+"/"+indname+".jsub","w") 
        f.write("#------ pjsub option --------# \n")
        f.write("#PJM -L rscgrp=regular-flat \n")
        f.write("#PJM -L node=1 \n")
        f.write("#PJM --mpi proc=256 \n")
        f.write("#PJM --omp thread=1 \n")
        f.write("#PJM -L elapse=48:00:00 \n")
        f.write("#PJM -g gk73 \n")
        f.write("#------- Program execution -------# \n")
        f.write("source /work/gk73/k73003/bashrc3 \n")
        f.write(HMpath+" -libdir "+HMlib+" -seqname "+indname+" -datadir "+wd+"/indname"+" -runstyle gnuparallel \n")
        f.close()
        #os.system("pjsub "+wd+"/"+indname+"/"+indname+".jsub")
    elif HPCtype==2:
        print("Going to write ISSP supercomputer style queue script")
        f=open(wd+"/"+indname+"/"+indname+".jsub","w") 
        f.write("#!/bin/bash \n")
        f.write("#QSUB -queue F4cpu \n")
        f.write("#QSUB -node 1 \n")
        f.write("#QSUB -mpi 1 \n")
        f.write("#QSUB -omp 24 \n")
        f.write("#QSUB -place distribute \n")
        f.write("#QSUB -over false \n")
        f.write("#PBS -l walltime=24:00:00 \n")
        f.write("#PBS -N "+indname+" \n")
        f.write("cd $PBS_O_WORKDIR \n")
        f.write(HMpath+" -libdir "+HMlib+" -seqname "+indname+" -datadir "+wd+"/indname"+" -runstyle gnuparallel \n")
        f.close()
        #os.system("pjsub "+wd+"/"+indname+"/"+indname+".jsub")
    elif HPCtype==3:
        print("Going to write Tsubame3 style queue script")
        f=open(wd+"/"+indname+"/"+indname+".jsub","w") 
        f.write("#------ pjsub option --------# \n")
        f.write("#PJM -L rscgrp=regular-flat \n")
        f.write("#PJM -L node=1 \n")
        f.write("#PJM --mpi proc=256 \n")
        f.write("#PJM --omp thread=1 \n")
        f.write("#PJM -L elapse=48:00:00 \n")
        f.write("#PJM -g gk73 \n")
        f.write("#------- Program execution -------# \n")
        f.write("source /work/gk73/k73003/bashrc3 \n")
        f.write(HMpath+" -libdir "+HMlib+" -seqname "+indname+" -datadir "+wd+"/indname"+" -runstyle gnuparallel \n")
        f.close()
        #os.system("pjsub "+wd+"/"+indname+"/"+indname+".jsub")
    elif HPCtype==4:
        print("Going to write Kcomputer style queue script")
        f=open(wd+"/"+indname+"/"+indname+".jsub","w") 
        f.write("#------ pjsub option --------# \n")
        f.write("#PJM -L rscgrp=regular-flat \n")
        f.write("#PJM -L node=1 \n")
        f.write("#PJM --mpi proc=256 \n")
        f.write("#PJM --omp thread=1 \n")
        f.write("#PJM -L elapse=48:00:00 \n")
        f.write("#PJM -g gk73 \n")
        f.write("#------- Program execution -------# \n")
        f.write("source /work/gk73/k73003/bashrc3 \n")
        f.write(HMpath+" -libdir "+HMlib+" -seqname "+indname+" -datadir "+wd+"/indname"+" -runstyle gnuparallel \n")
        f.close()
        #os.system("pjsub "+wd+"/"+indname+"/"+indname+".jsub")
    elif HPCtype==5:
        print("Going to write Kcomputer style queue script")
        f=open(wd+"/"+indname+"/"+indname+".jsub","w") 
        f.write("#------ pjsub option --------# \n")
        f.write("#PJM -L rscgrp=regular-flat \n")
        f.write("#PJM -L node=1 \n")
        f.write("#PJM --mpi proc=256 \n")
        f.write("#PJM --omp thread=1 \n")
        f.write("#PJM -L elapse=48:00:00 \n")
        f.write("#PJM -g gk73 \n")
        f.write("#------- Program execution -------# \n")
        f.write("source /work/gk73/k73003/bashrc3 \n")
        f.write(HMpath+" -libdir "+HMlib+" -seqname "+indname+" -datadir "+wd+"/indname"+" -runstyle gnuparallel \n")
        f.close()
        #os.system("pjsub "+wd+"/"+indname+"/"+indname+".jsub")
    else HPCtype==6:
        print("Going to write Kcomputer style queue script")
        f=open(wd+"/"+indname+"/"+indname+".jsub","w") 
        f.write("#------ pjsub option --------# \n")
        f.write("#PJM -L rscgrp=regular-flat \n")
        f.write("#PJM -L node=1 \n")
        f.write("#PJM --mpi proc=256 \n")
        f.write("#PJM --omp thread=1 \n")
        f.write("#PJM -L elapse=48:00:00 \n")
        f.write("#PJM -g gk73 \n")
        f.write("#------- Program execution -------# \n")
        f.write("source /work/gk73/k73003/bashrc3 \n")
        f.write(HMpath+" -libdir "+HMlib+" -seqname "+indname+" -datadir "+wd+"/indname"+" -runstyle gnuparallel \n")
        f.close()
        #os.system("pjsub "+wd+"/"+indname+"/"+indname+".jsub")
    return 

#if __name__ == "__main__":
def execHM(seqsel,wd,seqfn,HPC,HPCtype,HMpath,HMlib):  
    #wd="/work/gk73/k73003/AMP-design/pepmcts-8Feb2019/actor-critic/test3-0_999"
    f=open(wd+"/"+seqfn,"r")
    ln=f.readlines()
    for x in range(0,len(ln)):
        os.system("mkdir "+wd+"/"+"p3-"+str(x))
        writefasta(ln[x].strip(),"p3-"+str(x),wd)
        writesubmitoakforest("p3-"+str(x),wd,HPCtype)
        if HPC==True:
            currwd=os.getcwd()
            os.system("cd "+wd+"/"+indname+"/")
            os.system("pjsub "+wd+"/"+indname+"/"+indname+".jsub")
            os.system("cd "+currwd)
            print("The job "+idname+" has been dispatched to HPC.")
        else:
            print("No HPC is selected. Proceed the HM on this system.")
            currwd=os.getcwd()
            os.system("cd "+wd+"/"+indname+"/")
            os.system(HMpath+" -libdir "+HMlib+" -seqname "+indname+" -datadir "+wd+"/indname"+" -runstyle gnuparallel")
            os.system("cd "+currwd)
    return superkonjobid


    