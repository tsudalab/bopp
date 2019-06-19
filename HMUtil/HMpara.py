import os
#from mpi4py import MPI  

def writefasta(seq,indname,wd):
    f=open(wd+"/"+indname+"/seq.fasta","w") 
    f.write("> seq "+indname+" \n")
    f.write(seq+"\n")
    f.close()
    return 

def writesubmitoakforest(indname,wd):
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
    f.write("/work/gk73/k73003/software/I-TASSER5.1/I-TASSERmod/runI-TASSER.pl -libdir /work/gk73/k73003/software/I-TASSER5.1/libdir -seqname "+indname+" -datadir "+wd+"/indname"+" -runstyle gnuparallel \n")
    f.close()
    os.system("pjsub "+wd+"/"+indname+"/"+indname+".jsub")
    return 

if __name__ == "__main__":
    wd="/work/gk73/k73003/AMP-design/pepmcts-8Feb2019/actor-critic/test3-0_999"
    f=open(wd+"/genpep.txt","r")
    ln=f.readlines()
    for x in range(0,len(ln)):
        os.system("mkdir "+wd+"/"+"p3-"+str(x))
        writefasta(ln[x].strip(),"p3-"+str(x),wd)
        writesubmitoakforest("p3-"+str(x),wd)
