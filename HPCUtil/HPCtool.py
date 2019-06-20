from __future__ import print_function
import time,os

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

def writequeuescript(indname,wd,HPCtype):  #this should be recheck again
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
        f.close()
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
        f.close()
    elif HPCtype==3:
        print("Going to write Tsubame3 style queue script")
        f=open(wd+"/"+indname+"/"+indname+".jsub","w") 
        f.write("#!/bin/sh \n")
        f.write("#$-cwd \n")
        f.write("#$ -l f_node=1 \n")
        f.write("#$ -l h_rt=24:00:00 \n")
        f.write("#$ -N "+indname+" \n")
        f.write(". /etc/profile.d/modules.sh \n")
        f.write("module load cuda \n")
        f.write("module load intel \n")
        f.write("export OMP_NUM_THREADS=28 \n")
        f.close()
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
        #f.write("source /work/gk73/k73003/bashrc3 \n")
        f.close()
    elif HPCtype==5:
        print("Going to write IMS Molecular Simulator style queue script")
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
        f.close()
    else HPCtype==6:
        print("Going to write ITO supercomputer style queue script")
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
        f.close()
        #os.system("pjsub "+wd+"/"+indname+"/"+indname+".jsub")
    return 

def HPCsubmit(indname,wd,HPCtype,groupid):  #this should be recheck again 
    if HPCtype==1:
        print("Going to submit job to OakForest supercomputer")
        os.system("pjsub -g "+groupid+" "+wd+"/"+indname+"/"+indname+".jsub")
    elif HPCtype==2:
        print("Going to submit job to ISSP supercomputer")
        os.system("qsub "+wd+"/"+indname+"/"+indname+".jsub")
    elif HPCtype==3:
        print("Going to submit job to Tsubame3 supercomputer")
        os.system("qsub -g "+groupid+" "+wd+"/"+indname+"/"+indname+".jsub")
    elif HPCtype==4:
        print("Going to submit job to Kcomputer")
        os.system("pjsub -g "+groupid+" "+wd+"/"+indname+"/"+indname+".jsub")
    elif HPCtype==5:
        print("Going to submit job to IMS Molecular Simulator")
        os.system("jsub -q PN "+wd+"/"+indname+"/"+indname+".jsub")
    else HPCtype==6:
        print("Going to write ITO supercomputer style queue script")
        os.system("pjsub "+wd+"/"+indname+"/"+indname+".jsub")
    return 

def qstatcmd(HPCtype):  #this should be recheck again
    if HPCtype==1:
        print("Going to submit job to OakForest supercomputer")
        os.system("pjstat ")
    elif HPCtype==2:
        print("Going to submit job to ISSP supercomputer")
        os.system("qstat ")
    elif HPCtype==3:
        print("Going to submit job to Tsubame3 supercomputer")
        os.system("qstat ")
    elif HPCtype==4:
        print("Going to submit job to Kcomputer")
        os.system("pjstat ")
    elif HPCtype==5:
        print("Going to submit job to IMS Molecular Simulator")
        os.system("jobinfo -c -q PN -l ")
    else HPCtype==6:
        print("Going to write ITO supercomputer style queue script")
        os.system("pjstat ")
    return 
