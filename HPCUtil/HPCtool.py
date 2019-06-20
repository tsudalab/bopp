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
