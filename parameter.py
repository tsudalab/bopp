###############Initial parameters come from here#####################
#cutoffrate output from the second NN
cutoffrate=0.999 
#Number of for sampling within a node
numcore=8
#peptide length for sampling
peplength=10
#Generation loop within an epoch
geninter=1000
#Number of epoch
genepoch=1000
#working folder for sampling
wd="../test-len20"
#Using HPC or not. If set to True, HPC submission is activated.
HPC=True
#Choose the type of HPC submit script.
#1. OakForest (UTokyo)
#2. ISSP (UTokyo)
#3. Tsubame3 (TokyoTech)
#4. KComputer alike (IMS, NINS)
#5. IMS (NINS)
#6. ITO (UKyushu)
HPCtype=1
#Activate the user-defined function, if it set to False, MD-based evaluation is used
usereval=False
#Path to I-Tasser
HMpath="/work/gk73/k73003/software/I-TASSER5.1/I-TASSERmod/runI-TASSER.pl" 
#Path to I-Tasser library
HMlib="/work/gk73/k73003/software/I-TASSER5.1/libdir"
#Path to GROMACS executable file
GMXpath="/home/k0055/k005503/software/gromacs/bin/gmx_mpi"
#Calling GROMACS for running MPI jobs across the nodes in HPC"
GMXpathmpi="mpijob /home/k0055/k005503/software/gromacs/bin/gmx_mpi"
#Setting file for GROMACS
GMXconfig="source /home/k0055/k005503/software/gromacs/bin/GMXRC.bash"
#concentration of ions
conc=0.15 
#Setting files for AMBERTools
AMBERtoolconfig="source ~/software/amber16/amber.sh"
#needed for system run setting with OpenMP larger than 6 cores
ntomp=6 
#calling the mpi process
mpicall="mpijob -np 12 " 
#HPC group id, leave the groupid blank if you don't have
groupid="gk73" 
############End of parameters ##########################

############Don't change the below#####################
aalist=["B","A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","X"," "]
val=aalist
aalen=56 #56
dball=[] #this is for checking the database whether generated sequence is new or not 
genmod="model-1Feb2019-GRU256-64"
clasmod="AMPcls-GRU256-64"
