###############Initial parameters come from here#####################
cutoffrate=0.999
numcore=8
peplength=10
geninter=1000
genepoch=1000
wd="../test-len20"
HPC=True
HPCtype=1
usereval=False
cycle=10000
HMpath="/work/gk73/k73003/software/I-TASSER5.1/I-TASSERmod/runI-TASSER.pl"
HMlib="/work/gk73/k73003/software/I-TASSER5.1/libdir"
GMXpath="/home/k0055/k005503/software/gromacs/bin/gmx_mpi"
GMXpathmpi="mpijob /home/k0055/k005503/software/gromacs/bin/gmx_mpi"
GMXconfig="source /home/k0055/k005503/software/gromacs/bin/GMXRC.bash"
AMBERtoolconfig="source ~/software/amber16/amber.sh"
acpype="../acpype/scripts/acpype.py"
genmod="model-1Feb2019-GRU256-64"
clasmod="AMPcls-GRU256-64"
ntomp=6 #needed for system run setting with openmp larger than 6 cores
mpicall="mpijob -np 12 " #calling the mpi process
groupid="gk73" #leave the groupid blank if you don't have
############End of parameters ##########################

############Don't change the below#####################
aalist=["B","A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","X"," "]
val=aalist
aalen=56 #56
dball=[] #this is for checking the database whether generated sequence is new or not 
