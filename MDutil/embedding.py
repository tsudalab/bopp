from __future__ import print_function
import os
import mdtraj as md
import numpy as np
from math import *
from parameter import *
# checking protein size and protein center of mass
def protsize(path,name,wd):
    if not os.path.isfile(path+"/"+name):
        print("File "+path+"/"+name+" does not exist! Please check.")
        quit
    else:
        prot=md.load(path+"/"+name)
        minx=np.min(prot.xyz[0][:,0])
        miny=np.min(prot.xyz[0][:,1])
        minz=np.min(prot.xyz[0][:,2])
        maxx=np.max(prot.xyz[0][:,0])
        maxy=np.max(prot.xyz[0][:,1])
        maxz=np.max(prot.xyz[0][:,2])
        protcom=md.compute_center_of_mass(prot)
        #process pdb with AMBERTool
        #os.system("pdb4amber -i "+path+"/"+name+" -o "+wd+"/prot-reduce.pdb --reduce --dry")
        os.system("reduce -Trim "+path+"/"+name+" > "+wd+"/prot-reduce.pdb ")
    return minx,maxx,miny,maxy,minz,maxz,protcom[0]

# translating protein to fit to simulation box
def translate(prot,dx,dy,dz): 
    for x in range(0,len(prot.xyz[0])):
        prot.xyz[0][x,0]=prot.xyz[0][x,0]+dx 
        prot.xyz[0][x,1]=prot.xyz[0][x,1]+dy  
        prot.xyz[0][x,2]=prot.xyz[0][x,2]+dz
    return prot

######################remove overlapping with protein (this should be deprecated) ##################
def overlapremove(pdb,prot,wd):  
    # create index list 
    protid=prot.top.select("all")
    watid=pdb.top.select("(resname WAT) or (resname KI) or (resname CLI)")
    table, bonds = pdb.top.to_dataframe()
    print(table.head())
    overlap=[]
    # compute the distance with water
    dislist=[]
    for x in range(0,len(watid)):
        for y in range(0,len(protid)):
            dislist.append([watid[x],protid[y]])
    dist=md.compute_distances(pdb,dislist)  
    print(dist)
    for x in range(0,len(watid)):
        print(dislist[:][0])
        if dist[np.where(dislist[0][:]==watid[x])[0]].min()<3.5:
            overlap.extend(pdb.top.select("(chainid 1) and (index "+str(watid[x])+")"))
    print(overlap)
    overlap=list(set(overlap))
    nonoverlap=pdb.top.select("all")
    print(nonoverlap)
    for x in range(0,len(overlap)):
        index = np.argwhere(nonoverlap==overlap[x])
        nonoverlap=np.delete(nonoverlap,index)
    print(nonoverlap)
    pdbcut=pdb.atom_slice(nonoverlap)
    pdbcut.save_pdb(wd+"/lipprot_nonoverlap.pdb")
    return overlap
##################################(this should be deprecated) ####################################

# remove overlapping with protein
def addter(pdb,wd,conc,jobname):  
    solvateions=True   #add waters and ions  
    # simulation box dimensions
    xdim=92.5888
    ydim=89.6671
    zdim=118.2444
    # adding TER to the pdb file 
    f=open(wd+"/"+pdb,"r")
    fo=open(wd+"/"+pdb[:-4]+"-addter.pdb","w")
    ln=f.readlines()
    ind=0
    resind=0
    oldres=0
    #sorting the resid and index, and translating the coordinates.
    for x in range(0,len(ln)):  
        tmp=ln[x]
        if (tmp.find("ATOM") is not -1 ):
            ind=ind+1
            print(tmp[23:27])
            if (int(tmp[23:27])!=oldres):
                oldres=int(tmp[23:27])
                resind=resind+1
            if ((("H18T" == tmp[12:16]) and ("  OL" == tmp[16:20]) ) or ((" HO1" == tmp[12:16]) and (" CHL" == tmp[16:20]) )  ) :
                print(float(tmp[30:38]),float(tmp[38:46]),float(tmp[46:54]))
                fo.write(tmp[0:4]+"{:7d}".format(ind)+tmp[11:22]+"{:4d}".format(resind)+tmp[26:30]+"{:8.3f}{:8.3f}{:8.3f}".format(float(tmp[30:38])-xdim/2.0,float(tmp[38:46])-ydim/2.0,float(tmp[46:54])-zdim/2.0)+tmp[54:])
                fo.write("TER \n")
            else:
                print(float(tmp[30:38]),float(tmp[38:46]),float(tmp[46:54]))
                fo.write(tmp[0:4]+"{:7d}".format(ind)+tmp[11:22]+"{:4d}".format(resind)+tmp[26:30]+"{:8.3f}{:8.3f}{:8.3f}".format(float(tmp[30:38])-xdim/2.0,float(tmp[38:46])-ydim/2.0,float(tmp[46:54])-zdim/2.0)+tmp[54:])     
        else:
            fo.write(tmp) 
    f.close()
    fo.close()
    # create tleap file to generate the amber run file
    if solvateions:
        f=open(wd+"/"+pdb[:-4]+"tleap.script","w")
        f.write("logFile "+wd+"/"+"leap.log  \n")
        f.write("source leaprc.protein.ff14SB \n")
        f.write("source leaprc.lipid14 \n")
        f.write("source leaprc.water.tip3p \n")
        f.write("loadamberparams frcmod.ionsjc_tip3p \n")
        f.write("logFile "+wd+"/"+"leap.log  \n")
        f.write("complex = loadpdb "+wd+"/"+pdb[:-4]+"-addter.pdb"+" \n")
        f.write("set complex box { "+str(xdim)+" "+str(ydim)+" "+str(zdim)+" } \n")
        f.write("solvateBox complex TIP3PBOX 0.0 1.0 \n")  
        f.write("charge complex \n")
        #f.write("logFile "+wd+"/"+"leap.log  \n")
        f.write("quit \n")
        f.close()
        os.system("tleap -f "+wd+"/"+pdb[:-4]+"tleap.script ")
        f=open(wd+"/"+"leap.log","r")
        ln=f.readlines()
        print(ln)  
        f.close()
        for x in range(0,len(ln)):
            if (ln[x].find("Total perturbed charge:") is not -1):
                charge=float(ln[x].split()[3])   
                print(charge)
    f=open(wd+"/"+pdb[:-4]+"tleap.script","w")
    f.write("source leaprc.protein.ff14SB \n")
    f.write("source leaprc.lipid14 \n")
    f.write("source leaprc.water.tip3p \n")
    f.write("loadamberparams frcmod.ionsjc_tip3p \n")
    f.write("complex = loadpdb "+wd+"/"+pdb[:-4]+"-addter.pdb"+" \n")
    f.write("set complex box { "+str(xdim)+" "+str(ydim)+" "+str(zdim)+" } \n")
    if solvateions:
        f.write("solvateBox complex TIP3PBOX {0.0 0.0 20.0}  1.0 \n") 
        #f.write("solvateBox complex TIP3PBOX 20.0 \n") 
        vol=xdim*ydim*zdim
        ions=int(conc*6.023*0.0001*vol)
        print(str(ions)+" ions need to keep "+str(conc)+"M concentration.")    
        if ions<=0:
            f.write("addIonsRand complex NA "+str(ions+round(abs(charge)))+" CL "+str(ions)+" \n") 
        else:
            f.write("addIonsRand complex NA "+str(ions)+" CL "+str(ions+round(abs(charge)))+" \n")        
    f.write("saveamberparm complex "+wd+"/"+pdb[:-4]+"amber.prmtop "+wd+"/"+pdb[:-4]+"amber.inpcrd \n")  
    f.write("savepdb complex "+wd+"/"+pdb[:-4]+"amber.pdb  \n")
    f.write("quit")
    f.close()
    # call tleap to generate the file
    os.system("tleap -f "+wd+"/"+pdb[:-4]+"tleap.script")
    # call acpype to convert to GMX file
    os.system(acpype+" -r -p "+wd+"/"+pdb[:-4]+"amber.prmtop -x "+wd+"/"+pdb[:-4]+"amber.inpcrd -b  "+wd+"/"+pdb[:-4]+"amber ")   
    # create position restraint file 
    makendx(wd,pdb[:-4]+"amber.pdb",wd,pdb[:-4]+"amber.pdb",pdb[:-4])
    addrestrtopol(wd,pdb[:-4]+"amber_GMX.top",wd,pdb[:-4]+"amber_GMX.gro",wd,pdb[:-4]+"amber.pdb",pdb[:-4],jobname)
    # create the run file for GROMACS
    
    return

# insert protein to lipid simulation box
def insertion(path,name,lipidpath,lipidname,wd,conc,jobname): 
    if not os.path.isfile(lipidpath+"/"+lipidname):
        print("File "+lipidpath+"/"+lipidname+" does not exist! Please check.")
        quit
    else:
        minx,maxx,miny,maxy,minz,maxz,protcom=protsize(path,name,wd) 
        lipid=md.load(lipidpath+"/"+lipidname) 
        lowzliphead= lipid.xyz[0][lipid.topology.select("(resname CHL) or (resname OL) or (resname PA) or (resname PC)")][:,2].min()
        boxminx=np.min(lipid.xyz[0][:,0])
        boxminy=np.min(lipid.xyz[0][:,1])
        boxminz=np.min(lipid.xyz[0][:,2])
        boxmaxx=np.max(lipid.xyz[0][:,0])
        boxmaxy=np.max(lipid.xyz[0][:,1])
        boxmaxz=np.max(lipid.xyz[0][:,2])
        lipxdim=boxmaxx-boxminx
        lipydim=boxmaxy-boxminy
        lipzdim=boxmaxz-boxminz
        protxdim=maxx-minx
        protydim=maxy-miny
        protzdim=maxz-minz
        if ((protxdim>(lipxdim-2.0*1.0)) or (protydim>(lipydim-2.0*1.0)) or (protzdim>(lipzdim-2.0*1.0))):
            print("It seems protein size is larger than the good lipid box in term of simulation. Please replace!")
            quit
        else:
            #align the protein to the center of water region
            prot=md.load(wd+"/"+"prot-reduce.pdb")  
            prot=translate(prot,(lipxdim-protxdim)/2.0-maxx,(lipydim-protydim)/2.0-maxy,(lipzdim-protzdim)/2.0-maxz)
            prot.save_pdb(wd+"/"+"prot-trans.pdb")
            nonoverlap=[]
            rescheck=True
            protlip=prot.stack(lipid.atom_slice(lipid.top.select("(resname CHL) or (resname OL) or (resname PA) or (resname PC)")))
            #protlip=prot.stack(lipid)
            protlip.save_pdb(wd+"/"+name[:-4]+"protlip.pdb")
            #protlip=md.load("protlip.pdb")
            #overlapid=addter(wd+"/"+name[:-4]+"protlip.pdb",wd,conc,jobname)
            overlapid=addter(name[:-4]+"protlip.pdb",wd,conc,jobname)
    return

#make index file for GROMACS
#def makendx(fn,path,topfn,toppath):
def makendx(path,fn,toppath,topfn,prefix):
    #since gromacs index starts from 1, this should be +1  
    print(path+"/"+fn)
    print(toppath+"/"+topfn)
    trj=md.load(path+"/"+fn,top=toppath+"/"+topfn) 
    #create lipid index group 
    f=open(path+"/"+prefix+"index.ndx","w")
    indlist=["all","protein","protein and (not (name =~ 'H.*'))","name =~ 'CA'","backbone","sidechain","water or resname NA or resname CL","water","resname NA","resname CL","resname CHL","resname OL","resname PA","resname PC","(resname CHL) or (resname OL) or (resname PA) or (resname PC) ","(backbone or (resname CHL) or (resname OL) or (resname PA) or (resname PC)) and (not (name =~ 'H.*'))"]
    indname=["System","Protein","Protein-H","C-alpha","Backbone","Sidechain","Sol","Waters","NA","CL","CHL","OL","PA","PC","Lipid","ProtLip"]
    for y in range(0,len(indlist)):
        sel=trj.top.select(indlist[y])
        print(indlist[y])
        print(sel)  
        f.write("[ "+indname[y]+" ] \n")
        a=""
        for x in range(0,len(sel)):
            a=a+str(sel[x]+1)+" "
            if ((x%15==14) or (x==len(sel)-1)):
                f.write(a+" \n")
                a=""
    f.close()
    return

# write position restraint file
def posrewrite(fn,path,indatom,fcx,fcy,fcz):  
    f=open(path+"/"+fn,'w')
    f.write("[ position_restraints ] \n")
    f.write("; atom  type      fx      fy      fz \n")
    for x in indatom:
        f.write("{:10d}{:10d}{:10.3f}{:10.3f}{:10.3f}  \n".format(x+1,1,fcx,fcy,fcz) )    
    return

#write the energy minimization mdp file
def emwrite(path,fn):
    f=open(path+"/"+fn,'w')
    f.write("title           = Minimization \n")
    f.write("integrator	     = steep \n")
    f.write("emtol		     = 50.0  \n")
    f.write("emstep          = 0.001 \n")
    f.write("nsteps		     = 50000000 \n")
    f.write("energygrps	     = Protein Lipid Sol \n")
    f.write("nstlist	     = 1 \n")
    f.write("cutoff-scheme   = Verlet \n")
    f.write("ns_type		 = grid \n")
    f.write("rlist		     = 1.0 \n")
    f.write("coulombtype	 = PME \n")
    f.write("rcoulomb	     = 1.0 \n")
    f.write("rvdw		     = 1.0 \n")
    f.write("pbc		     = xyz \n")    
    return

#write the nvt mdp file
def nvtposrewrite(path,fn,temp):
    f=open(path+"/"+fn,'w')
    f.write("title		     = NVT simulation \n")
    f.write("define		     = -DPROT -DLIP \n")
    f.write("integrator	     = md \n")
    f.write("nsteps		     = 500000 \n")
    f.write("dt		         = 0.002 \n")
    f.write("nstxout    		 = 50000 \n")
    f.write("nstvout		     = 50000 \n")
    f.write("nstenergy	         = 5000 \n")
    f.write("nstlog		         = 5000 \n")
    f.write("nstxout-compressed  = 5000 \n")
    f.write("compressed-x-grps   = System \n")
    f.write("continuation	         = no\n")
    f.write("constraint_algorithm    = lincs\n")
    f.write("constraints	         = h-bonds \n")
    f.write("lincs_iter	             = 1	 \n") 
    f.write("lincs_order	         = 4	 \n")
    f.write("cutoff-scheme       = Verlet   \n")
    f.write("ns_type		     = grid	 \n") 
    f.write("nstlist		     = 10   \n")
    f.write("rcoulomb	         = 1.0	 \n") 
    f.write("rvdw		         = 1.0   \n")
    f.write("coulombtype	     = PME	 \n") 
    f.write("pme_order	         = 4   \n")
    f.write("fourierspacing	        = 0.16	 \n") 
    f.write("tcoupl		            = v-rescale   \n")
    f.write("tc-grps	            = Protein Lipid Sol		 \n") 
    f.write("tau_t		            = 0.4 0.4 0.4  \n")
    f.write("ref_t	                = {:8.3f} {:8.3f} {:8.3f}	 \n".format(temp,temp,temp)) 
    f.write("pbc	             = xyz	 \n") 
    f.write("DispCorr	         = EnerPres   \n")
    f.write("gen_vel		        = yes	 \n") 
    f.write("gen-temp               = {:8.3f}   \n".format(temp))
    f.write("gen-seed               = -1	 \n") 
    f.close()   
    return

#write the npt mdp file
def nptposrewrite(path,fn,temp,press):
    f=open(path+"/"+fn,"w")
    f.write("title		     = NPT simulation \n")
    f.write("define		     = -DPROT  \n")
    f.write("integrator	     = md \n")
    f.write("nsteps		     = 10000000 \n")
    f.write("dt		         = 0.002 \n")
    f.write("nstxout    		 = 50000 \n")
    f.write("nstvout		     = 50000 \n")
    f.write("nstenergy	         = 5000 \n")
    f.write("nstlog		         = 5000 \n")
    f.write("nstxout-compressed  = 5000 \n")
    f.write("compressed-x-grps   = System \n")
    f.write("continuation	         = yes \n")
    f.write("constraint_algorithm    = lincs \n")
    f.write("constraints	         = h-bonds \n")
    f.write("lincs_iter	             = 1	 \n") 
    f.write("lincs_order	         = 4	 \n")
    f.write("cutoff-scheme       = Verlet   \n")
    f.write("ns_type		     = grid	 \n") 
    f.write("nstlist		     = 10   \n")
    f.write("rcoulomb	         = 1.0	 \n") 
    f.write("rvdw		         = 1.0   \n")
    f.write("coulombtype	     = PME	 \n") 
    f.write("pme_order	         = 4   \n")
    f.write("fourierspacing	        = 0.16	 \n") 
    f.write("tcoupl		            = v-rescale   \n")
    f.write("tc-grps	            = Protein Lipid Sol	 \n") 
    f.write("tau_t		            = 0.4 0.4 0.4  \n")
    f.write("ref_t	                = {:8.3f} {:8.3f} {:8.3f}	 \n".format(temp,temp,temp)) 
    f.write("pcoupl		         = berendsen	 \n") 
    #f.write("pcoupltype	         = semiisotropic  \n")
    #f.write("tau_p		         = 2.0	 \n") 
    #f.write("compressibility		         = 4.5e-5  4.5e-5   \n")
    #f.write("ref_p    = {:8.3f} {:8.3f}  \n".format(press,press) )
    f.write("pcoupltype                 = isotropic  \n")
    f.write("tau_p                      = 2.0   \n")
    f.write("compressibility                    =  4.5e-5   \n")
    f.write("ref_p    = {:8.3f}  \n".format(press) )
    f.write("refcoord-scaling	 = com	 \n") 
    f.write("pbc	                = xyz	 \n") 
    f.write("DispCorr	            = EnerPres   \n")
    f.write("gen_vel		        = no	 \n") 
    f.close()   
    return

#write the npt mdp file
def smdposrewrite(path,fn,temp,press):
    f=open(path+"/"+fn,'w')
    f.write("title		     = SMD simulation \n")
    f.write("integrator	     = md \n")
    f.write("nsteps		     = 15000000 \n")
    f.write("dt		         = 0.002 \n")
    f.write("nstxout    		 = 50000 \n")
    f.write("nstvout		     = 50000 \n")
    f.write("nstenergy	         = 5000 \n")
    f.write("nstlog		         = 5000 \n")
    f.write("nstxout-compressed  = 5000 \n")
    f.write("compressed-x-grps   = System \n")
    f.write("continuation	         = yes \n")
    f.write("constraint_algorithm    = lincs \n")
    f.write("constraints	         = h-bonds \n")
    f.write("lincs_iter	             = 1	 \n") 
    f.write("lincs_order	         = 4	 \n")
    f.write("cutoff-scheme       = Verlet   \n")
    f.write("ns_type		     = grid	 \n") 
    f.write("nstlist		     = 10   \n")
    f.write("rcoulomb	         = 1.0	 \n") 
    f.write("rvdw		         = 1.0   \n")
    f.write("coulombtype	     = PME	 \n") 
    f.write("pme_order	         = 4   \n")
    f.write("fourierspacing	        = 0.16	 \n") 
    f.write("tcoupl		            = v-rescale   \n")
    f.write("tc-grps	            = Protein Lipid Sol	 \n") 
    f.write("tau_t		            = 0.4 0.4 0.4  \n")
    f.write("ref_t	                = {:8.3f} {:8.3f} {:8.3f}	 \n".format(temp,temp,temp)) 
    #f.write("pcoupl		         = berendsen	 \n") 
    #f.write("pcoupltype	         = semiisotropic  \n")
    #f.write("tau_p		         = 2.0	 \n") 
    #f.write("compressibility		         = 4.5e-5  4.5e-5    \n")
    #f.write("ref_p     = {:8.3f} {:8.3f}  \n".format(press,press) )
    #f.write("refcoord-scaling	 = com	 \n") 
    f.write("pbc	                = xyz	 \n") 
    f.write("DispCorr	            = EnerPres   \n")
    f.write("gen_vel		        = no	 \n") 
    f.write("pull                 = yes	 \n") 
    f.write("pull-ngroups         = 1	 \n") 
    f.write("pull-group1-name     = Protein	 \n") 
    f.write("pull-coord1-type     = umbrella	 \n") 
    f.write("pull-coord1-geometry = direction-periodic	 \n") 
    f.write("pull-coord1-vec      = 0 0 1	 \n") 
    f.write("pull-coord1-k        = 1000	 \n") 
    f.write("pull-coord1-rate     = 0.001	 \n") 
    f.write("pull-coord1-groups   = 0 1	 \n") 
    f.close()   
    return

def writerun(path,fn,prefix,jobname):
    f=open(path+"/"+fn,"w") 
    #f.write("#!/bin/bash \n")  
    #f.write("#QSUB -queue B18acc \n")
    #f.write("#QSUB -node 8 \n")
    #f.write("#QSUB -mpi 32 \n")
    #f.write("#QSUB -omp 6 \n")
    #f.write("#QSUB -place distribute \n")
    #f.write("#QSUB -over false \n")
    #f.write("#PBS -l walltime=5:59:00 \n")
    #f.write("#QSUB -place distribute \n")
    #f.write("#PBS -N "+jobname+" \n")
    #f.write("cd $PBS_O_WORKDIR  \n")  
    #f.write(". /etc/profile.d/modules.sh  \n")  
    #f.write("module unload intel/16.0.1.150 mpt/2.12 gnu/4.8.5 cuda/7.0  \n")  
    #f.write("export GMX_MAXBACKUP=-1  \n")  
    #f.write("module load intel/15.0.0.090 intel-mpi/5.0.3.048 intel-mkl/15.0.0.090 gnu/4.8.5 cuda/7.0  \n")
    f.write(GMXpath+" grompp -f "+path+"/"+prefix+"em.mdp -c "+path+"/"+prefix+"amber_GMX.gro -p "+path+"/"+prefix+"amber_GMX.top -n "+path+"/"+prefix+"index.ndx -o "+path+"/"+prefix+"em.tpr -maxwarn 10   \n")
    f.write(GMXpathmpi+" mdrun -deffnm "+path+"/"+prefix+"em -v  \n")
    f.write(GMXpath+" grompp -f "+path+"/"+prefix+"nvt.mdp -c "+path+"/"+prefix+"em.gro -p "+path+"/"+prefix+"amber_GMX.top -n "+path+"/"+prefix+"index.ndx -o "+path+"/"+prefix+"nvt.tpr -maxwarn 10   \n")
    f.write(GMXpathmpi+" mdrun -deffnm "+path+"/"+prefix+"nvt -v  \n")
    f.write(GMXpath+" grompp -f "+path+"/"+prefix+"npt.mdp -c "+path+"/"+prefix+"nvt.gro -p "+path+"/"+prefix+"amber_GMX.top -n "+path+"/"+prefix+"index.ndx -o "+path+"/"+prefix+"npt.tpr -maxwarn 1   \n")
    f.write(GMXpathmpi+" mdrun -deffnm "+path+"/"+prefix+"npt -v  \n")
    f.write(GMXpath+" grompp -f "+path+"/"+prefix+"smd.mdp -c "+path+"/"+prefix+"npt.gro -p "+path+"/"+prefix+"amber_GMX.top -n "+path+"/"+prefix+"index.ndx -o "+path+"/"+prefix+"smd.tpr -maxwarn 1   \n")
    f.write(GMXpathmpi+" mdrun -deffnm "+path+"/"+prefix+"smd -v  \n") 
    f.close()  
    return

def writeissprunrestrt(path,fn,prefix,jobname):
    f=open(path+"/"+fn,"w")
    f.write("#!/bin/bash \n")
    f.write("#QSUB -queue F18acc \n")
    f.write("#QSUB -node 8 \n")
    f.write("#QSUB -mpi 32 \n")
    f.write("#QSUB -omp 6 \n")
    f.write("#QSUB -place distribute \n")
    f.write("#QSUB -over false \n")
    f.write("#PBS -l walltime=23:59:00 \n")
    f.write("#QSUB -place distribute \n")
    f.write("#PBS -N "+jobname+" \n")
    f.write("cd $PBS_O_WORKDIR  \n")
    f.write(". /etc/profile.d/modules.sh  \n")
    f.write("module unload intel/16.0.1.150 mpt/2.12 gnu/4.8.5 cuda/7.0  \n")
    f.write("export GMX_MAXBACKUP=-1  \n")
    f.write("module load intel/15.0.0.090 intel-mpi/5.0.3.048 intel-mkl/15.0.0.090 gnu/4.8.5 cuda/7.0  \n")
    f.write("~/software/gromacs/bin/gmx_gpu5.1.2 grompp -f "+path+"/"+prefix+"em.mdp -c "+path+"/"+prefix+"amber_GMX.gro -p "+path+"/"+prefix+"amber_GMX.top -n "+path+"/"+prefix+"index.ndx -o "+path+"/"+prefix+"em.tpr -maxwarn 10   \n")
    f.write("mpijob ~/software/gromacs/bin/gmx_gpu5.1.2 mdrun -deffnm "+path+"/"+prefix+"em -v  \n")
    f.write("~/software/gromacs/bin/gmx_gpu5.1.2 grompp -f "+path+"/"+prefix+"nvt.mdp -c "+path+"/"+prefix+"em.gro -p "+path+"/"+prefix+"amber_GMX.top -n "+path+"/"+prefix+"index.ndx -o "+path+"/"+prefix+"nvt.tpr -maxwarn 10   \n")
    f.write("mpijob ~/software/gromacs/bin/gmx_gpu5.1.2 mdrun -deffnm "+path+"/"+prefix+"nvt -v  \n")
    f.write("~/software/gromacs/bin/gmx_gpu5.1.2 grompp -f "+path+"/"+prefix+"npt.mdp -c "+path+"/"+prefix+"nvt.gro -p "+path+"/"+prefix+"amber_GMX.top -n "+path+"/"+prefix+"index.ndx -o "+path+"/"+prefix+"npt.tpr -maxwarn 1   \n")
    f.write("mpijob ~/software/gromacs/bin/gmx_gpu5.1.2 mdrun -deffnm "+path+"/"+prefix+"npt -v  \n")
    f.write("~/software/gromacs/bin/gmx_gpu5.1.2 grompp -f "+path+"/"+prefix+"smd.mdp -c "+path+"/"+prefix+"npt.gro -p "+path+"/"+prefix+"amber_GMX.top -n "+path+"/"+prefix+"index.ndx -o "+path+"/"+prefix+"smd.tpr -maxwarn 10   \n")
    f.write("mpijob ~/software/gromacs/bin/gmx_gpu5.1.2 mdrun -deffnm "+path+"/"+prefix+"smd -v  \n")
    f.close() 
    return


def addrestrtopol(path,fn,trjpath,trjfn,toppath,topfn,prefix,jobname):
    #adding position restraint to the topology file
    trj=md.load(trjpath+"/"+trjfn,top=toppath+"/"+topfn)
    f=open(path+"/"+fn,"r")
    ln=f.readlines()
    f.close()
    f=open(path+"/"+fn,"w") 
    lenln=len(ln)
    for x in range(0,len(ln)):
        if ((x)<lenln-1):
            if ((ln[x+1].find("[ moleculetype ]") is not -1) and ((ln[x+3].find("WAT") is not -1))): 
                f.write("#ifdef PROT \n")
                f.write("#include \""+path+"/"+prefix+"protposre.itp\" \n")
                f.write("#endif \n")
                f.write("#ifdef LIP \n")
                f.write("#include \""+path+"/"+prefix+"protlipposre.itp\"  \n")
                f.write("#endif \n")
            else:
                f.write(ln[x])
        else:
            f.write(ln[x])   
    #create the position restraint file
    posrewrite(prefix+"protposre.itp",path,trj.top.select("backbone"),1000,1000,1000)
    posrewrite(prefix+"protlipposre.itp",path,trj.top.select("(backbone or (resname == 'CHL') or (resname == 'OL') or (resname == 'PA') or (resname == 'PC')) and (not (name =~ 'H.*'))"),1000,1000,1000)    
    emwrite(path,prefix+"em.mdp")
    nvtposrewrite(path,prefix+"nvt.mdp",300.0)
    nptposrewrite(path,prefix+"npt.mdp",300.0,1.0)
    smdposrewrite(path,prefix+"smd.mdp",300.0,1.0)
    print("#############################")
    print("Finishing submitting the script file "+path+"/"+prefix+".qsub")   
    print("#############################") 
    return  

#if __name__ == "__main__":
def execMD(peprange,wd,acpype,GMXpath,ntomp,mpicall,indname,HPCtype,groupid):
    conc=0.15
    #acpype="../acpype/scripts/acpype.py"
    gmxpathserial=GMXpath
    gmxpathparallel=GMXpath
    #ntomp=6 #needed for system run setting with openmp larger than 6 cores
    #mpicall="mpijob -np 12 " #calling the mpi process
    #os.system("source ~/software/amber16/amber.sh") 
    superkonjobid=[]
    for x in range(0,range(peprange)): 
        for y in range(1,6):
            if os.path.isfile(wd+"/p3-"+str(x)+"/"+"model"+str(y)+".pdb"):  
                jobname="amp"+str(x)+"-"+str(y)
                superkonjobid.append(jobname)
                insertion(wd+"/p3-"+str(x),"model"+str(y)+".pdb","MDutil","lipid.pdb",wd+"/p3-"+str(x),conc,jobname)  
                #os.system("qsub /work/k0055/k005503/AMP-design/pepmcts-8Feb2019/actor-critic/test3-0_999/p3-"+str(x)+"/model"+str(y)+"protlip.qsub") 
                #writeissprunrestrt("/work/k0055/k005503/AMP-design/pepmcts-8Feb2019/actor-critic/test3-0_999/p3-"+str(x),"model"+str(y)+"protlip.qsub","model"+str(y)+"protlip",jobname )   
                #os.system("qsub /work/k0055/k005503/AMP-design/pepmcts-8Feb2019/actor-critic/test3-0_999/p3-"+str(x)+"/model"+str(y)+"protlip.qsub")                 
                if HPC==True:
                    #write run script
                    writequeuescript(indname,wd,HPCtype,groupid,numnode,mpiproc,mpproc,timelim)
                    writerun(path,prefix+".qsub",prefix,jobname)  
                    #submit job to HPC
                    HPCsubmit(indname,wd,HPCtype,groupid)
                else:
                    writerun(path,prefix+".qsub",prefix,jobname) 
                    os.system("sh "+path+"/"+prefix+".qsub")
    return superkonjobid
    

