import unittest,os
from pathlib import Path

import pkg_resources

#_REQUIREMENTS_PATH = Path(__file__).parent.with_name("requirements.txt")


class TestRequirements(unittest.TestCase):
    def test_requirements(self):
        requirements = pkg_resources.parse_requirements("requirements.txt")
        for requirement in requirements:
            requirement = str(requirement)
            with self.subTest(requirement=requirement):
                pkg_resources.require(requirement)

#checking whether GROMACS is working or not.
GROtest=False
print("Testing whether GROMACS is installed")
for x in "gmx_mpi" "gmx" "gmx_mpi_d" "gmx_gpu":
    os.system(x+" -h > test")
    with open('test') as f:
        ln=f.read()
        if 'Executable' in ln:
            GROtest=True
            GROpath=ln.strip().split()[1]
            print("Modify the GROMACS command in parameter.py file with "+GROpath)
            break
        else:
            print(x+" : failed ...")

pwd=os.getcwd()
GROtest=False
if GROtest==False:
    os.system("wget http://ftp.gromacs.org/pub/gromacs/gromacs-2021.1.tar.gz ")
    os.system("tar -xvzf gromacs-2021.1.tar.gz ")
    os.system("mkdir gromacs-2021.1/build")
    os.system("mkdir soft")
    os.system("mkdir soft/gromacs")
    #os.chdir("gromacs-2019.6/build")
    os.system("cmake -S "+pwd+"/gromacs-2021.1  -B "+pwd+"/gromacs-2021.1/build -DCMAKE_INSTALL_PREFIX="+pwd+"/soft/gromacs -DGMX_FFT_LIBRARY=fftw3 -DGMX_BUILD_OWN_FFTW=ON -DGMX_MPI=ON -DGMX_GPU=CUDA; ")
    os.chdir(pwd+"/gromacs-2021.1/build")
    os.system("make -j4; make install")
    os.chdir(pwd)

#checking whether AMBER is working or not.
print("Testing whether tleap, amb2gro_top_gro.py in AMBERTools package is installed")
tleaptest=False
for x in "tleap":
    os.system(x+" -h > test")
    with open('test') as f:
        ln=f.read()
        if 'Usage' in ln:
            tleaptest=True
            tleappath=ln.strip().split()[1]
        else:
            print(x+" : failed ...")
#clean trash
os.system("rm -r test")

amb2grotest=False
for x in "tleap":
    os.system(x+" -h > test")
    with open('test') as f:
        ln=f.read()
        if 'Usage' in ln:
            amb2grotest=True
            amb2gropath=ln.strip().split()[1]
        else:
            print(x+" : failed ...")
#clean trash
os.system("rm -r test")

pwd=os.getcwd()
ambertooltest=amb2grotest or tleap
if ambertooltest==False:
    os.system("conda install -c conda-forge ambertools=20 ")
#clean trash
os.system("rm -r test")

