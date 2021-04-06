# bopp
Black-box optimization of peptides and proteins

This package is a part of black-box optimization of peptides and proteins using actor-critic model.
The algorithm will help to increase the confidency in sequence generating of neural network by feeding it with more and more sequence validated by MD simulations or other methods by user-defined.
The option for using MD simulation based validation is generally-talking time-consuming. We encourage users use it on HPC.


# Prerequisite:
- Tensorflow (currently using TF v2.+)
- Keras
- Numpy
- Multiprocessing

MD evaluation requires:
- I-TASSER
- GROMACS
- AMBERTOOLs


# Before running:
- Setting the parameters in parameter.py is required. There are explanations inside the file.
- We provide the setup file (setup.py) to help initialization such as installing GROMACS and AMBERTOOls. To run: ``` python setup.py ```
- The default settings of this tool is MD-simulation evaluation based. We provide the subroutine so-called `evaluate` in `sampling.py` for user-defined evaluating function.
- We have set up the queue format for some supercomputer such as: OakForest (UTokyo),, ISSP (UTokyo), Tsubame3 (TokyoTech), KComputer alike (IMS, NINS), ITO (UKyushu). Users can add additional submit script to HPCtype in the `HPCtool.py` files.

# To run:

just call: 
```
python sampling.py
```

The script will gentlely sample the sequence matching with the design goal.
Each running epoch: it will generate the selected sequences and report:

```
equence TCNSQSIAFD with probability at 0.88409495
Sequence MYTFTMWKNW with probability at 0.926888
Sequence MYTFTMWKNW with probability at 0.926888
Sequence NFEQFFYAQA with probability at 0.8918185
Sequence KNNFMWICYG with probability at 0.9766589
Sequence IRFHVSQSVH with probability at 0.96453977
Sequence LRFNKQAKPS with probability at 0.9339483
Sequence KYGIFNGGCF with probability at 0.97861457
Sequence NFEQFFYAQA with probability at 0.8918185
Sequence LRFNKQAKPS with probability at 0.9339483
Sequence KYGIFNGGCF with probability at 0.97861457
Sequence YGEIRNVKER with probability at 0.97028124
Sequence TCNSQSIAFD with probability at 0.88409495
Dealing with epoch 8
Sequence NPPDADSEQF with probability at 0.7778874
Sequence YGEIRNVKER with probability at 0.97028124
Sequence NPPDADSEQF with probability at 0.7778874
Sequence MWMWEFLNTY with probability at 0.9251394
Sequence TCNSQSIAFD with probability at 0.88409495
Sequence MYTFTMWKNW with probability at 0.926888
Sequence MWMWEFLNTY with probability at 0.9251394
Sequence KNNFMWICYG with probability at 0.9766589
```

After sampling the model, the selected sequence will then homology modeling by I-TASSER, you will find the generated pdb files (usually 5 models) in your selected working folder.

If MD-evaluated is used, the model peptides will embeded into a membrane and execute the Steered MD through out the predefined lipid bilayers `MDUtil/lipidnew.pdb`.

After running the MD simulation, the output time-course steering forces are reported in the corresponding folder.


