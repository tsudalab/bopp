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

