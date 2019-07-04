# bopp
Black-box optimization of peptides and proteins

This package is a part of black-box optimization of peptides and proteins using actor-critic model.
The algorithm will help to increase the confidency in sequence generating of neural network by feeding it with more and more sequence validated by MD simulations or other methods by user-defined.
The option for using MD simulations validated is generally-talking time-consuming. We encourage users use it on HPC.

Prerequisite:
- Tensorflow
- Keras
- Numpy
- Multiprocessing

MD evaluation requires:
- I-TASSER
- GROMACS
- AMBERTOOLs
- This option is time-consuming, considering using on HPC is saving-life.

Before running:
- Setting the parameters in parameter.py is required. There are explanation inside the file.

To run:

just call: 

`python sampling.py`

