# ODT
Code for "Optimal Decision Trees in Economics: a Framework to construct Binary Decision Trees with Categorical Data".

All software is written in Python3 and makes use of IBM ILOG CPLEX Optimisation studio (https://www.ibm.com/products/ilog-cplex-optimization-studio).

Collectively, the files in this repository are used to generate the numerical results for my bachelor thesis.

The main programme is 'ODT.py', which constructs decision trees and executes the optimisation using IBM CPLEX.\
This file needs 'DataProcessing.py' to read the data files and fetch all the parameters.\
All other files are data files.

In the ODT file, you can edit in the "settings" section certain specifications: the dataset, the tree topology and the options for the optimisation.

Kind regards,

Martijn Korf

2020-06-04
