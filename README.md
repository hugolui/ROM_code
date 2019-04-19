# DNN-ROM code (In construction)
The DNN-ROM code can be used for construction of reduced order models, ROMs, of fluid flows. The code employs a combination of flow modal decomposition and regression analysis. Spectral proper orthogonal decomposition, SPOD, is applied to reduce the dimensionality of the model and, at the same time, filter the POD temporal modes. The regression step is performed by a deep feedforward neural network, DNN, and the current framework is implemented in a context similar to the sparse identification of non-linear dynamics algorithm, SINDy. Test cases such as the compressible flow past a cylinder and the turbulent flow computed by a large eddy simulation of a plunging airfoil under dynamic stall are provided. 

# Required packages 
Your system will need the following packages to run the code:
1. Gfortran compiler
2. Python 3.6, including these modules
  - Numpy
  - Matplotlib
  - Scipy
  - Scikit-optimize
3. Tensorflow 
4. CGNS 
5. Surrogate modeling toolbox (SMT) 

# Installation 
Install the required packages and then clone the repository from github

# Example usage - Compressible flow past a cylinder



