# DNN-ROM code (In construction)
The DNN-ROM code can be used for construction of reduced order models, ROMs, of fluid flows. The code employs a combination of flow modal decomposition and regression analysis. Spectral proper orthogonal decomposition, SPOD, is applied to reduce the dimensionality of the model and, at the same time, filter the POD temporal modes. The regression step is performed by a deep feedforward neural network, DNN, and the current framework is implemented in a context similar to the sparse identification of non-linear dynamics algorithm, SINDy. Test cases such as the compressible flow past a cylinder and the turbulent flow computed by a large eddy simulation of a plunging airfoil under dynamic stall are provided. For more details, see https://arxiv.org/abs/1903.05206. 

# Required packages 
Your system will need the following packages to run the code:
1. GFortran compiler
2. Python 3.6, including these modules
    - Numpy
    - Matplotlib
    - Scipy
    - Scikit-optimize
    - Surrogate modeling toolbox (SMT) 
3. Tensorflow 
4. CGNS 

# Installation 
Install the required packages and then clone the repository from github.

# Example 1 - Compressible flow past a cylinder

First, we need to download the cylinder data: https://www.dropbox.com/sh/ji6i5u8valqyda8/AABtEYaZ7vG-62q6h4toQKBTa?dl=0. 

In the folder "code", we have a file called "inputs.inp". In this file, we specify the information required to construct the reduced order model. There are comments explaining the inputs list. In order to run this example, we just need to modify lines 1 and 3.

    - In line 1, we specify the path to the cylinder data. For example: /home/cfd/Desktop/hugo/cylinder_data/
    - In line 3, we specify the path to the folder "code". For example: /home/cfd/Desktop/hugo/ROM_code/code/
    
Now, we can run the code by open a terminal in the folder "code" and typing sh run.sh. The output file is a CGNS file containing the ROM solution for a determined number of snapshots (listed in the inputs.inp file). For all candidate models evaluated, the DNN parameters and hyperparameters can be found in the folder "line3/regression/deep_learning/results/". 

In the "inputs.inp" file, we can specify the training and validation data, the fluid region of interest for the construction of the ROM, the numerical scheme used to compute the derivative of the temporal modes, the number of POD modes, the SPOD size and type, the norm for the POD correlation matrix, the hyperparameters search space, the hyperparameter optimization strategy, the number of candidate models to be evaluated and the number of snapshots for reconstruction of the flowfield. So, there are many parameters to play with here to improve the accuracy of the reduced order model.
   
# Example 2 - Deep dynamic stall of plunging airfoil
