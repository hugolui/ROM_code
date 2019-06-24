# DNN-ROM code 
The DNN-ROM code can be used for construction of reduced-order models, ROMs, of fluid flows. The code employs a combination of flow modal decomposition and regression analysis. Spectral proper orthogonal decomposition, SPOD, is applied to reduce the dimensionality of the model and, at the same time, filter the POD temporal modes. The regression step is performed by a deep feedforward neural network, DNN, in order to obtain a system of ordinary differential equations (ODEs) which represents the dynamics of POD temporal modes. The current framework is implemented in a context similar to the sparse identification of non-linear dynamics algorithm, SINDy. Test cases such as the compressible flow past a cylinder and the turbulent flow computed by a large eddy simulation of a plunging airfoil under dynamic stall are provided. 

For more details, see

Lui, H., & Wolf, W. (2019). [Construction of reduced-order models for fluid flows using deep feedforward neural networks](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/construction-of-reducedorder-models-for-fluid-flows-using-deep-feedforward-neural-networks/ECEC52E32AEBBEA049CF26D6C79EE394). Journal of Fluid Mechanics, 872, 963-994. doi:10.1017/jfm.2019.358                                                               

# Required packages 
Your system will need the following packages to run the code:
1. GFortran compiler
2. Python 3.6, including these modules
    - Numpy
    - Matplotlib
    - Scipy
    - Scikit-optimize
    - Surrogate modeling toolbox (SMT) 
3. Tensorflow (GPU or CPU support)
4. CGNS 

# Installation 
See INSTALL.md 

# Example 1 - Compressible flow past a cylinder

In this example, the full order model (FOM) is obtained by solving the compressible Navier Stokes equations as detailed in [section 4.2](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/construction-of-reducedorder-models-for-fluid-flows-using-deep-feedforward-neural-networks/ECEC52E32AEBBEA049CF26D6C79EE394). The numerical simulations are conducted for Reynolds and Mach numbers Re = 150 and M = 0.4, respectively. The grid configuration consists of a body-fitted O-grid with 421 × 751 points in the streamwise and wall-normal directions, respectively. The training data comprises the first 280 snapshots of the FOM, and the validation data contains the next 145 snapshots. 

We are going to build a reduced order model capable to predict the flowfield beyond the training. The final form of the ROM is an ordinary differential equation (ODE) system. 

## 1 - Download the cylinder data
First, download the cylinder data: 

https://www.dropbox.com/sh/ji6i5u8valqyda8/AABtEYaZ7vG-62q6h4toQKBTa?dl=0

## 2 - Define the parameters in the inputs.inp file
In the folder "code", we have a file called "inputs.inp". In this file, we specify the information required to construct the reduced order model. There are comments explaining the inputs list. In order to run this example, we just need to modify lines 1 and 3.

    - In line 1, specify the path to the cylinder data. For example: /home/cfd/Desktop/hugo/cylinder_data/
    - In line 3, specify the path to the folder "code". For example: /home/cfd/Desktop/hugo/ROM_code/code/
    
## 3 - Run the code    
Now, we can run the code by open a terminal in the folder "code" and typing 
      
    sh run.sh 

The output file is a CGNS file containing the ROM solutions for a determined number of snapshots (listed in the inputs.inp file). For all candidate models evaluated, the DNN parameters and hyperparameters can be found in the folder ".../regression/deep_learning/results/". 

Videos comparing the FOM and ROM solutions:
    
https://www.youtube.com/watch?v=52n4CP-01a8<br/>
https://www.youtube.com/watch?v=e8vASdxMUg8
    
In the "inputs.inp" file, we can specify the:

    training and validation data;
    fluid region of interest for the construction of the ROM; 
    numerical scheme used to compute the derivative of the temporal modes;
    number of POD modes; 
    SPOD filter size and type;
    norm for the POD correlation matrix; 
    hyperparameters search space;
    hyperparameter optimization strategy; 
    number of candidate models to be evaluated;
    number of snapshots for reconstruction of the flowfield. 
 
 So, there are many parameters to play with here to improve the accuracy of the reduced order model.
   
# Example 2 - Deep dynamic stall of plunging airfoil

In this example, we are going to reconstruct a turbulent flow involving dynamic stall of a plunging airfoil. Further details about the FOM can be found in [section 4.2.2](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/construction-of-reducedorder-models-for-fluid-flows-using-deep-feedforward-neural-networks/ECEC52E32AEBBEA049CF26D6C79EE394). Insteady of training the neural network as in the example 1, it will be given the DNN parameters already trained as well the POD modes computed from the training data (the first two cycles).

## 1 - Download the dynamical stall data 
First, download the dynamical stall data to the directory ".../ROM_code/examples/dynamical_stall/data/"

https://www.dropbox.com/sh/6qjz6eyexp1m5yu/AABo4A0_btk-ciKDS2DRsoiTa?dl=0
    
## 2 - Predict the POD temporal modes
In the folder ".../ROM_code/examples/dynamic_stall/", there is a python code called "reconst_best_model.py" that reconstruct the temporal modes. As we have the DNN parameters and the initial conditions, we can solve the system of coupled ODEs that describe the dynamics of the temporal modes by runnning the python code

    python3 reconst_best_model.py
    
The output is the POD temporal modes. 

## 3 - Reconstruct the flowfield

Once we have the POD temporal and spatial modes, we can reconstruct the flowfield by running the following commands:

    sh compile.sh
    ./reconst.out
    
The output file is a CGNS file containing the ROM solutions for a determined number of snapshots. In the "reconst_best_model.py" code, we can set the time step and the number of snapshots for reconstruction, hence the reduced order model obtained allows to predict the flowfield beyond the training window with larger or smaller time increments than those used by the full order model.

Videos comparing the FOM and ROM solutions:
    
https://www.youtube.com/watch?v=LWS98UX-YqA<br/>
https://www.youtube.com/watch?v=sxgaJTJYCls<br/>
https://www.youtube.com/watch?v=eE3Lem1EX8w (3D case, not provided in this example)

