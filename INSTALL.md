# Installing the required packages 

Before installing DNN-ROM you'll need to install the following packages
1. GFortran compiler
2. Python 3.6, including these modules
    - Numpy
    - Matplotlib
    - Scipy
    - Scikit-optimize
    - Surrogate modeling toolbox (SMT) 
3. Tensorflow (GPU or CPU support)
4. CGNS 

## Installing GFortran 
	$ sudo apt-get install gfortran

## Installing numpy
	$ pip3 install numpy

## Installing matplotlib
	$ pip3 install matplotlib

## Installing scipy
	$ pip3 install scipy

## Installing scikit-Optimize
	$ pip3 install scikit-optimize

## Installing Surrogate Modeling Toolbox

	- Clone the repository from github (https://github.com/SMTorg/SMT) then run:

	$ pip3 install -e <smt_folder>

## Installing Tensorflow GPU support

	(a) Update your GPU driver if your version is below 410.x (CUDA 10.0 requires 410.x or higher)
	(b) Install the CUDA toolkit 10.0
		(b.1) Download the installer https://developer.nvidia.com/cuda-10.0-download-archive
		(b.2) Open a terminal in the CUDA toolkit directory and run the following commands
			$ sudo sh cuda_10.0.130_410.48_linux.run
			$ Follow the command-line prompts
			$ Download and install the patches if available
		(b.3) Update your PATH variable
			$ gedit ~/.bashrc
			Go to the last line and add the following lines
			export PATH=/usr/local/cuda-10.0/bin${PATH:+:${PATH}}
			export LD_LIBRARY_PATH=/usr/local/cuda-10.0/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
			export LD_LIBRARY_PATH=/usr/local/cuda-10.0/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
			
	(c) Install the cuDNN (version >= 7.4.1)
		(c.1) Go to https://developer.nvidia.com/cudnn and Download the cuDNN version for CUDA 10.0
		(c.2) Unzip the cuDNN package
			$ tar -xzvf cudnn-10.0-linux-x64-v7.6.1.34.tgz
		(c.3) Copy the following files into the CUDA Toolkit directory, and change the file permissions
			$ sudo cp cuda/include/cudnn.h /usr/local/cuda/include
			$ sudo cp cuda/lib64/libcudnn* /usr/local/cuda/lib64
			$ sudo chmod a+r /usr/local/cuda/include/cudnn.h /usr/local/cuda/lib64/libcudnn*
		(c.4) Update your PATH variable
			$ gedit ~/.bashrc
			Go to the last line and add the following lines
			export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/extras/CUPTI/lib64

	(d) Install tensorflow
		$ pip3 install tensorflow-gpu
		
## Installing CGNS

	1. Clone the repository from github (https://github.com/CGNS/CGNS) and then go to the CGNS directory
		(a) In the bin (../CGNS/bin) folder, there is a file called "config-cgns.sh". Copy this file to the ../CGNS directory 
		$ cp bin/config-cgns.sh .

	(b) Modify the config-cgns.sh to the following form:

		#!/bin/sh
		#
		# Configure CGNS for travis CI. 
		#
		set -e
		cd src
		OPTS=""
		if [ $TRAVIS_OS_NAME = "linux" ]; then
		  export FLIBS="-Wl,--no-as-needed -ldl"
		  export LIBS="-Wl,--no-as-needed -ldl"
		  OPTS="--enable-parallel"
		fi

		export CC="gcc"
		export FC="gfortran"

		./configure \
		--prefix=$PWD/cgns_build $OPTS \
		--with-hdf5=no \
		--with-fortran \
		--enable-lfs \
		--disable-shared \
		--enable-debug \
		--disable-cgnstools \

	(c) run the config-cgns.sh file
		$ sh config-cgns.sh 	
			
	(d) Build the CGNS libray using 
		$ cd src/
		$ make
		$ make install
		$ make test

	(e) Copy the CGNS build files to /usr/local/lib and /usr/local/lib
		$ cd cgns_build/
		$ sudo cp include/cgns* /usr/local/include/
		$ sudo cp lib/libcgns.a /usr/local/lib/
