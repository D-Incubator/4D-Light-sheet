# 4D-Light-sheet
4D image reconstruction and VR model generation for the paper "4D Light-Sheet Imaging and Interactive Analysis of Cardiac Contractility in Zebrafish Larvae".

## This package contains:
1. The MATLAB code for reconstructing 4D images using parallel computation.
2. The Python code for converting .tif images with mutiple cell labels into .obj files for VR analysis in Unity.

## Requirements
- MATLAB
- Jupyter Notebook (Python 3)
- 3D Slicer

## Usage
### Run the 4D image reconstruction code 
0. The platform is MATLAB(R). [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html) is required for parallel computation.
1. Download the [4D-Reconstruction](https://github.com/d-incubator1/4D-Light-sheet/tree/main/4D-Reconstruction) code.
2. Download the [example data](https://drive.google.com/file/d/1b6wrtxzxBQMfubPJYMH8rC0Q6WtakY-p/view?usp=share_link) on google drive and put the data in ```4D-Reconstruction/example```.
3. Run the ```test_Parallel.m``` file and check the output in the data folder.
### Run the python code for converting .tif labels to .obj models
0. Install Jupyter Notebook (using [Anaconda](https://www.anaconda.com/)) and [3D Slicer](https://www.slicer.org/) on the computer.
1. Install ```SlicerJupyter``` extension in 3D slicer following this [instruction](https://github.com/Slicer/SlicerJupyter).
2. Download the [example data](https://drive.google.com/file/d/1lTHBetnFLjHKUVjT5oqtXGFTF00-oI9r/view?usp=sharing) and the [code](https://github.com/d-incubator1/4D-Light-sheet/blob/main/CellLabelsToObj.ipynb).
3. Open the code in Jupyter Notebook. Change the kernel to ```Slicer 5.2``` (I used this but other versions may work as well). 
4. Set the ```dataFolder``` and ```saveFolder``` in the code for sample data and output data.
5. Run the code and check the output.
