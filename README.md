# newMCF

## clone:
`git clone --recurse-submodules git@github.com:meglstuart/newMCF.git`

## build:
```
mkdir build
cd build
ccmake ..
(ensure LIBIGL_WITH_OPENGL_GLFW_IMGUI is on)
make
```
## to run:
for example

`./mcf_srep_bin ../003_left_parotid.stl 0.001 0.0001 1000`

press space to step the algorithm
this will give you a radii.txt

## Now run the matlab script
First, install Medical Image Processing Toolkit: https://www.mathworks.com/matlabcentral/fileexchange/41594-medical-image-processing-toolbox

May have to add MedicalImageProcessingToolbox and subfolders to matlab path (right click on folder in view at left of screen)
