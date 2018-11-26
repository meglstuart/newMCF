# newMCF

clone:
git clone --recurse-submodules git@github.com:meglstuart/newMCF.git

build:
mkdir build
cd build
ccmake ..
(ensure LIBIGL_WITH_OPENGL_GLFW_IMGUI is on)
make
