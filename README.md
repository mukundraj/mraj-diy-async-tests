# Integration of DIY with VTK

A set of examples of vtk filters and data models written using DIY block-parallelism.

# Requirements

The latest version of DIY, and hence these examples, require a C++11 compiler. VTK 7.1 is required.

# Examples

- vtk-data-models: various VTK data models wrapped in DIY blocks
- particle-tracing: vtkStreamTracer with DIY.

# Installation

## Build dependencies

a. DIY

```
git clone https://github.com/diatomic/diy
```

b. VTK

```
wget http://www.vtk.org/files/release/7.1/VTK-7.1.0.tar.gz
tar -xvf VTK-7.1.0.tar.gz

// --or-- //

git clone https://gitlab.kitware.com/vtk/vtk.git

cmake /path/to/vtk \
-DCMAKE_INSTALL_PREFIX=/Users/tpeterka/software/vtk/install \
-DModule_vtkRenderingParallel=on \
-DVTK_Group_MPI=on \
-DModule_vtkParallelMPI=on \

make
make install
```

## Build examples

```
cmake /path/to/diy2-vtk7 \
-DCMAKE_CXX_COMPILER=/path/to/mpicxx \
-DCMAKE_C_COMPILER=/path/to/mpicc \
-DCMAKE_INSTALL_PREFIX=/path/to/diy2-vtk7/install \
-DDIY_INCLUDE_DIRS=/path/to/diy/include \
-DVTK_DIR=/path/to/vtk/build \

make install
```

# Execution

## Data model example
```
cd path/to/diy2-vtk7/install/examples/vtk-data-models
mpiexec -n <num_procs> ./diy-rgrid

```

## Particle tracing example
```
cd path/to/diy2-vtk7/install/examples/particle_tracing
./TORNADO_TEST
./PLUME_TEST

```
