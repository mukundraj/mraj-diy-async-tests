# Particle tracing test case for DIY Iexchange

Iexchange version can be built by setting the macro variable "IEXCHANGE" equal to 1 in examples/particle-tracing/ptrace.cpp. The block synchronous version can be built by setting "IEXCHANGE" equal to 0.

# Requirements

The latest version of DIY, and hence these examples, require a C++11 compiler. VTK 7.1 is required for visualization of streamlines only.

# Examples

- particle-tracing: particle tracing with DIY.

# Installation

## Build dependencies

a. DIY

```
git clone https://github.com/diatomic/diy
```

b. VTK (optional, for visualizing streamlines only)

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
cmake /path/to/mraj-diy-async-tests \
-DCMAKE_CXX_COMPILER=/path/to/mpicxx \
-DCMAKE_C_COMPILER=/path/to/mpicc \
-DCMAKE_INSTALL_PREFIX=/path/to/mraj-diy-async-tests/install \
-DDIY_INCLUDE_DIRS=/path/to/diy/include \
-DPNETCDF_DIR=/path/to/pnetcdf \
-DWITH_VTK=true \                   # optional
-DVTK_DIR=/path/to/vtk/build \      # optional, needed if -DWITH_VTK=true

make install
```

# Execution



```

## Particle tracing example
```
cd path/to/mraj-diy-async-tests/install/examples/particle_tracing
./TORNADO_TEST
./PLUME_TEST
./NEK_TEST1

```
