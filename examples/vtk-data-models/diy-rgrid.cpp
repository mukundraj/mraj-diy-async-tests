/*=========================================================================

  Program:   Visualization Toolkit
  Module:    RGrid.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/
// This example shows how to create a rectilinear grid.
//

#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkFloatArray.h"
#include "vtkRectilinearGrid.h"
#include "vtkRectilinearGridGeometryFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"

#include <diy/master.hpp>
#include <diy/reduce.hpp>
#include <diy/partners/merge.hpp>
#include <diy/decomposition.hpp>
#include <diy/assigner.hpp>

#include "../opts.h"

using namespace std;

typedef  diy::ContinuousBounds       Bounds;
typedef  diy::RegularContinuousLink  RCLink;

static float x[47] =
{
    -1.22396, -1.17188, -1.11979, -1.06771, -1.01562, -0.963542,
    -0.911458, -0.859375, -0.807292, -0.755208, -0.703125, -0.651042,
    -0.598958, -0.546875, -0.494792, -0.442708, -0.390625, -0.338542,
    -0.286458, -0.234375, -0.182292, -0.130209, -0.078125, -0.026042,
    0.0260415, 0.078125, 0.130208, 0.182291, 0.234375, 0.286458,
    0.338542, 0.390625, 0.442708, 0.494792, 0.546875, 0.598958,
    0.651042, 0.703125, 0.755208, 0.807292, 0.859375, 0.911458,
    0.963542, 1.01562, 1.06771, 1.11979, 1.17188
};
static float y[33] =
{
    -1.25, -1.17188, -1.09375, -1.01562, -0.9375, -0.859375,
    -0.78125, -0.703125, -0.625, -0.546875, -0.46875, -0.390625,
    -0.3125, -0.234375, -0.15625, -0.078125, 0, 0.078125,
    0.15625, 0.234375, 0.3125, 0.390625, 0.46875, 0.546875,
    0.625, 0.703125, 0.78125, 0.859375, 0.9375, 1.01562,
    1.09375, 1.17188, 1.25
};
static float z[44] =
{
    0, 0.1, 0.2, 0.3, 0.4, 0.5,
    0.6, 0.7, 0.75, 0.8, 0.9, 1,
    1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
    1.7, 1.75, 1.8, 1.9, 2, 2.1,
    2.2, 2.3, 2.4, 2.5, 2.6, 2.7,
    2.75, 2.8, 2.9, 3, 3.1, 3.2,
    3.3, 3.4, 3.5, 3.6, 3.7, 3.75,
    3.8, 3.9
};

// block
struct Block
{
    Block()                                                     {}
    static void*    create()                                    { return new Block; }
    static void     destroy(void* b)                            { delete static_cast<Block*>(b); }
    static void     save(const void* b, diy::BinaryBuffer& bb)
        { diy::save(bb, *static_cast<const Block*>(b)); }
    static void     load(void* b, diy::BinaryBuffer& bb)
        { diy::load(bb, *static_cast<Block*>(b)); }
    void generate_data()
        {
            // Create a rectilinear grid by defining three arrays specifying the
            // coordinates in the x-y-z directions.
            xCoords = vtkFloatArray::New();
            for (int i = 0; i < 47; i++)
                xCoords->InsertNextValue(x[i]);

            yCoords = vtkFloatArray::New();
            for (int i = 0; i < 33; i++)
                yCoords->InsertNextValue(y[i]);

            zCoords = vtkFloatArray::New();
            for (int i = 0; i < 44; i++)
                zCoords->InsertNextValue(z[i]);

            // The coordinates are assigned to the rectilinear grid. Make sure that
            // the number of values in each of the XCoordinates, YCoordinates,
            // and ZCoordinates is equal to what is defined in SetDimensions().
            //
            rgrid = vtkRectilinearGrid::New();
            rgrid->SetDimensions(47,33,44);
            rgrid->SetXCoordinates(xCoords);
            rgrid->SetYCoordinates(yCoords);
            rgrid->SetZCoordinates(zCoords);
        }

    int gid;
    Bounds bounds;                    // block bounding box
    // vtk data
    vtkFloatArray *xCoords;
    vtkFloatArray *yCoords;
    vtkFloatArray *zCoords;
    vtkRectilinearGrid *rgrid;
};

// add blocks to a master
struct AddBlock
{
    AddBlock(diy::Master& master_):
        master(master_) {}

    void operator()(int gid, const Bounds& core, const Bounds& bounds, const Bounds& domain,
                    const RCLink& link) const
        {
            Block*        b = new Block();
            RCLink*       l = new RCLink(link);
            diy::Master&  m = const_cast<diy::Master&>(master);
            m.add(gid, b, l);
            b->gid = gid;
            b->bounds = bounds;              // TODO: is this an actual deep copy (desired)?
            b->generate_data();
        }

    diy::Master&  master;
};

// debug: block foreach function to print the block
void print_block(Block* b,
                 const diy::Master::ProxyWithLink& cp)
{
    fprintf(stderr, "gid %d bounds min[%.3f %.3f %.3f] max[%.3f %.3f %.3f]\n", b->gid,
            b->bounds.min[0], b->bounds.min[1], b->bounds.min[2],
            b->bounds.max[0], b->bounds.max[1], b->bounds.max[2]);
}

// block foreach function to render the block
void render_block(Block* b,
                  const diy::Master::ProxyWithLink& cp)
{
    // Extract a plane from the grid to see what we've got.
    vtkRectilinearGridGeometryFilter *plane = vtkRectilinearGridGeometryFilter::New();
    plane->SetInputData(b->rgrid);
    plane->SetExtent(0,46, 16,16, 0,43);

    vtkPolyDataMapper *rgridMapper = vtkPolyDataMapper::New();
    rgridMapper->SetInputConnection(plane->GetOutputPort());

    vtkActor *wireActor = vtkActor::New();
    wireActor->SetMapper(rgridMapper);
    wireActor->GetProperty()->SetRepresentationToWireframe();
    wireActor->GetProperty()->SetColor(0,0,0);

    // Create the usual rendering stuff.
    vtkRenderer *renderer = vtkRenderer::New();
    vtkRenderWindow *renWin = vtkRenderWindow::New();
    renWin->AddRenderer(renderer);
    vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin);

    renderer->AddActor(wireActor);
    renderer->SetBackground(1,1,1);
    renderer->ResetCamera();
    renderer->GetActiveCamera()->Elevation(60.0);
    renderer->GetActiveCamera()->Azimuth(30.0);
    renderer->GetActiveCamera()->Zoom(1.0);

    renWin->SetSize(300,300);

    // interact with data
    renWin->Render();
    iren->Start();

    // cleanup
    plane->Delete();
    renderer->Delete();
    renWin->Delete();
    iren->Delete();
    b->xCoords->Delete();
    b->yCoords->Delete();
    b->zCoords->Delete();
    b->rgrid->Delete();
    rgridMapper->Delete();
    wireActor->Delete();
}

int main(int argc, char* argv[])
{
    int tot_blocks;                     // total number of blocks in the domain
    int num_threads;                    // number of threads diy can use
    int mem_blocks;                     // number of blocks to keep in memory
    int dim;                            // domain dimensionality
    int rank, size;                     // MPI usual

    diy::mpi::environment     env(argc, argv); // equivalent of MPI_Init
    diy::mpi::communicator    world;           // equivalent of MPI communicator

    rank = world.rank();
    size = world.size();

    Bounds domain;                           // not used in this example, just initialize to [0,1]
    domain.min[0] = domain.min[1] = domain.min[2] = 0.0;
    domain.max[0] = domain.max[1] = domain.max[2] = 1.0;

    // get/set options

    using namespace opts;

    tot_blocks    = size;
    num_threads   = 4;
    mem_blocks    = -1;
    dim           = 3;
    string prefix = "./DIY.XXXXXX";

    Options ops(argc, argv);

    ops
        >> Option('b', "blocks",    tot_blocks,   "Total number of blocks to use")
        >> Option('t', "threads",   num_threads,  "Number of threads to use")
        >> Option('m', "in-memory", mem_blocks,   "Number of blocks to keep in memory")
        >> Option('s', "storage",   prefix,       "Path for out-of-core storage")
        ;

    if ( ops >> Present('h', "help", "show help"))
    {
        if (rank == 0)
        {
            fprintf(stderr, "Usage: %s [OPTIONS]\n", argv[0]);
            std::cout << ops;
        }
        return 1;
    }

    // initialize DIY
    diy::FileStorage          storage("./DIY.XXXXXX");
    diy::Master               master(world,
                                     num_threads,
                                     mem_blocks,
                                     &Block::create,
                                     &Block::destroy,
                                     &storage,
                                     &Block::save,
                                     &Block::load);
    diy::ContiguousAssigner   assigner(world.size(), tot_blocks);
    AddBlock                  create(master);
    diy::decompose(dim, world.rank(), domain, assigner, create);

    // debug
    master.foreach(print_block);

    master.foreach(render_block);

    return 0;
}
