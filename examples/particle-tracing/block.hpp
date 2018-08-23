//---------------------------------------------------------------------------
//
// diy2-vtk7 parallel particle advection block class
//
// Tom Peterka
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// tpeterka@mcs.anl.gov
//
//--------------------------------------------------------------------------
#include <diy/mpi.hpp>
#include <diy/master.hpp>
#include <diy/reduce-operations.hpp>
#include <diy/decomposition.hpp>
#include <diy/mpi/datatypes.hpp>
#include <diy/io/bov.hpp>
#include <diy/pick.hpp>
#include <diy/reduce.hpp>
#include <diy/partners/merge.hpp>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>

typedef diy::DiscreteBounds            Bounds;
typedef diy::RegularGridLink           RGLink;
typedef diy::RegularDecomposer<Bounds> Decomposer;

using namespace std;

// the diy block
struct Block
{
    Block() : nvecs(0), init(0), done(0) {}
    ~Block()
        {
            if (nvecs)
            {
                for (int i = 0; i < 3; i++)
                    delete vel[i];
            }
        }

    static void* create()
        {
            return new Block;
        }
    static void destroy (void* b)
        {
            delete static_cast<Block*>(b);
        }
    static void save(const void* b_, diy::BinaryBuffer& bb)
        {
            const Block* b = static_cast<const Block*>(b_);
            diy::save(bb, b->nvecs);
            diy::save(bb, b->vel[0], b->nvecs);
            diy::save(bb, b->vel[1], b->nvecs);
            diy::save(bb, b->vel[2], b->nvecs);
            diy::save(bb, b->init);
            diy::save(bb, b->done);
            // TODO: serialize vtk structures
        }
    static void load(void* b_, diy::BinaryBuffer& bb)
        {
            Block* b = static_cast<Block*>(b_);
            diy::load(bb, b->nvecs);
            b->vel[0] = new float[b->nvecs];
            b->vel[1] = new float[b->nvecs];
            b->vel[2] = new float[b->nvecs];
            diy::load(bb, b->vel[0], b->nvecs);
            diy::load(bb, b->vel[1], b->nvecs);
            diy::load(bb, b->vel[2], b->nvecs);
            diy::load(bb, b->init);
            diy::load(bb, b->done);
            // TODO: serialize vtk structures
        }

    // debug
        void show_geometry(const diy::Master::ProxyWithLink& cp, void*)
            {
                diy::RegularLink<Bounds> *link = static_cast<diy::RegularLink<Bounds>*>(cp.link());

                fprintf(stderr, "rank=%d, gid=%d, "
                        "{%d, %d, %d}x{%d, %d, %d}, {%d, %d, %d}x{%d, %d, %d}, %d\n",
                        cp.master()->communicator().rank(),
                        cp.gid(),
                        link->core().min[0], link->core().min[1], link->core().min[2],
                        link->core().max[0], link->core().max[1], link->core().max[2],
                        link->bounds().min[0], link->bounds().min[1], link->bounds().min[2],
                        link->bounds().max[0], link->bounds().max[1], link->bounds().max[2],
                        link->size());
                for (size_t i = 0; i < segments.size(); i++)
                {
                    fprintf(stderr, "[pid %d num_pts %ld]\n", segments[i].pid, segments[i].pts.size());
                    for (size_t j = 0; j < segments[i].pts.size(); j++)
                        fprintf(stderr, "[%.3f %.3f %.3f] ",
                                segments[i].pts[j].coords[0], segments[i].pts[j].coords[1],
                                segments[i].pts[j].coords[2]);
                    fprintf(stderr, "\n");
                }
            }

        // convert vector of particle trajectories to vtk polylines and render them
        // TODO: use vtk data model from the outset?
        // copied from http://www.vtk.org/Wiki/VTK/Examples/Cxx/GeometricObjects/PolyLine
        void render()
            {
                printf("%ld", segments.size());
                // vtk points
                vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

                // vtk cells
                vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
                size_t n = 0;           // index of current point in all points in all traces

                // add points from each trace to one global list of vtkPoints
                for (size_t i = 0; i < segments.size(); i++)
                    for (size_t j = 0; j < segments[i].pts.size(); j++)
                        points->InsertNextPoint(segments[i].pts[j].coords); // deep copy, I assume?

                // create a polyline from each trace and a cell from each polyline
                for (size_t i = 0; i < segments.size(); i++)
                {
                    // vtk polyline
                    vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
                    polyLine->GetPointIds()->SetNumberOfIds(segments[i].pts.size());

                    for(unsigned int j = 0; j < segments[i].pts.size(); j++)
                        // map index of point in the streamline to point in points geometry
                        // setId(id of point in this streamline, id of point in all points)
                        polyLine->GetPointIds()->SetId(j, n++);

                    // store the polyline in a 1d cell
                    cells->InsertNextCell(polyLine);
                }

                // Create a polydata to store points and cells
                vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();

                // Add the points to the dataset
                polyData->SetPoints(points);

                // Add the lines to the dataset
                polyData->SetLines(cells);

                // Setup actor and mapper
                vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
                mapper->SetInputData(polyData);

                vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
                actor->SetMapper(mapper);

                // Setup render window, renderer, and interactor
                vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
                renderer->AddActor(actor);

                vtkSmartPointer<vtkRenderWindow> renderWindow =
                    vtkSmartPointer<vtkRenderWindow>::New();
                renderWindow->AddRenderer(renderer);
                renderWindow->SetSize(300,300);

                vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
                    vtkSmartPointer<vtkRenderWindowInteractor>::New();
                renderWindowInteractor->SetRenderWindow(renderWindow);

                // run the window
                renderWindow->Render();
                renderWindowInteractor->Start();
            }

        // debug: same render routine as above, but with a function signature so that it
        // can be called from a foreach function
        void render(const diy::Master::ProxyWithLink& cp, void*)
            {
                render();
            }

    float                *vel[3];            // pointers to vx, vy, vz arrays (v[0], v[1], v[2])
    size_t               nvecs;              // number of velocity vectors
    int                  init, done;         // initial and done flags
    vtkNew<vtkPoints>    points;             // points to be traced
    vtkNew<vtkPolyData>  all_polydata;       // finished streamlines
    vtkNew<vtkCellArray> all_cells;          // all streamlines as vtk cells
    vtkNew<vtkPoints>    all_points;         // all points in all streamlines

    vector<Segment> segments;                // finished segments of particle traces
};
