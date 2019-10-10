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

#ifdef WITH_VTK
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
#endif

#include <fstream>
#include <stdio.h>

#include <pnetcdf.h>
#include <iomanip>      // std::setprecision

typedef diy::DiscreteBounds            Bounds;
typedef diy::RegularGridLink           RGLink;
typedef diy::RegularDecomposer<Bounds> Decomposer;

using namespace std;

static void handle_error(int status, int lineno)
{
    fprintf(stderr, "Error at line %d: %s\n", lineno, ncmpi_strerror(status));
    MPI_Abort(MPI_COMM_WORLD, 1);
}

// the diy block
struct Block
{
    Block() : nvecs(0), init(0), done(0) {}
    ~Block()
    {
        if (nvecs)
        {
            for (int i = 0; i < 3; i++)
                delete[] vel[i];
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
    void show_geometry(const diy::Master::ProxyWithLink& cp)
    {
        diy::RegularLink<Bounds> *link = static_cast<diy::RegularLink<Bounds>*>(cp.link());

        fprintf(stderr, "rank=%d, gid=%d, core={%d, %d, %d} x {%d, %d, %d}, bounds={%d, %d, %d} x {%d, %d, %d}, %d\n",
                cp.master()->communicator().rank(),
                cp.gid(),
                link->core().min[0], link->core().min[1], link->core().min[2],
                link->core().max[0], link->core().max[1], link->core().max[2],
                link->bounds().min[0], link->bounds().min[1], link->bounds().min[2],
                link->bounds().max[0], link->bounds().max[1], link->bounds().max[2],
                link->size());

        fmt::print(stderr, "rank {} gid {} has {} segments\n", cp.master()->communicator().rank(), cp.gid(), segments.size());

        for (size_t i = 0; i < segments.size(); i++)
        {
            fprintf(stderr, "[pid %d  gid %d num_pts %ld]: ", segments[i].pid, segments[i].gid, segments[i].pts.size());
            if (segments[i].pts.size())     // print only first and last point in segment (all points are too much)
                fprintf(stderr, "[%.3f %.3f %.3f] ... [%.3f %.3f %.3f]\n ",segments[i].pts.front().coords[0], segments[i].pts.front().coords[1], segments[i].pts.front().coords[2],
                        segments[i].pts.back().coords[0], segments[i].pts.back().coords[1], segments[i].pts.back().coords[2]);
        }
    }

    void write_segments(std::string filename)
    {
        ofstream f;
        f.open(filename);

        // debug
//         fmt::print("writing {} segments\n", segments.size());

        for (size_t i = 0; i < segments.size(); i++)
        {
            for (size_t j = 0; j < segments[i].pts.size(); j++)
            {
                // debug
//                 fprintf(f, "%ld %f %f %f, ", segments[i].pts.size(), segments[i].pts[0].coords[0], segments[i].pts[0].coords[1], segments[i].pts[0].coords[2]);
                f << std::setprecision(8) << segments[i].pts[j].coords[0] << " " << segments[i].pts[j].coords[1] << " " << segments[i].pts[j].coords[2] << " ";
            }
            f << endl;
        }
        f.close();
    }

#ifdef WITH_VTK
    // convert vector of particle trajectories to vtk polylines and render them
    // TODO: use vtk data model from the outset?
    // copied from http://www.vtk.org/Wiki/VTK/Examples/Cxx/GeometricObjects/PolyLine
    void render()
    {
        // vtk points
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        // vtk cells
        vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
        size_t n = 0;           // index of current point in all points in all traces

        // add points from each trace to one global list of vtkPoints
        for (size_t i = 0; i < segments.size(); i++)
            for (size_t j = 0; j < segments[i].pts.size(); j++)
                points->InsertNextPoint(&(segments[i].pts[j].coords[0])); // deep copy, I assume?

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
    void render_block(const diy::Master::ProxyWithLink& cp)
    {
        render();
    }

#endif

    float                *vel[3];            // pointers to vx, vy, vz arrays (v[0], v[1], v[2])
    size_t               nvecs;              // number of velocity vectors
    int                  init, done;         // initial and done flags
    vector<Segment>      segments;           // finished segments of particle traces
    vector<EndPt>        particles;

#ifdef WITH_VTK
    vtkNew<vtkPoints>    points;             // points to be traced
    vtkNew<vtkPolyData>  all_polydata;       // finished streamlines
    vtkNew<vtkCellArray> all_cells;          // all streamlines as vtk cells
    vtkNew<vtkPoints>    all_points;         // all points in all streamlines
#endif

};

// add a block to the master
struct AddBlock
{
    AddBlock(diy::Master &master_) :
        master(master_)
    {}

    Block* operator()(int gid,
                      const Bounds& core,
                      const Bounds& bounds,
                      const Bounds& domain,
                      const RGLink& link) const
    {
        Block *b       = new Block;
        RGLink *l      = new RGLink(link);
        diy::Master &m = const_cast<diy::Master&>(master);
        m.add(gid, b, l);
        return b;
    }

    diy::Master& master;
};

// add a block to the master and read input data
struct AddAndRead : public AddBlock
{
    AddAndRead(diy::Master& m,
               const char*  infile_,
               diy::mpi::communicator& world_,
               const float vec_scale_,
               const int hdr_bytes_) :
        AddBlock(m),
        infile(infile_),
        world(world_),
        vec_scale(vec_scale_),
        hdr_bytes(hdr_bytes_) {}

    void operator()(int gid,
                    const Bounds& core,
                    const Bounds& bounds,
                    const Bounds& domain,
                    const RGLink& link) const
    {
        Block* b = AddBlock::operator()(gid, core, bounds, domain, link);
        MPI_Offset *start, *count;
        float *data_u=NULL, *data_v=NULL, *data_w=NULL;

        int ncfile, ndims, nvars, ngatts, unlimited;
        int ret;
        ret = ncmpi_open(world, infile, NC_NOWRITE, MPI_INFO_NULL,&ncfile);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        ret = ncmpi_inq(ncfile, &ndims, &nvars, &ngatts, &unlimited);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);


        // reversed order of shape and bounds needed because the sample data file
        // is linearized in row-major (C) order
        vector<int> shape(3);
        for (size_t i = 0; i < 3; i++)
            shape[2 - i] = domain.max[i] - domain.min[i] + 1;
        //        diy::io::BOV reader(in, shape, hdr_bytes);

        Bounds r_bounds { 3 };
        r_bounds.min[0] = bounds.min[2];
        r_bounds.max[0] = bounds.max[2];
        r_bounds.min[1] = bounds.min[1];
        r_bounds.max[1] = bounds.max[1];
        r_bounds.min[2] = bounds.min[0];
        r_bounds.max[2] = bounds.max[0];



        start = (MPI_Offset*) calloc(ndims, sizeof(MPI_Offset));
        count = (MPI_Offset*) calloc(ndims, sizeof(MPI_Offset));

        if (ndims==4){
            count[0] = 1;
            count[1] = r_bounds.max[0] - r_bounds.min[0]+1;
            count[2] = r_bounds.max[1] - r_bounds.min[1]+1;
            count[3] = r_bounds.max[2] - r_bounds.min[2]+1;

            start[0] =  0; start[1] = r_bounds.min[0]; start[2] = r_bounds.min[1]; start[3] = r_bounds.min[2];
        }else if(ndims==3){

            count[0] = r_bounds.max[0] - r_bounds.min[0]+1;
            count[1] = r_bounds.max[1] - r_bounds.min[1]+1;
            count[2] = r_bounds.max[2] - r_bounds.min[2]+1;

            start[0] = r_bounds.min[0]; start[1] = r_bounds.min[1]; start[2] = r_bounds.min[2];
        }

        //        std::cout<<"counts"<<count[0]<<" "<<count[1]<<" "<<count[2]<<"\n";
        //        std::cout<<"starts"<<start[0]<<" "<<start[1]<<" "<<start[2]<<"\n";

        size_t nvecs =
                (bounds.max[0] - bounds.min[0] + 1) *
                (bounds.max[1] - bounds.min[1] + 1) *
                (bounds.max[2] - bounds.min[2] + 1);
        vector<float> values(nvecs * 3); // temporary contiguous buffer of input vector values

        data_u = (float*) calloc(nvecs, sizeof(float));
        data_v = (float*) calloc(nvecs, sizeof(float));
        data_w = (float*) calloc(nvecs, sizeof(float));
        ret = ncmpi_get_vara_float_all(ncfile, 0, start, count, data_u);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);
        ret = ncmpi_get_vara_float_all(ncfile, 1, start, count, data_v);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);
        ret = ncmpi_get_vara_float_all(ncfile, 2, start, count, data_w);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);



        // copy from temp values into block
        b->vel[0] = new float[nvecs];
        b->vel[1] = new float[nvecs];
        b->vel[2] = new float[nvecs];
        b->nvecs = nvecs;
        for (size_t i = 0; i < nvecs; i++)
        {

            b->vel[0][i] = data_u[i] * vec_scale;
            b->vel[1][i] = data_v[i] * vec_scale;
            b->vel[2][i] = data_w[i] * vec_scale;

        }

        ret = ncmpi_close(ncfile);
        free(start);
        free(count);
        free(data_u);
        free(data_v);
        free(data_w);
    }


    const char*	infile;
    diy::mpi::communicator world;
    float vec_scale;
    int hdr_bytes;
};

// convert linear domain point index into (i,j,k,...) multidimensional index
// number of dimensions is the domain dimensionality
void idx2ijk(
        size_t                  idx,                // linear cell indx
        const vector<size_t>&   ds,                 // stride of input points
        const Bounds&           bounds,             // block bounds
        vector<size_t>&         ijk)                // i,j,k,... indices in all dimensions
{
    int dim = ds.size();
    if (dim == 1)
    {
        ijk[0] = idx;
        return;
    }

    for (auto i = 0; i < dim; i++)
    {
        if (i < dim - 1)
            ijk[i] = bounds.min[i] + (idx % ds[i + 1]) / ds[i];
        else
            ijk[i] = bounds.min[i] + idx / ds[i];
    }
}

// add a block to the master and set synthetic vector field
// block along diagonal of block coords is slower than the rest
struct AddSynthetic1 : public AddBlock
{
    AddSynthetic1(diy::Master&           m,
                 const float             slow_vel_,         // slow velocity
                 const float             fast_vel_,         // fast velocity
                 const Decomposer&       decomposer_) :
        AddBlock(m),
        slow_vel(slow_vel_),
        fast_vel(fast_vel_),
        decomposer(decomposer_) {}

    void operator()(int gid,
                    const Bounds& core,
                    const Bounds& bounds,
                    const Bounds& domain,
                    const RGLink& link) const
    {
        Block* b = AddBlock::operator()(gid, core, bounds, domain, link);

        b->nvecs =                                  // total number of vectors in the block
                (bounds.max[0] - bounds.min[0] + 1) *
                (bounds.max[1] - bounds.min[1] + 1) *
                (bounds.max[2] - bounds.min[2] + 1);

        b->vel[0] = new float[b->nvecs];
        b->vel[1] = new float[b->nvecs];
        b->vel[2] = new float[b->nvecs];

        // set synthetic velocity vectors
        std::vector<int> coords;                            // coordinates of block in each dimension
        decomposer.gid_to_coords(gid, coords);
        vector<int> divs(coords.size());                    // number of blocks in each dimension
        decomposer.fill_divisions(divs);

        // debug
//         if (gid == 0)
//             fmt::print(stderr, "divs = [{} {} {}]\n", divs[0], divs[1], divs[2]);

        // one slow diy block along diagonal
        bool slow_block = true;
        for (int i = 0; i < coords.size(); i++)
        {
            if (divs[i] == 1)
                continue;
            if (i > 0 && coords[i] != coords[i - 1])
            {
                slow_block =  false;
                break;
            }
        }

        // debug
//         if (slow_block)
//             fmt::print(stderr, "gid {} is slow\n", gid);

        // set the velocity
        for (size_t i = 0; i < b->nvecs; i++)
        {
            if (slow_block)
                b->vel[0][i] = slow_vel;
            else
                b->vel[0][i] = fast_vel;
            b->vel[1][i] = 0.0;
            b->vel[2][i] = 0.0;
        }
    }

    Decomposer  decomposer;
    float       slow_vel, fast_vel;
};

// add a block to the master and set synthetic vector field
// velocities vary along a gradient in y direction, with oppsite gradients in blocks adjacent in x direction
struct AddSynthetic2 : public AddBlock
{
    AddSynthetic2(diy::Master&           m,
                 const float             slow_vel_,
                 const float             fast_vel_,
                 const Decomposer&       decomposer_) :
        AddBlock(m),
        slow_vel(slow_vel_),
        fast_vel(fast_vel_),
        decomposer(decomposer_) {}

    void operator()(int gid,
                    const Bounds& core,
                    const Bounds& bounds,
                    const Bounds& domain,
                    const RGLink& link) const
    {
        Block* b = AddBlock::operator()(gid, core, bounds, domain, link);

        b->nvecs =                                  // total number of vectors in the block
                (bounds.max[0] - bounds.min[0] + 1) *
                (bounds.max[1] - bounds.min[1] + 1) *
                (bounds.max[2] - bounds.min[2] + 1);

        b->vel[0] = new float[b->nvecs];
        b->vel[1] = new float[b->nvecs];
        b->vel[2] = new float[b->nvecs];

        // set synthetic velocity vectors
        // blocks in adjacent x-columns have oppposite direction velocity gradients in y direction
        std::vector<int> coords;                            // coordinates of block in each dimension
        decomposer.gid_to_coords(gid, coords);

        for (size_t i = 0; i < b->nvecs; i++)
        {
            // we want to create a gradient in the x-velocity based on the y-z coordinates of the vector
            // velocity is constant across the x direction
            size_t nx = bounds.max[0] - bounds.min[0] + 1;          // number of vectors in x direction
            size_t j = i / nx;                                      // index of vector in y-z plane, ignoring x coordinate
            if (coords[0] % 2 == 0)
                b->vel[0][i] = slow_vel + (fast_vel - slow_vel) * (float)j / b->nvecs * nx;
            else
                b->vel[0][i] = fast_vel - (fast_vel - slow_vel) * (float)j / b->nvecs * nx;
            b->vel[1][i] = 0.0;
            b->vel[2][i] = 0.0;
        }
    }

    Decomposer  decomposer;
    float       slow_vel, fast_vel;
};


