#include "io.h"
#include <pnetcdf.h>
#include "../misc.h"

static void handle_error(int status, int lineno)
{
    fprintf(stderr, "Error at line %d: %s\n", lineno, ncmpi_strerror(status));
    MPI_Abort(MPI_COMM_WORLD, 1);
}

void read_data(diy::mpi::communicator& world, const char*  infile, std::map<int, std::vector<float>> &data, std::vector<int> &weights, int* C, int *cside, BBlock *b, bbounds &dom){

    for(size_t i=0; i<weights.size(); i++){
        weights[i] = 1;
    }


    MPI_Offset *start, *count;
    float *data_u=NULL, *data_v=NULL, *data_w=NULL;

    int ncfile, ndims, nvars, ngatts, unlimited;
    int ret;
    ret = ncmpi_open(world, infile, NC_NOWRITE, MPI_INFO_NULL,&ncfile);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    ret = ncmpi_inq(ncfile, &ndims, &nvars, &ngatts, &unlimited);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    std::map<int, std::vector<float>>::iterator it = data.begin();




    while (it != data.end()){

        int i = it->first;
        int ix = i%C[0];
        int iy = i%(C[0]*C[1])/C[0];
        int iz = i/(C[0]*C[1]);

        
       

        MPI_Offset count[3];
        // count[0] = bside[2]; // xa
        // count[1] = bside[1];  // y
        // count[2] = bside[0]; // z (slowest changing)
        MPI_Offset start [3];
        // start[0] = ix*bside[2];
        // start[1] = iy*bside[1];
        // start[2] = iz*bside[0];

        // count[0] = b->bounds[i].cside[2];
        // count[1] = b->bounds[i].cside[1];
        // count[2] = b->bounds[i].cside[0];

        // count[0] = b->bounds[i].max[2] - b->bounds[i].min[2];
        // count[1] = b->bounds[i].max[1] - b->bounds[i].min[1];
        // count[2] = b->bounds[i].max[0] - b->bounds[i].min[0]; 


        count[0] = b->bounds[i].max[2] - b->bounds[i].min[2];
        count[1] = b->bounds[i].max[1] - b->bounds[i].min[1];
        count[2] = b->bounds[i].max[0] - b->bounds[i].min[0]; 


        // start[0] = ix*dom.cside[2];
        // start[1] = iy*dom.cside[1];
        // start[2] = iz*dom.cside[0];
        start[0] = iz*dom.cside[2];
        start[1] = iy*dom.cside[1];
        start[2] = ix*dom.cside[0];

        if (world.rank()==7)
            dprint("start %d %d %d count %d %d %d, rank %d", start[0], start[1], start[2], count[0], count[1], count[2], world.rank());

        // int len = b->bounds[i].cside[2]*b->bounds[i].cside[1]*b->bounds[i].cside[0];
        int len = count[0]*count[1]*count[2];
        it->second.resize(3*len);

        //  int len = (1+b->bounds[i].cside[2])*(1+b->bounds[i].cside[1])*(1+b->bounds[i].cside[0]);
        // it->second.resize(3*len);
        ret = ncmpi_get_vara_float_all(ncfile, 0, start, count, &it->second[0]); // 
        if (ret != NC_NOERR){
            dprint("err start %d %d %d count %d %d %d, rank %d, cid %d", start[0], start[1], start[2], count[0], count[1], count[2], world.rank(), it->first);
        } //handle_error(ret, __LINE__);
         ret = ncmpi_get_vara_float_all(ncfile, 1, start, count, &it->second[len]); // 
        if (ret != NC_NOERR) handle_error(ret, __LINE__);
         ret = ncmpi_get_vara_float_all(ncfile, 2, start, count, &it->second[2*len]); // 
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        //  dprint("key %d,  (%d %d %d) | %f", i, ix*bside[0],iy*bside[1],iz*bside[2], it->second[0]);

        it++;
    }
    
    ret = ncmpi_close(ncfile);
}

