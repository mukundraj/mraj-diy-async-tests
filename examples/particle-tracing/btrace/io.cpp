#include "io.h"
#include <pnetcdf.h>
#include "../misc.h"

static void handle_error(int status, int lineno)
{
    fprintf(stderr, "Error at line %d: %s\n", lineno, ncmpi_strerror(status));
    MPI_Abort(MPI_COMM_WORLD, 1);
}

void read_data(diy::mpi::communicator& world, const char*  infile, std::map<int, std::vector<float>> &data, std::vector<int> &weights, int C, int *bside){

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
        int ix = i%C;
        int iy = i%(C*C)/C;
        int iz = i/(C*C);
       

        MPI_Offset count[3];
        count[0] = bside[2]; // x
        count[1] = bside[1];  // y
        count[2] = bside[0]; // z (slowest changing)
        MPI_Offset start [3];
        start[0] = ix*bside[2];
        start[1] = iy*bside[1];
        start[2] = iz*bside[0];

        ret = ncmpi_get_vara_float_all(ncfile, 0, start, count, &it->second[0]); // 
        if (ret != NC_NOERR) handle_error(ret, __LINE__);
         ret = ncmpi_get_vara_float_all(ncfile, 1, start, count, &it->second[0+bside[1]]); // 
        if (ret != NC_NOERR) handle_error(ret, __LINE__);
         ret = ncmpi_get_vara_float_all(ncfile, 2, start, count, &it->second[0+bside[1]+bside[2]]); // 
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

         dprint("key %d,  (%d %d %d) | %f", i, ix*bside[0],iy*bside[1],iz*bside[2], it->second[0]);

        it++;
    }


    ret = ncmpi_close(ncfile);
}

