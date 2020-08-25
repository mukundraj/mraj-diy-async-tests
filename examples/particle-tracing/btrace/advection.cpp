#include "advection.h"
#include "../misc.h"
#include <cmath>

inline float texel3d(std::vector<float> &vec, bbounds &bnd, int x, int y, int z, int v){
    size_t len = vec.size()/3;
    return vec[v*len + x + bnd.cside[0] * (y + bnd.cside[1] * z)];
}

// inline float texel3d(std::vector<float> &vec, bbounds &bnd, int x, int y, int z, int v){
//     size_t len = vec.size()/3;
//     return vec[v*len + x + (bnd.cside[0]+1) * (y + (bnd.cside[1]+1) * z)];
// }

bool badvect_rk1(BEndPt &pt, BBlock *b, bbounds &dom, int C, float h, BEndPt &next_p){

     int cid = in_block(&pt[0], b->data, b->data_ghost, dom, C);


     
    //  print_cellids_in_block(b->data);

    float p[3] = {pt[0], pt[1], pt[2]};
    float v[3];
     
     if (b->data.find(cid)!=b->data.end()){
         lerp(&p[0], b->bounds[cid], b->data[cid], &v[0]);
     }else if (b->data_ghost.find(cid) != b->data_ghost.end()){
        //  dprint("d->data_ghost[cid] size %ld", b->data_ghost[cid].size());
         lerp(&p[0], b->bounds_ghost[cid], b->data_ghost[cid], &v[0]);

     }
     else{
        next_p.pid = pt.pid;
        next_p[0] = pt[0];
        next_p[1] = pt[1];
        next_p[2] = pt[2];
        next_p.cid = pt.cid;
        next_p.nsteps = pt.nsteps;
        // dprint("particle's cell not in block. exiting.");
        return false;
     }

    
     next_p.pid = pt.pid;
     next_p[0] = p[0] + h*v[0];   
     next_p[1] = p[1] + h*v[1];   
     next_p[2] = p[2] + h*v[2]; 
     next_p.nsteps ++;
     int ccid = pos2cgid(next_p[0], next_p[1], next_p[2], dom, C);
     next_p.cid = ccid;
    //  dprint("advected cblok %d, nsteps %d, next_p(%f %f %f), ccid %d, pid %d", cid, next_p.nsteps, next_p[0], next_p[1], next_p[2], ccid, next_p.pid);


    return true;

}

void lerp(float *pt, bbounds &bnd,  std::vector<float> &vec, float *vars){

    float p[8];

    float x  = pt[0] - (float)(bnd.min[0]),       // physical coords of point relative to min. of block
          y  = pt[1] - (float)(bnd.min[1]),
          z  = pt[2] - (float)(bnd.min[2]);
    int   i  = floor(x),                     // indices of min corner of texel containing point
          j  = floor(y),
          k  = floor(z);
    int   i1 =i + 1,                         // indices of max corner of texel containing point
          j1 =j + 1,
          k1 =k + 1;
    float x0 = i, x1 = i1,                   // float version of min, max corner of texel
          y0 = j, y1 = j1,
          z0 = k, z1 = k1;
    int v;                                   // dimension
    // dprint("xyz %f %f %f, ijk [%d %d %d] [%d %d %d]", x, y, z, i, j, k, i1, j1, k1);
    for (v = 0; v < 3; v++)
    {
        p[0] = texel3d(vec, bnd, i, j, k, v);
        p[1] = texel3d(vec, bnd, i1, j, k, v);
        p[2] = texel3d(vec, bnd, i, j1, k, v);
        p[3] = texel3d(vec, bnd, i1, j1, k, v);
        p[4] = texel3d(vec, bnd, i, j, k1, v);
        p[5] = texel3d(vec, bnd, i1, j, k1, v);
        p[6] = texel3d(vec, bnd, i, j1, k1, v);
        p[7] = texel3d(vec, bnd, i1, j1, k1, v);
        // dprint("ps %f %f %f %f %f %f %f %f, vecsize %ld", p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], vec.size());
        
        vars[v] =
            p[0] * (x1 - x) * (y1 - y) * (z1 - z) +
            p[1] * (x - x0) * (y1 - y) * (z1 - z) +
            p[2] * (x1 - x) * (y - y0) * (z1 - z) +
            p[3] * (x - x0) * (y - y0) * (z1 - z) +
            p[4] * (x1 - x) * (y1 - y) * (z - z0) +
            p[5] * (x - x0) * (y1 - y) * (z - z0) +
            p[6] * (x1 - x) * (y - y0) * (z - z0) +
            p[7] * (x - x0) * (y - y0) * (z - z0);

    }


}

int in_block(float *pt, std::map<int, std::vector<float>> &data, std::map<int, std::vector<float>> &data_ghost, bbounds &dom, int C){

    int gid = pos2cgid(pt[0], pt[1], pt[2], dom, C);

    if (data.find(gid) == data.end() && data_ghost.find(gid) == data_ghost.end() )
        return -1;

    return gid;
}

int pos2cgid(float px, float py, float pz, bbounds &dom, int C){

    // int gid = C*C*(px/dom.cside[0]) + C*(py/dom.cside[1]) + pz/dom.cside[2]; 
    int gid = floor((px/(float)dom.cside[0])) + C*floor((py/(float)dom.cside[1])) + C*C*floor((pz/(float)dom.cside[2])); 
    // dprint("p %f %f %f, cid %d", px, py, pz, gid);

    return gid;
}

void gid2bounds(int gid, int *cside, int C, bbounds &dom, bbounds &bnd){


    int c2 = gid / (C*C); // slowest varying
    int c1 = (gid % (C*C)) / C;
    int c0 = (gid % (C*C)) % C; // fastest varying

    bnd.min[0] = c0*cside[0];
    bnd.max[0] = c0*cside[0] + cside[0];
    bnd.min[1] = c1*cside[1];
    bnd.max[1] = c1*cside[1] + cside[1];;
    bnd.min[2] = c2*cside[2];
    bnd.max[2] = c2*cside[2] + cside[2];

    bnd.cside[0] = cside[0];
    bnd.cside[1] = cside[1];
    bnd.cside[2] = cside[2];

    if (bnd.max[0] < dom.max[0]){
        bnd.max[0] += 1;
        // bnd.cside[0] += 1;
    }

    if (bnd.max[1] < dom.max[1]){
        bnd.max[1] += 1;
        // bnd.cside[1] += 1;
    }

    if (bnd.max[2] < dom.max[2]){
        bnd.max[2] += 1;
        // bnd.cside[2] += 1;
    }
        
}

bool inside(float *pt, bbounds &bnd){

      for (int i = 0; i < 3; i++)
        if (pt[i] < (float)(bnd.min[i]) || pt[i] >= (float)(bnd.max[i] + bnd.cside[i] - 1))
            return false;
    return true;
}

void seed(BBlock *b, bbounds &dom, int C, float sr, int rank){


    
    int count = 0;
    for (float k=dom.min[2]; k<dom.max[2]; k+=sr){
        for (float j=dom.min[1]; j<dom.max[1]; j+=sr){
            for (float i=dom.min[0]; i<dom.max[0]; i+=sr){

                // iterater over each cell in block and insert if seed in any cell

                    
                    float pt[3] = {i,j,k};
                    int gid = pos2cgid(pt[0], pt[1], pt[2], dom, C); 

                    // if (rank == 0){
                    //         dprint(" in rank 0, gid %d, (%f %f %f), dom [%f %f %f]", gid, i, j, k, i/(float)dom.cside[0], j/(float)dom.cside[1], k/(float)dom.cside[2]);
                    //   }

                    if (b->data.find(gid)!=b->data.end()){
                        BEndPt endp;
                        endp[0] = i;
                        endp[1] = j;
                        endp[2] = k;
                        endp.pid = (rank+1) * 100 + count;
                        endp.cid = gid;

                        // dprint("rank %d, seed %d| %f %f %f", rank, endp.cid, i,j,k);

                        b->particles[gid].push_back(endp);
                        count ++;
                    }

                    

                    // std::map<int, std::vector<float>>::iterator it = b->data.begin();
                    
                    // while (it != b->data.end()){

                        
                    //     float pt[3] = {i,j,k};
                    //     int gid = pos2cgid(pt[0], pt[1], pt[2], dom, C);
                    //     bbounds bnd;
                    //     gid2bounds(gid, &dom.cside[0], C, bnd);
                    //     if (inside(&pt[0], bnd)){
                    //         BEndPt endpt;
                    //         b->particles.push_back(endpt);
                    //         break;
                    //     }
                    //     it++;
                    // }


            }
        }
    }

    
}