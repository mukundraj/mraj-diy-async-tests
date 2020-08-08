#include "advection.h"
#include "../misc.h"
#include <cmath>

int pos2cgid(float px, float py, float pz, bbounds &dom, int C){

    // int gid = C*C*(px/dom.cside[0]) + C*(py/dom.cside[1]) + pz/dom.cside[2]; 
    int gid = floor((px/(float)dom.cside[0])) + C*floor((py/(float)dom.cside[1])) + C*C*floor((pz/(float)dom.cside[2])); 

    return gid;
}

void gid2bounds(int gid, int *cside, int C, bbounds &bnd){


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

}

bool inside(float *pt, bbounds &bnd){

      for (int i = 0; i < 3; i++)
        if (pt[i] < (float)(bnd.min[i]) || pt[i] >= (float)(bnd.max[i] + bnd.cside[i] - 1))
            return false;
    return true;
}

void seed(BBlock *b, bbounds &dom, int C, float sr, int rank){


  

    for (float k=dom.min[2]; k<dom.max[2]; k+=sr){
        for (float j=dom.min[1]; j<dom.max[1]; j+=sr){
            for (float i=dom.min[0]; i<dom.max[0]; i+=sr){

                // iterater over each cell in block and insert if seed in any cell

                    
                    float pt[3] = {i,j,k};
                    int gid = pos2cgid(pt[0], pt[1], pt[2], dom, C); 

                    if (rank == 0){
                            dprint(" in rank 0, gid %d, (%f %f %f), dom [%f %f %f]", gid, i, j, k, i/(float)dom.cside[0], j/(float)dom.cside[1], k/(float)dom.cside[2]);
                      }

                    if (b->data.find(gid)!=b->data.end()){
                        BEndPt endp;

                      

                        b->particles[gid].push_back(endp);
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