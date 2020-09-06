#include "advection.h"
#include "../misc.h"
#include <cmath>

bool in_global_dom(bbounds &dom, BEndPt &p){
    
    bool inside = true;
    for (int i=0; i<3; i++){
        if (p[i] < dom.min[i] || p[i] > dom.max[i])
            inside = false;
    }
    return inside;
}

inline float texel3d(std::vector<float> &vec, bbounds &bnd, int x, int y, int z, int v){
    size_t len = vec.size()/3;
    // dprint("bcs %d %d %d", bnd.cside[0], bnd.cside[1], bnd.cside[2]);
    return vec[v*len + x + bnd.cside[0] * (y + bnd.cside[1] * z)];
}

// inline float texel3d(std::vector<float> &vec, bbounds &bnd, int x, int y, int z, int v){
//     size_t len = vec.size()/3;
//     return vec[v*len + x + (bnd.cside[0]+1) * (y + (bnd.cside[1]+1) * z)];
// }

bool badvect_rk1(BEndPt &pt, BBlock *b, bbounds &dom, int *C, float h, BEndPt &next_p, int rank){

    // pt[0] = 100.3;
    // pt[1] = 100.4;
    // pt[2] = 100.5;

     int cid = in_block(&pt[0], b->data, b->data_ghost, dom, C);


     
    //  print_cellids_in_block(b->data);

    float p[3] = {pt[0], pt[1], pt[2]};
    float v[3];

     

     if (b->data.find(cid)!=b->data.end()){
         lerp(&p[0], b->bounds[cid], b->data[cid], &v[0]);
        //  dprint("in main d->data[cid] size %ld", b->data[cid].size());
     }else if (b->data_ghost.find(cid) != b->data_ghost.end()){
         if (b->data_ghost[cid].size()<1)
         dprint("d->data_ghost[cid] size %ld, cid %d, pid %d, rank %d", b->data_ghost[cid].size(), cid, pt.pid, rank);
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
    //  dprint("vel %f %f %f, cid %d", v[0], v[1], v[2], cid);
    //  exit(0);

    
     next_p.pid = pt.pid;
     next_p[0] = p[0] + h*v[0];   
     next_p[1] = p[1] + h*v[1];   
     next_p[2] = p[2] + h*v[2]; 
     next_p.nsteps = pt.nsteps + 1;
     int ccid = pos2cgid(next_p[0], next_p[1], next_p[2], dom, C);
     next_p.cid = ccid;
    //  if (pt.pid == 401)
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

int in_block(float *pt, std::map<int, std::vector<float>> &data, std::map<int, std::vector<float>> &data_ghost, bbounds &dom, int *C){

    int gid = pos2cgid(pt[0], pt[1], pt[2], dom, C);

    if (data.find(gid) == data.end() && data_ghost.find(gid) == data_ghost.end() )
        return -1;

    return gid;
}

int pos2cgid(float px, float py, float pz, bbounds &dom, int *C){

    // int gid = C*C*(px/dom.cside[0]) + C*(py/dom.cside[1]) + pz/dom.cside[2]; 
    int gid = floor((px/(float)dom.cside[0])) + C[0]*floor((py/(float)dom.cside[1])) + C[0]*C[1]*floor((pz/(float)dom.cside[2])); 
    // if (gid < 0)
    //     dprint("lessthanzero p %f %f %f, cid %d, domcs %d %d %d", px, py, pz, gid, dom.cside[0], dom.cside[1], dom.cside[2]);

    return gid;
}

void gid2bounds(int gid, int *cside, int *C, bbounds &dom, bbounds &bnd){


    int c2 = gid / (C[0]*C[1]); // slowest varying
    int c1 = (gid % (C[0]*C[1])) / C[0];
    int c0 = (gid % (C[0]*C[1])) % C[0]; // fastest varying

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
        bnd.cside[0] += 1;
    }

    if (bnd.max[1] < dom.max[1]){
        bnd.max[1] += 1;
        bnd.cside[1] += 1;
    }

    if (bnd.max[2] < dom.max[2]){
        bnd.max[2] += 1;
        bnd.cside[2] += 1;
    }
        
}

bool inside(float *pt, bbounds &bnd){

      for (int i = 0; i < 3; i++)
        if (pt[i] < (float)(bnd.min[i]) || pt[i] >= (float)(bnd.max[i] + bnd.cside[i] - 1))
            return false;
    return true;
}

int seed(BBlock *b, bbounds &dom, int *C, float sr, int rank){


    
    int count = 0;
    for (float k=dom.min[2]; k<dom.max[2]; k+=sr){
        for (float j=dom.min[1]; j<dom.max[1]; j+=sr){
            for (float i=dom.min[0]; i<dom.max[0]; i+=sr){

                // iterater over each cell in block and insert if seed in any cell

                    
                    float pt[3] = {i,j,k};
                    int gid = pos2cgid(pt[0], pt[1], pt[2], dom, &C[0]); 

                    // if (rank == 0){
                    //         dprint(" in rank 0, gid %d, (%f %f %f), dom [%f %f %f]", gid, i, j, k, i/(float)dom.cside[0], j/(float)dom.cside[1], k/(float)dom.cside[2]);
                    //   }

                    if (b->data.find(gid)!=b->data.end()){
                        BEndPt endp;
                        endp[0] = i;
                        endp[1] = j;
                        endp[2] = k;
                        // endp.pid = (rank+1) * 100 + count;
                        endp.pid = i*512*512 + j*512 + k;
                        endp.cid = gid;
                        // if (rank==0)
                        //     dprint("rank %d, seed %d| %f %f %f | count %d", rank, endp.pid, i,j,k, count);

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

    return count;
}


void enqueue_local_particles(std::map<int, std::vector<BEndPt>> &unfinished_local, 
std::map<int, std::vector<BEndPt>> &particles){

   
    // iterate over cells in block
    for (auto &p : particles){
        
        if (unfinished_local.find(p.first) != unfinished_local.end()){
            p.second = std::move(unfinished_local[p.first]);
            dprint("unfinised exist");
        }else{
            // p.second.clear();
            particles.erase(p.first);
        }

    }


}

inline void remote_enq_nonlocal(BBlock*,
        const diy::Master::ProxyWithLink&   cp,
        const diy::Assigner&                assigner, 
        datablock &db){

    int my_gid              = cp.gid();
    int dest_gid            = db.to_proc;
    int dest_proc           = assigner.rank(dest_gid);
    diy::BlockID dest_block = {dest_gid, dest_proc};

    // db.gid = my_gid;
    // if (my_gid == 4)
    //     dprint("enqing my_gid %d, dest_gid %d, dest_proc %d, size %ld", my_gid, dest_gid, dest_proc, db.cgid.size());

    cp.enqueue(dest_block, db);


}

inline void remote_deq_nonlocal(
       BBlock* b, const diy::Master::ProxyWithLink& cp
){
    std::vector<int> incoming_gids;
    cp.incoming(incoming_gids);
    for (size_t i = 0; i < incoming_gids.size(); i++)
        if (cp.incoming(incoming_gids[i]).size())
        {   
            datablock recvd_data;
            cp.dequeue(incoming_gids[i], recvd_data);
            // if (cp.gid()==0)
            //     dprint("recvd from %d in %d", recvd_data.from_proc, recvd_data.to_proc);
            for(int i=0; i< recvd_data.cgid.size(); i++){
                int cid = recvd_data.cgid[i]; 

                // if (cp.gid()==0)
                    // dprint("recvd cellid %d, fromrank %d", cid, recvd_data.from_proc);
                // dprint("received pid %d, nsteps %d", recvd_data.particles[i][0].pid, recvd_data.particles[i][0].nsteps);
                // enqueue received particles into b->particles[cid]
                b->particles[cid].insert(b->particles[cid].end(), recvd_data.particles[i].begin(), recvd_data.particles[i].end());
            }
        } 
}

void exchange_nonlocal_particles(diy::Master &master, const diy::Assigner& assigner){

    // enqueue nonlocal particles
     master.foreach ([&](BBlock *b, const diy::Master::ProxyWithLink &cp) {
        
           // iterating over each dest rank
           for (auto &plist : b->unfinished_nonlocal){
               datablock db;
               int dest_rank = plist.first;
               
               db.to_proc = dest_rank;
               db.from_proc = cp.gid(); // only for debugging
               
            //    if (db.from_proc==4)
            //         dprint("dest_rank %d, rank %d, unfin_nonloc %ld, to_proc %d", dest_rank , cp.gid(), b->unfinished_nonlocal.size(), db.to_proc);
                for (auto next_p: plist.second){
                    db.cgid.push_back(next_p.cid);
                    db.particles.push_back( std::move(plist.second));
                }

                remote_enq_nonlocal(b, cp, assigner, db);
           }

            
     });



    // rexchange, retrieve and pocess particles
    bool remote = true;
    master.exchange(remote);
    master.foreach(&remote_deq_nonlocal);


}