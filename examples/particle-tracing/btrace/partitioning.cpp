#include "partitioning.h"
#include "../misc.h"

void partition(bbounds &dom, int C, std::map<int, std::vector<float>> data, int rank, int worldsize, std::vector<int> &bid_to_rank){
    
   
    // int S =  C * cellsperblockside // cells per side;
    // int side[3] = {(dom.max[0] - dom.min[0]+1)/C, (dom.max[0] - dom.min[0]+1)/C, (dom.max[0] - dom.min[0]+1)/C}; 
    // int num_per_side[3] =  {(dom.max[0] - dom.min[0]+1)/side[0], (dom.max[0] - dom.min[0]+1)/side[1], (dom.max[0] - dom.min[0]+1)/side[2]}; 
    // int corner_cell_id = rank / (C*C) + rank / C + rank % C;
    // int corner_cell_id =  rank / (num_per_side[0]) +  rank / num_per_side[1] +  rank % num_per_side[2];

    int rpd = worldsize / C; //{worldsize/C, worldsize/C, worldsize/C}; // ranks per dimension

    int c2 = rank / (rpd*rpd); // slowest varying
    int c1 = (rank % (rpd*rpd)) / rpd;
    int c0 = (rank % (rpd*rpd)) % rpd; // fastest varying
    
    
    int side[3] = {(dom.max[0] - dom.min[0]+1)/rpd, (dom.max[0] - dom.min[0]+1)/rpd, (dom.max[0] - dom.min[0]+1)/rpd }; // cells per side of a rank
    int coords[3] = {c2, c1, c0};

    int ccoords[3] =  {coords[0]*side[0], coords[1]*side[1], coords[2]*side[2]};// side of block

    // dprint("rank %d, coords %d %d %d", rank, coords[0], coords[1], coords[2]);
    // dprint("rank %d, ccoords %d %d %d | rpd %d", rank, ccoords[0],  ccoords[1], ccoords[2], rpd);

    


    int bside[3] = {(dom.max[0] - dom.min[0]+1)/C, (dom.max[0] - dom.min[0]+1)/C, (dom.max[0] - dom.min[0]+1)/C }; // cells per side of a block
    

    // populate the block ids into the map

    // int bpd = C; // blocks per dimension
    // int T = C*C*C; // 

    // c2 = T / (bpd*bpd); // slowest varying
    // c1 = (T % (bpd*bpd)) / bpd;
    // c0 = (T % (bpd*bpd)) % bpd; // fastest varying

    int bpr = C/rpd;

    dprint("corner block ids rank %d | %d %d %d", rank, coords[0]*bpr, coords[1]*bpr, coords[2]*bpr);
    // dprint("rank block corner cells %d | bid %d %d %d", rank, coords[0]*bpr*side[0], coords[1]*bpr*side[0], coords[2]*bpr*side[0]);

    if (rank==0)
        dprint("sides %d %d %d | bpr %d", side[0], side[1], side[2], bpr);



    // add the cell ids to the map and set current rank in bid_to_rank
    
    for (int i=coords[0]*bpr; i<coords[0]*bpr+bpr; i++){
        for (int j=coords[1]*bpr; j< coords[1]*bpr+bpr; j++){
            for (int k=coords[2]*bpr; k < coords[2]*bpr+bpr; k++){
                int bid = i*C*C + j*C + k;
                // dprint ("rank %d, bid %d", rank, bid);

                std::vector<float> temp;
                std::pair<int,std::vector<float>> tpair (bid,temp);
                data.insert(tpair);
                bid_to_rank[bid] = rank;

            }
        }
    }
    dprint("data size %d, %ld", data.size(), bid_to_rank.size());
}
    // dprint("map size %ld", data.size());
// https://stackoverflow.com/questions/20834838/using-tuple-in-unordered-map