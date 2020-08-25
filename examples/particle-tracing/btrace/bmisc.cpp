#include "bmisc.h"
#include <diy/mpi.hpp>


void print_cellids_in_block(std::map<int, std::vector<float>> &data){

    std::map<int, std::vector<float>>::iterator it = data.begin();
    fprintf(stderr, "cellids :");
    while (it != data.end()){
        
        fprintf(stderr, "%d ", it->first );
        it++;
    }
    fprintf(stderr, "\n");


}

void update_weights(diy::mpi::communicator &world, std::map<int, std::vector<BEndPt>> &particles, std::vector<int> &weights){

    // reset weights
    std::fill(weights.begin(), weights.end(), 0);
    std::vector<int> tmp(weights.size());

    // populate weights
    for (auto const& plist : particles)
    {
        // std::cout << p.first << ' ' << p.second << '\n';
        tmp[plist.first] = plist.second.size();
    }

    // do all reduce
    diy::mpi::all_reduce (world, tmp, weights, std::plus<int>());

}