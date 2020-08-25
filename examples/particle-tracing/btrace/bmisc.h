#ifndef BMISC_HPP
#define BMISC_HPP

#include <map>
#include <vector>
#include <diy/mpi.hpp>

struct bbounds{

    // std::vector<int> max;
    // std::vector<int> min;
    // std::vector<int> cside; // elements per cell side along each dimension

    // bbounds(){
    //     max.resize(3);
    //     min.resize(3);
    //     cside.resize(3);
    // }

    int max[3];
    int min[3];
    int cside[3];

};

// namespace diy
// {
//     template<>
//     struct Serialization<bbounds>
//     {
//         static
//         void save(diy::BinaryBuffer& bb, const bbounds& x)
//             {
//                 diy::Serialization<int>::save(bb, x.max, 3);
//                 diy::Serialization<int>::save(bb, x.min, 3);
//                 diy::Serialization<int>::save(bb, x.cside, 3);
                
//             }
//         static
//         void load(diy::BinaryBuffer& bb, bbounds& x)
//             {
//                diy::Serialization<int>::load(bb, x.max, 3);
//                diy::Serialization<int>::load(bb, x.min, 3);
//                diy::Serialization<int>::load(bb, x.cside, 3); 
//             }
//     };
// }

struct BEndPt{

    int pid;
    int cid;
    int nsteps;
    // std::vector<float> pt(3);
    float pt[3];

    const float& operator [](int i) const { return pt[i]; }
    float& operator [](int i)             { return pt[i]; }

    BEndPt(){
        pid=-1;
        cid=-1;
        nsteps=0;
    }

};

void print_cellids_in_block(std::map<int, std::vector<float>> &data);

void update_weights(diy::mpi::communicator &world, std::map<int, std::vector<BEndPt>> &particles, std::vector<int> &weights);

#endif