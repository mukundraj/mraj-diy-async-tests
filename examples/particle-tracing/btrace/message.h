#ifndef BTRACE_HPP
#define BTRACE_HPP

#include <vector>
#include "bmisc.h"

// structure for requesting ghost data
struct data_req{
    int from_proc;
    int to_proc;
    std::vector<int> cids;
};

struct datablock{
    int to_proc;
    int from_proc;
    std::vector<int> cgid; // cell gid
    std::vector<std::vector<float>> data;
    std::vector<std::vector<BEndPt>>        particles;
    std::vector<bbounds> bounds;
};


namespace diy
{
     template<>
    struct Serialization<data_req>
    {
        static
        void save(diy::BinaryBuffer& bb, const data_req& x)
            {
                diy::Serialization<int>::save(bb, x.to_proc);
                diy::Serialization<int>::save(bb, x.from_proc);
                diy::Serialization<std::vector<int>>::save(bb, x.cids);
            }
        static
        void load(diy::BinaryBuffer& bb, data_req& x)
            {
                diy::Serialization<int>::load(bb, x.to_proc);
                diy::Serialization<int>::load(bb, x.from_proc);
                diy::Serialization<std::vector<int>>::load(bb, x.cids);
            }
    };
    
    template<>
    struct Serialization<datablock>
    {
        static
        void save(diy::BinaryBuffer& bb, const datablock& x)
            {
                diy::Serialization<int>::save(bb, x.to_proc);
                diy::Serialization<int>::save(bb, x.from_proc);
                diy::Serialization<std::vector<int>>::save(bb, x.cgid);
                diy::Serialization<std::vector<std::vector<float>>>::save(bb, x.data),
                diy::Serialization<std::vector<std::vector<BEndPt>>>::save(bb, x.particles),
                diy::Serialization<std::vector<bbounds>>::save(bb, x.bounds);
            }
        static
        void load(diy::BinaryBuffer& bb, datablock& x)
            {
                diy::Serialization<int>::load(bb, x.to_proc);
                diy::Serialization<int>::load(bb, x.from_proc);
                diy::Serialization<std::vector<int>>::load(bb, x.cgid);
                diy::Serialization<std::vector<std::vector<float>>>::load(bb, x.data),
                diy::Serialization<std::vector<std::vector<BEndPt>>>::load(bb, x.particles),
                diy::Serialization<std::vector<bbounds>>::load(bb, x.bounds);
            }
    };


    // template<>
    // struct Serialization<BEndPt>
    // {
    //     static
    //     void save(diy::BinaryBuffer& bb, const BEndPt& x)
    //         {
    //            diy::Serialization<int>::save(bb, x.pid); 
    //            diy::Serialization<int>::save(bb, x.gid); 
    //            diy::Serialization<int>::save(bb, x.nsteps); 
    //            diy::Serialization<float*>::save(bb, x.pt); 
    //         }
    //     static
    //     void load(diy::BinaryBuffer& bb, BEndPt& x)
    //         {
    //           diy::Serialization<int>::load(bb, x.pid); 
    //            diy::Serialization<int>::load(bb, x.gid); 
    //            diy::Serialization<int>::load(bb, x.nsteps); 
    //            diy::Serialization<std::vector<float>>::load(bb, x.pt);  
    //         }

    // }
}



#endif