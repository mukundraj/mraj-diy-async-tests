#ifndef BTRACE_HPP
#define BTRACE_HPP

#include <vector>
#include "bmisc.h"


struct datablock{
    int to_proc;
    int from_proc;
    std::vector<int> cgid; // cell gid
    std::vector<std::vector<float>> data;
    // std::vector<BEndPt>        particles;
};


namespace diy
{
    template<>
    struct Serialization<datablock>
    {
        static
        void save(diy::BinaryBuffer& bb, const datablock& x)
            {
                diy::Serialization<int>::save(bb, x.to_proc);
                diy::Serialization<int>::save(bb, x.from_proc);
                diy::Serialization<std::vector<int>>::save(bb, x.cgid);
                diy::Serialization<std::vector<std::vector<float>>>::save(bb, x.data);
            }
        static
        void load(diy::BinaryBuffer& bb, datablock& x)
            {
                diy::Serialization<int>::load(bb, x.to_proc);
                diy::Serialization<int>::load(bb, x.from_proc);
                diy::Serialization<std::vector<int>>::load(bb, x.cgid);
                diy::Serialization<std::vector<std::vector<float>>>::load(bb, x.data);
            }
    };
}



#endif