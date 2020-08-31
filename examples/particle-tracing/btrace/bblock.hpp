#ifndef BBLOCK_HPP
#define BBLOCK_HPP

#include "bmisc.h"
#include <diy/master.hpp>
#include "zoltan.h"
#include "message.h"
#include <diy/point.hpp>

struct BPt
{
    diy::Point<float, 3>    coords;
};

struct BSegment{
    int pid;        // 
    std::vector<BPt> pts; // points along trace
    int gid;       

    BSegment(BEndPt &cur_pt); 
    BSegment();
};

struct MESH_DATA{
  int numGlobalPoints;
  int numMyPoints;
  ZOLTAN_ID_PTR myGlobalIDs;
  float *x;
  float *y;
  float *z;
} ;

struct BBlock
{
  static void*    create()                                    { return new BBlock; }
  static void     destroy(void* b)                            { delete static_cast<BBlock*>(b); }
  static void     save(const void* b, diy::BinaryBuffer& bb)  { diy::save(bb, *static_cast<const BBlock*>(b)); }
  static void     load(void* b, diy::BinaryBuffer& bb)        { diy::load(bb, *static_cast<BBlock*>(b)); }
  // other functions and data members

   


  // store the data here

  std::map<int, std::vector<float>> data;
  std::map<int, std::vector<float>> data_ghost;
  std::vector<data_req> data_reqsts;
//   std::vector<int> bid_to_rank(C *C *C, 0);
  std::vector<int> cid_to_rank;
//   std::vector<int> weights(C *C *C);
  std::vector<int> weights;
  int cside[3]; // length of each cell side
  std::vector<int> partn; // map of gid to current partn
  MESH_DATA mesh_data;

  // for advection
  std::map<int, std::vector<BEndPt>> particles;
  std::map<int, std::vector<BEndPt>> unfinished_local;
  std::map<int, std::vector<BEndPt>> unfinished_nonlocal;

  std::map<int, bbounds> bounds;
  std::map<int, bbounds> bounds_ghost;

  std::vector<BSegment> segments;

};



#endif