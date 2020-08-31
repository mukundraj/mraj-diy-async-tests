#include "bblock.hpp"


BSegment::BSegment(BEndPt &cur_pt){

    
    pid = cur_pt.pid;
    gid = cur_pt.cid;

    BPt pt;
    pt.coords[0] = cur_pt[0];
    pt.coords[1] = cur_pt[1];
    pt.coords[2] = cur_pt[2];
    pts.push_back(pt);
}

BSegment::BSegment(){

}