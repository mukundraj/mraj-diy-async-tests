#ifndef PATITIONING_H
#define PATITIONING_H

#include "bmisc.h"
#include <map>
#include <vector>
#include <map>
#include <diy/mpi.hpp>
#include "zoltan.h"
#include "bblock.hpp"

/* Structure to hold mesh data */


void assign(diy::mpi::communicator &world, std::map<int, std::vector<float>> &data, std::map<int, std::vector<BEndPt>> particles,  std::vector<int> &weights, std::vector<int> &partn, MESH_DATA &myMesh, BBlock* b, const diy::Master::ProxyWithLink &cp, const diy::RoundRobinAssigner &assigner);

void partition(diy::mpi::communicator& world, bbounds &dom, int C, std::map<int, std::vector<float>> &data, int rank, int worldsize, std::vector<int> &bid_to_rank, int *bside, std::vector<int> &part, MESH_DATA &mesh_data); // determines the cell boundaries

void center_to_id(); 





#endif