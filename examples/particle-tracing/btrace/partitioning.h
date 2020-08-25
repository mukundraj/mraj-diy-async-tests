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


void assign(diy::mpi::communicator &world, std::map<int, std::vector<float>> &data, std::map<int, std::vector<BEndPt>> particles,  std::vector<int> &weights, std::vector<int> &partn, MESH_DATA &myMesh, BBlock* b, const diy::Master::ProxyWithLink &cp, const diy::RoundRobinAssigner &assigner, std::map<int,bbounds> &bounds);

void partition(diy::mpi::communicator& world, bbounds &dom, int C, std::map<int, std::vector<float>> &data, int rank, int worldsize, std::vector<int> &cid_to_rank, int *cside, std::vector<int> &part, MESH_DATA &mesh_data,  BBlock* b); // determines the cell boundaries

void center_to_id(); 

void update_cid_to_rank_vector(diy::mpi::communicator &world, std::map<int, std::vector<float>> &data, std::vector<int> &cid_to_rank);

// populates the b->data_ghosts based on depth of ghost region N
void get_ghost_cells(diy::Master &master, const diy::RoundRobinAssigner &assigner, int N, bbounds &dom, int C, int rank);

// identify cells in ghost region of current block
void id_ghost_cells(std::map<int, std::vector<float>> &data, std::map<int, std::vector<float>> &data_ghost, std::map<int, bbounds> &bounds, bbounds &dom, int C, int rank);

// get nbrs of cell (upto 27 nbrs)
void get_nbrs_of_cell(int cid, std::vector<int>& nbrs, bbounds &bnd, bbounds &dom, int C, int rank);



#endif