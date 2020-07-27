#ifndef IO_HPP
#define IO_HPP

#include <vector>
#include <map>
#include <diy/mpi.hpp>

void read_data(diy::mpi::communicator& world, const char*  infile_, std::map<int, std::vector<float>> &data, std::vector<int> &weights, int C, int *bside);

#endif