#ifndef PATITIONING_H
#define PATITIONING_H

#include "message.h"
#include <map>


void partition(bbounds &dom, int C, std::map<int, std::vector<float>> data, int rank, int worldsize, std::vector<int> &bid_to_rank); // determines the cell boundaries
void assign(); // assigns cells to processes

void center_to_id(); 





#endif