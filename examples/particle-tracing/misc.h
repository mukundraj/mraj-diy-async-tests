#ifndef MISC_H
#define MISC_H

#include <string>
#include <vector>
#include <set>
#include <stdint.h>
#include <limits.h>



std::string split_filename(std::string str);

#define dprint(format, ...) \
     fprintf(stderr, format " %s %d \n", ## __VA_ARGS__, split_filename(__FILE__).c_str(), __LINE__)

void pvi(std::vector<int> &v, int n=0);

void pvi(std::vector<double> &v, int n=0);

//#define dprint(format, ...) \
	if (0) fprintf(stderr, format " %s %d \n", ## __VA_ARGS__, split_filename(__FILE__).c_str(), __LINE__)
    

std::vector<std::vector<int>> read_csv(const std::string &filename);


#if SIZE_MAX == UCHAR_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
   #error "what is happening here?"
#endif





#endif
