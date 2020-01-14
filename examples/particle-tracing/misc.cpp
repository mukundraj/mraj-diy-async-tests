#include "misc.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

std::string split_filename(std::string str){
	std::size_t found = str.find_last_of("/\\");
	// std::string fullpath = str;
	return str.substr(found+1);
}


void pvi(std::vector<int> &v, int n) // print vector int
{   
    fprintf(stderr, "vector<int>: ");
    if (n>0){
        for (size_t i=0; i<n; i++){ 
                fprintf(stderr, " %d", v[i]); 
        }
    }else{
        for (size_t i=0; i<v.size(); i++){ 
                fprintf(stderr, " %d", v[i]); 
        }
    }

    fprintf(stderr, "\n");
} 

void pvi(std::vector<double> &v, int n) // print vector int
{   
    fprintf(stderr, "vector<double>: ");
    if (n>0){
        for (size_t i=0; i<n; i++){ 
                fprintf(stderr, " %f", v[i]); 
        }
    }else{
        for (size_t i=0; i<v.size(); i++){ 
                fprintf(stderr, " %f", v[i]); 
        }
    }

    fprintf(stderr, "\n");
} 



// std::vector<std::vector<int>> read_csv(const std::string &filename){

// 	std::vector<std::vector<int>> x;
// 	return x;
// }


std::vector<std::vector<int> > read_csv(const std::string &filename)
{
	std::ifstream file(filename.c_str());
 	
 	std::string line, field;
 
 	std::vector< std::vector<int> > array;  // the 2D array
    std::vector<int> v;                // array of values for one line only

	
    while ( getline(file,line) )    // get next line in file
    {
        v.clear();
        std::stringstream ss(line);

        while (getline(ss,field,','))  // break line into comma delimitted fields
        {
            v.push_back(std::stoi(field));  // add each field to the 1D array
        }


        // std::cout<<v.size()<<" ";
        
        array.push_back(v);  // add the 1D array to the 2D array
    }

	// Close the File
	file.close();
 
 	// std::cout<<filename;
	return array;
}


