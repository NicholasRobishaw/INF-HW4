#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <ctime>
#include <memory>
#include <math.h>
#include <cstdlib>

using namespace std;

#ifndef MAIN_H

    class Queries_NW{
        int neddleman_Wunsch( string oneString, string queryString);
        void search_NW(bool isQuery, long searchMax);
        string random_Fragment_Str( int size);
        long random_Index(bool isQuery);
        bool read_Qurey(const string& file_Name);
        void query_Constructor(const string newQuery);
        bool file_reader(const string& file_Name);
        void genome_Constructor(const string toAdd, long catIndex);
        void qurey_Deconstructor(char** queryPtr, int size);
        void genome_Deconstructor(char* subjectPtr);

    };

    class Queries_BL{

    };


#endif