#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <ctime>
#include <memory>
#include <math.h>
#include <cstdlib>
#include "main.h"

#ifndef QUERIES_NW_H
#define QUERIES_NW_H

using namespace std;


class Queries_NW{
    private:
    
    char**      queryData            = nullptr;

    char*       subjectData          = nullptr;

    long        allocatedQuerySize,
                querySize,
                allocatedSubjectSize,
                scaffoldCount,
                subjectSize,
                hits;

    const long  MAX_FRAGMENTS        = 1000000,
                MAX_SCAFFOLDS        = 607;

    const int   FRAGMENT_SIZE        = 33;

    int         mismatchMax,
                matchScore,
                mismatchScore,
                gapPenalty;

    
    public:
    Queries_NW();
    ~Queries_NW();
    int neddleman_Wunsch(string oneString, char queryString[]);
    void search_NW(bool isQuery, long searchMax);
    string random_Fragment_Str(int size);
    long random_Index(bool isQuery);
    bool read_Query(const string& file_Name);
    void query_Constructor(const string newQuery);
    bool file_reader(const string& file_Name);
    void genome_Constructor(const string toAdd, long catIndex);
    void query_Deconstructor(char** queryPtr, int size);
    void genome_Deconstructor(char* subjectPtr);
    int getHits() const;
};

#endif