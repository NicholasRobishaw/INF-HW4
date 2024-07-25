#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <ctime>
#include <memory>
#include <math.h>
#include <cstdlib>
#include "main.h"

using namespace std;

#ifndef QUERIES_BL_H
#define QUERIES_BL_H

// Subject Fragment string: CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA
// queryString : NCTCACACGAGT
// Test string :  CTAACCCTAAC

class Queries_BL{
    private:
    
    char**      queryData = nullptr;

    char*       subjectData = nullptr;

    long        allocatedQuerySize,
                querySize,
                allocatedSubjectSize,
                scaffoldCount,
                subjectSize,
                hits;

    const long  MAX_FRAGMENTS        = 125000000,
                MAX_SCAFFOLDS        = 607;

    const int   FRAGMENT_SIZE        = 33;

    int         mismatchMax,
                matchScore,
                mismatchScore,
                gapPenalty,
                initialSeed;
    
    public:

    Queries_BL();
    ~Queries_BL();
    void search_BL(bool isQuery, long searchMax);
    void blast(string subjectFragment, int size);
    int check_Expansion(int value, int limit, bool left);
    string copy_String(string fromString, int low, int high);
    long search_Seed(hash_Table** headPtr, string queryString, int size);
    hash_Table** initialize_Hash( int size);
    void hash_Constructor(hash_Table** headPtr, string fragment_String, int size, int subjectIndex);
    hash_Table* new_Node(hash_Table** headPtr, int hash_Index, string fragment_String, int subjectIndex);
    void hash_Deconstructor(hash_Table** headPtr, int size);
    void free_LL(hash_Table* current_Ptr);
    unsigned int radix_Notation(string fragment_Str, int size);
    int character_Value(char letter);
    int neddleman_Wunsch(string oneString, string queryString);
    string random_Fragment_Str(int size);
    long random_Index(int size);
    bool read_Query(const string& file_Name);
    void query_Constructor(const string newQuery);
    bool file_reader(const string& file_Name);
    void genome_Constructor(const string toAdd, long catIndex);
    void query_Deconstructor(char** queryPtr, int size);
    void genome_Deconstructor(char* subjectPtr);
    int getHits();
};

#endif