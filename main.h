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
#define MAIN_H

struct hash_Table{
    string seed;
    int subjectIndex;
    hash_Table* nextSeed;

};  

#endif