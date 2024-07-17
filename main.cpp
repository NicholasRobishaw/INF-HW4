#include "main.h"

using namespace std;

class Queries_NW{
    public:
    char**      queryData            = nullptr;

    char*       subjectData          = nullptr;

    long        allocatedQuerySize,
                querySize,
                allocatedSubjectSize,
                scaffoldCount,
                subjectSize,
                hits;

    const long  MAX_FRAGMENTS,
                MAX_SCAFFOLDS;

    const int   FRAGMENT_SIZE;

    int         mismatchMax,
                matchScore           =  2,
                mismatchScore        = -1,
                gapPenalty           = -1,
                hits                 =  0;

    // constructor function


    // NW function to compare two n-mer squences, returning similarity score of the best alignment
    int neddleman_Wunsch( string oneString, string queryString){
        // going to be using a scoring matrix to determine simularity score
        int rowSize = oneString.length() + 1, 
            colSize = queryString.length() + 1, 
            row, col, prevValue = 0,
            match, mismatch, gapLeft, gapUp, numOfMis = 0,
            gapPenalty =-1, matchScore = 2, mismatchScore=-1;

        // make a score matrix values (rows = oneString length + 1, cols = otherString length + 1), initialize all vlaues to 0
        int scoreMatrix[rowSize][colSize];


        // iterate through the rows
        for(row=0; row<rowSize; row++){
            // iterate through columns
            for(col=0; col<colSize; col++){
                // if in the first row/col fill in the gap scores
                if(row == 0 && col == 0){
                    scoreMatrix[row][col] = 0;
                }

                else if( row == 0 ){
                    scoreMatrix[row][col] = scoreMatrix[row][col-1] + gapPenalty;
                }

                else if( col == 0 ){
                    scoreMatrix[row][col] = scoreMatrix[row-1][col] + gapPenalty;
                }

                // otherwise check if the character matches
                else if( oneString[row-1] == queryString[col-1]){
                    // calulate the value for a match with the matchScore var
                    match = scoreMatrix[row-1][col-1] + matchScore;

                    // calculate the value for a gap with the gapScore var
                    gapLeft = scoreMatrix[row-1][col] + gapPenalty;
                    gapUp = scoreMatrix[row][col-1] + gapPenalty;

                    // set to which is higher
                    if( match >= gapLeft && match >= gapUp){
                        scoreMatrix[row][col] = match;
                    }
                    else if( gapLeft >= gapUp){
                        scoreMatrix[row][col] = gapLeft;
                    }
                    else{
                        scoreMatrix[row][col] = gapUp;
                    }

                }

                // otherwise assume mismatch
                else{
                    // calculate the value for a mismatch with the mismatchScore var
                    mismatch = scoreMatrix[row-1][col-1] + mismatchScore;

                    // calculate the value for a gap with the gapScore var
                    gapLeft = scoreMatrix[row-1][col] + gapPenalty;
                    gapUp = scoreMatrix[row][col-1] + gapPenalty;

                    // set to which is higher
                    if( mismatch >= gapLeft && mismatch >= gapUp){
                        scoreMatrix[row][col] = mismatch;

                        // increment mismatch counter
                        numOfMis++;
                    }
                    else if( gapLeft >= gapUp){
                        scoreMatrix[row][col] = gapLeft;
                    }
                    else{
                        scoreMatrix[row][col] = gapUp;
                    }

                    // if mismatch counter is > 2 return failure
                    if( numOfMis > 2){
                        // return failure
                        return 0;
                    }
                }
            }

        }

        // here is some code to display the finished matrix
        // int maxWidth = 0;
        // for (int i = 0; i < rowSize; i++) {
        //     for (int j = 0; j < colSize; j++) {
        //         int width = std::to_string(scoreMatrix[i][j]).length();
        //         if (width > maxWidth) {
        //             maxWidth = width;
        //         }
        //     }
        // }

        // // Print the matrix with proper alignment
        // for (int i = 0; i < rowSize; i++) {
        //     for (int j = 0; j < colSize; j++) {
        //         std::string element = std::to_string(scoreMatrix[i][j]);
        //         int padding = maxWidth - element.length();
        //         std::cout << std::string(padding + 1, ' ') << element;
        //     }
        //     std::cout << std::endl;
        // }

        return scoreMatrix[row-1][col-1];
    }


    // search function for a given n-mer function, returning the best possible match(based on similarity score)
        // pass in the random n-mers array from the subject data set
    void search_NW(bool isQuery, long searchMax){
        long index, counter=0, subjectIndex, queryIndex;
        string testStr="";

        // START TIMER HERE

        
        // loop through until hitting the search max
        while(counter < searchMax){
            // check for subproblem a
            if(isQuery){
                // randomly pick an index in the subject dataset
                index = random_Index(false);

                // create our test string
                for(subjectIndex=index; subjectIndex<index+FRAGMENT_SIZE-2; subjectIndex++){
                    testStr += subjectData[subjectIndex];
                }

            }

            // otherwise assume problem b
            else{
                // call to function to generate a radnom n-mer fragment, set at the teststring
                testStr = random_Fragment_Str(FRAGMENT_SIZE-1);
            }

            // loop through the query data array
            for(queryIndex=0; queryIndex<querySize; queryIndex++){
                // get the NW fuzzy match score ( NW function will return 0 if more than 2 mismatches were found)

                // NOTE POSSIBLE BUG HERE WITH THE QUERY FRAGMENT- THE FRAGMENT IN THE QUERY DATA INCLUDES THE '\0' AT THE END OF THE FRAGMENT
                hits += neddleman_Wunsch( testStr, queryData[queryIndex]);
            }

            // rest the test string
            testStr="";
        }


        // END TIMER HERE

    }


    // function to generate a random test string
    string random_Fragment_Str( int size){
        string randomStr= "";
        char alphabet[5] = {'A', 'C', 'G', 'N', 'T'};

        for(int index=0; index < size; index++){
            randomStr += alphabet[rand() % 5];
        }

        return randomStr;
    }


    // function to get random n-mers index (pass in desired dataset)
    long random_Index(bool isQuery){
        long datasetSize;

        // take the whole size of the desired dataset
        if(isQuery){
            datasetSize = querySize;
        }
        
        // NOTE: if its the subject dataset make sure to subtract the n-mer size to account for not going over array size
        else{
            datasetSize = subjectSize-FRAGMENT_SIZE;
        }

        // return a random number from 0 to the size of the desired dataset
        return rand() % datasetSize;
    }


    // function to read in the query data
    bool read_Qurey(const string& file_Name){
        ifstream file(file_Name);
        string currentLine;
        int queryNum = 0;
        time_t stop_Watch = 0;
        
        // check if the file can be opened
        if (!file.is_open()) {
            // display error and exit
            cerr << "Error: Unable to open file " << file_Name << "\n";
            return false;
        }
        
        
        // allocate the inital array storage, to save time complexity set to max frags
        allocatedQuerySize = MAX_FRAGMENTS;
        queryData = new char*[allocatedQuerySize];
        
        // iterate through file contents line by line
        while ( querySize < MAX_FRAGMENTS && getline(file, currentLine)) {
            // check if we are at a valid fragment
            if(currentLine[0] != '>' && !currentLine.empty()){
                // increment query count
                queryNum++;
                
                // resize query array and add new fragment to next free index
                query_Constructor(currentLine);
            }
        }

        // close the file
        file.close();

        // return sucess
        return true;

    }


   // custom constructor function for resizing the mulitdimesional array
    void query_Constructor(const string newQuery){
        // check if the query size is larger than the allocated space for the array
        if(querySize >= allocatedQuerySize){
            // increment the size by 1 to save space complexity
            long newSize = querySize + 1;
            
            // create new pointer array with new size
            char** newData = new char*[newSize];
            
            // copy pointers from old array to new array
            for(long index = 0; index < querySize; index++){
                newData[index] = queryData[index];
            }
            
            // free old array memory
            delete[] queryData;
            
            // set query_Data pointer to new array
            queryData = newData;
            
            // set new size for tracking
            allocatedQuerySize = newSize;
        }
        
        // create a new memory storage for the new fragment
        queryData[querySize] = new char[FRAGMENT_SIZE];
        
        // copy the data into the 33 character space array
        strncpy(queryData[querySize], newQuery.c_str(), FRAGMENT_SIZE);
        
        // add terminating character to the end of the fragment
        queryData[querySize][FRAGMENT_SIZE-1] = '\0';
        
        // incremnt fragment count
        querySize++;
    }


    // function to read in the subject data
    bool file_reader(const string& file_Name) {
        // variables
        ifstream file(file_Name);
        string currentScaffoldName;
        string currentLine;
        string tempGenome = "";
        long prevSize;
        unsigned long index;
        time_t stop_Watch = 0;

        // check if the file can be opened
        if (!file.is_open()) {
            // display error and exit
            cerr << "Error: Unable to open file " << file_Name << "\n";
            return false;
        }

        
        // initialize starting array, set to 3 billion to save on runtime complexity
        allocatedSubjectSize = 3000000000;
        subjectData = new char[allocatedSubjectSize];
        
        // iterate through file contents line by line
        while ( scaffoldCount < MAX_SCAFFOLDS && getline(file, currentLine)) {
            // check if we hit a header
            if(currentLine[0] == '>' && !tempGenome.empty()){
                
                // increment scaffold count by 1
                scaffoldCount++;
                
                // calculate new size for genome array
                prevSize = subjectSize;
            
                // resize the genome array to fit last scaffold and copy data into new array
                genome_Constructor(tempGenome, prevSize);

                // reset the temp genome string to empty
                tempGenome = "";
            }

            else if (!currentLine.empty() && currentLine[0] != '>') {
                // Iterate through current line and append characters to temp_Genome
                for (index = 0; index < currentLine.length() && index < 80; index++) {
                    if(  currentLine[index] == 'A' 
                        || currentLine[index] == 'C' 
                        || currentLine[index] == 'G' 
                        || currentLine[index] == 'T' 
                        || currentLine[index] == 'N'){
                    
                    
                        // Append character to temp_Genome
                        tempGenome += currentLine[index];
                    }
                }
            }
        }

        // post append information of last scaffold
        if (!tempGenome.empty()) {
            // calculate new genome array size
            prevSize = subjectSize;
            
            // resize the genome array to fit last scaffold and copy data into new array
            genome_Constructor(tempGenome, prevSize);
        }

        // close the file
        file.close();

        // return sucess
        return true;
    }// end of file reader function


    // helper function for resizing and adding a new scaffold into the subject array
    void genome_Constructor(const string toAdd, long catIndex) {
        long index;
        long length = toAdd.length();
        
        // check if the array needs to be resized
        if(subjectSize + length >= allocatedSubjectSize){
            // calculate the new size of the genome
            long newSize = subjectSize + length;
            
            // allocate new array size
            char* newSubjectArr = new char[newSize];
            
            // copy data from old array to new one
            for(index = 0; index < subjectSize; index++){
                newSubjectArr[index] = subjectData[index];
            }
            
            // deallocate memory from old array
            delete[] subjectData;
            
            // set genome pointer to new array
            subjectData = newSubjectArr;
            
            // set the new allocated size for genome array
            allocatedSubjectSize = newSize;
        }
        
            
        for(index = 0; index < length; index++){
            //add in the new data 
            subjectData[catIndex + index] = toAdd[index];
        }
        
        // set the new genome size
        subjectSize = catIndex + toAdd.length();
    }


    // deconstructor function for query
    void qurey_Deconstructor(char** queryPtr, int size){
        int index;

        // iterate through array of pointers and free each index
        for(index=0;index<size;index++){
            delete[] queryPtr[index];
        }

        // free main array
        delete[] queryPtr;
        queryPtr=nullptr;
    }


    // deconstructor for deallocating the genome character array
    void genome_Deconstructor(char* subjectPtr){
        delete[] subjectPtr;
        subjectPtr=nullptr;
    }

};

class Queries_BL{
    public:
    char** queryData = nullptr;
    char* subjectData = nullptr;
    long allocatedQuerySize,
         querySize,
         allocatedSubjectSize,
         scaffoldCount,
         subjectSize;
    const long MAX_FRAGMENTS,
               MAX_SCAFFOLDS;
    const int FRAGMENT_SIZE;

    // constructor function

    // deconstructor function

    // NW function to compare two n-mer sequences, returning the similarity score of the best alignment

    // Search function to search for a given n-mer within the Query_BL class, returning the best possible match (based on simialrity score)

    // function to read the query data in
    bool read_Qurey(const string& file_Name){
        ifstream file(file_Name);
        string currentLine;
        int queryNum = 0;
        time_t stop_Watch = 0;
        
        // check if the file can be opened
        if (!file.is_open()) {
            // display error and exit
            cerr << "Error: Unable to open file " << file_Name << "\n";
            return false;
        }
        
        
        // allocate the inital array storage, to save time complexity set to max frags
        allocatedQuerySize = MAX_FRAGMENTS;
        queryData = new char*[allocatedQuerySize];
        
        // iterate through file contents line by line
        while ( querySize < MAX_FRAGMENTS && getline(file, currentLine)) {
            // check if we are at a valid fragment
            if(currentLine[0] != '>' && !currentLine.empty()){
                // increment query count
                queryNum++;
                
                // resize query array and add new fragment to next free index
                query_Constructor(currentLine);
            }
        }

        // close the file
        file.close();

        // return sucess
        return true;

    }


    // custom constructor function for resizing the mulitdimesional array
    void query_Constructor(const string newQuery){
        // check if the query size is larger than the allocated space for the array
        if(querySize >= allocatedQuerySize){
            // increment the size by 1 to save space complexity
            long newSize = querySize + 1;
            
            // create new pointer array with new size
            char** newData = new char*[newSize];
            
            // copy pointers from old array to new array
            for(long index = 0; index < querySize; index++){
                newData[index] = queryData[index];
            }
            
            // free old array memory
            delete[] queryData;
            
            // set query_Data pointer to new array
            queryData = newData;
            
            // set new size for tracking
            allocatedQuerySize = newSize;
        }
        
        // create a new memory storage for the new fragment
        queryData[querySize] = new char[FRAGMENT_SIZE];
        
        // copy the data into the 33 character space array
        strncpy(queryData[querySize], newQuery.c_str(), FRAGMENT_SIZE);
        
        // add terminating character to the end of the fragment
        queryData[querySize][FRAGMENT_SIZE-1] = '\0';
        
        // incremnt fragment count
        querySize++;
    }


    // function to read in the subject data
    bool file_reader(const string& file_Name) {
        // variables
        ifstream file(file_Name);
        string currentScaffoldName;
        string currentLine;
        string tempGenome = "";
        long prevSize;
        unsigned long index;
        time_t stop_Watch = 0;

        // check if the file can be opened
        if (!file.is_open()) {
            // display error and exit
            cerr << "Error: Unable to open file " << file_Name << "\n";
            return false;
        }

        
        // initialize starting array, set to 3 billion to save on runtime complexity
        allocatedSubjectSize = 3000000000;
        subjectData = new char[allocatedSubjectSize];
        
        // iterate through file contents line by line
        while ( scaffoldCount < MAX_SCAFFOLDS && getline(file, currentLine)) {
            // check if we hit a header
            if(currentLine[0] == '>' && !tempGenome.empty()){
                
                // increment scaffold count by 1
                scaffoldCount++;
                
                // calculate new size for genome array
                prevSize = subjectSize;
            
                // resize the genome array to fit last scaffold and copy data into new array
                genome_Constructor(tempGenome, prevSize);

                // reset the temp genome string to empty
                tempGenome = "";
            }

            else if (!currentLine.empty() && currentLine[0] != '>') {
                // Iterate through current line and append characters to temp_Genome
                for (index = 0; index < currentLine.length() && index < 80; index++) {
                    if(  currentLine[index] == 'A' 
                        || currentLine[index] == 'C' 
                        || currentLine[index] == 'G' 
                        || currentLine[index] == 'T' 
                        || currentLine[index] == 'N'){
                    
                    
                        // Append character to temp_Genome
                        tempGenome += currentLine[index];
                    }
                }
            }
        }

        // post append information of last scaffold
        if (!tempGenome.empty()) {
            // calculate new genome array size
            prevSize = subjectSize;
            
            // resize the genome array to fit last scaffold and copy data into new array
            genome_Constructor(tempGenome, prevSize);
        }

        // close the file
        file.close();

        // return sucess
        return true;
    }// end of file reader function


    // helper function for resizing and adding a new scaffold into the subject array
    void genome_Constructor(const string toAdd, long catIndex) {
        long index;
        long length = toAdd.length();
        
        // check if the array needs to be resized
        if(subjectSize + length >= allocatedSubjectSize){
            // calculate the new size of the genome
            long newSize = subjectSize + length;
            
            // allocate new array size
            char* newSubjectArr = new char[newSize];
            
            // copy data from old array to new one
            for(index = 0; index < subjectSize; index++){
                newSubjectArr[index] = subjectData[index];
            }
            
            // deallocate memory from old array
            delete[] subjectData;
            
            // set genome pointer to new array
            subjectData = newSubjectArr;
            
            // set the new allocated size for genome array
            allocatedSubjectSize = newSize;
        }
        
            
        for(index = 0; index < length; index++){
            //add in the new data 
            subjectData[catIndex + index] = toAdd[index];
        }
        
        // set the new genome size
        subjectSize = catIndex + toAdd.length();
    }


    // deconstructor function for query
    void qurey_Deconstructor(char** queryPtr, int size){
        int index;

        // iterate through array of pointers and free each index
        for(index=0;index<size;index++){
            delete[] queryPtr[index];
        }

        // free main array
        delete[] queryPtr;
        queryPtr=nullptr;
    }


    // deconstructor for deallocating the genome character array
    void genome_Deconstructor(char* subjectPtr){
        delete[] subjectPtr;
        subjectPtr=nullptr;
    }


};




int main(int argc, char* argv[]){
    time_t stop_Watch = 0;
    Queries_NW problemA;

    time(&stop_Watch);
    cout << "Program Start at: " << ctime(&stop_Watch) << endl;

    //./my_program human.txt human_reads_125_32.fa <search size> <problem> <-Subpart>

    // check if the program has the appropriate number of parameter
    // create error handler for no input file argument
    if( argc != 5 ){
        cout << "Error please input the correct command in\n Program End";
        return 1;
    }


    // check for problem A or B


    // if problem A then read data and work in the Queries_NW class
    if(*argv[3] == 'A'){
        // read query data in
        if(problemA.read_Qurey(argv[2])){
            cout << "Successful read query completion" << endl;
        }

        else{
            cout << "Unsuccessful read query completion" << endl;
            return 1;
        }

        // read subject data in
        if(problemA.file_reader(argv[1])){
            cout << "Successful read subject completion" << endl;
        }

        else{
            cout << "Unsuccessful read subject completion" << endl;
            return 1;
        }

        time(&stop_Watch);
        cout << "Search start at: " << ctime(&stop_Watch) << endl;

        // check for search A
        if( *argv[5] == 1){
            problemA.search_NW(true, *argv[4]);
        }
            
        // otherwise assume search B
        else if( *argv[5] == 2){
            problemA.search_NW(false, *argv[4]);
        }
            
        time(&stop_Watch);
        cout << "Search Completion at: " << ctime(&stop_Watch) << endl;

        // display
            // how many hits with up to 2 mismatches did you find
        cout << "Number of hits: " << problemA.hits << endl;

            // how long did the search take

        // clean up memory
        problemA.qurey_Deconstructor(problemA.queryData, problemA.querySize);
        problemA.genome_Deconstructor(problemA.subjectData);
    }


    // otherwise read data and work in the Queries_BL class
    else if( *argv[3] == 'B'){

        // read query data in

        // read subject data in

        // use word size of 11 for initial seed matching within the blast algorithm and the NW to perform seed extension

        // check for search A
            // randomly pick (10k, 100k, 1M) character long segments from the subject dataset to conduct search with blast algorithm
            // Note: iterate through all possible n-mers of your character segment when you are conductnig your search

        // otherwise assume search B

            // generate completely random (10k, 100k, 1M) character long segment to conduct fuzzy searching within the query dataset using the blast algorithm 
            // Note: iterate through all possible n-mers of your character segment when you are conductnig your search


        // display

            // how many hits with up to 2 mismatches were found

            // how long did the search take
    }


    time(&stop_Watch);
    cout << "Program End at: " << ctime(&stop_Watch) << endl;

    return 0;
}           