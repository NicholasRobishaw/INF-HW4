#include "Queries_BL.h"

using namespace std;

#ifndef QUERIES_BL_CPP
#define QUERIES_BL_CPP

    // constructor function
    Queries_BL::Queries_BL(){
        queryData            = nullptr;

        subjectData          = nullptr;

        allocatedQuerySize   = 125000000,
        querySize            = 0,
        allocatedSubjectSize = 3000000000, // 3000000000
        scaffoldCount        = 0,
        subjectSize          = 0,
        hits                 = 0;

        mismatchMax          = 2,
        matchScore           = 2,
        mismatchScore        = -1,
        gapPenalty           = -1,
        initialSeed          = 7;

    }

    // deconstructor function
    Queries_BL::~Queries_BL(){
        query_Deconstructor(queryData, querySize);
        genome_Deconstructor(subjectData);
    }


    // Function to return back the hits value
        // inline meanes that the function is all on one line 
        // the first const states that the value returned will be a const value and cannot be changed outside
        // the int& means to send the reference to the value not send a copy - helps with efficiency
        // the const after function name states that this function will not change the variable
    int Queries_BL::getHits(){
        return hits;
        
    }
        

    // search loop
        // loop through for N amount and generate k-mers from subject or from random generation
    void Queries_BL::search_BL(bool isQuery, long searchMax){
        long index = 0, counter=0, subjectIndex;
        string testStr="";
        time_t stop_Watch = 0;

        // check for subproblem a
        if(isQuery){
            // randomly pick an index in the subject dataset
            index = rand() % (subjectSize - searchMax);

            // create our test string
            for(subjectIndex=index; subjectIndex<searchMax; subjectIndex++){
                testStr += subjectData[subjectIndex];
            }
        }

        // otherwise assume problem b
        else{
            // call to function to generate a radnom n-mer fragment, set at the teststring
            testStr = random_Fragment_Str(searchMax);
            // testStr = "NCTCACACGAGTGAAAGTTTTGTCAAGTTTTAATTAACCACTCTTAACCAAGGGAATTGTCTAA";
        }

        // pass the test fragment to the blast algorithm
        blast(testStr, searchMax);

        // reset the test string
        testStr="";
        counter++;
    }


    // blast
        // seed based search
        // hueristic -no guarantee of a match but will be much faster
        // break down the fragment into k-mers size ( from subject or random such that k < n)
        // store location in a hash table ( each possible fragment of seed size to store in the has table) indexing
        // break up the query into k-mers (k < n) and store location of each k-mer in the query
        // search for k-mers of query in k-mers of subject
        // if there was atleast one match from the hash table
            // extend the seeds using local or global alignment algorithm (NW)
            // check left (-n from seed location) if not at index 0
            // check right (+n from seed location) if not at the last value - seed
    void Queries_BL::blast(string subjectFragment, int size){
        int    index,
               count,
               checkIndex = 0,
               subjectLeft = 0,
               subjectRight = 0,
               queryLeft = 0,
               queryRight = 0,
               testStrLen,
               mismatchMax,
               loopCounter = 0;

        bool run = true;

        hash_Table** headPtr = nullptr;

        string queryFragment ="";
        string testFragment  ="";
        
        // create the hash table of seed sizes of the subject fragment string
        headPtr = initialize_Hash( size);
        
        // when storing the seed in the hash table store the index of where it is located in the genome
        for(index=0; index < size-initialSeed; index++){ 
            // creating seed
            testFragment = copy_String(subjectFragment, index, index + initialSeed);

            // adding to hash table
            hash_Constructor(headPtr, testFragment, size, index);

            // reset seed buffer
            testFragment="";
        }

        // loop through the query
        for(index=0; index<querySize; index++){
            // at the current query create, break down the query fragment into seed sizes(index to index+seedSize) - for loop possibly another hash table
            for(count =0; count< FRAGMENT_SIZE-1-initialSeed; count++){
                // check each break for a match in the hash table
                queryFragment = copy_String(queryData[index], count, initialSeed);

                // iterate through all seed matches
                checkIndex = search_Seed(headPtr, queryFragment, size);
                subjectLeft = checkIndex;
                subjectRight = checkIndex + initialSeed;
                queryLeft = count;
                queryRight = count + initialSeed;

                // check if the query seed was found in the hash table
                if(checkIndex >= 0){
                    // attempt to expand to the left and right if its in bounds in both strings
                    // append left and right indexed of the queryData[index] to queryFragment
                    while( run && loopCounter < FRAGMENT_SIZE){
                        // expand left and right in the subject
                        subjectLeft = check_Expansion(subjectLeft, 0, true);
                        subjectRight = check_Expansion(subjectRight, size, false);

                        // attempt to expand query fragment left and right
                        queryLeft = check_Expansion(queryLeft, 0, true);
                        queryRight = check_Expansion(queryRight, FRAGMENT_SIZE, false);

                        // remake the test strings with the updated
                        queryFragment = copy_String(queryData[index], queryLeft, queryRight );
                        testFragment = copy_String(subjectFragment, subjectLeft, subjectRight);

                        // get the current length of test subject string
                        testStrLen = testFragment.length();

                        // calc floor score ( allowing for up to 2 mismatches )
                        mismatchMax = (testStrLen * matchScore) - abs(3 * matchScore);

                        // get the similarity score from NW function
                        // if the similarity score returns less than the match score wanted then go then break and try next seed match if possible
                        if(neddleman_Wunsch(queryFragment, testFragment) < mismatchMax){
                            // go to the next loop 
                            run = false;
                        }
                        
                        // win condition
                        // otherwise check if at the fragment length
                        else if( testStrLen == FRAGMENT_SIZE - 2){
                            // increment the hit count and stop loop
                            hits++;
                            run = false;
                        }

                        loopCounter++;
                    }
                }

            // reset values
            loopCounter = 0;
            run = true;
            testFragment  = "";
            queryFragment = "";
            }       
        }    

        // clear the hash table
        hash_Deconstructor(headPtr, size);
    }


    // function to return updated value if it can expand
    int Queries_BL::check_Expansion(int value, int limit, bool left){
        // check if we need to expand to the left
        if( left ){
            // check if we can go -1 to the left
            if(value > limit){
                return value - 1;
            }

            // otherwise just return the limit back
            return limit;
        }

        // other wise assume right and check if we can expand to the right by 1
        if(value < limit){
            return value + 1;
        }

        // otherwise just return the right limit
        return limit;
    }


    // function to copy a substring from another function
    string Queries_BL::copy_String(string fromString, int low, int high){
        string output ="";

        // iterate through the input string and create a substring
        for(int index=low; index<high; index++){
            output += fromString[index];
        }

        // return the substring
        return output;
    }


    // function to search for the seed in the hash table
    long Queries_BL::search_Seed(hash_Table** headPtr, string queryString, int size){
        int hash_Index = radix_Notation(queryString, initialSeed) % size;
        hash_Table* tempPtr = headPtr[hash_Index];

        // iterate throught he hash table
        while(tempPtr!=nullptr){
            // check if the seeds match
            if(queryString == tempPtr->seed ){
                // return the index of the seed in the subject data
                return tempPtr->subjectIndex;
            }
            
            // iterate to next node
            tempPtr = tempPtr->nextSeed;
        }
    
        // otherwise return failure
        return -1;
    }


    // initialize hash table
    hash_Table** Queries_BL::initialize_Hash( int size){
        hash_Table** headPtr = new hash_Table*[size];

        // iterate through the hash table and set each index to NULL
        for(long index=0; index < size; index++){
            headPtr[index] = nullptr;
        }
    
        // return new hash table pointer
        return headPtr;
    }


    // function to populate the hash table ( Avoids adding duplicate fragments to help runtime)
    void Queries_BL::hash_Constructor(hash_Table** headPtr, string fragment_String, int size, int subjectIndex) {
        // get the radix notation of the fragment string
        unsigned int fragment_Radix = radix_Notation(fragment_String, initialSeed);
        
        // calculate the hash index
        int hash_Index = fragment_Radix % size;
        
        // Create new node and insert at hash_Index
        headPtr[hash_Index] = new_Node(headPtr, hash_Index, fragment_String, subjectIndex);
    }


    // function to add a new node that Link to the hash table
    hash_Table* Queries_BL::new_Node(hash_Table** headPtr, int hash_Index, string fragment_String, int subjectIndex){
        // create a new fragment node
        hash_Table* new_Node_Ptr = new hash_Table;

        // create fragment string
        new_Node_Ptr->seed = fragment_String;

        // add in the index the seed is located at in the subject data
        new_Node_Ptr->subjectIndex=subjectIndex;

        // set next pointers to previous head ptrs
        new_Node_Ptr->nextSeed = headPtr[hash_Index];

        // return ptr to new node
        return new_Node_Ptr;
    }


    // deconstructor
    void Queries_BL::hash_Deconstructor(hash_Table** headPtr, int size){
        // iterate through the hash table
        for(int index = 0; index < size; index++){
            // destroy the LL
            free_LL(headPtr[index]);
        }

        // destroy the main array
        delete[] headPtr;
    }


    // LL destroyer
    void Queries_BL::free_LL( hash_Table* current_Ptr){
        // check if we are not at the end of the LL
        if( current_Ptr != nullptr ){
            // check if we are at the end of the LL
            if( current_Ptr->nextSeed != nullptr){
                // recurse to next node
                free_LL(current_Ptr->nextSeed);
            }

            // free the current node
            delete current_Ptr;
        }
    }


    // function to convert the fragment string to a unsigned radix number
    unsigned int Queries_BL::radix_Notation(string fragment_Str, int size) {
        unsigned int radix_Conversion = 0;
        int radix_Base = 5;
        int character_Num;
        
        // iterate thorugh the string from 
        for (int index = 0; index < initialSeed; index++) {
            // get the character value
            character_Num = character_Value(fragment_Str[index]);
            
            //calculate the radix conversion for the character
            radix_Conversion += (character_Num * pow(radix_Base, size - 1 - index));
        }
    
        return radix_Conversion;
    }
    
    // returns the radix number for the letter
    int Queries_BL::character_Value(char letter){
        int index;
        
        switch (letter) {
            case 'A':
                index = 0;
                break;
            case 'C':
                index = 1;
                break;
            case 'G':
                index = 2;
                break;
            case 'T':
                index = 3;
                break;
            default:
                index = 4;
        }
        
        return index;
    }


    // NW function to compare two n-mer sequences, returning the similarity score of the best alignment
        // use location of seed in query to figure out cut location
    int Queries_BL::neddleman_Wunsch( string queryString, string testFragment){
        // going to be using a scoring matrix to determine simularity score
        int rowSize = queryString.length() + 1, 
            colSize = testFragment.length() + 1, 
            row, col,
            match, mismatch, gapLeft, gapUp,
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
                else if( queryString[row-1] == testFragment[col-1]){
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
                    }
                    else if( gapLeft >= gapUp){
                        scoreMatrix[row][col] = gapLeft;
                    }
                    else{
                        scoreMatrix[row][col] = gapUp;
                    }
                }
            }
        }

        // return back the similarity score at the bottom right
        return scoreMatrix[row-1][col-1];
    }


     // function to generate a random test string
    string Queries_BL::random_Fragment_Str( int size){
        string randomStr= "";
        char alphabet[5] = {'A', 'C', 'G', 'N', 'T'};

        // randomize the seed based off the time
        srand(time(NULL));

        // concatenate random letters from the alphabet together
        for(int index=0; index < size; index++){
            randomStr += alphabet[rand() % 5];
        }

        return randomStr;
    }


    // function to get random n-mers index (pass in desired dataset)
    long Queries_BL::random_Index(int size){
        return rand() % (subjectSize - size);
    }

    // Search function to search for a given n-mer within the Query_BL class, returning the best possible match (based on simialrity score)

    // function to read the query data in
    bool Queries_BL::read_Query(const string& file_Name){
        ifstream file(file_Name);
        string currentLine;
        int queryNum = 0;
        
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
    void Queries_BL::query_Constructor(const string newQuery){
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
    bool Queries_BL::file_reader(const string& file_Name) {
        // variables
        ifstream file(file_Name);
        string current_Scaffold_Name;
        string current_Line;
        string temp_Genome = "";
        long prev_genome_size;
        unsigned long index;

        // check if the file can be opened
        if (!file.is_open()) {
            // display error and exit
            cerr << "Error: Unable to open file " << file_Name << "\n";
            return false;
        }
    
        // initialize starting array, set to 3 billion to save on runtime complexity
        subjectData = new char[allocatedSubjectSize];
        
        // iterate through file contents line by line
        while ( scaffoldCount < MAX_SCAFFOLDS && getline(file, current_Line)) {
            // check if we hit a header
            if(current_Line[0] == '>' && !temp_Genome.empty()){
                
                // increment scaffold count by 1
                scaffoldCount++;
                
                // calculate new size for genome array
                prev_genome_size = subjectSize;
            
                // resize the genome array to fit last scaffold and copy data into new array
                genome_Constructor(temp_Genome, prev_genome_size);

                // reset the temp genome string to empty
                temp_Genome = "";
            }

            else if (!current_Line.empty() && current_Line[0] != '>') {
                // Iterate through current line and append characters to temp_Genome
                for (index = 0; index < current_Line.length() && index < 80; index++) {
                    if(                  current_Line[index] == 'A' 
                                      || current_Line[index] == 'C' 
                                      || current_Line[index] == 'G' 
                                      || current_Line[index] == 'T' 
                                      || current_Line[index] == 'N'){
                    
                        // Append character to temp_Genome
                        temp_Genome += current_Line[index];
                    }
                }
            }
        }

        // post append information of last scaffold
        if (!temp_Genome.empty()) {
            // calculate new genome array size
            prev_genome_size = subjectSize;
            
            // resize the genome array to fit last scaffold and copy data into new array
            genome_Constructor(temp_Genome, prev_genome_size);
        }

        // close the file
        file.close();

        // return sucess
        return true;
        
    }// end of file reader function


    // helper function for resizing and adding a new scaffold into the subject array
    void Queries_BL::genome_Constructor(const string to_Add, long cat_Index) {
        long index;
        long length = to_Add.length();
        
        // check if the array needs to be resized
        if(subjectSize + length >= allocatedSubjectSize){
            // calculate the new size of the genome
            long new_Size = subjectSize + length + 1;
            
            // allocate new array size
            char* new_Genome_Arr = new char[new_Size];
            
            // copy data from old array to new one
            for(index = 0; index < subjectSize; index++){
                new_Genome_Arr[index] = subjectData[index];
            }
            
            // deallocate memory from old array
            delete[] subjectData;
            
            // set genome pointer to new array
            subjectData = new_Genome_Arr;
            
            // set the new allocated size for genome array
            allocatedSubjectSize = new_Size;
        }
        
            
        for(index = 0; index < length; index++){
            //add in the new data 
            subjectData[cat_Index + index] = to_Add[index];
        }
        
        // set the new genome size
        subjectSize = cat_Index + to_Add.length();
    }


    // deconstructor function for query
    void Queries_BL::query_Deconstructor(char** queryPtr, int size){
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
    void Queries_BL::genome_Deconstructor(char* subjectPtr){
        delete[] subjectPtr;
        subjectPtr=nullptr;
    }

#endif