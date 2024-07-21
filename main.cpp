#include "main.h"

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

    const long  MAX_FRAGMENTS        = 125000000,
                MAX_SCAFFOLDS        = 607;

    const int   FRAGMENT_SIZE        = 33;

    int         mismatchMax,
                matchScore,
                mismatchScore,
                gapPenalty;

    public:
    // constructor function
    Queries_NW(){
        queryData            = nullptr;

        subjectData          = nullptr;

        allocatedQuerySize   = 100,
        querySize            = 0,
        allocatedSubjectSize = 100,
        scaffoldCount        = 0,
        subjectSize          = 0,
        hits                 = 0;

        mismatchMax          = 2,
        matchScore           = 2,
        mismatchScore        = -1,
        gapPenalty           = -1;

    }

    // deconstructor
    ~Queries_NW(){
        query_Deconstructor(queryData, querySize);
        genome_Deconstructor(subjectData);
    }

    // Accessors (getters)
    // ------------------------------------------------------------------------------------------------------ 

    // Function to return back the hits value
        // inline meanes that the function is all on one line 
        // the first const states that the value returned will be a const value and cannot be changed outside
        // the int& means to send the reference to the value not send a copy - helps with efficiency
        // the const after function name states that this function will not change the variable
    inline const int& Queries_NW::getHits() const {return hits;}

    const string Queries_NW::toString() const{
        return to_string(hits);

    }

    // ------------------------------------------------------------------------------------------------------ 



    // NW function to compare two n-mer squences, returning similarity score of the best alignment
    int Queries_NW::neddleman_Wunsch( string oneString, string queryString){
        
        // cout << "In NW function\n";
        // cout << "Test string (oneString) : " << oneString << endl;
        // cout << "queryString             : " << queryString << endl;
        
        
        // going to be using a scoring matrix to determine simularity score
        int rowSize = oneString.length() + 1, 
            colSize = queryString.length() + 1, 
            row, col,
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

        if( scoreMatrix[row-1][col-1] >= ((rowSize * matchScore) - abs(3 * matchScore))){
            return 1;
        }


        return 0;
    }


    // search function for a given n-mer function, returning the best possible match(based on similarity score)
        // pass in the random n-mers array from the subject data set
    void Queries_NW::search_NW(bool isQuery, long searchMax){
        long index, counter=0, subjectIndex, queryIndex;
        string testStr="";
        time_t stop_Watch = 0;


        cerr << "Entered search function\n";
        time(&stop_Watch);
        cerr << "Entered Search function at: " << ctime(&stop_Watch) << endl;
        // START TIMER HERE

        
        // loop through until hitting the search max
        while(counter < searchMax){
            
            cerr << "on iteration: " << counter << endl;
            time(&stop_Watch);
            cerr << "   at: " << ctime(&stop_Watch) << endl;
            
            
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
            
            counter++;
        }


        // END TIMER HERE

    }


    // function to generate a random test string
    string Queries_NW::random_Fragment_Str( int size){
        string randomStr= "";
        char alphabet[5] = {'A', 'C', 'G', 'N', 'T'};

        for(int index=0; index < size; index++){
            randomStr += alphabet[rand() % 5];
        }

        return randomStr;
    }


    // function to get random n-mers index (pass in desired dataset)
    long Queries_NW::random_Index(bool isQuery){
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
    bool Queries_NW::read_Qurey(const string& file_Name){
        ifstream file(file_Name);
        string currentLine;
        int queryNum = 0;
        //time_t stop_Watch = 0;
        
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
    void Queries_NW::query_Constructor(const string newQuery){
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
    bool Queries_NW::file_reader(const string& file_Name) {
        // variables
        ifstream file(file_Name);
        string currentScaffoldName;
        string currentLine;
        string tempGenome = "";
        long prevSize;
        unsigned long index;
        //time_t stop_Watch = 0;

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
    void Queries_NW::genome_Constructor(const string toAdd, long catIndex) {
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
    void Queries_NW::query_Deconstructor(char** queryPtr, int size){
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
    void Queries_NW::genome_Deconstructor(char* subjectPtr){
        delete[] subjectPtr;
        subjectPtr=nullptr;
    }

};

class Queries_BL{
    private:
    
    char**      queryData;

    char*       subjectData;

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

    // constructor function
    Queries_BL(){
        queryData            = nullptr;

        subjectData          = nullptr;

        allocatedQuerySize   = 100,
        querySize            = 0,
        allocatedSubjectSize = 100,
        scaffoldCount        = 0,
        subjectSize          = 0,
        hits                 = 0;

        mismatchMax          = 2,
        matchScore           = 2,
        mismatchScore        = -1,
        gapPenalty           = -1,
        initialSeed          = 11;

    }

    // deconstructor function
    ~Queries_BL(){
        query_Deconstructor(queryData, querySize);
        genome_Deconstructor(subjectData);
    }

    // Accessors (getters)
    // ------------------------------------------------------------------------------------------------------ 

    // Function to return back the hits value
        // inline meanes that the function is all on one line 
        // the first const states that the value returned will be a const value and cannot be changed outside
        // the int& means to send the reference to the value not send a copy - helps with efficiency
        // the const after function name states that this function will not change the variable
    inline const int& Queries_BL::getHits() const {return hits;}

    const string Queries_BL::toString() const{
        return to_string(hits);

    }

    // ------------------------------------------------------------------------------------------------------ 


    // Modifiers (setters)
    // ------------------------------------------------------------------------------------------------------ 

    void setHits(){

    }
    // ------------------------------------------------------------------------------------------------------ 


    // Member Functions
    // ------------------------------------------------------------------------------------------------------ 

    // search loop
        // loop through for N amount and generate k-mers from subject or from random generation
    void Queries_BL::search_NW(bool isQuery, long searchMax){
        long index, counter=0, subjectIndex, queryIndex;
        string testStr="";
        time_t stop_Watch = 0;


        cerr << "Entered search function\n";
        time(&stop_Watch);
        cerr << "Entered Search function at: " << ctime(&stop_Watch) << endl;
        // START TIMER HERE

            
        // check for subproblem a
        if(isQuery){
            // randomly pick an index in the subject dataset
            index = random_Index(false, searchMax);

            // create our test string
            for(subjectIndex=index; subjectIndex<searchMax; subjectIndex++){
                testStr += subjectData[subjectIndex];
            }

        }

        // otherwise assume problem b
        else{
            // call to function to generate a radnom n-mer fragment, set at the teststring
            testStr = random_Fragment_Str(searchMax);
        }

        // pass the test fragment to the blast algorithm
        blast(testStr, searchMax);
        

        // rest the test string
        testStr="";
        
        counter++;
        


        // END TIMER HERE

    }


    // blast


        // TODO: Fix the loop search process for incrementing the seed size left and right
        //       Did not do the seed expanstion for the query size of things
        //       Test the max mismatch score to see if it will work
        //       Add check for the see to match at the complete fragment size - possible change in loop
        //       optimize function 
        //       make test string creation its own modular function

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
        int    subjectIndex = 0, 
               queryIndex   = 0,
               currentSeed  = 0,
               index,
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
        initialize_Hash(headPtr, size);

            // when storing the seed in the hash table store the index of where it is located in the genome

        for(index=0; index < size-initialSeed; index++){
            
            testFragment = (subjectFragment, index, initialSeed);

            hash_Constructor(headPtr, testFragment, initialSeed, index);

            testFragment="";
        }

        // loop through the query
        for(index=0; index<querySize; index++){
            // at the current query create, break down the query fragment into seed sizes(index to index+seedSize) - for loop possibly another hash table
            for(count =0; count< FRAGMENT_SIZE-1-initialSeed; count++){
                // check each break for a match in the hash table
                // potential issue here with the 2d array pass in !!!!!!!!!!!!!!!!!!!!!
                queryFragment = copy_String(queryData[index], count, initialSeed);

                // iterate through all seed matches
                checkIndex = search_Seed(headPtr, queryFragment, size);
                subjectLeft = checkIndex;
                subjectRight = checkIndex += initialSeed;
                queryLeft = count;
                queryRight = count + initialSeed;

                    if(checkIndex != -1){

                        // at the index in the genome
                        // create the test string from the index at the genome
                        testFragment = copy_String(subjectFragment, subjectLeft, subjectRight+1);
                    
                    // attempt to expand to the left and right if its in bounds in both strings
                        // append left and right indexed of the queryData[index] to queryFragment
                        while( run && loopCounter < FRAGMENT_SIZE + 3){

                            // expand left and right in the subject
                            subjectLeft = check_Expansion(subjectLeft, 0, true);
                            subjectRight = check_Expansion(subjectRight, size-initialSeed, false);

                            // attempt to expand query fragment left and right
                            queryLeft = check_Expansion(queryLeft, 0, true);
                            queryRight = check_Expansion(queryRight, FRAGMENT_SIZE - initialSeed, false);

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
                            else if( testStrLen == FRAGMENT_SIZE - 1){
                                // increment the hit count and stop loop
                                hits++;
                                run = false;
                            }

                            loopCounter++;
                        }

                    }

                loopCounter = 0;
                run = true;
                testFragment  = "";
                queryFragment = "";
            }       
        }    

        // clear the hash tables?
        hash_Deconstructor(headPtr, size);
    }


    int Queries_BL::check_Expansion(int value, int limit, bool left){
        if( left ){
            if(value >= limit){
                return value - 1;
            }

            return value;
        }

        if(value <= limit){
            return value + 1;
        }

        return value;
    }


    // function to copy a substring from another function
    string Queries_BL::copy_String(string fromString, int low, int high){
        string output ="";

        for(int index=0; index<high; index++){
            output += fromString[index];
        }

        return output;
    }

    // function to search for the seed in the hash table
    long Queries_BL::search_Seed(hash_Table** headPtr, string queryString, int size){
        int hash_Index = radix_Noation(queryString) % size;
        hash_Table* tempPtr = headPtr[hash_Index];

        while(tempPtr!=nullptr){
            if(queryString == tempPtr->seed ){
                return tempPtr->subjectIndex;
            }
        }
    
        return -1;
    }


        // initialize hash table
    void Queries_BL::initialize_Hash(hash_Table** headPtr, int size){

        // iterate through the hash table and set each index to NULL
        for(long index=0; index < size; index++){
            headPtr[index] = nullptr;
        }

    }

    // function to populate the hash table ( Avoids adding duplicate fragments to help runtime)
    void Queries_BL::hash_Constructor(hash_Table** headPtr, string fragment_String, int size, int subjectIndex) {
        unsigned int fragment_Radix = radix_Noation(fragment_String);
        int hash_Index = fragment_Radix % size;
        
        // Traverse the linked list at hash_Index to check for existing fragment
        hash_Table* current = headPtr[hash_Index];
        
        // Create new node and insert at hash_Index
        headPtr[hash_Index] = new_Node(headPtr, hash_Index, fragment_String, subjectIndex);
    }


    // function to add a new node that Link to the hash table
    hash_Table* Queries_BL::new_Node(hash_Table** headPtr, int hash_Index, string fragment_String, int subjectIndex){
        // create a new fragment node
        hash_Table* new_Node_Ptr = new hash_Table;

        // create fragment string
        new_Node_Ptr->seed = fragment_String;

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
    unsigned int Queries_BL::radix_Noation(string fragment_Str) {
        unsigned int radix_Conversion = 0;
        int radix_Base = 5;
        int character_Num;
        
        // iterate thorugh the string from 
        for (int index = 0; index < initialSeed; index++) {
            // get the character value
            character_Num = character_Value(fragment_Str[index]);
            
            //calculate the radix conversion for the character
            radix_Conversion += (character_Num * pow(radix_Base, initialSeed - 1 - index));
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
    int Queries_BL::neddleman_Wunsch( string oneString, string queryString){
        
        // cout << "In NW function\n";
        // cout << "Test string (oneString) : " << oneString << endl;
        // cout << "queryString             : " << queryString << endl;
        
        
        // going to be using a scoring matrix to determine simularity score
        int rowSize = oneString.length() + 1, 
            colSize = queryString.length() + 1, 
            row, col,
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


     // function to generate a random test string
    string Queries_BL::random_Fragment_Str( int size){
        string randomStr= "";
        char alphabet[5] = {'A', 'C', 'G', 'N', 'T'};

        for(int index=0; index < size; index++){
            randomStr += alphabet[rand() % 5];
        }

        return randomStr;
    }


    // function to get random n-mers index (pass in desired dataset)
    long Queries_BL::random_Index(bool isQuery, int size){
        long datasetSize;

        // take the whole size of the desired dataset
        if(isQuery){
            datasetSize = querySize;
        }
        
        // NOTE: if its the subject dataset make sure to subtract the n-mer size to account for not going over array size
        else{
            datasetSize = subjectSize-size;
        }

        // return a random number from 0 to the size of the desired dataset
        return rand() % datasetSize;
    }

    // Search function to search for a given n-mer within the Query_BL class, returning the best possible match (based on simialrity score)

    // function to read the query data in
    bool Queries_BL::read_Qurey(const string& file_Name){
        ifstream file(file_Name);
        string currentLine;
        int queryNum = 0;
        //time_t stop_Watch = 0;
        
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
        string currentScaffoldName;
        string currentLine;
        string tempGenome = "";
        long prevSize;
        unsigned long index;
        //time_t stop_Watch = 0;

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
    void Queries_BL::genome_Constructor(const string toAdd, long catIndex) {
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
    void Queries_BL::querey_Deconstructor(char** queryPtr, int size){
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


};




int main(int argc, char* argv[]){
    time_t stop_Watch = 0;
    Queries_NW problemA;

    time(&stop_Watch);
    cout << "Program Start at: " << ctime(&stop_Watch) << endl;

    //./my_program human.txt human_reads_125_32.fa <search size> <problem> <-Subpart>

    // check if the program has the appropriate number of parameter
    // create error handler for no input file argument
    if( argc < 5 ){
        cout << "Error please input the correct command in\n Program End";
        return 1;
    }


    // check for problem A or B
    cout << "Argv 1 = " << argv[1] << endl;
    cout << "Argv 2 = " << argv[2] << endl;
    cout << "Argv 3 = " << argv[3] << endl;
    cout << "Argv 4 = " << argv[4] << endl;
    cout << "Argv 5 = " << argv[5] << endl;

    // if problem A then read data and work in the Queries_NW class
    if(strcmp(argv[4], "1") == 0){
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
        if( strcmp(argv[5], "A") == 0){
            
            cout << "Sub problem A selected\n";
            
            problemA.search_NW(true, stoi(argv[3]));
        }
            
        // otherwise assume search B
        else if( strcmp(argv[5], "B") == 0){
            
            cout << "Sub problem B slelected\n";
            
            problemA.search_NW(false, stoi(argv[3]));
        }
            
        time(&stop_Watch);
        cout << "Search Completion at: " << ctime(&stop_Watch) << endl;

        // display

            // how many hits with up to 2 mismatches did you find
        cout << "Number of hits: " << problemA.toString() << endl;
 
            // how long did the search take


    }


    // otherwise read data and work in the Queries_BL class
    else if( strcmp(argv[3], "2") == 0){

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
