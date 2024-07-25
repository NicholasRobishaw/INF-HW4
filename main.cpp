#include "main.h"
#include "Queries_NW.h"
#include "Queries_BL.h"

using namespace std;

int main(int argc, char* argv[]){
    time_t stop_Watch = 0;
    Queries_NW problemA;
    Queries_BL problemB;
    int hits = 0;

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
        if(problemA.read_Query(argv[2])){
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
        hits = problemA.getHits();

            // how many hits with up to 2 mismatches did you find
        cout << "Number of hits: " << hits << endl;
    }


    // otherwise read data and work in the Queries_BL class
    else if( strcmp(argv[4], "2") == 0){
        // read query data in
        if(problemB.read_Query(argv[2])){
            cout << "Successful read query completion" << endl;
        }

        else{
            cout << "Unsuccessful read query completion" << endl;
            return 1;
        }

        // check for search A
        if( strcmp(argv[5], "A") == 0){
            
            // read subject data in
            if(problemB.file_reader(argv[1]) ){
                cout << "Successful read subject completion" << endl;
                
            }
    
            else{
                cout << "Unsuccessful read subject completion" << endl;
                return 1;
            }
            
            cout << "Sub problem A selected\n";
        
            time(&stop_Watch);
            cout << "Search start at: " << ctime(&stop_Watch) << endl;
            
            // call to the search function for blast
            problemB.search_BL(true, stoi(argv[3]));
        }
            
        // otherwise assume search B
        else if( strcmp(argv[5], "B") == 0){
            
            cout << "Sub problem B slelected\n";
            
            time(&stop_Watch);
            cout << "Search start at: " << ctime(&stop_Watch) << endl;
            
            // call to the search function for blast
            problemB.search_BL(false, stoi(argv[3]));
        }
            
        time(&stop_Watch);
        cout << "Search Completion at: " << ctime(&stop_Watch) << endl;

        // display
        hits = problemB.getHits();

            // how many hits with up to 2 mismatches did you find
        cout << "Number of hits: " << hits << endl;
      
    }


    time(&stop_Watch);
    cout << "Program End at: " << ctime(&stop_Watch) << endl;

    return 0;
}           