#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <ctime>
#include <memory>
#include <math.h>
#include <cstdlib>

using namespace std;

int neddleman_Wunsch( string oneString, string otherString);
string random_Fragment_Str( int size);

int main(){

    int test = neddleman_Wunsch("CNGTCGACGTN", "CNGTCGACGCA");

    cout << "\nSimularity Score: " << test;

    // string testStr = random_Fragment_Str(32);

    // cout << "random string generated: " << testStr << endl;


    cout << "Program End";

    return 0;
}


string random_Fragment_Str( int size){
        string randomStr= "";
        char alphabet[5] = {'A', 'C', 'G', 'N', 'T'};

        for(int index=0; index < size; index++){
            randomStr += alphabet[rand() % 5];
        }

        return randomStr;
    }

// NW function to compare two n-mer squences, returning similarity score of the best alignment
int neddleman_Wunsch( string oneString, string otherString){
    // going to be using a scoring matrix to determine simularity score
    int rowSize = oneString.length() + 1, 
        colSize = otherString.length() + 1, 
        row, col, prevValue = 0,
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
            else if( oneString[row-1] == otherString[col-1]){
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

                // if set to mismatch - increment mismatch counter
                // if mismatch counter is > 2 return failure
            }
        }


    }

    
    int maxWidth = 0;
    for (int i = 0; i < rowSize; i++) {
        for (int j = 0; j < colSize; j++) {
            int width = std::to_string(scoreMatrix[i][j]).length();
            if (width > maxWidth) {
                maxWidth = width;
            }
        }
    }

    // Print the matrix with proper alignment
    for (int i = 0; i < rowSize; i++) {
        for (int j = 0; j < colSize; j++) {
            std::string element = std::to_string(scoreMatrix[i][j]);
            int padding = maxWidth - element.length();
            std::cout << std::string(padding + 1, ' ') << element;
        }
        std::cout << std::endl;
    }

    cout << "Mismatch floor for up to 2: " << (oneString.length() * matchScore) - abs(3 * matchScore) << endl;

    if( scoreMatrix[row-1][col-1] >= (oneString.length() * matchScore) - abs(3 * matchScore)){
        cout << "we have a hit\n";
    }

    else{
        cout << "We have a miss\n";
    }
    return scoreMatrix[row-1][col-1];
}