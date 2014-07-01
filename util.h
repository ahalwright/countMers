#include <vector>
#include <string>
#include <unordered_set>
using namespace std;

typedef unsigned int uint;
typedef unsigned long long ullong;

// Returns the number of one bits in a ullong 
uint oneCount( ullong n);

// Initializes the array ones16 which caches the number of one bits in a short int
void ones16init();

// Initializes complementArray[]
void initComplementArray( );

// Returns the number of one bits in a ullong (less efficiently)
uint oneCountByShift( ullong n );

// Hamming distance between two strings by the "obvious" algorithm
uint hammingString( string s1, string s2 );

// string representation of current date and time
const string currentDateTime();

// string representation of current date 
const string currentDate();

// Integer power
ullong intPower( uint base, uint exp );

// Returns the reverse complement of a DNA string (assumed to be upper case)
string reverseComplement( string dnaString );

// Returns the reverse complement of a DNA string (assumed to be upper case)
string reverseCompl( string dnaString );

// ullong to binary string representation
string printBinary( ullong n );

// convert unsigned int to a C++ string
string intToString( uint i );

// Convert an integer that represents a mer back to a string
string intToMer(ullong intToConvert, int len);

// Set up the array toInt which is used to convert characters to integers
void toIntSetup();

// Convert characters A, C, G, T to integers 0, 1, 2, 3
int convChar( char c );

// Select M random elements from the integers 0..N-1.
// Result is the unordered set S (which is cleared on entry).
// See http://stackoverflow.com/questions/2394246/algorithm-to-select-a-single-random-combination-of-values
// See Communications of the ACM, September 1987, Volume 30, Number 9 in Programming Pearls column of Jon Bentley
void selectRandomSubset( const uint M, const uint N, unordered_set<uint>& S );

// Splits the string s based on the delimiter delim and adds the results to vector elems
// Taken from http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c
vector<string> &split(const string &s, char delim, vector<string> &elems);

// Splits the string s based on the delimiter delim returns string of results
vector<string> split(const string &s, char delim);
