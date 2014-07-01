/*  
 *  Utilities for encodeFasta.cpp
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_set>
#include "util.h"
using namespace std;

static char complementArray[256];
static uint ones16[65536];
char toInt[256];
static const ullong m1  = 0x5555555555555555; //binary: 0101...
static const ullong m2  = 0x3333333333333333; //binary: 00110011..
static const ullong m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
static const ullong m8  = 0x00ff00ff00ff00ff; //binary:  8 zeros,  8 ones ...
static const ullong m16 = 0x0000ffff0000ffff; //binary: 16 zeros, 16 ones ...
static const ullong m32 = 0x00000000ffffffff; //binary: 32 zeros, 32 ones
static const ullong hff = 0xffffffffffffffff; //binary: all ones
static const ullong h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...

//extern char complementArray[];
//extern char toInt[];


unsigned long long intPower( uint base, uint exp ){
  unsigned long long result = 1;
  for( uint i = 0; i < exp; i++ ){
	 result *= base;
  }
  return result;
}

// Returns the number of one bits in a ullong (very efficiently)
// see http://en.wikipedia.org/wiki/Hamming_weight
uint oneCount(ullong x) {
	 x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
	 x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits 
	 x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits 
	 return (x * h01)>>56;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
}

// Initializes the array ones16 which caches the number of one bits in a short int
void ones16init(){
	for( uint i = 0; i < 65536; i++ ){
		 ones16[i] = oneCountByShift( (ullong) i );
	}
}

// Returns the number of one bits in a ullong by using 16 bit lookup table
uint oneCountLookup( ullong n){
  uint oneCount = 0;
  oneCount = ones16[ n & 0xFFFF ];
  n >>= 16;
  oneCount += ones16[ n & 0xFFFF ];
  n >>= 16;
  oneCount += ones16[ n & 0xFFFF ];
  n >>= 16;
  oneCount += ones16[ n & 0xFFFF ];
  return oneCount;
}

// Returns the number of one bits in a ullong (less efficiently)
uint oneCountByShift( ullong n ){
  uint oneCount = 0;
  for( uint i = 0; i < 64; i++ ){
	 oneCount += n&1;
	 n >>= 1;
  }
  return oneCount;
}

// returns a binary string representation of n
string printBinary( ullong n ){
  ostringstream convert;   // stream used for the conversion
  for(uint i = 0; i < 64; i++ ){
	 convert << ((n & 0x8000000000000000) >> 63);
	 n <<= 1;
  }
  return convert.str();
}

const string currentDateTime() {
	 time_t     now = time(0);
	 struct tm  tstruct;
	 char       buf[80];
	 tstruct = *localtime(&now);
	 strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
	 return buf;
}

const string currentDate() {
	 time_t     now = time(0);
	 struct tm  tstruct;
	 char       buf[80];
	 tstruct = *localtime(&now);
	 strftime(buf, sizeof(buf), "%Y-%m-%d", &tstruct);
	 return buf;
}

// Hamming distance between two strings by the "obvious" algorithm
uint hammingString( string s1, string s2 ){
  if( s1.length() != s2.length() ){
	 cerr << "strings s1 and s2 have different length in function hammingString\n";
  }
  uint count = 0;
  for( uint i = 0; i<s1.length(); i++ ){
	 count += (s1[i] != s2[i]) ? 1 : 0;
  }
  return count;
}

// convert unsigned int to a C++ string
string intToString( uint i ){
	 ostringstream stm;
	 stm << i;
	 // return resulting string
	 return  stm.str();
}

// convert a C++ string to an uint
uint stringToUint( string s){
	istringstream stm(s);
	uint result;
	stm >> result;
	return result;
}

// Initializes array complementArray[] which is used to convert a base to its complement
void initComplementArray( ){
	for( int i = 0; i < 256; i++ ){
		complementArray[i] = -1;
	}
	complementArray['A'] = 'T';
	complementArray['T'] = 'A';
	complementArray['C'] = 'G';
	complementArray['G'] = 'C';
	complementArray['a'] = 'T';
	complementArray['t'] = 'A';
	complementArray['c'] = 'G';
	complementArray['g'] = 'C';
}

// Returns the reverse complement of a DNA string (assumed to be upper case)
string reverseComplement( string dnaString ){
  string result("");
  result.reserve(dnaString.size());
  char c_ptr[2] = "X";  // C string of length 1 used to store result of convering a character
  for( int i = dnaString.size()-1; i >= 0; i-- ){
		c_ptr[0] = complementArray[(uint)dnaString[i]];
		result += c_ptr;
  }
  return string(result);
}

// Returns the reverse complement of a DNA string (assumed to be upper case)
string reverseCompl( string dnaString ){
	stringstream ss;
	for( int i = dnaString.size()-1; i >= 0; i-- ){
		ss << complementArray[(uint)dnaString[i]];
	}
	return ss.str();
}

// Set up the array toInt which is used to convert characters to integers
void toIntSetup(){
	for( int i = 0; i < 256; i++ ){
		toInt[i] = -1;
	}
	toInt[(int)'A'] = 0;
	toInt[(int)'C'] = 1;
	toInt[(int)'G'] = 2;
	toInt[(int)'T'] = 3;
}

// Convert an integer that represents a mer back to a string
string intToMer(ullong intToConvert, int len){
	char conv[] = "ACGT";
	string strResult = "";
	for(int i = 0; i < len; i++){
		strResult = conv[intToConvert & 3]+strResult;
		intToConvert = intToConvert>>2;
	}
	return strResult;
}


// Extract the genus and species from the name field of a fasta record
// The genus is the first word.
// The species is the second word unless the second word is "sp.", in which case it is the third word
vector<string> extractGenusSpecies( vector<string> result, string name ){
  uint GenusEnd = name.find(" ",0);
  //cout << "GenusEnd: " << GenusEnd << endl;
  result.push_back( name.substr(0,GenusEnd));
  uint SpeciesEnd = name.find(" ",GenusEnd+1);
  //cout << "SpeciesEnd: " << SpeciesEnd << endl;
  string species = name.substr(GenusEnd+1,SpeciesEnd-1-GenusEnd);
  if( species == "sp." ){
		GenusEnd = SpeciesEnd;
		SpeciesEnd = name.find(" ",GenusEnd+1);
		species = name.substr(GenusEnd+1,SpeciesEnd-1-GenusEnd);
  }
  /*
  bool isLowerTest = true;
  for( uint i = 0; (i<species.length()) && islower(species.c_str()[i]); i++  ){
	 if( i < species.length() ){
			isLowerTest = false;
	 }
  }
  */
  result.push_back( species );
  return result;
}

// Select M random elements from the unsigned integers 0..N-1.  
// If M > N, all integers in the range are selected.
// The unordered set S contains the selected integers.
// See http://stackoverflow.com/questions/2394246/algorithm-to-select-a-single-random-combination-of-values
// See Communications of the ACM, September 1987, Volume 30, Number 9 in Programming Pearls column of Jon Bentley
void selectRandomSubset( const uint M, const uint N, unordered_set<uint>& S ){
	S.clear();
	uint U = (M <= N) ? N-M : 0;  // U is the max of 0 and N-M
	for( uint j = U; j< N; j++ ){
		 uint T = j > 0 ? (rand() % j) : 0; // Avoid mod by zero
		 if( S.count(T) == 0 ){
			  S.insert(T);
		 }
		 else{
			  S.insert(j);
		 }
	}
}

vector<string> &split(const string &s, char delim, vector<string> &elems) {
	 stringstream ss(s);
	 string item;
	 while (getline(ss, item, delim)) {
			 elems.push_back(item);
	 }
	 return elems;
}

vector<string> split(const string &s, char delim) {
	 vector<string> elems;
	 split(s, delim, elems);
	 return elems;
}

/*
int main(){
  ifstream IN;
  string fileName = "cctgcagg.fasta";
  IN.open(fileName.c_str());
  vector<string> rrr;
  string hdrStr;
  string name;
  while( IN.good() ){
		getline(IN,hdrStr);
		if( hdrStr.length() < 5 ){
			 break;
		}
		uint refEnd = hdrStr.find("|",4);
		string gi = hdrStr.substr(4,refEnd-4);
		uint chromBegin = hdrStr.find("chromNum|",refEnd);
		uint nameEnd = hdrStr.find("|",chromBegin+11);
		name = hdrStr.substr(chromBegin+11,nameEnd-chromBegin-11);
		rrr.clear();
		rrr = extractGenusSpecies( rrr, name );
		cout << "Genus: " << rrr[0] << endl;
		cout << "Species: " << rrr[1] << endl;
		if( !IN.good() ){ break; };
		getline(IN,name);
	}
	return 0;
}
int main(){
  cout << currentDate() << endl;
  return 0;
}
*/

