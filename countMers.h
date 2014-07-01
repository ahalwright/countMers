/*  indexMers.h
 *  See the comments at the beginninf of indexMers.cpp
 */
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <cstring>
#include <algorithm>
using namespace std;
typedef unsigned int uint;
typedef unsigned short int usint;
typedef unsigned char uchar;
typedef long long llong;
typedef unsigned long long ullong;
//typedef unsigned int merLocType;   // Makes it easy to change this type
typedef unsigned long long merLocType;   // Makes it easy to change this type
const int merLength = 16;   // The value of K, i. e., the length of mers used.
const int queryLength = 31;   // The value of L, i. e., the length of queries.  K <= L <= 2*K.
const bool saveDnaString = false;
const bool useDnaString = false;  // use the stored dnaString in genomeChromList rather than doing a second pass through references files
const bool saveToBinaryFileFlag = true;    // Save the in-memory indices to a binary file
const bool readQueryFlag = false;     // Process a query file
const bool printAllMers = false;  // print all mers along with their location in the reference files.  Should be false if merLength > 3.
const bool checkKmerLocationsFlag = false;  // Check k-mer locations using dnaStrings
const uint maxHistSize = 200;
const uint maxNumberRefFiles = 10000;  // Maximum number of reference files to process in countMers.cpp
const uint min_display_size = 1;   // mer counts are not displayed if less than this number
const uint min_histogram_display_size = 3; // histogram entries are not displayed if less than this number
const string ACGTchars = "ACGT";
const string nonACGTchars = "URYKMSWBDHVNX";    // might be included in DNA seq

struct genChromLoc{
    usint genomeIndex;
    int  chromLoc;
};


struct db_location {  // Location of a k-mer within a contig of a genome
	uint genomeIndex;   // integer index of genome which corresponds to the order genomes are read
	/* These fields should be uncommented when the program is used to find extended matches  */
	//uint contigIndex;   // integer index of the contig within the genome
	//uint contigStart;   // start
	//uint dbStart;
	//uint dbLength;
};

typedef struct extendedMatch{
	uint length;
	string matchString;
	usint genomeIndex;
	int chromLoc;
	uint contigStartPos;
}extMatch;

typedef list<extendedMatch> extMatchList;
typedef map<string, extMatchList> extMatchMapList;

const uint goldId = 0;
const uint kingdom = 1;
const uint phylum = 2;
const uint cclass = 3;
const uint order = 4;
const uint family = 5;
const uint genus = 6;
const uint species = 7;
const uint strain = 8;
const uint taxonomic_level = phylum;

const string taxLevels[] = {"goldId","kingdom","phylum","class","order","family","genus","species","strain"};

// Used to store phylogenies are read from the file "phyloCrossRef.csv"
/*
struct phylogeny{
	string strain;
	string goldId;
	string kingdom;
	string taxLevel;
	string cclass;
	string order;
	string family;
	string genus;
	string species;
};
*/
//typedef vector<string> phylogeny;


// Indexed by genome
typedef struct genomeChrom{
	int index;
	string genome;
	string chromosome;
	uint dnaLength;     // length of dnaString
	string dnaString;  // contains the whole dna string for the genome.  Saved only if useDnaString is true
} genomeChrome;

// Integer power
unsigned long long intPower( uint base, uint exp );

// convert unsigned int to a C++ string
string intToString( uint i );

// Convert characters A, C, G, T to integers 0, 1, 2, 3
int convChar( char c );

// Convert an integer that represents a mer back to a string
//string intToMer(ullong intToConvert, int len);

// updates currentMer to reflect character c.
//  currentMer is shifted left by 2 bits, the integer representation of c is XOR'd
//  onto current mer, and then currentMer is reduced back to the appopriate range
//  by ANDing with the mask.
uint nextMerInt( char c );

// Given a K-mer (where K=merLength as set in countMers.h), returns a vector of locations of that K-mer.
// A location is a pair where the first element of the pair is the index of a genome, and the second
//     element of the pair is the location in the dna string (i. e., chromosome) of that genome.
//
vector <struct genChromLoc> lookup_K_mer( ullong mer, vector<struct genomeChrom> genomeChromList, 
	merLocType* merLocations, usint* genomeIndices, int* chromLocations, vector<struct genChromLoc>& 
	gchromlocList );

// Given a L-mer (where K<=L<=2*k, and K=merLength as set in countMers.h), returns a vector of locations of that L-mer.
// A location is a pair where the first element of the pair is the index of a genome, and the second
//     element of the pair is the location in the dna string (i. e., chromosome) of that genome.
vector <struct genChromLoc> lookup_L_mer( int L, ullong mer, vector<struct genomeChrom> genomeChromList, 
	merLocType* merLocations, usint* genomeIndices, int* chromLocations, 
	vector<struct genChromLoc>& gchromlocList );

void queryMer( ullong queryMer, string queryDescriptor, uint contigStartPos, extMatchMapList& extendeMatchMapList );

void processContigQuery( string& contigStr, uint contigStartPos, string& queryDescriptor, extMatchMapList& extendeMatchMapList );

void processContig( string& contigStr, string queryDescriptor, uint dnaLength );

void splitContigQuery( string contigStr, string queryDescriptor, extMatchMapList& extendeMatchMapList );

void splitContig( string contigStr, uint genomeCount, bool countMersFlag, uint dnaLength );

void processGenomeFromMemory( uint genomeCount, string queryDescriptor, extMatchMapList& extendeMatchMapList );

// read a fasta file as a source of queries
void readQueryFile(string fileName, extMatchMapList& extendeMatchMapList );

// For each k-mer in referenced by genomeIndices and chromLocations, check that the k-mer actually
//    occurs in the genome at the specified location.
// If the printAllMers flag is set in countMers.h, then print each such k-mer.
void checkKmerLocations( merLocType* merLocations, vector<struct genomeChrom> genomeChromList, 
	short unsigned int* genomeIndices, int* chromLocations ); 

// Reads the reference file and builds the in-memory database of mers in this file.
//void readReference(string fileName, int& dbseq, uint genomeCount );
void readReference(string fileName, uint genomeCount, vector<genomeChrom>& genomeChromList, bool countMersFlag  );

void readListOfReferenceFiles( string fileName, vector<genomeChrom>& genomeChromList, bool countMersFlag );

// Converts mer counts in merLocations to a cummulative distribution
ullong convertToCummulative( merLocType* merLocations );

// Prunes merLocations array to those elements within a certain size range
//void pruneMerArray( vector<db_location>* merLocations );

// Displays all non-empty locations in the mer array
void outputMerArray( merLocType* merLocations );

// Displays number of mers of size up to maxHistSize
void outputMerSizes( merLocType* merLocations );

// Check for zeros in cummulative merLocations array
void checkMerArray( merLocType* merLocations );

// Check differences in cummulative merLocations array
void checkCummMerArray( merLocType* merLocations );

void saveToBinaryFile( string binFileName, ullong countMers, merLocType* merLocations,
        vector<struct genomeChrom> genomeChromList, short unsigned int* genomeIndices, int* chromLocations );

void readBinaryFile( string binaryFileName, ullong& countMers );

// Displays a histogram of mer frequencies
//void merHistogram( vector<db_location>* merLocations );

