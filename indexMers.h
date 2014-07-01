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
using namespace std;
typedef unsigned int uint;
typedef unsigned short int usint;
const int merLength = 16;   // The value of k, i. e., the length of mers used.
const uint minOverlap = 22;  // The minimum overlap allowed in following predecessors and successors.
const bool useTaxonomicWeights = false;  // If false, taxLevel weights are set to 1.0 (equiv to no weights)
const bool sqrt_weights = true;  // If true, use sqrt inverse weights, if false, inverse weights
const uint minMatch = 26;
const uint noPredSuccMarker = 999999;  // Marker for no predecessor or no successor
const bool display_matches = false;
const bool displayVotes = false;
const uint minMerArraySize = 1;
const uint maxMerArraySize = 5000;
// These specify the input files
//string referenceFileList = "shortGenomes.txt";
//const string referenceFileList = "genomeFiles.txt";
//string referenceFileList = "almostAllGenomes.txt";

// The frag file list file can also be input on the command line
//string fragFileList = "fragFiles.txt";
//string fragFileList = "fastqFiles1_6.txt";
//string fragFileList = "allFastq.txt";


struct db_location {  // Location of a k-mer within a contig of a genome
  uint genomeIndex;   // integer index of genome which corresponds to the order genomes are read
  /* These fields should be uncommented when the program is used to find extended matches  */
  //uint contigIndex;   // integer index of the contig within the genome
  //uint contigStart;   // start
  //uint dbStart;
  //uint dbLength;
};


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
typedef vector<string> phylogeny;

// The structure that is used to store a k-mer match
// Note: Assumes that we know the frag since a separate struct is used for each frag
// Assume:  relOffset = fragStart - contigStart
// Note: all k-mer matches that make up a longer exact match will have same relOffset
struct match{
  int leftRelOffset;  // fragStart - contigStart
  uint genomeIndex;
  uint contigIndex;
  uint fragStart;
};

// Structure that is used to store an extended exact match between a db sequence and a frag
// Note: Assumes that we know the frag since a separate struct is used for each frag
struct extendedMatch{
  int leftRelOffset;  // fragStart - contigStart
  uint genomeIndex;
  uint contigIndex;
  uint fragStart;
  uint matchLength;
};

// Structure to store the predecessor and successor of a frag
//  and the lengths and rel offsets that correspond to them
typedef struct predsucc {
  uint pred;
  uint pred_match_length;
  int pred_rel_offset;
  uint succ;
  uint succ_match_length;
  int succ_rel_offset;
} pred_succ_struct;  

// Used for testing when the genome and offset within the genome is known for each frag
// These are properties of a frag that are stored in an array indexed by frags
struct fragLocation{
  uint genome;
  uint offset;
  uint length;
};

// Compare two matches.  Is m1 < m2 ?
bool match_lt( match m1, match m2 );

// Integer power
unsigned long long intPower( uint base, uint exp );

// convert unsigned int to a C++ string
string intToString( uint i );

// Set up the array toInt which is used to convert characters to integers
void toIntSetup();

// Convert characters A, C, G, T to integers 0, 1, 2, 3
int convChar( char c );

// Convert an integer that represents a mer back to a string
string intToMer(int intToConvert, int len);

// updates currentMer to reflect character c.
//  currentMer is shifted left by 2 bits, the integer representation of c is XOR'd
//  onto current mer, and then currentMer is reduced back to the appopriate range
//  by ANDing with the mask.
int nextMerInt( char c );

// Reads the reference file and builds the in-memory database of mers in this file.
void readReference(string fileName, int& dbseq, uint genomeCount );

void readListOfReferenceFiles( string fileName );

// Prunes merLocations array to those elements within a certain size range
void pruneMerArray( vector<db_location>* merLocations );

// Displays all non-empty locations in the mer array
void outputMerArray( uint* merLocations );

// Displays number of mers of size up to maxMerSize
void outputMerSizes( uint* merLocations );

// Displays a histogram of mer frequencies
void merHistogram( vector<db_location>* merLocations );

// Displays a header for a sequence of matches
void outputMatchHeader();

// Displays a match in tab-delimited form.
void outputMatch( string label, match mat, uint merLength, string line );

// Displays a header for a sequence of extended matches
void outputExtendedMatchHeader();

// Displays an extended match in tab-delimited form.
void outputExtendedMatch( string label, extendedMatch emat, string line );

// Displays a header for a sequence of pred_succ_structs
void outputPredSuccHeader();

// Displays a pred_succ_struct in tab-delimited form.
void outputPredSuccStruct( uint frag, pred_succ_struct ps );

void readFragFile(string fragFileName, uint merLength, uint& fragCount, uint genomeCount, map<string,double>taxLevelWeights);

void readListOfFragFiles( string fileName );

// Constructs a list of predecessors/successors of the given frag
// Uses the global array predSucc
vector<uint>* followPredSucc( uint frag );

void constructTaxLevelWeights(map<string,int>& taxLevelCount, map<string,double>& taxLevelWeights );

void countVotes( map<string,double>& taxLevelVoteFrag, map<string,double>& taxLevelVoteAll);

double computeMutualInfo( map<string,int> taxLevelCount, map<string,int> merCountOTU );

// extracts the Gold Id from a file name.  It is preceded by a period, and starts with a "G$.
string extractGoldId( string& fileName );

/*
#ifndef FRAGID
#define FRAGID
// compare two match structures based first on field relOffset, then on dbseq,
//    then on fragStart
bool match_lt( match m1, match m2 ) {
  if( (m1.leftRelOffset < m2.leftRelOffset) ){
    return true;
  }
  else if( (m1.leftRelOffset == m2.leftRelOffset) ){
     if( m1.dbseq < m2.dbseq ){
       return true;
     }
     else if( m1.dbseq == m2.dbseq ){
        if( m1.fragStart < m2.fragStart ){
           return true;
        }
        else{
           return false;
        }
     }
     else{
       return false;
     }
  }
  else{
    return false;
  }
}

const uint noPredSuccMarker = 999999;  // Marker for no predecessor or no successor
#endif
*/
