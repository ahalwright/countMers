/* This program creates an in-memory count of all of the k-mers in one collection of DNA
 * sequences.   These are called the reference sequences.  These DNA seqnences are read in 
 * from files listed in the file whose name is stored in the variable referenceFileList.
 * This is set below.  It can be easily adapted to create an in-memory index which 
 * refers back to the reference sequences.  
 *
 * There is an associated header file "countMers.h" and a file of utility functions "util.cpp".
 *
 * Note that each file of reference sequences (which is "usually" a genome file) may contain
 * multiple contigs.  Each contig is proceeded by a line starting with a ">$ character.
 * In addition, sequences may contain the N character in addition to the A, C, G, T characters.
 * The N character indicates that the sequencer could not determine the base at that position.
 * Since we do not know how to interpret k-mers containing the N character, we will assume
 * that any sequence of consecutive N characters separates the sequence into contigs.
 * In other words, the contigs of a reference files can be further divided into subcontigs
 * (which are just called contigs in this program) by consecutive sequences of N characters.
 *
 *  Author:  Alden Wright, University of Montana CS department, adapted from indexMers.cpp
 *  in Nov. 2013.
 *
 */

#include "countMers.h"
#include "util.h"
using namespace std;

// merLength is now set in countMers.h
//const int merLength = 14;   // The value of k, i. e., the length of mers used.

// Unfortunately, this program uses way too many global variables.  They are listed here.
unsigned long long MER_ARRAY_LENGTH;  //  MER_ARRAY_LENGTH is set to 4^merLength. 
merLocType mask; // mask is set to 2^merLength-1.
ullong queryMask; // queryMask is set to 2^queryLength-1.
merLocType* merLocations;  // This is the in-memory index of k-mers
short unsigned int* genomeIndices;
int* chromLocations;
ullong countMers;
vector<struct genomeChrom> genomeChromList;
//vector<struct genomeChrom> genomeChromList;
// currentMer is an integer in the range from 0 to 2^merLength-1 that stores the integer 
//   representation of the current mer in the reference sequence.
ullong currentMer=0LL;
// Used to convert the current mer to an integer via array lookup
//char toInt[256];
extern char toInt[];
char complementArray[256];



int main(int argc, char* argv[]){
	//extMatchList extendedMatchList;
	extMatchMapList extendedMatchMapList;
	//map<string, extendedMatch> extendedMatchList;
	// These specify the input files
	string referenceFileList = "genomeFiles.txt";  // The default if no command line argument is given
	//vector<genomeChrome> genomeChromList;
	string binaryFileName = "countMers.bin";
	string fastaFileName;
	toIntSetup();
	initComplementArray();
	if( argc > 2 ){
		if( saveToBinaryFileFlag ){
			binaryFileName = argv[2];
		}
		else if( readQueryFlag ){
			fastaFileName = argv[2];
		}
	}
	if( argc > 1 ){
		referenceFileList = argv[1];
	}
	else{
		if( saveToBinaryFileFlag ){
			cerr << "USAGE:  " << argv[0] << " <file name with list of genome files> <binary file name>" << endl;
		}
		else{
			cerr << "USAGE:  " << argv[0] << " <file name of file with list of genome files> " << endl;
		}
		exit(1);
	}
	//constructPhylogenyMap();
	MER_ARRAY_LENGTH = intPower(4,merLength);
	if( merLength == 16 && sizeof(merLocType) == 4 ){  // special case when the mask is the all-ones word
		mask = ~0;
	}
	else{
		mask = MER_ARRAY_LENGTH-1;
	}
	queryMask = intPower(4,queryLength)-1;
	//mask = MER_ARRAY_LENGTH-1;
	cout << dec << "merLength: "<<merLength<<endl;
	cout << dec << "sizeof merLocType: " << sizeof(merLocType) << endl;
	cout << dec << "queryLength: "<<queryLength<<endl;
	cout << "mask: " << hex << mask << endl;
	cout << "qeuryMask: " << hex << queryMask << endl;
	cout << "saveToBinaryFileFlag: " << saveToBinaryFileFlag << endl;
	cout << "readQueryFlag: " << readQueryFlag << endl;
	cout << "display size cutoff: " << dec <<  min_display_size << endl;
	//cout << dec;
	cout << dec << "array length dec: " << MER_ARRAY_LENGTH << endl;
	cout << hex << "array length hex: " << MER_ARRAY_LENGTH << endl;
	try{
		merLocations = new merLocType[MER_ARRAY_LENGTH];
	}
	catch( bad_alloc ){
		cerr << "could not allocate merLocations array\n";
	}
	cout << "merLocations:        " << hex << merLocations << endl;
	cout << "merLocations[0]:     " << hex << &merLocations[0] << endl;
	cout << "merLocations[MAL]::  " << hex << &merLocations[MER_ARRAY_LENGTH-1] << endl;
	for( ullong i = 0; i < MER_ARRAY_LENGTH; i++ ){
		merLocations[i] = 0;
	}
	bool countMersFlag = true;    // Phase 1 is to count the k-mers in the genome 
	readListOfReferenceFiles(referenceFileList, genomeChromList, countMersFlag );
	//outputMerArray( merLocations );
	outputMerSizes( merLocations );
	ullong sum = 0;
	for( ullong i = 0; i < MER_ARRAY_LENGTH; i++ ){
		sum += merLocations[i];
	}
	cout << "sum of merLocations: "<<sum<<endl;
	countMers = convertToCummulative(merLocations);
	cout << "countMers: " << dec << countMers << "  " << hex << countMers << endl;
	//outputMerArray( merLocations );
	//checkCummMerArray( merLocations );
	//cout << "after checkCumm countMers: "<< countMers << endl;
	//checkMerArray( merLocations );
	try{
		genomeIndices = new usint[countMers];
		chromLocations = new int[countMers];
	}
	catch( bad_alloc ){
		cerr << "could not allocate genomeIndices or chromLocations array\n";
	}
	/*
	cout << "genomeIndices: " << hex << genomeIndices << endl;
	cout << "genomeIndices[0]:     " << hex << &genomeIndices[0] << endl;
	cout << "genomeIndices[MAL]::  " << hex << &genomeIndices[MER_ARRAY_LENGTH-1] << endl;
	cout << "chromLocations: " << hex << chromLocations << endl;
	cout << "chromLocations[0]:     " << hex << &chromLocations[0] << endl;
	cout << "chromLocations[MAL]::  " << hex << &chromLocations[MER_ARRAY_LENGTH-1] << endl;
	cout << dec;
	*/
	for( ullong i = 0; i < countMers; i++ ){
		genomeIndices[i] = 0;
	}
	//checkMerArray( merLocations );
	for( ullong i = 0; i < countMers; i++ ){
		chromLocations[i] = 0;
	}
	//checkCummMerArray( merLocations );
	//checkMerArray( merLocations );
	for( vector<struct genomeChrom>::iterator gg = genomeChromList.begin(); gg != genomeChromList.end(); gg++ ){
		cout << dec << gg->index << "  genome: " << gg->genome << "  chrom: " << gg->chromosome << endl;
	}
	//checkCummMerArray( merLocations );
	//checkMerArray( merLocations );
	countMersFlag = false;    // Phase 2 is to build index from k-mers to their location in genome/chromosome.
	readListOfReferenceFiles(referenceFileList, genomeChromList, countMersFlag );
	//outputMerSizes( merLocations );  // Don't use this here.
	//outputMerArray( merLocations );
	/*
	for( vector<struct genomeChrom>::iterator gg = genomeChromList.begin(); gg != genomeChromList.end(); gg++ ){
		cout << gg->index << "  genome: " << gg->genome << "  chrom: " << gg->chromosome << " : "<< gg->dnaString << endl;
	}
	for( uint j = 0; j < countMers; j++ ){
		cout << j << ": " << genomeIndices[j] << "  " << chromLocations[j] << endl;
	}
   */
	//bool kmerOutputFlag = false;
	if( checkKmerLocationsFlag && saveDnaString ){
		checkKmerLocations( merLocations, genomeChromList, genomeIndices, chromLocations );
	}
	//checkMerArray( merLocations );
	
	if( saveToBinaryFileFlag ){
		saveToBinaryFile( binaryFileName, countMers, merLocations, genomeChromList, genomeIndices, chromLocations );	
	}
	if( readQueryFlag && !saveToBinaryFileFlag ){
		//readQueryFile(fastaFileName,extendedMatchMapList);
	}
	return 0;
}
