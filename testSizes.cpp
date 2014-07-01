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

#include <stdio.h>
#include "countMers.h"
#include "util.h"
using namespace std;

// merLength is now set in countMers.h
//const int merLength = 14;   // The value of k, i. e., the length of mers used.

// Unfortunately, this program uses way too many global variables.  They are listed here.
unsigned long long MER_ARRAY_LENGTH;  //  MER_ARRAY_LENGTH is set to 4^merLength. 
ullong mask; // mask is set to 2^merLength-1.
ullong queryMask; // queryMask is set to 2^queryLength-1.
merLocType* merLocations;  // This is the in-memory index of k-mers
short unsigned int* genomeIndices;
uint* chromLocations;
vector<struct genomeChrom> genomeChromList;
//vector<struct genomeChrom> genomeChromList;
// currentMer is an integer in the range from 0 to 2^merLength-1 that stores the integer 
//   representation of the current mer in the reference sequence.
ullong currentMer=0LL;
// Used to convert the current mer to an integer via array lookup
char toInt[256];
char complementArray[256];


int main(int argc, char* argv[]){
	MER_ARRAY_LENGTH = intPower(4,merLength);
	mask = MER_ARRAY_LENGTH-1;
	cout << dec << "merLength: "<<merLength<<endl;
	cout << dec << "sizeof merLocType: " << sizeof(merLocType) << endl;
	cout << "mask: " << hex << mask << endl;
	//cout << dec;
	cout << dec << "array length dec: " << MER_ARRAY_LENGTH << endl;
	cout << hex << "array length hex: " << MER_ARRAY_LENGTH << endl;
	try{
		merLocations = new merLocType[MER_ARRAY_LENGTH];
	}
	catch( bad_alloc ){
		cerr << "could not allocate merLocations array\n";
	}
	cout << "merLocations:      " << hex << merLocations << endl;
	cout << "merLocations[MAL]: " << hex << &merLocations[MER_ARRAY_LENGTH-1] << endl;
	cout << "diff:              " << hex << &merLocations[MER_ARRAY_LENGTH-1] -&merLocations[0] << endl;
	cout << "mask:              " << hex << mask << endl;
	long long diff = &merLocations[MER_ARRAY_LENGTH-1] -merLocations;
	printf("diff:               %#llx\n",diff);
	printf("mask:               %#llx\n",mask);
	for( ullong i = 0; i < MER_ARRAY_LENGTH; i++ ){
		merLocations[i] = i+5;
	}
	//bool countMersFlag = true;    // Phase 1 is to count the k-mers in the genome 
	//readListOfReferenceFiles(referenceFileList, genomeChromList, countMersFlag );
	//outputMerArray( merLocations );
	//outputMerSizes( merLocations );
	//uint countMers = convertToCummulative(merLocations);
	uint countMers = 972485348;
	cout << "countMers: " << dec << countMers << "  " << hex << countMers << endl;
	checkCummMerArray( merLocations );
	checkMerArrayForZeros( merLocations );
	cout << "after checkCumm countMers: "<< countMers << endl;
	checkMerArrayForZeros( merLocations );
	try{
		genomeIndices = new usint[countMers];
		chromLocations = new uint[countMers];
	}
	catch( bad_alloc ){
		cerr << "could not allocate genomeIndices or chromLocations array\n";
	}
	cout << "mask:                " << hex << mask << endl;
	cout << "genomeIndices[0]:    " << hex << &genomeIndices[0] << endl;
	cout << "genomeIndices[MAL]:  " << hex << &genomeIndices[MER_ARRAY_LENGTH-1] << endl;
	diff = &genomeIndices[MER_ARRAY_LENGTH-1] -genomeIndices;
	cout << "diff:                " << hex << (&genomeIndices[MER_ARRAY_LENGTH-1]-genomeIndices) << endl;
	cout << "chromLocations:      " << hex << chromLocations << endl;
	cout << "chromLocations[MAL]: " << hex << &chromLocations[MER_ARRAY_LENGTH-1] << endl;
	diff = &chromLocations[MER_ARRAY_LENGTH-1]-chromLocations;
	cout << "diff:                " << hex << diff << endl;
	cout << dec;
	int gicount = 0;
	for( ullong i = 0; i < countMers; i++ ){
		genomeIndices[i] = 0;
		if( (merLocations <= (ullong*)&genomeIndices[i]) && ((ullong*)&genomeIndices[i] < &merLocations[MER_ARRAY_LENGTH-1]) ){
			cout << hex << "i: " << i << "  @genomeIndices[i]: " << hex << &genomeIndices[i] << endl;
			if( gicount > 2 ){
				break;
			}
			gicount++;
		}
	}
	cout << "after genomeIndices init " << endl;
	checkMerArrayForZeros( merLocations );
	gicount = 0;
	for( ullong i = 0; i < countMers; i++ ){
		chromLocations[i] = 0;
		if( (merLocations <= (ullong*)&chromLocations[i]) && ((ullong*)&chromLocations[i] < &merLocations[MER_ARRAY_LENGTH-1]) ){
			cout << hex << "i: " << i << "  @chromLocations[i]: " << hex << &chromLocations[i] << endl;
			if( gicount > 6 ){
				break;
			}
			gicount++;
		}
	}
	cout << "after chromLocations init " << endl;
	checkMerArrayForZeros( merLocations );
	return 0;
}
