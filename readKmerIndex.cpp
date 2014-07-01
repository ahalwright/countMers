#include "countMers.h"
#include "util.h"
using namespace std;

merLocType* merLocations;
usint* genomeIndices;
int* chromLocations;
vector<struct genomeChrom> genomeChromList;
ullong MER_ARRAY_LENGTH;  // set to 4^merLength
ullong mask; // mask is set to 2^merLength-1.
ullong queryMask; // queryMask is set to 2^queryLength-1.
// currentMer is an integer in the range from 0 to 2^merLength-1 that stores the integer 
//    representation of the current mer in the sequence read from the query file.
ullong currentMer = 0;
ullong countMers = 0;

extern char toInt[256];


int main( int argc, char* argv[] ){
	extMatchMapList extendeMatchMapList;
	//map<string, extendedMatch> extendedMatchList;
	//merLocType* merLocations;
	//usint* genomeIndices;
	//int* chromLocations;
	//vector<struct genomeChrom> genomeChromList;
	//unsigned long long MER_ARRAY_LENGTH;
	toIntSetup();
	initComplementArray();
	string binaryFileName = "countMers.bin";
	string fastaFileName = "query.fasta";
	if( argc > 2 ){
		fastaFileName = argv[2];
	}
	if( argc > 1 ){
		binaryFileName = argv[1];
	}
	else{
		cerr << "USAGE:  " << argv[0] << " <binary file name> <fasta file name>" << endl; 
	}
	//readBinaryFile( binaryFileName, countMers, merLocations, genomeChromList, genomeIndices, chromLocations );
	readBinaryFile( binaryFileName, countMers );
	cout << "countMers: " << countMers << endl;
	checkMerArray( merLocations );
	//checkCummMerArray( merLocations );
	/*
	vector<int> queryList = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
	vector<int>::iterator i;
	for( i = queryList.begin(); i != queryList.end(); i++ ){
		uint genomeIndex = genomeIndices[*i];
		cout << intToMer(*i,merLength) << ": "<<genomeIndex << " : " << chromLocations[*i] << endl;
	}
	*/
	/*
	for( ullong i =0; i< MER_ARRAY_LENGTH; i++ ){
		vector<struct genChromLoc> gchromlocList;
		if( (i > 0 && merLocations[i-1] != merLocations[i]) || (i==0 && merLocations[0]>0) ){
			cout << "genomes and chroms for mer: " << intToMer(i,merLength) << endl;
			gchromlocList = lookup_K_mer(i,genomeChromList,merLocations,genomeIndices,chromLocations,gchromlocList);
			vector<struct genChromLoc>::iterator j;
			for( j = gchromlocList.begin(); j != gchromlocList.end(); j++ ){
				uint genomeInd = j->genomeIndex;
				cout << genomeChromList[genomeInd].genome <<" " << genomeChromList[genomeInd].chromosome<< ":" << j->chromLoc << endl;
			}
		}
	}
	*/
	/*
	ullong mer = 0;
	uint L=merLength+2;
	cout << "genomes and chroms for mer: " << intToMer(mer,L) << endl;
	vector<struct genChromLoc> gchromlocList;
	gchromlocList = lookup_L_mer(L,mer,genomeChromList,merLocations,genomeIndices,chromLocations,gchromlocList);
	vector<struct genChromLoc>::iterator j;
	for( j = gchromlocList.begin(); j != gchromlocList.end(); j++ ){
		uint genomeInd = j->genomeIndex;
		cout << genomeChromList[genomeInd].genome <<" " << genomeChromList[genomeInd].chromosome<< ":" << j->chromLoc << endl;
	}
	*/
	//cout << "mer: " << intToMer(1639927868,merLength)<< "  merLocations[1639927868]: " << merLocations[1639927868] << endl;
	readQueryFile(fastaFileName,extendeMatchMapList);
	if( checkKmerLocationsFlag && saveDnaString ){
		checkKmerLocations( merLocations, genomeChromList, genomeIndices, chromLocations );
	}
	extMatchMapList::iterator i;
	extMatchList::iterator j;
	for( i = extendeMatchMapList.begin(); i!=extendeMatchMapList.end(); i++ ){
		cout << "extended match list for query descriptor: " << i->first << endl;
		for( j = i->second.begin(); j != i->second.end(); j++ ){
			cout << "matchStr: " << j->matchString << "  contigStartPos: " << j->contigStartPos << endl;
			cout <<genomeChromList[j->genomeIndex].genome << " : " << genomeChromList[j->genomeIndex].chromosome << endl;
		}
	}
	return 0;
}
