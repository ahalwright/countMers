#include "countMers.h"
#include "util.h"
using namespace std;

extern merLocType* merLocations;
extern usint* genomeIndices;
extern int* chromLocations;
extern vector<struct genomeChrom> genomeChromList;
extern extMatchMapList extendedMatchMapList;
extern ullong mask; // mask is set to 2^merLength-1.
extern ullong queryMask; // queryMask is set to 2^queryLength-1.
extern ullong currentMer;
extern char toInt[256];
extern ullong MER_ARRAY_LENGTH;  // set to 4^merLength
extern ullong countMers;

// Given a K-mer (where K=merLength as set in countMers.h), returns a vector of locations of that K-mer.
// A location is a pair where the first element of the pair is the index of a genome, and the second
//     element of the pair is the location in the dna string (i. e., chromosome) of that genome.
vector <struct genChromLoc> lookup_K_mer( ullong mer, vector<struct genomeChrom> genomeChromList, merLocType* merLocations, 
		usint* genomeIndices, int* chromLocations, vector<struct genChromLoc>& gchromlocList ){
	merLocType startIndex;
	/*
	cout << "Entering lookup_K_mer with mer: " << mer << " : " << intToMer(mer,merLength) << endl;
	cout << "merLocations[mer]: " << merLocations[mer] << endl;
	*/
	gchromlocList.clear();
	if( mer == 0 ){
		startIndex = 0;
	}
	else{
		startIndex = merLocations[mer-1];
	}
	uint count = 0;
	for( merLocType j = startIndex; j < merLocations[mer]; j++ ){
		//cout << "j: " << j << "  genomeIndex: " << genomeIndex << "  chromLoc: " << chromLoc << endl;
		struct genChromLoc gc = {genomeIndices[j],chromLocations[j]};
		gchromlocList.push_back(gc);
		count++;
		/* remove when cleaning up code
		if( count > 200){
			cerr << "exiting \n";
			exit(1);
		}
		*/
	}
	return gchromlocList;
}

// Given a L-mer (where K<=L<=2*k, and K=merLength as set in countMers.h), returns a vector of locations of that L-mer.
// A location is a pair where the first element of the pair is the index of a genome, and the second
//     element of the pair is the location in the dna string (i. e., chromosome) of that genome.
vector <struct genChromLoc> lookup_L_mer( int L, ullong mer, vector<struct genomeChrom> genomeChromList, merLocType* merLocations, 
		usint* genomeIndices, int* chromLocations, vector<struct genChromLoc>& gchromlocList ){
	if( L > 2*merLength ){
		cerr << "L = "<<L<<" cannot be more than double merLength=K= "<<merLength<<endl;
		cerr << "exiting"<<endl;
		exit(1);
	}
	if( L == merLength ){
		return lookup_K_mer( mer, genomeChromList, merLocations, genomeIndices,chromLocations,gchromlocList);
	}
	gchromlocList.clear();
	// Make more efficient when cleaning up code
	ullong mer1 = mer >>2*(L-merLength);
	ullong mer2 = mer % (1ULL<<(2*merLength));  // 1<<(2*merLength) is equal to intPower(4,merLength)
	//cout << "mer1: " << mer1 << "  mer2: " << mer2 << "  mer2x: " << mer2x <<endl;
	//cout << "mer1: " << intToMer(mer1,merLength) << "  mer2: " << intToMer(mer2,merLength) << endl;
	merLocType startIndex1, startIndex2;
	if( mer1 == 0 ){
		startIndex1 = 0;
	}
	else{
		startIndex1 = merLocations[mer1-1];
	}
	if( mer2 == 0 ){
		startIndex2 = 0;
	}
	else{
		startIndex2 = merLocations[mer2-1];
	}
	merLocType j1 = startIndex1;
	merLocType j2 = startIndex2;
	while( j1 < merLocations[mer1] && j2 < merLocations[mer2] ){
		usint genomeIndex1 = genomeIndices[j1];
		int chromLoc1 = chromLocations[j1];
		int chromLoc1abs = abs(chromLoc1);   // use absolute value for comparison with chromLoc2
		//cout << "gI1: " << genomeIndex1 << "  cL1: " << chromLoc1 << endl;
		usint genomeIndex2 = genomeIndices[j2];
		int chromLoc2abs = abs(chromLocations[j2]);
		//cout << "gI2: " << genomeIndex2 << "  cL2: " << chromLoc2 << endl;
		if( genomeIndex1 == genomeIndex2 && chromLoc1abs+L-merLength == chromLoc2abs ){
			struct genChromLoc gc = {genomeIndex1, chromLoc1};  // Store negative value of chromLoc1 if found in rev complement
			gchromlocList.push_back(gc);
			j1++;  j2++;
		}
		else if( genomeIndex1 < genomeIndex2 ){
			j1++;
		}
		else if( genomeIndex1 > genomeIndex2 ){
			j2++;
		}
		else if( genomeIndex1 == genomeIndex2 && chromLoc1abs+L-merLength < chromLoc2abs ){
			j1++;
		}
		else if( genomeIndex1 == genomeIndex2 && chromLoc1abs+L-merLength > chromLoc2abs ){
			j2++;
		}
	}
	return gchromlocList;
}

// uses lookup_L_mer to get a list of matches for queryMer, then prints this list of matches
void queryMer( ullong queryMer, string queryDescriptor, uint contigStartPos, extMatchMapList& extendedMatchMapList ){
	vector <struct genChromLoc> queryLocList;
	extMatchMapList::iterator j;
	extMatchList::iterator k;
	queryLocList = lookup_L_mer( queryLength, queryMer, genomeChromList,  merLocations, 
		genomeIndices, chromLocations, queryLocList );
	if( !queryLocList.empty() ){
		vector <struct genChromLoc>::iterator i;
		string merStr = intToMer( queryMer, queryLength );
		cout << "mer: " << merStr << "  len: " << queryLocList.size();
		cout << "   descriptor: " << queryDescriptor << endl;
		for( i = queryLocList.begin(); i != queryLocList.end(); i++ ){
			cout << i->genomeIndex << " : " << i->chromLoc <<" : " << contigStartPos <<  "::";
			cout <<genomeChromList[i->genomeIndex].genome << " : " << genomeChromList[i->genomeIndex].chromosome << endl;
			j = extendedMatchMapList.find(queryDescriptor);
			bool found = false;
			if( j != extendedMatchMapList.end() ){
				int m = 0;
				for( k = j->second.begin(); k!=j->second.end(); k++ ){
					cout << "matchList iteration: " << m++ << endl;
					cout << "llabs(k->chromLoc): " << llabs(k->chromLoc) << "  llabs(i->chromLoc)-1: " << (llabs(i->chromLoc)-1) << endl;
					if( k->genomeIndex == i->genomeIndex && llabs(k->chromLoc) == llabs(i->chromLoc)-1 ){
						cout << "found ext match of length: " << k->length << "  contigStartPos: " << k->contigStartPos << endl;
						found = true;
						k->length++;
						k->matchString += merStr.at(queryLength-1);
						if( k->chromLoc > 0 ){
							k->chromLoc++;
						}
						else{
							k->chromLoc--;
						}
						k->contigStartPos = contigStartPos;
						cout << "len: " << k->length << "  genomeIndex: " << k->genomeIndex << "  chromLoc: " << k->chromLoc;
						cout << "  contigStartPos: " << k->contigStartPos;
						cout << "  matchStr: " << k->matchString << endl;
					}
				}
				if( !found ){
					extMatch extendedMatch = {queryLength,merStr,i->genomeIndex,i->chromLoc,contigStartPos};
					j->second.push_back(extendedMatch);
					cout << "added record to ext Match List for " << queryDescriptor << endl;
				}
			}
			else{
				cout << "creating extended match list: merStr: " << merStr << endl;
				extMatch extendedMatch = {queryLength,merStr,i->genomeIndex,i->chromLoc,contigStartPos};
				extMatchList* extendedMatchListPtr = new extMatchList;
				extendedMatchListPtr->push_back(extendedMatch);
				extendedMatchMapList.insert(pair<string,extMatchList>(queryDescriptor,*extendedMatchListPtr));
				j = extendedMatchMapList.find(queryDescriptor);
				cout << "matchString: " << j->second.begin()->matchString << endl;
			}
		}
	}
	else{
		//cout << "query unsuccessful for mer: " << intToMer( queryMer, queryLength ) << endl;
	}
}

// updates currentMer to reflect character c.  
//    currentMer is shifted left by 2 bits, the integer representation of c is XOR'd 
//    onto current mer, and then currentMer is reduced back to the appopriate range 
//    by ANDing with the mask.
/*
ullong nextMerInt( char c ){
	currentMer = ((currentMer<<2)^toInt[(int)c])&mask;
	return currentMer;
}
*/

// Extract the K-mers from the subcontigs produced by splitContigQuery.
// queryMer is called for each K-mer.
void processContigQuery( string& contigStr, uint contigStartPos, string& queryDescriptor, extMatchMapList& extendedMatchMapList ){
	//cout << "processContigQuery contig: " << contigStr << endl;
	currentMer = 0;
	int count=0;  // count of characters processed on a line
			if( contigStr.find_first_of(nonACGTchars) != string::npos ){  // Comment out for speed
				cerr << "contig contains a non-ACGT character which is not correct!!\n";
				cerr << contigStr << endl;
				exit(0);
			}
		for( uint i = 0; i < contigStr.length(); i++ ){
			//ullong mer = nextMerInt(contigStr[i]);
			//cout << (int)toInt[(int)contigStr[i]] << endl;
			currentMer = ((currentMer<<2)^toInt[(int)contigStr[i]])&queryMask;
			ullong mer = currentMer;
			if( count >= queryLength-1 ){
				//cout << "mer: " << intToMer(mer,queryLength) << endl;
				queryMer( mer, queryDescriptor, contigStartPos+count, extendedMatchMapList );
				//cout << "mer: " << mer << endl;
				/*
					if( !countMersFlag ){
					genomeIndices[merLocations[mer]] = genomeCount;
					chromLocations[merLocations[mer]] = contigStartPos + count - merLength + 1;
				}
				merLocations[mer]++;
				*/
			}
			count++;  // increment to next character
		}
}

// Split contigStr into multiple subcontigs separated by non-ACGT characters
// This version is called by readQuery.
void splitContigQuery( string contigStr, string queryDescriptor, extMatchMapList& extendedMatchMapList ){
	//cout <<"splitContigQuery: "<< contigStr<<endl;
	size_t nonACGTpos = 0;
	size_t ACGTpos = 0;
	while( (ACGTpos = contigStr.find_first_of(ACGTchars,nonACGTpos)) != string::npos ){
		//cout << "ACGTpos: " << ACGTpos << endl;
		nonACGTpos = contigStr.find_first_of(nonACGTchars,ACGTpos+1);
		if( nonACGTpos == string::npos ){
			nonACGTpos = contigStr.length();
		}
		//cout << "nonACGTpos: " << nonACGTpos << endl;
		string subContig = contigStr.substr(ACGTpos,nonACGTpos-ACGTpos);
		//cout << "subcontig: " << subContig  << endl;
		processContigQuery( subContig, ACGTpos, queryDescriptor, extendedMatchMapList );
        //processContig( subContig, genomeCount, ACGTpos, countMersFlag );
	}
}	

// Equivalent to functionReadReference below except that it uses the dnaString stored in
//    genomeChromList rather than rereading it from the file
void processGenomeFromMemory( uint genomeCount, string queryDescriptor, extMatchMapList& extendedMatchMapList ){
	splitContigQuery( genomeChromList[genomeCount].dnaString, queryDescriptor, extendedMatchMapList );
}

// read a fasta file as a source of queries
void readQueryFile(string fileName, extMatchMapList& extendedMatchMapList ){
	// A "contig" here is part or all of the DNA string for this file which is delimited by non-ACGT characters.
	uint contigCount = 0;  // the index of the contig.  Each genome may have multiple contigs
	ifstream IN;
	IN.open(fileName.c_str());
	if(IN.is_open()){
		cout << "Query file " << fileName << " opened." << endl;
	}else{
		cout << "couldn't open query file " << fileName << endl;
		exit(0);
	}
	// clear extendedMatchMapList which will store extended matches for each query descriptor
	extendedMatchMapList.clear();
	// fastq_flag is true if file is a fragment file.
	bool fastq_flag = (fileName.substr(fileName.length()-6,6)==".fastq") ? true : false ;
 
	//~ #~ read in dummy line. For a fasta file, all following lines are part of genome.
	string line;
	cout << dec;
	cout << "processing " << (fastq_flag?"fragments":"contigs") << " from reference file: " << fileName << endl;
	getline (IN,line);// descriptor line
	if( line[0] != '>' ){
		cout << "Error: the following line is supposed to be a descriptor line which starts with a > character." << endl;
		exit(1);
	}
	string DNA_string = "";
	string contigStr = "";
	string queryDescriptor = "";
	string newQueryDescriptor = line;
	//cout << "header line: " << line << "   genomeCount: " << genomeCount << endl;
	bool contigEnd = false;
	//size_t firstNpos, lastNpos;  // first and last positions of nonACGT character
	while ( IN.good() )  {
		contigStr = "";
		do{  // each line that starts with '>' terminates the contig
			getline (IN,line);  // may be a line of genome, or contig label
			//std::transform(line.begin(), line.end(),line.begin(), ::toupper);
			//cout << "line: " << line << endl;
			contigEnd = (line.substr(0,1)==">");
			//cout << "contigEnd: " << contigEnd << endl;
			if( !contigEnd ){
				 contigStr += line;
			}
			else{
				queryDescriptor = newQueryDescriptor;
				newQueryDescriptor = line;
			}
			contigCount++;
		}
		while( IN.good() && (!contigEnd) );
		//DNA_string += contigStr;  //  The whole DNA string for this file
		//cout << "desc str: " << descStr<< endl;
		if( !IN.good() ){
			queryDescriptor = newQueryDescriptor;
		}
		splitContigQuery( contigStr, queryDescriptor, extendedMatchMapList );
	}
	cout << "processed " << contigCount << (fastq_flag?" fragments":" contigs") << " from query file: " << fileName << endl;
	//cout << "DNAstr: " << DNA_string << endl;
}

// Duplicated from countMers.cpp
// Displays all non-empty locations in the mer array
void outputMerArray( merLocType* merLocations ){
	//const int min_display_size = 20;   // Now in countMers.h
	cout << "outputMerArray"<<endl;
	//vector<db_location>::iterator it;
	int count_mers = 0;
	merLocType count_large_mers = 0;
	for( ullong i =0; i< MER_ARRAY_LENGTH; i++ ){
		merLocType size = merLocations[i];
		count_mers += size;
		if( size >= min_histogram_display_size ){
			cout << intToMer(i,merLength) << " : " << (int)merLocations[i]; 
			cout << endl;
			count_large_mers += merLocations[i];
		}
	}
	cout << "Total mers: " << count_mers << endl;
	cout << "Large mers: " << count_large_mers << endl;
}

// For each k-mer in referenced by genomeIndices and chromLocations, check that the k-mer actually
//    occurs in the genome at the specified location.
// useDnaString must be true in countMers.h.
// If the printAllMers flag is set in countMers.h, then print each such k-mer.
void checkKmerLocations( 
		merLocType* merLocations, vector<struct genomeChrom> genomeChromList, short unsigned int* genomeIndices, uint* chromLocations ){
	if( !useDnaString ){
		cerr << "checkKmerLocations error:  useDnaString must be set to true in countMers.h to call checkKmerLocations.\n";
		return;
	}
	merLocType startIndex;
	for( ullong i =0; i< MER_ARRAY_LENGTH; i++ ){
		if( i == 0 ){
			startIndex = 0;
		}
		else{
			startIndex = merLocations[i-1];
		}
		for( uint j = startIndex; j < merLocations[i]; j++ ){
			uint genomeIndex = genomeIndices[j];
			if( printAllMers ){
				cout << i<< " " << j << " " << intToMer(i,merLength) << ": " << genomeIndex<< ": ";
				cout << genomeChromList[genomeIndex].genome << ":"<<genomeChromList[genomeIndex].chromosome ;
				cout << "  "<< chromLocations[j]<<"  " 
					<< genomeChromList[genomeIndex].dnaString.substr( chromLocations[j], merLength ) << endl; 
			}
			if( intToMer(i,merLength) != genomeChromList[genomeIndex].dnaString.substr( chromLocations[j], merLength ) ){
				cerr << "Error: incorrect mer at mer count " << j << endl;
				exit(1);
			}
		}
	}
}

// Check for zeros in cummulative merLocations array
void checkMerArray( merLocType* merLocations ){
	cout << "checkMerArray with countMers = "<< countMers <<endl;
	//vector<db_location>::iterator it;
	//int count_mers = 0;
	//int count_large_mers = 0;
	uint count = 0;
	for( ullong i =1; i< MER_ARRAY_LENGTH; i++ ){
		//cout << "i: "<<i<<endl;
		if( merLocations[i] == 00 || merLocations[i] > countMers ){
			cout <<"mer: " <<  intToMer(i,merLength) << "  merLocations[mer]: " << merLocations[i]; 
			cout << endl;
			count++;
			if( count > 5 ) break;
		}
	}
}

// Check differences in cummulative merLocations array
void checkCummMerArray( merLocType* merLocations ){
	//const int min_display_size = 20;   // Now in countMers.h
	cout << "checkCummMerArray"<<endl;
	//vector<db_location>::iterator it;
	//int count_mers = 0;
	//int count_large_mers = 0;
	uint errorCount = 0;
	for( ullong i =1; i< MER_ARRAY_LENGTH; i++ ){
		//cout << "i: "<<i<<endl;
		if( merLocations[i] > countMers ){
			cout << "Error  merLocations[i] greater than countmers  ";
			cout <<"i  : " << i<< " " <<  intToMer(i,merLength) << "  merLocations["<<i<<"]: " << merLocations[i]<< endl; 
			errorCount++;
		}
		if( (merLocations[i] < merLocations[i-1]) || (merLocations[i] - merLocations[i-1] > 10000) ){
			cout <<"i-1: " <<  intToMer(i-1,merLength) << "  merLocations["<<i-1<<"]: " << merLocations[i-1]; 
			cout << endl;
			cout <<"i  : " <<  intToMer(i,merLength) << "  merLocations["<<i<<"]: " << merLocations[i]; 
			cout << endl;
			errorCount++;
		}
		if( errorCount > 10 ) break;
	}
}

//void readBinaryFile( string binaryFileName, ullong& countMers, merLocType* merLocations,
//		vector<struct genomeChrom> genomeChromList, short unsigned int* genomeIndices, uint* chromLocations ){
void readBinaryFile( string binaryFileName, ullong& countMers ){
	ifstream IN (binaryFileName, ios::in|ios::binary|ios::ate);
	if (IN.is_open()) {
		// all is good
	}
	else{
		cerr << "unable to open binary input file " << binaryFileName << endl;
		exit(1);
	}
	ullong* lengths;
	lengths = new ullong[3];
	IN.seekg (0, ios::beg);
	IN.read(reinterpret_cast<char *>(lengths),3*sizeof(ullong));
	int myMerLength = lengths[0];
	countMers = lengths[1];
	uint genomeChromListLength = lengths[2];
	if( myMerLength != merLength ){
		cerr << "merLength from countMers.h: " << merLength << endl;
		cerr << "merLength from " << binaryFileName << ": " << myMerLength << endl;
		cerr << "These do not agree, so exiting!!" << endl;
		exit(1);
	}
	cout << "reading binary file: " << binaryFileName << endl;
	cout << "merLength: " << myMerLength << endl;
	cout << dec << "sizeof merLocType: " << sizeof(merLocType) << endl;
	cout << "queryLength: " << queryLength << endl;
	cout << "genomeChromListLength: " << genomeChromListLength << endl;
	MER_ARRAY_LENGTH = intPower(4,merLength);
	ullong QUERY_MER_ARRAY_LENGTH = intPower(4,queryLength);
	queryMask = QUERY_MER_ARRAY_LENGTH-1;
	cout << "countMers: " << countMers << endl;
	cout << "genomeChromListLength: "<<genomeChromListLength<<endl;
	cout << "genomeChromList: " << endl;
	for( uint ii = 0; ii < genomeChromListLength; ii++ ){
		size_t gsize=0, csize=0;
		IN.read(reinterpret_cast<char *>(&gsize),sizeof(size_t));
		char* genomeStr = new char[gsize+1];
		IN.read(genomeStr,gsize*sizeof(char));
		genomeStr[gsize] = '\0';
		IN.read(reinterpret_cast<char *>(&csize),sizeof(size_t));
		char* chromStr = new char[csize+1];
		IN.read(chromStr,csize*sizeof(char));
		chromStr[csize] = '\0';
		uint dnaLen;
		IN.read(reinterpret_cast<char *>(&dnaLen),sizeof(uint));
		struct genomeChrom gtemp = {ii,string(genomeStr), string(chromStr),dnaLen,string("")};
		genomeChromList.push_back(gtemp);
		//cout << "ii: " << ii << ": " << genomeChromList[ii].genome<<"  "<<genomeChromList[ii].chromosome<<"  len:"<<dnaLen<<endl;
	}
	try{
		merLocations = new merLocType[MER_ARRAY_LENGTH];
	}
	catch( bad_alloc ){
		cerr << "could not allocate merLocations array\n";
	}
	cout << "MER_ARRAY_LENGTH: "<<MER_ARRAY_LENGTH << endl;
	IN.read(reinterpret_cast<char *>(merLocations),MER_ARRAY_LENGTH*sizeof(merLocType));
	//outputMerArray( merLocations );
	/*
	for( uint i = 0; i < MER_ARRAY_LENGTH; i++ ){
		cout << i << ": "<<merLocations[i] << endl;
	}
	*/
	try{
		genomeIndices = new usint[countMers];
		chromLocations = new int[countMers];
	}
	catch( bad_alloc ){
		cerr << "could not allocate genomeIndices or chromLocations array\n";
	}
	IN.read(reinterpret_cast<char *>(genomeIndices),countMers*sizeof(usint));
	IN.read(reinterpret_cast<char *>(chromLocations),countMers*sizeof(uint));
	/*
	for( uint i = 0; i < countMers; i++ ){
		cout << i << ": "<<genomeIndices[i] << "  "<< chromLocations[i] << endl;
	}
	*/
	IN.close();
	//checkCummMerArray( merLocations );
}

// These specify the input files
string referenceFileList = "genomeFiles.txt";  // The default if no command line argument is given

/*  Moved to util.cpp
// Convert an integer that represents a mer back to a string
string intToMer(int intToConvert, int len){
	char conv[] = "ACGT";
	string strResult = "";
	for(int i = 0; i < len; i++){
		strResult = conv[intToConvert & 3]+strResult;
		intToConvert = intToConvert>>2;
	} 
	return strResult;
}
*/

// updates currentMer to reflect character c.  
//  currentMer is shifted left by 2 bits, the integer representation of c is XOR'd 
//  onto current mer, and then currentMer is reduced back to the appopriate range 
//  by ANDing with the mask.
/*
uint nextMerInt( char c ){
	currentMer = ((currentMer<<2)^toInt[(int)c])&mask;
	return currentMer;
}
*/

void processContig( string& contigStr, uint genomeCount, uint contigStartPos, bool countMersFlag, uint dnaLength ){
	if( !countMersFlag ){
		//cout << dec <<"processContig: "<<genomeCount<<endl;
		//cout <<contigStr<<endl;
	}
	string contigStrRevComplement = reverseComplement( contigStr );
	int count=0;  // count of characters processed on a line
	//db_location* locptr = NULL;
	    //cout << "contig strt: " << contigStartPos<< "  contig: " << contigStr << endl;
		if( contigStr.find_first_of(nonACGTchars) != string::npos ){  // Comment out for speed
			cerr << "contig contains a non-ACGT character which is not correct!!\n";
			cerr << contigStr << endl;
			exit(0);
		}
		for( uint i = 0; i < contigStr.length(); i++ ){
			//uint mer = nextMerInt(contigStr[i]);
			currentMer = ((currentMer<<2)^toInt[(int)contigStr[i]])&mask;
			ullong mer = currentMer;
			if( count >= merLength-1 ){
			 	if( !countMersFlag ){
					//cout << "pc i: " << i << "  mer: " << intToMer(mer,merLength) << "   ml[mer]: " << merLocations[mer] << endl;
					genomeIndices[merLocations[mer]] = genomeCount;
					chromLocations[merLocations[mer]] = contigStartPos + count - merLength + 1;
				}
				merLocations[mer]++;
			}
			count++;  // increment to next character
		}
		//cout << "processContig after first loop.  count: " << count << endl;
	    //cout << "contig strt: " << contigStartPos<< "  contig: " << contigStrRevComplement << endl;
	    count = 0;
		for( uint i = 0; i < contigStrRevComplement.length(); i++ ){
			//uint mer = nextMerInt(contigStrRevComplement[i]);
			currentMer = ((currentMer<<2)^toInt[(int)contigStrRevComplement[i]])&mask;
			ullong mer = currentMer;
			if( count >= merLength-1 ){
			 	if( !countMersFlag ){
					genomeIndices[merLocations[mer]] = genomeCount;
					//chromLocations[merLocations[mer]] = contigStartPos + count - merLength + 1;
					chromLocations[merLocations[mer]] = 
						-((dnaLength - contigStartPos - contigStr.length()) +count - merLength+1);
					//cout << "pc i: " << i <<"  count: "<<count<< "  mer: " << intToMer(mer,merLength) << "   ml[mer]: " << merLocations[mer] << "  cl[ml[mer]]: "<< chromLocations[merLocations[mer]]<<endl;
				}
				merLocations[mer]++;
			}
			//cout << "count: " << count << endl;
			count++;  // increment to next character
		}
}

void splitContig( string contigStr, uint genomeCount, bool countMersFlag, uint dnaLength ){
	if( !countMersFlag ){
		//cout << dec <<"splitContig: "<<genomeCount<<endl;
		//checkMerArrayForZeros( merLocations );
		//checkCummMerArray( merLocations );
		//cout <<contigStr<<endl;
	}
	size_t nonACGTpos = 0;
	size_t ACGTpos = 0;
	while( (ACGTpos = contigStr.find_first_of(ACGTchars,nonACGTpos)) != string::npos ){
		//cout << "ACGTpos: " << ACGTpos << endl;
		nonACGTpos = contigStr.find_first_of(nonACGTchars,ACGTpos+1);
		if( nonACGTpos == string::npos ){
			nonACGTpos = contigStr.length();
		}
		//cout << "nonACGTpos: " << nonACGTpos << endl;
		string subContig = contigStr.substr(ACGTpos,nonACGTpos-ACGTpos);
		//cout << "subcontig: " << subContig  << endl;
		processContig( subContig, genomeCount, ACGTpos, countMersFlag, dnaLength );
	}
}

// Equivalent to functionReadReference below except that it uses the dnaString stored in
//    genomeChromList rather than rereading it from the file
void processGenomeFromMemory( uint genomeCount, vector<genomeChrom>& genomeChromList, bool countMersFlag ){
	splitContig( genomeChromList[genomeCount].dnaString, genomeCount, countMersFlag, genomeChromList[genomeCount].dnaLength );
}

// Reads the reference file and adds the mers of this file to the the in-memory database 
void readReference(string fileName, uint genomeCount, vector<genomeChrom>& genomeChromList, bool countMersFlag ){
	// A "contig" here is part or all of the DNA string for this file which is delimited by non-ACGT characters.
	int contig = 0;  // the index of the contig.  Each genome may have multiple contigs
	ifstream IN;
	IN.open(fileName.c_str());
	if(IN.is_open()){
		//~ all is good
	}else{
		cout << "couldn't open genome file " << fileName << endl;
		exit(0);
	}
	// fastq_flag is true if file is a fragment file.
	bool fastq_flag = (fileName.substr(fileName.length()-6,6)==".fastq") ? true : false ;
 
	//~ #~ read in dummy line. For a fasta file, all following lines are part of genome.
	string line;
	//cout << "processing " << (fastq_flag?"fragments":"contigs") << " from reference file: " << fileName << endl;
	getline (IN,line);// descriptor line
	if( line[0] != '>' ){
		cout << "Error: the following line is supposed to be a descriptor line which starts with a > character." << endl;
		exit(1);
	}
	string DNA_string = "";
	string contigStr = "";
	//cout << "header line: " << line << "   genomeCount: " << genomeCount << endl;
	bool contigEnd = false;
	//size_t firstNpos, lastNpos;  // first and last positions of nonACGT character
	while ( IN.good() )  {
		do{  // each line that starts with '>' terminates the contig
			getline (IN,line);  // may be a line of genome, or contig label
			//std::transform(line.begin(), line.end(),line.begin(), ::toupper);
			//cout << "line: " << line << endl;
			contigEnd = (line.substr(0,1)==">");
			//cout << "contigEnd: " << contigEnd << endl;
			if( !contigEnd ){
			  contigStr += line;
			}
		}
		while( IN.good() && (!contigEnd) );
		DNA_string += contigStr;  //  The whole DNA string for this file
		//cout << "contigStartPos: " << contigStartPos << endl;
		//processContig( contig, contigStr, genomeCount, contigStartPos, countMersFlag );
		//processContig is called by splitContig
		splitContig( contigStr, genomeCount, countMersFlag, DNA_string.size() );
		contig++;
	}
	cout << "processed " << contig << (fastq_flag?" fragments":" contigs") << " from reference file: " << fileName << endl;
	//cout << "DNAstr: " << DNA_string << endl;
	genomeChromList[genomeCount].dnaLength = DNA_string.size();
	if( countMersFlag && saveDnaString ){
		genomeChromList[genomeCount].dnaString = DNA_string;
	}
}

void readListOfReferenceFiles( string fileName, vector<genomeChrom>& genomeChromList, bool countMersFlag ){
	uint genomeCount = 0;
	//uint pathOffset = 0;
	uint pathOffset = 6;   // AllGenomeFiles.txt on pardosa
	cout << "pathOffset: " << pathOffset << endl;
	ifstream IN;
	IN.open(fileName.c_str());
	if(IN.is_open()){
		//~ all is good
		cout << "File " << fileName << " of list of reference files opened."<<endl;
	}else{
		cout << "couldn't open file " << fileName << " with list of reference file names\n";
		exit(0);
	}
	string line;
	//int dbseq = 0;
	while ( IN.good() )  {
		string GoldId = "";
		getline (IN,line);
		//cout << "fline: " << line << endl;
		if( (line.length() > 0) && (line.substr(0,2) !=  "//") ){
			//string GoldId = extractGoldId( line );
			//cout << "GoldId: " << GoldId << endl;
			//genomeToGoldId[genomeCount] = GoldId;
			if( countMersFlag ){
				string genome, chromosome;
				vector<string> pathFields = split( line, '/' );
				//cout << "pathFields.size: " << pathFields.size() << endl;
				if( pathFields.size() > 2 ){
					genome = pathFields[pathOffset];
					chromosome = pathFields[pathOffset+1];
					cout << "genome: " << genome;
					cout <<	"    chrom: " << chromosome << endl;
				}
				else if( pathFields.size() == 2 ){
				  genome = pathFields[0];
				  chromosome = pathFields[1];
				}
				else{
				  genome = line;
				  chromosome = "";
				}
				genomeChrom tmp = { genomeCount, genome, chromosome, 0, "" };
				genomeChromList.push_back( tmp );
			}
			if( !countMersFlag && useDnaString ){
				processGenomeFromMemory( genomeCount, genomeChromList, countMersFlag );
			}
			else{
				readReference( line, genomeCount, genomeChromList, countMersFlag );
			}
			genomeCount++;
			if( genomeCount > maxNumberRefFiles ){
				break;
			}
		}
	}
}

// Output the number of mers of different sizes up to a maximum of maxHistSize
void outputMerSizes( merLocType* merLocations ){
	//static uint maxHistSize = 3;
	merLocType maxMerSize = 0;
	unsigned long long int count_mers = 0;
	unsigned long long int countMers[maxHistSize+1];
	for( uint i =0; i<=maxHistSize; i++ ){
		countMers[i] = 0;
	}
	for( ullong i =0; i< MER_ARRAY_LENGTH; i++ ){
		//count_mers += merLocations[i].size();
		merLocType size = merLocations[i];
		if( size > maxMerSize ){
			maxMerSize = size;
		}
		count_mers += size;
		if( size >= maxHistSize ){
			countMers[maxHistSize]++;
		}
		else{
			countMers[size]++;
		}
	}
	cout << "total mers: " << dec << count_mers << endl;
	cout << "max Mer Count: "<< maxMerSize << endl;
	cout << "size        count   percent\n";
	for( uint sizeInd = 0; sizeInd <= maxHistSize; sizeInd++ ){
		if( countMers[sizeInd] >= min_display_size ){
			cout.width(2);
			cout << sizeInd;
			cout.width(13);
			cout << countMers[sizeInd]; 
			cout.width(7);
			cout.precision(3);
			cout << "   "<< fixed << (100.0*(float)countMers[sizeInd])/MER_ARRAY_LENGTH << endl;
			//cout << "   "<< fixed << (100.0*(float)countMers[sizeInd])  << endl;
			//cout << endl;
			cout.width(0);
		}
	}
}

/*
void outputFragLocation( uint count, map<uint,fragLocation>::iterator i ){
	cout << count << "\t" << i->second.genome << "\t" << i->second.offset << "\t" << i->second.length << endl;
}
*/

// Converts mer counts in merLocations to a cummulative distribution
ullong convertToCummulative( merLocType* merLocations ){
	merLocType s = merLocations[0];
	merLocType countMers = s;
	merLocations[0] = 0;
	for( ullong i = 1; i < MER_ARRAY_LENGTH; i++ ){
		uint t = s;
		s = merLocations[i];
		countMers += s;
		merLocations[i] = merLocations[i-1]+t;
	}
	cout << "returning from convertToCummulative" << endl;
	cout << "countMers:                        " << countMers << endl;
	cout << "merLocations[MER_ARRAY_LENGTH-1]: " << merLocations[MER_ARRAY_LENGTH-1] << endl;
	return countMers;
}

// shift values in merLocations back to the cummulative distribution
// Not used as of 3/31/14.
void shiftMerLocations( merLocType* merLocations ){
	for( uint i = MER_ARRAY_LENGTH-1; i > 0; i-- ){
		merLocations[i] = merLocations[i-1];
	}
	merLocations[0] = 0;
}

void saveToBinaryFile( string binFileName, ullong countMers, merLocType* merLocations, 
		vector<struct genomeChrom> genomeChromList, short unsigned int* genomeIndices, int* chromLocations ){
	cout << "Saving to binary file: " << binFileName << endl;
	ofstream OUT (binFileName,ios::out | ios::binary);
   if( ! OUT.is_open()) {
      cerr << "Could not open binary output file\n";
      exit(1);
   }
	ullong lengths[3];     
	lengths[0] = merLength;
	lengths[1] = countMers;
	lengths[2] = genomeChromList.size();
	OUT.write(reinterpret_cast<char *>(lengths),3*sizeof(ullong));
	vector<struct genomeChrom>::iterator i;
	for( i = genomeChromList.begin(); i != genomeChromList.end(); i++ ){
		cout << "genome: " << i->genome << "  chromosome: " << i->chromosome << "  dnaLen: " << i->dnaLength<<endl;
		size_t gsz = i->genome.size();
		OUT.write(reinterpret_cast<char *>(&gsz),sizeof(size_t));
		OUT.write(reinterpret_cast<const char *>(i->genome.c_str()),(i->genome.size())*sizeof(char));
		size_t csz = i->chromosome.size();
		OUT.write(reinterpret_cast<char *>(&csz),sizeof(size_t));
		OUT.write(reinterpret_cast<const char *>(i->chromosome.c_str()),(i->chromosome.size())*sizeof(char));
		OUT.write(reinterpret_cast<const char *>(&(i->dnaLength)),sizeof(uint));
	}
	OUT.write(reinterpret_cast<char *>(merLocations),MER_ARRAY_LENGTH*sizeof(merLocType));
	OUT.write(reinterpret_cast<char *>(genomeIndices),countMers*sizeof(short unsigned int));
	OUT.write(reinterpret_cast<char *>(chromLocations),countMers*sizeof(uint));
	OUT.close();
}
