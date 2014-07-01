#include "countMers.h"
#include "util.h"
using namespace std;

uint* merLocations;
usint* genomeIndices;
uint* chromLocations;
vector<struct genomeChrom> genomeChromList;
ullong MER_ARRAY_LENGTH;  // set to 4^merLength
ullong mask; // mask is set to 2^merLength-1.
ullong queryMask; // queryMask is set to 2^queryLength-1.
// currentMer is an integer in the range from 0 to 2^merLength-1 that stores the integer 
//    representation of the current mer in the sequence read from the query file.
ullong currentMer = 0;
//string specialMer = "ACTTTAGATAATCCAA";
//ullong specialMerInt = 0x1fc8c350;   // look at  bioinf/cpp/stringToMer.cpp
string specialMer = "ACAATGTAGACCATAT";
ullong specialMerInt = 0x10ec8533;   // look at  bioinf/cpp/stringToMer.cpp

char toInt[256];
struct genChromLoc{
	usint genomeIndex;
	uint  chromLoc;
};

// Given a K-mer (where K=merLength as set in countMers.h), returns a vector of locations of that K-mer.
// A location is a pair where the first element of the pair is the index of a genome, and the second
//     element of the pair is the location in the dna string (i. e., chromosome) of that genome.
vector <struct genChromLoc> lookup_K_mer( uint mer, vector<struct genomeChrom> genomeChromList, uint* merLocations, 
		usint* genomeIndices, uint* chromLocations, vector<struct genChromLoc>& gchromlocList ){
	uint startIndex;
	gchromlocList.clear();
	if( mer == 0 ){
		startIndex = 0;
	}
	else{
		startIndex = merLocations[mer-1];
	}
	for( uint j = startIndex; j < merLocations[mer]; j++ ){
		usint genomeIndex = genomeIndices[j];
		uint chromLoc = chromLocations[j];
		struct genChromLoc gc = {genomeIndex,chromLoc};
		gchromlocList.push_back(gc);
	}
	return gchromlocList;
}

// Given a L-mer (where K<=L<=2*k, and K=merLength as set in countMers.h), returns a vector of locations of that L-mer.
// A location is a pair where the first element of the pair is the index of a genome, and the second
//     element of the pair is the location in the dna string (i. e., chromosome) of that genome.
vector <struct genChromLoc> lookup_L_mer( int L, uint mer, vector<struct genomeChrom> genomeChromList, uint* merLocations, 
		usint* genomeIndices, uint* chromLocations, vector<struct genChromLoc>& gchromlocList ){
	if( L > 2*merLength ){
		cerr << "L = "<<L<<" cannot be more than double merLength=K= "<<merLength<<endl;
		cerr << "exiting"<<endl;
		exit(1);
	}
	if( L == merLength ){
		return lookup_K_mer( mer, genomeChromList, merLocations, genomeIndices,chromLocations,gchromlocList);
	}
	gchromlocList.clear();
	uint mer1 = mer >>2*(L-merLength);
	uint mer2 = mer % intPower(4,merLength);  // Maybe do more efficiently by shifting or masking
	//cout << "mer1: " << mer1 << "  mer2: " << mer2 << endl;
	//cout << "mer1: " << intToMer(mer1,merLength) << "  mer2: " << intToMer(mer2,merLength) << endl;
	uint startIndex1, startIndex2;
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
	uint j1 = startIndex1;
	uint j2 = startIndex2;
	while( j1 < merLocations[mer1] && j2 < merLocations[mer2] ){
		usint genomeIndex1 = genomeIndices[j1];
		uint chromLoc1 = chromLocations[j1];
		//cout << "gI1: " << genomeIndex1 << "  cL1: " << chromLoc1 << endl;
		usint genomeIndex2 = genomeIndices[j2];
		uint chromLoc2 = chromLocations[j2];
		//cout << "gI2: " << genomeIndex2 << "  cL2: " << chromLoc2 << endl;
		if( genomeIndex1 == genomeIndex2 && chromLoc1+L-merLength == chromLoc2 ){
			struct genChromLoc gc = {genomeIndex1, chromLoc1};
			gchromlocList.push_back(gc);
			j1++;  j2++;
		}
		else if( genomeIndex1 < genomeIndex2 ){
			j1++;
		}
		else if( genomeIndex1 > genomeIndex2 ){
			j2++;
		}
		else if( genomeIndex1 == genomeIndex2 && chromLoc1+L-merLength < chromLoc2 ){
			j1++;
		}
		else if( genomeIndex1 == genomeIndex2 && chromLoc1+L-merLength > chromLoc2 ){
			j2++;
		}
	}
	return gchromlocList;
}


void queryMer( ullong queryMer, string queryDescriptor ){
	if( queryMer == specialMerInt ){
		cout << "special mer in query mer\n";
	}
	vector <struct genChromLoc> queryLocList;
	queryLocList = lookup_L_mer( queryLength, queryMer, genomeChromList,  merLocations, 
		genomeIndices, chromLocations, queryLocList );
	if( !queryLocList.empty() ){
		vector <struct genChromLoc>::iterator i;
		cout << "query results for mer: " << intToMer( queryMer, queryLength ) << "  len: " << queryLocList.size() <<endl;
		cout << "from query descriptor: " << queryDescriptor << endl;
		for( i = queryLocList.begin(); i != queryLocList.end(); i++ ){
			cout << i->genomeIndex << " : " << i->chromLoc << endl;
			cout <<genomeChromList[i->genomeIndex].genome << " : " << genomeChromList[i->genomeIndex].chromosome << endl;
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

void processContig( string& contigStr, string queryDescriptor ){
	cout << "contig: " << contigStr << endl;
	currentMer = 0;
	if( contigStr.find(specialMer) != string::npos ){
		cout << "special mer " << specialMer << " found in contig " << contigStr << endl;
	}
	int count=0;  // count of characters processed on a line
		 if( contigStr.find_first_of(nonACGTchars) != string::npos ){  // Comment out for speed
			 cerr << "contig contains a non-ACGT character which is not correct!!\n";
			 cerr << contigStr << endl;
			 exit(0);
		 }
		for( uint i = 0; i < contigStr.length(); i++ ){
			//uint mer = nextMerInt(contigStr[i]);
			//cout << (int)toInt[(int)contigStr[i]] << endl;
			currentMer = ((currentMer<<2)^toInt[(int)contigStr[i]])&queryMask;
		 	ullong mer = currentMer;
			if( count >= queryLength-1 ){
				//cout << "mer: " << intToMer(mer,queryLength) << endl;
				queryMer( mer, queryDescriptor );
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

void splitContig( string contigStr, string queryDescriptor ){
	//cout <<"splitContig: "<<genomeCount<<endl;
	//cout <<contigStr<<endl;
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
		processContig( subContig, queryDescriptor );
	}
}

// Equivalent to functionReadReference below except that it uses the dnaString stored in
//    genomeChromList rather than rereading it from the file
void processGenomeFromMemory( uint genomeCount, string queryDescriptor ){
	splitContig( genomeChromList[genomeCount].dnaString, queryDescriptor );
}

// read a fasta file as a source of queries
void readQueryFile(string fileName ){
	// A "contig" here is part or all of the DNA string for this file which is delimited by non-ACGT characters.
	uint contigCount = 0;  // the index of the contig.  Each genome may have multiple contigs
	ifstream IN;
	IN.open(fileName.c_str());
	if(IN.is_open()){
		//~ all is good
	}else{
		cout << "couldn't open query file " << fileName << endl;
		exit(0);
	}
	// fastq_flag is true if file is a fragment file.
	bool fastq_flag = (fileName.substr(fileName.length()-6,6)==".fastq") ? true : false ;
 
	//~ #~ read in dummy line. For a fasta file, all following lines are part of genome.
	string line;
	cout << "processing " << (fastq_flag?"fragments":"contigs") << " from reference file: " << fileName << endl;
	getline (IN,line);// descriptor line
	if( line[0] != '>' ){
		cout << "Error: the following line is supposed to be a descriptor line which starts with a > character." << endl;
		exit(1);
	}
	string DNA_string = "";
	string contigStr = "";
	string queryDescriptor = line;
	//cout << "header line: " << line << "   genomeCount: " << genomeCount << endl;
	bool contigEnd = false;
	//size_t firstNpos, lastNpos;  // first and last positions of nonACGT character
	while ( IN.good() )  {
		if( queryDescriptor.length() > 0 ){
			//cout << "desc str: " << descStr<< endl;
		}
		contigStr = "";
		queryDescriptor = "";
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
				queryDescriptor += line;
			}
			contigCount++;
		}
		while( IN.good() && (!contigEnd) );
		//DNA_string += contigStr;  //  The whole DNA string for this file
		//cout << "desc str: " << descStr<< endl;
		splitContig( contigStr, queryDescriptor );
	}
	cout << "processed " << contigCount << (fastq_flag?" fragments":" contigs") << " from query file: " << fileName << endl;
	//cout << "DNAstr: " << DNA_string << endl;
}

// For each k-mer in referenced by genomeIndices and chromLocations, check that the k-mer actually
//    occurs in the genome at the specified location.
// If the printAllMers flag is set in countMers.h, then print each such k-mer.
void checkKmerLocations( 
		uint* merLocations, vector<struct genomeChrom> genomeChromList, short unsigned int* genomeIndices, uint* chromLocations ){
	uint startIndex;
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

int main( int argc, char* argv[] ){
	//uint* merLocations;
	//usint* genomeIndices;
	//uint* chromLocations;
	//vector<struct genomeChrom> genomeChromList;
	//unsigned long long MER_ARRAY_LENGTH;
	cout << "special mer: " << specialMer << endl;
	toIntSetup();
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
	ifstream IN (binaryFileName, ios::in|ios::binary|ios::ate);
	if (IN.is_open()) {
		// all is good
	}
	else{
		cerr << "unable to open binary input file " << binaryFileName << endl;
		exit(1);
	}
	uint* lengths;
	lengths = new uint[3];
	IN.seekg (0, ios::beg);
	IN.read(reinterpret_cast<char *>(lengths),3*sizeof(uint));
	int myMerLength = lengths[0];
	uint countMers = lengths[1];
	uint genomeChromListLength = lengths[2];
	if( myMerLength != merLength ){
		cerr << "merLength from countMers.h: " << merLength << endl;
		cerr << "merLength from " << binaryFileName << ": " << myMerLength << endl;
		cerr << "These do not agree, so exiting!!" << endl;
		exit(1);
	}
	cout << "merLength: " << myMerLength << endl;
	cout << "queryLength: " << queryLength << endl;
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
		struct genomeChrom gtemp = {ii,string(genomeStr), string(chromStr),string("")};
		genomeChromList.push_back(gtemp);
		cout << genomeChromList[ii].genome<<"  "<<genomeChromList[ii].chromosome<<endl;
	}
   try{
		merLocations = new uint[MER_ARRAY_LENGTH];
   }
   catch( bad_alloc ){
		cerr << "could not allocate merLocations array\n";
   }
	IN.read(reinterpret_cast<char *>(merLocations),MER_ARRAY_LENGTH*sizeof(uint));
	/*
	for( uint i = 0; i < MER_ARRAY_LENGTH; i++ ){
		cout << i << ": "<<merLocations[i] << endl;
	}
	*/
	try{
		genomeIndices = new usint[countMers];
		chromLocations = new uint[countMers];
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

	readQueryFile(fastaFileName);
	if( checkKmerLocationsFlag && saveDnaString ){
		checkKmerLocations( merLocations, genomeChromList, genomeIndices, chromLocations );
	}
	return 0;
}
