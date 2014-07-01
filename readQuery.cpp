#include "countMers.h"
#include "util.h"
using namespace std;

//uint* merLocations;
//usint* genomeIndices;
//uint* chromLocations;
vector<struct genomeChrom> genomeChromList;
ullong MER_ARRAY_LENGTH;  // set to 4^merLength
ullong mask; // mask is set to 2^merLength-1.
// currentMer is an integer in the range from 0 to 2^merLength-1 that stores the integer 
//    representation of the current mer in the sequence read from the query file.
ullong currentMer = 0;

// Characters that are neither ACTG nor nonACGTchars as defined in countMers.h
const string badChars = "0123456789!@#$%^&*()_+-={}|[]:\";'<>?,./abcdefghijklmnopqrstuvxyzEFIJLOPQUZ";

char toInt[256];
struct genChromLoc{
	usint genomeIndex;
	uint  chromLoc;
};
void processContig( string& contigStr, uint genomeCount, uint contigStartPos, bool countMersFlag ){
	//cout << "contig: " << contigStr << endl;
	if( contigStr.find_first_of(badChars) != string::npos ){
		uint n = contigStr.find_first_of(badChars);
		cerr << "bad character " << contigStr[n] << " in contig at position " << n << endl;
		exit(1);
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
			currentMer = ((currentMer<<2)^toInt[(int)contigStr[i]])&mask;
		 	ullong mer = currentMer;
			if( count >= merLength-1 ){
				//cout << "mer: " << intToMer(mer,merLength) << endl;
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

void splitContig( string contigStr, uint genomeCount, bool countMersFlag ){
	//cout <<"splitContig: " <<contigStr<<endl;
	//cout <<genomeCount<<endl;
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
		processContig( subContig, genomeCount, ACGTpos, countMersFlag );
	}
}

// Equivalent to functionReadReference below except that it uses the dnaString stored in
//    genomeChromList rather than rereading it from the file
void processGenomeFromMemory( uint genomeCount, vector<genomeChrom>& genomeChromList, bool countMersFlag ){
	splitContig( genomeChromList[genomeCount].dnaString, genomeCount, countMersFlag );
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
	string descStr = line;
	//cout << "header line: " << line << "   genomeCount: " << genomeCount << endl;
	bool contigEnd = false;
	//size_t firstNpos, lastNpos;  // first and last positions of nonACGT character
	while ( IN.good() )  {
		if( descStr.length() > 0 ){
			//cout << "desc str: " << descStr<< endl;
		}
		contigStr = "";
		descStr = "";
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
				descStr += line;
			}
			contigCount++;
		}
		while( IN.good() && (!contigEnd) );
		//DNA_string += contigStr;  //  The whole DNA string for this file
		int genomeCount = 0;
		bool countMersFlag = true;
		//cout << "desc str: " << descStr<< endl;
		splitContig( contigStr, genomeCount, countMersFlag );
	}
	cout << "processed " << contigCount << (fastq_flag?" fragments":" contigs") << " from query file: " << fileName << endl;
	//cout << "DNAstr: " << DNA_string << endl;
}

int main( int argc, char* argv[] ){
	//uint* merLocations;
	//usint* genomeIndices;
	//uint* chromLocations;
	//vector<struct genomeChrom> genomeChromList;
	//unsigned long long MER_ARRAY_LENGTH;
	toIntSetup();
	//cout << "toInt: " << (int)toInt[(int)'T'] << endl;
	//string binaryFileName = "countMers.bin";
	string queryFileName = "query.fasta";
	MER_ARRAY_LENGTH = intPower(4,merLength);
	mask = MER_ARRAY_LENGTH-1;
   cout << dec << "merLength: "<<merLength<<endl;
   cout << "mask: " << hex << mask << endl;
	cout << dec;
	/*
	if( argc > 2 ){
		fastaFileName = argv[2];
   }
	*/
	if( argc > 1 ){
		queryFileName = argv[1];
   }
   else{
		//cerr << "USAGE:  " << argv[0] << " <binary file name> <fasta file name>" << endl; 
		cerr << "USAGE:  " << argv[0] << "  <query file name>" << endl; 
	}
	/*
	ifstream IN (queryFileName, ios::in);
	if (IN.is_open()) {
		// all is good
	}
	else{
		cerr << "unable to open binary input file " << binaryFileName << endl;
		exit(1);
	}
	*/
	cout << "merLength: " << merLength << endl;
	cout << "query file name: " << queryFileName << endl;
	readQueryFile(queryFileName );
	return 0;
}
