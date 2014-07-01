/* This program creates an in-memory index of all of the k-mers in one collection of DNA
 * sequences.   These are called the reference sequences.  These DNA seqnences are read in 
 * from files listed in the file whose name is stored in the variable referenceFileList.
 * This is set below.
 *
 * There is an associated header file "indexMers.h" and a file of utility functions "util.cpp".
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
 * There are two ways in which the index can be used.
 *
 * 1)  The in-memory index is used to find all exact matches of substrings of length
 * k or longer from a second collection of DNA sequences.  Currently, the first indexed
 * collection of sequences is called the reference and the second collection is called
 * a fragment collection, but it is anticipated that the program will be used in other ways.
 * You can think of the fragment collection as a metagenomic sample.
 * For example, both collections might be collections of fragments which are the ouput from
 * high-throughput sequencing.  Then the objective is to link fragments together to
 * create contigs of genomes from which the fragments are drawn.
 * To run this version, uncomment some fields of "struct db_location" in indexMers.h,
 * and compile with "readFragFile.cpp" as in 
 * "g++ -o inm indexMers.cpp readFragFile.cpp util.cpp".
 * Also, the functions called in the main program will need to be changed.
 *
 * 2)  In this case, the reference sequences are from a collection of known microbial 
 * genomes, and the fragment collection is anticipated to be the reads of a metagenomic 
 * sample.  The objective is to determine the taxonomic unit (such as the phylum) of each
 * fragment. Each reference genome has a unique identifier called the GoldID.  The file
 * "phyloCrossRef.csv" contains the known taxonomic classification for most of the
 * reference genomes.  The in-memory index is used to determine all matches of each fragment 
 * k-mer with a reference sequence.  The taxonomic unit (such as the phylum) of the reference 
 * sequence is determined, and this constitutes a "vote" for the taxonomic unit of the
 * fragment.
 * To run this version, comment out some fields of "struct db_location" in indexMers.h,
 * and compile with "processFrag.cpp" as in 
 * "g++ -o inm indexMers.cpp processFrag.cpp util.cpp".
 *
 *  Author:  Alden Wright, University of Montana CS department,  12/2010 throught 1/2012.
 *
 */

#include "indexMers.h"
using namespace std;

// merLength is now set in indexMers.h
//const int merLength = 14;   // The value of k, i. e., the length of mers used.

// Unfortunately, this program uses way too many global variables.  They are listed here.
unsigned long long MER_ARRAY_LENGTH;  //  MER_ARRAY_LENGTH is set to 4^merLength. 
uint mask; // mask is set to 2^merLength-1.
// currentMer is an integer in the range from 0 to 2^merLength-1 that stores the integer 
//   representation of the current mer in the reference sequence.
uint currentMer = 0;
//vector<db_location>* merLocations;  // This is the in-memory index of k-mers
uint* merLocations;  // This is the in-memory index of k-mers
// Used to convert the current mer to an integer via array lookup
char toInt[256];
char complementArray[256];

// The following are used in the construction of exact matches
map<uint,fragLocation> fragLoc;  // returns the genome and offset of a frag when that is known
pred_succ_struct* predSucc = NULL;

// The following are used for taxonomic classification
vector<string> genomeToGoldId = vector<string>(1000);  // converts genomeIndexes to GoldIds
map<string,phylogeny>* phylogenyMap = NULL;  // maps goldIds to phylogenies
map<string,int> taxLevelCount;  // maps OTU's to the number of k-mers that map to that OTU 
int* k_mer_size=NULL;

// These specify the input files
//string referenceFileList = "shortGenomes.txt";
string referenceFileList = "genomeFiles.txt";
//string referenceFileList = "almostAllGenomes.txt";

// The frag file list file can also be input on the command line
//string fragFileList = "fragFiles.txt";
//string fragFileList = "fastqFiles1_6.txt";
string fragFileList = "allFastq.txt";
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

// updates currentMer to reflect character c.  
//  currentMer is shifted left by 2 bits, the integer representation of c is XOR'd 
//  onto current mer, and then currentMer is reduced back to the appopriate range 
//  by ANDing with the mask.
int nextMerInt( char c ){
  currentMer = ((currentMer<<2)^toInt[(int)c])&mask;
  return currentMer;
}

// Extract GoldId from file name.
// The GoldId is preceded by a period, and then the first character is "G".
string extractGoldId( string& fileName ){
   //cout << "fname: " << fileName << endl;
   size_t pos = fileName.find( ".G", 0 );
   if( pos <= 0 ){
      return "";
   }
   size_t GidEnd = fileName.find( ".", pos+1 );
   if( GidEnd <= 0 ){
      return "";
   }
   string GiD = fileName.substr(pos+1,GidEnd-pos-1);
   return GiD;
}

phylogeny extractPhyloFields( string line, phylogeny phyl ){
  string  temp, strain;
  //vector <string> numbers;
  line = line + ",";

  int count = 0;
  while (line.find(",", 0) != string::npos)
  { //does the string a comma in it?
    size_t  pos = line.find(",", 0); //store the position of the delimiter
    temp = line.substr(0, pos);      //get the token
    line.erase(0, pos + 1);          //erase it from the source 
    //cout << "count: " << count << "   " << temp << endl;
    if( count > 0 ) {
       phyl[count-1] = temp;
       phyl.push_back(temp);                //and put it into the array
    }
    else{  // this string corresponds to the strain
       strain = temp;
    }
    /*
    switch( count ){
      case 0:
        phyl.strain = temp;
        break;
      case 1:
        phyl.goldId = temp;
        break;
      case 2:
        phyl.kingdom = temp;
        break;
      case 3:
        phyl.phylum = temp;
        break;
      case 4:
        phyl.cclass = temp;
        break;
      case 5:
        phyl.order = temp;
        break;
      case 6:
        phyl.family = temp;
        break;
      case 7:
        phyl.genus = temp;
        break;
      case 8:
        phyl.species = temp;
        break;
     }
     phyl.push_back(temp);                //and put it into the array
     */
     count++;
  }
  phyl.push_back(strain);    // add the strain at the end
  /*
  numbers.push_back(line);           //the last token is all alone 
  for( int i = 0; i < numbers.size(); i++ ){
     cout << "Number " << i << " is " << numbers[i] << endl;
  }
  */
  return phyl;
}

void printPhyloFields( phylogeny phyl ){
   /*
   cout << "strain:   " << phyl.strain << endl;
   cout << "goldId:   " << phyl.goldId << endl;
   cout << "kingdom:  " << phyl.kingdom << endl;
   cout << "phylum:   " << phyl.phylum << endl;
   cout << "class:   " << phyl.cclass << endl;
   cout << "order:   " << phyl.order << endl;
   cout << "family:   " << phyl.family << endl;
   cout << "genus:   " << phyl.genus << endl;
   cout << "species:   " << phyl.species << endl;
   */
};

void constructPhylogenyMap(){
  string line;
  try{
    phylogenyMap = new map<string,phylogeny>;
  }
  catch( bad_alloc ){
    cerr << "could not allocate phylogenyMap.\n";
  }
  phylogeny phyl;
  ifstream IN;
  IN.open("phyloCrossRef.csv");
  if(IN.is_open()){
      //~ all is good
  }else{
    cout << "couldn't open file " << "phyloCrossRef.csv" << "\n";
    exit(0);
  }
  getline(IN,line);
  //cout << line << endl;
  while( IN.good() ){
    getline(IN,line);
    //cout << line << endl;
    phyl = vector<string>(9); 
    phyl = extractPhyloFields( line, phyl );
    //taxLevelCount[phyl[taxonomic_level]] = 0;
    //printPhyloFields( phyl );
    phylogenyMap->insert(pair<string,phylogeny>(phyl[0],phyl));
  }
  /*
  map<string,phylogeny>:: iterator it;
  for( it =  phylogenyMap->begin(); it!=phylogenyMap->end(); it++){
     printPhyloFields( it->second );
  }
  */
}

void constructTaxLevelWeights( map<string,int>& taxLevelCount, map<string,double>& taxLevelWeights ){
    map<string,int>::iterator phylCount_iter;
    map<string,double>::iterator phylWeight_iter;
    double weightTotal = 0.0;
    for( phylCount_iter = taxLevelCount.begin(); phylCount_iter != taxLevelCount.end(); phylCount_iter++ ){
       if( useTaxonomicWeights ){
           if( (phylCount_iter->second > 0.005) ){
              if( sqrt_weights ){
                  taxLevelWeights[phylCount_iter->first] = sqrt(1.0/phylCount_iter->second);
              }
              else{  // use inverse weights 
                  taxLevelWeights[phylCount_iter->first] = 1.0/phylCount_iter->second;
              }
              weightTotal += taxLevelWeights[phylCount_iter->first];
           }
           else{
              taxLevelWeights[phylCount_iter->first] = 0.0;
           }
           //cout << phylCount_iter->first << ":: " << taxLevelWeights[phylCount_iter->first] << endl;
       }
       else{   // useTaxonomicWeights is false
          taxLevelWeights[phylCount_iter->first] = 1.0;
       }
    }
    cout << "taxonomic weights" << endl;
    for( phylWeight_iter = taxLevelWeights.begin(); phylWeight_iter != taxLevelWeights.end(); phylWeight_iter++ ){
       if( useTaxonomicWeights ){
          phylWeight_iter->second = phylWeight_iter->second/weightTotal;
          cout << phylWeight_iter->first << ": " << phylWeight_iter->second << endl;
       }
    }
}

/*  Count the k-mers from contigStr and store these in genomeIndex.
 */
void processContig( int contig, string& contigStr, uint genomeCount, string OTU_str ){
  int count=0;  // count of characters processed on a line
  db_location* locptr = NULL;
    //cout << "contig: " << contigStr << endl;
    for( uint i = 0; i < contigStr.length(); i++ ){
      if( contigStr[i] == 'N' ){  
         cerr << "contig contain an N which is not correct!!\n";
         exit(0);
      }
      uint mer = nextMerInt(contigStr[i]);
      if( count >= merLength-1 ){
		   // The next four lines commented out on 11/26/13
         //locptr = new db_location;
         //locptr->genomeIndex =genomeCount;
         //locptr->contigIndex = contig;
         //locptr->contigStart = count-merLength+1;
         //merLocations[mer].push_back(*locptr);  // commented out 11/26/13
         merLocations[mer]++;
         // The following has to do with taxonomic classification only
         taxLevelCount[OTU_str]++;
      }
      count++;  // increment to next character
    }
}

// Reads the reference file and builds the in-memory database of mers in this file.
void readReference(string fileName, uint genomeCount, string GoldId ){
  int contig = 0;  // the index of the contig.  Each genome may have multiple contigs
  string OTU_str;
  /*  Commented out 11/25/13
  map<string,phylogeny>::iterator phylogen_it;
  phylogen_it = phylogenyMap->find(GoldId);
  bool GoldId_found = (phylogen_it!=phylogenyMap->end())?true:false;
  if( ! GoldId_found ){
     cout << "GoldId: " << GoldId << " was not found in the phylogeny map database." << endl;
     return;
  }
  else {
     OTU_str = phylogenyMap->at(GoldId)[taxonomic_level];
     cout << taxLevels[taxonomic_level] << ": " << OTU_str << endl;
  }
  */
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
  cout << "processing " << (fastq_flag?"fragments":"contigs") << " from reference file: " << fileName << endl;
  getline (IN,line);// dummy header line
  string contigStr = "";
  cout << "header line: " << line << "   genomeCount: " << genomeCount << endl;
  bool contigEnd = false;
  size_t firstNpos, lastNpos;  // first and last positions of N character
  while ( IN.good() )  {
    do{  // each line that starts with '>' or contains an N finishes the contig
       getline (IN,line);  // may be a line of genome, or contig label
       //cout << "line: " << line << endl;
       firstNpos = line.find_first_of("N");
       contigEnd = (line.substr(0,1)==">")||( string::npos != firstNpos );
       if( contigEnd && (firstNpos > 0)){
          contigStr += line.substr(0,firstNpos);
       }
       if( !contigEnd ){
          contigStr += line;
       }
    }
    while( IN.good() && (!contigEnd) );
    processContig( contig, contigStr, genomeCount, OTU_str );
    contig++;
    //cout << "contig incremented to " << contig << endl;
    lastNpos = line.find_last_of("N");
    if( lastNpos != string::npos ){
       size_t firstACGT = line.find_first_of("ACGT",firstNpos);
       //cout << "lastNpos: " << lastNpos << "   firstACGT: " << firstACGT << endl;
       while( (firstACGT != string::npos) && (firstACGT < lastNpos) ){
          size_t nextNpos = line.find_first_of("N",firstACGT);
          //cout << "nextNpos: " << nextNpos << "   firstACGT: " << firstACGT << endl;
          contigStr = line.substr(firstACGT,nextNpos-firstACGT);
          //cout<< "cStr: " << contigStr << endl;
          if( nextNpos-firstACGT >= merLength-1 ){
             processContig( contig, contigStr, genomeCount, OTU_str );
          }
          contig++;
          //cout << "contig incremented to " << contig << endl;
          firstNpos = nextNpos;
          firstACGT = line.find_first_of("ACGT",firstNpos);
       }   
       contigStr = line.substr(lastNpos+1);
    }
    else{ 
       contigStr = "";
    }
  }
  cout << "processed " << contig << (fastq_flag?" fragments":" contigs") << " from reference file: " << fileName << endl;
}

void readListOfReferenceFiles( string fileName ){
  uint genomeCount = 0;
  ifstream IN;
  IN.open(fileName.c_str());
  if(IN.is_open()){
    //~ all is good
  }else{
    cout << "couldn't open file " << fileName << " with list of reference file names\n";
    exit(0);
  }
  string line;
  int dbseq = 0;
  while ( IN.good() )  {
     getline (IN,line);
     if( (line.length() > 0) && (line.substr(0,2) !=  "//") ){
       string GoldId = extractGoldId( line );
       //cout << "GoldId: " << GoldId << endl;
       genomeToGoldId[genomeCount] = GoldId;
       readReference( line, genomeCount, GoldId );
       genomeCount++;
     }
  }
  // allocate an array with one element per frag
  //predSucc = new pred_succ_struct[dbseq];
  cout << "exiting from readListOfReferenceFiles\n";
}
/* Commented out on 11/26/13
void readListOfFragFiles( string fileName, map<string,double>taxLevelWeights ){
  ifstream IN;
  IN.open(fileName.c_str());
  if(IN.is_open()){
    //~ all is good
  }else{
    cout << "couldn't open file " << fileName << " with list of frag files\n";
    exit(0);
  }
  string fragFileName;
  uint fragFileCount = 0;   // assume that each frag file corresponds to a genome
  uint fragCount = 0;
  while ( IN.good() )  {
     getline (IN,fragFileName);
     if( (fragFileName.length() > 0) && (fragFileName.substr(0,2) !=  "//") ){
       readFragFile( fragFileName, merLength, fragCount, fragFileCount, taxLevelWeights );
       fragFileCount++;
       //cout << "fragFileCount incremented to " << fragFileCount << endl;
     }
     fragCount=0;
  }
}
*/

// Prunes merLocations array to those elements within a certain size range
void pruneMerArray( vector<db_location>* merLocations ){
  cout << "pruning merLocations array of elemements whose size is greater than " << maxMerArraySize << endl;
  if( minMerArraySize > 0 ){
     cout << "pruning merLocations array of elemements whose size is less than " << minMerArraySize << endl;
  }
  for( uint i =0; i< MER_ARRAY_LENGTH; i++ ){
    if( (merLocations[i].size() > maxMerArraySize) || (merLocations[i].size() < minMerArraySize) ){
       merLocations[i].clear();
    }
  }
}

// Output the number of mers of different sizes up to a maximum of maxMerSize
void outputMerSizes( uint* merLocations ){
  static int maxMerSize = 3;
  int count_mers = 0;
  int countMers[maxMerSize+1];
  for( uint i =0; i<=maxMerSize; i++ ){
     countMers[i] = 0;
  }
  vector<db_location>::iterator it;
  for( uint i =0; i< MER_ARRAY_LENGTH; i++ ){
    //count_mers += merLocations[i].size();
    int size = merLocations[i];
    count_mers += size;
    if( size >= maxMerSize ){
	 countMers[maxMerSize]++;
    }
    else{
       countMers[size]++;
    }
  }
  cout << "total mers: " << count_mers << endl;
  cout << "size   count\n";
  for( uint sizeInd = 0; sizeInd <= maxMerSize; sizeInd++ ){
    cout << sizeInd<< "    " << countMers[sizeInd] << endl;
  }
}


// Displays all non-empty locations in the mer array
void outputMerArray( uint* merLocations ){
  const int min_display_size = 20;
  vector<db_location>::iterator it;
  int count_mers = 0;
  int count_large_mers = 0;
  for( uint i =0; i< MER_ARRAY_LENGTH; i++ ){
    int size = merLocations[i];
    count_mers += size;
    if( size >= min_display_size ){
 
       cout << intToMer(i,merLength) << " : " << merLocations[i]; 
        /*
       for( it = merLocations[i].begin(); it!=merLocations[i].end(); it++){
         cout << " : "<< it->genomeIndex << "." << it->contigIndex << "." << it->contigStart;
             //"." << it->contigIndex << "." << it->contigOffset << 
       }
       */
       cout << endl;
       count_large_mers += merLocations[i];
    }
  }
  cout << "Total mers: " << count_mers << endl;
  cout << "Large mers: " << count_large_mers << endl;
}

void outputFragLocation( uint count, map<uint,fragLocation>::iterator i ){
   cout << count << "\t" << i->second.genome << "\t" << i->second.offset << "\t" << i->second.length << endl;
}

void outputAllFragLocations( ){
   cout << "genome\toffset\tlength\n";
   map<uint,fragLocation>::iterator i;
   uint count = 0;
   for( i = fragLoc.begin(); i != fragLoc.end(); i++ ){
      outputFragLocation( count, i );
      count++;
   }
}

// Displays a histogram of mer frequencies
// Outdated--- not used
void merHistogram( vector<db_location>* merLocations ){
  cout << "MER_ARRAY_LENGTH: " <<  MER_ARRAY_LENGTH << endl;
  uint num = 1000; 
  cout << "num: " << num << endl;
  k_mer_size = (int *) realloc( k_mer_size, num*sizeof(int));
  for( uint j = 0; j < num; j++ ){
     k_mer_size[j]=0;
  }
  uint maxindex = 0;
  for( uint i =0; i< MER_ARRAY_LENGTH; i++ ){
     /*
     if( merLocations[i].size() > 0 ){
       cout << intToMer(i,merLength) << " : " << merLocations[i].size() << endl; 
     }
     */
     uint size = merLocations[i].size();
     if( size > maxindex ){
        maxindex = size;
     }
     if( size > num){
        uint old_num = num;
        num = (size >= 2*num)? size+1:2*num;
        cout << "num: " << num << endl;
        k_mer_size = (int *) realloc( k_mer_size, num*sizeof(int));
        for( uint j = old_num; j < num; j++ ){
           k_mer_size[j]=0;
        }
     }
     k_mer_size[size]++;
  }
  cout << "Histogram of mer count frequencies" << endl;
  cout << "freq\tcount\n";
  for( uint j = 0; j <= maxindex; j++ ){
     if( k_mer_size[j] > 0 ){
        cout << j << "\t" << k_mer_size[j] << endl;
     }
  } 
  cout << endl;
}

int main(int argc, char* argv[]){
  cout << "Taxonomic level: " << taxLevels[taxonomic_level] << endl;
  cout << "Not counting votes for k-mer matches between a frag and its generating genome (in processFrag.cpp)." << endl;
  if( useTaxonomicWeights ){
     if( sqrt_weights ){
        cout << "Taxonomic weights used with sqrt of inverse of frequency." << endl;
     }
     else{
        cout << "Taxonomic weights used with inverse of frequency." << endl;
     }
  }
  else{
     cout << "Taxonomic weights not used." << endl;
  }
  map<string,double> taxLevelWeights;
  if( argc > 1 ){
     fragFileList = argv[1];
  }
  //constructPhylogenyMap();
  MER_ARRAY_LENGTH = intPower(4,merLength);
  if( merLength == 16 ){  // 16 is the maximum for merLength
    mask = ~0;
  }
  else{
    mask = MER_ARRAY_LENGTH-1;
  }
  cout << dec << "merLength: "<<merLength<<endl;
  cout << "array length: " << MER_ARRAY_LENGTH << endl;
  cout << "mask: " << hex << mask << endl;
  //cout << dec;
  try{
    //merLocations = new vector<db_location>[MER_ARRAY_LENGTH];
    merLocations = new uint[MER_ARRAY_LENGTH];
  }
  catch( bad_alloc ){
    cerr << "could not allocate merLocations array\n";
  }
  for( uint i = 0; i < MER_ARRAY_LENGTH; i++ ){
     merLocations[i] = 0;
  }
  cout << dec << "array length: " << MER_ARRAY_LENGTH << endl;
  toIntSetup();
  /*
  cout << "toIntSetup" << endl;
  for( int i = 0; i < 256; i++ ){
     cout << "i: " << i << "  toInt: " << toInt[i] << endl;
  }
  */
  readListOfReferenceFiles(referenceFileList);
  cout << "After readListOfReferenceFiles\n";
  
  // Added 11/26/13 to only print histogram
  //outputMerArray( merLocations );
  outputMerSizes( merLocations );
  exit(1);

  /*  Commented out on 11/26/13

  map<string,int>:: iterator it;
  for( it = taxLevelCount.begin(); it != taxLevelCount.end(); it++ ){
     cout << it->first << ": " << it->second << endl;
  }
  
  constructTaxLevelWeights(taxLevelCount, taxLevelWeights );
  //merHistogram( merLocations );
  //outputMerArray( merLocations);  // Displays the nonempty elements of the reference array
  pruneMerArray( merLocations);  // clears the large and small elements of the reference array
  //outputAllFragLocations( );
  //merHistogram( merLocations );
  
  //for( int i =0; i< MER_ARRAY_LENGTH; i++ ){
  //  cout << merLocations[i].size() << ": "<< merLocations[i].capacity() << endl;
  //}
  
  //outputPredSuccHeader();
  readListOfFragFiles(fragFileList, taxLevelWeights );
  cout << "run completed successfully\n";
*/
  return 0;
}
