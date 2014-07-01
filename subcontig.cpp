//  program to break up a contig into subcontigs separated by non-ACGT characters
#include "countMers.h"
#include "util.h"
using namespace std;

typedef struct pos_string{
	size_t startPosition;
	string contigStr;
} pos_string;

string contigStr[] = {"ACCTGCAGTTACCGGGNNNTCTGGCATCAGNNTACCTAAAGGCAG", "AAANNAAAANAAN", "NNNGGGGNNGGNN", "NNNANANANA"};
// Used to convert the current mer to an integer via array lookup
char toInt[256];

vector<pos_string> splitContig( string contigStr, vector<pos_string>& subcontigList ){
	uint count = 0;
	size_t nonACGTpos = 0;
	size_t ACGTpos = 0;
	while( (ACGTpos = contigStr.find_first_of(ACGTchars,nonACGTpos)) != string::npos && count < 8 ){
		cout << "ACGTpos: " << ACGTpos << endl;
		nonACGTpos = contigStr.find_first_of(nonACGTchars,ACGTpos+1);
		if( nonACGTpos == string::npos ){
			nonACGTpos = contigStr.length();
		}
		cout << "nonACGTpos: " << nonACGTpos << endl;
		string subContig = contigStr.substr(ACGTpos,nonACGTpos-ACGTpos);
		cout << "subcontig: " << subContig  << endl;
		pos_string ps = {ACGTpos,subContig};
		subcontigList.push_back(ps);
		count++;
	}
	return subcontigList;
}

int main(){
	vector<pos_string> subcontigList;
	for( uint i=0; i<4; i++){
		cout << contigStr[i] << endl;
		subcontigList.clear();
		subcontigList = splitContig( contigStr[i], subcontigList );
		vector<pos_string>::iterator j;
		for( j = subcontigList.begin(); j!= subcontigList.end(); j++ ){
			cout << j->startPosition<< ":"<< j->contigStr << endl;
		}
		cout << "=============="<<endl;
	}
	return 0;
}
