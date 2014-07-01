//  Test allocation of a large array following by initializing it to zeros.
#include <iostream>

using namespace std;
typedef unsigned int uint;
typedef unsigned long long ullong;
uint* merLocations;  // This is the in-memory index of k-mers
const int merLength = 15;

unsigned long long intPower( uint base, uint exp ){
  unsigned long long result = 1;
  for( uint i = 0; i < exp; i++ ){
    result *= base;
  }
  return result;
}

int main(){
  ullong MER_ARRAY_LENGTH = intPower(4,merLength);
  cout << dec << "merLength: "<<merLength<<endl;
  try{
    merLocations = new uint[MER_ARRAY_LENGTH];
  }
  catch( bad_alloc ){
    cerr << "could not allocate merLocations array\n";
  }
  cout << "Allocated\n";
  for( ullong i = 0; i < MER_ARRAY_LENGTH; i++ ){
     merLocations[i] = 0;
  }
  cout << "Finished\n";
  return 0;
}
