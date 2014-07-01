CC = g++
#CC = icpc
DEBUG = -g
STATIC = -static

cms:  countMers.o   kmerInd.o util.o 
	$(CC) -std=c++0x $(DEBUG) -Wall -o cms countMers.o  kmerInd.o util.o 

readKmerInd:	readKmerIndex.o kmerInd.o util.o countMers.o
	$(CC) -std=c++0x $(DEBUG) -Wall -o readKmerInd readKmerIndex.o kmerInd.o  util.o

kmerInd.o:  kmerInd.cpp countMers.h  util.h
	$(CC) -std=c++0x $(DEBUG) -Wall -O -c kmerInd.cpp

readKmerIndDeb:	readKmerIndDebug.o util.o countMers.o
	$(CC) -std=c++0x $(DEBUG) -Wall -o readKmerIndDeb util.o readKmerIndDebug.o 

countMers.o:  countMers.cpp countMers.h  util.h
	$(CC) -std=c++0x $(DEBUG) -Wall -O -c countMers.cpp

util.o:  util.cpp util.h
	$(CC) -std=c++0x $(DEBUG) -Wall -O -c util.cpp

readKmerIndex.o:	readKmerIndex.cpp countMers.h util.h
	$(CC) -std=c++0x $(DEBUG) -Wall -O -c readKmerIndex.cpp

readKmerIndDebug.o:	readKmerIndDebug.cpp countMers.h util.h
	$(CC) -std=c++0x $(DEBUG) -Wall -O -c readKmerIndDebug.cpp

readQuery:	readQuery.o util.o
	$(CC) -std=c++0x $(DEBUG) -Wall -o readQuery readQuery.o  util.o

readQuery.o:	readQuery.cpp countMers.h util.h
	$(CC) -std=c++0x $(DEBUG) -Wall -O -c readQuery.cpp

subc:	subcontig.o	util.o
	$(CC) -std=c++0x $(DEBUG) -Wall -o subc subcontig.o  util.o 

subcontig.o:	subcontig.cpp countMers.h
	$(CC) -std=c++0x $(DEBUG) -Wall -O -c subcontig.cpp
