#include<iostream>
#include<string.h>
#include<vector>
#include<set>
#include<deque>
#include<set>
#include<fstream>
#include<cstdlib>
#include"union.h"
#include"pthread.h"
#define _DEBUG_ 
using namespace std ;


char* allReadConcatenated ; 
int k;
int L ;
int delta;


void open_file(ifstream & inFile, string filename) {
	inFile.open(filename.c_str());
	if (!inFile) {
		cerr <<  "Cannot open input file: " <<  filename  << endl;
		assert(0);
		exit(1);
	}
}
void open_file(ofstream & file, string filename) {
	file.open(filename.c_str());
	if (!file) {
		cerr <<  "Cannot open input file: " <<  filename  << endl;
		assert(0);
		exit(1);
	}
}



long GetNumberOfLine(ifstream& file)
{
	string line ;
	long counter = 0;
	while(getline(file,line))
	{
		counter++; 
	}

	file.seekg(0, ios::beg);
	return counter ;

}


void ConcatenateAllReads(char* array, ifstream& file)
{
	int bufsize = 200000000;
	//char buff[bufsize];
	char* buff = new char[bufsize];
	long counter = 0;
	cerr<<"begin read"<<endl;
	while(!file.eof())
	{
		file.read(buff,bufsize);
		long numberofCharRead = file.gcount();

		for(long i = 0 ; i < numberofCharRead ; i++)
		{
			if(buff[i] == 'A' || buff[i] == 'T' || buff[i] == 'G' || buff[i] == 'C' )
			{
				array[counter] = buff[i];
				counter++;
			}
		} 
	}
	delete buff;
	cerr<<"end filter "<<endl;
	return ;
}


int main(int argc, char* argv[])
{

	if(argc != 4)
	{
		cerr<<"generateKs base L k "<<endl;
		exit(1);
	}

	string base = argv[1];
	string inputReadsFileName(base);
	inputReadsFileName += ".readl";
	L = atoi(argv[2]) ;
	k = atoi(argv[3]);

	ifstream inputReadFile ;
	inputReadFile.open(inputReadsFileName.c_str());
	if(!inputReadFile.is_open())
	{
		cerr<<"Can not open file "<< inputReadsFileName << endl;
		exit(1);
	}
	long numberOfLines = GetNumberOfLine(inputReadFile);
	allReadConcatenated = new char[numberOfLines*2*L];
	inputReadFile.close();
	ifstream newInputReadFile ;
	open_file( newInputReadFile, inputReadsFileName);  
	ConcatenateAllReads(allReadConcatenated, newInputReadFile);
	newInputReadFile.close();


	ofstream outFile;
	outFile.open((base + ".readk").c_str());
	
	int line = 0;
	while (line < numberOfLines) {
		for (int offset = 0; offset <= L - k; offset++) {
			int pos = 2 * L * line + offset;
			outFile.write(&allReadConcatenated[pos], k);
			outFile << " " << pos << endl;
		}
		line++;
	}

	outFile.close();



	delete allReadConcatenated;
	return 0;

}

