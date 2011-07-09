#include<iostream>
#include<cstdlib>
#include<fstream>
#include<time.h>
using namespace std;
//this code generate single spectrum 


int main(int argc, char* argv[])
{


	if (argc != 5) {
		cerr << "Usage: " << argv[0] << " infile l k cov\nNote cov = 0 means spectrum.\n";
		exit(1);
	}
	cerr << "WARNING: This program has a bug that has not been corrected.  When the infile contains a fasta header line, the [acgt]s in that line are read into the genome.\n";
	char * inputFilename = argv[1];
	string outFileName(inputFilename);
	string outKFileName(inputFilename);
	outFileName += ".readl";
	outKFileName += ".readk";
	int l = atoi(argv[2]);
	int k = atoi(argv[3]);
	int cov = atoi(argv[4]);
	long numberofReads;

	//read and get one chromosome to an array
	ifstream::pos_type size;
	char* memblock;
	long counter = 0;
	ifstream file(inputFilename, ios::in|ios::binary|ios::ate);
	if(file.is_open())
	{
		size = file.tellg();
		memblock = new char[size];
		file.seekg(0,ios::beg);
		file.read(memblock,size);
		file.close();
		cout<<"The complete file is in memory"<<size<<endl;
	}
	long chromeLength = 0;
	for(long i = 0 ; i < size ; i++)
	{
		char c = toupper(memblock[i]);
		if(c == 'A' || c =='T' || c =='C' || c == 'G')
			chromeLength++;
	}
	char* chromosome = new char[chromeLength];
	counter =0;
	for(long i = 0 ; i < size ; i++)
	{
		char c = toupper(memblock[i]);
		if(c == 'A' || c =='T' || c =='C' || c == 'G') {
			chromosome[counter] = c;
			counter++;
		}
	}
	delete[] memblock;
	ofstream outFile ;
	outFile.open(outFileName.c_str());
	ofstream outKFile ; 
	outKFile.open(outKFileName.c_str());

	if (cov == 0) { //spectrum
		numberofReads = chromeLength - l + 1;
	} else {
		numberofReads = cov * chromeLength / l;
	}


	char* firstPointer = NULL ; 
	srand(time(NULL));
	for(long i = 0 ; i < numberofReads; i ++)
	{
		long position;
		if (cov == 0) {
			position = i; 
		} else {
			position  = rand() % (chromeLength - l +1 );
		}
		firstPointer = &chromosome[position]; 
		outFile.write(firstPointer, k+1 );
		outFile.put(' ');
		outFile << position;
		outFile.put('\n');
        for(long   j = 0 ; j < l - k +1; j++){
			outKFile.write(& chromosome[position+j],k);
			outKFile.put(' ');
			outKFile << j + 2*(l)*(i) ;
			outKFile.put('\n');
		}

	}
	outFile.close();
	outKFile.close();

}
