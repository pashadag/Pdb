#include<iostream>
#include<cstdlib>
#include<fstream>
#include<time.h>
#include<cassert>
using namespace std;

void open_file(ifstream & inFile, string filename) {
	inFile.open(filename.c_str());
	if (!inFile) {
		cerr <<  "Cannot open input file: " <<  filename  << endl;
		assert(0);
		exit(1);
	}
}

string read_genome(string filename) {
	string genome;
	ifstream inFile;
	open_file(inFile, filename);
	char buf[257];

	inFile.get(buf, 256);

	char c;
	while ((c = inFile.rdbuf()->sbumpc()) != EOF) {
		if (c != '\n') {
			c = toupper(c);
			if(c == 'A' || c =='T' || c =='C' || c == 'G') {
				genome.push_back(c);
			}
		}
	}
	inFile.close();
	return genome;
}

int main(int argc, char* argv[])
{
	if (argc != 7) {
		cerr << "Usage: " << argv[0] << " infile l k d delta cov\nNote cov = 0 means spectrum.\n";
		exit(1); } string outFileName = string(argv[1]) + ".readl";
	string outKFileName = string(argv[1]) + ".readk";
	int l = atoi(argv[2]);
	int k = atoi(argv[3]);
	int d = atoi(argv[4]);
	int delta = atoi(argv[5]);
	int cov = atoi(argv[6]);
	string genome = read_genome(argv[1]);
	long numberofReads;
	ofstream outFile ;
	outFile.open(outFileName.c_str());
	ofstream outKFile ; 
	outKFile.open(outKFileName.c_str());
	char* firstPointer = NULL ; 
	char* secondPointer = NULL;
	srand(time(NULL));
	int chromeLength = genome.length();
	if (cov == 0) { //spectrum
		numberofReads = chromeLength;
	} else {
		numberofReads = cov * chromeLength / (2*l);
	}
	cerr<<"Generating " << numberofReads << " reads.\n";
	genome+= genome; //in order to support circularity
	for(long i = 0 ; i < numberofReads  ; i ++)
	{
		long position;
		if (cov == 0) {
			position = i;
		} else {
			position  = rand() % (chromeLength);
		}
		firstPointer = &genome[position]; 
		long distance;
		if (delta == 0) {
			distance = d;
		} else {
			distance = (rand() % (2*delta)) - delta + d;
		}
		secondPointer = & genome[distance + position + l];
		outFile.write(firstPointer, l );
		outFile.put(' ');
		outFile.write(secondPointer,l);
		outFile << ' ' << position << ' ' << distance << endl;
		for(long   j = 0 ; j < l - k +1; j++){
			outKFile.write(& genome[position+j],k);
			outKFile.put(' ');
			//outKFile.write(& chromosome[position + j  + distance + l],k);
			//outKFile.put(' ');
			outKFile << j + 2*l*(i) ;
			outKFile.put('\n');
		}
	}
	outFile.close();
	outKFile.close();
}
