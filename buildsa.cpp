#include<cassert>
#include<cmath>
#include<string>
#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<map>
#include<list>
#include<queue>
#include<cstdarg>
#include<algorithm>
//#include<mysql++.h>
//#include <boost/foreach.hpp>

using namespace std;
ostringstream oMsg;
string sbuf;
string filename, baseFilename;


string genome;
int genSize;

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

	//read the file into genome
	inFile.get(buf, 256);

	char c;
	while ((c = inFile.rdbuf()->sbumpc()) != EOF) {
		if (c != '\n') {
			genome.push_back(c);
		}
	}
	inFile.close();
	return genome;
}

bool complt(int x, int y) {
	if (x==y) return false;
	int i =0;
	while ((i < genSize) && (genome[x+i] == genome[y+i])) i++;
	return genome[x+i] < genome[y+i];
}


void usage(int argc, char * argv[]) {
	cerr << "Usage: " << argv[0] << " <genome> " << endl;
	exit(1);
}


int main(int argc, char * argv[]) {

	if (argc != 2)  usage(argc, argv);

	vector<int> sa;
	genome = read_genome(argv[1]);
	genSize= genome.size();

	genome += genome; //to support circularity

	for (int i = 0; i < genSize; i++) {
		sa.push_back(i);
	}
	sort(sa.begin(), sa.end(), complt);
	for (int i = 0; i < sa.size() ; i++) cout << sa[i] << endl;
	return 0;


}

