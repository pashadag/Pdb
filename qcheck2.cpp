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

bool complt(int y, string s) {
	int i = 0;
	while ((i < s.length()) && (s[i] == genome[y+i])) i++;
	if (i == s.length()) return false;
	return  genome[y+i] < s[i];
}

bool compeq(int y, string s ) {
	int i = 0;
	while ((i < s.length()) && (s[i] == genome[y+i])) i++;
	if (i == s.length()) return true;
	return false ;
}

void usage(int argc, char * argv[]) {
	cout << "Usage: " << argv[0] << " <genome fasta file> <genome suffix array> " << endl;
	cout << "Standard in should be a fasta file of conigs.\n";
	exit(1);
}


int main(int argc, char * argv[]) {

	if (argc != 3)  usage(argc, argv);

	genome = read_genome(argv[1]);
	vector<int> sa;
	ifstream safile;
	
	open_file(safile, argv[2]);

	genome += genome; //to support circularity

	int ibuf;
	while (safile >> ibuf) sa.push_back(ibuf);
		
	while (getline(cin, sbuf)) {
		string contig(sbuf);
		getline(cin, sbuf);
		vector<int>::iterator res = lower_bound(sa.begin(), sa.end(), sbuf, complt);
		if ((res == sa.end()) || (!compeq(*res, sbuf))) {
			cout << "NO\t" << contig << "\t" << sbuf.length() << endl;
			//cout << "closest = " << *res << "\t" << genome.substr(*res, 100) << endl;
		} else {
			cout << "YES\t" << contig << "\t" << sbuf.length()<< "\t" << *res <<  endl;
		}
	}
	return 0;


}

