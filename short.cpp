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
long genSize;
long L;

void open_file(ifstream & inFile, string filename) {
	inFile.open(filename.c_str());
	if (!inFile) {
		cerr <<  "Cannot open input file: " <<  filename  << endl;
		assert(0);
		exit(1);
	}
}
void open_file(ofstream & inFile, string filename) {
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

bool complt(long x, long y) {
	if (x==y) return false;
	long i =0;
	while ((i < genSize) && (genome[x+i] == genome[y+i])) i++;
	return genome[x+i] < genome[y+i];
}


void usage(int argc, char * argv[]) {
	cerr << "Usage: " << argv[0] << " <genome> <l>" << endl;
	cerr << "This program is like singleGenerateReads, but it assumes k = l -1, and cov=0.\n";
	cerr << "However, it generates the readks file directly, using the fact the genome as a cheat in order to do a much faster sort than unix would.\n";

	exit(1);
}


int main(int argc, char * argv[]) {

/*	int x = 30000000;
	L = 1000;
	cout << fixed <<  2*L*(long(x)-1) + 1  << endl;
	return 0;
	*/
	if (argc != 3)  usage(argc, argv);

	vector<long> sa;
	string base = argv[1];
	genome = read_genome(base);
	L = atoi(argv[2]);

	genSize= genome.size();

	genome += genome; //to support circularity

	cerr << "Building suffix array...\n";
	for (long i = 0; i < genSize; i++) {
		sa.push_back(i);
	}
	sort(sa.begin(), sa.end(), complt);
	ofstream out;

	/*cerr << "Building readl...\n";
	open_file(out, base + ".readl");
	for (int i = 0; i < genSize; i++) out << genome.substr(i,L) << endl;
	out.close();
	*/

	cerr << "Building readks...\n";
	open_file(out, base + ".readks");
	for (long i = 0; i < sa.size() ; i++) {
		long x = sa[i];
		//long x = sa[i];
		if (x == 0) {
			x = genSize * 2 * long(L) + 1 - 2 * long(L);
			out << fixed <<  x  << endl;
		} else {
			out << fixed <<  2*L*(long (x)-1) + 1  << endl;
		}
		out << 2*L*long (sa[i]) << endl;
		//out << genome.substr(sa[i],L-1) << "\t" << 2*L*(sa[i]-1) + 1  << endl;
		//out << genome.substr(sa[i],L-1) << "\t" << 2*L*(sa[i]-1) + 1  << endl;
	}
}

