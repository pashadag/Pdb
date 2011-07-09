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


class position_id{
	public:
		position_id(long p, long n1, long n2) : position(p), nodeID(n1), nodeID2(n2) {};
		position_id(long p) : position(p), nodeID(-1), nodeID2(-2) {};
		position_id() : position(-1), nodeID(-1), nodeID2(-2) {};
		long position;
		long nodeID ;   // id of first mate
		long nodeID2;   // id of second mate  -2 means uninitialize, -1 means not found in DB graph
};




int k;
int L;
int delta;
int numthreads;
char* allReadConcatenated ; 

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
void ConcatenateAllReads(char* array, ifstream& file)
{
	int bufsize = 200000000;
	//char buff[bufsize];
	char* buff = new char[bufsize];
	long counter = 0;
	cout<<"begin read"<<endl;
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
	cout<<"end filter "<<endl;
	return ;
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

bool compltgen(int y, string s) {
    int i = 0;
    while ((i < s.length()) && (s[i] == genome[y+i])) i++;
    if (i == s.length()) return false;
    return  genome[y+i] < s[i];
}


void find_string(vector<int> & hits, string contig) {
	vector<int>::iterator res = lower_bound(sa.begin(), sa.end(), contig, compltgen);
	if ((res == sa.end()) || (!compeqgen(*res, contig))) {
		cerr << "string not found error.\n";
		for (int i = 0; i < info.size(); i++) cerr << info[i].position << "\t" << info[i].nodeID << endl;
		exit(1);
	} 
	while ( (res != sa.end()) && (compeqgen(*res, contig))) {
		hits.push_back(*res);
		res++;
	}
	return;
}

void validate(vector<position_id> & info) {

	//build groups
	vector<vector< position_id> > groups;
	position_id last = info[0];
	info.push_back();
	info[0].push_back(last)
	for (int i = 1; i < info.size; i++) {
		if (last.nodeID != info[i].nodeID) {
			info.push_back();
		}
		info.back().push_back(info[i]);
	}

	//get mappings
	vector<vector<int> > locs;
	for (int i = 0; i < groups.size(); i++) {
		vector<int> hits;
		for (int j = 0; j < groups[i].size; j++) {
			string contig = 
			find_string(hits, contig);
		}
		locs.push_back(hits);
	}
	//check if mappings match up

	return;

}



bool compeq(const position_id first , const position_id second) {
	return (0==strncmp(&allReadConcatenated[first.position], &allReadConcatenated[second.position], k));
}


vector<int> sa;
string genome;

bool compeqgen(int y, string s ) {
    int i = 0;
    while ((i < s.length()) && (s[i] == genome[y+i])) i++;
    if (i == s.length()) return true;
    return false ;
}


int main(int argc, char* argv[]) {

	if(argc < 6)
	{
		cerr<<"validate base L k delta genome" << endl;
		exit(1);
	}
	string base = argv[1];
	genome = read_genome(base);
	genome += genome;
	L = atoi(argv[2]);
	k = atoi(argv[3]);
	delta = atoi(argv[4]);
	string genome = read_genome(argv[5] + ".fa");

	ifstream safile;
	open_file(safile, argv[5] + ".sa");

	int ibuf;
	while (safile >> ibuf) sa.push_back(ibuf);
	safile.close();


	//numthreads = atoi(argv[5]);

	//assert (numthreads >=1 );

	ifstream inputReadFile ;
	open_file(inputReadFile, base + ".readl");
	long numberOfLines = GetNumberOfLine(inputReadFile);
	allReadConcatenated = new char[numberOfLines*2*L];
	inputReadFile.close();
	ifstream newInputReadFile ;
	open_file( newInputReadFile, base + ".readl");  
	ConcatenateAllReads(allReadConcatenated, newInputReadFile);
	newInputReadFile.close();


	

	ifstream file;
	open_file(file, base + ".info");
	position_id p, last;
	vector<position_id> info;

	file >> p.position >> p.nodeID >> p.nodeID2;
	last = p;
	info.push_back(p);

	while (file >> p.position >> p.nodeID >> p.nodeID2) {
		if (compeq(p, last)) {
			info.push_back(p);
		} else {
			cout << "read in a block with first = " << info.at(0).position << " and last = " << info.back().position << endl;
			validate(info); //process info
			info.clear();
			info.push_back(p);
		}
		last = p;
	}
		
	file.close();


	return 0;

}

