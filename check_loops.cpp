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

int main(int argc, char* argv[])
{

	int lineNum=0;
	while (getline(cin, sbuf)) {
		istringstream line(sbuf);
		int num; int dummy;
		while (line >> num >> dummy) {
			if (num == lineNum) {
				cout << "Found loop on line " << lineNum << ": " << sbuf << endl;
			}
		}
		lineNum++;
	}

	return 0;
}

