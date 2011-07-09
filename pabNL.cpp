#include<iostream>
#include<string.h>
#include<vector>
#include<set>
#include<deque>
#include<set>
#include<fstream>
#include<cstdlib>
#include"union.h"
#define _DEBUG_ 
using namespace std ;

typedef set< pair< long, int> > adjList;
typedef vector< adjList > adjType;

char* allReadConcatenated ; 
int k;
int L ;
int delta;
long numCollapsedVertices=0;

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
	cerr<<"begin read"<<endl;
    while(!file.eof())
    {
        file.read(buff,bufsize);
        int numberofCharRead = file.gcount();
       
       for(int i = 0 ; i < numberofCharRead ; i++)
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

class position_id{
public:
	position_id(long p, long n1, long n2) : position(p), nodeID(n1), nodeID2(n2) {};
	position_id(long p) : position(p), nodeID(-1), nodeID2(-1) {};
	position_id() : position(-1), nodeID(-1), nodeID2(-1) {};
    long position;
    long nodeID ;   // id of first mate
	long nodeID2;   // id of second mate
	/*bool operator< (const position_id & second) const {
		if( strncmp(&allReadConcatenated[position], &allReadConcatenated[second.position], k) < 0)
			return true;
		return false;
	}*/


};



class doffset {
public:
	long fwd_node;
	int fwd_offset;
	long bwd_node;
	int bwd_offset;
};

bool compeqk(const position_id first , const position_id second) {
	return (0==strncmp(&allReadConcatenated[first.position], &allReadConcatenated[second.position], k));
}

struct compareFunctor
{
    //This will find k-dmer sorted by the first part
    bool operator()(const position_id first , const position_id second) const
    {
        if( strncmp(&allReadConcatenated[first.position], &allReadConcatenated[second.position], k) < 0)
            return true;
        return false;

    }
};


void FillInKmers(deque<position_id> &info, const char* fileName, int k)
{
    ifstream file ;
    file.open(fileName);
    //long distance = 2* k +2;
    if(!file.is_open())
    {
        cerr<<"Error when openning file "<< fileName<<endl;
        exit(1);
    }
	long x;
	while (file >> x) {
		position_id kdmer(x);
		info.push_back(kdmer);
	}
	file.close();

}

void WriteInfo(deque<position_id> & info, const char* fileName)
{
	ofstream file ;
	file.open(fileName);
	if(!file.is_open())
	{
        cerr<<"Error when opening the file "<<fileName << endl;
        exit(1);
    }
    for(deque<position_id>::iterator it = info.begin() ; it!= info.end() ; ++it)
    {
        file << it->position << " " << it->nodeID <<endl;  

    }
    file.close();
}
struct compareFunctorByNodeID
{
    bool operator()(const position_id first , const position_id second) const
    {
        if(first.nodeID < second.nodeID)
            return true;
        return false ;

    }
};

deque<position_id>::iterator my_lower_bound(deque<position_id> &info, position_id &kmer)
{
    deque<position_id>::iterator info_iter = lower_bound(info.begin(), info.end(), kmer, compareFunctor());

    while(true)
    {
        if(strncmp(&allReadConcatenated[ info_iter->position + L ], &allReadConcatenated[kmer.position + L ],k) == 0)
            break;
        if(allReadConcatenated[info_iter->position] != allReadConcatenated[kmer.position])
        {
            cerr<<"Error in my_lower_bound \n";
			assert(0);
            exit(1);
        } 
        info_iter++;
    }
    return info_iter ; 

}
long FindNextNode(set<long> &positionSet, deque<position_id> & info )
{
    set<long> nextNodes ;
    for(set<long>::iterator it = positionSet.begin() ; it!= positionSet.end(); ++it)
    {
        //check the position to see whether it's the last kd mers of the ldmers
        long currentPosition = *it ; 
        if( ((currentPosition % L) - (L- k)) == 0 )
            continue;

        position_id kmer(currentPosition + 1); 
        deque<position_id>::iterator info_iter = my_lower_bound(info, kmer );
        if(info_iter == info.end() )
        {
            cerr<<"Some errors in binary search \n"<<endl;
            exit(1);
        }
        nextNodes.insert(info_iter->nodeID);
        if(nextNodes.size() > 1)
            return -1;

    }
    if(nextNodes.size() == 0)
        return -1;

    return *(nextNodes.begin());
}

long FindPrevNode(set<long> &positionSet, deque<position_id> & info )
{
	set<long> nodes ;
	for(set<long>::iterator it = positionSet.begin() ; it!= positionSet.end(); ++it)
	{
		//check the position to see whether it's the last kd mers of the ldmers
		long currentPosition = *it ; 
		if( (currentPosition % L) == 0 )
			continue;

		position_id kmer(currentPosition - 1);
		deque<position_id>::iterator info_iter = my_lower_bound(info, kmer );
		if(info_iter == info.end() )
		{
			cerr<<"Some errors in binary search \n"<<endl;
			exit(1);
		}
		nodes.insert(info_iter->nodeID);
		if(nodes.size() > 1)
			return -1;

	}
	if(nodes.size() == 0)
		return -1;

	return *(nodes.begin());
}



void GetSimpleGraphStructure(vector<long> &nextNodes, vector<long> &prevNodes, deque<position_id> &info)
{
    for(deque<position_id>::iterator it = info.begin() ; it != info.end(); it++)
    {
        set<long> positionSet;
        long currentNodeID = it->nodeID;
        while(it!= info.end() &&  it->nodeID == currentNodeID)
        {
            positionSet.insert(it->position);
            it++;
        }
        it--;
        long nextNode = FindNextNode(positionSet, info);
		long prevNode = FindPrevNode(positionSet, info);
        nextNodes[currentNodeID] = nextNode; 
		prevNodes[currentNodeID] = prevNode; 
    }


    return ; 
}

/*

//chains only version
set<long> GetSourcesNodes(vector<long>  & nextNodes, vector<long>  & prevNodes)
{
	set<long> sourceNodes ;
	for(long i = 0; i < nextNodes.size(); i++) { 
		if((nextNodes[i] != -1) && (prevNodes[i] == -1)) {
			sourceNodes.insert(i);
		}
	}
	return sourceNodes ;
}

//Son's version
set<long> GetSourcesNodes(vector<long>  & nextNodes)
{
	set<long> sourceNodes ;
	set<long> nodesOutGoing ;
	set<long> nodesInGoing;
	long counter = 0;
	for(vector<long>::iterator it = nextNodes.begin() ; it!= nextNodes.end() ; it++)
	{
		if((*it) != -1){
			nodesOutGoing.insert(counter);
			nodesInGoing.insert(*it);
		}

		counter++;
	}
	for(set<long>::iterator it = nodesOutGoing.begin() ; it!= nodesOutGoing.end(); it++)
	{
		if(nodesInGoing.find(*it) == nodesInGoing.end())
			sourceNodes.insert(*it);
	}
	return sourceNodes ;
}

 */

//Paul's version
set<long> GetSourcesNodes(deque<position_id> & info, vector<long>  & nextNodes)
{
    set<long> sourceNodes ;
	
	for(deque<position_id>::iterator it = info.begin() ; it != info.end(); it++)
	{
		long curNodeID = it->nodeID;
		if (nextNodes[curNodeID] == -1) {
			while(it!= info.end() &&  it->nodeID == curNodeID) {
				long curPos = it->position;
				if( ((curPos % L) - (L- k)) != 0 ) {
					position_id kmer(curPos + 1);
					deque<position_id>::iterator info_iter = my_lower_bound(info, kmer );
					if(info_iter == info.end() ) {
						cerr<<"Some errors in binary search \n"<<endl;
						exit(1);
					}
					long nextNode = info_iter->nodeID;
					sourceNodes.insert(nextNode);
				}
				it++;
			}
			it--;
		}
	}

	return sourceNodes ;
}

vector<long> GetNonLoop(vector<long> &path)
{
	set<long> nodeSet;
	vector<long> nonloopPath;
	for(long i = 0 ; i < path.size() ; i ++)
	{
		if(nodeSet.find(path[i]) != nodeSet.end())
		{
			return nonloopPath;
		}
		nodeSet.insert(path[i]);
		nonloopPath.push_back(path[i]);

	}
	return nonloopPath; 
}

/*
//stop at chains
vector<long> FindPath(vector<long> &nextNodes, vector<long> &prevNodes, long source)
{
	vector<long> path ;
	long currentNode = source;
	while(nextNodes[currentNode] != -1 && path.size() <= nextNodes.size()) 
	{
		path.push_back(currentNode);
		currentNode = nextNodes[currentNode];
		if (prevNodes[nextNodes[currentNode]] != -1 ) continue;
	} 
	path.push_back(currentNode);
	if(path.size() > nextNodes.size())
	{
		cerr << "Getting NonLoop path\n";
		return GetNonLoop(path);
	}
	return path;

}

*/

//   simple and naive algorithm. TODO put back chasing   
vector<long> FindPath(vector<long> &nextNodes, long source)
{
	vector<long> path ;
	long currentNode = source;
	while(nextNodes[currentNode] != -1 && path.size() <= nextNodes.size())
	{
		path.push_back(currentNode);
		currentNode = nextNodes[currentNode];
	} 
	path.push_back(currentNode);
	if(path.size() > nextNodes.size())
	{
		cerr << "Getting NonLoop path\n";
		return GetNonLoop(path);
	}
	return path;

}

void PrintVector(vector<long> &v, ostream &stream)
{
	for(long i = 0 ; i < v.size() ; i++)
		stream<<v[i]<<" ";
	stream<<endl;
}
void PrintSet (set<long> &s)
{
	for(set<long>::iterator it = s.begin() ; it!= s.end() ;it++)
	{
		cout<<*it<<" ";
	}
	cout<<endl;
}
void FindAndPrintContig(vector<long> &path,ostream &stream, deque<position_id> &info )
{
	ostringstream oMsg;
    //find the first path: 
    if(path.size() == 0)
        return ;
    position_id kmer(-1, path[0], -1);
    deque<position_id>::iterator it = lower_bound(info.begin(), info.end(),kmer, compareFunctorByNodeID());
    if(it == info.end())
    {
        cerr<<"Some error in binary search for vertexID \n";
    }
    oMsg.write(&allReadConcatenated[it->position] ,k);
    for(long i = 1 ; i < path.size() ; i++)
    {
        kmer.nodeID = path[i];
        it = lower_bound(info.begin(), info.end(), kmer, compareFunctorByNodeID());
		assert(kmer.nodeID == it->nodeID);
        oMsg.put( allReadConcatenated[(it->position) + k-1]);
    }
	stream << oMsg.str().substr(0, oMsg.str().length() - k);
    stream << endl;

}



void GetComplexGraphStructure(adjType & adj, deque<position_id> &info, const vector<doffset> &ep )
	
{
	for(deque<position_id>::iterator it = info.begin() ; it != info.end(); it++)
	{
		long curNodeID = it->nodeID;
		if (ep[curNodeID].fwd_offset != 0) //skip chains
			continue;
		
		adjList next;
		while(it!= info.end() &&  it->nodeID == curNodeID)
		{
			long curPos = it->position; 
			if( ((curPos % L) - (L- k)) != 0 ) {
				position_id kmer(curPos + 1);
				deque<position_id>::iterator info_iter = my_lower_bound(info, kmer );
				if(info_iter == info.end() ) {
					cerr<<"Some errors in binary search \n"<<endl;
					exit(1);
				}
				doffset x = ep[info_iter->nodeID];
				next.insert(make_pair(x.fwd_node, x.fwd_offset + 1));
			}
			it++;
		}
		it--;
		//list< pair< long, int> > nextList;
		//nextList.assign(next.begin(), next.end());
		adj.push_back(next);
	}

	return ; 
}

class sourceClass {
public:
	sourceClass() {};
	sourceClass(int f, int t, int o, int i, int l) : from(f), to(t), offset(o), id(i), length(l) {}
	int from;
	int to;
	int offset;
	int id;
	int length;
	bool operator<( const sourceClass & s2) const  {
		if (from != s2.from) return from < s2.from;
		if (to != s2.to) return to < s2.to;
		if (length != s2.length) return length < s2.length;
		if (offset != s2.offset) return offset < s2.offset;
		return id < s2.id;
	}

};



class KeyValPair {
	public:
		int key;
		int value;
		KeyValPair() : key(-1), value(-1) { }
		KeyValPair(int x, int y) : key(x), value(y) {}
		bool operator< (const KeyValPair & p) const { return value > p.value; }
};

void sp(vector<int> & dist, adjType & graph, priority_queue< KeyValPair> & heap, set<int> &found, int maxDist) {

	// all nodes whithin maxDist of start are stored in distances along with their distances
	// dist has to be initialized to -1 and sizez of num vertices of the graph

	while (!heap.empty()) {
		KeyValPair top = heap.top();
		heap.pop();
		top.value = top.value * -1;
		if (dist[top.key] != -1) continue;
		found.insert(top.key);
		dist[top.key] = top.value;
		//cout << "marked " << top.key << " " << top.value << endl;

		for (adjList::iterator i = graph[top.key].begin(); i != graph[top.key].end(); i++) {
			int adjNext = i->first;
			int adjLen  = i->second;
			if ((dist[adjNext] == -1) && (top.value + adjLen <= maxDist)) {
				heap.push(KeyValPair(adjNext, -1 *(top.value + adjLen)));
			}
		}
	}
	return;
}

void PartitionNL (deque< position_id> &info, vector<doffset> & ep, adjType & colAdj, ofstream & infofile) {

	long simpleCases=0;
	long complexCases=0;

	deque<position_id>::iterator info_iter = info.begin();
	
	long startBlockIndex = 0 ;
	long vertexID = 0 ; 
	long endBlockIndex =0;
	cout<<"size "<<info.size()<<endl;

	//startBlockIndex = 55348958; endBlockIndex = startBlockIndex; info_iter = info.begin() + startBlockIndex;

	while(info_iter != info.end())
	{
		//Identify block of kdmers with identitical left kmers (endBlockIndex)
		char* firstBlock = &allReadConcatenated[info_iter->position];
		char* current = &allReadConcatenated[info_iter->position];
		while(strncmp(firstBlock, current, k) == 0)
		{
			info_iter++;
			endBlockIndex++;
			if(info_iter == info.end())
				break;
			current = &allReadConcatenated[info_iter->position];
		}


		//initizlize unionfind 
		//and cluster things on the same collapsed edges
		if (startBlockIndex == 55348958) {
			cerr << "lets do it!\n";
		}
		unionFindClass unionclass(endBlockIndex - startBlockIndex); 
		vector< sourceClass> sources; //edge/offset/id triples 
		for(long i = startBlockIndex ; i < endBlockIndex  ; i++){
			unionclass.find_set(i-startBlockIndex);
			doffset o = ep[info[i].nodeID2];
			sources.push_back( sourceClass(o.bwd_node, o.fwd_node, o.bwd_offset, i - startBlockIndex, o.fwd_offset + o.bwd_offset));
		}
		sort(sources.begin(), sources.end());
		for (int i = 1; i < sources.size(); i++) {
			if ((sources[i].to == sources[i-1].to) && (sources[i].from == sources[i-1].from) && (sources[i].length == sources[i-1].length ) &&(sources[i].offset - sources[i-1].offset <= 2* delta)) 
				unionclass.unionn(sources[i-1].id, sources[i].id);

		}	

		vector< vector<int> > partitions ;
		int np = unionclass.num_classes();
		if (np > 1) { 
			complexCases++;
			partitions.clear();
			unionclass.get_classes(partitions);
			vector<int> dist(numCollapsedVertices,-1);
			vector<int> minBwdOff(np,100000000);
			vector<int> minFwdOff(np,100000000);
			for (int i = 0; i < np ; i++) { //calculate minOffsets for partitions
				for (int j = 0; j < partitions[i].size(); j++) { //find smallest offset
					int bwdOff = ep[info[startBlockIndex + partitions[i][j]].nodeID2].bwd_offset;
					int fwdOff = ep[info[startBlockIndex + partitions[i][j]].nodeID2].fwd_offset;
					minBwdOff[i] = min(minBwdOff[i], bwdOff);
					minFwdOff[i] = min(minFwdOff[i], fwdOff);
				}
			}
			for (int i = 0; i < np ; i++) {

				if (minFwdOff[i] <= 2* delta) {
					set<int> found;
					priority_queue< KeyValPair > starts;
					long ver1 = info[startBlockIndex + partitions[i][0]].nodeID2;
					starts.push(KeyValPair(ep[ver1].fwd_node, -1 * minFwdOff[i]));
					sp (dist, colAdj, starts, found, 2*delta);
					for (int j = 0; j < np ; j++) {
						if (i==j) continue;
						doffset o = ep[info[startBlockIndex + partitions[j][0]].nodeID2];
						if (found.find(o.bwd_node) != found.end()) {
							assert (minBwdOff[j] <= o.bwd_offset);
							if (dist[o.bwd_node] + minBwdOff[j] <= 2 * delta) {
								unionclass.unionn(partitions[i][0],partitions[j][0]);
							}
						}
					}
					for (set<int>::iterator it = found.begin(); it != found.end(); it++) dist[*it] = -1; //reset dist for future calls
				}
			}
		} else {
			simpleCases++;
		}


		partitions.clear();
		unionclass.get_classes(partitions);
		for(vector<vector<int> >::iterator p_iter = partitions.begin() ; p_iter != partitions.end() ; ++p_iter) {

			for(vector<int>::iterator it = p_iter->begin() ; it != p_iter->end() ; ++it) {
				info[(*it) + startBlockIndex ].nodeID = vertexID ; 
				//this invalidats the nodeID2, but we will not be needing to use them from now on.  so be carefuly to not actually use them
				infofile << info[(*it) + startBlockIndex ].position << "\t" <<  info[(*it) + startBlockIndex ].nodeID << "\t" << info[(*it) + startBlockIndex ].nodeID2 << endl;
			}
			vertexID++;
		}
		stable_sort(info.begin() + startBlockIndex, info.begin() + endBlockIndex, compareFunctorByNodeID());

		startBlockIndex = endBlockIndex ;

	}
	cout  << "Complex cases = " << complexCases << "\tSimple cases = " << simpleCases << endl;

}


bool fastmode = false;
bool savemode = false;
int main(int argc, char* argv[])
{

	if(argc < 5)
	{
		cerr<<"pabNL base L k delta "<<endl;
		exit(1);
	}
	if (argc >5) {
		cerr << "argv[5] = " << argv[5] << endl;
		string sbuf = argv[5];
		if (sbuf == "fast") {
			cerr << "fastmode\n";
			fastmode = true;
		} else if (sbuf == "save") {
			cerr << "savemode\n";
			savemode = true;
		}
	}

	string base = argv[1];
	string inputReadsFileName(base);
	string inputKmersFileName(base);
	inputReadsFileName += ".readl";
	inputKmersFileName += ".readks";
	L = atoi(argv[2]) ;
	k = atoi(argv[3]);
	delta = atoi(argv[4]);
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
	//read the sorted KmersFileName
	deque<position_id> info;
	vector<long> nextNodes;
	vector<long> prevNodes;
	vector<doffset> ep;
	adjType colAdj;

	if (!fastmode) {
		cout<<"Reading kmer file...\n";
		FillInKmers(info, inputKmersFileName.c_str(), k);
		cout<<"Info size "<<info.size()<<endl;

		cerr << "Fill in nodeIDs in info in the de Bruijn way.\n";
		vector<long> node2info;
		info.at(0).nodeID=0;
		long curId = 0;
		node2info.push_back(0);
		for (long i = 1; i < info.size(); i++) {
			if (!compeqk(info[i-1], info[i])) {
				curId++;
				node2info.push_back(i);
			}
			info[i].nodeID = curId;
		}

		
		cerr<< "Fill in nodeID2 in info, corresponding to the nodeId of the right mate.\n";
		int rightmersnotfound=0;
		for (long i = 0; i < info.size(); i++) {
			position_id kmer(info[i].position + L);
			deque<position_id>::iterator info_iter = lower_bound(info.begin(), info.end(), position_id(info[i].position + L), compareFunctor());
			if ((info_iter == info.end()) || !compeqk(kmer, *info_iter) ) {
				rightmersnotfound++;
				//assert(0); exit(1);
			} else {
				info[i].nodeID2 = info_iter->nodeID;
			}
		}
		cerr << rightmersnotfound << " right mates not found.\n";

		



		nextNodes.resize(info.back().nodeID + 1, -1) ;
		prevNodes.resize(nextNodes.size(), -1) ;
		cerr << "Filling in nextNodes and prevNodes.\n";
		GetSimpleGraphStructure(nextNodes, prevNodes, info);

		cerr << "Collapsing de Bruijn graph (filling in ep).\n";
		ep.resize(nextNodes.size());
		for (long i = 0; i < ep.size(); i++) {
			ep[i].bwd_node = i;
			ep[i].fwd_node = i;
			ep[i].bwd_offset = 0;
			ep[i].fwd_offset = 0;
		}
		for (long i = 0; i < ep.size(); i++) {
			if ((nextNodes[i] != -1) && (prevNodes[i] != -1) ) { //chain vertex
				if ((nextNodes[prevNodes[i]] == -1) || (prevNodes[prevNodes[i]] == -1))  { // prev is not a chain vertex, so start of a forward chain
					int startVer = prevNodes[i];
					int offset = 1;
					int curVer = i;
					do {
						ep[curVer].bwd_node = startVer;
						ep[curVer].bwd_offset = offset++;
						curVer = nextNodes[curVer];
					} while ((nextNodes[curVer] != -1) && (prevNodes[curVer] != -1)); //while curVer is still a chain

				} 
				if ((nextNodes[nextNodes[i]] == -1) || (prevNodes[nextNodes[i]] == -1))  { // next is not a chain vertex, so start of a backward chain
					int startVer = nextNodes[i];
					int offset = 1;
					int curVer = i;
					do {
						ep[curVer].fwd_node = startVer;
						ep[curVer].fwd_offset = offset++;
						curVer = prevNodes[curVer];
					} while ((nextNodes[curVer] != -1) && (prevNodes[curVer] != -1));
				}

			}
		}

		vector<long> epnew2old;

		cerr << "Renumber vertices in ep.\n";
		for (long i = 0; i < ep.size(); i++) { 
			if (ep[i].fwd_offset == 0) {
				if (ep[i].bwd_offset != 0) {
					cerr << "Trouble in paradise.\n";
				}
				epnew2old.push_back(i);
				//assert (ep[i].bwd_offset == 0);
				ep[i].fwd_node = numCollapsedVertices;
				ep[i].bwd_node = numCollapsedVertices;
				numCollapsedVertices++;
			}
		}
		for (long i = 0; i < ep.size(); i++) { 
			if (ep[i].fwd_offset != 0) {
				assert (ep[i].bwd_offset != 0);
				assert (ep[i].fwd_offset != 0);
				ep[i].fwd_node = ep[ep[i].fwd_node].fwd_node;
				ep[i].bwd_node = ep[ep[i].bwd_node].fwd_node;
			}
		}

		cerr << "Filling in colAdj...\n";
		GetComplexGraphStructure(colAdj, info, ep);
		assert(numCollapsedVertices == colAdj.size());


		if (savemode) {
			ofstream file ;

			open_file(file, base + ".adj");
			for (int i =0 ; i < colAdj.size(); i++) {
				for (adjList::iterator it = colAdj[i].begin(); it != colAdj[i].end(); it++) {
					file << it->first << "\t" << it->second << "\t";
				}
				file << endl;

			}
			file.close();

			open_file(file, base + ".dbinfo");
			for(deque<position_id>::iterator it = info.begin() ; it!= info.end() ; ++it) {
				file << it->position << "\t" << it->nodeID << "\t" << it->nodeID2 << endl;  

			}
			file.close();
			open_file(file, base + ".ep");
			for (int i = 0; i < ep.size(); i++) file << ep[i].bwd_node << "\t" << ep[i].bwd_offset << "\t" << ep[i].fwd_node << "\t" << ep[i].fwd_offset << endl;
			file.close();
		}
	}  else { //if (fastmode)
		cerr << "Reading in dbinfo...\n";
		ifstream file ;
		open_file(file, base + ".dbinfo");
		position_id p;
		while (file >> p.position >> p.nodeID >> p.nodeID2 ) info.push_back(p);
		file.close();

		cerr << "Reading in ep...\n";
		open_file(file, base + ".ep");
		doffset d;
		while (file >> d.bwd_node >> d.bwd_offset >> d.fwd_node >> d.fwd_offset ) ep.push_back(d);
		file.close();

		cerr << "Reading in adj...\n";
		open_file(file, base + ".adj");
		string sbuf;
		while (getline(file, sbuf)) {
			istringstream line(sbuf);
			pair<long,int> element;
			adjList cur; //set of long int pairs
			while (line >> element.first >> element.second) {
				cur.insert(element);
			}
			colAdj.push_back(cur);
		}
		file.close();
		numCollapsedVertices = colAdj.size();
	}

	cout<<"Start partioning...";
	ofstream infofile;
	open_file(infofile, base + ".info");
	PartitionNL(info, ep, colAdj, infofile);
	infofile.close();
	cout<<"end partioning.\nFinding chains again...";

	ep.clear(); //this is invalid, don't even think of using it anymore.  you changed all the vertexIDs in Partition, remember?
	nextNodes.clear(); prevNodes.clear();
	nextNodes.resize(info.back().nodeID + 1, -1) ;
	prevNodes.resize(nextNodes.size(), -1) ;
	GetSimpleGraphStructure(nextNodes, prevNodes, info);
	//set<long> sourceNodes = GetSourcesNodes(nextNodes, prevNodes);
	set<long> sourceNodes = GetSourcesNodes(info, nextNodes);
#ifdef _DEBUG_
	//PrintVector(nextNodes,cout);
	//PrintSet(sourceNodes);

#endif    
	cout<<"end finding chains.\nPrinting contigs...";
	ofstream outFile, detailOutFile, stringOutFile ;
	outFile.open((base + ".contigs").c_str());
	detailOutFile.open((base + ".contigsdetail").c_str());
	stringOutFile.open((base + ".contigsString").c_str());

	outFile<< "L="<< L <<  "\t" << "k=" << k << "\tDel=" << delta << endl;
	long counter = 0;
	for(set<long>::iterator it = sourceNodes.begin(); it!= sourceNodes.end(); it++)
	{
		//vector<long> path = FindPath(nextNodes, prevNodes, *it);
		vector<long> path = FindPath(nextNodes, *it);
		if (path.size() > k) {
			PrintVector(path, detailOutFile);
			outFile << (path.size() )<<endl;
			//outFile << (path.size() + k)<<endl;
			stringOutFile << ">" << counter++ << endl;
			FindAndPrintContig(path, stringOutFile, info);
		}
	}
	outFile.close();
	detailOutFile.close();
	stringOutFile.close();
	return 0;

}

