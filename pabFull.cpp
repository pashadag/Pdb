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

typedef set< pair< long, int> > adjList;
typedef vector< adjList > adjType;

class position_id{
	public:
		position_id(long p, long n1, long n2) : position(p), nodeID(n1), nodeID2(n2) {};
		position_id(long p) : position(p), nodeID(-1), nodeID2(-2) {};
		position_id() : position(-1), nodeID(-1), nodeID2(-2) {};
		long position;
		long nodeID ;   // id of first mate
		long nodeID2;   // id of second mate  -2 means uninitialize, -1 means not found in DB graph
};


class doffset {
	public:
		long fwd_node;
		int fwd_offset;
		long bwd_node;
		int bwd_offset;
};



int numthreads;
char* allReadConcatenated ; 
int k;
int L ;
int delta;
long numCollapsedVertices=0;
long numNextNodes;
position_id * info;
long infosize;
vector<doffset> ep;
adjType colAdj;
long * nextNodes;
long * prevNodes;

class block {
public:
	long start;
	long end;
};

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


void ThreadOverInfo(void * (*func)(void *), int nthreads) {

	pthread_t thread[nthreads];
	block    params[nthreads];
	long blocksize = (infosize / nthreads) + 1;
	long startblock = 0;
	for (int i = 0; i < nthreads; i++) {
		long endblock = min(blocksize*i + blocksize, long(infosize));
		while (endblock != infosize && info[endblock].nodeID == info[endblock-1].nodeID) endblock++; //align the block
		params[i].start = startblock;
		params[i].end   = endblock;
		pthread_create(&thread[i], NULL, func, (void *) &params[i]);
		startblock = endblock;
	}
	for (int i = 0; i < nthreads; i++) {
		pthread_join(thread[i], NULL);
	}
	return;

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


void FillInKmers(const char* fileName)
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
	long counter = 0;
	while (file >> x) {
		position_id kdmer(x);
		info[counter++] = kdmer;
	}
	assert (infosize==counter);
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

position_id * my_lower_bound(position_id &kmer)
{
    position_id * info_iter = lower_bound(&info[0], &info[infosize], kmer, compareFunctor());

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
long FindNextNode(set<long> &positionSet) 
{
    set<long> nodes ;
    for(set<long>::iterator it = positionSet.begin() ; it!= positionSet.end(); ++it)
    {
        //check the position to see whether it's the last kd mers of the ldmers
        long currentPosition = *it ; 
        if( ((currentPosition % L) - (L- k)) == 0 )
            continue;

        position_id kmer(currentPosition + 1); 
        position_id * info_iter = my_lower_bound(kmer );
        if(info_iter == &info[infosize] )
        {
            cerr<<"Some errors in binary search \n"<<endl;
            exit(1);
        }
		if (info_iter->nodeID2 == -1) 
			continue;
        nodes.insert(info_iter->nodeID);
        if(nodes.size() > 1)
            return -1;

    }
    if(nodes.size() == 0)
        return -1;

    return *(nodes.begin());
}

long FindPrevNode(set<long> &positionSet)
{
	set<long> nodes ;
	for(set<long>::iterator it = positionSet.begin() ; it!= positionSet.end(); ++it)
	{
		//check the position to see whether it's the last kd mers of the ldmers
		long currentPosition = *it ; 
		if( (currentPosition % L) == 0 )
			continue;

		position_id kmer(currentPosition - 1);
		position_id * info_iter = my_lower_bound(kmer);
		if(info_iter == &info[infosize] )
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




void *GetSimpleGraphStructure(void * params) {
	long start = ((block *) params)->start;
	long end   = ((block *) params)->end;
	for(long i = start; i < end; i++) {
		set<long> positionSet;
		long currentNodeID = info[i].nodeID;
		if (info[i].nodeID2 == -1) {
			assert (i + 1 == infosize || info[i+1].nodeID != currentNodeID);
			nextNodes[currentNodeID] = -1;
		} else {
			while(i != infosize && info[i].nodeID == currentNodeID) {
				positionSet.insert(info[i].position);
				i++;
			}
			i--;
			long nextNode = FindNextNode(positionSet);
			nextNodes[currentNodeID] = nextNode; 
			long prevNode = FindPrevNode(positionSet); 
			prevNodes[currentNodeID] = prevNode; 
		}
	}


	pthread_exit(NULL);
}


/*
//Paul's version
set<long> GetSourcesNodes()
{
	set<long> sourceNodes ;

	for(long it = 0; it < infosize; it++) {
		long curNodeID = info[it].nodeID;
		if (nextNodes[curNodeID] == -1) {
			if (info[it].nodeID2 == -1) {
				assert (it == infosize -1 || curNodeID != info[it+1].nodeID);
			} else {
				while(it!= infosize &&  info[it].nodeID == curNodeID) {
					long curPos = info[it].position;
					if( ((curPos % L) - (L- k)) != 0 ) {
						position_id kmer(curPos + 1);
						position_id * info_iter = my_lower_bound(kmer);
						if(info_iter == &info[infosize] ) {
							cerr<<"Some errors in binary search \n"<<endl;
							exit(1);
						}
						if (info_iter->nodeID2 != -1) {
							long nextNode = info_iter->nodeID;
							sourceNodes.insert(nextNode);
						}
					}
					it++;
				}
				it--;
			}
		}
	}

	return sourceNodes ;
}

*/

set<long> GetSourcesNodes()
{
	set<long> sourceNodes ;
	set<long> nodesOutGoing ;
	set<long> nodesInGoing;
	long counter = 0;
	for(int i = 0; i < numNextNodes; i++) {
		if(nextNodes[i] != -1){
			nodesOutGoing.insert(counter);
			nodesInGoing.insert(nextNodes[i]);
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


//   simple and naive algorithm. TODO put back chasing   
vector<long> FindPath(long source)
{
	vector<long> path ;
	long currentNode = source;
	while(nextNodes[currentNode] != -1 && path.size() <= numNextNodes)
	{
		path.push_back(currentNode);
		currentNode = nextNodes[currentNode];
	} 
	path.push_back(currentNode);
	if(path.size() > numNextNodes)
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
void FindAndPrintContig(vector<long> &path,ostream &stream )
{
	ostringstream oMsg;
	//find the first path: 
	if(path.size() == 0)
		return ;
	position_id kmer(-1, path[0], -2);
	position_id * it = lower_bound(&info[0], &info[infosize], kmer, compareFunctorByNodeID());
	if(it == &info[infosize])
	{
		cerr<<"Some error in binary search for vertexID \n";
	}
	oMsg.write(&allReadConcatenated[it->position] ,k);
	for(long i = 1 ; i < path.size() ; i++)
	{
		kmer.nodeID = path[i];
		it = lower_bound(&info[0], &info[infosize], kmer, compareFunctorByNodeID());
		assert(kmer.nodeID == it->nodeID);
		oMsg.put( allReadConcatenated[(it->position) + k-1]);
	}
	//stream << oMsg.str().substr(0, oMsg.str().length() - k);
	stream << oMsg.str();
	stream << endl;

}



void GetComplexGraphStructure(adjType & adj, const vector<doffset> &ep ) {
	for(long i = 0; i < infosize; i++) {
		long curNodeID = info[i].nodeID;
		if (ep[curNodeID].fwd_offset != 0) //skip chains
			continue;

		adjList next;
		while(i != infosize &&  info[i].nodeID == curNodeID)
		{
			long curPos = info[i].position; 
			if( ((curPos % L) - (L- k)) != 0 ) {
				position_id kmer(curPos + 1);
				position_id * info_iter = my_lower_bound(kmer );
				if(info_iter == &info[infosize] ) {
					cerr<<"Some errors in binary search \n"<<endl;
					exit(1);
				}
				doffset x = ep[info_iter->nodeID];
				next.insert(make_pair(x.fwd_node, x.fwd_offset + 1));
			}
			i++;
		}
		i--;
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
		bool operator< (const KeyValPair & p) const { return value < p.value; }
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

void *PartitionNL (void * params) {// startblock, long endblock) { 
	//the nodes whose right mates don't map anywhere in the de bruijn graph become their own partition.
	//however, we later pretend that these vertices don't actually exist

	long globalstart = ((block *) params)->start;
	long globalend = ((block *) params)->end;
	long startBlockIndex = globalstart;
	cerr << "PartitionNL called with start = " << startBlockIndex << " and end = " << globalend << endl;
	long endBlockIndex = startBlockIndex;

	long vertexID = 0 ; 
	while(startBlockIndex != globalend) {
		if ((vertexID % 100000) == 0) {
			cerr << (100*(startBlockIndex - globalstart)) / (globalend- globalstart+ 1) << "% of the way there.\tAt " <<  startBlockIndex << "\t of \t"<< globalend << endl; 
		}
		//Identify block of kdmers with identitical left kmers (endBlockIndex)
		char* firstBlock = &allReadConcatenated[info[endBlockIndex].position];
		char* current = &allReadConcatenated[info[endBlockIndex].position];
		while(strncmp(firstBlock, current, k) == 0)
		{
			endBlockIndex++;
			if (endBlockIndex == globalend) 
				break;
			current = &allReadConcatenated[info[endBlockIndex].position];
		}


		//initizlize unionfind 
		//and cluster things on the same collapsed edges
		unionFindClass unionclass(endBlockIndex - startBlockIndex); 
		vector< sourceClass> sources; //edge/offset/id triples 
		for(long i = startBlockIndex ; i < endBlockIndex  ; i++){
			unionclass.find_set(i-startBlockIndex);
			if (info[i].nodeID2 != -1) {
				doffset o = ep[info[i].nodeID2];
				sources.push_back( sourceClass(o.bwd_node, o.fwd_node, o.bwd_offset, i - startBlockIndex, o.fwd_offset + o.bwd_offset));
			}
		}
		sort(sources.begin(), sources.end());
		for (int i = 1; i < sources.size(); i++) {
			if ((sources[i].to == sources[i-1].to) && (sources[i].from == sources[i-1].from) && (sources[i].length == sources[i-1].length ) &&(sources[i].offset - sources[i-1].offset <= 2* delta)) 
				unionclass.unionn(sources[i-1].id, sources[i].id);

		}	

		vector< vector<int> > partitions ;
		int np = unionclass.num_classes();
		if (np > 1) { 
			partitions.clear();
			unionclass.get_classes(partitions);
			vector<int> dist(numCollapsedVertices,-1);
			vector<int> minBwdOff(np,100000000);
			vector<int> minFwdOff(np,100000000);
			for (int i = 0; i < np ; i++) { //calculate minOffsets for partitions
				if (info[startBlockIndex + partitions[i][0]].nodeID2 != -1) {
					for (int j = 0; j < partitions[i].size(); j++) { //find smallest offset
						int bwdOff = ep[info[startBlockIndex + partitions[i][j]].nodeID2].bwd_offset;
						int fwdOff = ep[info[startBlockIndex + partitions[i][j]].nodeID2].fwd_offset;
						minBwdOff[i] = min(minBwdOff[i], bwdOff);
						minFwdOff[i] = min(minFwdOff[i], fwdOff);
					}
				}
			}
			for (int i = 0; i < np ; i++) {
				if (info[startBlockIndex + partitions[i][0]].nodeID2 != -1) {
					assert (info[startBlockIndex + partitions[i][0]].nodeID2 >= 0);
					if (minFwdOff[i] <= 2* delta) {
						set<int> found;
						priority_queue< KeyValPair > starts;
						long ver1 = info[startBlockIndex + partitions[i][0]].nodeID2;
						starts.push(KeyValPair(ep[ver1].fwd_node, -1 * minFwdOff[i]));
						sp (dist, colAdj, starts, found, 2*delta);
						for (int j = 0; j < np ; j++) {
							if (i==j) continue;
							long nodeID2 = info[startBlockIndex + partitions[j][0]].nodeID2;
							if (nodeID2 != -1) {
								doffset o = ep[nodeID2];
								if (found.find(o.bwd_node) != found.end()) {
									assert (minBwdOff[j] <= o.bwd_offset);
									if (dist[o.bwd_node] + minBwdOff[j] <= 2 * delta) {
										unionclass.unionn(partitions[i][0],partitions[j][0]);
									}
								}
							}
						}
						for (set<int>::iterator it = found.begin(); it != found.end(); it++) dist[*it] = -1; //reset dist for future calls
					}
				}
			}
		} 

		partitions.clear();
		unionclass.get_classes(partitions);
		for(vector<vector<int> >::iterator p_iter = partitions.begin() ; p_iter != partitions.end() ; ++p_iter) {

			for(vector<int>::iterator it = p_iter->begin() ; it != p_iter->end() ; ++it) {
				info[*it + startBlockIndex ].nodeID = vertexID ; 
				//this invalidats the nodeID2, but we will not be needing to use them from now on.  so be carefuly to not actually use them
				//infofile << info[(*it) + startBlockIndex ].position << " " <<  info[(*it) + startBlockIndex ].nodeID << endl;
			}
			vertexID++;
		}

		startBlockIndex = endBlockIndex ;

	}
	stable_sort(&info[globalstart], &info[globalend], compareFunctorByNodeID());
	pthread_exit(NULL);

}

void PartitionNLpar() {

	cerr << "Info size = " << infosize << endl;;
	ThreadOverInfo(PartitionNL, numthreads);

	//clean up vertexids
	long curid = 0;
	long lastlocid = 0;
	for (long i = 1; i < infosize; i++) {
		long curlocid = info[i].nodeID;
		if (curlocid != lastlocid){
			curid++;
		} 
		info[i].nodeID = curid;
		lastlocid = curlocid;
	}

	return;

}


void* FindNodeID2(void * params) {
	long start = ((pair<long,long> *) params)->first;
	long end = ((pair<long,long> *) params)->second;
	long notfound = 0;
	for (long i = start; i < end; i++) {
		position_id kmer(info[i].position + L);
		position_id * info_iter = lower_bound(&info[0], &info[infosize], position_id(info[i].position + L), compareFunctor());
		if ((info_iter == &info[infosize]) || !compeqk(kmer, *info_iter) ) {
			notfound++;
			info[i].nodeID2 = -1;
		} else {
			info[i].nodeID2 = info_iter->nodeID;
		}
	}
	cerr << "notfound\t" << notfound << endl;

	pthread_exit(NULL);

}

void dumpGraph(string filename, string name) {
	ofstream out;
	open_file(out, filename);
	out << "digraph " << name << " {\n";


	for(long i = 0; i < infosize; i++) {
		long curNodeID = info[i].nodeID;


		while(i != infosize &&  info[i].nodeID == curNodeID)
		{
			long curPos = info[i].position; 
			if( ((curPos % L) - (L- k)) != 0 ) {
				position_id kmer(curPos + 1);
				position_id * info_iter = my_lower_bound(kmer );
				if(info_iter == &info[infosize] ) {
					cerr<<"Some errors in binary search \n"<<endl;
					exit(1);
				}
				//doffset x = ep[info_iter->nodeID]; next.insert(make_pair(x.fwd_node, x.fwd_offset + 1));
				string label;
				for (int offset = 0; offset < min(k+1, 5); offset++) label += allReadConcatenated[info[i].position + offset];
				out << curNodeID << " -> " << info_iter->nodeID << " [label=" << label << "];\n";
			}
			i++;
		}
		i--;
	}

	out << "}\n";
	out.close();
	return ; 
}



bool fastmode = false;
bool savemode = false;
bool superfastmode = false;
int main(int argc, char* argv[])
{

	if(argc < 7)
	{
		cerr<<"pabNLpal base L k delta threads infosize [save,fast, superfast]"<<endl;
		exit(1);
	}
	if (argc >7) {
		string sbuf = argv[7];
		if (sbuf == "fast") {
			cerr << "fastmode\n";
			fastmode = true;
		} else if (sbuf == "save") {
			cerr << "savemode\n";
			savemode = true;
		} else if (sbuf == "superfast") {
			cerr << "superfastmode\n";
			superfastmode = true;
		} else {
			cerr << "couldn't figure out option " << argv[7] << endl;
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
	numthreads = atoi(argv[5]);
	infosize = atol(argv[6]);
	info = new position_id[infosize];
	assert (numthreads >=1 );

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

	if (!(fastmode||superfastmode)) {
		cerr<<"Reading kmer file...\n";
		FillInKmers(inputKmersFileName.c_str());
		cerr<<"Info size "<<infosize<<endl;

		cerr << "Fill in nodeIDs in info in the de Bruijn way.\n";
		vector<long> node2info;
		info[0].nodeID = 0;
		long curId = 0;
		node2info.push_back(0);
		for (long i = 1; i < infosize; i++) {
			if (!compeqk(info[i-1], info[i])) {
				curId++;
				node2info.push_back(i);
			}
			info[i].nodeID = curId;
		}
		dumpGraph(base + ".db.dot", "DB");

		numNextNodes = info[infosize-1].nodeID + 1;
		nextNodes = new long[numNextNodes];
		prevNodes = new long[numNextNodes];
		for (long i = 0; i < numNextNodes; i++) {
			nextNodes[i] = -1;
			prevNodes[i] = -1;
		}

		cerr << "Filling in nextNodes and prevNodes.\n";
		ThreadOverInfo(GetSimpleGraphStructure, numthreads);

		cerr<< "Fill in nodeID2 in info, corresponding to the nodeId of the right mate.\n";
		ThreadOverInfo(FindNodeID2, numthreads);

		cerr << "Collapsing de Bruijn graph (filling in ep).\n";
		ep.resize(numNextNodes);
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
		GetComplexGraphStructure(colAdj, ep);
		assert(numCollapsedVertices == colAdj.size());


		if (savemode) {
			ofstream file ;
			cerr << "Saving files....\n";

			open_file(file, base + ".adj");
			for (int i =0 ; i < colAdj.size(); i++) {
				for (adjList::iterator it = colAdj[i].begin(); it != colAdj[i].end(); it++) {
					file << it->first << "\t" << it->second << "\t";
				}
				file << endl;

			}
			file.close();

			open_file(file, base + ".dbinfo");
			for(long i = 0 ; i < infosize; i++) {
				file << info[i].position << "\t" << info[i].nodeID << "\t" << info[i].nodeID2 << endl;  

			}
			file.close();
			open_file(file, base + ".ep");
			for (long i = 0; i < ep.size(); i++) file << ep[i].bwd_node << "\t" << ep[i].bwd_offset << "\t" << ep[i].fwd_node << "\t" << ep[i].fwd_offset << endl;
			file.close();
		}
	}  
	if (fastmode) {
		ifstream file ;
		cerr << "Reading in dbinfo...\n";
		open_file(file, base + ".dbinfo");
		position_id p;
		long counter = 0;
		while (file >> p.position >> p.nodeID >> p.nodeID2 ) info[counter++] = p;
		assert (infosize == counter);
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
	if (superfastmode) {
		ifstream file ;
		cerr << "Reading in info...\n";
		open_file(file, base + ".info");
		position_id p;
		long counter = 0;
		while (file >> p.position >> p.nodeID >> p.nodeID2 ) info[counter++] = p;
		assert (infosize == counter);
		file.close();

	} else {
		cerr<<"Start partioning...";
		PartitionNLpar();
		ofstream infofile;
		open_file(infofile, base + ".info");
		for(long i = 0; i < infosize; i++) infofile << info[i].position << "\t" << info[i].nodeID << "\t" << info[i].nodeID2 <<endl;  
		infofile.close();
		cerr<<" end partioning.\n";
		dumpGraph(base + ".pdb.dot", "DB");
		ep.clear(); //this is invalid, don't even think of using it anymore.  you changed all the vertexIDs in Partition, remember?
		colAdj.clear(); //we also don't use this anymore, so forget about it.
	}
	

	cerr << "Finding chains...";
	numNextNodes = info[infosize-1].nodeID + 1;
	delete nextNodes;
	delete prevNodes;
	nextNodes = new long[numNextNodes];
	prevNodes = new long[numNextNodes];
	for (long i = 0; i < numNextNodes; i++) {
		nextNodes[i] = -1;
		prevNodes[i] = -1;
	}
	ThreadOverInfo(GetSimpleGraphStructure,numthreads);
	set<long> sourceNodes = GetSourcesNodes();

	cerr<<"end finding chains.\nPrinting contigs...\n";
	ofstream outFile, detailOutFile, stringOutFile ;
	outFile.open((base + ".contigs").c_str());
	detailOutFile.open((base + ".contigsdetail").c_str());
	stringOutFile.open((base + ".contigsString").c_str());

	outFile<< "L="<< L <<  "\t" << "k=" << k << "\tDel=" << delta << endl;
	long counter = 0;
	for(set<long>::iterator it = sourceNodes.begin(); it!= sourceNodes.end(); it++)
	{
		vector<long> path = FindPath(*it);
		if (path.size() > k) {
			PrintVector(path, detailOutFile);
			outFile << (path.size() )<<endl;
			//outFile << (path.size() + k)<<endl;
			stringOutFile << ">" << counter++ << endl;
			FindAndPrintContig(path, stringOutFile);
		}
	}
	outFile.close();
	detailOutFile.close();
	stringOutFile.close();
	delete info;
	delete nextNodes;
	delete prevNodes;
	return 0;

}

