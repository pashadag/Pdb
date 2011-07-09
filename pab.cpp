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

char* allReadConcatenated ; 
int k;
int L ;
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

void ConcatenateAllReads(char* array, ifstream& file)
{
	int bufsize = 200000000;
    //char buff[bufsize];
	char* buff = new char[bufsize];
    long counter = 0;
    while(!file.eof())
    {
        cerr<<"begin read"<<endl;
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
       cerr<<"end filter "<<endl;
    }
	delete buff;
    return ;
}
struct position_id{
    long position;
    long nodeID ; 
};
typedef struct position_id position_id;
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
    if(!file.is_open())
    {
        cerr<<"Error when openning file "<< fileName<<endl;
        exit(1);
    }
	long x;
	while (file >> x) {
		struct position_id kdmer = {x,-1};
		info.push_back(kdmer);
	}
	/*file.seekg(0, ios::beg);
	  bool first  = true;
	  while(!file.fail())
	  {
	  long a=-1;
	  long pos = file.tellg();
	  if(first)
	  {

	  file.seekg(distance -1 + pos);
	  first = false ;
	  }else
	  file.seekg(distance+ pos);

	  file >> a;
	  if(a == -1)
	  break;
	  struct position_id kdmer = {a,-1};
	  info.push_back(kdmer);
	  }
	 */
	file.close();

}
bool Overlap(long first, long second, int overlapThreshold)
{
	//first before second 
	char* firstPointer = &allReadConcatenated[first+ k - overlapThreshold ] ;
	for ( int i = 0 ; i < k - overlapThreshold +1 ; i++)
	{
		if(strncmp( (firstPointer -i), &allReadConcatenated[second] , overlapThreshold +i ) == 0)
			return true;
	}
	//first after second 
	char* secondPointer = &allReadConcatenated[ second + k - overlapThreshold];
	for(int i = 0 ; i < k - overlapThreshold +1 ; i++)
	{
		if( strncmp(secondPointer -i, & allReadConcatenated[first], overlapThreshold +i) == 0)
			return true;
	}
	return false; 

}
void Partition(deque< position_id> &info, int overlapThreshold)
{
	deque<position_id>::iterator info_iter = info.begin();
	long startBlockIndex = 0 ;
	long endBlockIndex =0;
	long vertexID = 0 ; 
	cout<<"size "<<info.size()<<endl;
	while(info_iter != info.end())
	{
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
		unionFindClass unionclass(endBlockIndex - startBlockIndex); 
		for(long i = startBlockIndex ; i < endBlockIndex  ; i++)
			unionclass.find_set(i-startBlockIndex);

		for( long i = startBlockIndex ; i < endBlockIndex ; i++)
		{
			for(long j = startBlockIndex +1 ; j < endBlockIndex ; j++)
			{
				if (j == i)
					continue;
				if(Overlap(info[i].position + L, info[j].position + L , overlapThreshold))
					unionclass.unionn(i-startBlockIndex,j-startBlockIndex);

			}
		}
		vector<vector<int> > partitions ;
		unionclass.get_classes(partitions);
		for(vector<vector<int> >::iterator p_iter = partitions.begin() ; p_iter != partitions.end() ; ++p_iter)
		{
			for(vector<int>::iterator it = p_iter->begin() ; it != p_iter->end() ; ++it)
			{
				info[(*it) + startBlockIndex ].nodeID = vertexID ; 
			}
			vertexID++;
			if (vertexID == 4440) {
				cerr << "point of interest\n";
			}
		}

		startBlockIndex = endBlockIndex ;

	}

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

        position_id kmer = {currentPosition + 1 , -1};
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


void GetSimpleGraphStructure(vector<long> &nextNodes, deque<position_id> &info)
{
    vector<long> sources ;
    vector<long> previousNodes(nextNodes.size(), -1);
    vector<long> terminateNodes ;
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
        nextNodes[currentNodeID] = nextNode; 
    }

    return ; 
}


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
   simple and naive algorithm. TODO put back chasing    */
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
    //find the first path: 
    if(path.size() == 0)
        return ;
    position_id kmer  = {-1, path[0]};
    deque<position_id>::iterator it = lower_bound(info.begin(), info.end(),kmer, compareFunctorByNodeID());
    if(it == info.end())
    {
        cerr<<"Some error in binary search for vertexID \n";
    }
    stream.write(&allReadConcatenated[it->position] ,k);
    for(long i = 1 ; i < path.size() ; i++)
    {
        kmer.nodeID = path[i];
        it = lower_bound(info.begin(), info.end(), kmer, compareFunctorByNodeID());
        stream.put( allReadConcatenated[(it->position) + k-1]);
    }
    stream << endl;



}
int main(int argc, char* argv[])
{

    if(argc !=5)
    {
        cerr<<"pab base L k delta"<<endl;
        exit(1);
    }
	string base = argv[1];
    string inputReadsFileName(base);
    string inputKmersFileName(base);
	inputReadsFileName += ".readl";
	inputKmersFileName += ".readks";
	L = atoi(argv[2]) ;
    k = atoi(argv[3]);
	int delta = atoi(argv[4]);
	int overlapThreshold = k - 2*delta;
	if (overlapThreshold < 10) {
		cout << "Warning: overlapThreshold is kinda low (" << overlapThreshold << ").\n";
		assert (overlapThreshold > 6); //crash in this case, its just too much!
	}
	string infoFileName(base);
	infoFileName+= ".info";
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
	cout<<"Reading kmer file...\n";
    FillInKmers(info, inputKmersFileName.c_str(), k);
    cout<<"Info size "<<info.size()<<endl;
    Partition(info,overlapThreshold);
    cout<<"End partioning, start sort \n";
    stable_sort(info.begin(), info.end(), compareFunctorByNodeID());
    cout<<"End sort start write\n";
    WriteInfo( info, infoFileName.c_str());
    long maxVertexID = info[info.size()-1].nodeID;
    vector<long> nextNodes(maxVertexID +1, -1) ;

    cout<<"End write start get graph str\n";
    GetSimpleGraphStructure(nextNodes, info);
    set<long> sourceNodes = GetSourcesNodes(nextNodes);
#ifdef _DEBUG_
    //PrintVector(nextNodes,cout);
    //PrintSet(sourceNodes);

#endif    
    ofstream outFile, detailOutFile, stringOutFile ;
    outFile.open((base + ".contigs").c_str());
    detailOutFile.open((base + ".contigsdetail").c_str());
    stringOutFile.open((base + ".contigsString").c_str());
	int counter=0;

    outFile<< "L="<< L <<  "\t" << "k=" << k << "\tDel=" << delta << endl;
    for(set<long>::iterator it = sourceNodes.begin(); it!= sourceNodes.end(); it++)
    {
        vector<long> path = FindPath(nextNodes, *it);
        PrintVector(path, detailOutFile);
		stringOutFile << ">" << counter++ << endl;
        FindAndPrintContig(path, stringOutFile, info);
        outFile << (path.size() + k)<<endl;
    }
    outFile.close();
    detailOutFile.close();
    stringOutFile.close();
    return 0;

}

