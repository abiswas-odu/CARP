#include<map>
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include <vector>
#include <omp.h>

#define GAPLEN 200
#define GAPOVERLAP 20
#define SSTR( x ) dynamic_cast< std::ostringstream & >(( std::ostringstream() << std::dec << x ) ).str()

using std::istringstream;
using std::stringstream;
using std::ifstream;
using std::ofstream;
using std::string;
using std::cout;
using std::cin;
using std::cerr;
using std::map;
using std::multimap;

using namespace std;

vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
vector<string> split(const string &s, char delim) {
    vector<string> elems;
    return split(s, delim, elems);
}
void readContigFasta(ifstream& multiFile, map<string,string> &theMap)
{
	string str,id,seq; 
	getline(multiFile,str);
	id=str.substr(1); 
	while(getline(multiFile,str))
	{
		int found = str.find(">");
		if (found!=string::npos)
		{
			theMap.insert (pair<string,string>(id,seq) );
			seq="";
			id = str.substr(1);
		}
		else
		{
			seq=seq+str;
		}
	}
	theMap.insert(pair<string,string>(id,seq));
}

int main(int argc, char **argv)
{
	//Program to find genes/intergenic regions 
	//arg 1: IS Joined Contigs
	//arg 2: BLAST search vs Scaffold Sequences
	
	string blast_str,ctgID=""; 
	int maxlen=0;
	map<string, string> ctg_list;
	
	ifstream isJoinedCtgs(argv[1]);
	ifstream blastSearchSeqs(argv[2]);
			
	if(isJoinedCtgs.fail() || blastSearchSeqs.fail())
	{
		cerr<<"\nThe files could not be opened.";
		return -1;
	}
	readContigFasta(blastSearchSeqs, ctg_list);
	
	while(getline(isJoinedCtgs,blast_str))
	{
		std::vector<std::string> tok = split(blast_str, ',');
		if(ctgID==tok[0])
		{
			int mlen = atoi(tok[6].c_str());
			if(maxlen<mlen)
				maxlen=mlen;
		}
		else
		{
			if(ctgID!="")
			{
				map<string,string>::iterator itCtg;
				itCtg = ctg_list.find(ctgID);
				int diffLen = itCtg->second.size()-maxlen;
				if(abs(diffLen)<=100)
					cout<<ctgID<<"|T\n";
				else
					cout<<ctgID<<"|F\n";
			}
			ctgID=tok[0];
			maxlen = atoi(tok[6].c_str());
		}
			
	}
	
	blastSearchSeqs.close();
	isJoinedCtgs.close();
	return 0;	
}	
