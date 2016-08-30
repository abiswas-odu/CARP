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
int main(int argc, char **argv)
{
	//Program to find genes/intergenic regions 
	//arg 1: Mate Read ID map (gkpStore.fastqUIDmap)
	//arg 2: ASM File
	//arg 3: Output matched contigs 	
	string blast_str; 
	map<string,string> ctgReadMap;
	
	ifstream mateReads(argv[1]);
	ifstream asmFile(argv[2]);
	
	if(mateReads.fail() || asmFile.fail())
	{
		cerr<<"\nThe files could not be opened.";
		return -1;
	}
	int index = 0;
	while(getline(mateReads,blast_str))
	{
		std::vector<std::string> tok = split(blast_str, '\t');
		string readName=tok[2];
		index = readName.find("/");
		readName.replace(index, 1, ".");
		ctgReadMap.insert(pair<string,string>(tok[0],readName));
	}
	
	filebuf modifiedAsm;
	modifiedAsm.open(argv[3],ios::out);
	ostream os(&modifiedAsm);
	
	while(getline(asmFile,blast_str))
	{
		index = blast_str.find("acc:");
		string accStr = blast_str;
		if(index != string::npos) 
		{
			string iid = blast_str.substr(5,12);
			map<string,string>::iterator it;
			it = ctgReadMap.find(iid);
			if(it != ctgReadMap.end())
			{
				accStr = "acc:(" + it->second + blast_str.substr(blast_str.find(","));
			}
		}
		index = blast_str.find("mid:");
		if(index != string::npos) 
		{
			string iid = blast_str.substr(4);
			map<string,string>::iterator it;
			it = ctgReadMap.find(iid);
			if(it != ctgReadMap.end())
			{
				accStr = "mid:" + it->second;
			}
		}
		os<<accStr<<"\n";
	}
	modifiedAsm.close();
	asmFile.close();
	mateReads.close();
	return 0;	
}	
