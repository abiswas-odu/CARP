#include<map>
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include <vector>

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
	//arg 1: Adj Gene Matched Results
	//arg 2: WGS Matched Results
	string blast_str; 
	vector<string> adjCtg;
	vector<string> WGSCtg;
	
	ifstream adjRes(argv[1]);
	ifstream wgsRes(argv[2]);
	
	if(adjRes.fail() || wgsRes.fail())
	{
		cerr<<"\nThe files could not be opened.";
		return -1;
	}
	
	while(getline(adjRes,blast_str))
	{
		adjCtg.push_back(blast_str);
		cout<<blast_str<<"\n";
	}
	
	while(getline(wgsRes,blast_str))
	{
		std::vector<std::string> tok1 = split(blast_str, ',');
		int flg=1;
		for(int i=0;i<adjCtg.size();i++)
		{
			std::vector<std::string> tok2 = split(adjCtg[i], ',');
			if(tok1[0]==tok2[0])
				flg=0;
		}
		if(flg)
			cout<<blast_str<<"\n";
	}
	adjRes.close();
	wgsRes.close();
	return 0;	
}	
