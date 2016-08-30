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

void readFasta(ifstream& multiFile, vector<string> &theVec)
{
	string str,id,seq; 
	getline(multiFile,str);
	id=str.substr(1); 
	while(getline(multiFile,str))
	{
		int found = str.find(">");
		if (found!=string::npos)
		{
		    theVec.push_back(seq);
			id = str.substr(1);
		}
		else
		{
			seq=seq+str;
		}
	}
	theVec.push_back(seq);
}

int main(int argc, char **argv)
{
	//Program to find genes/intergenic regions 
	//arg 1: Milti Fasta File
	string len_str; 
	vector<string> seq_list;
		
	ifstream fstFile(argv[1]);
		
	if(fstFile.fail())
	{
		cerr<<"\nThe files could not be opened.";
		return -1;
	}
	
	readFasta(fstFile,seq_list);
	for(int i=0;i<seq_list.size();i++)
		len_str = len_str + "," + SSTR(seq_list[i].size());
	
	cout<<len_str<<"\n";
	return 0;
}
