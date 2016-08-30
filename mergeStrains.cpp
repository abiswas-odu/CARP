#include<map>
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include <vector>
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
	//arg 1: Insertion Loci
	//arg 2: Insertion Loci
		
	string locus_str; 
	map<string,string> interruptedList;
	int i=0;
	
	ifstream loci1(argv[1]);
	ifstream loci2(argv[2]);
	if(loci1.fail() || loci2.fail())
	{
		cerr<<"\nThe files could not be opened.";
		return -1;
	}
	while(getline(loci1,locus_str))
	{
		std::vector<std::string> tok = split(locus_str, ',');
		interruptedList.insert(pair<string,string>(tok[0],tok[1]+","+tok[2]));
	}
	while(getline(loci2,locus_str))
	{
		std::vector<std::string> tok = split(locus_str, ',');
		if(interruptedList.count(tok[0]) > 0)
		{
			map<string,string>::iterator it;
			it = interruptedList.find(tok[0]);
			it->second = it->second + "," +tok[1]+","+tok[2];
		}
		else
		{
			interruptedList.insert(pair<string,string>(tok[0],",,"+tok[1]+","+tok[2]));
		}
	}
	loci1.close();
	loci2.close();
	map<string,string>::iterator it;
	for (it=interruptedList.begin() ; it != interruptedList.end(); it++ )
	{
		cout<<(*it).first<<","<<(*it).second<<"\n";
	}
	return 0;	
}	
