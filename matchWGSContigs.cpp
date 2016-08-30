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

int getGeneStart(string gene) 
{
	string retStr = gene.substr(gene.find(":")+1,gene.find("-")-gene.find(":")+1);
	if(retStr[0]=='c')
		return atoi(retStr.substr(1).c_str());
	else
		return atoi(retStr.c_str());
}
		
int getGeneEnd(string gene)
{
	string retStr = gene.substr(gene.find("-")+1);
	return atoi(retStr.c_str());
}

bool isNotRepeatedCtg(vector<string> &joinedCtg, vector<string> &locA, vector<string> &locB, string ctgStr, string a, string b)
{
	for(int i=0;i<joinedCtg.size();i++)
	{
		if((joinedCtg[i]==ctgStr) && ((locA[i]==a && locB[i]==b) || (locA[i]==b && locB[i]==a)))
			return false;
	}
	return true;
}

int main(int argc, char **argv)
{
	//Program to find genes/intergenic regions 
	//arg 1: Gene Matched Results
	string blast_str; 
	multimap<string,string> blast_list;
	vector<string> joinedCtg;
	vector<string> locA;
	vector<string> locB;
		
	ifstream geneBlastRes(argv[1]);
		
	if(geneBlastRes.fail())
	{
		cerr<<"\nThe files could not be opened.";
		return -1;
	}
	
	while(getline(geneBlastRes,blast_str))
	{
		std::vector<std::string> tok = split(blast_str, ',');
		string locusTag = tok[4] + "," + tok[5];
		blast_list.insert(pair<string,string>(tok[0],locusTag));
	}
	for(std::multimap<string,string>::iterator it=blast_list.begin(); it!=blast_list.end(); ++it)
	{
		string key = (*it).first;
		string value = (*it).second;
		std::vector<std::string> tok = split(value, ',');
		int x1 = atoi(tok[0].c_str());
		int x2 = atoi(tok[1].c_str());
		char lastChar = key.at( key.length() - 1 );
		string newKey = "";
		if(lastChar=='a')
			newKey = key.substr(0,key.size()-1) + "b";
		else
			newKey = key.substr(0,key.size()-1) + "a";
		
		pair<multimap<string, string>::iterator, multimap<string, string>::iterator> ppp;
		ppp = blast_list.equal_range(newKey);
		for (multimap<string,string>::iterator it2 = ppp.first; it2 != ppp.second; ++it2)
		{
			string value2 = (*it2).second;
			std::vector<std::string> tok2 = split(value2, ',');
			int x3 = atoi(tok2[0].c_str());
			int x4 = atoi(tok2[1].c_str());
			int v1 = (x1+x2)/2;
			int v2 = (x3+x4)/2;
			
			if(abs(v1-v2)<=5000)
			{
				string printStr = newKey.substr(0,newKey.size()-2);
				if(isNotRepeatedCtg(joinedCtg,locA,locB,printStr,value,value2))
				{
					joinedCtg.push_back(printStr);
					locA.push_back(value);
					locB.push_back(value2);
					cout<<printStr<<"|("<<x1<<";"<<x2<<")|("<<x3<<";"<<x4<<")"<<"\n";
				}
			}
		}
	}
	geneBlastRes.close();
	return 0;	
}	
