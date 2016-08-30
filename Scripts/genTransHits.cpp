#include<map>
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include <vector>

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

void readFasta(ifstream& multiFile, map<string,string> &theMap)
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
	theMap.insert (pair<string,string>(id,seq) );
}
int main(int argc, char **argv)
{
	//Program to find genes/intergenic regions 
	//arg 1: Gene Matched Results
	string blast_str; 
	map<string,string> ctg_cuts;
	map<string,string> ctg_list;
	map<string,string> transHit_list;
		
	ifstream blastRes(argv[1]);
	ifstream ctgRes(argv[2]);
		
	if(blastRes.fail() || ctgRes.fail())
	{
		cerr<<"\nThe files could not be opened.";
		return -1;
	}
	readFasta(ctgRes,ctg_list);
	while(getline(blastRes,blast_str))
	{
		std::vector<std::string> tok = split(blast_str, ',');
		string ctgID =  tok[0];
		int endIndx = tok.size();
		int qstart = atoi(tok[endIndx-4].c_str());
		int qend = atoi(tok[endIndx-3].c_str());
		
		std::map<string,string>::iterator iter=ctg_cuts.find(ctgID);
		if(iter != ctg_cuts.end())
		{
			string cutRange = (*iter).second;
			std::vector<std::string> tok2 = split(cutRange, ',');
			int prev_qstart = atoi(tok2[0].c_str());
			int prev_qend = atoi(tok2[1].c_str());
			int cur_qsrart, cur_qend;
			if(prev_qstart < qstart)
				cur_qsrart = prev_qstart;
			else
				cur_qsrart = qstart;
				
			if(prev_qend < qend)
				cur_qend = qend;
			else
				cur_qend = prev_qend;
				
			cutRange =  SSTR(cur_qsrart) + "," + SSTR(cur_qend);
			iter->second = cutRange;
		}
		else
		{
			string cutRange =  SSTR(qstart) + "," + SSTR(qend);
			ctg_cuts.insert(pair<string,string>(ctgID,cutRange));
		}
	}
	int i=0;
	for(std::map<string,string>::iterator it=ctg_cuts.begin(); it!=ctg_cuts.end(); ++it)
	{
		string key = (*it).first;
		string cutRange = (*it).second;
		std::vector<std::string> tok = split(cutRange, ',');
		int qstart = atoi(tok[0].c_str())-1;
		int qend = atoi(tok[1].c_str())-1;
		
		std::map<string,string>::iterator iter=ctg_list.find(key);
		string seq = (*iter).second;
		string transHit = seq.substr(qstart,qend-qstart);
		
		if(transHit.size() > 100)
		{
			std::map<string,string>::iterator iter2=transHit_list.find(transHit);
			if(iter2 == transHit_list.end())
			{
				//cout<<">seq_"<<i++<<key<<"\n";
				cout<<">seq_"<<i++<<"\n";
				cout<<transHit<<"\n";
				transHit_list.insert(pair<string,string>(transHit,key));
			}
		}
	}
	blastRes.close();
	ctgRes.close();
	return 0;	
}	
