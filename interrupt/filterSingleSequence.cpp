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

class Contigs
{
	private: 
		string id;
		string seq;
		vector< vector<int> > matchPos;
	public:
		void setId(string s) { id=s;}
		string getId() {return id;}
		void setSeq(string s) {seq=s;}
		string getSeq() {return seq;}
		int getSeqLen()	{ return seq.size(); }
		void addMatch(int st,int en)
		{
			vector<int> row;
			row.push_back(st-1);
			row.push_back(en-1);
			matchPos.push_back(row);
		}
		int calculateGapLength()
		{
			int countCoverage = 0;
			if(matchPos.size()>0)
			{
				vector<int> gr;
				for(int i=0;i<seq.size();i++)
					gr.push_back(1);
				for(int i=0;i<matchPos.size();i++)
					for(int j=matchPos[i][0];j<=matchPos[i][1];j++)
						gr[j]=0;
				
				for(int i=0;i<seq.size();i++)
				{
					if(gr[i]==0)
						countCoverage++;
				}
			}
			return countCoverage;
		}
};

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

void readFasta(ifstream& multiFile, map<string,Contigs> &theMap)
{
	string str,id,seq; 
	getline(multiFile,str);
	id=str.substr(1); 
	while(getline(multiFile,str))
	{
		int found = str.find(">");
		if (found!=string::npos)
		{
			Contigs gen;
			gen.setId(id);
			gen.setSeq(seq);
			theMap.insert (pair<string,Contigs>(id,gen) );
			seq="";
			id = str.substr(1);
		}
		else
		{
			seq=seq+str;
		}
	}
	Contigs gen;
	gen.setId(id);
	gen.setSeq(seq);
	theMap.insert (pair<string,Contigs>(id,gen) );
}
int main(int argc, char **argv)
{
	//Program to find genes/intergenic regions 
	//arg 1: BLAST result file	
		
	string blast_str; 
	int i=0;
	
	ifstream isBlastRes(argv[1]);
		
	if(isBlastRes.fail())
	{
		cerr<<"\nThe files could not be opened.";
		return -1;
	}
	while(getline(isBlastRes,blast_str))
	{
		std::vector<std::string> tok = split(blast_str, ',');
		
		
		map<string,Contigs>::iterator it;
		it = ctg_list.begin();
		int x = atoi(tok[4].c_str());
		int y = atoi(tok[5].c_str());
		
		if(x<=y) {
			cout<<">"<<tok[0]<<"|"<<tok[1]<<"|"<<x<<"|"<<y<<"|"<<"F"<<"|"<<x-200<<"|"<<x<<"|Begin\n"<<(*it).second.getSeq().substr(x-200,200)<<"\n";
			cout<<">"<<tok[0]<<"|"<<tok[1]<<"|"<<x<<"|"<<y<<"|"<<"F"<<"|"<<y<<"|"<<y+200<<"|End\n"<<(*it).second.getSeq().substr(y,200)<<"\n";
		}
		else {
			cout<<">"<<tok[0]<<"|"<<tok[1]<<"|"<<x<<"|"<<y<<"|"<<"C"<<"|"<<y-200<<"|"<<y<<"|End\n"<<(*it).second.getSeq().substr(y-200,200)<<"\n";
			cout<<">"<<tok[0]<<"|"<<tok[1]<<"|"<<x<<"|"<<y<<"|"<<"C"<<"|"<<x<<"|"<<x+200<<"|Begin\n"<<(*it).second.getSeq().substr(x,200)<<"\n";
		}
	}
	isBlastRes.close();
	return 0;	
}	
