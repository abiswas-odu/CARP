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

class UnusedReads
{
	private: 
		string id;
		string seq;
		vector< vector<int> > matchPos;
		vector< vector<int> > gapPos;
	public:
		void setId(string s) { id=s;}
		string getId() {return id;}
		void setSeq(string s) {seq=s;}
		string getSeq() {return seq;}
		int getGapCtr() {return gapPos.size();}
		void addMatch(int st,int en)
		{
			vector<int> row;
			row.push_back(st-1);
			row.push_back(en-1);
			matchPos.push_back(row);
		}
		void calculateGaps()
		{
			if(matchPos.size()>0)
			{
				vector<int> gr;
				for(int i=0;i<seq.size();i++)
					gr.push_back(1);
				for(int i=0;i<matchPos.size();i++)
					for(int j=matchPos[i][0];j<=matchPos[i][1];j++)
						gr[j]=0;
				
				int st=0,ed=-1,inRun=0;
				for(int i=0;i<seq.size();i++)
				{
					if(gr[i]==1)
					{
						if(inRun==0)
						{
							st=ed=i;
							inRun=1;
						}
						else
							ed++;
					}
					else
					{
						if(inRun==1)
						{
							inRun=0;
							vector<int> grow;
							grow.push_back(st);
							grow.push_back(ed);
							gapPos.push_back(grow);
						}
					}
				}
				if(inRun==1)
				{
					vector<int> grow;
					grow.push_back(st);
					grow.push_back(ed);
					gapPos.push_back(grow);
				}
			}
		}
		void printMatchPos()
		{
			for(int i=0;i<matchPos.size();i++)
			{
				cout<<"("<<matchPos[i][0]+1<<","<<matchPos[i][1]+1<<"),";
			}
		}
		void printGaps()
		{
			for(int i=0;i<gapPos.size();i++)
			{
				cout<<"[("<<gapPos[i][0]+1<<","<<gapPos[i][1]+1<<");"<<(gapPos[i][1]-gapPos[i][0])<<"],";
			}
		}
		void printGapSeqs()
		{
			for(int i=0;i<gapPos.size();i++)
			{
				if((gapPos[i][1]-gapPos[i][0])>20 && ((gapPos[i][0]==0)||(gapPos[i][1]==(seq.length()-1))))
				{
					cout<<">"<<id<<":("<<gapPos[i][0]+1<<","<<gapPos[i][1]+1<<")\n"<<seq.substr(gapPos[i][0],(gapPos[i][1]-gapPos[i][0]))<<"\n";
				}
			}
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
void readFasta(ifstream& multiFile, map<string,UnusedReads> &theMap)
{
	string str,id,seq; 
	getline(multiFile,str);
	id=str.substr(1); 
	while(getline(multiFile,str))
	{
		int found = str.find(">");
		if (found!=string::npos)
		{
			UnusedReads gen;
			gen.setId(id);
			gen.setSeq(seq);
			theMap.insert (pair<string,UnusedReads>(id,gen) );
			seq="";
			id = str.substr(1);
		}
		else
		{
			seq=seq+str;
		}
	}
	UnusedReads gen;
	gen.setId(id);
	gen.setSeq(seq);
	theMap.insert (pair<string,UnusedReads>(id,gen) );
}

int main(int argc, char **argv)
{
	//Program to find genes/intergenic regions 
	//arg 1: Un-mapped Reads
	//arg 2: BLAST result
	
	string blast_str; 
	map<string,UnusedReads> read_list;
	int i=0;
	
	ifstream reads(argv[1]);
	ifstream blastRes(argv[2]);
	if(reads.fail() || blastRes.fail())
	{
		cerr<<"\nThe files could not be opened.";
		return -1;
	}
	
	readFasta(reads,read_list);
	
	//test read 
	/*int ctr=0;
	map<string,UnusedReads>::iterator it;
	for (it=read_list.begin() ; it != read_list.end(); it++ )
	{
		cout<<">"<<(*it).first<<'\n'<<(*it).second.getSeq()<<'\n';
		ctr++;
	}
	cout<<"\n Read Count"<<ctr;*/
	
	while(getline(blastRes,blast_str))
	{
		std::vector<std::string> tok = split(blast_str, ',');
		if(read_list.count(tok[1]) > 0)
		{
			map<string,UnusedReads>::iterator it;
			it = read_list.find(tok[1]);
			int s1 = atoi(tok[4].c_str());
			int s2 = atoi(tok[5].c_str());
			if(s1>s2)
				it->second.addMatch(s2,s1);
			else
				it->second.addMatch(s1,s2);
		}
	}
	map<string,UnusedReads>::iterator it;
	for (it=read_list.begin() ; it != read_list.end(); it++ )
		it->second.calculateGaps();
	
	//Print gapped sequences
	for (it=read_list.begin() ; it != read_list.end(); it++ )
	{
		if((*it).second.getGapCtr()>0)
		{
			(*it).second.printGapSeqs();
		}
	}
	
	//Print unmatched sequences
	/*for (it=read_list.begin() ; it != read_list.end(); it++ )
	{
		if((*it).second.getGapCtr()==0)
		{
			cout<<">"<<(*it).second.getId()<<"\n"<<(*it).second.getSeq()<<"\n";
		}
	}*/
	
}	
