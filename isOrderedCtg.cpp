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

class ISHits
{
	private: 
		string id;
		string seq;
		vector<string> contigIds;
		vector< vector<int> > QPos;
		vector< vector<int> > HPos;
		vector<int> begEndFlag;
		vector<int> fwdRevFlag;
	public:
		void setId(string s) { id=s;}
		string getId() {return id;}
		void setSeq(string s) {seq=s;}
		string getSeq() {return seq;}
		int getSeqLen() {return seq.size();}
		void addHit(string contigId,int qs,int qe, int hs, int he, int startEnd)
		{
			contigIds.push_back(contigId);
			vector<int> row;
			row.push_back(qs);
			row.push_back(qe);
			QPos.push_back(row);
			
			vector<int> row2;
			row2.push_back(hs);
			row2.push_back(he);
			HPos.push_back(row2);
			
			if(hs<=he)
				fwdRevFlag.push_back(1);
			else
				fwdRevFlag.push_back(0);
			
			begEndFlag.push_back(startEnd);
				
		}
		string getAdjSeq(string ctgSeq, int hs, int he)
		{
			if(hs<=10 || he<=10)
			{
				if(hs<=he)
					return ctgSeq.substr(he+1, 200);
				else
					return ctgSeq.substr(hs+1, 200);
			}
			else
			{
				if(hs<=he)
					return ctgSeq.substr(hs-201, 200);
				else
					return ctgSeq.substr(he-201, 200);
			}			
		}
		void printMatchedList(int i,int j, map<string,string> &theContigMap)
		{
			//cout<<contigIds[i]<<"|"<<HPos[i][0]<<"|"<<HPos[i][1]<<"|"<<contigIds[j]<<"|"<<HPos[j][0]<<"|"<<HPos[j][1]<<"|"<<id<<"|"<<QPos[j][0]<<"|"<<QPos[j][1]<<"|"<<QPos[i][0]<<"|"<<QPos[i][1]<<"\n";
			map<string,string>::iterator itCtg;
			itCtg = theContigMap.find(contigIds[i]);
			string contig1Substr = getAdjSeq(itCtg->second,HPos[i][0],HPos[i][1]); 
			cout<<">"<<contigIds[i]<<"|"<<HPos[i][0]<<"|"<<HPos[i][1]<<"|"<<contigIds[j]<<"|"<<HPos[j][0]<<"|"<<HPos[j][1]<<"|"<<id<<"|"<<QPos[j][0]<<"|"<<QPos[j][1]<<"|"<<QPos[i][0]<<"|"<<QPos[i][1]<<"|a\n"<<contig1Substr<<"\n";
			
			itCtg = theContigMap.find(contigIds[j]);
			string contig2Substr = getAdjSeq(itCtg->second,HPos[j][0],HPos[j][1]); 
			cout<<">"<<contigIds[i]<<"|"<<HPos[i][0]<<"|"<<HPos[i][1]<<"|"<<contigIds[j]<<"|"<<HPos[j][0]<<"|"<<HPos[j][1]<<"|"<<id<<"|"<<QPos[j][0]<<"|"<<QPos[j][1]<<"|"<<QPos[i][0]<<"|"<<QPos[i][1]<<"|b\n"<<contig2Substr<<"\n";
		}
		void joinedContigs(map<string,string> &theContigMap)
		{
			for(int i=0;i<QPos.size();i++)
			{
				for(int j=i+1;j<QPos.size();j++)
				{
					if(contigIds[i]!=contigIds[j])
					{
						if(QPos[i][1]<=QPos[j][0])
						{
							//BR/BF, BR/ER, EF/BF, EF/ER
							if(begEndFlag[i]==1 && fwdRevFlag[i]==0)
							{
								if(begEndFlag[j]==1 && fwdRevFlag[j]==1)
								{
									printMatchedList(i,j,theContigMap);
								}
								else if(begEndFlag[j]==0 && fwdRevFlag[j]==0)
								{
									printMatchedList(i,j,theContigMap);
								}
							}
							else if(begEndFlag[i]==0 && fwdRevFlag[i]==1)
							{
								if(begEndFlag[j]==1 && fwdRevFlag[j]==1)
								{
									printMatchedList(i,j,theContigMap);
								}
								else if(begEndFlag[j]==0 && fwdRevFlag[j]==0)
								{
									printMatchedList(i,j,theContigMap);
								}
							}
							//BF/BR, BF/EF, ER/BR, ER/EF
							if(begEndFlag[i]==1 && fwdRevFlag[i]==1)
							{
								if(begEndFlag[j]==1 && fwdRevFlag[j]==0)
								{
									printMatchedList(i,j,theContigMap);
								}
								else if(begEndFlag[j]==0 && fwdRevFlag[j]==1)
								{
									printMatchedList(i,j,theContigMap);
								}
							}
							else if(begEndFlag[i]==0 && fwdRevFlag[i]==0)
							{
								if(begEndFlag[j]==1 && fwdRevFlag[j]==0)
								{
									printMatchedList(i,j,theContigMap);
								}
								else if(begEndFlag[j]==0 && fwdRevFlag[j]==1)
								{
									printMatchedList(i,j,theContigMap);
								}
							}
						}
						else if(QPos[j][1]<=QPos[i][0])
						{
							//BR/BF, BR/ER, EF/BF, EF/ER
							if(begEndFlag[j]==1 && fwdRevFlag[j]==0)
							{
								if(begEndFlag[i]==1 && fwdRevFlag[i]==1)
								{
									printMatchedList(i,j,theContigMap);
								}
								else if(begEndFlag[j]==0 && fwdRevFlag[j]==0)
								{
									printMatchedList(i,j,theContigMap);
								}
							}
							else if(begEndFlag[j]==0 && fwdRevFlag[j]==1)
							{
								if(begEndFlag[i]==1 && fwdRevFlag[i]==1)
								{
									printMatchedList(i,j,theContigMap);
								}
								else if(begEndFlag[j]==0 && fwdRevFlag[j]==0)
								{
									printMatchedList(i,j,theContigMap);
								}
							}
							//BF/BR, BF/EF, ER/BR, ER/EF
							if(begEndFlag[i]==1 && fwdRevFlag[i]==1)
							{
								if(begEndFlag[j]==1 && fwdRevFlag[j]==0)
								{
									printMatchedList(i,j,theContigMap);
								}
								else if(begEndFlag[j]==0 && fwdRevFlag[j]==1)
								{
									printMatchedList(i,j,theContigMap);
								}
							}
							else if(begEndFlag[i]==0 && fwdRevFlag[i]==0)
							{
								if(begEndFlag[j]==1 && fwdRevFlag[j]==0)
								{
									printMatchedList(i,j,theContigMap);
								}
								else if(begEndFlag[j]==0 && fwdRevFlag[j]==1)
								{
									printMatchedList(i,j,theContigMap);
								}
							}
						}
					}
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
void readISFasta(ifstream& multiFile, map<string,ISHits> &theMap)
{
	string str,id,seq; 
	getline(multiFile,str);
	id=str.substr(1); 
	while(getline(multiFile,str))
	{
		int found = str.find(">");
		if (found!=string::npos)
		{
			ISHits gen;
			gen.setId(id);
			gen.setSeq(seq);
			theMap.insert (pair<string,ISHits>(id,gen) );
			seq="";
			id = str.substr(1);
		}
		else
		{
			seq=seq+str;
		}
	}
	ISHits gen;
	gen.setId(id);
	gen.setSeq(seq);
	theMap.insert (pair<string,ISHits>(id,gen) );
}

bool isAllowedContig(vector<string> &restricted_contigs, string ctg)
{
	for(int i=0;i<restricted_contigs.size();i++)
	{
		if(restricted_contigs[i]==ctg)
			return false;
	}
	return true; 
}

int main(int argc, char **argv)
{
	//Program to find genes/intergenic regions 
	//arg 1: Contig Sequences
	//arg 2: IS Sequences
	//arg 3: BLAST result
	//arg 4: Restricted Contigs
	
	string blast_str; 
	map<string,ISHits> is_list;
	map<string,string> contig_list;
	vector<string> restricted_contigs;
	int startEnd=0;
	
	ifstream contigs(argv[1]);
	ifstream insSeq(argv[2]);
	ifstream isBlastRes(argv[3]);
	ifstream resContigs(argv[4]);
		
	if(contigs.fail() || insSeq.fail() || isBlastRes.fail() || resContigs.fail())
	{
		cerr<<"\nThe files could not be opened.";
		return -1;
	}
	
	readContigFasta(contigs,contig_list);
	readISFasta(insSeq,is_list);
	
	while(getline(resContigs,blast_str))
	{
		restricted_contigs.push_back(blast_str);
	}
	
	while(getline(isBlastRes,blast_str))
	{
		std::vector<std::string> tok = split(blast_str, ',');
		map<string,string>::iterator itCtg;
		itCtg = contig_list.find(tok[1]);
		int len = itCtg->second.size();
		if((atoi(tok[4].c_str())<=250 || atoi(tok[5].c_str())<=250 || atoi(tok[4].c_str())>=(len-250) || atoi(tok[5].c_str())>=(len-250)) && isAllowedContig(restricted_contigs,tok[1]))
		{
			map<string,ISHits>::iterator it;
			it = is_list.find(tok[0]);
			startEnd=0;
			if(atoi(tok[4].c_str())<=10 || atoi(tok[5].c_str())<=10)
			{	startEnd=1;  }
			if(atoi(tok[2].c_str())<=50 || atoi(tok[3].c_str()) >= ((*it).second.getSeqLen()-50))
				(*it).second.addHit(tok[1],atoi(tok[2].c_str()),atoi(tok[3].c_str()),atoi(tok[4].c_str()),atoi(tok[5].c_str()),startEnd); 
		}
	}
	map<string,ISHits>::iterator it;
	for(it = is_list.begin(); it != is_list.end(); it++) {
		// iterator->first = key
		// iterator->second = value
		// Repeat if you also want to iterate through the second map.
		it->second.joinedContigs(contig_list);
	}
	contigs.close();
	insSeq.close();
	isBlastRes.close();
	return 0;	
}	
