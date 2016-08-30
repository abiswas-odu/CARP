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
string getFirstMatch(vector<int> &orList,vector<int> &ctgSig,int errorThresh)
{
	for(int i=0;i<orList.size();i++)
	{
		int accErr=0,k=0,j;
		for(j=0;j<ctgSig.size() && accErr<errorThresh;j++)
		{
			if((i+k)<orList.size())
				accErr += abs(ctgSig[j]-orList[i+k]);
			else
				accErr += abs(ctgSig[j]-orList[(i+k)-orList.size()]);
			k++;
		}
		if(j==ctgSig.size() && accErr<errorThresh)
		{
			if((i+k)<orList.size())
				return SSTR(i) + "," + SSTR(i+k) + "," + SSTR(accErr);
			else
				return SSTR(i) + "," + SSTR(i+k-orList.size()) + "," + SSTR(accErr);
		}
	}
	return "";
}
string getBestMatch(vector<int> &orList,vector<int> &ctgSig)
{
	int beg,end,minErr=(10000000-1); 
	string ori="";
	for(int i=0;i<orList.size();i++)
	{
		int accErr=0,k=0,j;
		for(j=0;j<ctgSig.size();j++)
		{
			if((i+k)<orList.size())
				accErr += abs(ctgSig[j]-orList[i+k]);
			else
				accErr += abs(ctgSig[j]-orList[(i+k)-orList.size()]);
			k++;
		}
		
		if(accErr<minErr)
		{
			ori="FWD";
			minErr=accErr;
			beg=i;
			if((i+k)<orList.size())
				end = (i+k);
			else
				end = (i+k-orList.size());
		}
	}
	//Check Reverse
	for(int i=0;i<orList.size();i++)
	{
		int accErr=0,k=0,j;
		for(j=ctgSig.size()-1;j>=0;j--)
		{
			if((i+k)<orList.size())
				accErr += abs(ctgSig[j]-orList[i+k]);
			else
				accErr += abs(ctgSig[j]-orList[(i+k)-orList.size()]);
			k++;
		}
		
		if(accErr<minErr)
		{
			ori="REV";
			minErr=accErr;
			beg=i;
			if((i+k)<orList.size())
				end = (i+k);
			else
				end = (i+k-orList.size());
		}
	}
	return SSTR(beg) + "," + SSTR(end) + "," + SSTR(minErr)+","+ori;
}
int main(int argc, char **argv)
{
	//Program to find genes/intergenic regions 
	//arg 1: Matched Contigs
	//arg 2: OR Data
	//arg 3: Output matched contigs 	
	string blast_str; 
	map<string, string> ctg_list;
	vector<int> orList;
	vector<string> mathedCtgList;
	vector<int> ctgStart;
	vector<int> ctgEnd;
	vector<int> ctgErr;
	vector<string> ctgOri;
	
	ifstream ctgFile(argv[1]);
	ifstream orFile(argv[2]);
	if(ctgFile.fail() || orFile.fail())
	{
		cerr<<"\nThe files could not be opened.";
		return -1;
	}
	readContigFasta(ctgFile, ctg_list);
	while(getline(orFile,blast_str))
	{
		std::vector<std::string> tok = split(blast_str, ',');
		orList.push_back(atoi(tok[5].c_str()));
	}
	filebuf matchedCtg;
	matchedCtg.open(argv[4],ios::out);
	ostream os(&matchedCtg);
	
	for(map<string, string>::iterator it = ctg_list.begin(); it != ctg_list.end(); it++) 
	{
			vector<int> ctgPos;
			vector<int> ctgSig;
			int found = it->second.find("GCTAGC");
			while(found != -1)
			{
				ctgPos.push_back(found);
				found = it->second.find("GCTAGC",found+1);
			}
			for(int i=0;i<(ctgPos.size()-1) && ctgPos.size()>0;i++)
			{
				if((ctgPos[i+1]-ctgPos[i])>=50)
					ctgSig.push_back(ctgPos[i+1]-ctgPos[i]);
				else
					ctgPos.erase(ctgPos.begin()+i+1);
				//cout<<ctgPos[i+1]-ctgPos[i]<<"\n";
			}
			
			if(ctgSig.size()>=5)
			{
				os<<"Processing:"<<it->first<<","<<it->second.size()<<"\nCuts Found:"<<ctgSig.size()<<"\n";
				for(int i=0;i<ctgSig.size();i++)
				{
					os<<ctgSig[i]<<"\n";
				}
				string bestMatch = getBestMatch(orList,ctgSig);
				
				std::vector<std::string> tok = split(bestMatch, ',');
				mathedCtgList.push_back(it->first);
				ctgStart.push_back(atoi(tok[0].c_str()));
				ctgEnd.push_back(atoi(tok[1].c_str()));
				ctgErr.push_back(atoi(tok[2].c_str()));
				ctgOri.push_back(tok[3]);
				os<<"Best Match:"<<bestMatch<<"\n";
			}
		
	}
	os<<"Contig Cut Bitmap\n";
	string cutPintCons(orList.size(),'0');
	for(int i=0;i<ctgStart.size();i++)
	{
		string dispStr(orList.size(),'0');
		os<<mathedCtgList[i]<<":";
		if(ctgStart[i]<ctgEnd[i])
		{
			for(int j=ctgStart[i];j<ctgEnd[i];j++)
			{
				dispStr[j] = '1';
				cutPintCons[j] = '1';
				/*if(cutPintCons[j]>0 && ctgErr[cutPintCons[j]]<) 
				{
					cutPintCons[j] = j;
				}
				else
				{
					cutPintCons[j] = j;
				}*/
			}
				
		}
		else
		{
			for(int j=ctgStart[i];j<orList.size();j++)
			{
				dispStr[j] = '1';
				cutPintCons[j] = '1';
			}
			for(int j=0;j<ctgEnd[i];j++)
			{
				dispStr[j] = '1';
				cutPintCons[j] = '1';
			}
		}
		os<<dispStr<<"\n";
	}
	
	
	
	
	os<<"All Matched Cuts:\n"<<cutPintCons<<"\n";
	int totalMatchedCuts=0;
	for(int j=0;j<orList.size();j++)
		totalMatchedCuts += (int)cutPintCons[j]-'0';
	os<<"Total Cuts Matched="<<totalMatchedCuts<<"\n";
	matchedCtg.close();
	orFile.close();
	ctgFile.close();
	return 0;	
}	
