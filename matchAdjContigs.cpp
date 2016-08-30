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

void readGenBankCDS(ifstream& cdsFile, map<string,string> &theMap)
{
	string id,locus="",prod="",str; 
	int wFlg=0;
	while(getline(cdsFile,str))
	{
		if(str.substr(5,3)=="CDS")
		{
			if(wFlg==1)
			{
				//write
				string temp = locus+","+prod;
				//cout<<id<<","<<temp<<'\n';
				theMap.insert (pair<string,string>(id,locus) );
			}
			wFlg=1;
			if(str[21]=='c')
			{
				id=str.substr(32);
				id=id.substr(0,id.length()-1);
			}	
			else
				id =  str.substr(21);
		}
		
		int found = str.find("locus_tag=");
		if (found != string::npos)
		{
			locus=str.substr(33);
			locus=locus.substr(0,locus.length()-1);
		}
		
		found = str.find("product=");
		if(found!=string::npos)
		{
			prod=str.substr(31);
			prod=prod.substr(0,prod.length()-1);
			std::replace(prod.begin(), prod.end(), ',', ';');
		}
	}
	string temp = locus+","+prod;
	theMap.insert (pair<string,string>(id,locus) );	
}

bool isNotRepeatedCtg(vector<string> &joinedCtg, vector<int> &locA, vector<int> &locB, string ctgStr, int a, int b)
{
	for(int i=0;i<joinedCtg.size();i++)
	{
		if((joinedCtg[i]==ctgStr) && ((locA[i]==a && locB[i]==b) || (locA[i]==b && locB[i]==a)))
			return false;
	}
	return true;
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

int main(int argc, char **argv)
{
	//Program to find genes/intergenic regions 
	//arg 1: Gene Matched Results
	//arg 2: Gene Bank
	string blast_str; 
	multimap<string,int> blast_list;
	map<string,string> genBank_rec;
	vector<string> joinedCtg;
	vector<int> locA;
	vector<int> locB;
	int i=0;
	
	ifstream geneBlastRes(argv[1]);
	ifstream genBank(argv[2]);
	
	if(geneBlastRes.fail() || genBank.fail())
	{
		cerr<<"\nThe files could not be opened.";
		return -1;
	}
	
	readGenBankCDS(genBank,genBank_rec);
	
	while(getline(geneBlastRes,blast_str))
	{
		std::vector<std::string> tok = split(blast_str, ',');
		
		string cdsKey="";
		if(getGeneStart(tok[1]) > getGeneEnd(tok[1]))
			cdsKey = SSTR(getGeneEnd(tok[1])) + ".." + SSTR(getGeneStart(tok[1]));
		else
			cdsKey = SSTR(getGeneStart(tok[1])) + ".." + SSTR(getGeneEnd(tok[1]));
		
		map<string,string>::iterator it;
		it = genBank_rec.find(cdsKey);
		string locusTag = (*it).second;
		int locusNum = atoi(locusTag.substr(locusTag.find("_")+1).c_str());
		blast_list.insert(pair<string,int>(tok[0],locusNum));
	}
	for(std::multimap<string,int>::iterator it=blast_list.begin(); it!=blast_list.end(); ++it)
	{
		string key = (*it).first;
		int value = (*it).second;
		char lastChar = key.at( key.length() - 1 );
		string newKey = "";
		if(lastChar=='a')
			newKey = key.substr(0,key.size()-1) + "b";
		else
			newKey = key.substr(0,key.size()-1) + "a";
		
		pair<multimap<string, int>::iterator, multimap<string, int>::iterator> ppp;
		ppp = blast_list.equal_range(newKey);
		for (multimap<string,int>::iterator it2 = ppp.first; it2 != ppp.second; ++it2)
		{
			int value2 = (*it2).second;
			if(abs(value-value2)<=2)
			{
				string printStr = newKey.substr(0,newKey.size()-2);
				if(isNotRepeatedCtg(joinedCtg,locA,locB,printStr,value,value2))
				{
					joinedCtg.push_back(printStr);
					locA.push_back(value);
					locB.push_back(value2);
					cout<<printStr<<"|MMAR_"<<value<<"|MMAR_"<<value2<<"\n";
				}
			}
		}
	}
	geneBlastRes.close();
	genBank.close();
	return 0;	
}	
