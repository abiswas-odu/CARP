#include<map>
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include <vector>
#include <algorithm>

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
				string temp = locus+","+prod;
				theMap.insert (pair<string,string>(id,temp));
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
	theMap.insert (pair<string,string>(id,temp) );	
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
	//arg 1: BLAST result file	
	//arg 2: Genbank File	
		
	string blast_str; 
	int i=0;
	
	ifstream isBlastRes(argv[1]);
	ifstream gbFile(argv[2]);
		
	if(isBlastRes.fail() || gbFile.fail())
	{
		cerr<<"\nThe files could not be opened.";
		return -1;
	}
	map<string,string> genBank_rec;
	readGenBankCDS(gbFile,genBank_rec);
	while(getline(isBlastRes,blast_str))
	{
		std::vector<std::string> tok = split(blast_str, ',');
		int isCoLocated=0,isMiddle=0; 
		
		int coLocX = atoi(tok[2].c_str());
		int coLocY = atoi(tok[3].c_str());
		
		if(coLocX<=5) {
			std::vector<std::string> ISBoundTok = split(tok[0], '|');
			int searchBoundCord = atoi(ISBoundTok[5].c_str());
			int ISBoundX= atoi(ISBoundTok[2].c_str());
			int ISBoundY= atoi(ISBoundTok[3].c_str());
			if(searchBoundCord==ISBoundX || searchBoundCord==ISBoundY)
				isCoLocated=1;
		}
		if(coLocY>=195) {
			std::vector<std::string> ISBoundTok = split(tok[0], '|');
			int searchBoundCord = atoi(ISBoundTok[6].c_str());
			int ISBoundX= atoi(ISBoundTok[2].c_str());
			int ISBoundY= atoi(ISBoundTok[3].c_str());
			if(searchBoundCord==ISBoundX || searchBoundCord==ISBoundY)
				isCoLocated=1;
		}
		
		std::vector<std::string> geneBoundTok = split(tok[1], ':');
		string geneBoundStr = geneBoundTok[1];
		if(geneBoundStr[0]=='c')
			geneBoundStr.erase(0,1);
		geneBoundTok = split(geneBoundStr, '-');
		int geneX = atoi(geneBoundTok[0].c_str());
		int geneY = atoi(geneBoundTok[1].c_str());
		
		int geneLen = abs(geneX-geneY);
		
		int geneBlastHitX = atoi(tok[4].c_str());
		int geneBlastHitY = atoi(tok[5].c_str());
		
		if((geneBlastHitX<(geneLen-3) && geneBlastHitX>3) || (geneBlastHitY<(geneLen-3) && geneBlastHitY>3))
			isMiddle=1;
		
		string cdsKey="";
		if(getGeneStart(tok[1]) > getGeneEnd(tok[1]))
			cdsKey = SSTR(getGeneEnd(tok[1])) + ".." + SSTR(getGeneStart(tok[1]));
		else
			cdsKey = SSTR(getGeneStart(tok[1])) + ".." + SSTR(getGeneEnd(tok[1]));
		
		map<string,string>::iterator it;
		it = genBank_rec.find(cdsKey);
		string infoTag = (*it).second;
		
		if(isMiddle==1 && isCoLocated==1) {
			cout<<blast_str<<",Interrupting,"<<infoTag<<"\n";
			
		}
		else {
			cout<<blast_str<<",Non-Interrupting,"<<infoTag<<"\n";
		}
	}
	gbFile.close();
	isBlastRes.close();
	return 0;	
}	
