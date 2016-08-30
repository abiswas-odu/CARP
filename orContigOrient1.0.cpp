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

class ContigCuts
{
	string seqId;
	string seq;
	vector<int> ctgCuts;
	
	public:
	static bool compareLenPredicate(ContigCuts lhs, ContigCuts rhs) 
	{ 
		return (lhs.ctgCuts.size() < rhs.ctgCuts.size()); 
	}
    void addCuts(string id, string seqStr, vector<int> &ctgSig)
	{
		seqId=id;
		seq=seqStr;
		ctgCuts = ctgSig;
	}
	int getContigLen()
	{
		return seq.size();
	}
	int getCutSize()
	{
		return ctgCuts.size();
	}
	int getCutPoint(int i)
	{
		return ctgCuts[i];
	}
	vector<int> getCutVector()
	{
		return ctgCuts;
	}
	string getId()
	{
		return seqId;
	}
	string getSeq()
	{
		return seq;
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
string getRevComp(string s)
{
	int ctgLen=s.length();
	char * cstr = new char [ctgLen+1];
	std::strcpy (cstr, s.c_str());
	char* retChars = new char[ctgLen+1];
	#pragma omp parallel for num_threads(2)
	for(int i=0; i<ctgLen; i++)
	{
		if(cstr[ctgLen-i-1]=='A')
			retChars[i] = 'T'; 
		else if(cstr[ctgLen-i-1]=='C')
			retChars[i] = 'G'; 
		else if(cstr[ctgLen-i-1]=='T')
			retChars[i] = 'A'; 
		else if(cstr[ctgLen-i-1]=='G')
			retChars[i] = 'C';
		else
			retChars[i] = cstr[ctgLen-i-1];
	}
	retChars[ctgLen]='\0';
	string retStr(retChars);
	delete[] retChars;
	delete[] cstr;
	return retStr;
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
void printLongSeq(string str, ostream &os)
{	
	for(int i=0;i<str.size();i+=80)
		os<<str.substr(i,80)<<'\n';
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
string getBestMatch(vector<int> &orList,vector<int> &orOccuList,vector<int> ctgSig)
{
	int beg=0,end=0,minErr=(10000000-1); 
	string ori="NA";
	for(int i=0;i<orList.size();i++)
	{
		//check run length is open
		int runOpen = 1, k=0;
		for(int j=0;j<ctgSig.size();j++)
		{
			if((i+k)<orList.size()){
				if(orOccuList[i+k]==1)
					runOpen=0;
			}
			else{
				if(orOccuList[(i+k)-orList.size()]==1)
					runOpen=0;
			}
			k++;
		}
		if(runOpen==1)
		{
			int accErr=0,j;
			k=0;
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
	}
	//Check Reverse
	for(int i=0;i<orList.size();i++)
	{
		int runOpen = 1, k=0;
		for(int j=0;j<ctgSig.size();j++)
		{
			if((i+k)<orList.size()){
				if(orOccuList[i+k]==1)
					runOpen=0;
			}
			else{
				if(orOccuList[(i+k)-orList.size()]==1)
					runOpen=0;
			}
			k++;
		}
		if(runOpen==1)
		{
			int accErr=0,j;
			k=0;
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
	}
	
	//Mark matched area
	if(beg<=end)
	{
		for(int j=beg;j<end;j++)
			orOccuList[j] = 1;
	}
	else
	{
		for(int j=beg;j<orList.size();j++)
			orOccuList[j] = 1;
		for(int j=0;j<end;j++)
			orOccuList[j] = 1;
	}
	minErr=minErr/ctgSig.size();
	return SSTR(beg) + "," + SSTR(end) + "," + SSTR(minErr)+","+ori;
}
int main(int argc, char **argv)
{
	//Program to find genes/intergenic regions 
	//arg 1: Matched Contigs
	//arg 2: OR Data
	//arg 3: Output matched contigs 	
	//arg 4: Output Draft Genome 	
	string blast_str; 
	map<string, string> ctg_list;
	vector<int> orList;
	vector<int> orOccuList;
	
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
		orOccuList.push_back(0);
	}
	filebuf matchedCtg;
	matchedCtg.open(argv[3],ios::out);
	ostream os(&matchedCtg);
	
	filebuf draftGen;
	draftGen.open(argv[4],ios::out);
	ostream os2(&draftGen);
	
	vector<ContigCuts> ctgSigAll;
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
			if((ctgPos[i+1]-ctgPos[i])>=1000)
				ctgSig.push_back(ctgPos[i+1]-ctgPos[i]);
			else
				ctgPos.erase(ctgPos.begin()+i+1);
			//cout<<ctgPos[i+1]-ctgPos[i]<<"\n";
		}
		ContigCuts temp;
		temp.addCuts(it->first,it->second,ctgSig);
		ctgSigAll.push_back(temp);
	}
	std::sort(ctgSigAll.rbegin(),ctgSigAll.rend(),ContigCuts::compareLenPredicate);
	
	for(int i=0;i<ctgSigAll.size();i++)
	{
		if(ctgSigAll[i].getCutSize()>=5)
		{
			os<<"Processing:"<<ctgSigAll[i].getId()<<"\nCuts Found:"<<ctgSigAll[i].getCutSize()<<"\n";
			for(int k=0;k<ctgSigAll[i].getCutSize();k++)
			{
				os<<ctgSigAll[i].getCutPoint(k)<<",";
			}
			os<<"\n";
			string bestMatch = getBestMatch(orList,orOccuList,ctgSigAll[i].getCutVector());
			os<<"Best Match:"<<bestMatch<<"\n";
			std::vector<std::string> tok = split(bestMatch, ',');
			ctgStart.push_back(atoi(tok[0].c_str()));
			ctgEnd.push_back(atoi(tok[1].c_str()));
			ctgErr.push_back(atoi(tok[2].c_str()));
			ctgOri.push_back(tok[3]);
		}
	}
	os<<"Contig Cut Bitmap\n";
	string cutPintCons(orList.size(),'0');
	vector<int> cutPintIndx(orList.size(),-1);
	for(int i=0;i<ctgStart.size();i++)
	{
		string dispStr(orList.size(),'0');
		os<<ctgSigAll[i].getId()<<":";
		if(ctgStart[i]<=ctgEnd[i])
		{
			for(int j=ctgStart[i];j<ctgEnd[i];j++)
			{
				dispStr[j] = '1';
				cutPintCons[j] = '1';
				cutPintIndx[j] = i;
			}
		}
		else
		{
			for(int j=ctgStart[i];j<orList.size();j++)
			{
				dispStr[j] = '1';
				cutPintCons[j] = '1';
				cutPintIndx[j] = i;
			}
			for(int j=0;j<ctgEnd[i];j++)
			{
				dispStr[j] = '1';
				cutPintCons[j] = '1';
				cutPintIndx[j] = i;
			}
		}
		os<<dispStr<<"\n";
	}
	
	os<<"All Matched Cuts:"<<cutPintCons<<"\n";
	int totalMatchedCuts=0;
	for(int j=0;j<orList.size();j++)
		totalMatchedCuts += (int)cutPintCons[j]-'0';
	os<<"Total Cuts Matched="<<totalMatchedCuts<<"\n";
	
	int coveredSeq=0;
	string prevStr="";
	vector<int> compressedIndx;
	for(int i=0;i<cutPintIndx.size();i++)
	{
		if(cutPintIndx[i]>=0)
		{
			if(ctgSigAll[cutPintIndx[i]].getId()!=prevStr)
			{
				prevStr=ctgSigAll[cutPintIndx[i]].getId();
				os<<prevStr<<"("<<ctgOri[cutPintIndx[i]]<<")\n";
				coveredSeq += ctgSigAll[cutPintIndx[i]].getContigLen();
				compressedIndx.push_back(cutPintIndx[i]);
			}
		}
	}
	
	
	
	string genomeSeq="";
	string gap(3000,'N');
	for(int i=0;i<compressedIndx.size()-1;i++)
	{
		string tempStr;
		if(ctgOri[compressedIndx[i]]=="REV")
			tempStr=getRevComp(ctgSigAll[compressedIndx[i]].getSeq());
		else
			tempStr=ctgSigAll[compressedIndx[i]].getSeq();
		
		genomeSeq += tempStr;
		
		if(i<(compressedIndx.size()-2))
			genomeSeq+=gap;
		 
	}
	os<<"Sequence Covered="<<genomeSeq.size()-(3000*(compressedIndx.size()-1))<<"\n";
	os2<<"> MSho Draft Genome|Length="<<genomeSeq.size()<<"\n";
	printLongSeq(genomeSeq,os2);
	
	matchedCtg.close();
	orFile.close();
	ctgFile.close();
	draftGen.close();
	
	return 0;	
}	
