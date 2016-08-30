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

class CtgJoins
{
	private:
	string ctgID1,ctgID2;
	int blastS1,blastS2,blastE1,blastE2;
	int joinTyp; //1 MP Only; 2 MP and IS 
	string isID;
	int isS1,isE1,isS2,isE2;
	string sy1,sy2; 
	double mpCount, meadDist, stdDev;
	
	public:
	CtgJoins()
	{
		ctgID1 = "";
		blastS1 = -1;
		blastE1 = -1;
		ctgID2 = "";
		blastS2 = -1;
		blastE2 = -1;
		isID = "";
		isS1 = -1;
		isE1 = -1;
		isS2 = -1;
		isE2 = -1;
		sy1 = "";
		sy2 = "";
		mpCount = -1;
		meadDist = -1;
		stdDev = -1;
	} 
	string getC1Id() {return ctgID1;}
	string getC2Id() {return ctgID2;}
	double getMates() {return mpCount;}
	
	void addJoinType1(string jstr, int firstCtg)
	{
		joinTyp = 1;
		if(firstCtg==1)
		{
			std::vector<std::string> tok = split(jstr, '|');
			ctgID1 = tok[0];
			blastS1 = atoi(tok[1].c_str());
			blastE1 = atoi(tok[2].c_str());
			ctgID2 = tok[3];
			blastS2 = atoi(tok[4].c_str());
			blastE2 = atoi(tok[5].c_str());
			mpCount = atof(tok[6].c_str());
			meadDist = atof(tok[7].c_str());
			stdDev = atof(tok[8].c_str());
		}
		else
		{
			std::vector<std::string> tok = split(jstr, '|');
			ctgID1 = tok[3];
			blastS1 = atoi(tok[4].c_str());
			blastE1 = atoi(tok[5].c_str());
			ctgID2 = tok[0];
			blastS2 = atoi(tok[1].c_str());
			blastE2 = atoi(tok[2].c_str());
			
			mpCount = atof(tok[6].c_str());
			meadDist = atof(tok[7].c_str());
			stdDev = atof(tok[8].c_str());
		}
		
	}
	void addJoinType2(string jstr,int firstCtg)
	{
		joinTyp = 2;
		if(firstCtg==1)
		{
			std::vector<std::string> tok = split(jstr, '|');
			ctgID1 = tok[0];
			blastS1 = atoi(tok[1].c_str());
			blastE1 = atoi(tok[2].c_str());
			ctgID2 = tok[3];
			blastS2 = atoi(tok[4].c_str());
			blastE2 = atoi(tok[5].c_str());
			isID = tok[6];
			isS1 = atoi(tok[9].c_str());
			isE1 = atoi(tok[10].c_str());
			isS2 = atoi(tok[7].c_str());
			isE2 = atoi(tok[8].c_str());
			mpCount = atof(tok[11].c_str());
			meadDist = atof(tok[12].c_str());
			stdDev = atof(tok[13].c_str());
		}
		else
		{
			std::vector<std::string> tok = split(jstr, '|');
			ctgID1 = tok[3];
			blastS1 = atoi(tok[4].c_str());
			blastE1 = atoi(tok[5].c_str());
			ctgID2 = tok[0];
			blastS2 = atoi(tok[1].c_str());
			blastE2 = atoi(tok[2].c_str());
			isID = tok[6];
			isS1 = atoi(tok[7].c_str());
			isE1 = atoi(tok[8].c_str());
			isS2 = atoi(tok[9].c_str());
			isE2 = atoi(tok[10].c_str());
			mpCount = atof(tok[11].c_str());
			meadDist = atof(tok[12].c_str());
			stdDev = atof(tok[13].c_str());
		}
	}
	/*void addJoinType3(string jstr,int firstCtg)
	{
		joinTyp = 3;
		if(firstCtg==1)
		{
			std::vector<std::string> tok = split(jstr, '|');
			ctgID1 = tok[0];
			blastS1 = atoi(tok[1].c_str());
			blastE1 = atoi(tok[2].c_str());
			ctgID2 = tok[3];
			blastS2 = atoi(tok[4].c_str());
			blastE2 = atoi(tok[5].c_str());
			isID = tok[6];
			isS1 = atoi(tok[9].c_str());
			isE1 = atoi(tok[10].c_str());
			isS2 = atoi(tok[7].c_str());
			isE2 = atoi(tok[8].c_str());
			sy1 = tok[11];
			sy2 = tok[12];
			mpCount = atof(tok[13].c_str());
			meadDist = atof(tok[14].c_str());
			stdDev = atof(tok[15].c_str());
		}
		else
		{
			std::vector<std::string> tok = split(jstr, '|');
			ctgID1 = tok[3];
			blastS1 = atoi(tok[4].c_str());
			blastE1 = atoi(tok[5].c_str());
			ctgID2 = tok[0];
			blastS2 = atoi(tok[1].c_str());
			blastE2 = atoi(tok[2].c_str());
			isID = tok[6];
			isS1 = atoi(tok[7].c_str());
			isE1 = atoi(tok[8].c_str());
			isS2 = atoi(tok[9].c_str());
			isE2 = atoi(tok[10].c_str());
			sy1 = tok[11];
			sy2 = tok[12];
			mpCount = atof(tok[13].c_str());
			meadDist = atof(tok[14].c_str());
			stdDev = atof(tok[15].c_str());
		}
	}*/
};

class Contigs
{
	private: 
		string seqId;
		string seq;
		map<string,CtgJoins> beginJoinCtg;
		map<string,CtgJoins> endJoinCtg;
	public:
		void setSeqId(string s) { seqId=s;}
		string getSeqId() {return seqId;}
		void setSeq(string s) { seq=s;}
		string getSeq() {return seq;}
		int getLength() { return seq.size(); };
		void addJoin(string bStr, int baseCtg, string jp,int jType)
		{
			string ctg="";
			std::vector<std::string> tok = split(bStr, '|');
			if(baseCtg==1)
				ctg=tok[3];
			else
				ctg=tok[0];
				
			if(jp=="B")
			{
				if(jType==1)
				{
					CtgJoins ctgJoin;
					ctgJoin.addJoinType1(bStr,baseCtg);
					beginJoinCtg.insert(pair<string,CtgJoins>(ctg,ctgJoin));
				}
				else if(jType==2)
				{
					CtgJoins ctgJoin;
					ctgJoin.addJoinType2(bStr,baseCtg);
					beginJoinCtg.insert(pair<string,CtgJoins>(ctg,ctgJoin));
				}
				else
				{
					CtgJoins ctgJoin;
					ctgJoin.addJoinType3(bStr,baseCtg);
					beginJoinCtg.insert(pair<string,CtgJoins>(ctg,ctgJoin));
				}
			}
			else
			{
				if(jType==1)
				{
					CtgJoins ctgJoin;
					ctgJoin.addJoinType1(bStr,baseCtg);
					endJoinCtg.insert(pair<string,CtgJoins>(ctg,ctgJoin));
				}
				else if(jType==2)
				{
					CtgJoins ctgJoin;
					ctgJoin.addJoinType2(bStr,baseCtg);
					endJoinCtg.insert(pair<string,CtgJoins>(ctg,ctgJoin));
				}
				else
				{
					CtgJoins ctgJoin;
					ctgJoin.addJoinType3(bStr,baseCtg);
					endJoinCtg.insert(pair<string,CtgJoins>(ctg,ctgJoin));
				}
			}
		}
		CtgJoins getBeginJoin(vector<string> &alreadyJoinCtg)
		{
			CtgJoins retJoin;
			for(map<string, CtgJoins>::iterator it = beginJoinCtg.begin(); it != beginJoinCtg.end(); it++) 
			{
				std::vector<string>::iterator itStr;
				itStr = find(alreadyJoinCtg.begin(), alreadyJoinCtg.end(), it->first);
				if(itStr == alreadyJoinCtg.end())
				{
					if(retJoin.getMates()<it->second.getMates())
						retJoin = it->second;
				}
			}
			return retJoin;
		}
		CtgJoins getEndJoin(vector<string> &alreadyJoinCtg)
		{
			CtgJoins retJoin;
			for(map<string, CtgJoins>::iterator it = endJoinCtg.begin(); it != endJoinCtg.end(); it++) 
			{
				std::vector<string>::iterator itStr;
				itStr = find(alreadyJoinCtg.begin(), alreadyJoinCtg.end(), it->first);
				if(itStr == alreadyJoinCtg.end())
				{
					if(retJoin.getMates()<it->second.getMates())
						retJoin = it->second;
				}
			}
			return retJoin;
		}
};
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

void readContigFasta(ifstream& multiFile, map<string,Contigs> &theMap)
{
	string str,id,seq; 
	getline(multiFile,str);
	id=str.substr(1); 
	while(getline(multiFile,str))
	{
		int found = str.find(">");
		if (found!=string::npos)
		{
			Contigs ctgs;
			ctgs.setSeqId(id);
			ctgs.setSeq(seq);			
			theMap.insert (pair<string,Contigs>(id,ctgs) );
			seq="";
			id = str.substr(1);
		}
		else
		{
			seq=seq+str;
		}
	}
	Contigs ctgs;
	ctgs.setSeqId(id);
	ctgs.setSeq(seq);	
	theMap.insert(pair<string,Contigs>(id,ctgs));
}

void readISFasta(ifstream& multiFile, map<string,string> &theMap)
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
	//arg 1: Mated and IS Matched Contigs
	//arg 2: All Mated Contigs
	//arg 3: IS Sequences
	//arg 4: Contig Sequences
	//arg 5: Output matched contigs 	
	string blast_str; 
	map<string, Contigs> ctg_list;
	map<string,string> is_list;
	map<string,int> allGoodCtgJoins;
	
	ifstream matedISCtgs(argv[1]);
	ifstream matedCtgs(argv[2]);
	ifstream isSeqs(argv[3]);
	ifstream ctgSeqs(argv[4]);
			
	if(matedISCtgs.fail() || matedCtgs.fail() || isSeqs.fail() || ctgSeqs.fail())
	{
		cerr<<"\nThe files could not be opened.";
		return -1;
	}
	readISFasta(isSeqs,is_list);
	readContigFasta(ctgSeqs, ctg_list);
	
	while(getline(matedISCtgs,blast_str))
	{
		std::vector<std::string> tok = split(blast_str, '|');
		
		map<string,Contigs>::iterator itCtg1;
		itCtg1 = ctg_list.find(tok[0]);
		string c1ORI = ((atoi(tok[1].c_str()) < 10) || (atoi(tok[2].c_str()) < 10))?"B":"E";
		
		itCtg1->second.addJoin(blast_str,1,c1ORI,2);
		
		map<string,Contigs>::iterator itCtg2;
		itCtg2 = ctg_list.find(tok[3]);
		string c2ORI = ((atoi(tok[4].c_str()) < 10) || (atoi(tok[5].c_str()) < 10))?"B":"E";
		itCtg2->second.addJoin(blast_str,2,c2ORI,2);
	}
	
	while(getline(matedCtgs,blast_str))
	{
		std::vector<std::string> tok = split(blast_str, '|');
		
		map<string,Contigs>::iterator itCtg1;
		itCtg1 = ctg_list.find(tok[0]);
		string c1ORI = ((atoi(tok[1].c_str()) < 10) || (atoi(tok[2].c_str()) < 10))?"B":"E";
		itCtg1->second.addJoin(blast_str,1,c1ORI,1);
		
		map<string,Contigs>::iterator itCtg2;
		itCtg2 = ctg_list.find(tok[3]);
		string c2ORI = ((atoi(tok[4].c_str()) < 10) || (atoi(tok[5].c_str()) < 10))?"B":"E";
		itCtg2->second.addJoin(blast_str,2,c2ORI,1);
	}
	

	/*while(getline(allMatchMates,blast_str))
	{
		std::vector<std::string> tok = split(blast_str, '|');
		
		string c1ORI = ((atoi(tok[1].c_str()) < 10) || (atoi(tok[2].c_str()) < 10))?"B":"E";
		string c2ORI = ((atoi(tok[4].c_str()) < 10) || (atoi(tok[5].c_str()) < 10))?"B":"E";
		
		allGoodCtgJoins.insert(pair<string,int>(tok[0]+c1ORI,atoi(tok[11].c_str())));
		allGoodCtgJoins.insert(pair<string,int>(tok[3]+c2ORI,atoi(tok[11].c_str())));
	}*/
	
	//double start = omp_get_wtime();
	//double end = omp_get_wtime( );
	//printf("start = %.16g\nend = %.16g\ndiff = %.16g\n",start, end, end - start);
	
	filebuf joinedCtg;
	joinedCtg.open(argv[5],ios::out);
	ostream os(&joinedCtg);
	
	vector<string> allJoinedCtg;
	for(map<string,Contigs>::iterator it = ctg_list.begin(); it != ctg_list.end(); it++)
	{
		std::vector<string>::iterator itStr;
		itStr = find(allJoinedCtg.begin(), allJoinedCtg.end(), it->first);
		if(itStr == allJoinedCtg.end())
		{
			string scfString = it->second.getSeq();
			string scfStringID = it->first;
			int scfLen = it->second.getLength();
			allJoinedCtg.push_back(it->first);
			CtgJoins begJoin = it->second.getBeginJoin(allJoinedCtg);
			CtgJoins endJoin = it->second.getEndJoin(allJoinedCtg);
			
			while(begJoin.getC1Id() != "")
			{
				allJoinedCtg.push_back(begJoin.getC2Id());
				scfStringID = begJoin.getC2Id() + "|" + scfStringID;
				
				map<string,Contigs>::iterator itCtg;
				itCtg = ctg_list.find(begJoin.getC2Id());
				
				scfLen += itCtg->second.getLength();
				begJoin = itCtg->second.getBeginJoin(allJoinedCtg);
			}
			while(endJoin.getC1Id() != "")
			{
				allJoinedCtg.push_back(endJoin.getC2Id());
				scfStringID =  scfStringID + "|" +  endJoin.getC2Id(); 
				
				map<string,Contigs>::iterator itCtg;
				itCtg = ctg_list.find(endJoin.getC2Id());
				scfLen += itCtg->second.getLength();
				endJoin = itCtg->second.getEndJoin(allJoinedCtg);
				
			}
			cout<<scfStringID<<"\n";
		}
	}
		
		
	/*while(getline(matchedCtgs,blast_str))
	{
		std::vector<std::string> tok = split(blast_str, '|');
		map<string,Contigs>::iterator itCtg1;
		itCtg1 = ctg_list.find(tok[0]);
		map<string,Contigs>::iterator itCtg2;
		itCtg2 = ctg_list.find(tok[3]);
		
		//Join Points Begin0 End1
		string c1ORI = ((atoi(tok[1].c_str()) < 10) || (atoi(tok[2].c_str()) < 10))?"B":"E";
		string c2ORI = ((atoi(tok[4].c_str()) < 10) || (atoi(tok[5].c_str()) < 10))?"B":"E";
		
		
		//Check mated and decide to further processing; It better match exists -- continue;
		if(atoi(tok[13].c_str()) == 0)
		{
			map<string,int>::iterator itAllG;
			itAllG = allGoodCtgJoins.find(tok[0]+c1ORI);
			if(itAllG != allGoodCtgJoins.end())
				continue;
			
			itAllG = allGoodCtgJoins.find(tok[3]+c2ORI);
			if(itAllG != allGoodCtgJoins.end())
				continue;
		}
		
		//Orientation FR
		c1ORI += (atoi(tok[1].c_str()) < atoi(tok[2].c_str()))?"F":"R";
		c2ORI += (atoi(tok[4].c_str()) < atoi(tok[5].c_str()))?"F":"R";
		
		
		map<string,string>::iterator itIS;
		itIS = is_list.find(tok[6]);
		
		string scfStr="";
				
		if(c1ORI=="EF" && c2ORI=="BF") //EF -- BF F-IS-F Tested
		{
			scfStr = itCtg1->second.getSeq();
			scfStr += itIS->second.substr(atoi(tok[10].c_str()),atoi(tok[7].c_str())-atoi(tok[10].c_str())-1);
			scfStr += itCtg2->second.getSeq();
		}
		else if(c1ORI=="BF" && c2ORI=="EF") //EF -- BF F-IS-F Tested
		{
			scfStr = itCtg2->second.getSeq();
			scfStr += itIS->second.substr(atoi(tok[8].c_str()),atoi(tok[9].c_str())-atoi(tok[8].c_str())-1);
			scfStr += itCtg1->second.getSeq();
		}
		
		else if(c1ORI=="EF" && c2ORI=="ER") //EF -- ER F-IS-R Tested
		{
			scfStr = itCtg1->second.getSeq();
			scfStr += itIS->second.substr(atoi(tok[10].c_str()),atoi(tok[7].c_str())-atoi(tok[10].c_str())-1);
			scfStr += getRevComp(itCtg2->second.getSeq());
		}
		else if(c1ORI=="ER" && c2ORI=="EF") //EF -- ER F-IS-R Tested
		{
			scfStr = itCtg2->second.getSeq();
			scfStr += itIS->second.substr(atoi(tok[8].c_str()),atoi(tok[9].c_str())-atoi(tok[8].c_str())-1);
			scfStr += getRevComp(itCtg1->second.getSeq());
		}
		
		else if(c1ORI=="BR" && c2ORI=="BF") //BR -- BF -> R-IS-F Tested
		{
			scfStr = getRevComp(itCtg1->second.getSeq());
			scfStr += itIS->second.substr(atoi(tok[10].c_str()),atoi(tok[7].c_str())-atoi(tok[10].c_str())-1);
			scfStr += itCtg2->second.getSeq();
		}
		else if(c1ORI=="BF" && c2ORI=="BR") //BR -- BF -> R-IS-F Tested
		{
			scfStr = getRevComp(itCtg2->second.getSeq());
			scfStr += itIS->second.substr(atoi(tok[8].c_str()),atoi(tok[9].c_str())-atoi(tok[8].c_str())-1);
			scfStr += itCtg1->second.getSeq();
		}
		
		else if(c1ORI=="BR" && c2ORI=="ER") //BR -- ER -> ER -- BR F-ISR-F Tested
		{
			scfStr = itCtg2->second.getSeq();
			scfStr += getRevComp(itIS->second.substr(atoi(tok[10].c_str()),atoi(tok[7].c_str())-atoi(tok[10].c_str())-1));
			scfStr += itCtg1->second.getSeq();
		}
		else if(c1ORI=="ER" && c2ORI=="BR") //BR -- ER -> ER -- BR F-ISR-F Tested
		{
			scfStr = itCtg1->second.getSeq();
			scfStr += getRevComp(itIS->second.substr(atoi(tok[8].c_str()),atoi(tok[9].c_str())-atoi(tok[8].c_str())-1));
			scfStr += itCtg2->second.getSeq();
		}
		os<<">"<<blast_str<<"\n"<<scfStr<<"\n";
	}*/
	joinedCtg.close();
	isSeqs.close();
	ctgSeqs.close();
	matedISCtgs.close();
	matedCtgs.close();
	return 0;	
}	
