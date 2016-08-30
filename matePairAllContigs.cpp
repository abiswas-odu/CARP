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

class Contigs
{
	private: 
		string seqId;
		int seqLen;
		vector<string> readId;
		vector<int> begPos;
		vector<int> endPos;
		vector<int> oriPos;
	public:
		void setSeqId(string s) { seqId=s;}
		string getSeqId() {return seqId;}
		void addReads(string readid,int st,int en,int ori)
		{
			readId.push_back(readid);
			begPos.push_back(st);
			endPos.push_back(en);
			oriPos.push_back(ori);
		}
		int getLength() { return seqLen; };
		void setLength(map<string,int> &ctgLen)
		{
			map<string,int>::iterator itCTG;
			itCTG = ctgLen.find("ctg"+seqId);
			seqLen = itCTG->second;
		}
		int getReadDist(int joinPt, int indx)
		{
			int dist=0;
			
			if(joinPt==0)
				dist = begPos[indx];
			else
				dist = seqLen - endPos[indx];
				
			return dist;
		}
		friend double mateContigs(Contigs c1, int c1FR, int c1jp, Contigs c2, int c2FR, int c2jp, map<string,string> &ctgReadMap, int midLen, double &sumDist, double &sumSqDist);
};

double mateContigs(Contigs c1, int c1FR, int c1jp, Contigs c2, int c2FR, int c2jp, map<string,string> &ctgReadMap, int midLen, double& sumDist, double& sumSqDist)
{
	int countMates=0;
	int sumMateDist=0;
	int sumSqMateDist=0;
	#pragma omp parallel for
	for(int i=0;i<c1.readId.size();i++)
	{
		//cout<<"Entering Parallel Code:"<<omp_get_num_threads()<<" threads and "<<omp_get_thread_num()<<" thread."<<'\n';
		map<string,string>::iterator itCtg1;
		itCtg1 = ctgReadMap.find(c1.readId[i]);
		for(int j=0;j<c2.readId.size();j++)
		{
			map<string,string>::iterator itCtg2;
			itCtg2 = ctgReadMap.find(c2.readId[j]);
			
			if(itCtg1 != ctgReadMap.end() && itCtg2 != ctgReadMap.end())
			{
				if(itCtg1->second == itCtg2->second)// && )
				{
					//F-F or R-R orientation
					if((c1FR==0 && c2FR==0) || (c1FR==1 && c2FR==1))
					{
						if(c1.oriPos[i]==c2.oriPos[j])
						{
							int mateDist = midLen + c1.getReadDist(c1jp,i) + c2.getReadDist(c2jp,j);
							if(mateDist<=6000)
							{
								#pragma omp atomic 
									sumSqMateDist+=(mateDist*mateDist);
								#pragma omp atomic 
									sumMateDist+=mateDist;
								#pragma omp atomic 
									countMates++;
							}
						}
					}
					else
					{
						if(c1.oriPos[i] != c2.oriPos[j])
						{
							int mateDist = midLen + c1.getReadDist(c1jp,i) + c2.getReadDist(c2jp,j);
							if(mateDist<=6000)
							{
								#pragma omp atomic 
									sumSqMateDist+=(mateDist*mateDist);
								#pragma omp atomic 
									sumMateDist+=mateDist;
								#pragma omp atomic 
									countMates++;
							}
						}
					}
				}
			}
		}
	}
	sumDist = sumMateDist;
	sumSqDist = sumSqMateDist;
	return countMates;
}
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

void readContigPosmatch(ifstream& multiFile, map<string, Contigs> &ctgs, map<string,string> &ctgReadMap, map<string,int> &ctgLen)
{
	string str,ctgId; 
	getline(multiFile,str);
	std::vector<std::string> tok = split(str, '\t');
	ctgId=tok[1]; 
	Contigs* gen = new Contigs();
	gen->setSeqId(ctgId);
	
	map<string,string>::iterator itCtg;
	itCtg = ctgReadMap.find(tok[0]);
	if(itCtg != ctgReadMap.end())
	{
		int ori = tok[4]=="f"?0:1;
		gen->addReads(tok[0],atoi(tok[2].c_str()),atoi(tok[3].c_str()),ori);
	}
	while(getline(multiFile,str))
	{
		tok = split(str, '\t');
		if (ctgId==tok[1])
		{
			itCtg = ctgReadMap.find(tok[0]);
			if(itCtg != ctgReadMap.end())
			{
				int ori = tok[4]=="f"?0:1;
				gen->addReads(tok[0],atoi(tok[2].c_str()),atoi(tok[3].c_str()),ori);
			}
		}
		else
		{
			(*gen).setLength(ctgLen);
			ctgs.insert(pair<string, Contigs>(ctgId,(*gen)));
			delete gen;
			ctgId=tok[1];
			gen = new Contigs();
			gen->setSeqId(ctgId);
			itCtg = ctgReadMap.find(tok[0]);
			if(itCtg != ctgReadMap.end())
			{
				int ori = tok[4]=="f"?0:1;
				gen->addReads(tok[0],atoi(tok[2].c_str()),atoi(tok[3].c_str()),ori);
			}
		}
	}
	(*gen).setLength(ctgLen);
	ctgs.insert(pair<string, Contigs>(ctgId,(*gen)));
	delete gen;
}
void readContigFasta(ifstream& multiFile, map<string,int> &theMap)
{
	string str,id,seq; 
	getline(multiFile,str);
	id=str.substr(1); 
	while(getline(multiFile,str))
	{
		int found = str.find(">");
		if (found!=string::npos)
		{
			theMap.insert (pair<string,int>(id,seq.size()) );
			seq="";
			id = str.substr(1);
		}
		else
		{
			seq=seq+str;
		}
	}
	theMap.insert(pair<string,int>(id,seq.size()));
}

void readISFasta(ifstream& multiFile, map<string,int> &theMap)
{
	string str,id,seq; 
	getline(multiFile,str);
	id=str.substr(1); 
	while(getline(multiFile,str))
	{
		int found = str.find(">");
		if (found!=string::npos)
		{
			theMap.insert (pair<string,int>(id,seq.size()) );
			seq="";
			id = str.substr(1);
		}
		else
		{
			seq=seq+str;
		}
	}
	theMap.insert (pair<string,int>(id,seq.size()) );
}

int main(int argc, char **argv)
{
	//Program to find genes/intergenic regions 
	//arg 1: Reads mapped contigs (posmatch.frgctg)
	//arg 2: Mate Read ID map (gkpStore.fastqUIDmap)
	//arg 3: List of Contigs
	//arg 4: Contig Sequences
	//arg 5: Output matched contigs 	
	string blast_str; 
	map<string, Contigs> ctg_list;
	map<string,string> ctgReadMap;
	vector<string> matchedCtg_list;
	map<string,int> ctgLen;

	ifstream ctgFile(argv[1]);
	ifstream mateReads(argv[2]);
	ifstream ctgList(argv[3]);
	ifstream ctgSeqs(argv[4]);
	
	if(ctgFile.fail() || mateReads.fail() || ctgList.fail() || ctgSeqs.fail())
	{
		cerr<<"\nThe files could not be opened.";
		return -1;
	}
	readContigFasta(ctgSeqs, ctgLen);
	while(getline(ctgList,blast_str))
	{
		matchedCtg_list.push_back(blast_str);
	}
	
	while(getline(mateReads,blast_str))
	{
		std::vector<std::string> tok = split(blast_str, '\t');
		ctgReadMap.insert(pair<string,string>(tok[0],tok[2].substr(0,tok[2].size()-2)));
	}
	readContigPosmatch(ctgFile,ctg_list,ctgReadMap, ctgLen);
	
	filebuf matedCtg;
	matedCtg.open(argv[5],ios::out);
	ostream os(&matedCtg);
	
	double start = omp_get_wtime();
	for(int i=0;i<matchedCtg_list.size();i++)
	{
		for(int j=i+1;j<matchedCtg_list.size();j++)
		{
			//Join Orient FWD0 REV1
			//Join Points Begin0 End1
			
			map<string,Contigs>::iterator itCtg1;
			itCtg1 = ctg_list.find(matchedCtg_list[i].substr(3));
			map<string,Contigs>::iterator itCtg2;
			itCtg2 = ctg_list.find(matchedCtg_list[j].substr(3));
						
			//Join (a,b)(c,d) --> (b,c) --(0,1)(0,0)
			double sumDist = 0.0;
			double sumSqDist = 0.0;
			double n = mateContigs(itCtg1->second,0,1,itCtg2->second,0,0,ctgReadMap,0,sumDist,sumSqDist);
			double avg = 0.0;
			double stdDev = 0.0;
			if(n>0.0)
			{
				avg = sumDist/n;
				stdDev = sqrt((sumSqDist-((sumDist*sumDist)/n))/n);
				os<<matchedCtg_list[i]<<"|0|1|"<<matchedCtg_list[j]<<"|0|0|"<<n<<"|"<<avg<<"|"<<stdDev<<"\n";
			}
			//Join (a,b)(c,d) --> (b,d) --(0,1)(1,1)
			
			sumDist = 0.0;
			sumSqDist = 0.0;
			n = mateContigs(itCtg1->second,0,1,itCtg2->second,1,1,ctgReadMap,0,sumDist,sumSqDist);
			avg = 0.0;
			stdDev = 0.0;
			if(n>0.0)
			{
				avg = sumDist/n;
				stdDev = sqrt((sumSqDist-((sumDist*sumDist)/n))/n);
				os<<matchedCtg_list[i]<<"|0|1|"<<matchedCtg_list[j]<<"|1|1|"<<n<<"|"<<avg<<"|"<<stdDev<<"\n";
			}
			
			//Join (a,b)(c,d) --> (a,c) --(1,0)(0,0)
			
			sumDist = 0.0;
			sumSqDist = 0.0;
			n = mateContigs(itCtg1->second,1,0,itCtg2->second,0,0,ctgReadMap,0,sumDist,sumSqDist);
			avg = 0.0;
			stdDev = 0.0;
			if(n>0.0)
			{
				avg = sumDist/n;
				stdDev = sqrt((sumSqDist-((sumDist*sumDist)/n))/n);
				os<<matchedCtg_list[i]<<"|0|1|"<<matchedCtg_list[j]<<"|1|1|"<<n<<"|"<<avg<<"|"<<stdDev<<"\n";
			}
			
			//Join (a,b)(c,d) --> (a,d) --(0,0)(0,1)
			
			sumDist = 0.0;
			sumSqDist = 0.0;
			n = mateContigs(itCtg1->second,0,0,itCtg2->second,0,1,ctgReadMap,0,sumDist,sumSqDist);
			avg = 0.0;
			stdDev = 0.0;
			if(n>1.0)
			{
				avg = sumDist/n;
				stdDev = sqrt((sumSqDist-((sumDist*sumDist)/n))/n);
				os<<matchedCtg_list[i]<<"|0|0|"<<matchedCtg_list[j]<<"|0|1|"<<n<<"|"<<avg<<"|"<<stdDev<<"\n";
			}
		}
	}
	double end = omp_get_wtime( );
	
	printf("start = %.16g\nend = %.16g\ndiff = %.16g\n",start, end, end - start);
	
	matedCtg.close();
	ctgFile.close();
	mateReads.close();
	ctgList.close();
	ctgSeqs.close();
	return 0;	
}	
