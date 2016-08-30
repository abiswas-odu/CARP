#include<map>
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<vector>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<omp.h>
#include<cstring>
#include<algorithm>

#define MAX_TRANSHIT_LEN 100
#define MAX_SEED_LENGTH 3000
#define MAX_DIFF_HAMMING_DIST 5

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
string getRevComp(string s)
{
	int ctgLen=s.length();
	char * cstr = new char [ctgLen+1];
	std::strcpy(cstr, s.c_str());
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
bool isCutOverlap(int s,int e,int s1, int e1)
{
	if(s >= s1 && s <= e1)
     return true;
	if(s1 >= s && s1 <= e)
     return true;
	return false;
}
map<string,string> getTransHitSeqs(vector<string>& isRelatedCtgBLASTResult, ifstream& CtgFile)
{
    multimap<string,string> ctg_cuts;
	map<string,string> ctg_list;
	map<string,string> transHit_list;
    map<string,string> retSeqs;
	map<string,string> transHitContigCoordMap;
    readFasta(CtgFile,ctg_list);
	for(int i=0;i<isRelatedCtgBLASTResult.size();i++)
	{
		std::vector<std::string> tok = split(isRelatedCtgBLASTResult[i], ',');
		string ctgID =  tok[0];
		int endIndx = tok.size();
		int qstart = atoi(tok[endIndx-4].c_str());
		int qend = atoi(tok[endIndx-3].c_str());

		std::multimap<string,string>::iterator iter=ctg_cuts.find(ctgID);
		
		if(iter != ctg_cuts.end())
		{
			pair<multimap<string, string>::iterator, multimap<string, string>::iterator> ppp;
			ppp = ctg_cuts.equal_range(ctgID);
			int overlapFound=0;
			for (multimap<string,string>::iterator it2 = ppp.first; it2 != ppp.second; ++it2)
			{
				string cutRange = (*it2).second;
				std::vector<std::string> tok2 = split(cutRange, ',');
				int prev_qstart = atoi(tok2[0].c_str());
				int prev_qend = atoi(tok2[1].c_str());
				//Check if there is overlap with prev contig cut ranges
				if(isCutOverlap(qstart,qend,prev_qstart,prev_qend))
				{
					overlapFound=1;
					cout<<"Overlap Found..."<<ctgID<<","<<cutRange<<"\n";
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
					it2->second = cutRange;
					break;
				}
			}
			//If the range is a another new one for this contig...
			if(overlapFound == 0)
			{
				string cutRange =  SSTR(qstart) + "," + SSTR(qend);
				ctg_cuts.insert(pair<string,string>(ctgID,cutRange));
			}
		}
		else
		{
			string cutRange =  SSTR(qstart) + "," + SSTR(qend);
			ctg_cuts.insert(pair<string,string>(ctgID,cutRange));
		}
	}
	int hitSeqCtr=0;
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

		if(transHit.size() > MAX_TRANSHIT_LEN)
		{
			//avoid duplicate tranposase hits even from different contigs...
			std::map<string,string>::iterator iter2=transHit_list.find(transHit);
			if(iter2 == transHit_list.end())
			{
				retSeqs.insert(pair<string,string>("transposaseHit_"+SSTR(hitSeqCtr),transHit));
                transHit_list.insert(pair<string,string>(transHit,key));
				transHitContigCoordMap.insert(pair<string,string>("transposaseHit_"+SSTR(hitSeqCtr),key+","+ SSTR(qstart) +","+SSTR(qend)));
				hitSeqCtr++;
			}
		}
    }
	//Write out the mapping of tranposase hit to contig ID
	ofstream coordFile;
	string coordFileName = "IS_Contig_Map.txt";
	coordFile.open(coordFileName.c_str());
	for(std::map<string,string>::iterator it=transHitContigCoordMap.begin(); it!=transHitContigCoordMap.end(); ++it)
	{
		coordFile<<it->first<<","<<it->second<<"\n";
	}
	coordFile.close();
    return retSeqs;
}

void getRepeatSeqs(string hitSeqName,string hitSeq, map<string,string> &rawReads, multimap<string,string> &blastHits,vector<string> &beginMatches,vector<string> &endMatches)
{
    pair<multimap<string, string>::iterator, multimap<string, string>::iterator> ppp;
	ppp = blastHits.equal_range(hitSeqName);
	for (multimap<string,string>::iterator it2 = ppp.first; it2 != ppp.second; ++it2)
	{
		string blast_str = (*it2).second;
		std::vector<std::string> tok = split(blast_str, ',');
		string readID =  tok[0];
		int qstart = atoi(tok[2].c_str());
		int qend = atoi(tok[3].c_str());

		int hstart = atoi(tok[4].c_str());
		int hend = atoi(tok[5].c_str());

		std::map<string,string>::iterator iter=rawReads.find(readID);
		string readSeq = (*iter).second;
		
		//Typical begin hit...
		if(hstart == 1){
			if(qstart>1 && qend==readSeq.size()){
				//isBorderHit = isFwd = isBeginBorder = 1;
				beginMatches.push_back(readSeq.substr(0,qstart));
			}
		}
		//Reverse hit at the rear
		else if(hstart==hitSeq.size()){
			if(qstart==1 && qend<readSeq.size()){
				//isBorderHit = 1;
				endMatches.push_back(getRevComp(readSeq.substr(qend)));
			}
		}
		//
		else if(hend==1){
			if(qstart>1 && qend==readSeq.size()){
				//isBorderHit = isBeginBorder = 1;
				beginMatches.push_back(getRevComp(readSeq.substr(0,qstart)));
			}
		}
		//Typical end hit...
		else if(hend==hitSeq.size()){
			if(qstart==1 && qend<readSeq.size()){
				//isBorderHit = isFwd = 1;
				endMatches.push_back(readSeq.substr(qend));
			}
		}
    }
}

void generateFASTA(string fastaFileName, vector<string> &seqs,string baseName)
{
	ofstream fastaFile;
	fastaFile.open(fastaFileName.c_str());
	for(int i=0;i<seqs.size();i++)
	{
		fastaFile<<">"<<baseName<<"_"<<i<<"\n";
		fastaFile<<seqs[i]<<"\n";
	}
	fastaFile.close();
}

void generateSAM(string samFileName,string baseName, vector<string> &beginSeqs, vector<string> &endSeqs,string refSeqID,string refSeq)
{
	ofstream samFile;
	samFile.open(samFileName.c_str());
	samFile<<"@RG\tID:"<<refSeqID<<"\n";
	for(int i=0;i<beginSeqs.size();i++)
	{
		string revStr = beginSeqs[i];
		std::reverse(revStr.begin(),revStr.end());
		string fullSeq = revStr+refSeq;
		string cigarStr = SSTR(revStr.size()) + std::string("I") + SSTR(refSeq.size()) + std::string("M");
		string retStr = baseName +"Begin_" + SSTR(i) + std::string("\t") + "0" + std::string("\t") + refSeqID + std::string("\t") + "1" + std::string("\t255\t") + cigarStr + std::string("\t*\t0\t0\t") + fullSeq + std::string("\t*\tRG:Z:Unpaired reads assembled against ") + refSeqID + std::string(" \n");
		samFile<<retStr;
	}
	for(int i=0;i<endSeqs.size();i++)
	{
		string fullSeq = refSeq + endSeqs[i];
		string cigarStr = SSTR(refSeq.size()) + std::string("M") + SSTR(endSeqs[i].size()) + std::string("I");
		string retStr = baseName +"End_" + SSTR(i) + std::string("\t") + "0" + std::string("\t") + refSeqID + std::string("\t") + "1" + std::string("\t255\t") + cigarStr + std::string("\t*\t0\t0\t") + fullSeq + std::string("\t*\tRG:Z:Unpaired reads assembled against ") + refSeqID + std::string(" \n");
		samFile<<retStr;
	}
	samFile.close();
}

int hammingDist(string str1,string str2)
{
	string shortString,longString;
	if(str1.length()>str2.length())
	{
		shortString=str2;
		longString=str1;
	}
	else
	{
		shortString=str1;
		longString=str2;
	}
	int diffCtr=0;
	for(int i=0;i<shortString.length();i++)
	{
		if(shortString[i]!=longString[i])
			diffCtr++;
	}
	return diffCtr;
}
bool greatLength(const string& s1, const string& s2)
{
    return s1.length() > s2.length();
}
bool lessLength (const string& s1, const string& s2)
{
    return s1.length() < s2.length();
}
vector<string> myriadClasses(vector<string> seqs)
{
	vector<string> classes;
	sort(seqs.begin(),seqs.end(), greatLength);
	if(seqs.size()>0)
	{
		classes.push_back(seqs[0]);
		for(int i=1;i<seqs.size();i++)
		{
			int belongsToClass=0;
			for(int j=0;j<classes.size();j++)
			{
				int dist = hammingDist(seqs[i],classes[j]);
				if(dist <= MAX_DIFF_HAMMING_DIST)
					belongsToClass=1; 
			}
			if(belongsToClass==0)
				classes.push_back(seqs[i]);
		}
	}
	return classes;
}
string longestLeftSeq(vector<string> seqs)
{
	std::sort(seqs.begin(),seqs.end(),lessLength);
	string commonBeginSeq="";
	if(seqs.size()>0)
	{
		for(int j=0;j<seqs[0].size();j++)
		{
			char columNT=seqs[0][j];
			int flgMismatch=0;
			for(int k=1;k<seqs.size();k++)
			{
				if(seqs[k][j]!=columNT)
				{
					flgMismatch=1;
					break;
				}
			}
			if(flgMismatch==0)
				commonBeginSeq+=columNT;
			else
				break;

		}
	}
	return commonBeginSeq;
}
string longestRightSeq(vector<string> seqs)
{
	std::sort(seqs.begin(),seqs.end(),lessLength);
	string commonEndSeq="";
	if(seqs.size()>0)
	{
		for(int j=0;j<seqs[0].size();j++)
		{
			char columNT=seqs[0][j];
			int flgMismatch=0;
			for(int k=1;k<seqs.size();k++)
			{
				if(seqs[k][j]!=columNT)
				{
					flgMismatch=1;
					break;
				}
			}
			if(flgMismatch==0)
				commonEndSeq+=columNT;
			else
				break;
		}
	}
	return commonEndSeq;
}

map<string,string> expandTransposaseSeqs(map<string,string> &hitSeqs, map<string,string> &rawReads,string readFilePath, string blast_path, int iterCtr, string baseName)
{
	map<string,string> newHitSeqs;
	ofstream transHitFile;
    transHitFile.open ("transHits.fasta");

	for(std::map<string,string>::iterator iter = hitSeqs.begin(); iter != hitSeqs.end(); iter++) {
		transHitFile<<">"<<(*iter).first<<"\n"<<(*iter).second<<"\n";
    }
    transHitFile.close();

    string cmdStr = blast_path + "makeblastdb -in transHits.fasta -dbtype nucl -title transHits -out transHits";
    int ret = system(cmdStr.c_str());

    cmdStr = blast_path + "blastn -task megablast -db ./transHits -evalue 1e-01 -ungapped -num_threads 16 -query " + readFilePath + " -out transHitsBLASTResult.txt -outfmt \"10 qseqid sseqid qstart qend sstart send\" ";
    ret = system(cmdStr.c_str());
	
    string fname = "./transHitsBLASTResult.txt";
    ifstream transHitsBlastResult(fname.c_str());

	multimap<string,string> blastHits;
	string blast_str;
    while(getline(transHitsBlastResult,blast_str))
	{
		std::vector<std::string> tok = split(blast_str, ',');
		string hitSeqName =  tok[1];
		blastHits.insert(pair<string,string>(hitSeqName,blast_str));
	}
	transHitsBlastResult.close();
	int processedSeqCtr=0;
	for(std::map<string,string>::iterator iter = hitSeqs.begin(); iter != hitSeqs.end(); iter++) {
		string thisHitSeq = (*iter).second;
		string thisHitSeqName = (*iter).first;
		vector<string> beginMatches, endMatches;
		getRepeatSeqs(thisHitSeqName,thisHitSeq, rawReads, blastHits,beginMatches,endMatches);
		cout<<"\n\n\n\nTransposase Hit Sequence: "<<processedSeqCtr++<<"\n";
		cout<<"Begin Matches:"<<beginMatches.size()<<"\n";
		cout<<"End Matches:"<<endMatches.size()<<"\n";

		for(int j=0;j<beginMatches.size();j++)
		{
			string revStr = beginMatches[j];
			std::reverse(revStr.begin(),revStr.end());
			beginMatches[j] = revStr;
		}
		vector<string> uniqBeginMatches;
		for(int j=0;j<beginMatches.size();j++)
		{
			int found=0;
			for(int k=0;k<beginMatches.size();k++)
			{
				if(j!=k){
					if (beginMatches[k].find(beginMatches[j]) != std::string::npos) {
						found=1;
						break;
					}
				}
			}
			if(found==0)
				uniqBeginMatches.push_back(beginMatches[j]);
		}
		std::sort(uniqBeginMatches.begin(),uniqBeginMatches.end(),lessLength);
		//Get Longest Common Left Match
		cout<<"Longest Common Left Sequence:"<<'\n';
		string commonBeginSeq = longestLeftSeq(uniqBeginMatches);
		cout<<commonBeginSeq<<'\n';
		cout<<"\nLeft Match Options:"<<"\n";
		for(int j=0;j<uniqBeginMatches.size();j++)
		{
			int copyCtr=0;
			for(int k=0;k<beginMatches.size();k++)
			{
				if (uniqBeginMatches[j].find(beginMatches[k]) != std::string::npos &&
				     commonBeginSeq.find(beginMatches[k]) == std::string::npos) {
					copyCtr++;
				}
			}
			cout<<uniqBeginMatches[j]<<",ReadSupport:"<<copyCtr<<"\n";
		}
		vector<string> beginClasses = myriadClasses(uniqBeginMatches);
		cout<<"Class Count:"<<beginClasses.size()<<"\n";
		
		//The other side...
		vector<string> uniqEndMatches;
		for(int j=0;j<endMatches.size();j++)
		{
			int found=0;
			for(int k=0;k<endMatches.size();k++)
			{
				if(j!=k){
					if (endMatches[k].find(endMatches[j]) != std::string::npos) {
						found=1;
						break;
					}
				}
			}
			if(found==0)
				uniqEndMatches.push_back(endMatches[j]);
		}

		//Get Longest Common Right Match
		std::sort(uniqEndMatches.begin(),uniqEndMatches.end(),lessLength);
		cout<<"\n\nLongest Common Right Sequence:"<<'\n';
		string commonEndSeq = longestRightSeq(uniqEndMatches);
		cout<<commonEndSeq<<'\n';
		cout<<"\nRight Match Options:"<<"\n";
		for(int j=0;j<uniqEndMatches.size();j++)
		{
			int copyCtr=0;
			for(int k=0;k<endMatches.size();k++)
			{
				if (uniqEndMatches[j].find(endMatches[k]) != std::string::npos &&
				     commonEndSeq.find(endMatches[k]) == std::string::npos) {
					copyCtr++;
				}
			}
			cout<<uniqEndMatches[j]<<",ReadSupport:"<<copyCtr<<"\n";
		}
		vector<string> endClasses = myriadClasses(uniqEndMatches);
		cout<<"Class Count:"<<endClasses.size()<<"\n";
		
		if(uniqBeginMatches.size()>0 && uniqEndMatches.size()>0)
		{
			//Write fasta sequence for sam file
			string fastaFileName = baseName+"_"+SSTR(iterCtr)+"_"+thisHitSeqName+".fasta";
			ofstream fastaFile;
			fastaFile.open(fastaFileName.c_str());
			fastaFile<<">"<<thisHitSeqName<<"\n";
			fastaFile<<thisHitSeq<<"\n";
			fastaFile.close();
			//write SAM file
			generateSAM(baseName+"_"+SSTR(iterCtr)+"_"+ thisHitSeqName +".sam",thisHitSeqName, uniqBeginMatches, uniqEndMatches,thisHitSeqName,thisHitSeq);
		}
		int childSeqID=0;
		if(thisHitSeq.size() <= MAX_SEED_LENGTH)
		{
			for(int j=0;j<beginClasses.size();j++)
			{
				for(int k=0;k<endClasses.size();k++)
				{
					string revStr = beginClasses[j];
					std::reverse(revStr.begin(),revStr.end());
					if(beginClasses.size()<3 && endClasses.size()<3)
					{
						newHitSeqs.insert(pair<string,string>(thisHitSeqName + "_" + SSTR(childSeqID),revStr + thisHitSeq + endClasses[k]));
						childSeqID++;
					}
					else if(beginClasses.size()<3) 
					{
						newHitSeqs.insert(pair<string,string>(thisHitSeqName + "_" + SSTR(childSeqID),revStr + thisHitSeq));
						childSeqID++;
					}
					else if(endClasses.size()<3) 
					{
						newHitSeqs.insert(pair<string,string>(thisHitSeqName + "_" + SSTR(childSeqID),thisHitSeq + endClasses[k]));
						childSeqID++;
					}
				}
			}
		}
	}	
	cout<<"Working ISeqs:"<<newHitSeqs.size()<<"\n";
	return newHitSeqs;
}

int main(int argc, char **argv)
{
	//ISQuest
	//arg 1: Conf File
	//arg 2: Contigs
    //arg 3: Raw Reads
	//arg 4: Out Name
	string blast_str,conf_str,blast_path="",blast_db_path="";
	int maxLoopCtr=0;
    ifstream ctgRes(argv[2]);
	ifstream confFile(argv[1]);
	ifstream readFile(argv[3]);

	if(confFile.fail() || ctgRes.fail() || readFile.fail())
	{
		cerr<<"\nThe files could not be opened.";
		return -1;
	}

	while(getline(confFile,conf_str))
	{
	    if(conf_str.find("BLASTPATH") != std::string::npos)
	    {
            blast_path = conf_str.substr(conf_str.find("=")+1);
	    }
		if(conf_str.find("BLASTDB") != std::string::npos)
	    {
            blast_db_path = conf_str.substr(conf_str.find("=")+1);
	    }
		if(conf_str.find("LOOPCOUNTER") != std::string::npos)
	    {
            maxLoopCtr = atoi(conf_str.substr(conf_str.find("=")+1).c_str());
	    }
	}
	map<string,string> rawReads;
	readFasta(readFile,rawReads);
    
	//string cmdStr = blast_path + "blastx -db "+ blast_db_path +"nr -evalue 10 -max_target_seqs 500 -threshold 99 -window_size 4 -gapopen 12 -gapextend 2 -num_threads 16 -query " + string(argv[2]) + " -out ctgBLASTResult.txt -outfmt \"10 qseqid sseqid stitle salltitles qstart qend sstart send\"";
	string cmdStr = blast_path + "blastn -task megablast -db "+ blast_db_path +"nt -evalue 1e-01 -query " + string(argv[2]) + " -out ctgBLASTResult.txt -outfmt \"10 qseqid sseqid stitle salltitles qstart qend sstart send\"";
    int ret = system(cmdStr.c_str());
    
    string fname = "./ctgBLASTResult.txt";
    ifstream ctgBlastFile(fname.c_str());
	vector<string> isRelatedCtgBLASTResult;
	while(getline(ctgBlastFile,blast_str))
	{
		if (blast_str.find("transposase") != std::string::npos)
			isRelatedCtgBLASTResult.push_back(blast_str);
	}

    map<string,string> hitSeqs = getTransHitSeqs(isRelatedCtgBLASTResult,ctgRes);
	
	int loopCtr=0;
	do{
		hitSeqs = expandTransposaseSeqs(hitSeqs,rawReads, string(argv[3]), blast_path,loopCtr, string(argv[4]));
		loopCtr++;
	}while(hitSeqs.size()>0 && loopCtr < maxLoopCtr);
	
	readFile.close();
    ctgBlastFile.close();
	confFile.close();
	ctgRes.close();
	return 0;
}
