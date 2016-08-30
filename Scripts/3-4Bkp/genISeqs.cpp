#include "genISeqs.h"
#include "swaln.h"
extern "C" {
	#include "gbfp.h"
}

void IRepeat::setScore(double s){
	score=s;}
void IRepeat::setIrlstart(int s){
	irlstart=s;}
void IRepeat::setIrlend(int s){
	irlend=s;}
void IRepeat::setIrrstart(int s){
	irrstart=s;}
void IRepeat::setIrrend(int s){
	irrend=s;}
void IRepeat::setSeqName(string name){
	seqName=name;}
void IRepeat::setSeq(string s){
	seq=s;}
double IRepeat::getScore(){
	return score;}
int IRepeat::getIrlstart(){
	return irlstart;}
int IRepeat::getIrlend(){
	return irlend;}
int IRepeat::getIrrstart(){
	return irrstart;}
int IRepeat::getIrrend(){
	return irrend;}
string IRepeat::getSeq(){
	return seq;}
string IRepeat::getSeqName(){
	return seqName;}

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
const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%d-%m-%Y", &tstruct);
    return buf;
}
void create_directory(string dirName)
{
	string cmdStr = "mkdir " + dirName;
    int ret = system(cmdStr.c_str());
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

void generateGKBFile(string fname, string seqName,string seq, IRepeat &ir)
{
	ofstream gkbFile;
	std::transform(seq.begin(), seq.end(), seq.begin(), ::tolower);
	gkbFile.open(fname.c_str());
	gkbFile<<"LOCUS       "<<seqName<<" "<<seq.size()<<" bp "<<"DNA linear UNA "<<currentDateTime()<<"\n";
	gkbFile<<"DEFINITION  .\n";
	gkbFile<<"ACCESSION   "<<seqName<<"\n";
	gkbFile<<"VERSION     "<<seqName<<"\n";
	gkbFile<<"KEYWORDS    .\n";
	gkbFile<<"SOURCE      \n";
	gkbFile<<"  ORGANISM  \n";
	gkbFile<<"            .\n";
	gkbFile<<"FEATURES             Location/Qualifiers\n";
	gkbFile<<"     repeat_unit     "<<ir.getIrlstart()<<".."<<ir.getIrlend()<<"\n";
	gkbFile<<"                     /label=\"IRL\"\n";
	gkbFile<<"     repeat_unit     "<<ir.getIrrstart()<<".."<<ir.getIrrend()<<"\n";
	gkbFile<<"                     /label=\"IRR\"\n";	
	gkbFile<<"ORIGIN      \n";
	int isLineCtr=1, usedStr=0;
	for(int i=0;i<seq.size();i+=60)
	{
		gkbFile<<setw(9);
		gkbFile<<isLineCtr;
	
		isLineCtr+=60;
		if(usedStr<seq.size())
			gkbFile<<" "<<seq.substr(usedStr,10);
			
		usedStr+=10;
		if(usedStr<seq.size())
			gkbFile<<" "<<seq.substr(usedStr,10);
			
		usedStr+=10;
		if(usedStr<seq.size())
			gkbFile<<" "<<seq.substr(usedStr,10);
			
		usedStr+=10;
		if(usedStr<seq.size())
			gkbFile<<" "<<seq.substr(usedStr,10);
			
		usedStr+=10;
		if(usedStr<seq.size())
			gkbFile<<" "<<seq.substr(usedStr,10);
			
		usedStr+=10;
		if(usedStr<seq.size())
			gkbFile<<" "<<seq.substr(usedStr,10);
			
		usedStr+=10;
		gkbFile<<"\n";
	}
	gkbFile<<"//\n";
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
bool isInvertedRepeatFoundBLAST(string seq, string seqName, string blast_path, IRepeat &ir)
{
	if(seq.size()>1000)
	{
		string upStream = seq.substr(0,500);
		string downStream = seq.substr(seq.size()-500);
		
		ofstream irHitFile;
		irHitFile.open("upstream.fasta");
		irHitFile<<">upstream"<<"\n"<<upStream<<"\n";
		irHitFile.close();
	
		string cmdStr = blast_path + "makeblastdb -in upstream.fasta -dbtype nucl -title upstream -out upstream";
		int ret = system(cmdStr.c_str());
		
		irHitFile.open("downstream.fasta");
		irHitFile<<">downstream"<<"\n"<<downStream<<"\n";
		irHitFile.close();

		cmdStr = blast_path + "blastn -task blastn -db ./upstream -query downstream.fasta -gapopen 4 -gapextend 4 -penalty -3 -reward 1 -out irHitResult.txt -outfmt \"10 qseqid sseqid score qstart qend sstart send\" ";
		ret = system(cmdStr.c_str());
		
		ifstream readResultFile("irHitResult.txt");
		if(readResultFile.fail())
		{
			cerr<<"\nThe irHitResult.txt file could not be opened.";
			return false;
		}
		string blast_str;
		/*Get best IR Hit*/
		int max_score = MIN_IR_SCORE, max_qstart = 0, max_qend = 0,max_sstart = 0, max_send = 0, isHitFound = 0;
		int max_transLen = 0;
		while(getline(readResultFile,blast_str))
		{
			std::vector<std::string> tok = split(blast_str, ',');
			int score = atoi(tok[2].c_str());
			int qstart = atoi(tok[3].c_str());
			int qend = atoi(tok[4].c_str());
			int sstart = atoi(tok[5].c_str());
			int send = atoi(tok[6].c_str());
			int transLen = seq.size()-500 + send - qend;
			if((qend-qstart) >= MIN_IR_LEN && (qend-qstart) <= MAX_IR_LEN && score >= max_score && transLen > max_transLen && (sstart>send))
			{	
				isHitFound=1; 
				max_score = score;
				max_qstart = qstart;
				max_qend = qend;
				max_sstart = sstart;
				max_send = send;
				max_transLen = seq.size()-500 + send - qend;
			}
		} 
		if(isHitFound)
		{
			cout<<seqName<<"\n"<<seq<<"\n";
			cout<<"Max Score:"<<max_score<<",Downstream:("<<max_qstart<<","<<max_qend<<"),Upstream("<<seq.size()-500+max_sstart<<","<<seq.size()-500+max_send<<")\n";
			ir.setScore(max_score);
			ir.setSeqName(seqName);
			ir.setSeq(seq);
			ir.setIrlstart(max_qstart);
			ir.setIrlend(max_qend);
			ir.setIrrstart(seq.size()-500+max_send);
			ir.setIrrend(seq.size()-500+max_sstart);
			return true;
		}
	}
	return false;
}
bool isInvertedRepeatFoundSW(string seq, string seqName, IRepeat &ir)
{
	if(seq.size()>1200)
	{
		string upStream = seq.substr(0,500);
		string downStream = seq.substr(seq.size()-500);
		int max_qstart = 0, max_qend = 0,max_sstart = 0, max_send = 0, irLen=0;
		downStream=getRevComp(downStream);
		double max_score = swalign(upStream,downStream,max_qstart,max_qend,max_sstart,max_send,irLen);
		if(irLen >= MIN_IR_LEN && irLen <= MAX_IR_LEN && max_score>=MIN_IR_SCORE)
		{
			cout<<seqName<<"\n"<<seq<<"\n";
			cout<<"Max Score:"<<max_score<<",Downstream:("<<max_qstart<<","<<max_qend<<"),Upstream("<<seq.size()-max_sstart<<","<<seq.size()-max_send<<")\n";
			ir.setScore(max_score);
			ir.setSeqName(seqName);
			ir.setSeq(seq);
			ir.setIrlstart(max_qstart);
			ir.setIrlend(max_qend);
			ir.setIrrstart(seq.size()-max_send);
			ir.setIrrend(seq.size()-max_sstart);
			return true;
		}
	}
	return false;
}
map<string,string> expandTransposaseSeqs(map<string,string> &hitSeqs, map<string,string> &rawReads,string readFilePath, string blast_path, int iterCtr, string baseName,string writeIntermideateFiles)
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
		//cout<<"Longest Common Left Sequence:"<<'\n';
		string commonBeginSeq = longestLeftSeq(uniqBeginMatches);
		//cout<<commonBeginSeq<<'\n';
		//cout<<"\nLeft Match Options:"<<"\n";
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
			//cout<<uniqBeginMatches[j]<<",ReadSupport:"<<copyCtr<<"\n";
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
		//cout<<"\n\nLongest Common Right Sequence:"<<'\n';
		string commonEndSeq = longestRightSeq(uniqEndMatches);
		//cout<<commonEndSeq<<'\n';
		//cout<<"\nRight Match Options:"<<"\n";
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
			//cout<<uniqEndMatches[j]<<",ReadSupport:"<<copyCtr<<"\n";
		}
		vector<string> endClasses = myriadClasses(uniqEndMatches);
		cout<<"Class Count:"<<endClasses.size()<<"\n";
		IRepeat ir;
		bool foundIR = isInvertedRepeatFoundSW(thisHitSeq,thisHitSeqName, ir);
		if(foundIR)
		{
			//Write fasta sequence for sam file
			cout<<"Writing Result File\n";
			string fastaFileName = "FinalResult/" + baseName+"_"+SSTR(iterCtr)+"_"+thisHitSeqName+".fasta";
			ofstream fastaFile;
			fastaFile.open(fastaFileName.c_str());
			fastaFile<<">"<<thisHitSeqName<<"\n";
			fastaFile<<thisHitSeq<<"\n";
			fastaFile.close();
			//write SAM file
			generateSAM("FinalResult/" + baseName+"_"+SSTR(iterCtr)+"_"+ thisHitSeqName +".sam",thisHitSeqName, uniqBeginMatches, uniqEndMatches,thisHitSeqName,thisHitSeq);
			generateGKBFile("FinalResult/" + baseName+"_"+SSTR(iterCtr)+"_"+ thisHitSeqName +".gkb",thisHitSeqName, thisHitSeq, ir);
		}
		if(uniqBeginMatches.size()>0 && uniqEndMatches.size()>0 && writeIntermideateFiles == "YES")
		{
			//Write fasta sequence for sam file
			string fastaFileName = "level"+SSTR(iterCtr) + "/" + baseName+"_"+SSTR(iterCtr)+"_"+thisHitSeqName+".fasta";
			ofstream fastaFile;
			fastaFile.open(fastaFileName.c_str());
			fastaFile<<">"<<thisHitSeqName<<"\n";
			fastaFile<<thisHitSeq<<"\n";
			fastaFile.close();
			//write SAM file
			generateSAM("level"+SSTR(iterCtr) + "/" + baseName+"_"+SSTR(iterCtr)+"_"+ thisHitSeqName +".sam",thisHitSeqName, uniqBeginMatches, uniqEndMatches,thisHitSeqName,thisHitSeq);

		}
		int childSeqID=0;
		if(thisHitSeq.size() <= MAX_SEED_LENGTH && !foundIR)
		{
			for(int j=0;j<beginClasses.size();j++)
			{
				for(int k=0;k<endClasses.size();k++)
				{
					string revStr = beginClasses[j];
					std::reverse(revStr.begin(),revStr.end());
					if(beginClasses.size()<3 && endClasses.size()<3)
					{
						newHitSeqs.insert(pair<string,string>(revStr + thisHitSeq + endClasses[k],thisHitSeqName + "_" + SSTR(childSeqID)));
						childSeqID++;
					}
					else if(beginClasses.size()<3) 
					{
						newHitSeqs.insert(pair<string,string>(revStr + thisHitSeq,thisHitSeqName + "_" + SSTR(childSeqID)));
						childSeqID++;
					}
					else if(endClasses.size()<3) 
					{
						newHitSeqs.insert(pair<string,string>(thisHitSeq + endClasses[k],thisHitSeqName + "_" + SSTR(childSeqID)));
						childSeqID++;
					}
				}
			}
		}
	}	
	cout<<"Working ISeqs:"<<newHitSeqs.size()<<"\n";
	map<string,string> newReturnHitSeqs;
	for(std::map<string,string>::iterator iter = newHitSeqs.begin(); iter != newHitSeqs.end(); iter++) {
		string thisHitSeq = (*iter).first;
		string thisHitSeqName = (*iter).second;
		newReturnHitSeqs.insert(pair<string,string>(thisHitSeqName,thisHitSeq));
	}
	return newReturnHitSeqs;
}
static char *searchQualValue(const char *sWord, gb_feature *ptFeature) {
    gb_qualifier *i;
    for (i = ptFeature->ptQualifier; (i - ptFeature->ptQualifier) < ptFeature->iQualifierNum; i++) 
        if (strstr(i->sValue,sWord) != NULL)
            return i->sValue;
    return NULL;
}
bool isAnnotationHit(vector<string> sWord, gb_feature *ptFeature)
{
	char *sValue = NULL;
	for(int i=0;i<sWord.size();i++)
	{
		if ((sValue = searchQualValue(sWord[i].c_str(), ptFeature)) != NULL) {
			return true;
		}
	}
	return false;
}
bool isTransposaseHit(string blast_str,string keywords,string genbank_path, string wget_path)
{
	if (blast_str.find("transposase") != std::string::npos)
	{
		return true;
	}
	/* Some accession files cause errors. Black list them*/
	vector<string> blackListAcc;
	blackListAcc.push_back("CP005080");
	
	std::vector<std::string> tok = split(blast_str, ',');
	std::vector<std::string> searchKeys = split(keywords, ',');
	int endIndx = tok.size();
	string accCode =  tok[1];
	int sstart = atoi(tok[endIndx-2].c_str());
	int send = atoi(tok[endIndx-1].c_str());
    string accFilePath = genbank_path + accCode + ".gbk";
	ifstream gbkFile(accFilePath.c_str());
	
	if(std::find(blackListAcc.begin(), blackListAcc.end(), accCode)==blackListAcc.end())
	{
		if(gbkFile.fail())
		{
			cout<<"\nThe genebank files could not be found:"<<accFilePath<<"\n";
			cout<<"Downloading..."<<"\n";
			string urlStr = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&id=" + accCode;
			string cmdStr = wget_path + "wget \"" + urlStr +"\" -O " + accFilePath;
			int ret = system(cmdStr.c_str());
			if(ret==0)
			{
				char * cstr = new char [accFilePath.length()+1];
				std::strcpy(cstr, accFilePath.c_str());
				gb_data **pptSeqData, *ptSeqData;
				gb_feature *ptFeature;
				cout<<"\nThe genebank file downloaded:"<<accFilePath<<"\n";
				pptSeqData = parseGBFF(cstr);
				for (int i = 0; (ptSeqData = *(pptSeqData + i)) != NULL; i++) {
					cout<<"File ptr:"<<ptSeqData->iFeatureNum<<"\n";
					for (int j = 0; j < ptSeqData->iFeatureNum; j++) {
						ptFeature = (ptSeqData->ptFeatures + j);
						if(((sstart>=ptFeature->lStart && sstart<=ptFeature->lEnd) || 
							(send>=ptFeature->lStart && send<=ptFeature->lEnd))	&&
								isAnnotationHit(searchKeys,ptFeature))
						{
							cout<<blast_str<<"\n"<<ptFeature->sFeature<<","<<ptFeature->lStart<<","<<ptFeature->cDirection<<","<<ptFeature->lEnd<<","<<"\n";
							freeGBData(pptSeqData);
							delete[] cstr;
							return true;
						}
					}
				}
				freeGBData(pptSeqData);
				delete[] cstr;
			}
			else
				cout<<"\nThe genebank file download error:"<<accCode<<"\n";
		}
		else
		{
			char * cstr = new char [accFilePath.length()+1];
			std::strcpy(cstr, accFilePath.c_str());
			gb_data **pptSeqData, *ptSeqData;
			gb_feature *ptFeature;
			cout<<"\nThe genebank file found:"<<accFilePath<<"\n";
			pptSeqData = parseGBFF(cstr);
			for (int i = 0; (ptSeqData = *(pptSeqData + i)) != NULL; i++) {
				cout<<"File ptr:"<<ptSeqData->iFeatureNum<<"\n";
				for (int j = 0; j < ptSeqData->iFeatureNum; j++) {
					ptFeature = (ptSeqData->ptFeatures + j);
					if(((sstart>=ptFeature->lStart && sstart<=ptFeature->lEnd) || 
						(send>=ptFeature->lStart && send<=ptFeature->lEnd))	&&
							isAnnotationHit(searchKeys,ptFeature))
					{
						cout<<blast_str<<"\n"<<ptFeature->sFeature<<","<<ptFeature->lStart<<","<<ptFeature->cDirection<<","<<ptFeature->lEnd<<","<<"\n";
						freeGBData(pptSeqData);
						delete[] cstr;
						return true;
					}
				}
			}
			freeGBData(pptSeqData);
			delete[] cstr;
		}
	}
	else
        cout<<"\nThe genebank file blacklisted:"<<accFilePath<<"\n";
	return false;
}

int main(int argc, char **argv)
{
	//ISQuest
	//arg 1: Conf File
	//arg 2: Contigs
    //arg 3: Raw Reads
	//arg 4: Out Name
	string blast_str,conf_str,blast_path="",blast_db_path="",genbank_path="",wget_path="",writeIntermideateFiles="NO";
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
		if(conf_str.find("GENBANKPATH") != std::string::npos)
	    {
            genbank_path = conf_str.substr(conf_str.find("=")+1);
	    }
		if(conf_str.find("WGETPATH") != std::string::npos)
	    {
            wget_path = conf_str.substr(conf_str.find("=")+1);
	    }
		if(conf_str.find("WRITE_INTERMIDEATE_FILES") != std::string::npos)
	    {
            writeIntermideateFiles = conf_str.substr(conf_str.find("=")+1);
	    }
	}
	map<string,string> rawReads;
	readFasta(readFile,rawReads);
    
	//string cmdStr = blast_path + "blastx -db "+ blast_db_path +"nr -evalue 10 -max_target_seqs 500 -threshold 99 -window_size 4 -gapopen 12 -gapextend 2 -num_threads 16 -query " + string(argv[2]) + " -out ctgBLASTResult.txt -outfmt \"10 qseqid sseqid stitle salltitles qstart qend sstart send\"";
	//string cmdStr = blast_path + "blastn -task megablast -db "+ blast_db_path +"nt -evalue 1e-01 -query " + string(argv[2]) + " -out ctgBLASTResult.txt -outfmt \"10 qseqid sacc stitle salltitles qstart qend sstart send\"";
    //int ret = system(cmdStr.c_str());
    
    string fname = "./ctgBLASTResult.txt";
    ifstream ctgBlastFile(fname.c_str());
	vector<string> isRelatedCtgBLASTResult;
	string annotationSearchKeywords="transposase,intergrase,insertion";
	while(getline(ctgBlastFile,blast_str))
	{
		if(isTransposaseHit(blast_str,annotationSearchKeywords,genbank_path,wget_path))
			isRelatedCtgBLASTResult.push_back(blast_str);
	}

    map<string,string> transHitSeqs = getTransHitSeqs(isRelatedCtgBLASTResult,ctgRes);
	map<string,string> hitSeqs;
	for(std::map<string,string>::iterator iter = transHitSeqs.begin(); iter != transHitSeqs.end(); iter++) {
		string thisHitSeqName = (*iter).first;
		string thisHitSeq = (*iter).second;
		bool isSubString = false;
		for(std::map<string,string>::iterator iter2 = transHitSeqs.begin(); iter2 != transHitSeqs.end(); iter2++) {
			if((*iter2).second.find(thisHitSeq) != std::string::npos && thisHitSeqName != (*iter2).first) {
				isSubString=true;
				break;
			}
		}
		if(!isSubString)
			hitSeqs.insert(pair<string,string>(thisHitSeqName,thisHitSeq));
	}
	int loopCtr=0;
	create_directory("FinalResult");
	do{
		if(writeIntermideateFiles == "YES")
			create_directory("level"+SSTR(loopCtr)); 
		hitSeqs = expandTransposaseSeqs(hitSeqs,rawReads, string(argv[3]), blast_path,loopCtr, string(argv[4]),writeIntermideateFiles);
		loopCtr++;
	}while(hitSeqs.size()>0 && loopCtr < maxLoopCtr);
	
	readFile.close();
    ctgBlastFile.close();
	confFile.close();
	ctgRes.close();
	return 0;
}
