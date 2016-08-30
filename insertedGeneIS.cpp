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
				if((gapPos[i][1]-gapPos[i][0])>20)
				{
					cout<<">"<<id<<":("<<gapPos[i][0]+1<<","<<gapPos[i][1]+1<<")\n"<<seq.substr(gapPos[i][0],(gapPos[i][1]-gapPos[i][0]))<<"\n";
				}
			}
		}
};

class GeneAlign
{
	private: 
		string id;
		string seq;
		vector<string> alignedReads;
		vector< vector<int> > readPos; 
		vector< vector<int> > genePos;
		vector< vector<int> > gapPos;
		
	public:
		void setId(string s) { id=s;}
		string getId() {return id;}
		void setSeq(string s) {seq=s;}
		string getSeq() {return seq;}
		int getReadCtr() {return alignedReads.size();}
		int getGapCtr() {return gapPos.size();}
		void addRead(string read,int qs,int qe,int ss,int se)
		{
			alignedReads.push_back(read);
			
			vector<int> qrow;
			qrow.push_back(qs-1);
			qrow.push_back(qe-1);
			genePos.push_back(qrow);
			
			vector<int> srow;
			srow.push_back(ss-1);
			srow.push_back(se-1);
			readPos.push_back(srow);
		}
		string getRead(string read)
		{	
			for(int i=0;i<alignedReads.size();i++)
			{
				if(alignedReads[i]==read)
				{
					string ret_str = read + std::string(",") + SSTR(genePos[i][0]) + std::string(",") + SSTR(genePos[i][1]) + std::string(",") + SSTR(readPos[i][0]) + std::string(",") + SSTR(readPos[i][1]);
					return ret_str; 
				}
			}
		}
		void calculateGaps()
		{
			vector<int> gr;
			for(int i=0;i<seq.size();i++)
				gr.push_back(1);
			for(int i=0;i<genePos.size();i++)
				for(int j=genePos[i][0];j<=genePos[i][1];j++)
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
		
		void printGenePos(map<string,UnusedReads> &read_list)
		{
			for(int i=0;i<genePos.size();i++)
			{
				map<string,UnusedReads>::iterator it;
				it = read_list.find(alignedReads[i]);
				cout<<genePos[i][0]+1<<","<<genePos[i][1]+1<<";"<<readPos[i][0]+1<<","<<readPos[i][1]+1<<";"<<(*it).second.getSeq().size()<<"\n";
			}
		}
		void printGaps(map<string,UnusedReads> &read_list)
		{
			for(int i=0;i<gapPos.size();i++)
			{
				cout<<"[("<<gapPos[i][0]+1<<","<<gapPos[i][1]+1<<");"<<(gapPos[i][1]-gapPos[i][0])<<"],";
			}
		}
		
		bool isInsertionGap(map<string,UnusedReads> &read_list)
		{
			for(int i=0;i<gapPos.size();i++)
			{
				int sf=0,ef=0;
				if(gapPos[i][0]>0 && gapPos[i][1]<seq.size() && (gapPos[i][1]-gapPos[i][0])<=GAPLEN)
				{
					vector<string> gapStartReads;
					vector<string> gapEndReads;
					vector< vector<int> > startReadPos; 
					vector< vector<int> > endReadPos; 
					for(int j=0;j<genePos.size();j++)
					{
						//hit nust be within 5bp of gap start
						if((genePos[j][1]+1)<=gapPos[i][0] && (genePos[j][1]+1)>=(gapPos[i][0]-5))
						{
							gapStartReads.push_back(alignedReads[j]);
							map<string,UnusedReads>::iterator it;
							it = read_list.find(alignedReads[j]);
							
							if(readPos[j][1]>=readPos[j][0])
							{
								vector<int> rrow;
								rrow.push_back(readPos[j][1]+1);
								rrow.push_back((*it).second.getSeq().size()-1);
								startReadPos.push_back(rrow);
							}
							else
							{
								vector<int> rrow;
								rrow.push_back(0);
								rrow.push_back(readPos[j][1]-1);
								startReadPos.push_back(rrow);
							}
							if(readPos[j][1]>=GAPOVERLAP && readPos[j][1]<=((*it).second.getSeq().size()-GAPOVERLAP))
								sf=1;
						}
					}
					for(int j=0;j<genePos.size();j++)
					{
						//hit nust be within 5bp of gap end
						if((genePos[j][0]-1)>=gapPos[i][1] && (genePos[j][0]-1)<=(gapPos[i][1]+5))
						{
							gapEndReads.push_back(alignedReads[j]);
							map<string,UnusedReads>::iterator it;
							it = read_list.find(alignedReads[j]);
							if(readPos[j][0]>=readPos[j][1])
							{
								vector<int> rrow;
								rrow.push_back(readPos[j][0]+1);
								rrow.push_back((*it).second.getSeq().size()-1);
								endReadPos.push_back(rrow);
							}
							else
							{
								vector<int> rrow;
								rrow.push_back(0);
								rrow.push_back(readPos[j][0]-1);
								endReadPos.push_back(rrow);
							}
							if(readPos[j][0]>=GAPOVERLAP && readPos[j][0]<=((*it).second.getSeq().size()-GAPOVERLAP))
								ef=1;
						}
					}
					if(sf==1 && ef==1)
					{
						//check different reads hit at both ends
						int returnFlag=1;
						for(int x=0;x<gapStartReads.size();x++)
							for(int y=0;y<gapEndReads.size();y++)
								if(gapStartReads[x]==gapEndReads[y])
									returnFlag=0;
						//if no insersection
						if(returnFlag)
						{
							cout<<">"<<id<<"[("<<gapPos[i][0]+1<<","<<gapPos[i][1]+1<<");"<<(gapPos[i][1]-gapPos[i][0])<<"]\n"<<seq<<"\n";
							filebuf fb_gen_alignment;
							fb_gen_alignment.open((id + std::string("_reads.fasta")).c_str(),ios::out);
							ostream os(&fb_gen_alignment);
							for(int j=0;j<gapStartReads.size();j++)
							{
								map<string,UnusedReads>::iterator it;
								it = read_list.find(gapStartReads[j]);
								os<<">"<<gapStartReads[j]<<"S\n"<<(*it).second.getSeq().substr(startReadPos[j][0],(startReadPos[j][1]-startReadPos[j][0])+1)<<"\n";
							}
							for(int j=0;j<gapEndReads.size();j++)
							{
								map<string,UnusedReads>::iterator it;
								it = read_list.find(gapEndReads[j]);
								os<<">"<<gapEndReads[j]<<"E\n"<<(*it).second.getSeq().substr(endReadPos[j][0],(endReadPos[j][1]-endReadPos[j][0])+1)<<"\n";
							}
							fb_gen_alignment.close();
							return true; 
						}
					}
				}
			}
			return false;
		}
		void printInsertionGaps(map<string,UnusedReads> &read_list,ostream &os)
		{
			for(int i=0;i<gapPos.size();i++)
			{
				int sf=0,ef=0;
				if(gapPos[i][0]>0 && gapPos[i][1]<seq.size() && (gapPos[i][1]-gapPos[i][0])<=GAPLEN)
				{
					for(int j=0;j<genePos.size();j++)
					{
						if(gapPos[i][0]==genePos[j][1]+1)
						{
							map<string,UnusedReads>::iterator it;
							it = read_list.find(alignedReads[j]);
							if(readPos[j][1]>=GAPOVERLAP && readPos[j][1]<=((*it).second.getSeq().size()-GAPOVERLAP))
								sf=1;
						}
						if(gapPos[i][1]==genePos[j][0]-1)
						{
							map<string,UnusedReads>::iterator it;
							it = read_list.find(alignedReads[j]);
							if(readPos[j][0]>=GAPOVERLAP && readPos[j][0]<=((*it).second.getSeq().size()-GAPOVERLAP))
								ef=1;
						}
					}
					if(sf==1 && ef==1)
					{
						//Print gene seqs in Fasta format
						//os<<">"<<id<<"[("<<gapPos[i][0]+1<<","<<gapPos[i][1]+1<<");"<<(gapPos[i][1]-gapPos[i][0])<<"]\n"<<seq<<"\n";
						//Print a csv list of interrupted genes
						//cout<<"@RG\tID:"<<id<<","<<gapPos[i][0]+1<<","<<gapPos[i][1]+1<<"\n";
						//os<<"@RG\tID:"<<id<<","<<gapPos[i][0]+1<<","<<gapPos[i][1]+1<<"\n";
						for(int j=0;j<alignedReads.size();j++)
						{
							map<string,UnusedReads>::iterator it;
							it = read_list.find(alignedReads[j]);
							//cout<<alignedReads[j]<<","<<id<<","<<gapPos[j][0]+1<<","<<gapPos[j][1]+1<<","<<readPos[j][0]+1<<","<<readPos[j][1]+1<<","<<alignedReads[j]<<"\n";
							//string anotherSAMLine=createSAMEntry(alignedReads[j],id,genePos[j][0]+1,genePos[j][1]+1,readPos[j][0]+1,readPos[j][1]+1,(*it).second.getSeq());
							//os<<anotherSAMLine;
							os<<">"<<alignedReads[j]<<"\n"<<(*it).second.getSeq()<<"\n";
						}
						return;
					}
				}
			}
		}
		string createSAMEntry(string read, string ref, int refMatchStart, int refMatchEnd, int readMatchStart, int readMatchEnd, string readSeq)
		{
			string f_r="0";
			string cigarStr="";
			if(readMatchStart>readMatchEnd)
			{
				f_r="16";
				int startMisMatch = readSeq.size() - readMatchStart;
				if(startMisMatch>0)
					cigarStr= cigarStr + convertInt(startMisMatch) + std::string("I");
				
				int matchLen = readMatchStart-readMatchEnd+1;
				if(matchLen>0)
					cigarStr= cigarStr + convertInt(matchLen) + std::string("M");
				
				int endMisMatch = readMatchEnd - 1;
				if(endMisMatch>0)
					cigarStr= cigarStr + convertInt(endMisMatch) + std::string("I");
			}
			else
			{
				int startMisMatch = readMatchStart-1;
				if(startMisMatch>0)
					cigarStr= cigarStr + convertInt(startMisMatch) + std::string("I");
				
				int matchLen = readMatchEnd-readMatchStart+1;
				if(matchLen>0)
					cigarStr= cigarStr + convertInt(matchLen) + std::string("M");
				
				int endMisMatch = readSeq.size() - readMatchEnd;
				if(endMisMatch>0)
					cigarStr=  cigarStr + convertInt(endMisMatch) + std::string("I");
			}
			string retStr = read + std::string("\t") + f_r + std::string("\t") + ref + std::string("\t") + convertInt(refMatchStart) + std::string("\t255\t") + cigarStr + std::string("\t*\t0\t0\t") + readSeq + std::string("\t*\tRG:Z:Unpaired reads assembled against ") + id + std::string(" \n");
			return retStr;

		}
		string convertInt(int number)
		{
		   stringstream ss;//create a stringstream
		   ss << number;//add number to the stream
		   return ss.str();//return a string with the contents of the stream
		}
};

class ISInsertion
{
	private: 
		string gene;
		int ins_start;
		int ins_end;
		string read;
		int isMatch_start;
		int isMatch_end;
		int readMatch_start;
		int readMatch_end;
	public:
		void setInsertion(string g, int inss, int inse, string r, int isMs, int isMe, int rMs,int rMe)
		{
			gene=g;
			ins_start=inss;
			ins_end=inse;
			read=r;
			isMatch_start=isMs;
			isMatch_end=isMe;
			readMatch_start=rMs;
			readMatch_end=rMe;
		}
		string getG()
		{ return gene;}
		string getR()
		{ return read;}
		int getInss()
		{ return ins_start;}
		int getInse()
		{ return ins_end;}
		int getIsMs()
		{ return isMatch_start;}
		int getIsMe()
		{ return isMatch_end;}
		int getRMs()
		{ return readMatch_start;}
		int getRMe()
		{ return readMatch_end;}	

		int getGeneStart() const
		{
			string retStr = gene.substr(gene.find(":")+1,gene.find("-")-gene.find(":")+1);
			if(retStr[0]=='c')
				return atoi(retStr.substr(1).c_str());
			else
				return atoi(retStr.c_str());
		}
		
		int getGeneEnd() const
		{
			string retStr = gene.substr(gene.find("-")+1);
			return atoi(retStr.c_str());
		}
		
		int getGeneLength() const
		{
			return getGeneEnd()-getGeneStart()+1;
		}
		
		bool operator<( const ISInsertion& val ) const { 
			if(getGeneStart() == val.getGeneStart())
				return ins_start < val.ins_start;
			else
				return getGeneStart() < val.getGeneStart(); 
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
				theMap.insert (pair<string,string>(id,temp) );
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

void readGeneFasta(ifstream& multiFile, vector<GeneAlign> &theVector)
{
	string str,id,seq; 
	getline(multiFile,str);
	id=str.substr(1,str.find(" ")-1); 
	while(getline(multiFile,str))
	{
		int found = str.find(">");
		if (found!=string::npos)
		{
			GeneAlign gen;
			gen.setId(id);
			gen.setSeq(seq);
			theVector.push_back(gen);
			seq="";
			id = str.substr(1,str.find(" ")-1);
		}
		else
		{
			seq=seq+str;
		}
	}
	GeneAlign gen;
	gen.setId(id);
	gen.setSeq(seq);
	theVector.push_back(gen);
}

int main(int argc, char **argv)
{
	//Program to find genes/intergenic regions 
	//arg 1: Gene Sequences
	//arg 2: Reads
	//arg 3: Gene BLAST result
	//arg 4: IS BLAST result
	//arg 5: GenBank CDS Records
	
	string blast_str; 
	vector<GeneAlign> gene_list;
	map<string,UnusedReads> read_list;
	multimap<string,string> gene_read_map;
	map<string,string> genBank_rec;
	int i=0;
	
	ifstream genes(argv[1]);
	ifstream reads(argv[2]);
	ifstream geneBlastRes(argv[3]);
	ifstream isBlastRes(argv[4]);
	ifstream genBank(argv[5]);
	
	if(genes.fail() || reads.fail() || geneBlastRes.fail() || isBlastRes.fail() || genBank.fail())
	{
		cerr<<"\nThe files could not be opened.";
		return -1;
	}
	
	readGeneFasta(genes,gene_list);
	readFasta(reads,read_list);
	readGenBankCDS(genBank,genBank_rec);
	
	while(getline(geneBlastRes,blast_str))
	{
		std::vector<std::string> tok = split(blast_str, ',');
		if(tok[0] == gene_list[i].getId())
		{
			gene_list[i].addRead(tok[1],atoi(tok[2].c_str()),atoi(tok[3].c_str()),atoi(tok[4].c_str()),atoi(tok[5].c_str()));
			gene_read_map.insert(pair<string,string>(tok[1],gene_list[i].getId()));
			
		}
		else
		{
			i=0;
			while(tok[0] != gene_list[i].getId())
				i++;
			gene_list[i].addRead(tok[1],atoi(tok[2].c_str()),atoi(tok[3].c_str()),atoi(tok[4].c_str()),atoi(tok[5].c_str()));
			gene_read_map.insert(pair<string,string>(tok[1],gene_list[i].getId()));
		}
	}
	
	vector<ISInsertion> isIns_list;
	while(getline(isBlastRes,blast_str))
	{
		std::vector<std::string> tok = split(blast_str, ',');
		map<string,UnusedReads>::iterator it;
		it = read_list.find(tok[0]);
		if(atoi(tok[2].c_str())>50 || atoi(tok[3].c_str())<((*it).second.getSeq().size()-50) )
		{
			pair<multimap<string, string>::iterator, multimap<string, string>::iterator> ppp;
			ppp = gene_read_map.equal_range(tok[0]);
			for (multimap<string,string>::iterator it = ppp.first; it != ppp.second; ++it)
			{
				for(int i=0;i<gene_list.size();i++)
				{
					if(gene_list[i].getId() == (*it).second)
					{
						ISInsertion isIns;
						std::vector<std::string> isCoords = split(gene_list[i].getRead(tok[0]),',');
						isIns.setInsertion((*it).second,atoi(isCoords[1].c_str()),atoi(isCoords[2].c_str()),isCoords[0],atoi(tok[2].c_str()),atoi(tok[3].c_str()),atoi(isCoords[3].c_str()),atoi(isCoords[4].c_str()));
						isIns_list.push_back(isIns);
						//cout<<(*it).second<<","<<tok[2]<<","<<tok[3]<<","<<gene_list[i].getRead(tok[0])<<"\n";
						break;
					}
				}
			}
		}
	}
	genes.close();
	reads.close();
	geneBlastRes.close();
	isBlastRes.close();
	std::sort(isIns_list.begin(),isIns_list.end());
	//Header
	cout<<"CDS Start,CDS End,Length,Insertion Start,Insertion End,Type,Read,IS Match Start,IS Match End,Read Match Start,Read Match End,Unique,Locus,Product,Gene\n";
	int covCounter=0, igCounter=0, endsCounter=0;
	for(int i=0;i<isIns_list.size();i++)
	{
		string igClass;
		covCounter++;
		if((abs(isIns_list[i].getGeneLength()-1) == isIns_list[i].getInse()) || isIns_list[i].getInss()==0)
		{
			igClass="Intergenic";
			igCounter++;
		}
		else if(((abs(isIns_list[i].getGeneLength()-1) - isIns_list[i].getInse()) <=3) || isIns_list[i].getInss()<=3)
		{
			igClass="Ends";
			endsCounter++;
		}
		else
			igClass = SSTR(abs(isIns_list[i].getGeneLength()) - isIns_list[i].getInse()) + ";" + SSTR(isIns_list[i].getInss());
			
		cout<<isIns_list[i].getGeneStart()<<","<<isIns_list[i].getGeneEnd()<<","<<isIns_list[i].getGeneLength()-1<<","<<isIns_list[i].getInss()<<","<<isIns_list[i].getInse()<<","<<igClass<<","<<isIns_list[i].getR()<<","<<isIns_list[i].getIsMs()<<","<<isIns_list[i].getIsMe()<<","<<isIns_list[i].getRMs()<<","<<isIns_list[i].getRMe();
		
		if(i==isIns_list.size()-1)
		{
			string gene_class="";
			if(covCounter<4)
				gene_class="Low Coverage";
			else
			{	
				int instCounter = covCounter - igCounter - endsCounter;
				gene_class=igCounter <= instCounter?(endsCounter<=instCounter?"Good":"Ends"):(endsCounter<=igCounter?"IG":"Ends");
			}
		
			string cdsKey="";
			if(isIns_list[i].getGeneStart() > isIns_list[i].getGeneEnd())
				cdsKey = SSTR(isIns_list[i].getGeneEnd()) + ".." + SSTR(isIns_list[i].getGeneStart());
			else
				cdsKey = SSTR(isIns_list[i].getGeneStart()) + ".." + SSTR(isIns_list[i].getGeneEnd());
			map<string,string>::iterator it;
			it = genBank_rec.find(cdsKey);
			cout<<","<<isIns_list[i].getGeneEnd()<<","<<(*it).second<<","<<gene_class<<"\n";
			endsCounter=igCounter=covCounter=0;
			//cout<<","<<isIns_list[i].getGeneEnd()<<","<<cdsKey<<"\n";
		}
		else
		{
			if(isIns_list[i].getG() == isIns_list[i+1].getG())
				cout<<",0"<<"\n";
			else
			{
			
				string gene_class="";
				if(covCounter<4)
				gene_class="Low Coverage";
				else
				{	
					int instCounter = covCounter - igCounter - endsCounter;
					gene_class=igCounter <= instCounter?(endsCounter<=instCounter?"Good":"Ends"):(endsCounter<=igCounter?"IG":"Ends");
				}
			
				string cdsKey="";
				if(isIns_list[i].getGeneStart() > isIns_list[i].getGeneEnd())
					cdsKey = SSTR(isIns_list[i].getGeneEnd()) + ".." + SSTR(isIns_list[i].getGeneStart());
				else
					cdsKey = SSTR(isIns_list[i].getGeneStart()) + ".." + SSTR(isIns_list[i].getGeneEnd());
				map<string,string>::iterator it;
				it = genBank_rec.find(cdsKey);
				cout<<","<<isIns_list[i].getGeneEnd()<<","<<(*it).second<<","<<gene_class<<"\n";
				//cout<<","<<isIns_list[i].getGeneEnd()<<","<<cdsKey<<"\n";
				endsCounter=igCounter=covCounter=0;
			}
		}
	}
	return 0;	
}	
