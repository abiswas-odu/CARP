#include<map>
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include <vector>

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

int readSAMFile(ifstream& multiFile, vector<string> &theVec)
{
	string str,id,seq; 
	int maxLen=0;
	while(getline(multiFile,str))
	{
		if (str[0] == '@')
			continue;
		std::vector<std::string> tok = split(str, '\t');
		std::size_t found = tok[5].find("I");
		
		if (found==std::string::npos)
			return maxLen; 
			
		//Incorporate P later
		std::size_t foundP = tok[5].find("P");
		if (foundP != std::string::npos)
			continue;
		
		int insertLen = atoi(tok[5].substr(0,found).c_str());
		string insertSeq = tok[9].substr(0,insertLen);
		theVec.push_back(insertSeq);
		if(insertLen>maxLen)
			maxLen=insertLen;
	}
	return maxLen;
}
int main(int argc, char **argv)
{
	//Program to find genes/intergenic regions 
	//arg 1: Gene Matched Results
	map<string,int> class_map;
	vector<string> samReads; 
			
	ifstream samFile(argv[1]);
	
	if(samFile.fail())
	{
		cerr<<"\nThe files could not be opened.";
		return -1;
	}
	
	int maxLen = readSAMFile(samFile, samReads);
	
	for(int i=0;i<samReads.size();i++)
	{
		string revStr = samReads[i];
		std::reverse(revStr.begin(),revStr.end());
		samReads[i] = revStr;
	}
	
	for(int i=1;i<=maxLen;i++)
	{
		for(int j=0;j<samReads.size();j++)
		{
			string readMer = samReads[j].substr(0,i);
			if(readMer.size()==i)
			{
				std::map<string,int>::iterator iter=class_map.find(readMer);
				if(iter != class_map.end())
				{
					int class_ctr = (*iter).second;
					class_ctr++;
					iter->second = class_ctr;
				}
				else
				{
					class_map.insert(pair<string,int>(readMer,1));
				}
			}
		}
		//Column Index and map size
		cout<<maxLen-i+1<<","<<class_map.size()<<",";
		vector<int> class_sizes;  
		vector<string> class_seq;  
		for(std::map<string,int>::iterator it=class_map.begin(); it!=class_map.end(); ++it)
		{
			int readsInClass = (*it).second;
			if (readsInClass > 0)
			{
				class_sizes.push_back(readsInClass);
				class_seq.push_back((*it).first);
			}
		}
		cout<<class_sizes.size()<<",";
		sort(class_sizes.begin(),class_sizes.end(),std::greater<int>());
		for(int k=0;k<class_sizes.size();k++)
		{
			cout<<class_sizes[k]<<",";
		}
		cout<<"\n";
		/*for(int k=0;k<class_sizes.size();k++)
		{
			cout<<"repeat_variant"<<k<<"\n";
			cout<<class_seq[k]<<"\n";
		}
		cout<<"\n";
		cout<<"\n";
		cout<<"\n";
		cout<<"\n";*/
		class_map.clear();
	}
	samFile.close();
	return 0;	
}	
