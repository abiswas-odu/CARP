#include<map>
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<vector>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<cstring>
#include<algorithm>
#include<iomanip>
#include<omp.h>

#define MAX_TRANSHIT_LEN 100
#define MAX_SEED_LENGTH 5000
#define MAX_DIFF_HAMMING_DIST 10
#define MIN_CLASS_SIZE 4
#define MIN_IR_SCORE 13
#define MAX_IR_LEN 40
#define MIN_IR_LEN 20
#define MIN_GLOBAL_MATCHES_THRESH 0.950

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

class IRepeat
{
	private:
	double score;
	int irlstart,irlend,irrstart,irrend;
	string seqName,seq;
	public:
	void setScore(double s);
	void setIrlstart(int s);
	void setIrlend(int s);
	void setIrrstart(int s);
	void setIrrend(int s);
	void setSeqName(string name);
	void setSeq(string s);
	double getScore();
	int getIrlstart();
	int getIrlend();
	int getIrrstart();
	int getIrrend();
	string getSeq();
	string getSeqName();
};
