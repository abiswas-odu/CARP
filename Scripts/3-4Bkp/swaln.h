#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <sys/time.h>

#define GAP 2.0
#define MATCH 1.0
#define MISMATCH 2.0

using namespace std;

double similarity_score(char a,char b);
double find_array_max(double array[],int length, int &ind);
double swalign(string seq_a, string seq_b, int &seq_a_start, int &seq_a_end,int &seq_b_start,int &seq_b_end, int &hitLen);


