#ifndef MBBP_BGRAPH_H
#define MBBP_BGRAPH_H
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <map>


using namespace std;

#define SZ(x) ((int)x.size())

//ofstream __HEAD_H_FOUT;
#define ASSERT(v) {if (!(v)) {cerr<<"ASSERT FAIL @ "<<__FILE__<<":"<<__LINE__<<endl; exit(1);}}
#define INFO(...) do { ___debug( #__VA_ARGS__,  __VA_ARGS__,cout); cout<<endl; if(__HEAD_H_FOUT.is_open()){___debug( #__VA_ARGS__,  __VA_ARGS__,__HEAD_H_FOUT); __HEAD_H_FOUT<<endl;}  } while(0)

typedef unsigned int uint;
typedef unsigned long ulong;

using namespace std;
class bgraph{
public:
    uint des_count=0;
    uint novL;
    uint novR;
    uint noE;
    int max_degree_L{};
    int max_degree_R{};
    int offset;
    bool isSwitch=false;

    vector<vector<uint>> adjL;
    vector<uint> cohesionL; // todo
    vector<unordered_set<uint>> adjHL;
    vector<vector<uint>> adjR;
    vector<uint> cohesionR;  // todo
    vector<unordered_set<uint>> adjHR;

    string datasets;



    bgraph();
    void dispalyAdj();
    void sortAdj();
    static bool checkNumber(string &ss);
    static void Stringsplit(const string& str, const string& split, vector<string>& res);
    void readEdgeList(const string&);

    void switchLR();
    uint maxDegree(vector<vector<uint>>&);
    void calculate_number_of_edges();

    ~bgraph();

};


#endif //MBBP_BGRAPH_H

