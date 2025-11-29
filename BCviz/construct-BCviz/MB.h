#ifndef RELEASEBCVIZ_MB_H
#define RELEASEBCVIZ_MB_H

#include <utility>
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <map>
#include <chrono>

#define ROW_SIZE 2
#define SZ(x) ((int)x.size())
#define ASSERT(v) {if (!(v)) {cerr<<"ASSERT FAIL @ "<<__FILE__<<":"<<__LINE__<<endl; exit(1);}}

using namespace std;
using namespace std::chrono;

typedef unsigned uint;
typedef unsigned long ulong;

class bgraph {
public:
    uint novL;
    uint novR;
    uint noE;
    uint max_degree_L{};
    uint max_degree_R{};
    vector<vector<uint>> adjL, adjR;

    bgraph(uint& e, uint& dl, uint& dr, vector<vector<uint>>& L, vector<vector<uint>>& R):
        noE(e), max_degree_L(dl), max_degree_R(dr), adjL(L), adjR(R)
        {
            novL = L.size();
            novR = R.size();
        }
    void sortAdj();
    static uint maxDegree(vector<vector<uint>>&);
    void calculate_number_of_edges();
    ~bgraph();
};

struct Biclique{
    vector<uint> subU;
    vector<uint> subV;
    ulong MBsize{};
};

class MB {
public:
    bgraph &g;
    string problem_type = "MEB";  // MVB MBB
    uint t_U, t_V;
    Biclique _C;
    string method = "iMBEA";
    bool init_MEB = false;

    MB(bgraph &bg, uint &tu, uint &tv, string mt, string pt):g(bg),
    t_U(tu), t_V(tv), method(std::move(mt)), problem_type(std::move(pt)) {  _C.MBsize = 0; };

    void iMBEA_init(bgraph &Gk,
                    vector<uint> &U,
                    vector<uint> &V,
                    vector<uint> &CV,
                    vector<uint> &XV);
    void iMBEA_Gk(bgraph &Gk,
                  uint &tkU,
                  uint &tkV,
                  vector<uint> &U,
                  vector<uint> &V,
                  vector<uint> &CV,
                  vector<uint> &XV);
    // TODO: MBC-algorithm
    void MBC_Gk(bgraph &Gk, uint tkU, uint tkV);
    void initMBC();
    ulong MBC_improved();

    bgraph Reduce(bgraph &Gt, uint tkU, uint tkV);
    static void Reduce1hop(bgraph &Gt, uint &tku, uint &tkv);
    static void Reduce2hop(vector<vector<uint>> &U, vector<vector<uint>> &V, uint &tku, uint &tkv);
};


#endif //RELEASEBCVIZ_MB_H
