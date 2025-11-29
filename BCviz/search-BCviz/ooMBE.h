
#ifndef MBBP_OOMBE_H
#define MBBP_OOMBE_H

#include "bgraph.h"
#include <utility>
#include <functional>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <set>

using namespace std::chrono;
typedef unsigned int uint;
typedef unsigned long ulong;


struct Biclique{
    vector<uint> subU;
    vector<uint> subV;
    ulong MBsize;
};

struct SortedBySize
{
    unsigned row;
    unsigned vex;
    unsigned oneNodeSize;
    unsigned twoNodeSize;

    SortedBySize(unsigned r=0, unsigned v=0, unsigned o=0, unsigned t=0): row(r), vex(v), oneNodeSize(o), twoNodeSize(t) {};
};

class ooMBE {
public:
    bgraph &g;
    int *counts= 0;

    explicit ooMBE(bgraph &bg):g(bg){  _C.MBsize = 0; t_U = 1; t_V = 1; };
    ooMBE(bgraph &bg, uint &tu, uint &tv, const string mt, const string pt):g(bg),
            t_U(tu), t_V(tv), method(mt), problem_type(pt) {  _C.MBsize = 0; };
    unsigned long int nomb = 0;

    /*** baseline* */
    bgraph tmpG; // for BCviz to reduce graph
    string problem_type = "MEB";
    uint t_U, t_V;
    Biclique _C;
    string method = "ooMBEA";
    bool init_MEB = false;
    uint sizeL, sizeR;
    vector<vector<SortedBySize>> ordered_vex; // record the ordered vertex. (row, vex, onevexsize, twovexsize)
    uint *posMapidL = 0, *posMapidR = 0, *idMapposL = 0, *idMapposR = 0;
    bool *visitL = 0, *visitR = 0;
    bool switchLR = false;
    uint *L_order = 0, *R_order = 0;
    string method_name;

    /// tool
    void printMapVec(vector<uint> &vec, uint row);
    bgraph Reduce(bgraph &Gt, uint tkU, uint tkV);
    static void Reduce1hop(bgraph &Gt, uint &tku, uint &tkv);
    static void Reduce2hop(vector<vector<uint>> &U, vector<vector<uint>> &V, uint &tku, uint &tkv);
    void calculate_sizeLR(bgraph &Gt);
    static void printMEB(Biclique &bc);
    void printMapMEB(Biclique &bc);

    // TODO: iMBEA-algorithm
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
    void MBC_improved(uint tkU, uint tkV);
    /** baseline end */

    /** Index search */
    void readBCviz(string filename);
    void printBCviz();
    void initMBC_Gk(bgraph &Gt);
    ulong BicliqueSize(uint l, uint r);
    void MBC_improved_Gk(bgraph &Gt, uint tkU, uint tkV);
    void initBiclique2Vertex(bgraph& Gk, uint &vv);
    void initBiclique(bgraph& Gk);
    void initGraph();
    void getCandidateVertices(uint &max_ts, uint &pos_l, uint &pos_r, uint interval_size,
                              uint &num_L, uint &numR, uint pos_i,
            bool *reduce_L, bool *reduce_R, const bool *core_L, const bool *core_R);

    void init_BCviz();
    void init_BCviz_opt();
    void baseline_initMEB_BCviz();
    void MEB_BCvizReduce(bgraph &Gk);
    void MVB_BCvizReduce(bgraph &Gk);
    void MBB_BCvizReduce(bgraph &Gk);
    void BCviz_MBC(string& LSIfile);
    /** Index search end */

    virtual ~ooMBE();
};



#endif //MBBP_OOMBE_H
