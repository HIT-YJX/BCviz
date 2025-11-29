#ifndef OOBCVIZ_BCVIZ_H
#define OOBCVIZ_BCVIZ_H

#include "Bigraph.h"
#include "MB.h"

#include <list>
#include <unordered_map>
#include <map>
#include <queue>
#include <utility>
#include <vector>
#include <unordered_set>
#include <bitset>
#include <cstring>
#include <cmath>
#include <string>
#include <stack>
#include <set>
#include <list>
#include <chrono>

//#define DEBUG
//#define EST_DEBUG
//#define HYBRID

using namespace std::chrono;

typedef unsigned uint;
typedef unsigned long ulong;


struct SortedVertexByDegree
{
    unsigned vertex;
    unsigned degree;
    SortedVertexByDegree(unsigned v, unsigned d): vertex(v), degree(d) {}
    SortedVertexByDegree() = default;

    inline bool operator < (const SortedVertexByDegree &other) const {
        return degree < other.degree || (degree == other.degree && vertex < other.vertex); // non-descending by degree
    }
};

struct SortedBySize_BCviz
{
    unsigned row;
    unsigned vex;
    unsigned oneNodeSize;
    unsigned twoNodeSize;

    SortedBySize_BCviz(unsigned r=0, unsigned v=0, unsigned o=0, unsigned t=0): row(r), vex(v), oneNodeSize(o), twoNodeSize(t) {};
};

struct cmp
{
    bool operator () (const SortedBySize_BCviz &a, const SortedBySize_BCviz &b) const{
        return a.twoNodeSize > b.twoNodeSize || (a.twoNodeSize == b.twoNodeSize && a.oneNodeSize > b.oneNodeSize) || \
        (a.twoNodeSize == b.twoNodeSize && a.oneNodeSize == b.oneNodeSize && a.vex < b.vex) || \
        (a.twoNodeSize == b.twoNodeSize && a.oneNodeSize == b.oneNodeSize && a.vex == b.vex && a.row < b.row);
    }
};


class BCVIZ {
public:
    Bigraph bigraph;
    unsigned main_row;
    unsigned sub_row;
    string dataset, pro_name;
    int sizetype;
    int method_type;

/// for ordering vetex
    vector<ulong> vec_size[ROW_SIZE];  /// record onesize
    vector<ulong> vec_tsize[ROW_SIZE];  /// record twosize
    vector<vector<SortedBySize_BCviz>> ordered_vex; // record the ordered vertex. (row, vex, onevexsize, twovexsize)
    vector<bool> vec_vis[ROW_SIZE]; //record visited vex
    vector<bool> is_in_O[ROW_SIZE]; //record visited vex
    map<pair<uint, uint>, ulong> edgeCoh;

/// for new framework
    vector<vector<unsigned>> THnei; // no use
    uint t_U, t_V;


    BCVIZ(string d, int& s, int& m, string p, uint &tu, uint &tv):
    dataset(std::move(d)), sizetype(s), method_type(m), pro_name(std::move(p)), t_U(tu), t_V(tv){};
    ~BCVIZ();

/// read graph
    void read_graph(string &inputfile); //, uint scale

/*** Estimate Methods ***/
    /// TODO:  1.structure design  2.prepare_onesize  3.calc_twosize  4.print_ordering_vertex
    static void reduce1hop(Bigraph &Gt, unsigned &tku, unsigned &tkv);

    void BCviz1_prepare_graph();
    void offline_BCviz_prepare_graph();
    void BCviz2_prepare_graph();
    unsigned BCviz_EdgeCohesivenessShrinkGraph(const unsigned uu, const unsigned vv, const unsigned uu_row);
    unsigned BCviz_EdgeCohesivenessExact(const unsigned uu, const unsigned vv, const unsigned uu_row);


/// ordering framework
    void BCviz_Construct();
    void BCviz1_Construct();
    void BCviz2_Construct();

/// tools
    static void printVec(vector<uint> &vec);
    static void printGraph(Bigraph &PG);
    static void printIDGraph(Bigraph &PG);
    void BCviz_out_result(const string& tail);
};


#endif //OOBCVIZ_BCVIZ_H
