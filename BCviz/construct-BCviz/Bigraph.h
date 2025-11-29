#ifndef OOBCVIZ_BIGRAPH_H
#define OOBCVIZ_BIGRAPH_H

#include <algorithm>
#include <iostream>
#include <vector>
#include <algorithm>
#include <queue>
#include <cstdio>
#include <list>
#include <string>
#include <unordered_map>
#include <set>
#include <sstream>


#define ROW_SIZE 2
#define SZ(x) ((int)x.size())
#define ASSERT(v) {if (!(v)) {cerr<<"ASSERT FAIL @ "<<__FILE__<<":"<<__LINE__<<endl; exit(1);}}
#define INFO(...) do { ___debug( #__VA_ARGS__,  __VA_ARGS__,cout); cout<<endl; if(__HEAD_H_FOUT.is_open()){___debug( #__VA_ARGS__,  __VA_ARGS__,__HEAD_H_FOUT); __HEAD_H_FOUT<<endl;}  } while(0)


using namespace std;

typedef unsigned uint;
typedef unsigned long long lint;

struct Vertex
{
    unsigned id;
    unsigned rownum;
    vector<unsigned> nei;

    unsigned deg; // the number of neighbor

    Vertex();
    Vertex(unsigned i, unsigned r, unsigned d);
    ~Vertex();
};


class Bigraph{
public:
    unsigned num_edges;
    unsigned num_vertices;

    vector<Vertex> vertices[ROW_SIZE];
    unordered_map<unsigned, unsigned> id_map_pos[ROW_SIZE];
    set< pair<unsigned, unsigned> > edges;
    vector<unsigned> deg_num;
    unsigned max_deg;

    void read_graph(const string& filename); //, uint scale
    void add_undirected_edge(unsigned u, unsigned v);
    int compute_deg(); // return row that includes vex with max degree

    ~Bigraph();


    //tools
    static bool all_num(string &s);
    void print_edges();
    void print_vertices();
    void record_deg_num(string filename);


};


#endif //OOBCVIZ_BIGRAPH_H
