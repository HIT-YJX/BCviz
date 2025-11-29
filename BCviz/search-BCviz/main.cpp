#include <iostream>
#include "bgraph.h"
#include "ooMBE.h"
#include <fstream>
#include <sstream>
#include <string>

/**
run:
./MBS [dataset] [problem_type] [algo_type] [s] [t] ([BCviz_file])

[dataset]: dataset filename
[problem_type]:  1. MEB (size=pxq); 2. MVB (size = p+q); 3. MBB (size=pxq and p=q)
[algo_type]:  1. BCviz; 2. BCviz+; 3. BCviz-
[s]: the threshold on the first vertex set of bipartite graph
[t]: the threshold on the second vertex set of bipartite graph
([BCviz_file]): the file name of BCviz index,
               [BCviz_file]="[dataset]_[problem_type]_[algo_type].txt" is default
*/


using namespace std;

void help();
void identify_file(const string& problem_type, const string& algo_type, string& file_tail);

int main(int argc, char  **argv) {
    string path, dataset, algo_type, BCviz_file, head_file, file_tail, problem_type;
    uint s=4, t=4;
    head_file = "../datasets/";
    BCviz_file = "../Index-results/";

    if (argc < 6) help();
    dataset = argv[1];
    problem_type = argv[2];
    algo_type = argv[3];
    s = atoi(argv[4]);
    t = atoi(argv[5]);
/*    head_file = "D:/Source Code/BCviz/datasets/";
    BCviz_file = "D:/Source Code/BCviz/Index-results/";
    cout<<"[datasets] [problem_type] [algo_name] [t_U] [t_V]: ";
    cin>>dataset>>problem_type>>algo_type>>s>>t;*/

    
    identify_file(problem_type, algo_type, file_tail);

    path = head_file + dataset;
    bgraph bg;
    bg.readEdgeList(path);


    BCviz_file += dataset + file_tail;
    if(argc == 7) BCviz_file = argv[6];
//if(algo_type != "MBC") cout<<"BCviz path: "<<BCviz_file<<endl;
//cout<<"algo_type: "<<algo_type<<", problem type: "<<problem_type<<endl;
    ooMBE alg(bg, s, t, algo_type, problem_type);
    alg.method_name = algo_type;
    if(algo_type == "MBC") alg.MBC_improved(s, t); // baseline: progressive reduction
    else alg.BCviz_MBC(BCviz_file);

    return 0;
}

void help() {
    cout<<"Options are error! Please dataset:"<<endl;
    cout<<"[algorithm] [datasets] [t_U] [t_V] [problem_type] [algo_name]"<<endl;
    cout<<"[problem_type] options: MEB, MVB, MBB"<<endl;
    cout<<"[algo_name] options: MBC, BCviz, BCviz+, BCviz-"<<endl;
    exit(0);
}

void identify_file(const string& problem_type, const string& algo_type, string& file_tail) {
    if(!(problem_type == "MEB"||problem_type == "MVB"||problem_type == "MBB")){
        cout<<"problem_type is error! (problem_type: "<<problem_type<<")"<<endl;
        exit(0);
    }
    if(!(algo_type == "BCviz"||algo_type == "BCviz+"||algo_type == "BCviz-"||algo_type == "MBC")) {
        cout<<"algo_type name is error! (algo_type: "<<algo_type<<")"<<endl;
        exit(0);
    }

    file_tail += "_" + problem_type + "_" + algo_type + ".txt";
}
