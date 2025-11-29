#include <iostream>
#include "BCVIZ.h"
#include "Bigraph.h"

#include <iostream>
#include <cstdio>
#include <cmath>
#include <random>
#include <chrono>
#include <bitset>

/**
run:
./BCviz [dataset] [problem_type] [algo_type] ([s_min]) ([t_min])

[dataset]: dataset filename
[problem_type]:  1: MEB (size=pxq); 2: MVB (size = p+q); 3. MBB (size=pxq and p=q)
[algo_type]:  0: BCviz; 1: BCviz+; 2: BCviz-
([s_min]): the minimum threshold on the first vertex set of bipartite graph, s_min=3 by default
([t_min]): the minimum threshold on the second vertex set of bipartite graph, t_min=3 by default
*/

using namespace std;
using namespace std::chrono;


void help();
void setting(int &pro, int &algo, const string& pro_t, const string& algo_t);

int main(int argc, char *argv[])
{
    uint s_min=3, t_min=3;
    string dataset, path, algo_type, problem_type, result_tail;
    int problem_ID, algo_ID, method_type = 3;
    path = "../datasets/";

    if(argc < 4) help();
    dataset = argv[1];
    problem_type = argv[2];
    algo_type = argv[3];
    path += dataset;
    if(argc == 6) {
        s_min = atoi(argv[4]);
        t_min = atoi(argv[5]);
    }

    setting(problem_ID, algo_ID, problem_type, algo_type);
    BCVIZ example(dataset, problem_ID, method_type, problem_type, s_min, t_min);
    example.read_graph(path);


    // prepare graph and order vertices
    auto start = high_resolution_clock::now();

    BCVIZ::reduce1hop(example.bigraph, s_min, t_min);
    if(algo_type == "BCviz") {
        example.offline_BCviz_prepare_graph();
        example.BCviz_Construct();
    }
    else if(algo_type == "BCviz+") {
        example.BCviz1_prepare_graph();
        example.BCviz1_Construct();
    }
    else if (algo_type == "BCviz-") {
        example.BCviz2_prepare_graph();
        example.BCviz2_Construct();
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    double order_time = duration.count() * 1.0 / 1000000.0;


    /// output: [problem_type] [algo_type] [dataset] [s_min] [t_min] [total_time]
    cout<<problem_type<<" "<<algo_type<<" "<<dataset
        <<" "<<s_min<<" "<<t_min<<" "<<order_time<<endl;

    // save index in file
    result_tail = algo_type;
    if(s_min != 3 || t_min != 3)
        result_tail += "_" + to_string(s_min) + "_" + to_string(t_min);
    example.BCviz_out_result(result_tail);

    return 0;
}

void help() {
    printf("parameter number is less than 4!\n");
    printf("Please input: [execution] [dataset] [problem_type] [algo_type] ");
    printf(" ([s_min]) ([t_min])\n");
    exit(0);
}

void setting(int &pro, int &algo, const string& pro_t, const string& algo_t) {
    if(pro_t == "MEB") {
        pro = 1;
    }
    else if(pro_t == "MVB") {
        pro = 2;
    }
    else if(pro_t == "MBB") {
        pro = 3;
    }
    else {
        cout<<"problem_type is error (problem_type="<<pro_t<<")! \n Please choose: MEB, MVB, MBB."<<endl;
        exit(0);
    }

    if(algo_t == "BCviz") {
        algo = 0;
    }
    else if(algo_t == "BCviz+") {
        algo = 1;
    }
    else if (algo_t == "BCviz-") {
        algo = 2;
    }
    else {
        cout<<"algo_type is error (algo_type="<<algo_t<<")! \n Please choose: BCviz, BCviz+, BCviz-"<<endl;
        exit(0);
    }
}
