#include "bgraph.h"

using namespace std;

bgraph::bgraph() {
    novL = 0;
    novR = 0;
    noE = 0;
    offset = 0;
    max_degree_R = 0;
    max_degree_L = 0;
}

bgraph::~bgraph() {
    for(auto itvec: adjL) {
        itvec.clear();
    }
    adjL.clear();
    for(auto itvec: adjR) {
        itvec.clear();
    }
    adjR.clear();

    cohesionL.clear();
    cohesionR.clear();

    for(auto itset: adjHL) {
        itset.clear();
    }
    adjHL.clear();
    for(auto itset: adjHR) {
        itset.clear();
    }
    adjHR.clear();
}

void bgraph::dispalyAdj() {
    if(adjL.size()!=0||adjR.size()!=0){
        cout<<"The adj list for L"<<endl;
        for(int i=0; i<adjL.size();i++){
            cout<<i<<":";
            for(int j=0; j<adjL[i].size(); j++){
                cout<<adjL[i][j]<<",";
            }
            cout<<endl;
        }
        cout<<"The adj list for R"<<endl;
        for(int i=0; i<adjR.size();i++){
            cout<<i<<":";
            for(int j=0; j<adjR[i].size(); j++){
                cout<<adjR[i][j]<<",";
            }
            cout<<endl;
        }
    }else{
        cout<<"This bipartite graph is an empty set!";
    }

}

void bgraph::sortAdj() {
    for(auto &list:adjL){
        sort(list.begin(),list.end());
    }
    for(auto &list:adjR){
        sort(list.begin(),list.end());
    }
}

bool bgraph::checkNumber(string &ss) {
    for (int i = 0; i < SZ(ss); i++) if (!(ss[i] >= '0' && ss[i] <= '9')) return false;
    return true;
}

void bgraph::Stringsplit(const string& str, const string& split, vector<string>& res)
{
    char* strc = new char[str.size() + 1];
    char* pt;
    strcpy(strc, str.c_str());
    char* temp = strtok_r(strc, split.c_str(), &pt);
    while (temp != nullptr)
    {
        res.emplace_back(temp);
        temp = strtok_r(nullptr, split.c_str(), &pt);
    }
    delete[] strc;
}

void bgraph::readEdgeList(const string& file) {
    ifstream myfile (file.c_str());
    string line,u,v;
    string split = " \t";

    vector<pair<uint,uint>> edge;
    uint numL = 0;
    uint numR = 0;
    int num_line = 0;

    if (myfile.is_open()){
        while ( getline(myfile, line) ) {
            num_line++;
            vector<string> vecstr;
            Stringsplit(line, split, vecstr);
            if(vecstr.size() >= 2) {
                u = vecstr[0];
                v = vecstr[1];
                if(!(checkNumber(u) && checkNumber(v))) {
                    if(num_line >= 4) {
                        cout<<"The dataset is not the file format: u1(integer) v1(integer). "<<endl;
                        exit(0);
                    }
                    continue;
                }

                int ul = stoi(u)-1;  //todo: vid = pos + 1
                int vr = stoi(v)-1;  //todo: vid = pos + 1
                if(ul>numL){
                    numL = ul;
                }
                if(vr>numR){
                    numR=vr;
                }
                edge.emplace_back(make_pair(ul,vr));

            }
            else if(num_line>=4){
                cout<<"The dataset is not the file format: u1 v1. "<<endl;
                exit(0);
            }
            else {
                continue;
            }


        }
    }

    vector<vector<uint>> adjL_prime(numL+1);
    vector<vector<uint>> adjR_prime(numR+1);
    /// reduce multiple edges;
    map<pair<uint, uint>, bool> reduceMap;
    for(auto &e:edge){
        if(reduceMap[e]) continue;
        reduceMap[e] = true;
        adjL_prime[e.first].emplace_back(e.second);
        adjR_prime[e.second].emplace_back(e.first);
    }
    //reduceMap.clear();
    adjL = adjL_prime;  //it is copying but this is the loading time
    adjR = adjR_prime;
    novL = adjL.size();
    novR = adjR.size();
    noE = edge.size();
    max_degree_L = maxDegree(adjL);
    max_degree_R = maxDegree(adjR);
    myfile.close();
    datasets = file;

    cohesionL.resize(novL);
    cohesionR.resize(novR);


    /*cout<<"--------------read graph--------------"<<endl;
    cout<<"number of vertices in U: "<<novL<<endl;
    cout<<"number of vertices in V: "<<novR<<endl;
    cout<<"number of edges: "<<noE<<endl;
    cout<<"max degree of U: "<<max_degree_L<<endl;
    cout<<"max degree of V: "<<max_degree_R<<endl<<endl;*/
}

void bgraph::switchLR() {
    adjL.swap(adjR);
    int tmp = max_degree_L;
    max_degree_L = max_degree_R;
    max_degree_R = tmp;
    tmp = novL;
    novL = novR;
    novR = tmp;
    isSwitch=!isSwitch;
    /*cout<<"--------------Switch graph--------------"<<endl;
    cout<<"number of vertices in U: "<<novL<<endl;
    cout<<"number of vertices in V: "<<novR<<endl;
    cout<<"number of edges: "<<noE<<endl;
    cout<<"max degree of U: "<<max_degree_L<<endl;
    cout<<"max degree of V: "<<max_degree_R<<endl;*/
}

uint bgraph::maxDegree(vector<vector<uint>> &adj) {
    uint max=0;
    for(auto &neighbour_list:adj){
        if(neighbour_list.size()>max){
            max = neighbour_list.size();
        }
    }
    return max;
}

void bgraph::calculate_number_of_edges() {
    if(noE==0){
        if(novL<novR){
            uint m=0;
            for(auto &neighbours:adjL){
                m=m+neighbours.size();
            }

            noE=m;
        }else{
            uint m=0;
            for(auto& neighbours:adjR){
                m+=neighbours.size();
            }
            noE = m;
        }

    }
}



