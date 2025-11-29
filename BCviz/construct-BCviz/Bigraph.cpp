#include"Bigraph.h"

Vertex::Vertex(): id(-1), rownum(-1), deg(-1){}

Vertex::Vertex(unsigned i, unsigned r, unsigned d): id(i), rownum(r), deg(d) {}

Vertex::~Vertex(){
    nei.clear();
}

Bigraph::~Bigraph() {
    for(int i=0;i<ROW_SIZE;i++){
        id_map_pos[i].clear();
        vertices[i].clear();
    }
    deg_num.clear();
}


bool Bigraph::all_num(string &s) {
    for (int i = 0; i < SZ(s); i++) if (!(s[i] >= '0' && s[i] <= '9')) return false;
    return true;
}

void Bigraph::add_undirected_edge(unsigned u, unsigned v)
{
    if(edges.count({u,v})) return;
    edges.insert({u,v});

    if(!id_map_pos[0].count(u)){
        id_map_pos[0][u] = vertices[0].size();
        Vertex uu(u, 0, 0);
        vertices[0].emplace_back(uu);
    }
    if(!id_map_pos[1].count(v)){
        id_map_pos[1][v] = vertices[1].size();
        Vertex vv(v, 1, 0);
        vertices[1].emplace_back(vv);
    }

    ///add edge;
    unsigned pos_u, pos_v;
    pos_u = id_map_pos[0][u];
    pos_v = id_map_pos[1][v];
    vertices[0][pos_u].nei.emplace_back(pos_v);
    vertices[1][pos_v].nei.emplace_back(pos_u);

}

void Bigraph::read_graph(const string& filename) //, uint scale
{
    FILE* fp = nullptr;
    fp = freopen(filename.c_str(), "r", stdin);
    if(!fp){
        cout<<"can't open dataset file!"<<endl;
        return;
    }

    string s;
    while(getline(cin, s)){
        stringstream ss;
        ss << s;
        vector <string> vec_str;
        for(string z; ss >> z; vec_str.emplace_back(z));
        if(SZ(vec_str) >= 2){
            bool is_all_number = true;
            for(unsigned i = 0; i < min(2, SZ(vec_str)); i++)
                is_all_number &= all_num(vec_str[i]);
            if(is_all_number){
                unsigned u, v;
                ss.clear(); ss << vec_str[0]; ss >> u;
                ss.clear(); ss << vec_str[1]; ss >> v;
                add_undirected_edge(u, v);
            }
        }
    }

    num_edges = edges.size();
    num_vertices = vertices[0].size() + vertices[1].size();

/*cout << "end of read graph" << endl;
cout << "#vertices = " << num_vertices << endl;
cout << "#edges = " << num_edges << endl;*/

    fclose(stdin);
    cin.clear();
    edges.clear();
}

int Bigraph::compute_deg()
{
    max_deg = 0;
    int sum_deg;
    double average_deg = 0;
    int row = 0;
    for(int i=0; i<ROW_SIZE; i++){
        sum_deg = 0;
        for(auto &it: vertices[i]){
            it.deg = it.nei.size();
            //it.com_deg = it.deg;  ///used for peel 1 degree
            sum_deg += it.deg;
            if(max_deg < it.deg) {
                max_deg = it.deg;
                row = i;
            }
        }
        average_deg += sum_deg;
        ///cout<<"the average degree of "<<i<<"th row is "<<double(sum_deg) / vertices[i].size() <<endl;
    }
    average_deg = average_deg / num_vertices;
    ///cout<<"the max degree is "<<max_deg<<endl;
    ///cout<<"the average degree is "<<average_deg<<endl;

    for(int i=0;i<=max_deg;i++){
        deg_num.emplace_back(unsigned(0));
    }
    for(int i=0;i<ROW_SIZE;i++){
        for(auto &it: vertices[i]){
            deg_num[it.deg]++;
        }
    }
    return row;
}

void Bigraph::record_deg_num(string filename)
{
    string degreefile = filename + "_Degree.infor";
    FILE* fp = nullptr;
    fp = freopen(degreefile.c_str(), "w", stdout);
    if(!fp){
        cout<<"can't open dataset file!"<<endl;
        return;
    }


    cout<<"number of vertices:"<<num_vertices<<", number of edges:"<<num_edges<<endl;
    cout<<"degree:num"<<endl;
    for(int i=1;i<=max_deg;i++){
        if(deg_num[i])
            cout<<i<<": "<<deg_num[i]<<endl;
    }
    fclose(stdout);
}

/*void Bigraph::print_biclique()
{
    for(int i=0; i<ROW_SIZE; i++)
        for(auto &ve: vertices[i]){
            for(auto &pq: ve.pxq){
                cout<<pq.first<<"x"<<pq.second<<", ";
            }
            cout<<"\n\n";
        }
}*/

void Bigraph::print_edges()
{
    for(auto &it: edges){
        cout<<it.first<<" "<<it.second<<endl;
    }
}

void Bigraph::print_vertices()
{
    for(int i=0; i<ROW_SIZE; i++){
        cout<<"the "<<i<<"th row:"<<endl;
        for(auto &it: vertices[i]){
            cout<<it.id<<endl;
        }
    }
}


