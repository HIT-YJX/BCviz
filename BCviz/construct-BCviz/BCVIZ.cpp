#include "BCVIZ.h"

BCVIZ::~BCVIZ() {
    for(int i=0;i<ROW_SIZE;i++){
        vec_size[i].clear();
        vec_tsize[i].clear();
        vec_vis[i].clear();
    }

    for(auto & i : ordered_vex){
        i.clear();
    }

    if(!edgeCoh.empty()) edgeCoh.clear();
    /*for(auto & j : THnei)
        j.clear();*/
}

void BCVIZ::read_graph(string &inputfile)
{
    bigraph.read_graph(inputfile);
//cout<<"read graph successfully!"<<endl;
}

void BCVIZ::reduce1hop(Bigraph &Gt, unsigned &tku, unsigned &tkv) {
    if(tku<=1 && tkv<=1) return;
    uint num_u = Gt.vertices[0].size();
    uint num_v = Gt.vertices[1].size();
    uint *degree_u = new uint[num_u];
    uint *degree_v = new uint[num_v];
    uint *to_remove_u = new uint[num_u];
    uint *to_remove_v = new uint[num_v];
    bool *removed_u = new bool[num_u];
    bool *removed_v = new bool[num_v];
    for(uint i = 0; i < num_u; i++) {
        degree_u[i] = Gt.vertices[0][i].nei.size();
    }
    for(uint i = 0; i < num_v; i++) {
        degree_v[i] = Gt.vertices[1][i].nei.size();
    }
    for(uint i = 0; i < num_u; i++){
        to_remove_u[i] = 0;
        removed_u[i] = false;
    }
    for(uint i = 0; i < num_v; i++){
        to_remove_v[i] = 0;
        removed_v[i] = false;
    }

    uint start_idx_u, end_idx_u, start_idx_v, end_idx_v;
    start_idx_u = end_idx_u = start_idx_v = end_idx_v = 0;

    // step 1: collect the initial vertices with degree less than q in left and p in right
    for(uint i = 0; i < num_u; i++){
        if(degree_u[i] < tkv){
            to_remove_u[end_idx_u++] = i;
            removed_u[i] = true;
        }
    }
    for(uint i = 0; i < num_v; i++){
        if(degree_v[i] < tku){
            to_remove_v[end_idx_v++] = i;
            removed_v[i] = true;
        }
    }

    // step 2: recursively remove all vertices with degree less than q in left and p in right, i.e., (q,p)-core
    while(start_idx_u != end_idx_u || start_idx_v != end_idx_v){

        while(start_idx_u != end_idx_u) {
            uint uu = to_remove_u[start_idx_u++];

            //cout << "remove : " << uu << endl;
            for(uint i = 0; i < Gt.vertices[0][uu].nei.size(); i++){
                uint other = Gt.vertices[0][uu].nei[i];
                if(!removed_v[other]){
                    degree_v[other]--;

                    if(degree_v[other] < tku){
                        to_remove_v[end_idx_v++] = other;
                        removed_v[other] = true;
                    }
                }
            }
            Gt.vertices[0][uu].nei.clear();
            degree_u[uu] = 0;
        }

        while(start_idx_v != end_idx_v) {
            uint vv = to_remove_v[start_idx_v++];

            //cout << "remove : " << vv << endl;
            for(uint i = 0; i < Gt.vertices[1][vv].nei.size(); i++){
                uint other = Gt.vertices[1][vv].nei[i];
                if(!removed_u[other]){
                    degree_u[other]--;

                    if(degree_u[other] < tkv){
                        to_remove_u[end_idx_u++] = other;
                        removed_u[other] = true;
                    }
                }
            }
            Gt.vertices[1][vv].nei.clear();
            degree_v[vv] = 0;
        }

    }

    // update adj
    for(uint i = 0; i < num_u; i++){
        if(!removed_u[i]){
            vector<uint> new_vec(degree_u[i]);
            uint idx = 0;
            for(uint j = 0; j < Gt.vertices[0][i].nei.size(); j++){
                if(!removed_v[Gt.vertices[0][i].nei[j]]){
                    new_vec[idx++] = Gt.vertices[0][i].nei[j];
                }
            }
            Gt.vertices[0][i].nei.clear();
            Gt.vertices[0][i].nei = new_vec;
        }
    }
    for(uint i = 0; i < num_v; i++){
        if(!removed_v[i]){
            vector<uint> new_vec(degree_v[i]);
            uint idx = 0;
            for(uint j = 0; j < Gt.vertices[1][i].nei.size(); j++){
                if(!removed_u[Gt.vertices[1][i].nei[j]]){
                    new_vec[idx++] = Gt.vertices[1][i].nei[j];
                }
            }
            Gt.vertices[1][i].nei.clear();
            Gt.vertices[1][i].nei = new_vec;
        }
    }

    delete[] degree_u;
    delete[] degree_v;
    delete[] to_remove_u;
    delete[] to_remove_v;
    delete[] removed_u;
    delete[] removed_v;
}


void BCVIZ::BCviz1_prepare_graph() {
    /// initial row
    bigraph.compute_deg();
    main_row = 0;
    sub_row = 1;

    /// initial data structure
    vec_tsize[0].resize(bigraph.vertices[0].size(), 0);
    vec_tsize[1].resize(bigraph.vertices[1].size(), 0);
    vec_vis[0].resize(bigraph.vertices[0].size(), false);
    vec_vis[1].resize(bigraph.vertices[1].size(), false);
    is_in_O[0].resize(bigraph.vertices[0].size(), false);
    is_in_O[1].resize(bigraph.vertices[1].size(), false);

    ///cout<<"BCviz prepare graph successfully!"<<endl;
}

void BCVIZ::offline_BCviz_prepare_graph() {
    /// initial row
    bigraph.compute_deg();
    main_row = 0;
    sub_row = 1;

    /// initial data structure
    vec_tsize[0].resize(bigraph.vertices[0].size(), 0);
    vec_tsize[1].resize(bigraph.vertices[1].size(), 0);
    vec_vis[0].resize(bigraph.vertices[0].size(), false);
    vec_vis[1].resize(bigraph.vertices[1].size(), false);
    is_in_O[0].resize(bigraph.vertices[0].size(), false);
    is_in_O[1].resize(bigraph.vertices[1].size(), false);

    /// compute \sigma(u,v)
    ulong coh = 0;
    for(uint u = 0; u < bigraph.vertices[0].size(); u++) {
        for(auto v: bigraph.vertices[0][u].nei) {
            coh = BCviz_EdgeCohesivenessExact(u, v, 0);
            pair<uint, uint> pr = make_pair(u,v);
            edgeCoh[pr] = coh;
        }
    }

    ///cout<<"BCviz prepare graph successfully!"<<endl;
}

void BCVIZ::BCviz2_prepare_graph() {
    /// initial row
    bigraph.compute_deg();
    main_row = 0;
    sub_row = 1;

    /// initial data structure
    vec_tsize[0].resize(bigraph.vertices[0].size(), 0);
    vec_tsize[1].resize(bigraph.vertices[1].size(), 0);
    vec_vis[0].resize(bigraph.vertices[0].size(), false);
    vec_vis[1].resize(bigraph.vertices[1].size(), false);
    is_in_O[0].resize(bigraph.vertices[0].size(), false);
    is_in_O[1].resize(bigraph.vertices[1].size(), false);

    ///cout<<"BCviz prepare graph successfully!"<<endl;
}


/// shrinking graph O(m)
unsigned BCVIZ::BCviz_EdgeCohesivenessShrinkGraph(const unsigned uu, const unsigned vv, const unsigned uu_row)
{
    unsigned res, edgen, vertexn, u_num, v_num;
    unsigned lower[2]={t_V, t_U};
    vector<Vertex>& urow_vertices = bigraph.vertices[uu_row];
    vector<unsigned>& uu_nei = urow_vertices[uu].nei;
    unsigned vv_row = 1-uu_row;
    vector<Vertex>& vrow_vertices = bigraph.vertices[vv_row];
    vector<unsigned>& vv_nei = vrow_vertices[vv].nei;
    res = edgen = vertexn = u_num = v_num = 0;
#ifdef EST_DEBUG
    cout<<"uu_row="<<uu_row<<", uu_id="<<bigraph.vertices[uu_row][uu].id<<", vv_id=";
    cout<<bigraph.vertices[vv_row][vv].id<<endl;
#endif

    /// special situation1: vex has just one degree.
    //if(uu_nei.size() == 0 || vv_nei.size() == 0) return (unsigned)0;
    if(uu_nei.size() < lower[uu_row] || vv_nei.size() < lower[vv_row]) return (unsigned)0;


    /// 1.get center graph O(m)
    Bigraph centerG;
    unsigned max_v = *max_element(uu_nei.begin(), uu_nei.end());
    unsigned max_u = *max_element(vv_nei.begin(), vv_nei.end());
    vector<unsigned> u_map_id(max_u+1, 0);
    vector<unsigned> v_map_id(max_v+1, 0);
    for(auto ui: vv_nei) u_map_id[ui] = ++u_num;
    for(auto vj: uu_nei) v_map_id[vj] = ++v_num;
    centerG.vertices[uu_row].resize(u_num+1);
    centerG.vertices[vv_row].resize(v_num+1);
    for(auto ui: vv_nei) {
        for(auto vj: bigraph.vertices[uu_row][ui].nei) {
            if(ui > max_u || vj > max_v) continue; /// O(1) map
            uint utmp = u_map_id[ui];
            uint vtmp = v_map_id[vj];
            if(utmp==0 || vtmp==0) continue;
            centerG.vertices[uu_row][ utmp ].nei.emplace_back( vtmp );
            //centerG.vertices[uu_row][ utmp ].id = ui;
            centerG.vertices[vv_row][ vtmp ].nei.emplace_back( utmp );
            //centerG.vertices[vv_row][ vtmp ].id = vj;

            //edgen++;
        }
    }

#ifdef EST_DEBUG
    cout<<"get center graph is ok!"<<endl;
    //printGraph(centerG);
    /*cout<<"End of getting center graph"<<endl;
    cout << "#vertices = " << u_num + v_num << endl;
    cout << "#edges = " << edgen << endl;*/
#endif

    /// 2.core Reduce: get (\beta, \alpha)-core O(m)
    reduce1hop(centerG, lower[vv_row], lower[uu_row]);
#ifdef EST_DEBUG
    cout<<"core reduce is ok!"<<endl;
    //printGraph(centerG);
#endif

    /// 3.swap O(n)
    uint core_u_num, core_v_num;
    core_u_num = core_v_num = 0;
    for(auto &uit : centerG.vertices[uu_row]) /// pos from 1
        if(!uit.nei.empty()) core_u_num++;
    for(auto &vit : centerG.vertices[vv_row])
        if(!vit.nei.empty()) core_v_num++;
    //if(u_num == 0 || v_num == 0) return 0;
    if(core_v_num < lower[uu_row] || core_u_num < lower[vv_row]) return (unsigned)0;


    unsigned u_row, v_row;
    u_row = uu_row;
    v_row = vv_row;

    if(core_u_num > core_v_num) {
        u_row = vv_row;
        v_row = uu_row;

        uint tmp = core_u_num;
        core_u_num = core_v_num;
        core_v_num = tmp;

        tmp = u_num; // todo: new code
        u_num = v_num;
        v_num = tmp;

#ifdef EST_DEBUG
        cout<<"swap U and V"<<endl;
#endif
    }

#ifdef EST_DEBUG
    cout<<"core_u_num="<<core_u_num<<", core_v_num="<<core_v_num<<endl;
#endif


    /// 4. prepare data structure of degree O(n)
    vertexn = core_u_num + core_v_num;
    edgen = 0;

    // sort degree
    auto *u_degree = new SortedVertexByDegree[core_u_num+1];
    auto *v_degree = new SortedVertexByDegree[core_v_num+1];
    unsigned u_max_degree, v_max_degree, u_degree_st, u_degree_ed, v_degree_st, v_degree_ed;
    u_max_degree = v_max_degree = u_degree_st = u_degree_ed = v_degree_st = v_degree_ed = 0;
    for(uint uid = 1; uid < centerG.vertices[u_row].size(); uid++) {
        uint ud = centerG.vertices[u_row][uid].nei.size(); // degree of u
        if(ud > 0) { /// todo: ud >= lower[u_row] ?
            SortedVertexByDegree node(uid, ud);
            u_degree[u_degree_ed++] = node;
            edgen += ud;
        }
    }
#ifdef HYBRID
    //cout<<"l="<<core_u_num<<", r="<<core_v_num<<", m="<<edgen<<" ";
#endif
#ifdef EST_DEBUG
    cout<<"edgen="<<edgen<<endl;
#endif
    for(uint vid = 1; vid < centerG.vertices[v_row].size(); vid++) {
        uint vd = centerG.vertices[v_row][vid].nei.size(); // degree of v
        if(vd > 0) { /// todo: vd >= lower[v_row] ?
            SortedVertexByDegree node(vid, vd);
            v_degree[v_degree_ed++] = node;
        }
    }
#ifdef EST_DEBUG
    cout<<"u_degree_ed="<<u_degree_ed<<", v_degree_ed="<<v_degree_ed<<endl;
#endif
    if(u_degree_ed < lower[u_row] || v_degree_ed < lower[v_row]) return res; /// check u_degree[u_degree_ed-1]
    sort(u_degree, u_degree + u_degree_ed); // non-descending order by degree
    sort(v_degree, v_degree + v_degree_ed); // non-descending order by degree
#ifdef EST_DEBUG
    cout<<"u_degree:"<<endl;
    for(int i = 0; i < u_degree_ed; i++){
        cout<<u_degree[i].degree<<",";
    }
    cout<<endl;
    cout<<"v_degree:"<<endl;
    for(int i = 0; i < v_degree_ed; i++){
        cout<<v_degree[i].degree<<",";
    }
    cout<<endl;
#endif

    // get pointer of degree array
    // get array of vertex mapping pos
    auto *u_pos = new uint[u_num+1]; /// uid form 1 to u_num
    auto *v_pos = new uint[v_num+1];

    u_max_degree = u_degree[u_degree_ed-1].degree;
    v_max_degree = v_degree[v_degree_ed-1].degree;
    double density = edgen*1.0/core_u_num/core_v_num;
#ifdef HYBRID
/*    cout<<"du="<<u_max_degree<<", dv="<<v_max_degree<<endl;
    cout<<"density="<<edgen*1.0/core_u_num/core_v_num<<", p="
        <<edgen*1.0/sqrt(core_u_num*1.0)/sqrt(core_v_num*1.0)<<", cost_est="
        <<(u_max_degree+lower[v_row])*lower[u_row]*1.0/2<<endl;*/
#endif
    if(density == 1) {
        if(pro_name=="MEB") return core_u_num * core_v_num;
        if(pro_name=="MVB") return core_u_num + core_v_num;
        if(pro_name=="MBB") return min(core_u_num * core_u_num, core_v_num * core_v_num);
    }

    uint u_pointer[u_max_degree+2]; /// u_pointer[0]~u_pointer[1] are the pointer of vertices with degree=0
    uint v_pointer[v_max_degree+2];
    for(uint i=0; i<=u_max_degree; i++) {
        u_pointer[i] = 0;
    }
    for(uint i=0; i<=v_max_degree; i++) {
        v_pointer[i] = 0;
    }

    uint ite_d = u_degree[u_degree_st].degree;
    u_pointer[ite_d] = u_degree_st;
    for(uint pos = u_degree_st; pos < u_degree_ed; pos++) {
#ifdef EST_DEBUG
        /*if(u_degree[pos].vertex > u_num) {
            cout<<"u_degree[pos].vertex("<<u_degree[pos].vertex<< ") > u_num("<<u_num<<")"<<endl;
            exit(0);
        }*/
#endif
        u_pos[u_degree[pos].vertex] = pos; // get array of vertex mapping pos

        if(u_degree[pos].degree > ite_d) {
            uint before_d = ite_d;
            ite_d = u_degree[pos].degree;
            for(uint tt = before_d+1; tt <= ite_d; tt++)
                u_pointer[tt] = pos;
        }
    }
    u_pointer[ite_d+1] = u_degree_ed;
#ifdef EST_DEBUG
    cout<<"u_pointer:"<<endl;
    for(int i = 0; i <= u_max_degree+1; i++){
        cout<<u_pointer[i]<<",";
    }
    cout<<endl;
#endif


    ite_d = v_degree[v_degree_st].degree;
    v_pointer[ite_d] = v_degree_st;
    for(uint pos = v_degree_st; pos < v_degree_ed; pos++) {
#ifdef EST_DEBUG
        /*if(v_degree[pos].vertex > v_num) {
            cout<<"v_degree[pos].vertex("<<v_degree[pos].vertex<< ") > v_num("<<v_num<<")"<<endl;
            exit(0);
        }*/
#endif
        v_pos[v_degree[pos].vertex] = pos; // get array of vertex mapping pos

        if(v_degree[pos].degree > ite_d) {
            uint before_d = ite_d;
            ite_d = v_degree[pos].degree;
            for(uint tt = before_d+1; tt <= ite_d; tt++)
                v_pointer[tt] = pos;
        }
    }
    v_pointer[ite_d+1] = v_degree_ed;
#ifdef EST_DEBUG
    cout<<"v_pointer:"<<endl;
    for(int i = 0; i <= v_max_degree+1; i++){
        cout<<v_pointer[i]<<",";
    }
    cout<<endl;
#endif


    /// 5.peel + calculate cohesion
    /// O(|centerG|) <= O(d(uu)*d(vv))
    unsigned pp, qq, max_pp, max_qq;
    pp = core_u_num;
    //todo: get max_pp and max_qq (method check)
    ASSERT(v_degree_ed >= lower[u_row])
    max_pp = v_degree[v_degree_ed - lower[u_row]].degree;
    ASSERT(u_degree_ed >= lower[v_row])
    max_qq = u_degree[u_degree_ed - lower[v_row]].degree;
#ifdef EST_DEBUG
    cout<<"u_num="<<u_num<<", v_num="<<v_num<<endl;
    cout<<"max_pp="<<max_pp<<", max_qq="<<max_qq<<endl;
#endif
    //SortedVertexByDegree tmpnode(0,0);
    u_pointer[lower[u_row]-1] = 0;
    v_pointer[lower[v_row]-1] = 0;
    if(sizetype == 1){
        /// MEB
        //res = max(vv_nei.size(), uu_nei.size());
        //res = max(core_u_num * lower[u_row], core_v_num * lower[v_row]);
        res = max(max_pp * lower[u_row], max_qq * lower[v_row]); // todo: check
#ifdef EST_DEBUG
        cout<<"init res=max()="<<res<<endl;
#endif
        if(res >= edgen) return edgen;

        for(qq = lower[u_row] + 1; qq <= u_max_degree; qq++){
            // qq: u_degree, vertex number of Biclique_v
            // pp: vertex number of Core_u
#ifdef EST_DEBUG
            cout<<"qq="<<qq<<endl;
#endif
            bool finish = false;
            bool over = false;
            while(!finish) {
                finish = true;
                for(uint upos = u_pointer[qq-1]; upos < u_pointer[qq]; upos++) {
                    uint uid = u_degree[upos].vertex;
                    u_degree[upos].degree = 0;
                    u_pointer[qq-1] ++;

                    pp--; // delete u
                    finish = false;

                    for(auto vid: centerG.vertices[u_row][uid].nei){
                        uint vpos = v_pos[vid];
                        if(v_degree[vpos].degree >= lower[v_row]) {
                            uint vd = v_degree[vpos].degree;
                            v_degree[vpos].degree--;

                            // swap node in v_degree
                            uint sid = v_degree[v_pointer[vd]].vertex;
                            swap(v_pos[sid], v_pos[vid]);
                            SortedVertexByDegree tmpnode = v_degree[vpos];
                            v_degree[vpos] = v_degree[v_pointer[vd]];
                            v_degree[v_pointer[vd]] = tmpnode;
                            v_pointer[vd] ++;
                        }
                    }
                }
#ifdef EST_DEBUG
                cout<<"after peel u_degree="<<qq<<", v_pointer:"<<endl;
                for(int i = 0; i <= v_max_degree+1; i++){
                    cout<<v_pointer[i]<<",";
                }
                cout<<endl;
                cout<<"u_degree:"<<endl;
                for(int i = 0; i < u_degree_ed; i++){
                    cout<<u_degree[i].degree<<",";
                }
                cout<<endl;
                cout<<"v_degree:"<<endl;
                for(int i = 0; i < v_degree_ed; i++){
                    cout<<v_degree[i].degree<<",";
                }
                cout<<endl;
#endif

                /// pruning: p < \alpha
                if(pp < lower[v_row]) {
                    //return min(res, edgen);
#ifdef EST_DEBUG
                    cout<<"pp="<<pp<<" < lower[v_row]: res="<<res<<endl;
#endif
                    res = min(res, edgen);
                    over = true;
                    break;
                }

                for(uint vpos = v_pointer[lower[v_row]-1]; vpos < v_pointer[lower[v_row]]; vpos++) {
                    uint vid = v_degree[vpos].vertex;
                    v_degree[vpos].degree = 0;
                    v_pointer[lower[v_row]-1]++;

                    finish = false;

                    for(auto uid: centerG.vertices[v_row][vid].nei) {
                        uint upos = u_pos[uid];
                        if(u_degree[upos].degree >= qq) {
                            uint ud = u_degree[upos].degree;
                            //if(ud == qq) pp--; // delete u
                            u_degree[upos].degree--;

                            // swap node in u_degree
                            uint sid = u_degree[u_pointer[ud]].vertex;
                            swap(u_pos[sid], u_pos[uid]);
                            SortedVertexByDegree tmpnode = u_degree[upos];
                            u_degree[upos] = u_degree[u_pointer[ud]];
                            u_degree[u_pointer[ud]] = tmpnode;
                            u_pointer[ud] ++;
                        }
                    }
                }
#ifdef EST_DEBUG
                cout<<"after peel v_degree="<<lower[v_row]<<", u_pointer:"<<endl;
                for(int i = 0; i <= u_max_degree+1; i++){
                    cout<<u_pointer[i]<<",";
                }
                cout<<endl;
                cout<<"u_degree:"<<endl;
                for(int i = 0; i < u_degree_ed; i++){
                    cout<<u_degree[i].degree<<",";
                }
                cout<<endl;
                cout<<"v_degree:"<<endl;
                for(int i = 0; i < v_degree_ed; i++){
                    cout<<v_degree[i].degree<<",";
                }
                cout<<endl;
#endif

                //todo: check
                ASSERT(v_degree_ed >= qq)
                max_pp = v_degree[v_degree_ed - qq].degree;
                if(max_pp < lower[v_row]) {
                    //return min(res, edgen);
#ifdef EST_DEBUG
                    cout<<"max_pp="<<max_pp<<" < lower[v_row]: res="<<res<<endl;
#endif
                    res = min(res, edgen);
                    over = true;
                    break;
                }

            }
#ifdef EST_DEBUG
            cout<<"qq="<<qq<<", max_pp="<<max_pp<<", pp="<<pp<<endl;
#endif
            if(over) break;

            //res = max(res, (pp * qq));
            res = max(res, (max_pp * qq)); //todo: check method
#ifdef EST_DEBUG
            cout<<"res="<<res<<endl;
#endif

            if(res >= edgen) {
                //return edgen;
                res = edgen;
                break;
            }

        }

        //res = min(edgen, res);
    }
    else if(sizetype == 2){
        /// MVB
        res = max(max_pp + lower[u_row], max_qq + lower[v_row]); // todo: check
        if(res >= vertexn) return vertexn;

        for(qq = lower[u_row]+1; qq <= u_max_degree; qq++){
            // qq: u_degree, vertex number of Biclique_v
            // pp: vertex number of Core_u
            bool over = false;
            bool finish = false;
            while(!finish) {
                finish = true;
                for(uint upos = u_pointer[qq-1]; upos < u_pointer[qq]; upos++) {
                    uint uid = u_degree[upos].vertex;
                    u_degree[upos].degree = 0;
                    u_pointer[qq-1] ++;

                    pp--; // delete u
                    finish = false;

                    for(auto vid: centerG.vertices[u_row][uid].nei){
                        uint vpos = v_pos[vid];
                        if(v_degree[vpos].degree >= lower[v_row]) {
                            uint vd = v_degree[vpos].degree;
                            v_degree[vpos].degree--;

                            // swap node in v_degree
                            uint sid = v_degree[v_pointer[vd]].vertex; /// todo: revise in ooBCviz !!!
                            swap(v_pos[vid], v_pos[sid]);
                            SortedVertexByDegree tmpnode = v_degree[vpos];
                            v_degree[vpos] = v_degree[v_pointer[vd]];
                            v_degree[v_pointer[vd]] = tmpnode;
                            v_pointer[vd] ++;
                        }
                    }
                }

                /// pruning: p < \alpha
                if(pp < lower[v_row]) {
                    //return min(res, edgen);
                    //cout<<"pp < lower[v_row]: res="<<res<<endl;
                    res = min(res, vertexn);
                    over = true;
                    break;
                }

                for(uint vpos = v_pointer[lower[v_row]-1]; vpos < v_pointer[lower[v_row]]; vpos++) {
                    uint vid = v_degree[vpos].vertex;
                    v_degree[vpos].degree = 0;
                    v_pointer[lower[v_row]-1]++;

                    finish = false;

                    for(auto uid: centerG.vertices[v_row][vid].nei) {
                        uint upos = u_pos[uid];
                        if(u_degree[upos].degree >= qq) {
                            uint ud = u_degree[upos].degree;
                            //if(ud == qq) pp--; // delete u
                            u_degree[upos].degree--;

                            // swap node in u_degree
                            uint sid = u_degree[u_pointer[ud]].vertex; /// todo: revise in ooBCviz !!!
                            swap(u_pos[uid], u_pos[sid]);
                            SortedVertexByDegree tmpnode = u_degree[upos];
                            u_degree[upos] = u_degree[u_pointer[ud]];
                            u_degree[u_pointer[ud]] = tmpnode;
                            u_pointer[ud] ++;
                        }
                    }
                }

                //todo: check
                ASSERT(v_degree_ed >= qq)
                max_pp = v_degree[v_degree_ed - qq].degree;
                if(max_pp < lower[v_row]) {
                    //return min(res, edgen);
#ifdef EST_DEBUG
                    cout<<"max_pp="<<max_pp<<" < lower[v_row]: res="<<res<<endl;
#endif
                    res = min(res, edgen);
                    over = true;
                    break;
                }
            }

            if(over) break;

            //res = max(res, (pp + qq));
            res = max(res, (max_pp + qq)); //todo: check method
            //cout<<"res="<<res<<endl;

            if(res >= vertexn) {
                //return vertexn;
                res = vertexn;
                break;
            }

        }

    }
    else if(sizetype == 3){
        /// MBB
        //res = max(vv_nei.size(), uu_nei.size());
        //res = max(res, (qq * qq));
        res = max(lower[u_row] * lower[u_row], lower[v_row] * lower[v_row]);
        if(res >= edgen) return edgen;

        for(qq = lower[u_row]+1; qq <= u_max_degree; qq++){
            // qq: u_degree, vertex number of Biclique_v
            // pp: vertex number of Core_u
            bool over = false;
            bool finish = false;
            while(!finish) {
                finish = true;
                for(uint upos = u_pointer[qq-1]; upos < u_pointer[qq]; upos++) {
                    uint uid = u_degree[upos].vertex;
                    u_degree[upos].degree = 0;
                    u_pointer[qq-1] ++;

                    pp--; // delete u
                    finish = false;

                    for(auto vid: centerG.vertices[u_row][uid].nei){
                        uint vpos = v_pos[vid];
                        if(v_degree[vpos].degree >= lower[v_row]) {
                            uint vd = v_degree[vpos].degree;
                            v_degree[vpos].degree--;

                            // swap node in v_degree
                            uint sid = v_degree[v_pointer[vd]].vertex; /// todo: revise in ooBCviz !!!
                            swap(v_pos[vid], v_pos[sid]);
                            SortedVertexByDegree tmpnode = v_degree[vpos];
                            v_degree[vpos] = v_degree[v_pointer[vd]];
                            v_degree[v_pointer[vd]] = tmpnode;
                            v_pointer[vd] ++;
                        }
                    }
                }

                /// pruning: p < \alpha
                if(pp < lower[v_row] || pp < qq) {
                    //return min(res, edgen);
                    //cout<<"pp < lower[v_row]: res="<<res<<endl;
                    res = min(res, edgen);
                    over = true;
                    break;
                }

                for(uint vpos = v_pointer[lower[v_row]-1]; vpos < v_pointer[lower[v_row]]; vpos++) {
                    uint vid = v_degree[vpos].vertex;
                    v_degree[vpos].degree = 0;
                    v_pointer[lower[v_row]-1]++;

                    finish = false;

                    for(auto uid: centerG.vertices[v_row][vid].nei) {
                        uint upos = u_pos[uid];
                        if(u_degree[upos].degree >= qq) {
                            uint ud = u_degree[upos].degree;
                            //if(ud == qq) pp--; // delete u
                            u_degree[upos].degree--;

                            // swap node in u_degree
                            uint sid = u_degree[u_pointer[ud]].vertex; /// todo: revise in ooBCviz !!!
                            swap(u_pos[uid], u_pos[sid]);
                            SortedVertexByDegree tmpnode = u_degree[upos];
                            u_degree[upos] = u_degree[u_pointer[ud]];
                            u_degree[u_pointer[ud]] = tmpnode;
                            u_pointer[ud] ++;
                        }
                    }
                }

                //todo: check
                ASSERT(v_degree_ed >= qq)
                max_pp = v_degree[v_degree_ed - qq].degree;
                if(max_pp < lower[v_row] || max_pp < qq) {
                    //return min(res, edgen);
#ifdef EST_DEBUG
                    cout<<"max_pp="<<max_pp<<" < lower[v_row]: res="<<res<<endl;
#endif
                    res = min(res, edgen);
                    over = true;
                    break;
                }
            }

            if(over) break;

            res = max(res, (qq * qq));
            //cout<<"res="<<res<<endl;

            if(res >= edgen) {
                //return edgen;
                res = edgen;
                break;
            }

        }
    }

#ifdef DEBUG
    cout<<"return edge.res="<<res<<endl;
#endif // DEBUG

    delete[] u_degree;
    delete[] v_degree;
    delete[] u_pos;
    delete[] v_pos;

#ifdef EST_DEBUG
    cout<<"calculate cohesion estimation is ok!"<<endl;
    cout<<"res="<<res<<endl<<endl<<endl;
#endif

    return res;
}

/// exact O(R2^(dmax))
unsigned BCVIZ::BCviz_EdgeCohesivenessExact(const unsigned uu, const unsigned vv, const unsigned uu_row)
{
    unsigned res, edgen, vertexn, u_num, v_num;
    unsigned lower[2]={t_V, t_U};
    vector<Vertex>& urow_vertices = bigraph.vertices[uu_row];
    vector<unsigned>& uu_nei = urow_vertices[uu].nei;
    unsigned vv_row = 1-uu_row;
    vector<Vertex>& vrow_vertices = bigraph.vertices[vv_row];
    vector<unsigned>& vv_nei = vrow_vertices[vv].nei;
    res = edgen = vertexn = u_num = v_num = 0;
#ifdef EST_DEBUG
    cout<<"uu_row="<<uu_row<<", uu_id="<<bigraph.vertices[uu_row][uu].id<<", vv_id=";
    cout<<bigraph.vertices[vv_row][vv].id<<endl;
#endif

    /// special situation1: vex has just one degree.
    if(uu_nei.size() < lower[uu_row] || vv_nei.size() < lower[vv_row]) return (unsigned)0;


    /// 1.get center graph O(m)
    Bigraph centerG;
    unsigned max_v = *max_element(uu_nei.begin(), uu_nei.end());
    unsigned max_u = *max_element(vv_nei.begin(), vv_nei.end());
    vector<unsigned> u_map_id(max_u+1, 0);
    vector<unsigned> v_map_id(max_v+1, 0);
    for(auto ui: vv_nei) u_map_id[ui] = ++u_num;
    for(auto vj: uu_nei) v_map_id[vj] = ++v_num;
    centerG.vertices[uu_row].resize(u_num+1);
    centerG.vertices[vv_row].resize(v_num+1);
    for(auto ui: vv_nei) {
        for(auto vj: bigraph.vertices[uu_row][ui].nei) {
            if(ui > max_u || vj > max_v) continue; /// O(1) map
            uint utmp = u_map_id[ui];
            uint vtmp = v_map_id[vj];
            if(utmp==0 || vtmp==0) continue;
            centerG.vertices[uu_row][ utmp ].nei.emplace_back( vtmp );
            //centerG.vertices[uu_row][ utmp ].id = ui;
            centerG.vertices[vv_row][ vtmp ].nei.emplace_back( utmp );
            //centerG.vertices[vv_row][ vtmp ].id = vj;

            //edgen++;
        }
    }

#ifdef EST_DEBUG
    cout<<"get center graph is ok!"<<endl;
    //printGraph(centerG);
    /*cout<<"End of getting center graph"<<endl;
    cout << "#vertices = " << u_num + v_num << endl;
    cout << "#edges = " << edgen << endl;*/
#endif

    /// 2.core Reduce: get (\beta, \alpha)-core O(m)
    reduce1hop(centerG, lower[vv_row], lower[uu_row]);
#ifdef EST_DEBUG
    cout<<"core reduce is ok!"<<endl;
    //printGraph(centerG);
#endif

    /// 3. swap O(n)
    uint core_u_num, core_v_num;
    core_u_num = core_v_num = 0;
    for(auto &uit : centerG.vertices[uu_row]) /// pos from 1
        if(!uit.nei.empty()) core_u_num++;
    for(auto &vit : centerG.vertices[vv_row])
        if(!vit.nei.empty()) core_v_num++;
    if(core_v_num < lower[uu_row] || core_u_num < lower[vv_row]) return (unsigned)0;

    unsigned u_row, v_row;
    u_row = uu_row;
    v_row = vv_row;
    if(core_u_num > core_v_num) {
        u_row = vv_row;
        v_row = uu_row;

        uint tmp = core_u_num;
        core_u_num = core_v_num;
        core_v_num = tmp;

        tmp = u_num; // todo: new code
        u_num = v_num;
        v_num = tmp;

#ifdef EST_DEBUG
        cout<<"swap U and V"<<endl;
#endif
    }

#ifdef EST_DEBUG
    cout<<"core_u_num="<<core_u_num<<", core_v_num="<<core_v_num<<endl;
#endif
    for(auto &uit : centerG.vertices[u_row]){
        edgen += uit.nei.size();
    }

    double density = edgen*1.0/core_u_num/core_v_num;
    if(density == 1) {
        if(pro_name=="MEB") return core_u_num * core_v_num;
        if(pro_name=="MVB") return core_u_num + core_v_num;
        if(pro_name=="MBB") return min(core_u_num * core_u_num, core_v_num * core_v_num);
    }

    /// 5. compute MBsize exactly
    vector<vector<uint>> L, R;
    uint u_max_degree = 0;
    uint v_max_degree = 0;
    for(auto & i : centerG.vertices[u_row]) {
        L.emplace_back(i.nei);
        if(u_max_degree < i.nei.size()) u_max_degree = i.nei.size();
    }
    for(auto & i : centerG.vertices[v_row]) {
        R.emplace_back(i.nei);
        if(v_max_degree < i.nei.size()) v_max_degree = i.nei.size();
    }
    bgraph tmpG(edgen, u_max_degree, v_max_degree, L, R);

    MB exactMB(tmpG, lower[v_row], lower[u_row], "exact", pro_name);
    res = exactMB.MBC_improved();

#ifdef EST_DEBUG
    cout<<"calculate cohesion estimation is ok!"<<endl;
    cout<<"res="<<res<<endl<<endl<<endl;
#endif

    return res;
}


/**exact ordering framework
 * time: O(MlogN + Mm)*/
void BCVIZ::BCviz1_Construct() {
    ///cout<<"####### BCviz1_Construct #######"<<endl;
    //priority_queue<SortedBySize_BCviz, vector<SortedBySize_BCviz>, cmp> que;
    set<SortedBySize_BCviz, cmp> heap;

    unsigned row = 0;
    unsigned iipos = 0;
    if(vec_tsize[1].size() < vec_tsize[0].size()) row = 1;
    unsigned allnum = bigraph.vertices[row].size();

    while(iipos < allnum){
        /**optimization1: choose the min size vex for the first vex is better*/
        unsigned ii=iipos;
        for(; ii < vec_tsize[row].size(); ii++){
            if(!vec_vis[row][ii]){
                SortedBySize_BCviz stvex(row, ii, 0, 0);
                //que.push(stvex);
                heap.insert(stvex);
                vec_vis[row][ii] = true;
                iipos = ii + 1;
                break;
            }
        }
        if(ii == allnum) break;

        vector<SortedBySize_BCviz> node_ordered_vex;
        while(!heap.empty()){
            SortedBySize_BCviz topvex = *heap.begin();
            unsigned urow = topvex.row;
            unsigned ui = topvex.vex;
            unsigned usize = topvex.twoNodeSize;
#ifdef DEBUG
            cout<<"row="<<urow<<", vex="<<bigraph.vertices[urow][ui].id<<", oneNodeSize="<<topvex.oneNodeSize<<", twoNodeSize="<<topvex.twoNodeSize<<endl;
#endif // DEBUG

            node_ordered_vex.emplace_back(topvex);
            is_in_O[urow][ui] = true;
            //que.pop();
            heap.erase(*heap.begin());

            ///1. state[v] = ORDERED -> continue
            ///2. state[v] = UNPROCESSED -> calc cohesion
            ///3. state[v] = PENDING -> update cohesion, reorder heap
            ASSERT(ui < bigraph.vertices[urow].size())
            vector<unsigned>& ui_nei = bigraph.vertices[urow][ui].nei;
            unsigned vrow = 1 - urow;
#ifdef DEBUG
            cout<<"start check vex neibor"<<endl;
#endif // DEBUG
            for(auto &vj: ui_nei){
                ASSERT(vj < vec_vis[vrow].size())
                if(vec_vis[vrow][vj]){
                    if(is_in_O[vrow][vj]) continue;
                    //unsigned os = vec_size[vrow][vj];
                    unsigned ts = vec_tsize[vrow][vj];
#ifdef DEBUG
                    cout<<"update: row="<<vrow<<", vex="<<bigraph.vertices[vrow][vj].id<<": ";
#endif // DEBUG
                    unsigned newts = 0;
                    newts = BCviz_EdgeCohesivenessShrinkGraph(ui, vj, urow);


                    if(ts < newts) {
                        SortedBySize_BCviz tmpvjnode(vrow, vj, 0, ts);
                        auto iter = heap.find(tmpvjnode); //log(n) --> O(n) use Q_visit
                        heap.erase(*iter);

                        vec_tsize[vrow][vj] = newts;
                        SortedBySize_BCviz newvjnode(vrow, vj, 0, newts);
                        heap.insert(newvjnode);
                    }
#ifdef DEBUG
                    //cout<<"oneNodeSize="<<vec_size[vrow][vj]<<", twoNodeSize="<<ts<<endl;
                    cout<<"update twoNodeSize="<<ts<<endl;
                    //print_heap(heap);
#endif // DEBUG
                }
                else{
                    vec_vis[vrow][vj] = true;
#ifdef DEBUG
                    cout<<"add ui_nei: row="<<vrow<<", vex="<<bigraph.vertices[vrow][vj].id<<": ";
#endif // DEBUG
                    unsigned ts = 0;
                    ts = BCviz_EdgeCohesivenessShrinkGraph(ui, vj, urow);
                    vec_tsize[vrow][vj] = ts;
                    SortedBySize_BCviz vjnode(vrow, vj, 0, ts);
                    heap.insert(vjnode);
#ifdef DEBUG
                    //cout<<"oneNodeSize="<<vec_size[vrow][vj]<<", twoNodeSize="<<ts<<endl;
                    cout<<"calc twoNodeSize="<<ts<<endl;
#endif // DEBUG
                }
            }
#ifdef DEBUG
            cout<<"End"<<endl;
            //cout<<" of row="<<urow<<", vex="<<bigraph.vertices[urow][ui].id<<", oneNodeSize="<<topvex.oneNodeSize<<", twoNodeSize="<<topvex.twoNodeSize<<endl;
            //cout<<"calculating. "<<endl<<endl;
#endif // DEBUG
        }
        if(node_ordered_vex.size() >= t_U * t_V)
            ordered_vex.emplace_back(node_ordered_vex);
    }//while visnum < allnum

}

/**approximate ordering framework
 * time: O(NlogN + Nm + Tupdate)*/
void BCVIZ::BCviz2_Construct() {
    ///cout<<"####### BCviz2_Construct #######"<<endl;
    //priority_queue<SortedBySize_BCviz, vector<SortedBySize_BCviz>, cmp> que;
    set<SortedBySize_BCviz, cmp> heap;

    unsigned row = 0;
    unsigned iipos = 0;
    if(vec_tsize[1].size() < vec_tsize[0].size()) row = 1;
    unsigned allnum = bigraph.vertices[row].size();
    while(iipos < allnum){
        /**optimization1: choose the min size vex for the first vex is better*/
        unsigned ii=iipos;
        for(; ii < vec_tsize[row].size(); ii++){
            if(!vec_vis[row][ii]){
                SortedBySize_BCviz stvex(row, ii, 0, 0);
                //que.push(stvex);
                heap.insert(stvex);
                vec_vis[row][ii] = true;
                iipos = ii + 1;
                break;
            }
        }
        if(ii == allnum) break;

        vector<SortedBySize_BCviz> node_ordered_vex;
        uint last_tsize = 0; /// todo:climbing update 1
        while(!heap.empty()){
            SortedBySize_BCviz topvex = *heap.begin();
            unsigned urow = topvex.row;
            unsigned ui = topvex.vex;
            unsigned usize = topvex.twoNodeSize;
#ifdef DEBUG
            cout<<"row="<<urow<<", vex="<<bigraph.vertices[urow][ui].id<<", oneNodeSize="<<topvex.oneNodeSize<<", twoNodeSize="<<topvex.twoNodeSize<<endl;
#endif // DEBUG

            node_ordered_vex.emplace_back(topvex);
            is_in_O[urow][ui] = true;
            //que.pop();
            heap.erase(*heap.begin());

            ///1. state[v] = ORDERED -> continue
            ///2. state[v] = UNPROCESSED -> calc cohesion
            ///3. state[v] = PENDING -> update cohesion, reorder heap
            ASSERT(ui < bigraph.vertices[urow].size())
            vector<unsigned>& ui_nei = bigraph.vertices[urow][ui].nei;
            unsigned vrow = 1 - urow;
            bool update = false; /// todo:climbing update 2
#ifdef DEBUG
            cout<<"start check vex neibor"<<endl;
#endif // DEBUG
            for(auto &vj: ui_nei){
                ASSERT(vj < vec_vis[vrow].size())
                if(!vec_vis[vrow][vj]){
                    vec_vis[vrow][vj] = true;
#ifdef DEBUG
                    cout<<"add ui_nei: row="<<vrow<<", vex="<<bigraph.vertices[vrow][vj].id<<": ";
#endif // DEBUG
                    unsigned ts = 0;
                    ts = BCviz_EdgeCohesivenessShrinkGraph(ui, vj, urow);
                    vec_tsize[vrow][vj] = ts;
                    SortedBySize_BCviz vjnode(vrow, vj, 0, ts);
                    heap.insert(vjnode);

                    if(ts > usize) update = true; /// todo:climbing update 2
#ifdef DEBUG
                    //cout<<"oneNodeSize="<<vec_size[vrow][vj]<<", twoNodeSize="<<ts<<endl;
                    cout<<"twoNodeSize="<<ts<<endl;
#endif // DEBUG
                }
                else {
                    /// todo:climbing update 1
                    if(usize > last_tsize && (!is_in_O[vrow][vj])) {
                        unsigned ts = vec_tsize[vrow][vj];
                        unsigned newts = 0;
#ifdef DEBUG
                        cout<<"update: row="<<vrow<<", vex="<<bigraph.vertices[vrow][vj].id<<": ";
#endif // DEBUG

                        newts = BCviz_EdgeCohesivenessShrinkGraph(ui, vj, urow);
                        if(ts < newts) {
                            SortedBySize_BCviz tmpvjnode(vrow, vj, 0, ts);
                            auto iter = heap.find(tmpvjnode); //log(n) --> O(n) use Q_visit
                            if(iter!=heap.end()) heap.erase(*iter);

                            //ts = newts;
                            vec_tsize[vrow][vj] = newts;
                            SortedBySize_BCviz newvjnode(vrow, vj, 0, newts);
                            heap.insert(newvjnode);
                        }

#ifdef DEBUG
                        //cout<<"oneNodeSize="<<vec_size[vrow][vj]<<", twoNodeSize="<<ts<<endl;
                        cout<<"update twoNodeSize="<<ts<<endl;
                        //print_heap(heap);
#endif // DEBUG
                    }
                }


            }

            /// todo:climbing update 2
            if(update && usize <= last_tsize) {
                for(auto &vj: ui_nei) {
                    if(vec_vis[vrow][vj] && (!is_in_O[vrow][vj])) {
                        unsigned ts = vec_tsize[vrow][vj];
                        unsigned newts = 0;
                        newts = BCviz_EdgeCohesivenessShrinkGraph(ui, vj, urow);
                        if(ts < newts) {
                            SortedBySize_BCviz tmpvjnode(vrow, vj, 0, ts);
                            auto iter = heap.find(tmpvjnode); //log(n) --> O(n) use Q_visit
                            if(iter!=heap.end()) heap.erase(*iter);

                            //ts = newts;
                            vec_tsize[vrow][vj] = newts;
                            SortedBySize_BCviz newvjnode(vrow, vj, 0, newts);
                            heap.insert(newvjnode);
                        }
                    }
                }
            }

            /// todo:climbing update 1
            last_tsize = usize;

#ifdef DEBUG
            cout<<"End"<<endl;
            //cout<<" of row="<<urow<<", vex="<<bigraph.vertices[urow][ui].id<<", oneNodeSize="<<topvex.oneNodeSize<<", twoNodeSize="<<topvex.twoNodeSize<<endl;
            //cout<<"calculating. "<<endl<<endl;
#endif // DEBUG
        }
        if(node_ordered_vex.size() >= t_U * t_V)
            ordered_vex.emplace_back(node_ordered_vex);
    }//while visnum < allnum

}

/**exact ordering framework with offline cohesion
 * time: O(MlogN + T_MB)*/
void BCVIZ::BCviz_Construct() {
    ///cout<<"####### BCviz_Construct #######"<<endl;
    //priority_queue<SortedBySize_BCviz, vector<SortedBySize_BCviz>, cmp> que;
    set<SortedBySize_BCviz, cmp> heap;

    unsigned row = 0;
    unsigned iipos = 0;
    if(vec_tsize[1].size() < vec_tsize[0].size()) row = 1;
    unsigned allnum = bigraph.vertices[row].size();
    pair<uint, uint> e;
    //unsigned id_u, id_v;

    while(iipos < allnum){
        /**optimization1: choose the min size vex for the first vex is better*/
        unsigned ii=iipos;
        for(; ii < vec_tsize[row].size(); ii++){
            if(!vec_vis[row][ii]){
                SortedBySize_BCviz stvex(row, ii, 0, 0);
                //que.push(stvex);
                heap.insert(stvex);
                vec_vis[row][ii] = true;
                iipos = ii + 1;
                break;
            }
        }
        if(ii == allnum) break;

        vector<SortedBySize_BCviz> node_ordered_vex;
        while(!heap.empty()){
            SortedBySize_BCviz topvex = *heap.begin();

            unsigned urow = topvex.row;
            unsigned ui = topvex.vex;
            //unsigned usize = topvex.twoNodeSize;
            //unsigned uosize = topvex.oneNodeSize;
#ifdef DEBUG
            cout<<"row="<<urow<<", vex="<<bigraph.vertices[urow][ui].id<<", oneNodeSize="<<topvex.oneNodeSize<<", twoNodeSize="<<topvex.twoNodeSize<<endl;
#endif // DEBUG
            node_ordered_vex.emplace_back(topvex);
            is_in_O[urow][ui] = true;
            //que.pop();
            heap.erase(*heap.begin());

            ///1. state[v] = ORDERED -> continue
            ///2. state[v] = UNPROCESSED -> calc cohesion
            ///3. state[v] = PENDING -> update cohesion, reorder heap
            ASSERT(ui < bigraph.vertices[urow].size())
            vector<unsigned>& ui_nei = bigraph.vertices[urow][ui].nei;
            unsigned vrow = 1 - urow;
#ifdef DEBUG
            cout<<"start check vex neibor"<<endl;
#endif // DEBUG
            for(auto &vj: ui_nei){
                ASSERT(vj < vec_vis[vrow].size())
                if(vec_vis[vrow][vj]){
                    if(is_in_O[vrow][vj]) continue;
                    unsigned ts = vec_tsize[vrow][vj];
#ifdef DEBUG
                        cout<<"update: row="<<vrow<<", vex="<<bigraph.vertices[vrow][vj].id<<": ";
#endif // DEBUG
                    ulong newts = 0;
                    //id_u = bigraph.vertices[urow][ui].id;
                    //id_v = bigraph.vertices[vrow][vj].id;
                    if(urow == 0) {e = make_pair(ui, vj);}
                    else {e = make_pair(vj, ui);}
                    ASSERT(edgeCoh.count(e))
                    newts = edgeCoh[e];
                    if(ts < newts) {
                        SortedBySize_BCviz tmpvjnode(vrow, vj, 0, ts);
                        auto iter = heap.find(tmpvjnode); //log(n) --> O(n) use Q_visit
                        heap.erase(*iter);

                        vec_tsize[vrow][vj] = newts;
                        SortedBySize_BCviz newvjnode(vrow, vj, 0, newts);
                        heap.insert(newvjnode);
                    }
#ifdef DEBUG
                    //cout<<"oneNodeSize="<<vec_size[vrow][vj]<<", twoNodeSize="<<ts<<endl;
                    cout<<"update twoNodeSize="<<ts<<endl;
                    //print_heap(heap);
#endif // DEBUG
                }
                else{
#ifdef DEBUG
                    cout<<"add ui_nei: row="<<vrow<<", vex="<<bigraph.vertices[vrow][vj].id<<": ";
#endif // DEBUG
                    vec_vis[vrow][vj] = true;
                    ulong ts = 0;
                    //id_u = bigraph.vertices[urow][ui].id;
                    //id_v = bigraph.vertices[vrow][vj].id;
                    if(urow == 0) {e = make_pair(ui, vj);}
                    else {e = make_pair(vj, ui);}
                    ASSERT(edgeCoh.count(e))
                    ts = edgeCoh[e];

                    vec_tsize[vrow][vj] = ts;
                    SortedBySize_BCviz vjnode(vrow, vj, 0, ts);
                    heap.insert(vjnode);
#ifdef DEBUG
                    //cout<<"oneNodeSize="<<vec_size[vrow][vj]<<", twoNodeSize="<<ts<<endl;
                    cout<<"calc twoNodeSize="<<ts<<endl;
#endif // DEBUG
                }
            }
#ifdef DEBUG
            cout<<"End"<<endl;
            //cout<<" of row="<<urow<<", vex="<<bigraph.vertices[urow][ui].id<<", oneNodeSize="<<topvex.oneNodeSize<<", twoNodeSize="<<topvex.twoNodeSize<<endl;
            //cout<<"calculating. "<<endl<<endl;
#endif // DEBUG
        }
        if(node_ordered_vex.size() >= t_U * t_V)
            ordered_vex.emplace_back(node_ordered_vex);
    }//while visnum < allnum

}

/**
    tool
*/
/*
long MSIF::get_memory_used(){
    int who = RUSAGE_SELF;
    struct rusage usage;
    getrusage(who, &usage);
    return usage.ru_maxrss;
}
*/

void BCVIZ::printVec(vector<uint> &vec) {
    for(auto it: vec) {
        cout<<it<<" ";
    }
    cout<<endl;
}

void BCVIZ::printGraph(Bigraph &PG) {
    for(int i=0; i<ROW_SIZE; i++) {
        for(int vit=0; vit < PG.vertices[i].size(); vit++) {
            cout<<vit<<": ";
            printVec(PG.vertices[i][vit].nei);
        }
        cout<<endl;
    }

}

void BCVIZ::printIDGraph(Bigraph &PG) {
    for(int i=0; i<ROW_SIZE; i++) {
        for(auto vit: PG.vertices[i]) {
            cout<<vit.id<<": ";
            printVec(vit.nei);
        }
        cout<<endl;
    }

}

void BCVIZ::BCviz_out_result(const string& tail) {
    string resultfile = "../Index-results/" + dataset + "_" + pro_name + "_" + tail + ".txt";
    //string resultfile = "D:/Source Code/BCviz/Index-results/" + dataset + "_" + pro_name + "_" + tail + ".txt";

//cout<<"store the result in path: "<<resultfile<<endl;

    FILE* fp = nullptr;
    fp = freopen(resultfile.c_str(), "w", stdout);
    if(!fp){
        cout<<"can't open XXX_"<<tail<<".txt !"<<endl;
        return ;
    }

    // save vertex and cohesion
/*    cout<<ordered_vex.size()<<endl;
    for(int i = 0; i < ordered_vex.size(); i++){
        cout<<ordered_vex[i].size()<<endl;
        for(auto &it: ordered_vex[i]){
            cout<<it.row<<" "<<bigraph.vertices[it.row][it.vex].id<<" "<<it.twoNodeSize<<endl;
        }
    }
    fclose(stdout);*/

    // reduce some useless vertices
    int Cnum = 0;
    unsigned zero = 0;
    for(uint i = 0; i < ordered_vex.size(); i++) {
        int num = 0;
        int j = 0;
        for (j = 0; j < ordered_vex[i].size() - 1; j++) {
            if (ordered_vex[i][j].twoNodeSize == zero && ordered_vex[i][j + 1].twoNodeSize == zero) { continue; }
            num++;
        }
        if (ordered_vex[i][j].twoNodeSize == zero) {}
        else num++;

        if ((sizetype == 1 || sizetype == 3) && num < t_U * t_V) continue;
        if ((sizetype == 2) && num < t_U + t_V) continue;

        Cnum++;
    }

    cout<<Cnum<<endl;
    for(int i = 0; i < ordered_vex.size(); i++) {
        int num = 0;
        int j = 0;
        for (j = 0; j < ordered_vex[i].size() - 1; j++) {
            if (ordered_vex[i][j].twoNodeSize == zero && ordered_vex[i][j + 1].twoNodeSize == zero) { continue; }
            num++;
        }
        if (ordered_vex[i][j].twoNodeSize == zero) {}
        else num++;

        if ((sizetype == 1 || sizetype == 3) && num < t_U * t_V) continue;
        if ((sizetype == 2) && num < t_U + t_V) continue;

        cout<<num<<endl;
        for(j = 0; j < ordered_vex[i].size() - 1; j++){
            if(ordered_vex[i][j].twoNodeSize == zero && ordered_vex[i][j+1].twoNodeSize == zero) {continue;}
            auto it = ordered_vex[i][j];
            cout<<it.row<<" "<<bigraph.vertices[it.row][it.vex].id<<" "<<it.twoNodeSize<<endl;
        }

        if(ordered_vex[i][j].twoNodeSize == zero) {}
        else {
            auto it = ordered_vex[i][j];
            cout<<it.row<<" "<<bigraph.vertices[it.row][it.vex].id<<" "<<it.twoNodeSize<<endl;
        }
    }

    fclose(stdout);

    //cout<<"storing the result is ok!\n\n";
}

