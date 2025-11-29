#include "MB.h"

bgraph::~bgraph() {
    for(auto itvec: adjL) {
        itvec.clear();
    }
    adjL.clear();
    for(auto itvec: adjR) {
        itvec.clear();
    }
    adjR.clear();
}

void bgraph::sortAdj() {
    for(auto &list:adjL){
        sort(list.begin(),list.end());
    }
    for(auto &list:adjR){
        sort(list.begin(),list.end());
    }
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
            int m=0;
            for(auto& neighbours:adjR){
                m+=neighbours.size();
            }
            noE = m;
        }

    }
}

void MB::iMBEA_Gk(bgraph &Gk,
                     uint &tkU,
                     uint &tkV,
                     vector<uint> &U,
                     vector<uint> &V,
                     vector<uint> &CV,
                     vector<uint> &XV) {
    // TODO: assume nei_id is in increasing order
    if(U.empty()) {return;}

    /// update maximum edge Biclique
    if((U.size() >= tkU) && (V.size() >= tkV)) {
        if(problem_type == "MEB" && V.size() * U.size() > _C.MBsize) {
            //_C.subU = U;
            //_C.subV = V;
            _C.MBsize = V.size() * U.size();
        }
        else if (problem_type == "MVB" && V.size() + U.size() > _C.MBsize) {
            //_C.subU = U;
            //_C.subV = V;
            _C.MBsize = V.size() + U.size();
        }
        else if (problem_type == "MBB" && min(V.size(),U.size())*min(V.size(),U.size()) > _C.MBsize) {
            //_C.subU = U;
            //_C.subV = V;
            vector<uint> tmpU = U.size() < V.size() ? U : V;
            _C.MBsize = tmpU.size() * tmpU.size();
        }
    }

    /// update _U, _V, _CV, _XV, where select v in CV and put it into V
    vector<uint> _U, _V, _CV, _XV;
    const uint CV_size = CV.size();
    bool reduce_CV_pos[CV_size];
    const uint U_size = U[U.size()-1]; /// assume nei is in increasing order

    bool *U_hash = new bool[U_size+1];
    for(uint i = 0; i <= U_size; i++) U_hash[i] = false;
    for(uint u : U) {
        //ASSERT(u <= U_size)
        U_hash[u] = true;
    }

    for(uint i = 0; i < CV.size(); i++){
        uint v = CV[i];

        /// 1. update _U; O(degree_R_max)
        _U.clear();
        for(auto u : Gk.adjR[v]) {
            if(u > U_size) break;
            if(U_hash[u])
                _U.emplace_back(u);
        }


        // TODO-pruning: don't exist a vertex xv in XV, _U is a subset of N(v,G)
        /// assume nei is in increasing order
        uint us = _U.size();
        bool is_not_exist = true;
        /// function:  _U is a subset of N(v,G) where v is in XV
        /// O(|XV|*(degree_V_max + num_xv_nei)) <= O(|XV|*degree_V_max + |E+|)
        for(auto xv : XV) {
            vector<uint> &nxv = Gk.adjR[xv];
            uint num_xv_nei = nxv.size();
            int pos_v = 0;
            if(num_xv_nei >= us) {
                bool is_subset = true;
                // O(degree_R_max + degree_R_max)
                for(int pos_u = 0; pos_u < _U.size(); pos_u++) {
                    uint unode = _U[pos_u];
                    while(pos_v < num_xv_nei && nxv[pos_v] < unode) pos_v++;
                    if(pos_v >= num_xv_nei || (pos_v < num_xv_nei && nxv[pos_v] > unode)) {
                        is_subset = false;
                        break;
                    }
                }
                if(is_subset) {
                    is_not_exist = false;
                    break;
                }
            }
        }

        /// 2.update _V, _CV, _XV;
        /// assume nei is in increasing order
        if(is_not_exist){
            /// update _V <=O(|E+|*degree_V_max)
            /// function:  find all vertex v in subCV that _U is a subset of N(v,G)
            ///            compute subCV and put subCV to _V
            /// O(|CV|*(|_U|+num_nei_cv)) <= O(|CV|*|_U|+|E-|) |_U| = degree_V_max
            _V.clear();
            _V = V;
            _V.emplace_back(v);
            memset(reduce_CV_pos, false, sizeof(reduce_CV_pos));
            for(uint j = i+1; j < CV.size(); j++) { ///test
                //if(j == i) continue;
                auto cv = CV[j];
                vector<uint> &ncv = Gk.adjR[cv];
                uint num_cv_nei = ncv.size();
                int pos_v = 0;
                if(num_cv_nei >= us) {
                    bool is_subset = true;
                    // O(|_U|*num_nei_cv)
                    for(int pos_u = 0; pos_u < _U.size(); pos_u++){
                        uint unode = _U[pos_u];
                        while(pos_v < num_cv_nei && ncv[pos_v] < unode){
                            pos_v++;
                        }
                        if(pos_v >= num_cv_nei || (pos_v < num_cv_nei && ncv[pos_v] > unode)) {
                            is_subset = false;
                            break;
                        }
                    }
                    if(is_subset){
                        _V.emplace_back(cv);
                        reduce_CV_pos[j] = true;
                    }

                }
            }

            /// update _CV O(|CV|*|_U|*|num_cv_nei|) = O(|E-|*degree_V_max)
            _CV.clear();
            for(uint j = i+1; j < CV.size(); j++) {
                if(reduce_CV_pos[j]) continue;
                auto cv = CV[j];
                vector<uint> &ncv = Gk.adjR[cv];
                uint num_cv_nei = ncv.size();
                uint num_joint = 0;
                int pos_u = 0;
                int pos_v = 0;
                // O(|_U|*|num_cv_nei|)
                while(pos_u < us && pos_v < num_cv_nei) {
                    while(pos_u < us && pos_v < num_cv_nei && _U[pos_u] < ncv[pos_v]) pos_u++;
                    while(pos_u < us && pos_v < num_cv_nei && _U[pos_u] > ncv[pos_v]) pos_v++;
                    while(pos_u < us && pos_v < num_cv_nei && _U[pos_u] == ncv[pos_v]) {
                        num_joint++;
                        pos_u++;
                        pos_v++;
                    }
                }
                if(num_joint >= tkU) {
                    _CV.emplace_back(cv);
                }
            }

            /// update _XV O(|XV|*|_U|*|num_xv_nei|)=O(|E-|*degree_V_max)
            _XV.clear();
            for(unsigned int xv : XV) {
                vector<uint> &nxv = Gk.adjR[xv];
                uint num_xv_nei = nxv.size();
                uint num_joint = 0;
                int pos_u = 0;
                int pos_v = 0;
                while(pos_u < us && pos_v < num_xv_nei) {
                    while(pos_u < us && pos_v < num_xv_nei && _U[pos_u] < nxv[pos_v]) pos_u++;
                    while(pos_u < us && pos_v < num_xv_nei && _U[pos_u] > nxv[pos_v]) pos_v++;
                    while(pos_u < us && pos_v < num_xv_nei && _U[pos_u] == nxv[pos_v]){
                        num_joint++;
                        pos_u++;
                        pos_v++;
                    }
                }
                if(num_joint >= tkU) {
                    _XV.emplace_back(xv);
                }
            }

            /// branch and bound: recall iMBEA itself
            uint vs = _V.size();
            uint cvs = _CV.size();
            if((us >= tkU) && (vs + cvs >= tkV) ) {
                if(problem_type == "MEB" && (us * (vs+cvs) > _C.MBsize))
                    iMBEA_Gk(Gk, tkU, tkV, _U, _V, _CV, _XV);
                if(problem_type == "MVB" && (us + (vs+cvs) > _C.MBsize))
                    iMBEA_Gk(Gk, tkU, tkV, _U, _V, _CV, _XV);
                if(problem_type == "MBB" && (min(us,(vs+cvs))*min(us,(vs+cvs)) > _C.MBsize))
                    iMBEA_Gk(Gk, tkU, tkV, _U, _V, _CV, _XV);
            }

        }

        /// record N(v,G) result in XV
        XV.emplace_back(v);

    }

    delete[] U_hash;

}


void MB::MBC_Gk(bgraph &Gk, uint tkU, uint tkV) {
    if(Gk.noE == 0) return;

    // initial U, V, CV, XV
    vector<uint> U, V, CV, XV;
    for(uint i = 0; i < Gk.adjL.size(); i++) {
        if(Gk.adjL[i].size() >= tkV)
            U.emplace_back(i);
    }
    for(uint i = 0; i < Gk.adjR.size(); i++) {
        if(Gk.adjR[i].size() >= tkU)
            CV.emplace_back(i);
    }

    if(U.size() >= tkU && CV.size() >= tkV)
        iMBEA_Gk(Gk, tkU, tkV, U, V, CV, XV);

    U.clear();
    V.clear();
    CV.clear();
    XV.clear();

}

void MB::iMBEA_init(bgraph &Gk,
                       vector<uint> &U,
                       vector<uint> &V,
                       vector<uint> &CV,
                       vector<uint> &XV) {
    if(U.empty() || init_MEB) {return;}

    /// update maximum edge Biclique
    if(U.size() >= t_U && V.size() >= t_V) {
        if(problem_type == "MEB" && V.size() * U.size() > _C.MBsize) {
            //_C.subU = U;
            //_C.subV = V;
            _C.MBsize = V.size() * U.size();
        }
        else if (problem_type == "MVB" && V.size() + U.size() > _C.MBsize) {
            //_C.subU = U;
            //_C.subV = V;
            _C.MBsize = V.size() + U.size();
        }
        else if (problem_type == "MBB" && min(V.size(),U.size())*min(V.size(),U.size()) > _C.MBsize) {
            //_C.subU = U;
            //_C.subV = V;
            vector<uint> tmpU = U.size() < V.size() ? U : V;
            _C.MBsize = tmpU.size() * tmpU.size();
        }
        init_MEB = true;
        return;
        //printMEB(_C);
    }

    /// update _U, _V, _CV, _XV, where select v in CV and put it into V
    vector<uint> _U, _V, _CV, _XV;
    const uint CV_size = CV.size();
    bool *reduce_CV_pos = new bool[CV_size];

    const uint U_size = U[U.size()-1]; /// assume nei is in increasing order
    bool *U_hash = new bool[U_size+1];
    for(uint i = 0; i <= U_size; i++) U_hash[i] = false;
    for(uint u : U) {
        U_hash[u] = true;
    }

    for(uint i = 0; i < CV.size(); i++){
        if(init_MEB) return;
        uint v = CV[i];

        /// 1. update _U; O(degree_R_max)
        _U.clear();
        for(uint u : Gk.adjR[v]) {
            if(u > U_size) break;
            if(u <= U_size && U_hash[u])
                _U.emplace_back(u);
        }
        uint us = _U.size();

        /// checking whether _U is a subset of v's neighbors in XV
        /** TODO-pruning: don't exist a vertex xv in XV, _U is a subset of N(v,G)
            Assume nei is in increasing order
            function:  _U is a subset of N(v,G) where v is in XV
            O(|XV|*(degree_V_max + num_xv_nei)) <= O(|XV|*degree_V_max + |E+|) */
        bool is_not_exist = true;
        for(auto xv : XV) {
            ASSERT(xv < Gk.adjR.size())
            vector<uint> &nxv = Gk.adjR[xv];
            uint num_xv_nei = nxv.size();
            int pos_v = 0;
            if(num_xv_nei >= us) {
                bool is_subset = true;
                // O(degree_R_max + degree_R_max)
                for(int pos_u = 0; pos_u < _U.size(); pos_u++) {
                    uint unode = _U[pos_u];
                    while(pos_v < num_xv_nei && nxv[pos_v] < unode) pos_v++;
                    if(pos_v >= num_xv_nei || (pos_v < num_xv_nei && nxv[pos_v] > unode)) {
                        is_subset = false;
                        break;
                    }
                }
                if(is_subset) {
                    is_not_exist = false;
                    break;
                }
            }
        }


        /// 2.update _V, _CV, _XV;
        /// assume nei is in increasing order
        if(is_not_exist){
            /// update _V  O(|E+|*degree_V_max)
            /** function:  find all vertex v in subCV that _U is a subset of N(v,G)
              *            compute subCV and put subCV to _V
              * O(|CV|*(|_U|+num_nei_cv)) <= O(|CV|*|_U|+|E-|) |_U| = degree_V_max */
            _V.clear();
            _V = V;
            _V.emplace_back(v);
            for(uint i = 0; i < CV_size; i++) reduce_CV_pos[i] =false;
            for(uint j = i+1; j < CV.size(); j++) {
                auto cv = CV[j];
                ASSERT(cv < Gk.adjR.size())
                vector<uint> &ncv = Gk.adjR[cv];
                uint num_cv_nei = ncv.size();
                int pos_v = 0;
                if(num_cv_nei >= us) {
                    bool is_subset = true;
                    // O(|_U|*num_nei_cv)
                    for(int pos_u = 0; pos_u < _U.size(); pos_u++){
                        uint unode = _U[pos_u];
                        while(pos_v < num_cv_nei && ncv[pos_v] < unode){
                            pos_v++;
                        }
                        if(pos_v >= num_cv_nei || (pos_v < num_cv_nei && ncv[pos_v] > unode)) {
                            is_subset = false;
                            break;
                        }
                    }
                    if(is_subset){
                        _V.emplace_back(cv);
                        //ASSERT(j < CV.size())
                        reduce_CV_pos[j] = true;
                    }

                }
            }

            /// update _CV O(|CV|*|_U|*|num_cv_nei|) = O(|E-|*degree_V_max)
            _CV.clear();
            for(uint j = i+1; j < CV.size(); j++) {
                if(reduce_CV_pos[j]) continue;
                auto cv = CV[j];
                vector<uint> &ncv = Gk.adjR[cv];
                uint num_cv_nei = ncv.size();
                uint num_joint = 0;
                int pos_u = 0;
                int pos_v = 0;
                // O(|_U|*|num_cv_nei|)
                while(pos_u < us && pos_v < num_cv_nei) {
                    while(pos_u < us && pos_v < num_cv_nei && _U[pos_u] < ncv[pos_v]) pos_u++;
                    while(pos_u < us && pos_v < num_cv_nei && _U[pos_u] > ncv[pos_v]) pos_v++;
                    while(pos_u < us && pos_v < num_cv_nei && _U[pos_u] == ncv[pos_v]) {
                        num_joint++;
                        pos_u++;
                        pos_v++;
                    }
                }
                if(num_joint >= t_U) {
                    _CV.emplace_back(cv);
                }
            }

            /// update _XV O(|XV|*|_U|*|num_xv_nei|)=O(|E-|*degree_V_max)
            _XV.clear();
            for(unsigned int xv : XV) {
                vector<uint> &nxv = Gk.adjR[xv];
                uint num_xv_nei = nxv.size();
                uint num_joint = 0;
                int pos_u = 0;
                int pos_v = 0;
                while(pos_u < us && pos_v < num_xv_nei) {
                    while(pos_u < us && pos_v < num_xv_nei && _U[pos_u] < nxv[pos_v]) pos_u++;
                    while(pos_u < us && pos_v < num_xv_nei && _U[pos_u] > nxv[pos_v]) pos_v++;
                    while(pos_u < us && pos_v < num_xv_nei && _U[pos_u] == nxv[pos_v]){
                        num_joint++;
                        pos_u++;
                        pos_v++;
                    }
                }
                if(num_joint >= t_U) {
                    _XV.emplace_back(xv);
                }
            }

            /// branch and bound: recall iMBEA itself
            uint vs = _V.size();
            uint cvs = _CV.size();
            if(us >= t_U && (vs + cvs >= t_V)) {
                if(problem_type == "MEB" && (us * (vs+cvs) > _C.MBsize))
                    iMBEA_init(Gk, _U, _V, _CV, _XV);
                if(problem_type == "MVB" && (us + (vs+cvs) > _C.MBsize))
                    iMBEA_init(Gk, _U, _V, _CV, _XV);
                if(problem_type == "MBB" && (min(us,(vs+cvs))*min(us,(vs+cvs)) > _C.MBsize))
                    iMBEA_init(Gk, _U, _V, _CV, _XV);
            }
        }

        /// record N(v,G) result in XV
        XV.emplace_back(v);

    } // enumerate v in CV

    delete[] U_hash;
    delete[] reduce_CV_pos;
}
void MB::initMBC() {
    if(g.noE == 0) return;
    vector<uint> U, V, CV, XV;

    //TODO: method 2. choose a vertex with max degree (made it)

    bool *visit_adjR = new bool[g.adjR.size()];
    for(uint i = 0; i < g.adjR.size(); i++) {visit_adjR[i] = true;}
    int res = 0;

    do{

        uint vv = 0;
        uint maxD = 0;
        for(uint i = 0; i < g.adjR.size(); i++) {
            if(g.adjR[i].size() > maxD && visit_adjR[i]) {
                vv = i; // find the first vv with max degree
                maxD = g.adjR[i].size();
            }
            //if(visit_adjR[i]) vv = i; // randomly
        }
        visit_adjR[vv] = false;
        //cout<<"maxD="<<maxD<<", vv="<<vv<<endl;

        //vector<uint> U, V, CV, XV;
        const uint map_size = g.adjR.size();
        bool *map_2hop = new bool[map_size];
        for(uint i = 0; i < map_size; i++) map_2hop[i] = true;

        V.emplace_back(vv);
        for(uint it : g.adjR[vv]) {
            U.emplace_back(it);
        }

        // get 2-hop nei;
        vector<uint> THnei;
        for(int i=0; i<g.adjR[vv].size(); i++) {
            uint itu = g.adjR[vv][i];

            for(int j=0; j<g.adjL[itu].size(); j++){
                uint itv = g.adjL[itu][j];
                if(itv==vv) continue;
                // itv is 2-hop nei of vv
                if(map_2hop[itv]) {
                    map_2hop[itv] = false;
                    THnei.emplace_back(itv);
                }

            }
        }

        // update V and CV; O(|THnei|*(|U|+|num_cv_nei|))
        for(auto itv: THnei) {
            vector<uint> &ncv = g.adjR[itv];
            uint num_cv_nei = ncv.size();
            uint num_joint = 0;
            uint us = U.size();
            uint pos_u = 0;
            uint pos_v = 0;
            // O(|U|+|num_cv_nei|)
            while(pos_u < us && pos_v < num_cv_nei) {
                while(pos_u < us && pos_v < num_cv_nei && U[pos_u] < ncv[pos_v]) pos_u++;
                while(pos_u < us && pos_v < num_cv_nei && U[pos_u] > ncv[pos_v]) pos_v++;
                while(pos_u < us && pos_v < num_cv_nei && U[pos_u] == ncv[pos_v]) {
                    num_joint++;
                    pos_u++;
                    pos_v++;
                }
            }

            if(num_joint==U.size()){
                V.emplace_back(itv);
            }
            else if(num_joint >= t_U) {
                CV.emplace_back(itv);
            }
        }


        if(U.size() >= t_U && V.size() >= t_V){
            if(problem_type == "MEB" && V.size() * U.size() > _C.MBsize) {
                //_C.subU = U;
                //_C.subV = V;
                _C.MBsize = V.size() * U.size();
            }
            else if (problem_type == "MVB" && V.size() + U.size() > _C.MBsize) {
                //_C.subU = U;
                //_C.subV = V;
                _C.MBsize = V.size() + U.size();
            }
            else if (problem_type == "MBB" && min(V.size(),U.size())*min(V.size(),U.size()) > _C.MBsize) {
                //_C.subU = U;
                //_C.subV = V;
                vector<uint> tmpU = U.size() < V.size() ? U : V;
                _C.MBsize = tmpU.size() * tmpU.size();
            }
            init_MEB = true;
        }
        else{
            if(U.size() >= t_U && (V.size()+CV.size())>=t_V)
                iMBEA_init(g, U, V, CV, XV);
        }

        res++;
        V.clear();
        U.clear();
        CV.clear();
        XV.clear();
        delete[] map_2hop;

    }while((!init_MEB) && res < 10);

    delete[] visit_adjR;
}

bgraph MB::Reduce(bgraph &Gk, uint tkU, uint tkV) {
    bgraph Gt = Gk;

    /// one-hop pruning (maintain ordering) O(|E|)
    Reduce1hop(Gt, tkU, tkV);

    if(!(tkU==1||tkV==1)) {
        /// two-hop pruning O(d(u,G)^2 + d(v,G)^2)
        Reduce2hop(Gt.adjL, Gt.adjR, tkU, tkV);
        Reduce2hop(Gt.adjR, Gt.adjL, tkV, tkU);
        /// one-hop pruning again
        Reduce1hop(Gt, tkU, tkV);
    }

    /// update G information
    Gt.max_degree_L = Gt.maxDegree(Gt.adjL);
    Gt.max_degree_R = Gt.maxDegree(Gt.adjR);
    Gt.noE = 0;
    Gt.calculate_number_of_edges();

    return Gt;
}
void MB::Reduce1hop(bgraph &Gt, uint &tku, uint &tkv) {
    uint num_u = Gt.adjL.size();
    uint num_v = Gt.adjR.size();
    uint *degree_u = new uint[num_u];
    uint *degree_v = new uint[num_v];
    uint *to_remove_u = new uint[num_u];
    uint *to_remove_v = new uint[num_v];
    bool *removed_u = new bool[num_u];
    bool *removed_v = new bool[num_v];
    for(uint i = 0; i < num_u; i++) {
        degree_u[i] = Gt.adjL[i].size();
    }
    for(uint i = 0; i < num_v; i++) {
        degree_v[i] = Gt.adjR[i].size();
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
            for(uint i = 0; i < Gt.adjL[uu].size(); i++){
                uint other = Gt.adjL[uu][i];
                if(!removed_v[other]){
                    degree_v[other]--;

                    if(degree_v[other] < tku){
                        to_remove_v[end_idx_v++] = other;
                        removed_v[other] = true;
                    }
                }
            }
            Gt.adjL[uu].clear();
            degree_u[uu] = 0;
        }

        while(start_idx_v != end_idx_v) {
            uint vv = to_remove_v[start_idx_v++];

            //cout << "remove : " << vv << endl;
            for(uint i = 0; i < Gt.adjR[vv].size(); i++){
                uint other = Gt.adjR[vv][i];
                if(!removed_u[other]){
                    degree_u[other]--;

                    if(degree_u[other] < tkv){
                        to_remove_u[end_idx_u++] = other;
                        removed_u[other] = true;
                    }
                }
            }
            Gt.adjR[vv].clear();
            degree_v[vv] = 0;
        }

    }

    // update adj
    for(uint i = 0; i < num_u; i++){
        if(!removed_u[i]){
            vector<uint> new_vec(degree_u[i]);
            uint idx = 0;
            for(uint j = 0; j < Gt.adjL[i].size(); j++){
                if(!removed_v[Gt.adjL[i][j]]){
                    new_vec[idx++] = Gt.adjL[i][j];
                }
            }
            Gt.adjL[i].clear();
            Gt.adjL[i] = new_vec;
        }
    }
    for(uint i = 0; i < num_v; i++){
        if(!removed_v[i]){
            vector<uint> new_vec(degree_v[i]);
            uint idx = 0;
            for(uint j = 0; j < Gt.adjR[i].size(); j++){
                if(!removed_u[Gt.adjR[i][j]]){
                    new_vec[idx++] = Gt.adjR[i][j];
                }
            }
            Gt.adjR[i].clear();
            Gt.adjR[i] = new_vec;
        }
    }

    delete[] degree_u;
    delete[] degree_v;
    delete[] to_remove_u;
    delete[] to_remove_v;
    delete[] removed_u;
    delete[] removed_v;
}
void MB::Reduce2hop(vector<vector<uint>> &U, vector<vector<uint>> &V, uint &tku, uint &tkv) {
    if(tkv == 1 || tku == 1) return; // this special case may don't have 2-hop nei
    uint num_u = U.size();
    uint num_v = V.size();
    map<uint, uint> cnt_common_nei; // degree is small
    bool *removed_u = new bool[num_u];
    for(uint i = 0; i < num_u; i++) removed_u[i] = false;

    // delete u0 whose the number of common nei with its number lager than tkv is smaller than tku
    // O(sum(deg(u)^2)) or O(sum(deg(v)^2)) in paper
    // <= O(|E|deg(v)_max log(|U|)) or O(|E|deg(u)_max) I think
    for(uint i = 0; i < num_u; i++){
        if(U[i].size() < tku) continue;
        uint u0 = i;

        for(uint j = 0; j < U[i].size(); j++) {
            uint vv = U[u0][j];
            for(uint k = 0; k < V[vv].size(); k++){
                uint uu = V[vv][k];
                if(cnt_common_nei.count(uu)==0){
                    cnt_common_nei[uu] = 1;
                }
                else{
                    cnt_common_nei[uu]++;
                }
            }
        }


        uint cnt2hop = 0;
        for(auto it : cnt_common_nei) {
            if(it.second >= tkv) cnt2hop++;
        }
        if(cnt2hop < tku) {
            removed_u[u0] = true;
            U[u0].clear();
        }
        cnt_common_nei.clear();
    }

    // update V adj O(|E|)
    for(uint i = 0; i < num_v; i++) {
        vector<uint> new_vec;
        for(uint j = 0; j < V[i].size(); j++) {
            uint uu = V[i][j];
            if(!removed_u[uu]) {
                new_vec.emplace_back(uu);
            }
        }
        V[i].clear();
        V[i] = new_vec;
    }

    delete[] removed_u;

}

ulong MB::MBC_improved() {
    if(g.noE == 0) {
        cout<<"The graph is empty."<<endl;
        return 0;
    }

    /*g.max_degree_L = g.maxDegree(g.adjL);
    g.max_degree_R = g.maxDegree(g.adjR);
    g.noE = 0;
    g.calculate_number_of_edges();
    t_U = tkU;
    t_V = tkV;*/

    // sort neighbors with non-decreasing order by degrees of neighbors
    g.sortAdj();

    // find a large biclique
    initMBC();

    // find the size of maximum biclique
    uint t_0V, t_1U, t_1V;
    t_0V = g.max_degree_L;
    if(problem_type == "MEB") {
        while(t_0V > t_V) {
            t_1U = max(_C.MBsize/t_0V, (ulong)t_U);
            t_1V = max(t_0V/2, t_V);  /// t_1V = max(t_0V-1, t_V);  faster method
            bgraph Gk = Reduce(g, t_1U, t_1V); // one-hop and two-hop reduce
            if(!Gk.adjL.empty()) MBC_Gk(Gk, t_1U, t_1V);
            t_0V = t_1V;  /// t_0V = t_1V/2;
        }
    }
    else if(problem_type == "MVB") {
        while(t_0V > t_V) {
            ulong tmpS = _C.MBsize > t_0V ? (_C.MBsize - t_0V) : 0;
            t_1U = max(tmpS, (ulong)t_U);
            t_1V = max(t_0V/2, t_V);  /// t_1V = max(t_0V-1, t_V);  faster method
            bgraph Gk = Reduce(g, t_1U, t_1V); // one-hop and two-hop reduce
            if(!Gk.adjL.empty()) MBC_Gk(Gk, t_1U, t_1V);
            t_0V = t_1V;  /// t_0V = t_1V/2;
        }
    }
    else if (problem_type == "MBB") {
        uint mean1, tt;
        tt = max(t_U, t_V);
        if(_C.MBsize > 0 && tt < sqrt(_C.MBsize)) tt = sqrt(_C.MBsize);
        mean1 = g.max_degree_L < g.max_degree_R ? g.max_degree_L : g.max_degree_R;
        ulong last_MBB = _C.MBsize;
        while (mean1 > tt && last_MBB == _C.MBsize) {
            if(mean1 / tt > 10) mean1 = max(tt, mean1 / 2);
            else mean1 = max(tt, (mean1+tt) / 2);
            //mean1 = max(tt, mean1 / 2);
//cout<<"mean="<<mean1<<", MBsize="<<_C.MBsize<<endl;
            bgraph Gk = Reduce(g, mean1, mean1); // one-hop and two-hop reduce
            if(Gk.noE != 0) MBC_Gk(Gk, mean1, mean1);
        }
    }

    return _C.MBsize;
}