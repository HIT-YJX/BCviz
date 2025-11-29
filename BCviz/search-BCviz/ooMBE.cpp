
#include "ooMBE.h"

ooMBE::~ooMBE() {
    _C.subU.clear();
    _C.subV.clear();

    for(auto & i : ordered_vex) i.clear();

    if(L_order!= nullptr) {delete [] L_order;L_order=nullptr;}
    if(R_order!= nullptr) {delete [] R_order;R_order=nullptr;}

    if(idMapposL != nullptr) {delete[] idMapposL; idMapposL = nullptr;}
    if(idMapposR != nullptr) {delete[] idMapposR; idMapposR = nullptr;}
    if(posMapidL != nullptr) {delete[] posMapidL; posMapidL = nullptr;}
    if(posMapidR != nullptr) {delete[] posMapidR; posMapidR = nullptr;}
    if(visitL != nullptr) {delete[] visitL; visitL = nullptr;}
    if(visitR != nullptr) {delete[] visitR; visitR = nullptr;}
}

/****** tool algorithm ******/
void ooMBE::printMapVec(vector<uint> &vec, uint row) {
    if(row == 0) {
        for(auto &it: vec) {
            cout<<posMapidL[it]<<" ";
        }
    }
    else if (row == 1){
        for(auto &it: vec) {
            cout<<posMapidR[it]<<" ";
        }
    }
    cout<<endl;
}

void ooMBE::printMEB(Biclique &bc) {
    cout<<"print MEB: ";
    cout<<"size=p*q="<<bc.MBsize<<endl;
    cout<<"U: ";
    for(auto &it: bc.subU) {
        cout<<it<<" ";
    }
    cout<<endl;
    cout<<"V: ";
    for(auto &it: bc.subV) {
        cout<<it<<" ";
    }
    cout<<endl;
}
void ooMBE::printMapMEB(Biclique &bc) {
    cout<<"print posMapid MEB: ";
    cout<<"size=p*q="<<bc.MBsize<<endl;
    cout<<"U: ";
    for(auto &it: bc.subU) {
        cout<<posMapidL[it]<<" ";
    }
    cout<<endl;
    cout<<"V: ";
    for(auto &it: bc.subV) {
        cout<<posMapidR[it]<<" ";
    }
    cout<<endl;
}

void ooMBE::calculate_sizeLR(bgraph &Gt) {
    sizeL = 0;
    sizeR = 0;
    for(auto &it: Gt.adjL) {
        if(!it.empty()) sizeL++;
    }
    for(auto &it: Gt.adjR) {
        if(!it.empty()) sizeR++;
    }
}
/****** tool algorithm ******/


/****** baseline improved ******/
void ooMBE::iMBEA_Gk(bgraph &Gk,
                  uint &tkU,
                  uint &tkV,
                  vector<uint> &U,
                  vector<uint> &V,
                  vector<uint> &CV,
                  vector<uint> &XV) {
    // TODO: nei_id is in increasing order
    if(U.empty()) {return;}
    nomb++;

    /// update maximum edge Biclique
    if((U.size() >= tkU) && (V.size() >= tkV)) {
        if(problem_type == "MEB" && V.size() * U.size() > _C.MBsize) {
            _C.subU = U;
            _C.subV = V;
            _C.MBsize = V.size() * U.size();
        }
        else if (problem_type == "MVB" && V.size() + U.size() > _C.MBsize) {
            _C.subU = U;
            _C.subV = V;
            _C.MBsize = V.size() + U.size();
        }
        else if (problem_type == "MBB" && min(V.size(),U.size())*min(V.size(),U.size()) > _C.MBsize) {
            _C.subU = U;
            _C.subV = V;
            vector<uint> tmpU = U.size() < V.size() ? U : V;
/*            _C.subU = tmpU;
            _C.subV = tmpU;*/
            _C.MBsize = tmpU.size() * tmpU.size();
        }
    }



    /*cout<<"--------------information start--------------\n";
    cout<<"U array:"<<endl;
    printVector(U);
    cout<<"V array:"<<endl;
    printVector(V);
    cout<<"CV array:"<<endl;
    printVector(CV);
    cout<<"XV array:"<<endl;
    printVector(XV);
    cout<<"--------------information end--------------\n";*/

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
        /// nei is in increasing order
        if(is_not_exist){

            /// update _V <=O(|E+|*degree_V_max)
            /// function:  find all vertex v in subCV that _U is a subset of N(v,G)
            ///            compute subCV and put subCV to _V
            /// O(|CV|*(|_U|+num_nei_cv)) <= O(|CV|*|_U|+|E-|) |_U| = degree_V_max
            _V.clear();
            _V = V;
            _V.emplace_back(v);
            memset(reduce_CV_pos, false, sizeof(reduce_CV_pos));
            for(uint j = i+1; j < CV.size(); j++) {
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


void ooMBE::MBC_Gk(bgraph &Gk, uint tkU, uint tkV) {
    if(Gk.adjL.empty() || Gk.adjR.empty()) return;

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


    /*cout<<"--------------information start--------------\n";
    cout<<"U array:"<<endl;
    printVector(U);
    cout<<"V array:"<<endl;
    printVector(V);
    cout<<"CV array:"<<endl;
    printVector(CV);
    cout<<"XV array:"<<endl;
    printVector(XV);
    cout<<"--------------information end--------------\n";*/



    // call iMBEA(Gk)
    if(U.size() >= tkU && CV.size() >= tkV)
        iMBEA_Gk(Gk, tkU, tkV, U, V, CV, XV);

    U.clear();
    V.clear();
    CV.clear();
    XV.clear();

}

void ooMBE::iMBEA_init(bgraph &Gk,
                       vector<uint> &U,
                       vector<uint> &V,
                       vector<uint> &CV,
                       vector<uint> &XV) {
    if(U.empty() || init_MEB) {return;}
    nomb++;

    /// update maximum edge Biclique
    if(U.size() >= t_U && V.size() >= t_V) {
        if(problem_type == "MEB" && V.size() * U.size() > _C.MBsize) {
            _C.subU = U;
            _C.subV = V;
            _C.MBsize = V.size() * U.size();
        }
        else if (problem_type == "MVB" && V.size() + U.size() > _C.MBsize) {
            _C.subU = U;
            _C.subV = V;
            _C.MBsize = V.size() + U.size();
        }
        else if (problem_type == "MBB" && min(V.size(),U.size())*min(V.size(),U.size()) > _C.MBsize) {
            _C.subU = U;
            _C.subV = V;
            vector<uint> tmpU = U.size() < V.size() ? U : V;
/*            _C.subU = tmpU;
            _C.subV = tmpU;*/
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

    const uint U_size = U[U.size()-1]; /// nei is in increasing order
    bool *U_hash = new bool[U_size+1];
    for(uint i = 0; i <= U_size; i++) U_hash[i] = false;
    for(uint u : U) {
        //ASSERT(u <= U_size)
        U_hash[u] = true;
    }

    for(uint i = 0; i < CV.size(); i++){
        if(init_MEB) return;
        uint v = CV[i];

        /// 1. update _U; O(degree_R_max)
        _U.clear();
        ASSERT(v < Gk.adjR.size())
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
                    // TODO: check! ------right
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
                        ASSERT(j < CV.size())
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
void ooMBE::initMBC() {
    if(g.adjR.empty()) return;
    vector<uint> U, V, CV, XV;

    // choose a vertex with max degree
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
        /*cout<<"--------------information start--------------\n";
        cout<<"U array:"<<endl;
        printVector(U);
        cout<<"V array:"<<endl;
        printVector(V);
        cout<<"CV array:"<<endl;
        printVector(CV);
        cout<<"XV array:"<<endl;
        printVector(XV);
        cout<<"--------------information end--------------\n";*/


        if(U.size() >= t_U && V.size() >= t_V){
            if(problem_type == "MEB" && V.size() * U.size() > _C.MBsize) {
                _C.subU = U;
                _C.subV = V;
                _C.MBsize = V.size() * U.size();
            }
            else if (problem_type == "MVB" && V.size() + U.size() > _C.MBsize) {
                _C.subU = U;
                _C.subV = V;
                _C.MBsize = V.size() + U.size();
            }
            else if (problem_type == "MBB" && min(V.size(),U.size())*min(V.size(),U.size()) > _C.MBsize) {
                _C.subU = U;
                _C.subV = V;
                vector<uint> tmpU = U.size() < V.size() ? U : V;
/*                _C.subU = tmpU;
                _C.subV = tmpU;*/
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

bgraph ooMBE::Reduce(bgraph &Gk, uint tkU, uint tkV) {
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
void ooMBE::Reduce1hop(bgraph &Gt, uint &tku, uint &tkv) {
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
void ooMBE::Reduce2hop(vector<vector<uint>> &U,  vector<vector<uint>> &V, uint &tku, uint &tkv) {
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

void ooMBE::MBC_improved(uint tkU, uint tkV) {
    if(g.adjL.empty()) {
        //cout<<"The graph is empty."<<endl;
        return;
    }
    uint originU = tkU;
    uint originV = tkV;
    if(problem_type=="MBB") {
        uint tt = max(tkU, tkV);
        tkU = tt;
        tkV = tt;
    }

    ///cout<<"--------------reduced graph--------------"<<endl;
    auto start = high_resolution_clock::now();
    // reduce original graph
    if(tkU >= 2 || tkV >= 2){
        Reduce1hop(g, tkU, tkV);
    }

    if(g.adjL.size() < tkU || g.adjR.size() < tkV) {
        cout<<"|adjL| < tkU or |adjR| < tkV."<<endl;
        return;
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    double reduce_time = duration.count() * 1.0 / 1000.0;
    calculate_sizeLR(g);
    g.max_degree_L = g.maxDegree(g.adjL);
    g.max_degree_R = g.maxDegree(g.adjR);
    g.noE = 0;
    g.calculate_number_of_edges();

/*cout<<"------1. (4,4) degree-based reduced graph-------"<<endl;
cout<<"number of vertices in U: "<<sizeL<<endl;
cout<<"number of vertices in V: "<<sizeR<<endl;
cout<<"number of edges in E: "<<g.noE<<endl;
cout<<"max degree of U: "<<g.max_degree_L<<endl;
cout<<"max degree of V: "<<g.max_degree_R<<endl;
cout<<"reduce original graph time: "<<reduce_time<<" ms"<<endl<<endl;*/
    //g.dispalyAdj();

    ///cout<<"--------------prepare--------------"<<endl;
    start = high_resolution_clock::now();
    // initMBC
    /*if(g.max_degree_R < g.max_degree_L) {
        // github good, team bad
        g.switchLR();
        uint tmp = tkU;
        tkU = tkV;
        tkV = tmp;
        switchLR = !switchLR;
        cout<<"switchLR!"<<endl;
        //g.dispalyAdj();
    }*/
    t_U = tkU;
    t_V = tkV;

    // todo: sort neighbors with non-decreasing order by degrees of neighbors
    g.sortAdj();
    //cout<<"sorting graph"<<endl;
    //g.dispalyAdj();

    // find a large biclique
    initMBC();
    uint initSize = _C.MBsize;
    /*if(_C.MBsize == 0) {
        cout<<"Can't find a initMBC."<<endl;
        exit(0);
    }*/
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    double initMEB_time = duration.count() * 1.0 / 1000.0;
    //printMEB(_C);
//cout<<"init MBsizze="<<_C.MBsize<<", find time: "<<initMEB_time<<" ms"<<endl<<endl;


    /// test of Cohesion-based reduction (MEB)
/*    bgraph testG = g;
    uint _tkU = _C.MBsize / tkV;
    //uint _tkV = _C.MBsize / tkU;

    start = high_resolution_clock::now();
    Reduce1hop(testG, _tkU, tkV);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    double cohesion_reduce_time = duration.count() * 1.0 / 1000.0;

    uint numL, numR, degreeL, degreeR, numE;
    numL = numR = degreeL = degreeR = numE =0;
    for(const auto& itu: testG.adjL) {
        if(itu.size() > 0) {
            numL++;
            numE += itu.size();
            degreeL = max(degreeL, (uint)itu.size());
        }
    }
    for(const auto& itv: testG.adjR) {
        if(itv.size() > 0) {
            numR++;
            degreeR = max(degreeR, (uint)itv.size());
        }
    }*/
/*cout<<"------2. (4,4) cohesion-based reduced graph (1)-------"<<endl;
cout<<"s="<<_tkU<<", t="<<tkV<<endl;
cout<<"number of vertices in U: "<<numL<<endl;
cout<<"number of vertices in V: "<<numR<<endl;
cout<<"number of edges in E: "<<numE<<endl;
cout<<"max degree of U: "<<degreeL<<endl;
cout<<"max degree of V: "<<degreeR<<endl;
cout<<"reduce original graph time: "<<cohesion_reduce_time+reduce_time<<" ms"<<endl<<endl;*/

/*    testG = g;
    uint _tkV = _C.MBsize / tkU;

    start = high_resolution_clock::now();
    Reduce1hop(testG, tkU, _tkV);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    cohesion_reduce_time = duration.count() * 1.0 / 1000.0;

    numL = numR = degreeL = degreeR = numE =0;
    for(const auto& itu: testG.adjL) {
        if(itu.size() > 0) {
            numL++;
            numE += itu.size();
            degreeL = max(degreeL, (uint)itu.size());
        }
    }
    for(const auto& itv: testG.adjR) {
        if(itv.size() > 0) {
            numR++;
            degreeR = max(degreeR, (uint)itv.size());
        }
    }*/
/*    cout<<"------2. (4,4) cohesion-based reduced graph (2)-------"<<endl;
    cout<<"s="<<tkU<<", t="<<_tkV<<endl;
    cout<<"number of vertices in U: "<<numL<<endl;
    cout<<"number of vertices in V: "<<numR<<endl;
    cout<<"number of edges in E: "<<numE<<endl;
    cout<<"max degree of U: "<<degreeL<<endl;
    cout<<"max degree of V: "<<degreeR<<endl;
    cout<<"reduce original graph time: "<<cohesion_reduce_time+reduce_time<<" ms"<<endl<<endl;*/



    ///cout<<"--------------find MEB--------------"<<endl;
    auto allstart = high_resolution_clock::now();

    uint t_0V, t_1U, t_1V;
    t_0V = g.max_degree_L;
    /// main idea: bounding
    if(problem_type == "MEB") {
        uint ite = 1;
        while(t_0V > t_V) {
            t_1U = max(_C.MBsize/t_0V, (ulong)t_U);
            t_1V = max(t_0V/2, t_V);  /// t_1V = max(t_0V-1, t_V);  faster method
            bgraph Gk = Reduce(g, t_1U, t_1V); // one-hop and two-hop reduce
/*{
    numL = numR = degreeL = degreeR = numE =0;
    for(const auto& itu: Gk.adjL) {
        if(itu.size() > 0) {
            numL++;
            numE += itu.size();
            degreeL = max(degreeL, (uint)itu.size());
        }
    }
    if(numE == 0) {t_0V = t_1V; continue;}
    for(const auto& itv: Gk.adjR) {
        if(itv.size() > 0) {
            numR++;
            degreeR = max(degreeR, (uint)itv.size());
        }
    }
    cout<<"------3. (4,4) progressive reduced graph ("<<ite++<<")-------"<<endl;
    cout<<"s="<<t_1U<<", t="<<t_1V<<endl;
    cout<<"number of vertices in U: "<<numL<<endl;
    cout<<"number of vertices in V: "<<numR<<endl;
    cout<<"number of edges in E: "<<numE<<endl;
    cout<<"max degree of U: "<<degreeL<<endl;
    cout<<"max degree of V: "<<degreeR<<endl;
    //cout<<"reduce original graph time: "<<cohesion_reduce_time+reduce_time<<" ms"<<endl<<endl;
}*/
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
            mean1 = max(tt, mean1 / 2);
//cout<<"mean="<<mean1<<", MBsize="<<_C.MBsize<<endl;
            bgraph Gk = Reduce(g, mean1, mean1); // one-hop and two-hop reduce
            if(Gk.noE != 0) MBC_Gk(Gk, mean1, mean1);
        }
    }
    auto allstop = high_resolution_clock::now();
    auto allduration = duration_cast<microseconds>(allstop - allstart);
    double allfindtime = allduration.count() * 1.0 / 1000.0;
    //double findtime = alltime1 + alltime2;
    //cout<<"Finding MEB time: "<<findtime<<" ms. (where reduce time is "<<alltime1<<" ms, MBC time is "<<alltime2<<" ms)"<<endl<<endl;



    /** print result
     * [algorithm] [dataset] [tu] [tv] [time (m seconds)] [size of MEB(p*q)]  [initMBCsize] */
    double exetime = allfindtime;
    string dataset_name;
    const char split = '/';
    istringstream iss(g.datasets);
    string token;
    while (getline(iss, token, split))
    {
        dataset_name = token;
    }

    cout<<problem_type<<" MBC "<<dataset_name<<" "<<originU<<" "<<originV<<" "<<exetime<<" "
    <<_C.MBsize<<" "<<initSize<<" "<<initMEB_time<<endl;

    /*cout<<"\n------------MBC result-------------"<<endl;
    printMEB(_C);*/
}
/****** baseline improved ******/



/****** index search ******/
void ooMBE::readBCviz(string filename) {
    FILE* fp2 = nullptr;
    fp2 = freopen(filename.c_str(), "r", stdin);
    if(!fp2){
        printf("Can't read %s file! \n", filename.c_str());
        exit(0);
    }

    int c_size, num_v;
    uint r, vid, os, ts, deg;
    os = 0;
    L_order = new uint[g.adjL.size()];
    R_order = new uint[g.adjR.size()];

    cin>>c_size;
    uint order_pos = 0;
    for(int i = 0; i < c_size; i++){
        cin>>num_v;
        vector<SortedBySize> component;
        for(int j = 0; j < num_v; j++){
            //cin>>r>>vid>>ts>>os>>deg;
            cin>>r>>vid>>ts;
            SortedBySize node(r, (vid-1), os, ts);
            component.emplace_back(node);

            if(r == 0) {
                g.cohesionL[vid-1] = ts;
                L_order[vid-1] = order_pos;
            }
            else{
                g.cohesionR[vid-1] = ts;
                R_order[vid-1] = order_pos;
            }
            order_pos++;
        }
        ordered_vex.emplace_back(component);
    }

    fclose(stdin);
    cin.clear();
//cout<<"reading the path of BCviz index: "<<filename<<endl;
}
void ooMBE::printBCviz() {
    cout<<"ordered_vex.num="<<ordered_vex.size()<<endl;
    for(int i = 0; i < ordered_vex.size(); i++){
        cout<<ordered_vex[i].size()<<endl;
        cout<<"twoNodeSize, oneNodeSize, row, vex"<<endl;
        for(auto &it: ordered_vex[i]){
            cout<<it.twoNodeSize<<","<<it.oneNodeSize<<","<<it.row<<","<<it.vex<<endl;
        }
        cout<<endl;
    }
    cout<<endl;
}
ulong ooMBE::BicliqueSize(uint l, uint r){
    if(problem_type == "MEB") return l * r;
    if(problem_type == "MVB") return l + r;
    if(problem_type == "MBB") return min(l*l, r*r);
    return 0;
}
void ooMBE::initMBC_Gk(bgraph &Gt) {
    if(Gt.adjR.empty()) return;
    vector<uint> U, V, CV, XV;

    //TODO: method 1. choose a vertex randomly

    //TODO: method 2. choose a vertex with max degree (made it)

    bool *visit_adjR = new bool[Gt.adjR.size()];
    //memset(visit_adjR, true, sizeof visit_adjR);
    for(uint i = 0; i < Gt.adjR.size(); i++) {visit_adjR[i] = true;}
    int res = 0;

    do{

        uint vv = 0;
        uint maxD = 0;
        for(uint i = 0; i < Gt.adjR.size(); i++) {
            if(Gt.adjR[i].size() > maxD && visit_adjR[i]) {
                vv = i; // find the first vv with max degree
                maxD = Gt.adjR[i].size();
            }
        }
        visit_adjR[vv] = false;
        //cout<<"maxD="<<maxD<<", vv="<<vv<<endl;

        //vector<uint> U, V, CV, XV;
        const uint map_size = Gt.adjR.size();
        bool *map_2hop = new bool[map_size];
        //memset(map_2hop, true, sizeof map_2hop);
        for(uint i = 0; i < map_size; i++) map_2hop[i] = true;

        V.emplace_back(vv);
        for(uint it : Gt.adjR[vv]) {
            U.emplace_back(it);
        }

        // get 2-hop nei;
        vector<uint> THnei;
        for(int i=0; i<Gt.adjR[vv].size(); i++) {
            uint itu = Gt.adjR[vv][i];

            for(int j=0; j<Gt.adjL[itu].size(); j++){
                uint itv = Gt.adjL[itu][j];
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
            vector<uint> &ncv = Gt.adjR[itv];
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
        /*cout<<"--------------information start--------------\n";
        cout<<"U array:"<<endl;
        printVector(U);
        cout<<"V array:"<<endl;
        printVector(V);
        cout<<"CV array:"<<endl;
        printVector(CV);
        cout<<"XV array:"<<endl;
        printVector(XV);
        cout<<"--------------information end--------------\n";*/


        if(U.size() >= t_U && V.size() >= t_V){
            if(problem_type == "MEB" && V.size() * U.size() > _C.MBsize) {
                _C.subU = U;
                _C.subV = V;
                _C.MBsize = V.size() * U.size();
            }
            else if (problem_type == "MVB" && V.size() + U.size() > _C.MBsize) {
                _C.subU = U;
                _C.subV = V;
                _C.MBsize = V.size() + U.size();
            }
            else if (problem_type == "MBB" && min(V.size(),U.size())*min(V.size(),U.size()) > _C.MBsize) {
                _C.subU = U;
                _C.subV = V;
                vector<uint> tmpU = U.size() < V.size() ? U : V;
/*                _C.subU = tmpU;
                _C.subV = tmpU;*/
                _C.MBsize = tmpU.size() * tmpU.size();
            }
            init_MEB = true;
        }
        else{
            if(U.size() >= t_U && (V.size()+CV.size())>=t_V)
                iMBEA_init(Gt, U, V, CV, XV);
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
void ooMBE::MBC_improved_Gk(bgraph &Gt, uint tkU, uint tkV) {
    if(Gt.adjL.empty()) {
        cout<<"The graph is empty."<<endl;
        return;
    }

    uint t_0V, t_1U, t_1V;
    t_0V = Gt.max_degree_L;

    if(problem_type == "MEB") {
        while(t_0V > t_V) {
            t_1U = max(_C.MBsize / t_0V, (ulong) t_U);
            t_1V = max(t_0V / 2, t_V);  /// t_1V = max(t_0V-1, t_V);  faster method
//cout<<"MEBsize="<<_C.MBsize<<", t_1U="<<t_1U<<", t_1V="<<t_1V<<endl;
            bgraph Gk = Reduce(Gt, t_1U, t_1V); // one-hop and two-hop reduce
            if (!Gk.adjL.empty()) MBC_Gk(Gk, t_1U, t_1V);
            t_0V = t_1V;  /// t_0V = t_1V/2;
        }
    }
    else if (problem_type == "MVB") {
        while(t_0V > t_V) {
            ulong tmpS = _C.MBsize > t_0V ? _C.MBsize - t_0V : 0;
            t_1U = max(tmpS, (ulong) t_U);
            t_1V = max(t_0V / 2, t_V);  /// t_1V = max(t_0V-1, t_V);  faster method
//cout<<"MEBsize="<<_C.MBsize<<", t_1U="<<t_1U<<", t_1V="<<t_1V<<endl;
            bgraph Gk = Reduce(Gt, t_1U, t_1V); // one-hop and two-hop reduce
            if (!Gk.adjL.empty()) MBC_Gk(Gk, t_1U, t_1V);
            t_0V = t_1V;  /// t_0V = t_1V/2;
        }
    }
    else if (problem_type == "MBB") {
        uint mean1, tt;
        tt = max(t_U, t_V);
        if(_C.MBsize > 0 && tt < sqrt(_C.MBsize)) tt = sqrt(_C.MBsize);
        mean1 = Gt.max_degree_L < Gt.max_degree_R ? Gt.max_degree_L : Gt.max_degree_R;
        ulong last_MBB = _C.MBsize;
        while (mean1 > tt && last_MBB == _C.MBsize) {
            if(mean1 / tt > 10) mean1 = max(tt, mean1 / 2);
            else mean1 = max(tt, (mean1+tt) / 2);
            //mean1 = max(tt, mean1 / 2);
//cout<<"mean="<<mean1<<", MBsize="<<_C.MBsize<<endl;
            bgraph Gk = Reduce(Gt, mean1, mean1); // one-hop and two-hop reduce
            if(Gk.noE != 0) MBC_Gk(Gk, mean1, mean1);
        }
    }

}
void ooMBE::initBiclique2Vertex(bgraph& Gk, uint &vv) {
    vector<uint> U, V, CV, XV;

    const uint map_size = Gk.adjR.size();
    bool *map_2hop = new bool[map_size];
    //memset(map_2hop, true, sizeof map_2hop);
    for(uint i = 0; i < map_size; i++) map_2hop[i] = true;

    V.emplace_back(vv);
    for(uint it : Gk.adjR[vv]) {
        U.emplace_back(it);
    }

    // get 2-hop nei;
    vector<uint> THnei;
    for(int i=0; i<Gk.adjR[vv].size(); i++) {
        uint itu = Gk.adjR[vv][i];

        for(int j=0; j<Gk.adjL[itu].size(); j++){
            uint itv = Gk.adjL[itu][j];
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
        vector<uint> &ncv = Gk.adjR[itv];
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
/*cout<<"--------------information start--------------\n";
cout<<"U array:"<<endl;
printVector(U);
cout<<"V array:"<<endl;
printVector(V);
cout<<"CV array:"<<endl;
printVector(CV);
cout<<"XV array:"<<endl;
printVector(XV);
cout<<"--------------information end--------------\n";*/

//cout<<"U.size="<<U.size()<<"CV.size="<<CV.size()<<endl;

    if(U.size() >= t_U && V.size() >= t_V){
        if(problem_type == "MEB"){
            _C.MBsize = V.size() * U.size();
            _C.subV = V;
            _C.subU = U;
        }
        else if(problem_type == "MVB"){
            _C.MBsize = V.size() + U.size();
            _C.subV = V;
            _C.subU = U;
        }
        else if(problem_type == "MBB"){
            _C.subV = V;
            _C.subU = U;
            vector<uint> tmpU = U.size() < V.size() ? U : V;
/*            _C.subU = tmpU;
            _C.subV = tmpU;*/
            _C.MBsize = tmpU.size() * tmpU.size();
        }

        //init_MEB = true;
    }

    if(U.size() >= t_U && (V.size()+CV.size())>=t_V) {
        /// initMEB is slow
        if(min(CV.size(),U.size()) <= t_V * 3) iMBEA_Gk(Gk, t_U, t_V, U, V, CV, XV); //exact MEB
        else if(min(CV.size(),U.size()) <= t_V * 1000) {
            iMBEA_init(Gk, U, V, CV, XV); //approximate MEB
            MBC_improved_Gk(Gk, t_U, t_V); // fast exact MEB
        }
        else iMBEA_init(Gk, U, V, CV, XV); //approximate MEB
    }


    V.clear();
    U.clear();
    CV.clear();
    XV.clear();
    THnei.clear();
    delete[] map_2hop;
}

void ooMBE::initBiclique(bgraph& Gk) {
    vector<uint> U, V, CV, XV;
    for(uint i = 0; i < Gk.adjL.size(); i++){
        if(!Gk.adjL.empty())
            U.push_back(i);
    }
    for(uint i = 0; i < Gk.adjR.size(); i++) {
        if(!Gk.adjR.empty())
            CV.push_back(i);
    }
//cout<<"U.size="<<U.size()<<"CV.size="<<CV.size()<<endl;
    if(U.size() >= t_U && CV.size()>=t_V) {
        /// initMEB is slow
        if(min(CV.size(),U.size()) <= t_V * 3) iMBEA_Gk(Gk, t_U, t_V, U, V, CV, XV); //exact MEB
        else if(min(CV.size(),U.size()) <= t_V * 1000) {
            iMBEA_init(Gk, U, V, CV, XV); //approximate MEB
            MBC_improved_Gk(Gk, t_U, t_V); // fast exact MEB
        }
        else iMBEA_init(Gk, U, V, CV, XV); //approximate MEB
    }

    V.clear();
    U.clear();
    CV.clear();
    XV.clear();
}

void ooMBE::init_BCviz() {
    g.sortAdj();
    uint Lsize = g.adjL.size();
    uint Rsize = g.adjR.size();
    bool *reduce_L = new bool[Lsize];
    bool *reduce_R = new bool[Rsize];
    for(uint i = 0; i < Lsize; i++) {
        reduce_L[i] = true;
    }
    for(uint i = 0; i < Rsize; i++) {
        reduce_R[i] = true;
    }

    uint v0, max_ts, size_component, num_L, num_R, lsize, rsize, resL, resR, num_E;
    long int pos_i, pos_j, pos_l, pos_r, i, j, k;
    v0 = 0;
    max_ts = 0;
    pos_i = pos_j = 0;

    // get v0 with largest local maximum cohesion
    for(i = 0; i < ordered_vex.size(); i++) {
        for(j = 0; j < ordered_vex[i].size(); j++) {
            auto node = ordered_vex[i][j];
            if( (node.row) != switchLR && node.twoNodeSize > max_ts) {
                //(node.row == 1 && !switchLR) || (node.row == 0 && switchLR)
                max_ts = node.twoNodeSize;
                v0 = node.vex;
                pos_i = i;
                pos_j = j;
            }
        }
    }


    num_L = 0;  // v0
    reduce_R[v0] = false;
    num_R = 1;
    rsize = size_component = max_ts;
    lsize = 0;
    if(pos_j != 0) lsize = ordered_vex[pos_i][pos_j-1].twoNodeSize;
    pos_r = pos_j + 1;
    pos_l = pos_j - 1;

    // reduce graph: find candidate vertex set in BCviz
    while(BicliqueSize(num_L, num_R) < size_component || num_R < t_V || num_L < t_U) {
        j = pos_r;
        while( (BicliqueSize(num_L, num_R) < size_component || num_R < t_V || num_L < t_U) && j < ordered_vex[pos_i].size()) {
            auto &node = ordered_vex[pos_i][j];
            if( node.twoNodeSize < size_component) {
                size_component = node.twoNodeSize;
            }

            if((node.row) == switchLR) {
                //(node.row == 0 && !switchLR) || (node.row == 1 && switchLR)
                reduce_L[node.vex] = false;
                num_L++;
            }
            else{
                reduce_R[node.vex] = false;
                num_R++;
            }
            j++;

            if(size_component < lsize) {
                rsize = size_component;
                break;
            }
        }
        pos_r = j;

        k = pos_l;
        while((BicliqueSize(num_L, num_R) < size_component || num_R < t_V || num_L < t_U) && k >= 0) {
            auto &node = ordered_vex[pos_i][k];
            if( node.twoNodeSize < size_component) {
                size_component = node.twoNodeSize;
            }

            if((node.row) == switchLR) {
                //(node.row == 0 && !switchLR) || (node.row == 1 && switchLR)
                reduce_L[node.vex] = false;
                num_L++;
            }
            else{
                reduce_R[node.vex] = false;
                num_R++;
            }
            k--;

            if(size_component < rsize) {
                lsize = size_component;
                break;
            }
        }
        pos_l = k;
    }

    // get the last vetices from right line
    for(i = pos_r; i < ordered_vex[pos_i].size(); i++) {
        auto node = ordered_vex[pos_i][i];
        if( node.twoNodeSize < size_component ) {
            break;
        }
        if((node.row) == switchLR) {
            //(node.row == 0 && !switchLR) || (node.row == 1 && switchLR)
            reduce_L[node.vex] = false;
            num_L++;
        }
        else{
            reduce_R[node.vex] = false;
            num_R++;
        }

    }

    // get the last vetices from left line
    for(i = pos_l; i >= 0; i--) {
        auto node = ordered_vex[pos_i][i];
        if( node.twoNodeSize < size_component ) {
            break;
        }
        if((node.row) == switchLR) {
            //(node.row == 0 && !switchLR) || (node.row == 1 && switchLR)
            reduce_L[node.vex] = false;
            num_L++;
        }
        else{
            reduce_R[node.vex] = false;
            num_R++;
        }

    }
    /// attention: the tsize of start node is smaller than filter size!
    if(i >= 0) {
        auto node = ordered_vex[pos_i][i];
        if((node.row) == switchLR) {
            //(node.row == 0 && !switchLR) || (node.row == 1 && switchLR)
            reduce_L[node.vex] = false;
            num_L++;
        }
        else{
            reduce_R[node.vex] = false;
            num_R++;
        }
    }


    // reduce Graph, save in tmpG
    idMapposL = new uint[Lsize];
    idMapposR = new uint[Rsize];
    posMapidL = new uint[Lsize];
    posMapidR = new uint[Rsize];
    visitL = new bool[Lsize];
    visitR = new bool[Rsize];
    for(uint t = 0; t < Lsize; t++) visitL[t] = false;
    for(uint t = 0; t < Rsize; t++) visitR[t] = false;

    // based on R
    vector<vector<uint>> &Rlist = g.adjR;
    resL = resR = 0;
    num_E = 0;
    for(i = 0; i < Rlist.size(); i++) {
        if(!reduce_R[i]) {
            uint vv = i;
            visitR[vv] = true;
            idMapposR[vv] = resR;
            posMapidR[resR++] = vv;

            vector<uint> vex;
            for(j = 0; j < Rlist[i].size(); j++) {
                uint uu = Rlist[i][j];
                if(!visitL[uu]) {
                    visitL[uu] = true;
                    idMapposL[uu] = resL;
                    posMapidL[resL++] = uu;
                }
                vex.emplace_back(idMapposL[uu]);
                num_E++;
            }

            tmpG.adjR.emplace_back(vex);
        }
    }
//cout<<"resR="<<resR<<", resL="<<resL<<endl;
//tmpG.dispalyAdj();

    vector<vector<uint>> &tmpR = tmpG.adjR;
    tmpG.adjL.resize(resL);
    for(i = 0; i < tmpR.size(); i++) {
        uint vvpos = i;
        for(j = 0; j < tmpR[i].size(); j++){
            uint uupos = tmpR[i][j];
            //ASSERT(uupos < tmpG.adjL.size())
            tmpG.adjL[uupos].emplace_back(vvpos);
        }
    }

    // sort adj and reduce graph
    // Reduce1hop(tmpG, t_U, t_V);
    tmpG.sortAdj();
//cout<<"sortAdj is ok"<<endl;


    tmpG.noE = num_E;
    tmpG.novL = resL;
    tmpG.novR = resR;
    tmpG.max_degree_L = tmpG.maxDegree(tmpG.adjL);
    tmpG.max_degree_R = tmpG.maxDegree(tmpG.adjR);
/*cout<<"------BCviz initial--------"<<endl;
cout<<"Information of tmpG:"<<endl;
cout<<"number of vertices in U: "<<tmpG.novL<<endl;
cout<<"number of vertices in V: "<<tmpG.novR<<endl;
cout<<"number of edges: "<<tmpG.noE<<endl;
cout<<"max degree of U: "<<tmpG.max_degree_L<<endl;
cout<<"max degree of V: "<<tmpG.max_degree_R<<endl<<endl;*/


    // initBiclique with v0
//cout<<"initBiclique: v0 = "<<v0<<", idMapposR[v0]="<<idMapposR[v0]<<endl;
    if(t_U == 3 && t_V == 3){
        initBiclique2Vertex(tmpG, idMapposR[v0]);
    }
    else{
        Reduce1hop(tmpG, t_U, t_V);
        if(tmpG.noE > 0) {
            initBiclique(tmpG);
        }
    }
//cout<<"initMEB time: "<<elapsed_time<<endl;

    delete[] reduce_L;
    delete[] reduce_R;

}

void ooMBE::getCandidateVertices(uint &size_component, uint &pos_l, uint &pos_r, uint interval_size,
    uint &num_L, uint &num_R, uint pos_i,
    bool *reduce_L, bool *reduce_R, const bool *core_L, const bool *core_R)
{
    uint l, r;
    // get the last vetices from right line
    for(r = pos_r; r < interval_size; r++) {
        auto &node = ordered_vex[pos_i][r];
        if( node.twoNodeSize < size_component ) {
            break;
        }
        if((node.row) == switchLR) {
            //(node.row == 0 && !switchLR) || (node.row == 1 && switchLR)
            if(!core_L[node.vex]) {
                reduce_L[node.vex] = false;
                num_L++;
            }

        }
        else if(!core_R[node.vex]){
            reduce_R[node.vex] = false;
            num_R++;
        }

    }
    pos_r = r-1;

    // get the last vetices from left line
    for(l = pos_l; l >= 0; l--) {
        auto node = ordered_vex[pos_i][l];
        if( node.twoNodeSize < size_component ) {
            break;
        }
        if((node.row) == switchLR) {
            //(node.row == 0 && !switchLR) || (node.row == 1 && switchLR)
            if(!core_L[node.vex]) {
                reduce_L[node.vex] = false;
                num_L++;
            }
        }
        else if(!core_R[node.vex]){
            reduce_R[node.vex] = false;
            num_R++;
        }

    }
    pos_l = l+1;

    /// attention: the tsize of start node is smaller than filter size!
    if(l >= 0) {
        auto node = ordered_vex[pos_i][l];
        if((node.row) == switchLR) {
            //(node.row == 0 && !switchLR) || (node.row == 1 && switchLR)
            if(!core_L[node.vex]) {
                reduce_L[node.vex] = false;
                num_L++;
            }

        }
        else if(!core_R[node.vex]){
            reduce_R[node.vex] = false;
            num_R++;
        }
    }

}

void ooMBE::init_BCviz_opt() {
//cout<<"initial feasible biclique"<<endl;
    g.sortAdj();
    uint Lsize = g.adjL.size();
    uint Rsize = g.adjR.size();
    bool *reduce_L = new bool[Lsize];
    bool *reduce_R = new bool[Rsize];
    bool *core_L = new bool[Lsize];
    bool *core_R = new bool[Rsize];
    for(uint i = 0; i < Lsize; i++) reduce_L[i] = core_L[i] = true;
    for(uint i = 0; i < Rsize; i++) reduce_R[i] = core_R[i] = true;
    for(uint i = 0; i < g.adjL.size(); i++)
        if(!g.adjL[i].empty()) core_L[i] = false;
    for(uint i = 0; i < g.adjR.size(); i++)
        if(!g.adjR[i].empty()) core_R[i] = false;


    uint v0, max_ts, size_component, num_L, num_R, resL, resR, num_E;
    uint pos_i, pos_j, pos_l, pos_r, i, j, k, l, r;

    v0 = max_ts = pos_i = pos_j = 0;

    // get v0 with largest local maximum cohesion in R
    for(i = 0; i < ordered_vex.size(); i++) {
        for(j = 0; j < ordered_vex[i].size(); j++) {
            auto node = ordered_vex[i][j];
            if( (node.row) != switchLR && node.twoNodeSize > max_ts) {
                //(node.row == 1 && !switchLR) || (node.row == 0 && switchLR)
                max_ts = node.twoNodeSize;
                v0 = node.vex;
                pos_i = i;
                pos_j = j;
            }
        }
    }
//cout<<"v0="<<v0<<", max_ts="<<max_ts<<endl;

    num_L = 0;
    num_R = 0;
    size_component = max_ts;
    pos_r = pos_j + 1;
    pos_l = pos_j;
    uint interval_size = ordered_vex[pos_i].size();
    l = pos_l, r = pos_r;
    uint smaller_size = size_component;
    // reduce graph: find candidate vertex set in BCviz
    while(BicliqueSize(num_L, num_R) < size_component || num_R < t_V || num_L < t_U) {
        if(pos_l==0 && pos_r==interval_size) break; /// new code
        size_component = smaller_size;
        // get the last vetices from right line
        for(r = pos_r; r < interval_size; r++) {
            auto &node = ordered_vex[pos_i][r];
            if( node.twoNodeSize < size_component ) {
                smaller_size = node.twoNodeSize;
                break;
            }
            if((node.row) == switchLR) {
                //(node.row == 0 && !switchLR) || (node.row == 1 && switchLR)
                if(!core_L[node.vex]) {
                    reduce_L[node.vex] = false;
                    num_L++;
                }

            }
            else if(!core_R[node.vex]){
                reduce_R[node.vex] = false;
                num_R++;
            }

        }
        pos_r = r;
        // get the last vetices from left line
        for(l = pos_l; l > 0; l--) {
            auto node = ordered_vex[pos_i][l];
            uint size = node.twoNodeSize;
            if( size < size_component ) {
                smaller_size = min(smaller_size, size);
                break;
            }
            if((node.row) == switchLR) {
                //(node.row == 0 && !switchLR) || (node.row == 1 && switchLR)
                if(!core_L[node.vex]) {
                    reduce_L[node.vex] = false;
                    num_L++;
                }
            }
            else if(!core_R[node.vex]){
                reduce_R[node.vex] = false;
                num_R++;
            }

        }
        pos_l = l;
        /// attention: the tsize of start node is smaller than filter size!
        if(l >= 0) {
            auto node = ordered_vex[pos_i][l];
            if((node.row) == switchLR) {
                //(node.row == 0 && !switchLR) || (node.row == 1 && switchLR)
                if(!core_L[node.vex]) {
                    reduce_L[node.vex] = false;
                    num_L++;
                }

            }
            else if(!core_R[node.vex]){
                reduce_R[node.vex] = false;
                num_R++;
            }
        }
    }
//cout<<"size_component="<<size_component<<", num_L="<<num_L<<", num_R="<<num_R<<endl;

    if(pos_r+t_U*t_V < interval_size) {
        size_component = ordered_vex[pos_i][pos_r+t_U*t_V].twoNodeSize;
//cout<<"another size_component="<<size_component<<endl;
        for(r = pos_r; r < interval_size; r++) {
            auto &node = ordered_vex[pos_i][r];
            if( node.twoNodeSize < size_component ) {
                break;
            }
            if((node.row) == switchLR) {
                //(node.row == 0 && !switchLR) || (node.row == 1 && switchLR)
                if(!core_L[node.vex]) {
                    reduce_L[node.vex] = false;
                    //num_L++;
                }

            }
            else if(!core_R[node.vex]){
                reduce_R[node.vex] = false;
                //num_R++;
            }

        }
        pos_r = r;
        for(l = pos_l; l > 0; l--) {
            auto node = ordered_vex[pos_i][l];
            uint size = node.twoNodeSize;
            if( size < size_component ) {
                break;
            }
            if((node.row) == switchLR) {
                //(node.row == 0 && !switchLR) || (node.row == 1 && switchLR)
                if(!core_L[node.vex]) {
                    reduce_L[node.vex] = false;
                    //num_L++;
                }
            }
            else if(!core_R[node.vex]){
                reduce_R[node.vex] = false;
                //num_R++;
            }

        }
        pos_l = l;
        if(l >= 0) {
            auto &node = ordered_vex[pos_i][l];
            if((node.row) == switchLR) {
                //(node.row == 0 && !switchLR) || (node.row == 1 && switchLR)
                if(!core_L[node.vex]) {
                    reduce_L[node.vex] = false;
                    //num_L++;
                }

            }
            else if(!core_R[node.vex]){
                reduce_R[node.vex] = false;
                //num_R++;
            }
        }
    }

    // reduce Graph, save in tmpG
    idMapposL = new uint[Lsize];
    idMapposR = new uint[Rsize];
    posMapidL = new uint[Lsize];
    posMapidR = new uint[Rsize];
    visitL = new bool[Lsize];
    visitR = new bool[Rsize];
    for(uint t = 0; t < Lsize; t++) visitL[t] = false;
    for(uint t = 0; t < Rsize; t++) visitR[t] = false;

    // get R
    vector<vector<uint>> &Rlist = g.adjR;
    resL = resR = 0;
    num_E = 0;
    for(i = 0; i < Rlist.size(); i++) {
        if(!reduce_R[i]) {
            uint vv = i;
            visitR[vv] = true;
            idMapposR[vv] = resR;
            posMapidR[resR++] = vv;

            vector<uint> vex;
            for(j = 0; j < Rlist[i].size(); j++) {
                uint uu = Rlist[i][j];
                if(!visitL[uu]) {
                    visitL[uu] = true;
                    idMapposL[uu] = resL;
                    posMapidL[resL++] = uu;
                }
                vex.emplace_back(idMapposL[uu]);
                num_E++;
            }

            tmpG.adjR.emplace_back(vex);
        }
    }

    // get L
    vector<vector<uint>> &tmpR = tmpG.adjR;
    tmpG.adjL.resize(resL);
    for(i = 0; i < tmpR.size(); i++) {
        uint vvpos = i;
        for(j = 0; j < tmpR[i].size(); j++){
            uint uupos = tmpR[i][j];
            //ASSERT(uupos < tmpG.adjL.size())
            tmpG.adjL[uupos].emplace_back(vvpos);
        }
    }

    // sort adj and reduce graph
    Reduce1hop(tmpG, t_U, t_V);
    tmpG.sortAdj();
    tmpG.noE = 0;
    tmpG.calculate_number_of_edges();
    tmpG.novL = resL;
    tmpG.novR = resR;
    tmpG.max_degree_L = tmpG.maxDegree(tmpG.adjL);
    tmpG.max_degree_R = tmpG.maxDegree(tmpG.adjR);
/*cout<<"------BCviz initial--------"<<endl;
cout<<"Information of tmpG:"<<endl;
cout<<"number of vertices in U: "<<tmpG.novL<<endl;
cout<<"number of vertices in V: "<<tmpG.novR<<endl;
cout<<"number of edges: "<<tmpG.noE<<endl;
cout<<"max degree of U: "<<tmpG.max_degree_L<<endl;
cout<<"max degree of V: "<<tmpG.max_degree_R<<endl<<endl;*/
    if(tmpG.noE > 0) initBiclique(tmpG);

    /// another step
    uint ite = 0;
    while(_C.MBsize==0 && ite < 10) {
        if(pos_r+t_U*t_V < interval_size) size_component = ordered_vex[pos_i][pos_r+t_U*t_V].twoNodeSize;
        else break;
        /// find more candidate vertex set in BCviz

        // get the last vetices from right line
        for(r = pos_r; r < interval_size; r++) {
            auto &node = ordered_vex[pos_i][r];
            if( node.twoNodeSize < size_component ) {
                break;
            }
            if((node.row) == switchLR) {
                //(node.row == 0 && !switchLR) || (node.row == 1 && switchLR)
                if(!core_L[node.vex]) {
                    reduce_L[node.vex] = false;
                    //num_L++;
                }

            }
            else if(!core_R[node.vex]){
                reduce_R[node.vex] = false;
                //num_R++;
            }

        }
        pos_r = r;

        // get the last vetices from left line
        for(l = pos_l; l > 0; l--) {
            auto &node = ordered_vex[pos_i][l];
            if( node.twoNodeSize < size_component ) {
                break;
            }
            if((node.row) == switchLR) {
                //(node.row == 0 && !switchLR) || (node.row == 1 && switchLR)
                if(!core_L[node.vex]) {
                    reduce_L[node.vex] = false;
                    //num_L++;
                }
            }
            else if(!core_R[node.vex]){
                reduce_R[node.vex] = false;
                //num_R++;
            }

        }
        pos_l = l;

        /// attention: the tsize of start node is smaller than filter size!
        if(l >= 0) {
            auto node = ordered_vex[pos_i][l];
            if((node.row) == switchLR) {
                //(node.row == 0 && !switchLR) || (node.row == 1 && switchLR)
                if(!core_L[node.vex]) {
                    reduce_L[node.vex] = false;
                    num_L++;
                }

            }
            else if(!core_R[node.vex]){
                reduce_R[node.vex] = false;
                num_R++;
            }
        }


        // get R
        //tmpG.adjR.resize(resR);
        for(i = 0; i < Rlist.size(); i++) {
            if(!reduce_R[i]) {
                uint vv = i;
                vector<uint> vex;
                for(j = 0; j < Rlist[i].size(); j++) {
                    uint uu = Rlist[i][j];
                    if(!visitL[uu]) {
                        visitL[uu] = true;
                        idMapposL[uu] = resL;
                        posMapidL[resL++] = uu;
                    }
                    vex.emplace_back(idMapposL[uu]);
                    //num_E++;
                }
                if(!visitR[vv]) {
                    visitR[vv] = true;
                    idMapposR[vv] = resR;
                    posMapidR[resR++] = vv;
                    tmpG.adjR.emplace_back(vex);
                }
                else
                    tmpG.adjR[idMapposR[vv]] = vex;
            }
        }

        // get L
        tmpG.adjL.clear();
        tmpG.adjL.resize(resL);
        for(i = 0; i < tmpG.adjR.size(); i++) {
            uint vvpos = i;
            for(j = 0; j < tmpR[i].size(); j++){
                uint uupos = tmpR[i][j];
                tmpG.adjL[uupos].emplace_back(vvpos);
            }
        }

        // find MEB
        Reduce1hop(tmpG, t_U, t_V);
        tmpG.sortAdj();
        tmpG.noE = 0;
        tmpG.calculate_number_of_edges();
        tmpG.novL = resL;
        tmpG.novR = resR;
        tmpG.max_degree_L = tmpG.maxDegree(tmpG.adjL);
        tmpG.max_degree_R = tmpG.maxDegree(tmpG.adjR);
        if(tmpG.noE > 0) initBiclique(tmpG);


        ite++;
    }

    if(_C.MBsize==0) initMBC();

    delete[] reduce_L;
    delete[] reduce_R;
    delete[] core_L;
    delete[] core_R;
}

void ooMBE::baseline_initMEB_BCviz() {
    g.sortAdj();
    uint Lsize = g.adjL.size();
    uint Rsize = g.adjR.size();

    // reduce Graph, save in tmpG
    idMapposL = new uint[Lsize];
    idMapposR = new uint[Rsize];
    posMapidL = new uint[Lsize];
    posMapidR = new uint[Rsize];
    visitL = new bool[Lsize];
    visitR = new bool[Rsize];
    for(uint t = 0; t < Lsize; t++) visitL[t] = false;
    for(uint t = 0; t < Rsize; t++) visitR[t] = false;

    initMBC();
}

void ooMBE::MEB_BCvizReduce(bgraph &Gk) {
//cout<<"************ReduceMEB***********"<<endl;
    // linear reduce in O(|E|)
    uint Lsize = Gk.adjL.size();
    uint Rsize = Gk.adjR.size();
    vector<vector<uint>> &Rlist = Gk.adjR;
    vector<vector<uint>> &Llist = Gk.adjL;
    bool *reduce_L = new bool[Lsize];
    bool *reduce_R = new bool[Rsize];
    for(uint i = 0; i < Lsize; i++) {reduce_L[i] = true;}
    for(uint i = 0; i < Rsize; i++) {reduce_R[i] = true;}

    uint i, j, k, l, num_L, num_R, pos_j, resL, resR, num_E;
    uint t_0V, t_1U, t_1V;
    double time1, time2, alltime1, alltime2;
    alltime1 = alltime2 = 0;
    uint ite=1; uint subit=1;
    for(i = 0; i < ordered_vex.size(); i++) {
        if((ordered_vex[i].size()/2) * (ordered_vex[i].size()-ordered_vex[i].size()/2) < _C.MBsize) continue;

//cout<<"ordered_vex coponent: "<<i<<endl;
        // todo: pruning - if the cohesion of peak vertex is small than _C.MBsize
        pos_j = 0;
        uint start_j = 0;
        vector<uint> tmpL, tmpR;
        uint max_coh = 0;
        do {

            // find a candidate component
            /// attention: the tsize of start node is smaller than filter size!
            //cout<<"++++++ num_group = "<<num_group<<", pos_j = "<<pos_j<<" ++++++"<<endl;
            for(j = pos_j; j < ordered_vex[i].size(); j++) {
                auto &node = ordered_vex[i][j];
                if(node.twoNodeSize > _C.MBsize) {
                    start_j = j; /// get the start node
                    break;
                }
            }
            if(j == ordered_vex[i].size()) {break;}

            /// get the start node
            vector<uint> canL, canR;
            if(start_j > 0) {
                auto &node = ordered_vex[i][--start_j];
                if((node.row) == switchLR) {
                    //(ordered_vex[i][start_j].row == 0 && !switchLR) || (ordered_vex[i][start_j].row == 1 && switchLR)
                    if(Llist[node.vex].size() > 0)
                        canL.emplace_back(node.vex);
                }
                else{
                    if(Rlist[node.vex].size() > 0)
                        canR.emplace_back(node.vex);
                }
            }

            // get candidate node
            while(j < ordered_vex[i].size() && ordered_vex[i][j].twoNodeSize > _C.MBsize) {
                auto &node = ordered_vex[i][j];
                if((node.row) == switchLR) {
                    //(node.row == 0 && !switchLR) || (node.row == 1 && switchLR)
                    if(Llist[node.vex].size() > 0)
                        canL.emplace_back(node.vex);
                }
                else {
                    if(Rlist[node.vex].size() > 0)
                        canR.emplace_back(node.vex);
                }
                j++;
            }

            pos_j = j;
            //num_group++;
            if(canL.size() >= t_U && canR.size() >= t_V && canL.size() * canR.size() > _C.MBsize) {
                for(auto it: canL) tmpL.push_back(it);
                for(auto it: canR) tmpR.push_back(it);
/*{
    /// output reduced graph
    for(auto v: canL) {reduce_L[v] = false;}
    for(auto v: canR) {reduce_R[v] = false;}
    uint numE = 0;
    for(uint uit = 0; uit < Llist.size(); uit++){
        if(!reduce_L[uit])
        for(auto vit: Llist[uit]) {
            if(!reduce_R[vit])
                numE++;
        }
    }
    cout<<"-----BCviz (subgraph "<<subit++<<")------"<<endl;
    cout<<"number of vertices in U: "<<canL.size()<<endl;
    cout<<"number of vertices in V: "<<canR.size()<<endl;
    cout<<"number of edges: "<<numE<<endl;
    for(auto v: canL) {reduce_L[v] = true;}
    for(auto v: canR) {reduce_R[v] = true;}
}*/

            }

        }while(pos_j < ordered_vex[i].size()); // a lot of component

        /// update _C using MBC_improved
        // todo: new code
        if(tmpL.size() >= t_U && tmpR.size() >= t_V && tmpL.size() * tmpR.size() > _C.MBsize) {

            /// 1. update tmpG
            for(auto &it : tmpL) reduce_L[it] = false;
            for(auto &it : tmpR) reduce_R[it] = false;

            // reduce Graph, save in tmpG
            // based on R, get tmpG.adjR

            for(auto &it : tmpG.adjL) it.clear();
            for(auto &it : tmpG.adjR) it.clear();
            resL = tmpG.adjL.size();
            resR = tmpG.adjR.size();
            num_E = 0;
//cout<<"init: resR="<<resR<<", resL="<<resL<<", num_E="<<num_E<<endl;
            for(k = 0; k < Rlist.size(); k++) {
                if(!reduce_R[k]) {
                    uint vv = k;

                    vector<uint> vex;
                    for(l = 0; l < Rlist[vv].size(); l++) {
                        uint uu = Rlist[vv][l];
                        if(!reduce_L[uu]) {
                            if(!visitL[uu]) {
                                visitL[uu] = true;
                                idMapposL[uu] = resL;
                                posMapidL[resL++] = uu;
                            }

                            vex.emplace_back(idMapposL[uu]);
                            num_E++;
                        }

                    }

                    if(!visitR[vv]) {
                        visitR[vv] = true;
                        idMapposR[vv] = resR;
                        posMapidR[resR++] = vv;
                        tmpG.adjR.emplace_back(vex);
                    }
                    else{

                        tmpG.adjR[idMapposR[vv]] = vex;
                    }
                }
            }
//cout<<"after: resR="<<resR<<", resL="<<resL<<", num_E="<<num_E<<endl;
//tmpG.dispalyAdj();

            if(num_E > _C.MBsize) {
                // get tmpG.adjL
                vector<vector<uint>> &tmpAdjR = tmpG.adjR;
                tmpG.adjL.resize(resL);
                for(k = 0; k < tmpAdjR.size(); k++) {
                    uint vvpos = k;
                    for(l = 0; l < tmpAdjR[vvpos].size(); l++){
                        uint uupos = tmpAdjR[vvpos][l];
                        //ASSERT(uupos < tmpG.adjL.size())
                        tmpG.adjL[uupos].emplace_back(vvpos);
                    }
                }

                // sort adj and reduce graph
                //Reduce1hop(tmpG, t_U, t_V);
                tmpG.sortAdj();

                //cout<<"Information of tmpG:"<<endl;
                tmpG.noE = num_E;
                tmpG.novL = resL;
                tmpG.novR = resR;
                tmpG.max_degree_L = tmpG.maxDegree(tmpG.adjL);
                tmpG.max_degree_R = tmpG.maxDegree(tmpG.adjR);
                if(tmpG.noE == 0) continue;
/*cout<<"-----BCviz (component "<<ite++<<")------"<<endl;
cout<<"number of vertices in U: "<<tmpG.novL<<endl;
cout<<"number of vertices in V: "<<tmpG.novR<<endl;
cout<<"number of edges: "<<tmpG.noE<<endl;
cout<<"max degree of U: "<<tmpG.max_degree_L<<endl;
cout<<"max degree of V: "<<tmpG.max_degree_R<<endl<<endl;*/

                /// 2. call MBC_improved
//cout<<"******MBC_improved******"<<endl;
                t_0V = g.max_degree_L;
                while(t_0V > t_V) {
                    t_1U = max(_C.MBsize/t_0V, (ulong)t_U);
                    t_1V = max(t_0V/2, t_V);  /// t_1V = max(t_0V-1, t_V);  faster method
//cout<<"MEBsize="<<_C.MBsize<<", t_1U="<<t_1U<<", t_1V="<<t_1V<<endl;
                    bgraph Gt = Reduce(tmpG, t_1U, t_1V); // one-hop and two-hop reduce
                    if(!Gt.adjL.empty()) MBC_Gk(Gt, t_1U, t_1V);
                    t_0V = t_1V;  /// t_0V = t_1V/2;
                }
            }

            /// 3. reback state
            for(auto &it : tmpL) reduce_L[it] = true;
            for(auto &it : tmpR) reduce_R[it] = true;

        }


    } // connected component

//cout<<"all reduce time: "<<alltime1<<" ms, all finding MEB time: "<<alltime2<<" ms."<<endl<<endl;


    delete[] reduce_L;
    delete[] reduce_R;

}

void ooMBE::MVB_BCvizReduce(bgraph &Gk) {
//cout<<"************ReduceBBVS_all***********"<<endl;
    // linear reduce in O(|E|)
    uint Lsize = Gk.adjL.size();
    uint Rsize = Gk.adjR.size();
    bool *reduce_L = new bool[Lsize];
    bool *reduce_R = new bool[Rsize];
    for(uint i = 0; i < Lsize; i++) reduce_L[i] = true;
    for(uint i = 0; i < Rsize; i++) reduce_R[i] = true;

    uint i, j, k, l, num_L, num_R, pos_j, resL, resR, num_E;
    uint t_0V, t_1U, t_1V;
    double time1, time2, alltime1, alltime2;
    alltime1 = alltime2 = 0;

    for(i = 0; i < ordered_vex.size(); i++) {
        if( problem_type == "MVB" && ordered_vex[i].size() < _C.MBsize) continue; // MVB
//cout<<"ordered_vex coponent: "<<i<<endl;

        pos_j = 0;
        //int num_group = 1;
        uint start_j = 0;
        vector<uint> tmpL, tmpR;
        do {

            // find a candidate component
            /// attention: the tsize of start node is smaller than filter size!
            //cout<<"++++++ num_group = "<<num_group<<", pos_j = "<<pos_j<<" ++++++"<<endl;
            for(j = pos_j; j < ordered_vex[i].size(); j++) {
                auto &node = ordered_vex[i][j];
                if(node.twoNodeSize > _C.MBsize) {
                    start_j = j; /// get the start node
                    break;
                }
            }
            if(j == ordered_vex[i].size()) {break;}

            /// get the start node
            if(start_j > 0) {
                auto &node = ordered_vex[i][--start_j];
                if((node.row) == switchLR) {
                    //(ordered_vex[i][start_j].row == 0 && !switchLR) || (ordered_vex[i][start_j].row == 1 && switchLR)
                    tmpL.emplace_back(node.vex);
                }
                else{
                    tmpR.emplace_back(node.vex);
                }
            }

            // get candidate node
            while(j < ordered_vex[i].size() && ordered_vex[i][j].twoNodeSize > _C.MBsize) {
                auto &node = ordered_vex[i][j];
                if((node.row) == switchLR) {
                    //(node.row == 0 && !switchLR) || (node.row == 1 && switchLR)
                    tmpL.emplace_back(node.vex);
                }
                else {
                    tmpR.emplace_back(node.vex);
                }
                j++;
            }

            pos_j = j;
            //num_group++;


        }while(pos_j < ordered_vex[i].size()); // a lot of component

        /// update _C using MBC_improved
        // todo: new code
        if(tmpL.size() > t_U && tmpR.size() > t_V && tmpL.size() + tmpR.size() > _C.MBsize) {

            /// 1. update tmpG
            for(auto &it : tmpL) reduce_L[it] = false;
            for(auto &it : tmpR) reduce_R[it] = false;

            // reduce Graph, save in tmpG
            // based on R, get tmpG.adjR
            vector<vector<uint>> &Rlist = Gk.adjR;
            for(auto &it : tmpG.adjL) it.clear();
            for(auto &it : tmpG.adjR) it.clear();
            resL = tmpG.adjL.size();
            resR = tmpG.adjR.size();
            num_E = 0;
//cout<<"init: resR="<<resR<<", resL="<<resL<<", num_E="<<num_E<<endl;
            for(k = 0; k < Rlist.size(); k++) {
                if(!reduce_R[k]) {
                    uint vv = k;

                    vector<uint> vex;
                    for(l = 0; l < Rlist[vv].size(); l++) {
                        uint uu = Rlist[vv][l];
                        if(!reduce_L[uu]) {
                            if(!visitL[uu]) {
                                visitL[uu] = true;
                                idMapposL[uu] = resL;
                                posMapidL[resL++] = uu;
                            }

                            vex.emplace_back(idMapposL[uu]);
                            num_E++;
                        }

                    }

                    if(!visitR[vv]) {
                        visitR[vv] = true;
                        idMapposR[vv] = resR;
                        posMapidR[resR++] = vv;
                        tmpG.adjR.emplace_back(vex);
                    }
                    else{

                        tmpG.adjR[idMapposR[vv]] = vex;
                    }
                }
            }
//cout<<"after: resR="<<resR<<", resL="<<resL<<", num_E="<<num_E<<endl;
//tmpG.dispalyAdj();

            // get tmpG.adjL
            vector<vector<uint>> &tmpAdjR = tmpG.adjR;
            tmpG.adjL.resize(resL);
            for(k = 0; k < tmpAdjR.size(); k++) {
                uint vvpos = k;
                for(l = 0; l < tmpAdjR[vvpos].size(); l++){
                    uint uupos = tmpAdjR[vvpos][l];
                    //ASSERT(uupos < tmpG.adjL.size())
                    tmpG.adjL[uupos].emplace_back(vvpos);
                }
            }

            // sort adj and reduce graph
            //Reduce1hop(tmpG, t_U, t_V);
            tmpG.sortAdj();

            //cout<<"Information of tmpG:"<<endl;
            tmpG.noE = num_E;
            tmpG.novL = resL;
            tmpG.novR = resR;
            tmpG.max_degree_L = tmpG.maxDegree(tmpG.adjL);
            tmpG.max_degree_R = tmpG.maxDegree(tmpG.adjR);
/*cout<<"number of vertices in U: "<<tmpG.novL<<endl;
cout<<"number of vertices in V: "<<tmpG.novR<<endl;
cout<<"number of edges: "<<tmpG.noE<<endl;
cout<<"max degree of U: "<<tmpG.max_degree_L<<endl;
cout<<"max degree of V: "<<tmpG.max_degree_R<<endl<<endl;*/

            /// 2. call MBC_improved
//cout<<"******MBC_improved******"<<endl;
            t_0V = g.max_degree_L;
            while(t_0V > t_V) {
                ulong tmpS = _C.MBsize > t_0V ? _C.MBsize - t_0V : 0;
                t_1U = max(tmpS, (ulong)t_U);
                t_1V = max(t_0V/2, t_V);  /// t_1V = max(t_0V-1, t_V);  faster method
//cout<<"MEBsize="<<_C.MBsize<<", t_1U="<<t_1U<<", t_1V="<<t_1V<<endl;
                bgraph Gt = Reduce(tmpG, t_1U, t_1V); // one-hop and two-hop reduce
                if(!Gt.adjL.empty()) MBC_Gk(Gt, t_1U, t_1V);
                t_0V = t_1V;  /// t_0V = t_1V/2;
            }

            /// 3. reback state
            for(auto &it : tmpL) reduce_L[it] = true;
            for(auto &it : tmpR) reduce_R[it] = true;

        }


    } // connected component in LSI

//cout<<"all reduce time: "<<alltime1<<" ms, all finding MEB time: "<<alltime2<<" ms."<<endl<<endl;


    delete[] reduce_L;
    delete[] reduce_R;

}

void ooMBE::MBB_BCvizReduce(bgraph &Gk) {
//cout<<"************ReduceBBVS_all***********"<<endl;
    // linear reduce in O(|E|)
    uint Lsize = Gk.adjL.size();
    uint Rsize = Gk.adjR.size();
    bool *reduce_L = new bool[Lsize];
    bool *reduce_R = new bool[Rsize];
    for(uint i = 0; i < Lsize; i++) reduce_L[i] = true;
    for(uint i = 0; i < Rsize; i++) reduce_R[i] = true;

    uint i, j, k, l, num_L, num_R, pos_j, resL, resR, num_E;
    uint t_0V, t_1U, t_1V;
    double time1, time2, alltime1, alltime2;
    alltime1 = alltime2 = 0;

    for(i = 0; i < ordered_vex.size(); i++) {
        if((ordered_vex[i].size()/2) * (ordered_vex[i].size()/2) < _C.MBsize) continue; // MBB
        // todo: pruning - max size p q
//cout<<"ordered_vex coponent: "<<i<<endl;

        pos_j = 0;
        //int num_group = 1;
        uint start_j = 0;
        vector<uint> tmpL, tmpR;
        do {

            // find a candidate component
            /// attention: the tsize of start node is smaller than filter size!
            //cout<<"++++++ num_group = "<<num_group<<", pos_j = "<<pos_j<<" ++++++"<<endl;
            for(j = pos_j; j < ordered_vex[i].size(); j++) {
                auto &node = ordered_vex[i][j];
                if(node.twoNodeSize > _C.MBsize) {
                    start_j = j; /// get the start node
                    break;
                }
            }
            if(j == ordered_vex[i].size()) {break;}

            /// get the start node
            if(start_j > 0) {
                auto &node = ordered_vex[i][--start_j];
                if((node.row) == switchLR) {
                    //(ordered_vex[i][start_j].row == 0 && !switchLR) || (ordered_vex[i][start_j].row == 1 && switchLR)
                    tmpL.emplace_back(node.vex);
                }
                else{
                    tmpR.emplace_back(node.vex);
                }
            }

            // get candidate node
            while(j < ordered_vex[i].size() && ordered_vex[i][j].twoNodeSize > _C.MBsize) {
                auto &node = ordered_vex[i][j];
                if((node.row) == switchLR) {
                    //(node.row == 0 && !switchLR) || (node.row == 1 && switchLR)
                    tmpL.emplace_back(node.vex);
                }
                else {
                    tmpR.emplace_back(node.vex);
                }
                j++;
            }

            pos_j = j;
            //num_group++;


        }while(pos_j < ordered_vex[i].size()); // a lot of component

        /// update _C using MBC_improved
        // todo: new code
        //if(min(tmpL.size(), tmpR.size()) > t_U && min(tmpL.size(), tmpR.size()) * min(tmpL.size(), tmpR.size()) > _C.MBsize) {
        if(tmpL.size() > t_U && tmpR.size() > t_V && min(tmpL.size(), tmpR.size()) * min(tmpL.size(), tmpR.size()) > _C.MBsize) {

            /// 1. update tmpG
            for(auto &it : tmpL) reduce_L[it] = false;
            for(auto &it : tmpR) reduce_R[it] = false;

            // reduce Graph, save in tmpG
            // based on R, get tmpG.adjR
            vector<vector<uint>> &Rlist = Gk.adjR;
            for(auto &it : tmpG.adjL) it.clear();
            for(auto &it : tmpG.adjR) it.clear();
            resL = tmpG.adjL.size();
            resR = tmpG.adjR.size();
            num_E = 0;
//cout<<"init: resR="<<resR<<", resL="<<resL<<", num_E="<<num_E<<endl;
            for(k = 0; k < Rlist.size(); k++) {
                if(!reduce_R[k]) {
                    uint vv = k;

                    vector<uint> vex;
                    for(l = 0; l < Rlist[vv].size(); l++) {
                        uint uu = Rlist[vv][l];
                        if(!reduce_L[uu]) {
                            if(!visitL[uu]) {
                                visitL[uu] = true;
                                idMapposL[uu] = resL;
                                posMapidL[resL++] = uu;
                            }

                            vex.emplace_back(idMapposL[uu]);
                            num_E++;
                        }

                    }

                    if(!visitR[vv]) {
                        visitR[vv] = true;
                        idMapposR[vv] = resR;
                        posMapidR[resR++] = vv;
                        tmpG.adjR.emplace_back(vex);
                    }
                    else{

                        tmpG.adjR[idMapposR[vv]] = vex;
                    }
                }
            }
//cout<<"after: resR="<<resR<<", resL="<<resL<<", num_E="<<num_E<<endl;
//tmpG.dispalyAdj();

            // get tmpG.adjL
            vector<vector<uint>> &tmpAdjR = tmpG.adjR;
            tmpG.adjL.resize(resL);
            for(k = 0; k < tmpAdjR.size(); k++) {
                uint vvpos = k;
                for(l = 0; l < tmpAdjR[vvpos].size(); l++){
                    uint uupos = tmpAdjR[vvpos][l];
                    //ASSERT(uupos < tmpG.adjL.size())
                    tmpG.adjL[uupos].emplace_back(vvpos);
                }
            }

            // sort adj and reduce graph
            //Reduce1hop(tmpG, t_U, t_V);
            tmpG.sortAdj();

            //cout<<"Information of tmpG:"<<endl;
            tmpG.noE = num_E;
            tmpG.novL = resL;
            tmpG.novR = resR;
            tmpG.max_degree_L = tmpG.maxDegree(tmpG.adjL);
            tmpG.max_degree_R = tmpG.maxDegree(tmpG.adjR);
/*cout<<"number of vertices in U: "<<tmpG.novL<<endl;
cout<<"number of vertices in V: "<<tmpG.novR<<endl;
cout<<"number of edges: "<<tmpG.noE<<endl;
cout<<"max degree of U: "<<tmpG.max_degree_L<<endl;
cout<<"max degree of V: "<<tmpG.max_degree_R<<endl<<endl;*/

            /// 2. call MBC_improved
//cout<<"******MBC_improved******"<<endl;
            uint mean1, tt;
            tt = max(t_U, t_V);
            if(_C.MBsize > 0 && tt < sqrt(_C.MBsize)) tt = sqrt(_C.MBsize);
            mean1 = tmpG.max_degree_L < tmpG.max_degree_R ? tmpG.max_degree_L : tmpG.max_degree_R;
            ulong last_MBB = _C.MBsize;
            while (mean1 > tt && last_MBB == _C.MBsize) {
                if(mean1 / tt > 10) mean1 = max(tt, mean1 / 2);
                else mean1 = max(tt, (mean1+tt) / 2);
//cout<<"mean="<<mean1<<", MBsize="<<_C.MBsize<<endl;
                bgraph Gt = Reduce(tmpG, mean1, mean1); // one-hop and two-hop reduce
                if(Gt.noE != 0) MBC_Gk(Gt, mean1, mean1);
            }

            /// 3. reback state
            for(auto &it : tmpL) reduce_L[it] = true;
            for(auto &it : tmpR) reduce_R[it] = true;

        }


    } // connected component in LSI

//cout<<"all reduce time: "<<alltime1<<" ms, all finding MEB time: "<<alltime2<<" ms."<<endl<<endl;


    delete[] reduce_L;
    delete[] reduce_R;

}

void ooMBE::initGraph() {
    Reduce1hop(g, t_U, t_V); // some vertices have been removed in BCviz
/*    vector<pair<uint, uint>> new_edges;
    for(uint u = 0; u < g.adjL.size(); u++) {
        for(auto v: g.adjL[u])
            new_edges.emplace_back(u,v);
    }*/

    g.max_degree_L = g.maxDegree(g.adjL);
    g.max_degree_R = g.maxDegree(g.adjR);
    g.noE = 0;
    g.calculate_number_of_edges();
//cout<<"After 1hopreduce, Information of g:"<<endl;
//g.dispalyAdj();

//cout<<"number of vertices in U: "<<g.novL<<endl;
//cout<<"number of vertices in V: "<<g.novR<<endl;
/*cout<<"number of edges: "<<g.noE<<endl;
cout<<"max degree of U: "<<g.max_degree_L<<endl;
cout<<"max degree of V: "<<g.max_degree_R<<endl<<endl;*/
}

void ooMBE::BCviz_MBC(string& BCvizfile) {
///cout<<"--------------read BCviz--------------"<<endl;
    auto start = high_resolution_clock::now();
    readBCviz(BCvizfile);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    uint originU = t_U;
    uint originV = t_V;
    if(problem_type == "MBB") {
        uint tt = max(t_U, t_V);
        t_U = tt;
        t_V = tt;
    }
    initGraph();

///cout<<"--------------prepare--------------"<<endl;
    start = high_resolution_clock::now();

    // init_BCviz();
    if(g.noE > 0) init_BCviz_opt();
    // baseline_initMEB_BCviz();

    stop = high_resolution_clock::now();
//cout<<"find a initial large biclique: size="<<_C.MBsize<<endl;
    //printMapMEB(_C);
    unsigned initCsize = _C.MBsize;
    duration = duration_cast<microseconds>(stop - start);
    double preptime = duration.count() * 1.0 / 1000.0;
//cout<<"Preparing and init MB time: "<<preptime<<" ms"<<endl<<endl;


///cout<<"--------------Reducing graph using BCviz and using MBC_improved to find MEB--------------"<<endl;
    start = high_resolution_clock::now();

    if(g.noE > 0) {
        if(problem_type == "MEB") MEB_BCvizReduce(g);
        if(problem_type == "MVB") MVB_BCvizReduce(g);
        if(problem_type == "MBB") MBB_BCvizReduce(g);
    }

    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    double reducetime = duration.count() * 1.0 / 1000.0;
//cout<<"find time: "<<reducetime<<" ms"<<endl;


    /** print result
     * [problem type] [algorithm] [dataset] [tu] [tv] [time (m seconds)] [size of MEB(p*q)] [init size of MBC]*/
    double exetime = reducetime; // remove initMBC time
    string dataset_name;
    const char split = '/';
    istringstream iss(g.datasets);
    string token;
    while (getline(iss, token, split)) dataset_name = token;

    //cout<<"------------MBC result-------------"<<endl;
    //printMapMEB(_C);
    /// output: [problem_type] [algo_type] [dataset] [s] [t] [search time]
    /// [size of maximum biclique] [size of initial biclique]
    cout<<problem_type<<" "<<method<<" "<<dataset_name<<" "<<originU<<" "<<originV<<" ";
    cout<<exetime<<endl;
    //cout<<"MB |U|="<<_C.subU.size()<<", |V|="<<_C.subV.size()<<endl<<endl;
}
/****** index search ******/



