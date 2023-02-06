#include "GrM.cpp"

#include "eigenHelper.cpp"
#include "cs.h"
#include <chrono>
#include<thread>
#include <Eigen/StdVector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include<vector>
#include <utility> 
#include <cmath>
#include <Eigen/SparseCore>
#include <unordered_map>
#include <unordered_set>
#include <set>
inline bool checkNeigb(Graph *query_graph, Graph *data_graph, vector<pair<ui, int>> QueryNlabel, VertexID vid)
{
    ui u_nbrs_countD;

    for (int i; i < QueryNlabel.size(); i++)
    {
        const VertexID *u_nbrsD = data_graph->getNeighborsByLabel(vid, QueryNlabel[i].first, u_nbrs_countD);
        if (u_nbrs_countD < QueryNlabel[i].second)
            return false;
    }
    return true;
}
void makeLD(SparseMatrix<double> &M)
{
    int countDD = 0;
    for (int ka = 0; ka < M.outerSize(); ++ka)
    {
        countDD = 0;
        for (SparseMatrix<double>::InnerIterator it(M, ka); it; ++it)
        {
            if (it.value() == -1)
                countDD++;
        }
        if (countDD > 0)
            M.insert(ka, ka) = countDD;

        countDD = 0;
    }
}
inline void removeVertexAndEgjes(vector<vector<CSV>> &FCS, CSV cvertex, int i, int deli)
{
    VertexID vx1;
    CSV cvertexpair;
    for (int j = 0; j < cvertex.edges.size(); j++)
    {
        // BSCHange
        vx1 = findIndBS(FCS, cvertex.edges[j].second, cvertex.edges[j].first);
        // vx1=findInd(FCS,cvertex.edges[j].second,cvertex.edges[j].first);
        cvertexpair = FCS[cvertex.edges[j].first][vx1];
        for (int k = 0; k < cvertexpair.edges.size(); k++)
        {

            if (cvertexpair.edges[k].first == i && cvertexpair.edges[k].second == cvertex.ID)
            {
                FCS[cvertex.edges[j].first][vx1].edges.erase(FCS[cvertex.edges[j].first][vx1].edges.begin() + k);
                break;
            }
        }
    }
    // while (!FCS[i][deli].edges.empty())
    //     FCS[i][deli].edges.pop_back();
    // or
    FCS[i][deli].edges.clear();
    FCS[i].erase(FCS[i].begin() + deli);
}
inline VertexID checkAN(vector<VertexID> ID, VertexID CID)
{
    for (int i = 0; i < ID.size(); i++)
        if (ID[i] == CID)
            return i;
    return 1000000;
}
bool comparisonFunction(const T &firstObject, const T &secondObject)
{
    if (firstObject.row() < secondObject.row())
        return true;
    if (firstObject.row() < secondObject.row())
        return false;
    if (firstObject.row() == secondObject.row())
        if (firstObject.col() < secondObject.col())
            return true;
    return false;
    // Here is an example how how you can compare.
    // return firstObject.firstParameter() < secondObject.firstParameter();
}
bool degcomp(Graph *data_graph, Graph *query_graph, int di, int qi)
{
    if (data_graph->getVertexDegree(di) >= query_graph->getVertexDegree(qi))
        return true;
    return false;
}
void getAllEdges(CSV &cvertex, vector<VertexID> &temp2, vector<pair<VertexID, VertexID>> &q_curr, vector<vector<CSV>> &FCS)
{
    // Add all edges of cvertex ID cvertex.edges + add if there exist a repetition
    // of the Edge ID from other candidate vertex. If it does add the extra edges
    // and remove also from q_curr the candidate vartex.
    CSV cvertex2;
    VertexID vx1;
    temp2.emplace_back(cvertex.edges[0].second);
    for (int i = 1; i < cvertex.edges.size(); i++)
    {
        vx1 = checkAN(temp2, cvertex.edges[i].second);
        if (vx1 == 1000000)
        {
            temp2.emplace_back(cvertex.edges[i].second);
        }
    }
    bool continue1 = true;
    // cout<<"first loop good "<< temp2.size()<<endl;
    int j = 0;
    VertexID test = 0;
    while (j < q_curr.size())
    {
        if (q_curr[j].second == cvertex.ID)
        {
            // cout<<"in baby"<<endl;
            // test=findInd(FCS,q_curr[j].second,q_curr[j].first);
            test = findIndBS(FCS, q_curr[j].second, q_curr[j].first);
            cvertex2 = FCS[q_curr[j].first][test];
            for (int i = 0; i < cvertex2.edges.size(); i++)
            {
                vx1 = checkAN(temp2, cvertex.edges[i].second);
                if (vx1 == 1000000)
                {
                    temp2.emplace_back(cvertex.edges[i].second);
                }
            }
            q_curr.erase(q_curr.begin() + j);
            // cout<<"inside list"<< q_curr.size()<<endl;
            // cout<<"in baby size"<<temp2.size()<<endl;
        }
        else
        {
            j++;
        }
    }
}
inline VertexID findInd(vector<vector<CSV>> &FCS, VertexID IDC, VertexID IDQ)
{
    // cout<<"s..."<<FCS[IDQ].size()<<endl;
    VertexID eleos;
    for (VertexID i = 0; i < FCS[IDQ].size(); i++)
    {
        if (FCS[IDQ][i].ID == IDC)
        {
            // cout<<"i..."<<i<<endl;
            eleos = i;
            return i;
        }
        // return i;
    }
    // cout<<"malaka..."<<eleos<<endl;
    // cout<<"pos ginete"<<endl;
    // cout<<FCS[IDQ].size()<<endl;

    return eleos;
}
inline bool allcomp(Graph *data_graph, Graph *query_graph, int di, int qi)
{
    return (degcomp(data_graph, query_graph, di, qi) && Labelcomp(data_graph, query_graph, di, qi));
}
bool Labelcomp(Graph *data_graph, Graph *query_graph, int di, int qi)
{
    if (data_graph->getVertexLabel(di) == query_graph->getVertexLabel(qi))
        return true;
    return false;
}
inline VertexID checkANSBS(VertexID *ID, VertexID CID, int sizA)
{
    int lo = 0, hi = sizA - 1;
    int mid;
    // This below check covers all cases , so need to check
    // for mid=lo-(hi-lo)/2
    while (hi - lo > 1)
    {
        int mid = (hi + lo) / 2;
        if (ID[mid] > CID)
        {
            lo = mid + 1;
        }
        else
        {
            hi = mid - 1;
        }
    }
    if (ID[lo] == CID)
    {
        return lo;
    }
    else if (ID[hi] == CID)
    {
        return hi;
    }
    // cout<<"error Prob"<<endl;
    return 1000000;
}
inline VertexID checkANS(VertexID *ID, VertexID CID, int sizA)
{
    for (int i = 0; i < sizA; i++)
        if (ID[i] == CID)
            return i;
    return 1000000;
}
inline bool OneHopEigenSFPN(CSV &cvertex, VertexID *ID, VertexID *IDD, int *IDDLC,
                            vector<T> &tripletList, vector<pair<ui, int>> &EvalNeigb, vector<pair<ui, int>> &EvalNeigb2, Graph *query_graph, int countD, vector<pair<VertexID, VertexID>> &q_curr)
{
    ID[0] = cvertex.ID;

    IDDLC[0] = 1;
    IDDLC[1] = 1;
    IDDLC[2] = 1;
    VertexID vx1 = 0;
    VertexID vx2 = 0;
    int count = 1;
    // cout<<"sh1"<<endl;
    for (int k = 0; k < cvertex.edges.size(); k++)
    {
        if (cvertex.edges[k].first == 1000000)
            continue;
        // vx1=checkAN(ID,cvertex.edges[k].second);
        // vx2=checkAN(IDD,cvertex.edges[k].first);
        vx1 = checkANS(ID, cvertex.edges[k].second, count);
        vx2 = checkANS(IDD, cvertex.edges[k].first, IDDLC[1]);

        // cout<<k<<"sh1R"<<endl;
        if (vx1 == 1000000)
        {

            // cout<<count<<"sh4.1"<<endl;
            ID[count] = cvertex.edges[k].second;

            tripletList.emplace_back(T(0, count, -1));
            tripletList.emplace_back(T(count, 0, -1));
            // testcor=true;
            count++;
            // cout<<"sh4.15"<<endl;
            // cout<<EvalNeigb[0].first<<endl;
            for (int t = 0; t < EvalNeigb.size(); t++)
            {
                // cout<<EvalNeigb[t].first<<endl;
                // cout<<cvertex.edges[k].first;
                // cout<<query_graph->getVertexLabel(cvertex.edges[k].first)<<endl;
                if (EvalNeigb[t].first == query_graph->getVertexLabel(cvertex.edges[k].first))
                {
                    EvalNeigb[t].second--;
                    EvalNeigb2[t].second--;
                    // testcor=false;
                    // break;
                }
            }
            // cout<<"sh4.2"<<endl;
            /*if(testcor){
                    cout<<"prob_mate"<<endl;
                    cout<<i<<","<<j<<","<<k<<endl;
                    cout<<query_graph->getVertexLabel(cvertex.edges[k].first)<<endl;
                }  */
        }
        if (vx2 == 1000000)
        {
            countD--;
            IDD[IDDLC[1]] = cvertex.edges[k].first;
            IDDLC[1]++;
            // cout<<"sh4.3"<<endl;
        }
        q_curr.emplace_back(cvertex.edges[k]);
    }
    IDDLC[0] = count;
    IDDLC[2] = count;
    // cout<<count<<"shF"<<endl;
    return (checkA(EvalNeigb) && checkB(countD));
}
void InitPrunCS(vector<vector<CSV>> &FCS, int qsiz, Graph *query_graph)
{
    int jj = 0;
    ui VDP;
    // degree -edges pruning-
    for (VertexID kk = 0; kk < qsiz; kk++)
    {
        jj = 0;
        VDP = query_graph->getVertexDegree(kk);
        while (jj < FCS[kk].size())
        {
            if (FCS[kk][jj].edges.size() == 0)
                FCS[kk].erase(FCS[kk].begin() + jj);
            else if (FCS[kk][jj].edges.size() < VDP)
            {
                VertexID i = 0;
                while (i < FCS[kk][jj].edges.size())
                {
                    VertexID rev = findIndBS(FCS, FCS[kk][jj].edges[i].second, FCS[kk][jj].edges[i].first); // vertex to remove ID?
                    FCS[FCS[kk][jj].edges[i].first][rev].Ichange = true;
                    for (int dd = 0; dd < FCS[FCS[kk][jj].edges[i].first][rev].edges.size(); dd++)
                        if (FCS[FCS[kk][jj].edges[i].first][rev].edges[dd].first == kk && FCS[FCS[kk][jj].edges[i].first][rev].edges[dd].second == FCS[kk][jj].ID)
                        {
                            FCS[FCS[kk][jj].edges[i].first][rev].edges.erase(FCS[FCS[kk][jj].edges[i].first][rev].edges.begin() + dd);
                            i++;
                            break;
                        }
                }
                FCS[kk].erase(FCS[kk].begin() + jj);
            }
            else
                jj++;
        }
    }
}
void edgesHelperCSN(CSV **FCSN, ui *count, const VertexID u_nbrs[], const VertexID u_nbrsD[], ui u_nbrs_count, ui u_nbrs_countD, int a, int b)
{
    for (VertexID c = 0; c < u_nbrs_count; c++)
    { // for every neigboor of the query
        // u_nbrs[c] gives the vertex ID of the queris current neigborhood

        for (VertexID d = 0; d < count[c]; d++)
        { // for every candidate of the query vertex
            int countx = 0;
            for (VertexID e = 0; e < u_nbrs_countD; e++)
            { // for every neigboor of the candidate vertex of the real graph
                if (u_nbrsD[e] == FCSN[u_nbrs[c]][d].ID)
                {
                    FCSN[a][b].edgesQ[countx] = u_nbrs[c];
                    FCSN[a][b].edgesV[countx] = u_nbrsD[c];
                    countx++;
                    // cout<<"hi"<<endl;
                    // FCS[a][b].edges.emplace_back(make_pair(u_nbrs[c],FCS[u_nbrs[c]][d].ID));
                }
                else if (u_nbrsD[e] > FCSN[u_nbrs[c]][d].ID)
                { // optimization
                    //;
                    FCSN[a][b].edgesV[countx + 1] = 10000;
                    FCSN[a][b].edgesQ[countx + 1] = 10000;
                    break;
                }
            }
            FCSN[a][b].edgesV[countx + 1] = 10000;
            FCSN[a][b].edgesQ[countx + 1] = 10000;
        }
    }
}
inline bool OneHopEigenNOTESFPN(CSV &cvertex, VertexID *ID, VertexID *IDD, int *IDDLC,
                                vector<pair<ui, int>> &EvalNeigb, vector<pair<ui, int>> &EvalNeigb2, Graph *query_graph, int countD, vector<pair<VertexID, VertexID>> &q_curr)
{
    ID[0] = cvertex.ID;

    IDDLC[0] = 1;
    IDDLC[1] = 1;
    VertexID vx1 = 0;
    VertexID vx2 = 0;
    int count = 1;
    for (int k = 0; k < cvertex.edges.size(); k++)
    {
        if (cvertex.edges[k].first == 1000000)
            continue;

        // vx1=checkANSBSADD(ID,cvertex.edges[k].second,count);

        // vx2=checkANSBSADD(IDD,cvertex.edges[k].first,IDDLC[1]);
        vx1 = checkANS(ID, cvertex.edges[k].second, count);

        vx2 = checkANS(IDD, cvertex.edges[k].first, IDDLC[1]);
        // cout<<cvertex.edges[k].second<<endl;
        // for (int d=1;d<=count;d++)
        //   if(ID[d-1]>ID[d-1])
        // cout<<ID[d]<<",";
        //   cout<<" "<<endl;
        if (vx1 == 1000000)
        {
            ID[count] = cvertex.edges[k].second;

            count++;
            // cout<<count<<endl;
            for (int t = 0; t < EvalNeigb.size(); t++)
            {
                if (EvalNeigb[t].first == query_graph->getVertexLabel(cvertex.edges[k].first))
                {
                    EvalNeigb[t].second--;
                    EvalNeigb2[t].second--;
                }
            }
        }
        if (vx2 == 1000000)
        {
            countD--;
            IDD[IDDLC[1]] = cvertex.edges[k].first;
            IDDLC[1]++;
        }
        q_curr.emplace_back(cvertex.edges[k]);
    }
    IDDLC[0] = count;
    // IDDLC[2]=count;

    return (checkA(EvalNeigb) && checkB(countD));
}
inline bool OneHopEigenPN(CSV &cvertex, vector<VertexID> &ID, vector<VertexID> &IDL, vector<VertexID> &IDD, vector<T> &tripletList, vector<pair<ui, int>> &EvalNeigb, vector<pair<ui, int>> &EvalNeigb2, Graph *query_graph, int countD, vector<pair<VertexID, VertexID>> &q_curr)
{
    ID.emplace_back(cvertex.ID);
    IDL.emplace_back(1);

    VertexID vx1;
    VertexID vx2;
    int count = 0;

    for (int k = 0; k < cvertex.edges.size(); k++)
    {
        vx1 = checkAN(ID, cvertex.edges[k].second);
        vx2 = checkAN(IDD, cvertex.edges[k].first);
        if (cvertex.edges[k].first == 1000000)
            continue;
        if (vx1 == 1000000)
        {
            count++;
            ID.emplace_back(cvertex.edges[k].second);
            IDL.emplace_back(1);
            tripletList.emplace_back(T(0, count, -1));
            tripletList.emplace_back(T(count, 0, -1));
            // testcor=true;
            for (int t = 0; t < EvalNeigb.size(); t++)
            {
                if (EvalNeigb[t].first == query_graph->getVertexLabel(cvertex.edges[k].first))
                {
                    EvalNeigb[t].second--;
                    EvalNeigb2[t].second--;
                    // testcor=false;
                    // break;
                }

            } /*if(testcor){
                     cout<<"prob_mate"<<endl;
                     cout<<i<<","<<j<<","<<k<<endl;
                     cout<<query_graph->getVertexLabel(cvertex.edges[k].first)<<endl;
                 }  */
        }
        if (vx2 == 1000000)
        {
            countD--;
            IDD.emplace_back(cvertex.edges[k].first);
        }
        q_curr.emplace_back(cvertex.edges[k]);
    }
    return (checkA(EvalNeigb) && checkB(countD));
}

inline bool OneHopEigenFPN(CSV &cvertex, vector<VertexID> &ID, vector<VertexID> &IDL, vector<VertexID> &IDD, vector<T> &tripletList, vector<pair<ui, int>> &EvalNeigb, vector<pair<ui, int>> &EvalNeigb2, Graph *query_graph, int countD, vector<pair<VertexID, VertexID>> &q_curr)
{
    ID.emplace_back(cvertex.ID);
    IDL.emplace_back(1);

    VertexID vx1;
    VertexID vx2;
    int count = 0;

    for (int k = 0; k < cvertex.edges.size(); k++)
    {
        if (cvertex.edges[k].first == 100000)
            continue;
        vx1 = checkAN(ID, cvertex.edges[k].second);
        vx2 = checkAN(IDD, cvertex.edges[k].first);
        if (vx1 == 1000000)
        {
            count++;
            ID.emplace_back(cvertex.edges[k].second);
            IDL.emplace_back(1);
            tripletList.emplace_back(T(0, count, -1));
            tripletList.emplace_back(T(count, 0, -1));
            // testcor=true;
            for (int t = 0; t < EvalNeigb.size(); t++)
            {
                if (EvalNeigb[t].first == query_graph->getVertexLabel(cvertex.edges[k].first))
                {
                    EvalNeigb[t].second--;
                    EvalNeigb2[t].second--;
                    // testcor=false;
                    // break;
                }

            } /*if(testcor){
                     cout<<"prob_mate"<<endl;
                     cout<<i<<","<<j<<","<<k<<endl;
                     cout<<query_graph->getVertexLabel(cvertex.edges[k].first)<<endl;
                 }  */
        }
        if (vx2 == 1000000)
        {
            countD--;
            IDD.emplace_back(cvertex.edges[k].first);
        }
        q_curr.emplace_back(cvertex.edges[k]);
    }
    return (checkA(EvalNeigb) && checkB(countD));
}

inline bool SecHopEigenNOTEPNSS(vector<pair<VertexID, VertexID>> &q_curr, unordered_set<int> &ID, unordered_set<int> &IDD, int *IDDLC, vector<vector<CSV>> &FCS, vector<pair<ui, int>> &EvalNeigb, Graph *query_graph, int countD)
{
    int count2 = 0;
    pair<VertexID, VertexID> temp1;
    VertexID tempxx = 0;
    CSV cvertexpair;
    vector<VertexID> temp2;
    bool vxD = 0;
    int count = IDDLC[0];
    bool checkNeigb = true;
    VertexID QN = 0;
    VertexID DN = 0;
    int EvalSize = EvalNeigb.size();
    int ti = q_curr.size();
    int telos = 0;
    // countD=countD-IDD.size();
    // while(!q_curr.empty()){
    while (telos < ti)
    {
        // cout<<telos<<endl;
        // count2=0;
        temp1 = q_curr[telos];
        telos++;
        // q_curr.pop_back();
        QN = temp1.first;
        // DN=temp1.second;
        if (QN == 1000000)
            continue;
        tempxx = findIndBS(FCS, temp1.second, QN);
        cvertexpair.ID = FCS[QN][tempxx].ID;
        // cvertexpair.edges.clear();
        // cvertexpair.edges=FCS[QN][tempxx].edges;
        // temp2.clear();
        for (int i = 0; i < FCS[QN][tempxx].edges.size(); i++)
        {
            if (FCS[QN][tempxx].edges[i].first == 1000000)
                continue;
            DN = FCS[QN][tempxx].edges[i].first;
            if (countD > 0)
            {
                if (IDD.insert(DN).second)
                {
                    //
                    // IDD.insert(DN);
                    // IDD[IDDLC[1]]=QN;
                    // IDDLC[1]++;
                    countD--;
                }
            }
            // cout<<"hm"<<endl;
            if (checkNeigb)
            {
                VertexID vtemp = FCS[QN][tempxx].edges[i].second;
                if ((ID.insert(vtemp).second))
                {
                    // if((ID.count(vtemp)==0)){

                    // ID.insert(vtemp);
                    count++;

                    for (int t = 0; t < EvalSize; t++)
                    {
                        if (EvalNeigb[t].first == query_graph->getVertexLabel(DN))
                        {
                            EvalNeigb[t].second--;
                            break;
                        }
                    }
                }
            }
        }
        if (checkNeigb == false || checkA(EvalNeigb))
        {
            checkNeigb = false;
            if (checkB(countD))
                return true;
        }
    }
    // temp2.clear();

    return (checkA(EvalNeigb) && checkB(countD));
}
inline bool SecHopEigenNOTEPNS(vector<pair<VertexID, VertexID>> &q_curr, VertexID *ID, VertexID *IDD, int *IDDLC, vector<vector<CSV>> &FCS, vector<pair<ui, int>> &EvalNeigb, Graph *query_graph, int countD)
{
    int count2 = 0;
    pair<VertexID, VertexID> temp1;
    VertexID tempxx = 0;
    CSV cvertexpair;
    vector<VertexID> temp2;
    VertexID vx1 = 0;
    VertexID vx2 = 0;
    VertexID vxD = 0;
    int count = IDDLC[0];

    VertexID QN = 0;
    VertexID DN = 0;
    int EvalSize = EvalNeigb.size();
    bool checkNeigb = true;
    // countD=countD-IDD.size();
    int ti = q_curr.size();
    int telos = 0;
    while (telos < ti)
    {
        // checkNeigb=true;
        // count2=0;
        // temp1=q_curr[q_curr.size()-1];
        temp1 = q_curr[telos];
        telos++;
        // q_curr.pop_back();
        QN = temp1.first;
        // DN=temp1.second;
        if (QN == 1000000)
            continue;
        tempxx = findIndBS(FCS, temp1.second, QN);
        // if(FCS[QN][tempxx].deleted)
        // cout<<"me paratas"<<endl;
        cvertexpair.ID = FCS[QN][tempxx].ID;
        // cvertexpair.edges.clear();
        // cvertexpair.edges=FCS[QN][tempxx].edges;
        // temp2.clear();

        // vx1=checkANS(ID,cvertexpair.ID,count);
        for (int i = 0; i < FCS[QN][tempxx].edges.size(); i++)
        {
            if (FCS[QN][tempxx].edges[i].first == 1000000)
                continue;
            DN = FCS[QN][tempxx].edges[i].first;
            if (countD > 0)
            {
                vxD = checkANS(IDD, DN, IDDLC[1]);
                // cout<<"you are an angel"<<endl;
                // vxD=checkANSBSADD(IDD,DN,IDDLC[1]);
                // cout<<"takes you back"<<endl;
                if (vxD == 1000000)
                {

                    IDD[IDDLC[1]] = DN;
                    IDDLC[1]++;
                    countD--;
                }
            }
            if (checkNeigb)
            {
                VertexID vtemp = FCS[QN][tempxx].edges[i].second;
                vx2 = checkANS(ID, vtemp, count);
                // vx2=checkANSBSADD(ID,vtemp,count);
                if (vx2 == 1000000)
                {

                    ID[count] = vtemp;
                    // vx2=count;
                    count++;

                    for (int t = 0; t < EvalSize; t++)
                    {
                        if (EvalNeigb[t].first == query_graph->getVertexLabel(DN))
                        {
                            EvalNeigb[t].second--;
                            break;
                        }

                    } // cout<<"tear from the side of my face"<<endl;
                }
            }
        }
        if (checkNeigb == false || checkA(EvalNeigb))
        {
            checkNeigb = false;
            if (checkB(countD))
                return true;
        }
    }
    // temp2.clear();
    IDDLC[0] = count;
    // IDDLC[2]=count;
    q_curr.clear();
    return (checkA(EvalNeigb) && checkB(countD));
}

inline bool SecHopEigenPNS(vector<pair<VertexID, VertexID>> &q_curr, VertexID *ID, VertexID *IDD, int *IDDLC, vector<vector<CSV>> &FCS, vector<T> &tripletList, vector<pair<ui, int>> &EvalNeigb, Graph *query_graph, int countD)
{
    int count2 = 0;
    pair<VertexID, VertexID> temp1;
    VertexID tempxx = 0;
    CSV cvertexpair;
    vector<VertexID> temp2;
    VertexID vx1 = 0;
    VertexID vx2 = 0;
    VertexID vxD = 0;
    int count = IDDLC[0];

    VertexID QN = 0;
    VertexID DN = 0;
    int EvalSize = EvalNeigb.size();
    // countD=countD-IDD.size();
    while (!q_curr.empty())
    {
        count2 = 0;
        temp1 = q_curr[q_curr.size() - 1];
        q_curr.pop_back();
        QN = temp1.first;
        // DN=temp1.second;
        if (QN == 1000000)
            continue;
        tempxx = findIndBS(FCS, temp1.second, QN);
        cvertexpair.ID = FCS[QN][tempxx].ID;
        cvertexpair.edges.clear();
        cvertexpair.edges = FCS[QN][tempxx].edges;

        temp2.clear();

        vx1 = checkANS(ID, cvertexpair.ID, count);
        // getAllEdges(cvertexpair,temp2,q_curr,FCS);
        for (int i = 0; i < cvertexpair.edges.size(); i++)
        {
            if (cvertexpair.edges[i].first == 1000000)
                continue;
            QN = cvertexpair.edges[i].first;
            vxD = checkANS(IDD, QN, IDDLC[1]);
            if (vxD == 1000000)
            {

                IDD[IDDLC[1]] = QN;
                IDDLC[1]++;
                countD--;
            }
            VertexID vtemp = cvertexpair.edges[i].second;
            vx2 = checkANS(ID, vtemp, count);
            if (vx2 == 1000000)
            {

                // ID.emplace_back(count);
                // ID.emplace_back(vtemp);
                ID[count] = vtemp;
                vx2 = count;
                count++;

                for (int t = 0; t < EvalSize; t++)
                {
                    if (EvalNeigb[t].first == query_graph->getVertexLabel(QN))
                    {
                        EvalNeigb[t].second--;
                        // testcor=false;
                        // break;
                    }
                }
                // q_temp.emplace_back(vx2);
                if (vx2 == vx1)
                    continue;

                tripletList.emplace_back(T(vx1, vx2, -1));
                tripletList.emplace_back(T(vx2, vx1, -1));
            }
            else
            {
                if (vx2 == vx1)
                {
                    cout << "hola" << endl;
                    continue;
                }

                tripletList.emplace_back(T(vx2, vx1, -1));
                tripletList.emplace_back(T(vx1, vx2, -1));
            }
        }
    }
    temp2.clear();
    IDDLC[0] = count;

    return (checkA(EvalNeigb) && checkB(countD));
}
inline void SecHopEigen(vector<pair<VertexID, VertexID>> &q_curr, vector<VertexID> &ID, vector<vector<CSV>> &FCS, vector<T> &tripletList)
{
    int count2 = 0;
    pair<VertexID, VertexID> temp1;
    VertexID tempxx = 0;
    CSV cvertexpair;
    vector<VertexID> temp2;
    VertexID vx1 = 0;
    VertexID vx2 = 0;
    int count = ID.size() - 1;
    while (!q_curr.empty())
    {
        count2 = 0;
        temp1 = q_curr[q_curr.size() - 1];

        tempxx = findIndBS(FCS, temp1.second, temp1.first);
        cvertexpair.ID = FCS[temp1.first][tempxx].ID;
        cvertexpair.edges.clear();
        cvertexpair.edges = FCS[temp1.first][tempxx].edges;
        temp2.clear();
        q_curr.pop_back();
        // getAllEdges(cvertexpair,temp2,q_curr,FCS);
        for (int i = 0; i < cvertexpair.edges.size(); i++)
        {
            if (cvertexpair.edges[i].first == 100000)
                continue;
            temp2.emplace_back(cvertexpair.edges[i].second);
        }
        vx1 = checkAN(ID, cvertexpair.ID);
        for (int e = 0; e < temp2.size(); e++)
        {
            VertexID vtemp = temp2[e];
            vx2 = checkAN(ID, vtemp);
            if (vx2 == 1000000)
            {
                count++;
                // ID.emplace_back(count);
                ID.emplace_back(vtemp);
                // IDL.emplace_back(1);
                vx2 = count;
                // q_temp.emplace_back(vx2);
                if (vx2 == vx1)
                    continue;
                tripletList.emplace_back(T(vx1, vx2, -1));
                tripletList.emplace_back(T(vx2, vx1, -1));
            }
            else
            {
                if (vx2 == vx1)
                {
                    cout << "hola" << endl;
                    continue;
                }

                // IDL[vx2]++;//can be count based.
                tripletList.emplace_back(T(vx2, vx1, -1));
                tripletList.emplace_back(T(vx1, vx2, -1));
            }
        }
    }
    temp2.clear();
}
inline bool SecHopEigenPNSS(vector<pair<VertexID, VertexID>> &q_curr, unordered_map<int, int> &ID, unordered_set<int> &IDD,
                            int *IDDLC, vector<vector<CSV>> &FCS, vector<T> &tripletList, vector<pair<ui, int>> &EvalNeigb, Graph *query_graph, int countD)
{
    int count2 = 0;
    pair<VertexID, VertexID> temp1;
    VertexID tempxx = 0;
    CSV cvertexpair;
    vector<VertexID> temp2;
    VertexID vx1 = 0;
    VertexID vx2 = 0;
    VertexID vxD = 0;
    int count = IDDLC[0];

    VertexID QN = 0;
    VertexID DN = 0;
    int EvalSize = EvalNeigb.size();
    // countD=countD-IDD.size();
    while (!q_curr.empty())
    {
        count2 = 0;
        temp1 = q_curr[q_curr.size() - 1];
        q_curr.pop_back();
        QN = temp1.first;
        // DN=temp1.second;
        if (QN == 1000000)
            continue;
        tempxx = findIndBS(FCS, temp1.second, QN);
        cvertexpair.ID = FCS[QN][tempxx].ID;
        cvertexpair.edges.clear();
        cvertexpair.edges = FCS[QN][tempxx].edges;

        temp2.clear();
        auto it = ID.find(cvertexpair.ID);
        vx1 = it->second;
        // vx1=checkANS(ID,cvertexpair.ID,count);
        // getAllEdges(cvertexpair,temp2,q_curr,FCS);
        for (int i = 0; i < cvertexpair.edges.size(); i++)
        {
            if (cvertexpair.edges[i].first == 1000000)
                continue;
            QN = cvertexpair.edges[i].first;
            if (countD > 0)
            {
                auto it3 = IDD.insert(QN);
                if (it3.second)
                {

                    // IDD[IDDLC[1]]=QN;
                    IDDLC[1]++;
                    countD--;
                }
            }
            VertexID vtemp = cvertexpair.edges[i].second;
            auto it1 = ID.insert({vtemp, count});
            // vx2=checkANS(ID,vtemp,count);
            if (it1.second)
            {

                // ID.emplace_back(count);
                // ID.emplace_back(vtemp);
                // ID[count]=vtemp;
                vx2 = count;
                count++;

                for (int t = 0; t < EvalSize; t++)
                {
                    if (EvalNeigb[t].first == query_graph->getVertexLabel(QN))
                    {
                        EvalNeigb[t].second--;
                        // testcor=false;
                        // break;
                    }
                }
                // q_temp.emplace_back(vx2);
                if (vx2 == vx1)
                    continue;

                tripletList.emplace_back(T(vx1, vx2, -1));
                tripletList.emplace_back(T(vx2, vx1, -1));
            }
            else
            {
                vx2 = it1.first->second;
                if (vx2 == vx1)
                {
                    cout << "hola" << endl;
                    continue;
                }

                tripletList.emplace_back(T(vx2, vx1, -1));
                tripletList.emplace_back(T(vx1, vx2, -1));
            }
        }
        if (checkB(countD) && checkA(EvalNeigb) && count > 101)
        {
            IDDLC[0] = count;
            return true;
        }
    }
    temp2.clear();
    IDDLC[0] = count;

    return (checkA(EvalNeigb) && checkB(countD));
}
inline VertexID checkANSBSADD(VertexID *&ID, VertexID CID, int sizA)
{
    int lo = 0, hi = sizA - 1;
    int mid;
    // This below check covers all cases , so need to check
    // for mid=lo-(hi-lo)/2
    while (hi - lo > 1)
    {
        mid = (hi + lo) / 2;
        if (ID[mid] == CID)
            return mid;
        if (ID[mid] > CID)
        {
            lo = mid;
        }
        else
        {
            hi = mid;
        }
    }
    if (ID[lo] == CID)
    {
        return lo;
    }
    else if (ID[hi] == CID)
    {
        return hi;
    }
    if (sizA == 1)
    {
        if (ID[0] < CID)
            ID[1] = CID;
        else
        {
            ID[1] = ID[0];
            ID[0] = CID;
        }
        return 1000000;
    }
    if (ID[sizA - 1] < CID)
    {
        ID[sizA] = CID;
        return 1000000;
    }
    if (ID[0] > CID)
    {
        for (int i = sizA; i > 0; i--)
            ID[i] = ID[i - 1];

        ID[0] = CID;
        return 1000000;
    }
    if (ID[hi] < CID)
        hi++;
    if (ID[lo] > CID)
        hi = lo;
    if (ID[lo - 1] > CID)
        hi--;
    for (int i = sizA; i > hi; i--)
        ID[i] = ID[i - 1];
    // cout<<"error Prob"<<endl;
    ID[hi] = CID;

    return 1000000;
}
inline bool eigcomp(VectorXd &qevalues, VectorXd &devalues, int comp)
{
    // return true;
    // comp=10;
    for (int i = 0; i < comp; i++)
        if (devalues[i] < qevalues[i])
            // if(qevalues[i]-devalues[i]>0.001&&(qevalues[i]>-1)&&(devalues[i]>-1))
            if (qevalues[i] - devalues[i] > 0.0001)
                return false;
    return true;
}
void printEigenValues(MatrixXd &eigenVD1, int sizq, int ec)
{
    for (int ko = 0; ko < sizq; ko++)
    {
        for (int jo = 0; jo < ec; jo++)
            cout << eigenVD1(ko, jo) << ",";
        cout << "VD: " << ko << endl;
    }
}
void printPrunRes(MatrixXd &eigenVD1, MatrixXd &eigenVq1, Graph *data_graph, Graph *query_graph, int sizq, int sizd)
{
    int prun = 0;
    bool pruned = false;
    int prunC = 0;
    bool prunedC = false;
    int Dprun = 0;
    bool pruned1 = false;
    int check = 0;
    int check1 = 0;
    bool pruned2;
    int VL = 0;
    int magkas = 0;

    for (int aaa = 1; aaa <= 10; aaa++)
    {
        check = 0;
        Dprun = 0;
        prunC = 0;
        prun = 0;
        check1 = 0;
        VL = 0;
        magkas = 0;
        pruned1 = false;
        prunedC = false;
        pruned = false;
        cout << aaa << "Eigen ";

        for (int i = 0; i < sizd; i++)
        {
            // cout<<"vertex "<<i<<" Eigs: ";
            for (int k = 0; k < sizq; k++)
            {
                pruned2 = true;
                prunedC = false;
                for (int j = 0; j < aaa; j++)
                {
                    // if(!(eigenVD1(i,j)<0||eigenVq1(k,j)<0))
                    if (eigenVD1(i, j) < eigenVq1(k, j))
                    {
                        pruned = true;
                        if (eigenVq1(k, j) - eigenVD1(i, j) > 0.0001)
                            prunedC = true;
                        if (eigenVD1(i, j) < -100000)
                        {
                            cout << "error";
                            cout << eigenVD1(i, j);
                            cout << i << "," << j << endl;
                            cout << data_graph->getVertexDegree(i) << endl;
                        }
                        // break;
                    }

                    if (eigenVD1(i, j) != eigenVq1(k, j))
                        pruned2 = false;
                    // cout<<eigenVq1(k,j)<<",";
                }
                if (pruned == true)
                {
                    prun++;
                    if (prunedC == true)
                        prunC++;
                }

                if (query_graph->getVertexDegree(k) > data_graph->getVertexDegree(i))
                    Dprun = Dprun + 1;
                if (pruned == true & query_graph->getVertexDegree(k) > data_graph->getVertexDegree(i))
                    check++;
                if (pruned == true || query_graph->getVertexDegree(k) > data_graph->getVertexDegree(i))
                    check1++;
                pruned = false;
                prunedC = false;
                if (query_graph->getVertexLabel(k) != data_graph->getVertexLabel(i))
                    VL++;
                if (pruned == true || query_graph->getVertexDegree(k) > data_graph->getVertexDegree(i) || query_graph->getVertexLabel(k) != data_graph->getVertexLabel(i))
                    magkas++;
            }

            // cout<<" "<<endl;
        }

        // cout<<"x"<<x<<endl;
        cout << "pruned" << (prun) << "prunedC" << prunC / (sizd * sizq) << endl;
        cout << "prunedD" << Dprun / (sizd * sizq) << endl;
        cout << "VertexL" << VL << endl;
        cout << "checkTomi" << check / (sizd * sizq) << endl;
        cout << "checkUnion " << check1 / (sizd * sizq) << endl;
        cout << "Possible" << sizd * sizq << endl;
        // cout<<"magkas?"<<magkas<<endl;
        // cout<<"Possible -Pruned"<<sizd*sizq<<"-"<<Dprun+prun-check<<endl;
        // cout<<"Remaining"<<(sizd*sizq)-magkas<<endl;
    }
}
inline bool SecHopEigenPN(vector<pair<VertexID, VertexID>> &q_curr, vector<VertexID> &ID, vector<VertexID> &IDD, vector<vector<CSV>> &FCS, vector<T> &tripletList, vector<pair<ui, int>> &EvalNeigb, Graph *query_graph, int countD)
{
    int count2 = 0;
    pair<VertexID, VertexID> temp1;
    VertexID tempxx = 0;
    CSV cvertexpair;
    vector<VertexID> temp2;
    VertexID vx1 = 0;
    VertexID vx2 = 0;
    VertexID vxD = 0;
    int count = ID.size() - 1;
    VertexID QN = 0;
    VertexID DN = 0;
    int EvalSize = EvalNeigb.size();
    // countD=countD-IDD.size();
    while (!q_curr.empty())
    {
        count2 = 0;
        temp1 = q_curr[q_curr.size() - 1];
        q_curr.pop_back();
        QN = temp1.first;
        // DN=temp1.second;
        if (QN == 1000000)
            continue;
        tempxx = findInd(FCS, temp1.second, QN);
        cvertexpair.ID = FCS[QN][tempxx].ID;
        cvertexpair.edges.clear();
        cvertexpair.edges = FCS[QN][tempxx].edges;
        temp2.clear();

        vx1 = checkAN(ID, cvertexpair.ID);
        // getAllEdges(cvertexpair,temp2,q_curr,FCS);
        for (int i = 0; i < cvertexpair.edges.size(); i++)
        {
            if (cvertexpair.edges[i].first == 1000000)
                continue;
            QN = cvertexpair.edges[i].first;
            vxD = checkAN(IDD, QN);
            if (vxD == 1000000)
            {
                IDD.emplace_back(QN);
                countD--;
            }
            VertexID vtemp = cvertexpair.edges[i].second;
            vx2 = checkAN(ID, vtemp);
            if (vx2 == 1000000)
            {
                count++;
                // ID.emplace_back(count);
                ID.emplace_back(vtemp);
                // IDL.emplace_back(1);
                vx2 = count;
                for (int t = 0; t < EvalSize; t++)
                {
                    if (EvalNeigb[t].first == query_graph->getVertexLabel(QN))
                    {
                        EvalNeigb[t].second--;
                        // testcor=false;
                        // break;
                    }
                }
                // q_temp.emplace_back(vx2);
                if (vx2 == vx1)
                    continue;

                tripletList.emplace_back(T(vx1, vx2, -1));
                tripletList.emplace_back(T(vx2, vx1, -1));
            }
            else
            {
                if (vx2 == vx1)
                {
                    cout << "hola" << endl;
                    continue;
                }

                // IDL[vx2]++;//can be count based.
                tripletList.emplace_back(T(vx2, vx1, -1));
                tripletList.emplace_back(T(vx1, vx2, -1));
            }
        }
    }
    temp2.clear();
    return (checkA(EvalNeigb) && checkB(countD));
}
inline bool SecHopEigenDPN(vector<pair<VertexID, VertexID>> &q_curr, vector<VertexID> &ID, vector<VertexID> &IDD, vector<vector<CSV>> &FCS, vector<T> &tripletList,
                           vector<pair<ui, int>> &EvalNeigb, Graph *query_graph, int countD, MatrixXd &M)
{
    int count2 = 0;
    pair<VertexID, VertexID> temp1;
    VertexID tempxx = 0;
    CSV cvertexpair;
    vector<VertexID> temp2;
    VertexID vx1 = 0;
    VertexID vx2 = 0;
    VertexID vxD = 0;
    int count = ID.size() - 1;
    VertexID QN = 0;
    VertexID DN = 0;
    int EvalSize = EvalNeigb.size();
    // countD=countD-IDD.size();
    while (!q_curr.empty())
    {
        count2 = 0;
        temp1 = q_curr[q_curr.size() - 1];
        QN = temp1.first;
        // DN=temp1.second;
        tempxx = findIndBS(FCS, temp1.second, QN);
        cvertexpair.ID = FCS[QN][tempxx].ID;
        cvertexpair.edges.clear();
        cvertexpair.edges = FCS[QN][tempxx].edges;
        temp2.clear();
        q_curr.pop_back();
        vx1 = checkAN(ID, cvertexpair.ID);
        // getAllEdges(cvertexpair,temp2,q_curr,FCS);
        for (int i = 0; i < cvertexpair.edges.size(); i++)
        {
            QN = cvertexpair.edges[i].first;
            vxD = checkAN(IDD, QN);
            if (vxD == 1000000)
            {
                IDD.emplace_back(QN);
                countD--;
            }
            VertexID vtemp = cvertexpair.edges[i].second;
            vx2 = checkAN(ID, vtemp);
            if (vx2 == 1000000)
            {
                count++;
                // ID.emplace_back(count);
                ID.emplace_back(vtemp);
                // IDL.emplace_back(1);
                vx2 = count;
                for (int t = 0; t < EvalSize; t++)
                {
                    if (EvalNeigb[t].first == query_graph->getVertexLabel(QN))
                    {
                        EvalNeigb[t].second--;
                        // testcor=false;
                        // break;
                    }
                }
                // q_temp.emplace_back(vx2);
                if (vx2 == vx1)
                    continue;
                M(vx1, vx2) = -1;
                M(vx2, vx1) = -1;
            }
            else
            {
                if (vx2 == vx1)
                {
                    cout << "hola" << endl;
                    continue;
                }

                // IDL[vx2]++;//can be count based.
                M(vx1, vx2) = -1;
                M(vx2, vx1) = -1;
            }
        }
    }
    temp2.clear();
    return (checkA(EvalNeigb) && checkB(countD));
}

inline bool OneHopEigen(CSV &cvertex, vector<VertexID> &ID, vector<VertexID> &IDL, vector<VertexID> &IDD, vector<T> &tripletList, vector<pair<ui, int>> &EvalNeigb, Graph *query_graph, int countD, vector<pair<VertexID, VertexID>> &q_curr)
{
    ID.emplace_back(cvertex.ID);
    IDL.emplace_back(1);

    VertexID vx1;
    VertexID vx2;
    int count = 0;

    for (int k = 0; k < cvertex.edges.size(); k++)
    {
        vx1 = checkAN(ID, cvertex.edges[k].second);
        vx2 = checkAN(IDD, cvertex.edges[k].first);
        if (vx1 == 1000000)
        {
            count++;
            ID.emplace_back(cvertex.edges[k].second);
            IDL.emplace_back(1);
            tripletList.emplace_back(T(0, count, -1));
            tripletList.emplace_back(T(count, 0, -1));
            // testcor=true;
            for (int t = 0; t < EvalNeigb.size(); t++)
            {
                if (EvalNeigb[t].first == query_graph->getVertexLabel(cvertex.edges[k].first))
                {
                    EvalNeigb[t].second--;
                    // testcor=false;
                    // break;
                }

            } /*if(testcor){
                     cout<<"prob_mate"<<endl;
                     cout<<i<<","<<j<<","<<k<<endl;
                     cout<<query_graph->getVertexLabel(cvertex.edges[k].first)<<endl;
                 }  */
        }
        if (vx2 == 1000000)
        {
            countD--;
            IDD.emplace_back(cvertex.edges[k].first);
        }
        q_curr.emplace_back(cvertex.edges[k]);
    }
    return (checkA(EvalNeigb) && checkB(countD));
}
void ExtractNI2H(vector<vector<pair<ui, int>>> &QueryNlabel, vector<vector<pair<ui, int>>> &QueryNlabel2, vector<int> &Dnew, Graph *query_graph, int qsiz)
{
    // If I want to calculate second hop I need to create a QUEUE for BFS +ID count as i will find repeated.
    const VertexID *u_nbrs;
    const VertexID *u_nbrs1;
    bool ad = true;
    int countkk = 0;
    ui u_nbrs_count;
    ui u_nbrs_count1;
    vector<VertexID> ID;
    vector<VertexID> q_curr;
    ui temp1;
    VertexID vx1;
    for (int i = 0; i < qsiz; i++)
    {
        ID.clear();
        ID.emplace_back(i);
        u_nbrs = query_graph->getVertexNeighbors(i, u_nbrs_count);
        vector<pair<ui, int>> QueryVec;
        vector<pair<ui, int>> QueryVec2;

        for (int j = 0; j < u_nbrs_count; j++)
        {
            ui labela = query_graph->getVertexLabel(u_nbrs[j]);
            q_curr.emplace_back(u_nbrs[j]);
            ID.emplace_back(u_nbrs[j]);
            ad = true;

            for (int k = 0; k < countkk; k++)
            {

                if (QueryVec[k].first == labela)
                {
                    QueryVec[k].second++;
                    QueryVec2[k].second++;
                    ad = false;
                }
            }
            if (ad)
            {
                QueryVec.emplace_back(make_pair(labela, 1));
                QueryVec2.emplace_back(make_pair(labela, 1));
                countkk++;
            }
        }
        QueryNlabel.emplace_back(QueryVec);
        while (!q_curr.empty())
        {
            temp1 = q_curr[q_curr.size() - 1];
            q_curr.pop_back();
            u_nbrs1 = query_graph->getVertexNeighbors(temp1, u_nbrs_count1);
            for (int b = 0; b < u_nbrs_count1; b++)
            {
                vx1 = checkAN(ID, u_nbrs1[b]);
                if (vx1 == 1000000)
                {
                    ID.emplace_back(u_nbrs1[b]);
                    ui labela = query_graph->getVertexLabel(u_nbrs1[b]);

                    // ID.emplace_back(u_nbrs1[b]);
                    ad = true;

                    for (int c = 0; c < countkk; c++)
                    {

                        if (QueryVec2[c].first == labela)
                        {
                            QueryVec2[c].second++;
                            ad = false;
                        }
                    }
                    if (ad)
                    {
                        QueryVec2.emplace_back(make_pair(labela, 1));
                        countkk++;
                    }
                }
            }
        }
        QueryNlabel2.emplace_back(QueryVec2);
        QueryVec2.clear();
        QueryVec.clear();
        Dnew.emplace_back((ID.size() - 1) - query_graph->getVertexDegree(i));
        ID.clear();
        countkk = 0;
    }
}
void PrunRes(Graph *query_graph, Graph *data_graph, MatrixXd &eigenVD1, MatrixXd &eigenVq1, int sizd, int sizq, int EN, float avgPruned[], float avgPrunedC[],
             float avgDprun[], float avgcheck[], float avgcheck1[], float avgVL[], float avgmagkas[],
             float avgDL[], float avgEL[])
{
    int x = 0;
    for (int aaa = 1; aaa <= EN; aaa++)
    {
        int check = 0;
        int Dprun = 0;
        int prunC = 0;
        int prun = 0;
        int check1 = 0;
        int VL = 0;
        int magkas = 0;
        bool pruned1 = false;
        bool prunedC = false;
        bool pruned = false;
        int DL = 0;
        int EL = 0;
        // cout<<aaa<<"Eigen ";

        for (int i = 0; i < sizd; i++)
        {
            // cout<<"vertex "<<i<<" Eigs: ";
            for (int k = 0; k < sizq; k++)
            {
                // pruned2=true;
                prunedC = false;
                for (int j = 0; j < aaa; j++)
                {
                    // if(!(eigenVD1(i,j)<0||eigenVq1(k,j)<0))
                    if (eigenVD1(i, j) < eigenVq1(k, j))
                    {
                        pruned = true;
                        if (eigenVq1(k, j) - eigenVD1(i, j) > 0.0001)
                            prunedC = true;
                        if (eigenVD1(i, j) < -100000)
                        {
                            cout << "error E";
                            cout << eigenVD1(i, j);
                            cout << i << "," << j << endl;
                            cout << data_graph->getVertexDegree(i) << endl;
                        }
                        // break;
                    }

                    if (pruned == true)
                    {
                        prun++;
                        if (prunedC == true)
                            prunC++;
                    }

                    if (query_graph->getVertexDegree(k) > data_graph->getVertexDegree(i))
                        Dprun = Dprun + 1;
                    if (prunedC == true & query_graph->getVertexDegree(k) > data_graph->getVertexDegree(i))
                        check++;
                    if (prunedC == true || query_graph->getVertexDegree(k) > data_graph->getVertexDegree(i))
                        check1++;

                    if (query_graph->getVertexLabel(k) != data_graph->getVertexLabel(i))
                        VL++;
                    if (prunedC == true || query_graph->getVertexLabel(k) != data_graph->getVertexLabel(i))
                        EL++;
                    if (query_graph->getVertexDegree(k) > data_graph->getVertexDegree(i) || query_graph->getVertexLabel(k) != data_graph->getVertexLabel(i))
                        DL++;
                    if (prunedC == true || query_graph->getVertexDegree(k) > data_graph->getVertexDegree(i) || query_graph->getVertexLabel(k) != data_graph->getVertexLabel(i))
                        magkas++;
                    pruned = false;
                    prunedC = false;
                    // CSV candidate=new CSV(aaa,i,query_graph->getVertexDegree)
                }

                // cout<<" "<<endl;
            }
            // cout<<"pruned"<<prun/(sizd*sizq)<<"prunedC"<<prunC/(sizd*sizq)<<endl;
            // cout<<"prunedD"<<Dprun/(sizd*sizq)<<endl;
            // cout<<"VertexL"<<VL/(sizd*sizq)<<endl;
            // cout<<"checkTomi"<<check/(sizd*sizq)<<endl;
            // cout<<"checkUnion "<<check1/(sizd*sizq)<<endl;
            // cout<<"Possible"<<sizd*sizq<<endl;
            x = sizd * sizq;

            avgPruned[aaa - 1] += prun / (x);
            avgPrunedC[aaa - 1] += prunC / (x);
            avgDprun[aaa - 1] += Dprun / (x);
            avgcheck[aaa - 1] += check / (x);
            avgcheck1[aaa - 1] += check1 / (x);
            avgmagkas[aaa - 1] += magkas / (x);
            avgVL[aaa - 1] += VL / (x);
            avgDL[aaa - 1] += DL / (x);
            avgEL[aaa - 1] += EL / (x);
            // cout<<"magkas?"<<magkas<<endl;
            // cout<<"Possible -Pruned"<<sizd*sizq<<"-"<<Dprun+prun-check<<endl;
            // cout<<"Remaining"<<(sizd*sizq)-magkas<<endl;
        }
    }
}
bool Refinement(vector<vector<pair<ui, int>>> QueryNlabel, vector<vector<pair<ui, int>>> QueryNlabel2, vector<int> &Qnew, vector<vector<CSV>> &FCS, int qsiz, Graph *query_graph, MatrixXd &eigenVq1)
{
    vector<T> tripletList;
    vector<VertexID> ID; // add size then resize
    vector<VertexID> IDD;
    vector<VertexID> IDL;
    VertexID vx1 = 0;
    VertexID vx2 = 0;
    VertexID vertexpair = 0;
    bool returnhere = false;
    VertexID vertex = 0;
    CSV cvertex;
    // CSV temp;
    vector<VertexID> temp2;
    vector<pair<VertexID, VertexID>> q_curr;

    vector<pair<ui, int>> EvalNeigb;
    vector<pair<ui, int>> EvalNeigb2;
    // bool eigencacl=true;
    VectorXd evalues;
    VectorXd qevalues;
    // pair<VertexID,VertexID> temp1;
    VertexID tempxx = 0;
    // CSV cvertexpair;
    // hop 1
    VertexID vertexDegree = 0;
    VertexID vertexlabel = 0;

    int prunA = 0;
    int prunAa = 0;
    int prunAb = 0;
    int prunB = 0;
    int prunC = 0;
    // test
    bool testcor;
    bool continueEL = false;
    for (int i = 0; i < qsiz; i++)
    {
        int j = 0;
        qevalues = eigenVq1.row(i);
        vertexDegree = query_graph->getVertexDegree(i);
        vertexlabel = query_graph->getVertexLabel(i);
        while (j < FCS[i].size())
        {
            continueEL = false;
            // EdgesCheck(FCS, qsiz);
            EvalNeigb.clear();
            EvalNeigb = QueryNlabel[i];
            EvalNeigb2 = QueryNlabel2[i];
            tripletList.clear();
            ID.clear();
            IDD.clear();
            IDL.clear();
            q_curr.clear();
            temp2.clear();

            cvertex = FCS[i][j];
            IDD.emplace_back(i);
            // bool continueE=OneHopEigen(cvertex,ID,IDL,IDD,tripletList,EvalNeigb,query_graph,vertexDegree,q_curr);
            bool continueE = OneHopEigenPN(cvertex, ID, IDL, IDD, tripletList, EvalNeigb, EvalNeigb2,
                                           query_graph, vertexDegree, q_curr);
            // BFS STEP1
            // go to hop2 and eigens
            if (continueE)
            {
                // SecHopEigen(q_curr,ID,FCS,tripletList);
                continueEL = SecHopEigenPN(q_curr, ID, IDD, FCS, tripletList, EvalNeigb2, query_graph, Qnew[i]);
                continueEL = true;
            }
            if (continueEL)
            {
                // cout<<"i"<<i<<"j"<<j<<endl;
                SparseMatrix<double> M(ID.size() + 1, ID.size() + 1);
                M.setFromTriplets(tripletList.begin(), tripletList.end(), [](double a, double b)
                                  { return b; });
                int countDD;
                for (int ka = 0; ka < M.outerSize(); ++ka)
                {
                    countDD = 0;
                    for (SparseMatrix<double>::InnerIterator it(M, ka); it; ++it)
                    {
                        if (it.value() == -1)
                            countDD++;
                    }
                    if (countDD > 0)
                        M.insert(ka, ka) = countDD;
                    countDD = 0;
                }
                calcEigens1(M, 10, evalues);

                /*Laplacian Check keep in comments
                if(checkM(M)){
                    EdgesCheck(FCS, qsiz);
                    cout<<"wrong Laplacian Matrix";
                    //printSM(M);
                    //if(i==0&&j==22)
                    //printSM(M);
                }
    */

                if (eigcomp(qevalues, evalues, 10) == false)
                {
                    prunC++;
                    removeVertexAndEgjes(FCS, cvertex, i, j);
                    returnhere = true;
                }
                else
                {
                    // cvertex.eigenvalues=evalues;
                    j++;
                }
            }
            else
            { // return either a degree or a neigborhood error.
                // if(checkA(EvalNeigb)==false)
                prunAa++;
                // if(checkB(countD)==false)
                prunAb++;
                // if(checkB(countD)==true && checkA(EvalNeigb)==false)
                prunB++;
                removeVertexAndEgjes(FCS, cvertex, i, j);
                returnhere = true;
                // q_curr.clear();
                prunA++;
            }
        }
    }
    // cout<<prunA<<"prunA"<<endl;
    // cout<<prunAa<<"prunAa"<<endl;
    // cout<<prunAb<<"prunAb"<<endl;
    // cout<<prunB<<"prunB"<<endl;
    // cout<<prunC<<"prunC"<<endl;
    // extract neigborhood level 1
    return returnhere;
}
bool ReverseRefinementOCS(VertexID *ID, VertexID *IDD, vector<vector<pair<ui, int>>> &QueryNlabel,
                          vector<vector<pair<ui, int>>> &QueryNlabel2, vector<int> &Qnew, vector<vector<CSV>> &FCS,
                          int qsiz, Graph *query_graph, MatrixXd &eigenVq1)
{
    vector<T> tripletList;
    tripletList.reserve(1000);
    int IDDLC[3] = {0, 0, 0};
    VertexID vx1 = 0;
    VertexID vx2 = 0;
    VertexID vertexpair = 0;
    bool returnhere = false;
    VertexID vertex = 0;
    CSV cvertex;
    vector<VertexID> temp2;
    vector<pair<VertexID, VertexID>> q_curr;
    vector<pair<ui, int>> EvalNeigb;
    vector<pair<ui, int>> EvalNeigb2;
    VectorXd evalues;
    VectorXd qevalues;
    VertexID tempxx = 0;
    int metrwiters = 0;
    // VectorXd eigens<<500,500,500,500,500,500,500;
    // hop 1
    VertexID vertexDegree = 0;
    VertexID vertexlabel = 0;
    int prunC = 0;
    // test
    bool testcor;
    bool continueEL = false;
    for (int i = qsiz - 1; i >= 0; i--)
    {

        int j = FCS[i].size();
        qevalues = eigenVq1.row(i);
        vertexDegree = query_graph->getVertexDegree(i);
        vertexlabel = query_graph->getVertexLabel(i);

        while (j > 0)
        {

            j--;
            cvertex = FCS[i][j];

            if (cvertex.deleted == true || cvertex.change == false)
                continue;
            continueEL = false;
            // metrwiters++;

            EvalNeigb.clear();
            EvalNeigb = QueryNlabel[i];
            EvalNeigb2 = QueryNlabel2[i];
            tripletList.clear();
            IDDLC[0] = 0;
            IDDLC[1] = 1;
            // IDDLC[2]=0;
            q_curr.clear();
            IDD[0] = i;
            // IDDLC[1]=1;
            // cout<<"A"<<endl;
            bool continueE = OneHopEigenSFPN(cvertex, ID, IDD, IDDLC, tripletList, EvalNeigb, EvalNeigb2,
                                             query_graph, vertexDegree, q_curr);
            // BFS STEP1
            // go to hop2 and eigens
            if (continueE)
            {
                // cout<<"B"<<endl;
                continueEL = SecHopEigenPNS(q_curr, ID, IDD, IDDLC, FCS, tripletList, EvalNeigb2, query_graph, Qnew[i]);
                continueEL = true;
            }
            if (continueEL)
            {
                // endtst
                // tripletList=(tripletList.begin(), tripletList.end(), [] (double a,double b));
                if (IDDLC[0] < 50)
                {

                    // cout<<"C"<<endl;

                    // triplet with sort etc
                    /*
                     sort(tripletList.begin(), tripletList.end(), tripletCompare);
            tripletList.erase(unique(tripletList.begin(), tripletList.end(), tripletUnique), tripletList.end());

            unordered_map<int, int> row_count;
            for (auto& t : tripletList) {
                row_count[t.row()]++;
            }
            int N = tripletList.size();

            for (auto &row_element : row_count) {
                tripletList.emplace_back(row_element.first, row_element.first, row_element.second);
            }
                SparseMatrix<double> M(IDDLC[0],IDDLC[0]);
          M.setFromTriplets(tripletList.begin(), tripletList.end());

          */
                    /*
                    SparseMatrix<double> M1(IDDLC[0],IDDLC[0]);

                    M1.setFromTriplets(tripletList.begin(), tripletList.end(), [] (double a,double b) { return b; });
                    //

                    int countDD;
                    for (int ka=0; ka<M1.outerSize(); ++ka){
                        countDD=0;
                    for (SparseMatrix<double>::InnerIterator it(M1,ka); it; ++it)
                        {if (it.value()==-1)
                            countDD++;
                     }if(countDD>0)
                        //M.insert(ka,ka) = countDD;
                    tripletList.emplace_back(T(ka,ka,countDD));

                     countDD=0;

          }
                            //makeLD(M);


            // create sparse matrix










          SparseMatrix<double> M(IDDLC[0],IDDLC[0]);
                 M.setFromTriplets(tripletList.begin(), tripletList.end(), [] (double a,double b) { return b; });
                */
                    // Matrix<double> M(ID.size()+1,ID.size()+1);
                    // checkM(M);

                    std::map<int, int> count_uniques;
                    for (auto t : tripletList)
                    {
                        count_uniques[t.row()]++;
                    }
                    // add the count to the triplets
                    for (auto it : count_uniques)
                    {
                        tripletList.push_back(Triplet<double>(it.first, it.first, it.second));
                    }
                    SparseMatrix<double> M(IDDLC[0], IDDLC[0]);
                    M.setFromTriplets(tripletList.begin(), tripletList.end(), [](double a, double b)
                                      { return b; });

                    calcEigens1(M, 10, evalues);
                    // eigcomp1(qevalues,evalues,10);
                    if (eigcomp(qevalues, evalues, 10) == false)
                    {
                        prunC++;
                        removeVertexAndEgjesFK(FCS, cvertex, i, j);
                        returnhere = true;
                    }
                    else
                    {
                        // cvertex.eigenvalues=evalues;
                        cvertex.change = false;
                        // j++;
                    }
                }
                else
                {
                    cvertex.change = false;
                }
            }

            else
            {
                removeVertexAndEgjesFK(FCS, cvertex, i, j);
                returnhere = true;
            }
        }
    } // cout<<"D"<<endl;
      // cout<<"metrw iterts"<<metrwiters<<endl;
    return returnhere;
}

bool ReverseRefinementOC(vector<vector<pair<ui, int>>> &QueryNlabel, vector<vector<pair<ui, int>>> &QueryNlabel2, vector<int> &Qnew, vector<vector<CSV>> &FCS, int qsiz, Graph *query_graph, MatrixXd &eigenVq1)
{
    vector<T> tripletList;
    tripletList.reserve(1000);
    vector<VertexID> ID; // add size then resize
    vector<VertexID> IDD;
    vector<VertexID> IDL;
    VertexID vx1 = 0;
    VertexID vx2 = 0;
    VertexID vertexpair = 0;
    bool returnhere = false;
    VertexID vertex = 0;
    CSV cvertex;
    vector<VertexID> temp2;
    vector<pair<VertexID, VertexID>> q_curr;
    vector<pair<ui, int>> EvalNeigb;
    vector<pair<ui, int>> EvalNeigb2;
    VectorXd evalues;
    VectorXd qevalues;
    VertexID tempxx = 0;

    // hop 1
    VertexID vertexDegree = 0;
    VertexID vertexlabel = 0;
    int prunC = 0;
    // test
    bool testcor;
    bool continueEL = false;
    for (int i = qsiz - 1; i >= 0; i--)
    {

        int j = FCS[i].size();
        qevalues = eigenVq1.row(i);
        vertexDegree = query_graph->getVertexDegree(i);
        vertexlabel = query_graph->getVertexLabel(i);
        while (j > 0)
        {
            j--;
            cvertex = FCS[i][j];
            if (cvertex.ID == 1000000 || cvertex.change == false)
                continue;
            continueEL = false;
            EvalNeigb.clear();
            EvalNeigb = QueryNlabel[i];
            EvalNeigb2 = QueryNlabel2[i];
            tripletList.clear();
            ID.clear();
            IDD.clear();
            IDL.clear();
            q_curr.clear();
            IDD.emplace_back(i);
            bool continueE = OneHopEigenFPN(cvertex, ID, IDL, IDD, tripletList, EvalNeigb, EvalNeigb2,
                                            query_graph, vertexDegree, q_curr);
            // BFS STEP1
            // go to hop2 and eigens
            if (continueE)
            {
                continueEL = SecHopEigenPN(q_curr, ID, IDD, FCS, tripletList, EvalNeigb2, query_graph, Qnew[i]);
                continueEL = true;
            }
            if (continueEL)
            {
                // tes

                // Accessing structure members using their
                // names.

                // endtst
                // tripletList=(tripletList.begin(), tripletList.end(), [] (double a,double b));
                SparseMatrix<double> M(ID.size(), ID.size());
                M.setFromTriplets(tripletList.begin(), tripletList.end(), [](double a, double b)
                                  { return b; });
                makeLD(M);
                /*int countDD;
                for (int ka=0; ka<M.outerSize(); ++ka){
                    countDD=0;
                for (SparseMatrix<double>::InnerIterator it(M,ka); it; ++it)
                    {if (it.value()==-1)
                        countDD++;
                 }if(countDD>0)
                    M.insert(ka,ka) = countDD;

                 countDD=0;

      }*/
                // Matrix<double> M(ID.size()+1,ID.size()+1);
                calcEigens1(M, 10, evalues);

                if (eigcomp(qevalues, evalues, 10) == false)
                {
                    prunC++;
                    removeVertexAndEgjesFK(FCS, cvertex, i, j);
                    returnhere = true;
                }
                else
                {
                    // cvertex.eigenvalues=evalues;
                    cvertex.change = false;
                    // j++;
                }
            }
            else
            {
                removeVertexAndEgjesFK(FCS, cvertex, i, j);
                returnhere = true;
            }
        }
    }

    return returnhere;
}
bool ReverseRefinement(vector<vector<pair<ui, int>>> &QueryNlabel, vector<vector<pair<ui, int>>> &QueryNlabel2, vector<int> &Qnew, vector<vector<CSV>> &FCS, int qsiz, Graph *query_graph, MatrixXd &eigenVq1)
{
    vector<T> tripletList;
    vector<VertexID> ID; // add size then resize
    vector<VertexID> IDD;
    vector<VertexID> IDL;
    VertexID vx1 = 0;
    VertexID vx2 = 0;
    VertexID vertexpair = 0;
    bool returnhere = false;
    VertexID vertex = 0;
    CSV cvertex;
    vector<VertexID> temp2;
    vector<pair<VertexID, VertexID>> q_curr;

    vector<pair<ui, int>> EvalNeigb;
    vector<pair<ui, int>> EvalNeigb2;
    // bool eigencacl=true;
    VectorXd evalues;
    VectorXd qevalues;
    // pair<VertexID,VertexID> temp1;
    VertexID tempxx = 0;
    // CSV cvertexpair;
    // hop 1
    VertexID vertexDegree = 0;
    VertexID vertexlabel = 0;

    int prunA = 0;
    int prunAa = 0;
    int prunAb = 0;
    int prunB = 0;
    int prunC = 0;
    // test
    bool testcor;
    bool continueEL = false;
    for (int i = qsiz - 1; i >= 0; i--)
    {

        int j = FCS[i].size();
        qevalues = eigenVq1.row(i);
        vertexDegree = query_graph->getVertexDegree(i);
        vertexlabel = query_graph->getVertexLabel(i);
        while (j > 0)
        {
            j--;
            cvertex = FCS[i][j];
            if (cvertex.ID == 1000000)
                continue;
            continueEL = false;
            EvalNeigb.clear();
            EvalNeigb = QueryNlabel[i];
            EvalNeigb2 = QueryNlabel2[i];
            tripletList.clear();
            ID.clear();
            IDD.clear();
            IDL.clear();
            q_curr.clear();
            temp2.clear();

            IDD.emplace_back(i);
            bool continueE = OneHopEigenFPN(cvertex, ID, IDL, IDD, tripletList, EvalNeigb, EvalNeigb2,
                                            query_graph, vertexDegree, q_curr);
            // BFS STEP1
            // go to hop2 and eigens
            if (continueE)
            {
                continueEL = SecHopEigenPN(q_curr, ID, IDD, FCS, tripletList, EvalNeigb2, query_graph, Qnew[i]);
                continueEL = true;
            }
            if (continueEL)
            {

                SparseMatrix<double> M(ID.size(), ID.size());
                M.setFromTriplets(tripletList.begin(), tripletList.end(), [](double a, double b)
                                  { return b; });
                int countDD;
                for (int ka = 0; ka < M.outerSize(); ++ka)
                {
                    countDD = 0;
                    for (SparseMatrix<double>::InnerIterator it(M, ka); it; ++it)
                    {
                        if (it.value() == -1)
                            countDD++;
                    }
                    if (countDD > 0)
                        M.insert(ka, ka) = countDD;

                    countDD = 0;
                }
                // Matrix<double> M(ID.size()+1,ID.size()+1);
                calcEigens1(M, 10, evalues);

                if (eigcomp(qevalues, evalues, 10) == false)
                {
                    prunC++;
                    removeVertexAndEgjesFK(FCS, cvertex, i, j);
                    returnhere = true;
                }
                else
                {
                    // cvertex.eigenvalues=evalues;
                    cvertex.change = false;
                    // j++;
                }
            }
            else
            {
                prunAa++;

                prunAb++;

                prunB++;
                removeVertexAndEgjesFK(FCS, cvertex, i, j);
                returnhere = true;

                prunA++;
            }
        }
    }

    return returnhere;
}
void MTEdgesCS(vector<vector<CSV>> &FCS,int qsiz,int dsiz,Graph *data_graph,Graph *query_graph){
 auto EHMT = [](vector<vector<CSV>> &FCS,int qsiz,int dsiz,Graph *data_graph,Graph *query_graph,int qstart,int qend) {
    ui u_nbrs_count=0;
    ui u_nbrs_countD=0;
    for (VertexID a=qstart;a<qend;a++){ //for every node of the query
   const VertexID* u_nbrs = query_graph->getVertexNeighbors(a, u_nbrs_count);
   // VertexID* u_nbrs = query_graph->getVertexNeighbors(a, u_nbrs_count);
   int siza=FCS[a].size()-1;
    for (VertexID b=0;b<siza;b++){ //for every candidate vertex of every node of the query
//candidate vertex for a query a is FCS[a][b].ID
    
    const VertexID* u_nbrsD = data_graph->getVertexNeighbors(FCS[a][b].ID, u_nbrs_countD); //real neigboors of the candidate vertex
     //VertexID* u_nbrsD = data_graph->getVertexNeighbors(FCS[a][b].ID, u_nbrs_countD);     
         //neigboors of query A node.
        //edgesHelperCS(FCS,u_nbrs,u_nbrsD,u_nbrs_count,u_nbrs_countD,a,b);
        VertexID vID=FCS[a][b].ID;
        edgesHelperCSPD(FCS,u_nbrs,u_nbrsD,u_nbrs_count,u_nbrs_countD,a,b,vID);
    }
}
    
    
    
 };
int Tnum=1;
int div=(int)(qsiz/Tnum);
thread th[Tnum];
int slotS[Tnum];
//const VertexID* u_nbrs;
 ui u_nbrs_count=0;
//const VertexID* u_nbrsD;
 ui u_nbrs_countD=0;
    int count=0;
    
    
    //for (int i=0;i<qsiz;i++){
     //   count=count+FCS[i].size();
    //}
    int range= int(count*0.1);
    int slotsX=int(count/Tnum);
    int pos=1;
    /*
    slotS[0]=0;
    slotS[1]=0;
    slotS[2]=0;
    slotS[3]=0;
    int qstart=0;
    int qend=0;
    int countN=0;
    cout<<"slots"<<slotsX<<endl;
    cout<<"Count"<<count<<endl;
    for(int i=0;i<qsiz;i++){
        countN=countN+FCS[i].size();
        if(countN>slotsX){
            slotS[pos]=i+1;
            countN=0;
            pos++;
            if(pos==Tnum)
            break;
        }
    }*/
    //for (int j=0;j<Tnum;j++){
     //   cout<<slotS[j]<<endl;
    //}


 for (int d=0;d<Tnum-1;d++){
    //cout<<div*d<<endl;
    //cout<<div*(d+1)<<endl;
    //qstart=slotS[d];
    //qend=slotS[d+1];
        th[d]=thread (EHMT,ref(FCS), qsiz, dsiz,data_graph,query_graph,div*d,div*(d+1));
        //th[d]=thread (EHMT,ref(FCS), qsiz, dsiz,data_graph,query_graph,qstart,qend);
    }
    th[Tnum-1]=thread (EHMT,ref(FCS), qsiz, dsiz,data_graph,query_graph, div*(Tnum-1),qsiz);
    //th[Tnum-1]=thread (EHMT,ref(FCS), qsiz, dsiz,data_graph,query_graph, qend,qsiz);
    //th[Tnum-1]=thread (AdJAdl,data_graph,degree,siz-2,siz,ref(eigenVD));
        for (int d=0;d<Tnum;d++)
                th[d].join();
                //cout<<"hi"<<endl;
        //std::cout << "Eigenvalues found:\n" <<  eigenVD.row(0) << std::endl;
 

 //ad edges to candidate space

}
void edgesHelperCSPD(vector<vector<CSV>> &FCS,const VertexID u_nbrs[],const VertexID u_nbrsD[],ui u_nbrs_count,ui u_nbrs_countD,int a,int b,VertexID VID){
for (VertexID c=0;c<u_nbrs_count;c++){ //for every neigboor of the query
            VertexID cne=u_nbrs[c]; //to cne gives the vertex ID of the queris current neigborhood 
        //for (const  VertexID cne:u_nbrs){
            //if (u_nbrs[c]<a)
            if (cne<a)
                continue;
            int sizqn=FCS[cne].size();
            for (VertexID d=0;d<sizqn;d++){ //for every candidate of the query vertex
                      VertexID de=FCS[cne][d].ID;
            //for (VertexID de:FCS[cne]){
                for (VertexID e=0;e<u_nbrs_countD;e++){//for every neigboor of the candidate vertex of the real graph
                    if(u_nbrsD[e]==de){
                        FCS[a][b].edges.emplace_back(make_pair(cne,de));
                        FCS[cne][d].edges.emplace_back(make_pair(a,VID));
                    }else if(u_nbrsD[e]>de){//optimization
                        //;
                        break;
                    }
                }
            }
        }
}
void edgesHelperCS(vector<vector<CSV>> &FCS,const VertexID u_nbrs[],const VertexID u_nbrsD[],ui u_nbrs_count,ui u_nbrs_countD,int a,int b){
for (VertexID c=0;c<u_nbrs_count;c++){ //for every neigboor of the query
            //u_nbrs[c] gives the vertex ID of the queris current neigborhood 
            
            for (VertexID d=0;d<FCS[u_nbrs[c]].size();d++){ //for every candidate of the query vertex
                       
                for (VertexID e=0;e<u_nbrs_countD;e++){//for every neigboor of the candidate vertex of the real graph
                    if(u_nbrsD[e]==FCS[u_nbrs[c]][d].ID){
                        FCS[a][b].edges.emplace_back(make_pair(u_nbrs[c],FCS[u_nbrs[c]][d].ID));
                    }else if(u_nbrsD[e]>FCS[u_nbrs[c]][d].ID){//optimization
                        //;
                        break;
                    }
                }
            }
        }
}
void testtriplets1(){
       vector<T> tripletList;
    SparseMatrix<double> M(120,120);
       tripletList.emplace_back(T(0,1,-1));tripletList.emplace_back(T(0,2,-1));
    tripletList.emplace_back(T(0,3,-1));tripletList.emplace_back(T(0,5,-1));
    tripletList.emplace_back(T(0,8,-1));tripletList.emplace_back(T(0,0,5));
    tripletList.emplace_back(T(1,0,-1));tripletList.emplace_back(T(2,0,-1));
    tripletList.emplace_back(T(3,0,-1));tripletList.emplace_back(T(5,0,-1));
    tripletList.emplace_back(T(8,0,-1));

    tripletList.emplace_back(T(3,1,-1));tripletList.emplace_back(T(8,1,-1));

    tripletList.emplace_back(T(1,3,-1));tripletList.emplace_back(T(1,8,-1));
    tripletList.emplace_back(T(1,1,3));

    tripletList.emplace_back(T(2,7,-1));tripletList.emplace_back(T(2,9,-1));
    tripletList.emplace_back(T(2,10,-1)); tripletList.emplace_back(T(2,2,4));

     tripletList.emplace_back(T(7,2,-1));tripletList.emplace_back(T(9,2,-1));
    tripletList.emplace_back(T(10,2,-1));
tripletList.emplace_back(T(4,3,-1));
    tripletList.emplace_back(T(3,4,-1));tripletList.emplace_back(T(3,3,3));

     tripletList.emplace_back(T(4,4,-1));

 tripletList.emplace_back(T(6,5,-1)); tripletList.emplace_back(T(7,5,-1));
      tripletList.emplace_back(T(8,5,-1)); tripletList.emplace_back(T(9,5,-1));
     tripletList.emplace_back(T(5,6,-1)); tripletList.emplace_back(T(5,7,-1));
      tripletList.emplace_back(T(5,8,-1)); tripletList.emplace_back(T(5,9,-1));
      tripletList.emplace_back(T(5,5,6));
tripletList.emplace_back(T(6,6,-1));tripletList.emplace_back(T(7,7,2));
tripletList.emplace_back(T(9,8,-1));tripletList.emplace_back(T(10,8,-1));
tripletList.emplace_back(T(8,9,-1));tripletList.emplace_back(T(8,10,-1));tripletList.emplace_back(T(8,8,4));
tripletList.emplace_back(T(9,9,2));tripletList.emplace_back(T(10,10,2));
    M.setFromTriplets(tripletList.begin(), tripletList.end(), [] (double a,double b) { return b; });
    VectorXd evalues;
    calcEigens1(M,100,evalues);
    for (int i=0;i<100;i++)
    cout<<evalues[i]<<endl;
}
void testtriplets(){
       vector<T> tripletList;
    SparseMatrix<double> M(12,12);
    tripletList.emplace_back(T(0,1,-1));tripletList.emplace_back(T(0,2,-1));
    tripletList.emplace_back(T(0,3,-1));tripletList.emplace_back(T(0,5,-1));
    tripletList.emplace_back(T(0,8,-1));tripletList.emplace_back(T(0,0,5));
    tripletList.emplace_back(T(1,0,-1));tripletList.emplace_back(T(2,0,-1));
    tripletList.emplace_back(T(3,0,-1));tripletList.emplace_back(T(5,0,-1));
    tripletList.emplace_back(T(8,0,-1));

    tripletList.emplace_back(T(3,1,-1));tripletList.emplace_back(T(8,1,-1));

    tripletList.emplace_back(T(1,3,-1));tripletList.emplace_back(T(1,8,-1));
    tripletList.emplace_back(T(1,1,3));

    tripletList.emplace_back(T(2,7,-1));tripletList.emplace_back(T(2,9,-1));
    tripletList.emplace_back(T(2,11,-1)); tripletList.emplace_back(T(2,2,4));

     tripletList.emplace_back(T(7,2,-1));tripletList.emplace_back(T(9,2,-1));
    tripletList.emplace_back(T(11,2,-1));
tripletList.emplace_back(T(4,3,-1));
    tripletList.emplace_back(T(3,4,-1));tripletList.emplace_back(T(3,3,3));

     tripletList.emplace_back(T(4,4,-1));

 tripletList.emplace_back(T(6,5,-1)); tripletList.emplace_back(T(7,5,-1));
      tripletList.emplace_back(T(8,5,-1)); tripletList.emplace_back(T(9,5,-1));
     tripletList.emplace_back(T(5,6,-1)); tripletList.emplace_back(T(5,7,-1));
      tripletList.emplace_back(T(5,8,-1)); tripletList.emplace_back(T(5,9,-1));
      tripletList.emplace_back(T(5,5,6));
tripletList.emplace_back(T(6,6,-1));tripletList.emplace_back(T(7,7,2));
tripletList.emplace_back(T(9,8,-1));tripletList.emplace_back(T(11,8,-1));
tripletList.emplace_back(T(8,9,-1));tripletList.emplace_back(T(8,11,-1));tripletList.emplace_back(T(8,8,4));
tripletList.emplace_back(T(9,9,2));tripletList.emplace_back(T(11,11,2));
    M.setFromTriplets(tripletList.begin(), tripletList.end(), [] (double a,double b) { return b; });
    VectorXd evalues;
    calcEigens1(M,10,evalues);
    for (int i=0;i<10;i++)
    cout<<evalues[i]<<endl;
}
void EdgesCS(vector<vector<CSV>> &FCS,int qsiz,int dsiz,Graph *data_graph,Graph *query_graph){
//const VertexID* u_nbrs;
 ui u_nbrs_count=0;
//const VertexID* u_nbrsD;
 ui u_nbrs_countD=0;
 //ad edges to candidate space
for (VertexID a=0;a<qsiz;a++){ //for every node of the query
   const VertexID* u_nbrs = query_graph->getVertexNeighbors(a, u_nbrs_count);
   // VertexID* u_nbrs = query_graph->getVertexNeighbors(a, u_nbrs_count);
    for (VertexID b=0;b<FCS[a].size();b++){ //for every candidate vertex of every node of the query
//candidate vertex for a query a is FCS[a][b].ID

    const VertexID* u_nbrsD = data_graph->getVertexNeighbors(FCS[a][b].ID, u_nbrs_countD); //real neigboors of the candidate vertex
     //VertexID* u_nbrsD = data_graph->getVertexNeighbors(FCS[a][b].ID, u_nbrs_countD);     
         //neigboors of query A node.
        edgesHelperCS(FCS,u_nbrs,u_nbrsD,u_nbrs_count,u_nbrs_countD,a,b);
    }
}
}
void EdgesCSBasicJ(vector<vector<CSV>> &FCS,int qsiz,int dsiz,Graph *data_graph,Graph *query_graph){
ui u_nbrs_count=0;
//const VertexID* u_nbrsD;
 ui u_nbrs_countD=0;
 //ad edges to candidate space
 int sizA=0;
 int sizC;
 VertexID VID=0;
for (VertexID a=0;a<qsiz;a++){ //for every node of the query
   const VertexID* u_nbrs = query_graph->getVertexNeighbors(a, u_nbrs_count);
   // VertexID* u_nbrs = query_graph->getVertexNeighbors(a, u_nbrs_count);
   sizA=FCS[a].size();
   for (VertexID c=0;c<u_nbrs_count;c++){
    VertexID cne=u_nbrs[c];
    if (cne<a)
        continue;
    sizC=FCS[cne].size();
    for (VertexID b=0;b<sizA;b++){ //for every candidate vertex of every node of the query
//candidate vertex for a query a is FCS[a][b].ID
    CSV VIA=FCS[a][b]; 
    VID=VIA.ID;
    for (VertexID d=0;d<sizC;d++){ //for every candidate of the query vertex
        CSV VIN=FCS[cne][d];
        if (data_graph->checkEdgeExistence(VIA.ID,VIN.ID)){
                    FCS[a][b].edges.emplace_back(make_pair(cne,VIN.ID));
                    FCS[cne][d].edges.emplace_back(make_pair(a,VID));}
                }
            }
        }
    }

}
void EdgesCSN(CSV **FCSN,ui* candidates_count,int qsiz,int dsiz,Graph *data_graph,Graph *query_graph){
//const VertexID* u_nbrs;
 ui u_nbrs_count=0;
//const VertexID* u_nbrsD;
 ui u_nbrs_countD=0;
 //ad edges to candidate space
for (VertexID a=0;a<qsiz;a++){ //for every node of the query
   const VertexID* u_nbrs = query_graph->getVertexNeighbors(a, u_nbrs_count);
   // VertexID* u_nbrs = query_graph->getVertexNeighbors(a, u_nbrs_count);
    for (VertexID b=0;b<candidates_count[a];b++){ //for every candidate vertex of every node of the query
//candidate vertex for a query a is FCS[a][b].ID

    const VertexID* u_nbrsD = data_graph->getVertexNeighbors(FCSN[a][b].ID, u_nbrs_countD); //real neigboors of the candidate vertex
     //VertexID* u_nbrsD = data_graph->getVertexNeighbors(FCS[a][b].ID, u_nbrs_countD);     
         //neigboors of query A node.
        edgesHelperCSN(FCSN,candidates_count,u_nbrs,u_nbrsD,u_nbrs_count,u_nbrs_countD,a,b);
    }
}
}
void EdgesCSNX(CSV **FCSN,ui* candidates_count,int qsiz,int dsiz,Graph *data_graph,Graph *query_graph){
//const VertexID* u_nbrs;
 ui u_nbrs_count=0;
//const VertexID* u_nbrsD;
 ui u_nbrs_countD=0;
 int countx=0;
  VertexID uneibC=0;
 //ad edges to candidate space
for (VertexID a=0;a<qsiz;a++){ //for every node of the query
   const VertexID* u_nbrs = query_graph->getVertexNeighbors(a, u_nbrs_count);
   // VertexID* u_nbrs = query_graph->getVertexNeighbors(a, u_nbrs_count);
    ui countSA=candidates_count[a];
    for (VertexID c=0;c<u_nbrs_count;c++){
        uneibC=u_nbrs[c];
        if (uneibC<a)
            continue;
        
        ui countSC=candidates_count[c];
       
    for (VertexID b=0;b<countSA;b++){ //for every candidate vertex of every node of the query
//candidate vertex for a query a is FCS[a][b].ID

    const VertexID* u_nbrsD = data_graph->getVertexNeighbors(FCSN[a][b].ID, u_nbrs_countD); //real neigboors of the candidate vertex
     //VertexID* u_nbrsD = data_graph->getVertexNeighbors(FCS[a][b].ID, u_nbrs_countD);     
         //neigboors of query A node.
         //for every neigboor of the query
            //u_nbrs[c] gives the vertex ID of the queris current neigborhood 
            
            for (VertexID d=0;d<countSC;d++){ //for every candidate of the query vertex
                        countx=0;
                        VertexID add=FCSN[uneibC][d].ID;
                for (VertexID e=0;e<u_nbrs_countD;e++){//for every neigboor of the candidate vertex of the real graph
                    if(u_nbrsD[e]==add){
                        FCSN[a][b].edgesQ[countx]=uneibC;
                        FCSN[a][b].edgesV[countx]=add;
                        //FCSN[b][a].edgesQ[countx]=uneibC; prob half?
                        //FCSN[b][a].edgesV[countx]=add;
                        countx++;
                        //cout<<"hi"<<endl;
                        //FCS[a][b].edges.emplace_back(make_pair(u_nbrs[c],FCS[u_nbrs[c]][d].ID));
                    }else if(u_nbrsD[e]>add){//optimization
                        //;
                        FCSN[a][b].edgesV[countx+1]=100000;
                        FCSN[a][b].edgesQ[countx+1]=100000;
                        break;
                    }
                }FCSN[a][b].edgesV[countx+1]=100000;
                FCSN[a][b].edgesQ[countx+1]=100000;

            }
        }
    }
}
}
void VerticesCS(vector<vector<CSV>> &FCS,int qsiz,int dsiz,Graph *data_graph,Graph *query_graph,MatrixXd &eigenVD1,MatrixXd &eigenVq1){
    VectorXd devalues;
    VectorXd qevalues;
    int maxdeg=data_graph->getGraphMaxDegree();
    int vc[qsiz];
    ui u_nbrs_count=0;
    int play=0;
    for (int i=0;i<qsiz;i++){
        play=0;
        const VertexID* u_nbrs = query_graph->getVertexNeighbors(i, u_nbrs_count);
        for (int j=0;j<u_nbrs_count;j++){
            play=play+FCS[u_nbrs[j]].size();
        }
        vc[i]=play;
        play=0;
    } 
    for (int i=0;i<qsiz;i++){
        qevalues=eigenVq1.row(i);
        int qmaxdeg=query_graph->getVertexDegree(i);
        vector<CSV> CS;
        for (int j=0;j<dsiz;j++){
            devalues=eigenVD1.row(j);
            
            if (allcomp(data_graph,query_graph,j,i) && eigcomp(qevalues,devalues,10)){//faster?
            //if (allcomp(data_graph,query_graph,j,i) && eigcomp(qevalues,devalues,10)){
                CSV cat(10,j,play);
                CS.emplace_back(cat);          
            }
        }
        FCS.emplace_back(CS);
        CS.clear();
    }
}

void VerticesCSNT(vector<CSV> FCS[],int qsiz,int dsiz,Graph *data_graph,Graph *query_graph,MatrixXd &eigenVD1,MatrixXd &eigenVq1,vector<vector<pair<ui,int>>> QueryNlabel){
 
     VectorXd devalues;
    VectorXd qevalues;
    ui com=data_graph->getGraphMaxLabelFrequency();
    LabelID label=0;
    vector<CSV> CS;
    int k=0;
    ui data_vertex_num;
    for (int i=0;i<qsiz;i++){
         label= query_graph->getVertexLabel(i);
        ui degree = query_graph->getVertexDegree(i);
        qevalues=eigenVq1.row(i);
        data_vertex_num=0;
        ui reserveS=com*degree;
        const ui* data_vertices = data_graph->getVerticesByLabel(label, data_vertex_num);
            for (ui j = 0; j < data_vertex_num; ++j) {           
            ui data_vertex = data_vertices[j];
            //devalues=eigenVD1.row(data_vertex);
                //
                if(data_graph->getVertexDegree(data_vertex) >= degree){
                 //if (data_graph->getVertexDegree(data_vertex) >= degree&&eigcomp(qevalues,devalues,10))  {
                    //if (data_graph->getVertexDegree(data_vertex) >= degree&&eigcomp(qevalues,devalues,10)&&checkNeigb(query_graph,data_graph,QueryNlabel[i],data_vertex))  {
                    devalues=eigenVD1.row(data_vertex);
                    if(eigcomp(qevalues,devalues,10)){
                        
                        //if(checkNeigb(query_graph,data_graph,QueryNlabel[i],data_vertex)){
                                

                         //checkNeigb function
                            ui u_nbrs_countD;
                        
                for (k=0;k<QueryNlabel[i].size()-1;k++){
                const VertexID* u_nbrsD = data_graph->getNeighborsByLabel(data_vertex,QueryNlabel[i][k].first, u_nbrs_countD);
                    if(u_nbrs_countD<QueryNlabel[i][k].second)
                            break;
                        }    
                    if (k==QueryNlabel[i].size()-1){
                    const VertexID* u_nbrsD = data_graph->getNeighborsByLabel(data_vertex,QueryNlabel[i][k].first, u_nbrs_countD);
                    if(u_nbrs_countD>=QueryNlabel[i][k].second){
                         //until here 
                         CSV cat(10,data_vertex,reserveS);
                        CS.emplace_back(cat);      
                        }
                        }
                         //CSV cat(10,data_vertex,reserveS);
                        //CS.emplace_back(cat); 
                    //}
                    }
                    }
                //CSV cat(10,data_vertex,reserveS);
                //CS.emplace_back(cat);          
            //}
        }
        FCS[i]=CS;
        CS.clear();
    }
}
void MTVerticesCS(vector<vector<CSV>> &FCS,int qsiz,int dsiz,Graph *data_graph,Graph *query_graph,MatrixXd &eigenVD1,MatrixXd &eigenVq1){
auto VHMT = [](vector<vector<CSV>> &FCS,int qsiz,int dsiz,Graph *data_graph,Graph *query_graph,MatrixXd &eigenVD1,MatrixXd &eigenVq1,int qstart,int qend) {
    VectorXd devalues;
    VectorXd qevalues;
    for (int i=qstart;i<qend;i++){
        qevalues=eigenVq1.row(i);
        
        vector<CSV> CS;
        for (int j=0;j<dsiz;j++){
            devalues=eigenVD1.row(j);
            if (allcomp(data_graph,query_graph,j,i) && eigcomp(qevalues,devalues,10)){
                CSV cat(10,j);
                CS.emplace_back(cat);          
            }
        }
        FCS.emplace_back(CS);
        CS.clear();
    }
 };

 int Tnum=4;
int div=(int)(qsiz/Tnum);
thread th[Tnum];

//const VertexID* u_nbrs;
 ui u_nbrs_count=0;
//const VertexID* u_nbrsD;
 ui u_nbrs_countD=0;

 for (int d=0;d<Tnum-1;d++){
        th[d]=thread (VHMT,ref(FCS), qsiz, dsiz,data_graph,query_graph,ref(eigenVD1),ref(eigenVq1),div*d,div*(d+1));
    }
    th[Tnum-1]=thread (VHMT,ref(FCS), qsiz, dsiz,data_graph,query_graph,ref(eigenVD1),ref(eigenVq1),div*(Tnum-1),qsiz);
    //th[Tnum-1]=thread (AdJAdl,data_graph,degree,siz-2,siz,ref(eigenVD));
        for (int d=0;d<Tnum;d++)
                th[d].join();
}
void VerticesCSMN(ui **&candidates, ui *&candidates_count,int qsiz,int dsiz,Graph *data_graph,Graph *query_graph,MatrixXd &eigenVD1,MatrixXd &eigenVq1){

    int count=0;
    VectorXd devalues;
    VectorXd qevalues;
    //int maxdeg=data_graph->getGraphMaxDegree();
    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        count=0;
        LabelID label = query_graph->getVertexLabel(i);
        ui degree = query_graph->getVertexDegree(i);
        qevalues=eigenVq1.row(i);
        ui data_vertex_num;
        const ui* data_vertices = data_graph->getVerticesByLabel(label, data_vertex_num);

        for (ui j = 0; j < data_vertex_num; ++j) {
            ui data_vertex = data_vertices[j];
            devalues=eigenVD1.row(data_vertex);
            if (data_graph->getVertexDegree(data_vertex) >= degree&&eigcomp(qevalues,devalues,10)) {
                candidates[i][count] = data_vertex;
                count++;
            }
        }candidates_count[i] =count;
    }
}
void VerticesCSM(ui **&candidates, ui *&candidates_count,int qsiz,int dsiz,Graph *data_graph,Graph *query_graph,MatrixXd &eigenVD1,MatrixXd &eigenVq1){
    VectorXd devalues;
    VectorXd qevalues;
    int count=0;
    //int maxdeg=data_graph->getGraphMaxDegree();
    for (int i=0;i<qsiz;i++){
        count=0;
        qevalues=eigenVq1.row(i);
        int qmaxdeg=query_graph->getVertexDegree(i);
        //vector<CSV> CS;
        for (int j=0;j<dsiz;j++){
            devalues=eigenVD1.row(j);
            
            if (allcomp(data_graph,query_graph,j,i) && eigcomp(qevalues,devalues,10)){//faster?
            //if (allcomp(data_graph,query_graph,j,i) && eigcomp(qevalues,devalues,10)){   
            candidates[i][count]=j;
            count++;     
            }
        }candidates_count[i]=count;

    }
}
void calcres(int sizd,MatrixXd eigenVD1,Graph *data_graph){
    string input_query_graph_file;
    //string input_qpartA="dataset\\wordnet\\query_graph\\query_dense_";
    string input_qpartA="dataset\\yeast\\query_graph\\query_dense_";
    //input_qpartA="test_case_1.graph";
    string input_qpartW="hprd_dense_4_L";
     //input_qpartW="res\\wordnet\\L\\wordnet_dense_";
      input_qpartW="res\\yeas\\L\\Testyeast_dense_";
    string input_end=".graph";
    string txtend="E.txt";
    ofstream myfile;
    
    
    string qr [5]={"4_","8_","16_","24_","32_"};
    //string qr [5]={"4_","8_","12_","16_","20_"};
    int prun=0;bool pruned=false;
    int prunC=0;bool prunedC=false;int Dprun=0;
    bool pruned1=false;int check=0;int check1=0;bool pruned2;int VL=0;
    int magkas=0;int DL=0;int EL=0;
    float avgPruned [10];float avgPrunedC[10] ;
    float avgDprun [10];float avgcheck[10];
    float avgcheck1 [10];float avgVL [10];float avgmagkas [10];
    float avgDL [10];float avgEL [10];
    int testeigen=10;
    for (int xx=1;xx<2;xx++){
    myfile.open (input_qpartW+(qr[xx])+txtend);
    for (int avg=0;avg<10;avg++){
            avgPruned [avg]=0;
            avgPrunedC[avg]=0;
            avgDprun [avg]=0;
            avgcheck[avg]=0;
            avgcheck1 [avg]=0;
            avgVL [avg]=0;
            avgmagkas [avg]=0;
            avgDL [avg]=0;
            avgEL [avg]=0;
        }    myfile <<"Prunned";
            myfile <<";";
            myfile <<"PrunnedT";
            myfile <<";";
           myfile <<"DegreeP";
           myfile <<";";
           //cout<<"VertexL"<<VL<<endl;
           myfile <<"IntersectionDE";
           myfile <<";";
           myfile <<"unionDE";
           myfile <<";";
            myfile <<"label";
            myfile <<";";
            myfile <<"Dlabel";
            myfile <<";";
            myfile <<"Elabel";
            myfile <<";";
           myfile <<"All";
           int expsizestart=5;int expsizeend=10;int expsize=expsizeend-expsizestart;
    for (int kk=expsizestart;kk<expsizeend;kk++){
    Graph* query_graph = new Graph(true);
    input_query_graph_file=input_qpartA+(qr[xx])+to_string(kk)+input_end;
    //input_query_graph_file="test_case_1.graph";
    //input_query_graph_file="query_graph\\query_dense_16_1.graph";
    query_graph->loadGraphFromFile(input_query_graph_file);
    int sizq=query_graph->getVerticesCount();
    MatrixXd eigenVq1(sizq,10);
    MTcalc12(query_graph, query_graph->getGraphMaxDegree(),eigenVq1,true,30);
    double x;
    cout<<"before";
    
    PrunRes(query_graph,data_graph,eigenVD1,eigenVq1,sizd, sizq, testeigen,avgPruned , avgPrunedC,
     avgDprun , avgcheck, avgcheck1 , avgVL , avgmagkas ,avgDL,avgEL );
     //CSInit(data_graph,query_graph,eigenVD1,eigenVq1);
    cout<<"after";
        

    }
        for (int avg=1;avg<10;avg++){
            cout<<" "<<endl;
            cout<<"Number of Eigenvalues "<< avg<<endl;
            cout<<"pruned   "<<avgPruned[avg]/expsize<<endl;
            cout<<"prunedC  "<<avgPrunedC[avg]/expsize<<endl;
           cout<<"prunedD   "<<avgDprun[avg]/expsize<<endl;
           //cout<<"VertexL"<<VL<<endl;
           cout<<"checkIntersection "<<avgcheck[avg]/expsize<<endl;
           cout<<"checkUnion    "<<avgcheck1[avg]/expsize<<endl;
            cout<<"VL "<<avgVL[avg]/expsize<<endl;
           cout<<"all    "<<avgmagkas[avg]/expsize<<endl;
        }
        for (int avg=1;avg<10;avg++){
        
         myfile <<" \n";
         //myfile <<avg;
         //myfile <<" : Number of Eigenvalues used is \n";
            myfile <<avgPruned[avg]/expsize;
            myfile <<";";
            myfile <<avgPrunedC[avg]/expsize;
            myfile <<";";
           myfile <<avgDprun[avg]/expsize;
           myfile <<";";
           //cout<<"VertexL"<<VL<<endl;
           myfile <<avgcheck[avg]/expsize;
           myfile <<";";
           myfile <<avgcheck1[avg]/expsize;
           myfile <<";";
            myfile <<avgVL[avg]/expsize;
            myfile <<";";
            myfile <<avgDL[avg]/expsize;
            myfile <<";";
            myfile <<avgEL[avg]/expsize;
            myfile <<";";
           myfile <<avgmagkas[avg]/expsize;
           
        }

  
  myfile.close();
    }
}