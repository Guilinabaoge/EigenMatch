#include "GrM.cpp"
#include "eigenHelper.cpp"
#include "cs.h"
#include <chrono>
#include <thread>
#include <Eigen/StdVector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <utility>
#include <cmath>
#include <Eigen/SparseCore>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include "IO.cpp"
// using namespace cs;
using namespace std;
using namespace Eigen;
// using namespace Spectra;
using namespace std::chrono;
int CSSizeReal(vector<vector<CSV>> &FCS, int qsiz);
bool ReverseRefinementNOTES(vector<vector<pair<ui, int>>> &QueryNlabel, vector<vector<CSV>> &FCS, int qsiz, Graph *query_graph);
int GQL(Graph *data_graph, Graph *query_graph);

bool ReverseRefinementOCS1(vector<vector<pair<ui, int>>> &QueryNlabel,
                           vector<vector<pair<ui, int>>> &QueryNlabel2, vector<int> &Qnew, vector<vector<CSV>> &FCS,
                           int qsiz, Graph *query_graph, float **&eigenVq1eigenVq1);
int CSInit(Graph *data_graph, Graph *query_graph, float **&eigenVD1, float **&eigenVq1);
void EdgesCheck(vector<vector<CSV>> &FCS, int qsiz);
inline VertexID findIndBS(vector<vector<CSV>> &FCS, VertexID IDC, VertexID IDQ);
inline bool checkB(int countD);
inline bool checkA(vector<pair<ui, int>> &EvalNeigb);
inline void removeVertexAndEgjesFK(vector<vector<CSV>> &FCS, CSV cvertex, int i, int deli);
void ExtractNI(vector<vector<pair<ui, int>>> &QueryNlabel, Graph *query_graph, int qsiz);

inline bool OneHopEigenSFPNS(CSV &cvertex, unordered_map<int, int> &ID, unordered_set<int> &IDD, int *IDDLC,
                             vector<T> &tripletList, vector<pair<ui, int>> &EvalNeigb, Graph *query_graph, int countD, vector<pair<VertexID, VertexID>> &q_curr)
{

    ID.insert({cvertex.ID, 0});
    IDDLC[0] = 1;
    IDDLC[1] = 1;
    IDDLC[2] = 1;
    VertexID vx1 = 0;
    VertexID vx2 = 0;
    int count = 1;
    for (int k = 0; k < cvertex.edges.size(); k++)
    {
        if (cvertex.edges[k].first == 1000000)
            continue;

        auto result = ID.insert({cvertex.edges[k].second, count});

        if (result.second)
        {
            tripletList.emplace_back(T(0, count, -1));
            tripletList.emplace_back(T(count, 0, -1));
            count++;
            for (int t = 0; t < EvalNeigb.size(); t++)
            {
                if (EvalNeigb[t].first == query_graph->getVertexLabel(cvertex.edges[k].first))
                {
                    EvalNeigb[t].second--;
                    break;
                }
            }
        }
        auto result1 = IDD.insert(cvertex.edges[k].first);
        if (result1.second)
        {
            countD--;
            IDDLC[1]++;
        }

        q_curr.emplace_back(cvertex.edges[k]);
    }
    IDDLC[0] = count;

    return (checkA(EvalNeigb) && checkB(countD));
}

void combineEdgesSet(int kk, CSV cvertexpair, unordered_set<ui> edgesc, vector<vector<CSV>> &FCS, vector<pair<VertexID, VertexID>> &q_curr)
{
    CSV cvertexpair2;
    pair<VertexID, VertexID> temp1;
    VertexID QN = 0;
    VertexID tempxx = 0;
    int counter = kk;
    vector<pair<LabelID, ui>> tempadd;
    bool stop = true;
    for (int i = 0; i < cvertexpair.edges.size(); i++)
    {
        if (cvertexpair.edges[i].first == 1000000)
            continue;
        if (edgesc.count(cvertexpair.edges[i].second) == 0)
            edgesc.insert(cvertexpair.edges[i].second);
    }

    while (counter < q_curr.size() - 1)
    {
        counter++;
        if (cvertexpair.ID == (q_curr[counter].second))
        {
            temp1 = q_curr[counter];
            QN = temp1.first;
            if (QN == 1000000)
                continue;

            tempxx = findIndBS(FCS, temp1.second, QN);
            bool add = true;
            for (int i = 0; i < FCS[QN][tempxx].edges.size(); i++)
            {
                if (FCS[QN][tempxx].edges[i].first == 1000000)
                    continue;
                if (edgesc.count(FCS[QN][tempxx].edges[i].second) == 0)
                    edgesc.insert(FCS[QN][tempxx].edges[i].second);
            }
            q_curr[counter].first = 1000000;
        }
    }
}

void combineEdges(int kk, CSV &cvertexpair, vector<vector<CSV>> &FCS,
                  vector<pair<VertexID, VertexID>> &q_curr, int qID)
{
    CSV cvertexpair2;
    pair<VertexID, VertexID> temp1;
    VertexID QN = 0;
    VertexID tempxx = 0;
    int counter = kk - 1;
    vector<pair<LabelID, ui>> tempadd;
    ui size1 = 0;
    ui size2 = 0;
    while (counter < q_curr.size() - 1)
    {
        counter++;
        if (cvertexpair.ID == (q_curr[counter].second))
        {

            temp1 = q_curr[counter];
            QN = temp1.first;
            if (QN == 1000000 || QN == qID)
                continue;
            tempxx = findIndBS(FCS, temp1.second, QN);
            size1 = FCS[QN][tempxx].edges.size();
            bool add = true;
            size2 = cvertexpair.edges.size();
            for (int i = 0; i < size1; i++)
            {
                if (FCS[QN][tempxx].edges[i].first == 1000000 || FCS[QN][tempxx].edges[i].first == qID)
                    // if (FCS[QN][tempxx].edges[i].first == 1000000 )
                    continue;
                add = true;
                for (int j = 0; j < size2; j++)
                {
                    if (cvertexpair.edges[j].first == 1000000 || cvertexpair.edges[j].first == qID)
                        continue;
                    if (FCS[QN][tempxx].edges[i].second == cvertexpair.edges[j].second)
                    {
                        add = false;
                        break;
                    }
                }
                if (add)
                    tempadd.push_back(make_pair(FCS[QN][tempxx].edges[i].first, FCS[QN][tempxx].edges[i].second));
            }
            q_curr[counter].first = 1000000;
            for (int a = 0; a < tempadd.size(); a++)
            {
                cvertexpair.edges.push_back(tempadd[a]);
            }
            tempadd.clear();
        }
    }
}

inline bool SecHopEigenPNSSJE(vector<pair<VertexID, VertexID>> &q_curr, unordered_map<int, int> &ID,
                              int *IDDLC, vector<vector<CSV>> &FCS, vector<T> &tripletList, Graph *query_graph, ui Omax, int qID)
{
    pair<VertexID, VertexID> temp1;
    VertexID tempxx = 0;
    CSV cvertexpair;
    vector<VertexID> temp2;
    VertexID vx1 = 0;
    VertexID vx2 = 0;
    int count = IDDLC[0];
    VertexID QN = 0;
    VertexID DN = 0;
    unordered_set<ui> EdgeF;
    EdgeF.reserve((query_graph->getGraphMaxDegree() + 2) * 50);
    int kk = 0;
    while (kk < q_curr.size())
    {
        temp1 = q_curr[kk];
        kk++;

        QN = temp1.first;
        if (QN == 1000000)
            continue;

        if (count > Omax)
        {
            IDDLC[0] = count;
            tripletList.clear();
            return true;
        }
        tempxx = findIndBS(FCS, temp1.second, QN);
        cvertexpair.ID = FCS[QN][tempxx].ID;
        cvertexpair.edges.clear();
        cvertexpair.edges = FCS[QN][tempxx].edges;
        EdgeF.clear();
        auto it = ID.find(cvertexpair.ID);
        vx1 = it->second;
        combineEdges(kk, cvertexpair, FCS, q_curr, qID);

        // combineEdgesSet(kk,cvertexpair,EdgeF,FCS,q_curr);
        for (int i = 0; i < cvertexpair.edges.size(); i++)
        {
            // if (cvertexpair.edges[i].first == 1000000)
            if (cvertexpair.edges[i].first == 1000000 || cvertexpair.edges[i].first == qID)
                continue;

            VertexID vtemp = cvertexpair.edges[i].second;
            auto [it1, success] = ID.try_emplace(vtemp, count);

            if (success)
            {
                vx2 = count;
                count++;
                IDDLC[0] = count;

                if (count > Omax)
                    break;
            }
            else
            {
                vx2 = ID[vtemp];
            }
            tripletList.emplace_back(T(vx2, vx1, -1));
            tripletList.emplace_back(T(vx1, vx2, -1));
        }
        if (count > Omax)
            break;
    }
    q_curr.clear();
    IDDLC[0] = count;
    return (true);
}

inline bool OneHopEigenNOTESFPNS(unordered_set<int> &IDN, unordered_set<int> &IDD, CSV &cvertex, vector<pair<ui, int>> &EvalNeigb, Graph *query_graph,
                                 int countD)
{

    // unordered_set<int> IDN;
    // unordered_set<int> IDD;

    IDN.clear();
    IDD.clear();
    IDN.insert(cvertex.ID);
    int count = 1;
    for (int k = 0; k < cvertex.edges.size(); k++)
    {
        if (cvertex.edges[k].first == 1000000)
            continue;

        // if (IDN.insert(cvertex.edges[k].second).second)
        //{
        if (IDN.count(cvertex.edges[k].second) == 0)
        {
            IDN.insert(cvertex.edges[k].second);
            count++;
            for (int t = 0; t < EvalNeigb.size(); t++)
            {
                if (EvalNeigb[t].first == query_graph->getVertexLabel(cvertex.edges[k].first))
                {
                    EvalNeigb[t].second--;
                }
            }
        }
        // if (IDD.insert(cvertex.edges[k].first).second)
        //{
        if (IDD.count(cvertex.edges[k].first) == 0)
        {
            IDD.insert(cvertex.edges[k].first);
            countD--;
        }
    }

    IDD.clear();
    IDN.clear();

    return (checkA(EvalNeigb) && checkB(countD));
}

void ExtractNI(vector<vector<pair<ui, int>>> &QueryNlabel, Graph *query_graph, int qsiz)
{
    const VertexID *u_nbrs;
    bool ad = true;
    int countkk = 0;
    ui u_nbrs_count;
    for (int i = 0; i < qsiz; i++)
    {
        u_nbrs = query_graph->getVertexNeighbors(i, u_nbrs_count);
        vector<pair<ui, int>> QueryVec;

        for (int j = 0; j < u_nbrs_count; j++)
        {
            ui labela = query_graph->getVertexLabel(u_nbrs[j]);

            ad = true;

            for (int k = 0; k < countkk; k++)
            // for (int k = 0; k < QueryVec.size(); k++)
            {

                if (QueryVec[k].first == labela)
                {
                    QueryVec[k].second++;
                    ad = false;
                }
            }
            if (ad)
            {
                QueryVec.emplace_back(make_pair(labela, 1));
                countkk++;
            }
        }
        QueryNlabel.emplace_back(QueryVec);
        QueryVec.clear();
        countkk = 0;
    }
}

bool InitPrunTCSR(vector<vector<CSV>> &FCS, int qsiz, Graph *query_graph)
{
    int jj = 0;
    ui VDP;
    VertexID rev;
    bool ret = false;
    for (VertexID kk = 0; kk < qsiz; kk++)
    {
        jj = FCS[kk].size();
        VDP = query_graph->getVertexDegree(kk);
        while (jj > 0)
        {
            jj--;
            if (FCS[kk][jj].Ichange == true)
            {
                if (FCS[kk][jj].edges.size() == 0)
                {
                    ret = true;
                    FCS[kk].erase(FCS[kk].begin() + jj);
                }

                else if (FCS[kk][jj].edges.size() < VDP)
                {
                    VertexID i = 0;

                    while (i < FCS[kk][jj].edges.size())
                    {
                        rev = findIndBS(FCS, FCS[kk][jj].edges[i].second, FCS[kk][jj].edges[i].first); // vertex to remove ID?
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
                    ret = true;
                }

                FCS[kk][jj].Ichange = false;
            }
        }
    }
    return ret;
}

void EdgesCSBasic(vector<vector<CSV>> &FCS, int qsiz, int dsiz, Graph *data_graph, Graph *query_graph)
{
    ui u_nbrs_count = 0;
    ui u_nbrs_countD = 0;
    int sizA = 0;
    int sizC;
    VertexID VID = 0;
    VertexID de = 0;
    VertexID cne = 0;
    VertexID labela = 0;
    for (VertexID a = 0; a < qsiz; a++)
    { // for every node of the query
        const VertexID *u_nbrs = query_graph->getVertexNeighbors(a, u_nbrs_count);
        // VertexID* u_nbrs = query_graph->getVertexNeighbors(a, u_nbrs_count);
        sizA = FCS[a].size();
        for (VertexID c = 0; c < u_nbrs_count; c++)
        { // can do binary search to find if an elemnt i in FCS
            if (u_nbrs[c] < a)
                continue;
            cne = u_nbrs[c];
            labela = query_graph->getVertexLabel(cne);
            sizC = FCS[cne].size();
            for (VertexID b = 0; b < sizA; b++)
            { // for every candidate vertex of every node of the query
                // candidate vertex for a query a is FCS[a][b].ID
                VID = FCS[a][b].ID;
                // const VertexID* u_nbrsD = data_graph->getVertexNeighbors(FCS[a][b].ID, u_nbrs_countD); //real neigboors of the candidate vertex
                const VertexID *u_nbrsD = data_graph->getNeighborsByLabel(FCS[a][b].ID, labela, u_nbrs_countD); // real neigboors of the candidate vertex
                                                                                                                // neigboors of query A node.
                                                                                                                // for every neigboor of the query
                // u_nbrs[c] gives the vertex ID of the queris current neigborhood

                for (VertexID d = 0; d < sizC; d++)
                { // for every candidate of the query vertex
                    // de=FCS[cne][d].ID;
                    for (VertexID e = 0; e < u_nbrs_countD; e++)
                    { // for every neigboor of the candidate vertex of the real graph
                        if (u_nbrsD[e] == FCS[cne][d].ID)
                        {
                            // FCS[a][b].edges.emplace_back(make_pair(cne,de));
                            // FCS[cne][d].edges.emplace_back(make_pair(a,VID));
                            FCS[a][b].edges.emplace_back(make_pair(cne, FCS[cne][d].ID));
                            FCS[cne][d].edges.emplace_back(make_pair(a, VID));
                        }
                        else if (u_nbrsD[e] > FCS[cne][d].ID)
                        { // optimization
                            //;
                            break;
                        }
                    }
                }
            }
        }
    }
}

void EdgesCSBasicSet(vector<vector<CSV>> &FCS, int qsiz, int dsiz, Graph *data_graph, Graph *query_graph)
{
    ui u_nbrs_count = 0;
    ui u_nbrs_countD = 0;
    int sizA = 0;
    int sizC;
    VertexID VID = 0;
    VertexID de = 0;
    VertexID cne = 0;
    VertexID labela = 0;
    unordered_map<ui, ui> s[qsiz];
    for (int i = 0; i < qsiz; i++)
    {
        s[i].reserve(FCS[i].size());
        for (int j = 0; j < FCS[i].size(); j++)
            s[i].insert({FCS[i][j].ID, j});
    }

    for (VertexID a = 0; a < qsiz; a++)
    { // for every node of the query
        const VertexID *u_nbrs = query_graph->getVertexNeighbors(a, u_nbrs_count);
        // VertexID* u_nbrs = query_graph->getVertexNeighbors(a, u_nbrs_count);
        sizA = FCS[a].size();
        for (VertexID c = 0; c < u_nbrs_count; c++)
        {
            if (u_nbrs[c] < a)
                continue;
            cne = u_nbrs[c];
            labela = query_graph->getVertexLabel(cne);
            sizC = FCS[cne].size();
            for (VertexID b = 0; b < sizA; b++)
            { // for every candidate vertex of every node of the query
                // candidate vertex for a query a is FCS[a][b].ID
                VID = FCS[a][b].ID;
                // const VertexID* u_nbrsD = data_graph->getVertexNeighbors(FCS[a][b].ID, u_nbrs_countD); //real neigboors of the candidate vertex
                const VertexID *u_nbrsD = data_graph->getNeighborsByLabel(FCS[a][b].ID, labela, u_nbrs_countD); // real neigboors of the candidate vertex
                for (VertexID e = 0; e < u_nbrs_countD; e++)
                { // for every neigboor of the candidate vertex of the real graph
                    // if(u_nbrsD[e]== FCS[cne][d].ID){
                    auto got = s[cne].find(u_nbrsD[e]);
                    if (got != s[cne].end())
                    {
                        FCS[a][b].edges.emplace_back(make_pair(cne, FCS[cne][got->second].ID));
                        FCS[cne][got->second].edges.emplace_back(make_pair(a, VID));
                    }
                }
            }
        }
    }
    for (int i = 0; i < qsiz; i++)
    {
        s[i].clear();
    }
}

void VerticesCSN(vector<vector<CSV>> &FCS, int qsiz, int dsiz, Graph *data_graph, Graph *query_graph, float **&eigenVD1, float **&eigenVq1, vector<vector<pair<ui, int>>> &QueryNlabel)
{

    VectorXd devalues;
    VectorXd qevalues;
    ui com = data_graph->getGraphMaxLabelFrequency();
    ui copies = query_graph->getLabelsCount();
    ui labelsNum[copies];
    ui reverseLab[310];
    const ui *labelData[copies];
    LabelID label = 0;
    for (int i = 0; i < 310; i++)
        reverseLab[i] = 310;
    int pos = 0;
    for (int i = 0; i < qsiz; i++)
    {
        label = query_graph->getVertexLabel(i);
        ui vdata_vertex_num = 0;
        if (reverseLab[label] == 310)
        {
            reverseLab[label] = pos;
            labelData[pos] = data_graph->getVerticesByLabel(label, vdata_vertex_num);
            labelsNum[pos] = vdata_vertex_num;
            pos++;
        }
    }

    vector<CSV> CS;
    int k = 0;
    ui data_vertex_num;
    int prunES = qsiz - 1;
    if (prunES > 16)
        prunES = 20;
    if (prunES < 10)
        prunES = 10;
    prunES = 10;
    for (int i = 0; i < qsiz; i++)
    {
        label = query_graph->getVertexLabel(i);
        ui degree = query_graph->getVertexDegree(i);
        data_vertex_num = 0;
        ui reserveS = com * degree;
        for (ui j = 0; j < labelsNum[reverseLab[label]]; ++j)
        {
            VertexID data_vertex = labelData[reverseLab[label]][j];
            //
            if (data_graph->getVertexDegree(data_vertex) >= degree)
            {

                bool con = true;
                for (int kk = 0; kk < prunES; kk++)
                {
                    if (eigenVD1[data_vertex][kk] < eigenVq1[i][kk])
                        if (eigenVq1[i][kk] - eigenVD1[data_vertex][kk] > 0.1)
                        {
                            con = false;
                            break;
                        }
                } // con=true;
                if (con)
                {
                    ui u_nbrs_countD;

                    for (k = 0; k < QueryNlabel[i].size() - 1; k++)
                    {
                        const VertexID *u_nbrsD = data_graph->getNeighborsByLabel(data_vertex, QueryNlabel[i][k].first, u_nbrs_countD);

                        if (u_nbrs_countD < QueryNlabel[i][k].second)

                            break;
                    } // k = QueryNlabel[i].size() - 1;
                    if (k == QueryNlabel[i].size() - 1)
                    {
                        const VertexID *u_nbrsD = data_graph->getNeighborsByLabel(data_vertex, QueryNlabel[i][k].first, u_nbrs_countD);
                        if (u_nbrs_countD >= QueryNlabel[i][k].second)

                        {

                            CSV cat(10, data_vertex, reserveS);
                            CS.emplace_back(cat);
                        }
                    }
                }
            }
        }
        FCS.emplace_back(CS);
        CS.clear();
    }
}

void EdgesCheck(vector<vector<CSV>> &FCS, int qsiz)
{
    CSV cat;
    CSV dog;
    bool KAt = false;
    VertexID vertexIDDD;
    for (VertexID a = 0; a < qsiz; a++)
    { // for every node of the query
        for (VertexID b = 0; b < FCS[a].size(); b++)
        {
            cat = FCS[a][b];
            for (int i = 0; i < cat.edges.size(); i++)
            {
                KAt = false;
                VertexID qID = cat.edges[i].first;
                VertexID vID = cat.edges[i].second;
                VertexID cvID = findIndBS(FCS, cat.edges[i].second, cat.edges[i].first);
                dog = FCS[qID][cvID];
                for (int j = 0; j < dog.edges.size(); j++)
                {
                    if (dog.edges[j].first == a && dog.edges[j].second == cat.ID)
                        KAt = true;
                }
                if (KAt == false)
                {
                    cout << "Wrong" << endl;
                }
            }
        }
    }
}
inline bool checkA(vector<pair<ui, int>> &EvalNeigb)
{
    for (int i = 0; i < EvalNeigb.size(); i++)
    {
        if (EvalNeigb[i].second > 0)
            return false;
    }
    return true;
}
inline bool checkB(int countD)
{
    if (countD > 0)
        return false;
    return true;
}

inline void removeVertexAndEgjesFK(vector<vector<CSV>> &FCS, CSV cvertex, int i, int deli)
{
    VertexID vx1;
    CSV cvertexpair;
    for (int j = 0; j < cvertex.edges.size(); j++)
    {
        // BSCHange
        if (cvertex.edges[j].first == 1000000)
            continue;
        vx1 = findIndBS(FCS, cvertex.edges[j].second, cvertex.edges[j].first);
        // vx1=findInd(FCS,cvertex.edges[j].second,cvertex.edges[j].first);
        cvertexpair = FCS[cvertex.edges[j].first][vx1];
        cvertexpair.change = true;
        FCS[cvertex.edges[j].first][vx1].change = true;
        for (int k = 0; k < cvertexpair.edges.size(); k++)
        {
            if (cvertexpair.edges[k].first == 1000000)
                continue;
            if (cvertexpair.edges[k].first == i && cvertexpair.edges[k].second == cvertex.ID)
            {
                FCS[cvertex.edges[j].first][vx1].edges[k].first = 1000000;
                break;
            }
        }
    }
    // while (!FCS[i][deli].edges.empty())
    //   FCS[i][deli].edges.pop_back();
    FCS[i][deli].edges.clear();
    FCS[i][deli].deleted = true;
    // FCS.erase( remove (FCS[i].begin(), FCS[i].end(), 1000000), FCS[i].end() );

    // FCS[i].erase(FCS[i].begin()+deli);
}

inline void removeVertexAndEgjesFKNP(vector<vector<CSV>> &FCS, CSV cvertex, int i, int deli)
{
    VertexID vx1;
    CSV cvertexpair;
    for (int j = 0; j < cvertex.edges.size(); j++)
    {
        // BSCHange
        if (cvertex.edges[j].first == 1000000)
            continue;
        vx1 = findIndBS(FCS, cvertex.edges[j].second, cvertex.edges[j].first);
        // vx1=findInd(FCS,cvertex.edges[j].second,cvertex.edges[j].first);
        cvertexpair = FCS[cvertex.edges[j].first][vx1];
        cvertexpair.IPchange = true;
        FCS[cvertex.edges[j].first][vx1].IPchange = true;
        for (int k = 0; k < cvertexpair.edges.size(); k++)
        {
            if (cvertexpair.edges[k].first == 1000000)
                continue;
            if (cvertexpair.edges[k].first == i && cvertexpair.edges[k].second == cvertex.ID)
            {
                FCS[cvertex.edges[j].first][vx1].edges[k].first = 1000000;
                break;
            }
        }
    }
    // while (!FCS[i][deli].edges.empty())
    //   FCS[i][deli].edges.pop_back();
    FCS[i][deli].edges.clear();
    FCS[i][deli].deleted = true;
    // FCS.erase( remove (FCS[i].begin(), FCS[i].end(), 1000000), FCS[i].end() );

    // FCS[i].erase(FCS[i].begin()+deli);
}

void SpectralMatchingIter(int sizd, MatrixXd eigenVD1, Graph *data_graph, string input_query_graph_file, int expsizestart, int expsizeend)
{
    string input_end = ".graph";

    for (int kk = expsizestart; kk < expsizeend; kk++)
    {
        Graph *query_graph = new Graph(true);
        query_graph->loadGraphFromFile(input_query_graph_file + to_string(kk) + input_end);
        int sizq = query_graph->getVerticesCount();
        MatrixXd eigenVq1(sizq, 10);
        MTcalc12(query_graph, query_graph->getGraphMaxDegree(), eigenVq1, true, 30);
        // CSInit(data_graph,query_graph,eigenVD1,eigenVq1);
    }
}
int SpectralMatching(int sizd, float **&eigenVD1, Graph *data_graph, string input_query_graph_file)
{
    // auto start2 = high_resolution_clock::now();

    Graph *query_graph = new Graph(true);
    query_graph->loadGraphFromFile(input_query_graph_file);
    int sizq = query_graph->getVerticesCount();
    // int sizq=query_graph->getVerticesCount();
    int Eprun = sizq - 1;
    if (Eprun > 16)
        Eprun = 20;
    else if (Eprun < 10)
        Eprun = 10;
    Eprun = 10;
    MatrixXd eigenVq1(sizq, Eprun);

    MTcalc12(query_graph, query_graph->getGraphMaxDegree(), eigenVq1, true, Eprun);
    float **eigenQ = NULL;
    eigenQ = new float *[sizq];

    for (ui i = 0; i < sizq; ++i)
    {
        eigenQ[i] = new float[Eprun];
        for (ui j = 0; j < Eprun; j++)
        {
            eigenQ[i][j] = eigenVq1(i, j);
            // cout << eigenQ[i][j] << ",";
        }
        // cout << i << endl;
    }
    // auto stop2 = high_resolution_clock::now();
    // auto duration2 = duration_cast<milliseconds>(stop2 - start2);
    // cout<<"First Step"<<duration2.count()<<endl;
    return CSInit(data_graph, query_graph, eigenVD1, eigenQ);
}

inline void eigcomp1(VectorXd &qevalues, VectorXd &devalues, int comp)
{

    bool print1 = true;
    double lala = 0;
    for (int i = 0; i < comp; i++)
    {
        if (abs(devalues[i] - qevalues[i]) > 0.001)
            print1 = false;
        lala = lala + abs(devalues[i] - qevalues[i]);
    }
    if (print1)
        cout << "WoW" << endl;
}

inline VertexID findIndBS(vector<vector<CSV>> &FCS, VertexID IDC, VertexID IDQ)
{
    int lo = 0, hi = FCS[IDQ].size() - 1;
    int mid;
    // This below check covers all cases , so need to check
    // for mid=lo-(hi-lo)/2
    while (hi - lo > 1)
    {
        int mid = (hi + lo) / 2;
        if (FCS[IDQ][mid].ID < IDC)
        {
            lo = mid + 1;
        }
        else
        {
            hi = mid;
        }
    }
    if (FCS[IDQ][lo].ID == IDC)
    {
        return lo;
    }
    else if (FCS[IDQ][hi].ID == IDC)
    {
        return hi;
    }
    cout << "error Prob" << endl;
    return -10000;
}

int CSInit(Graph *data_graph, Graph *query_graph, float **&eigenVD1, float **&eigenVq1)
{
    int qsiz = query_graph->getVerticesCount();
    int dsiz = data_graph->getVerticesCount();
    vector<vector<CSV>> FCS;

    FCS.reserve(qsiz);

    vector<vector<pair<ui, int>>> QueryNlabel;
    vector<vector<pair<ui, int>>> QueryNlabel2;
    vector<int> Qnew(qsiz);

    ExtractNI(QueryNlabel, query_graph, qsiz);

    // auto VCSN1 = high_resolution_clock::now();
    // auto VCS1 = high_resolution_clock::now();

    VerticesCSN(FCS, qsiz, dsiz, data_graph, query_graph, eigenVD1, eigenVq1, QueryNlabel); // parallel-Qsiz?
    // auto VCS2 = high_resolution_clock::now();
    // auto VCSdur = duration_cast<milliseconds>(VCS2 - VCS1);

    // auto ES1 = high_resolution_clock::now();
    // EdgesCSBasic(FCS,qsiz,dsiz,data_graph,query_graph);
    EdgesCSBasicSet(FCS, qsiz, dsiz, data_graph, query_graph);
    // auto ES2 = high_resolution_clock::now();
    // auto ESdur = duration_cast<milliseconds>(ES2 - ES1);
    // EdgesCheck(FCS, qsiz);
    // cout<<"EdgesCS"<<ESdur.count() <<endl;
    //

    auto IP1 = high_resolution_clock::now();

    while (InitPrunTCSR(FCS, qsiz, query_graph))
        ;
    // EdgesCheck(FCS, qsiz);
    //  printCSSize(FCS,qsiz);

    // auto IP2 = high_resolution_clock::now();
    // auto IPdur = duration_cast<milliseconds>(IP2 - IP1);

    // auto RE1 = high_resolution_clock::now();
    int GDegree = query_graph->getGraphMaxDegree();
    int maxDK = 0;
    int stop = 0;
    while (ReverseRefinementNOTES(QueryNlabel,
                                  FCS, qsiz, query_graph))
        stop++;
    // CSSizeReal(FCS, qsiz);
    // while (ReverseRefinementOCS1(QueryNlabel, QueryNlabel2, Qnew, FCS, qsiz, query_graph, eigenVq1))
    //    ;
    ReverseRefinementOCS1(QueryNlabel, QueryNlabel2, Qnew, FCS, qsiz, query_graph, eigenVq1);
    return CSSizeReal(FCS, qsiz);
}

int CSSizeReal(vector<vector<CSV>> &FCS, int qsiz)
{
    int count = 0;
    int Tcount = 0;
    for (int kk = 0; kk < qsiz; kk++)
    {
        for (int zz = 0; zz < FCS[kk].size(); zz++)
            // if(FCS[kk][zz].ID!=1000000)
            if (FCS[kk][zz].deleted == false)
            {
                count++;
                // cout<<", "<<FCS[kk][zz].ID;
            } // cout<<" "<<endl;

        Tcount = count + Tcount;
        // cout<<"i= "<<kk<<" size "<<count<<endl;
        count = 0;
    }
    return Tcount;
}

bool ReverseRefinementNOTES(vector<vector<pair<ui, int>>> &QueryNlabel, vector<vector<CSV>> &FCS, int qsiz, Graph *query_graph)
{
    unordered_set<int> SID;
    unordered_set<int> SIDD;
    SID.reserve(500);
    SIDD.reserve(qsiz);
    int IDDLC[3] = {0, 0, 0}; // IFD=0,IDD=1 IDL=2;
    VertexID vx1 = 0;
    VertexID vx2 = 0;
    VertexID vertexpair = 0;
    bool returnhere = false;
    VertexID vertex = 0;
    CSV cvertex;
    vector<pair<VertexID, VertexID>> q_curr;

    vector<pair<ui, int>> EvalNeigb;
    VectorXd evalues;
    VectorXd qevalues;

    // hop 1
    VertexID vertexDegree = 0;
    VertexID vertexlabel = 0;

    bool testcor;
    bool continueEL = false;
    // for (int i=qsiz-1;i>=0;i--){
    for (int i = 0; i < qsiz; i++)
    {

        vertexDegree = query_graph->getVertexDegree(i);
        // vertexlabel = query_graph->getVertexLabel(i);
        for (int j = 0; j < FCS[i].size(); j++)
        {
            cvertex = FCS[i][j];
            if (cvertex.deleted == true || cvertex.IPchange == false)
                continue;
            EvalNeigb.clear();
            EvalNeigb = QueryNlabel[i];
            // bool continueE = OneHopEigenNOTESFPNS(cvertex, IDDLC, EvalNeigb,
            //   query_graph, vertexDegree, q_curr);
            if (!OneHopEigenNOTESFPNS(SID, SIDD, cvertex, EvalNeigb,
                                      query_graph, vertexDegree))
            {
                removeVertexAndEgjesFKNP(FCS, cvertex, i, j);

                returnhere = true;
            }

            cvertex.IPchange = false;
        }
    }

    return returnhere;
}

bool ReverseRefinementOCS1(vector<vector<pair<ui, int>>> &QueryNlabel,
                           vector<vector<pair<ui, int>>> &QueryNlabel2, vector<int> &Qnew, vector<vector<CSV>> &FCS,
                           int qsiz, Graph *query_graph, float **&eigenVq1)
{
    vector<T> tripletList;
    tripletList.reserve(50 * query_graph->getGraphMaxDegree() * 2);
    unordered_map<int, int> SID;
    unordered_set<int> SIDD;
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

    VectorXd evalues;
    VectorXd qevalues;
    VertexID tempxx = 0;
    int metrwiters = 0;
    VertexID vertexDegree = 0;
    VertexID vertexlabel = 0;
    int prunC = 0;
    bool testcor;
    bool continueEL = false;
    ui oMax = qsiz * 1.5;
    for (int i = 0; i < qsiz; i++)
    {

        vertexDegree = query_graph->getVertexDegree(i);
        vertexlabel = query_graph->getVertexLabel(i);

        for (int j = 0; j < FCS[i].size(); j++)
        {
            cvertex = FCS[i][j];
            if (cvertex.deleted == true || cvertex.change == false)
                continue;

            EvalNeigb.clear();
            EvalNeigb = QueryNlabel[i];

            tripletList.clear();
            IDDLC[0] = 0;
            IDDLC[1] = 1;
            q_curr.clear();
            SID.clear();
            SIDD.clear();
            SIDD.insert(i);

            bool continueE = OneHopEigenSFPNS(cvertex, SID, SIDD, IDDLC, tripletList, EvalNeigb,
                                              query_graph, vertexDegree - 1, q_curr);

            if (continueE)
            {
                continueEL = SecHopEigenPNSSJE(q_curr, SID, IDDLC, FCS, tripletList, query_graph, oMax, i);
                if (IDDLC[0] <= oMax)
                {
                    std::map<int, int> count_uniques;
                    std::set<std::pair<int, int>> seen;
                    std::vector<Triplet<double>> unique_triplets;
                    for (auto t : tripletList)
                    {
                        if (seen.count({t.row(), t.col()}) == 0)
                        {
                            unique_triplets.push_back(Triplet<double>(t.row(), t.col(), -1));
                            seen.insert({t.row(), t.col()});
                            count_uniques[t.row()]++;
                        }
                    }
                    for (auto it : count_uniques)
                    {
                        unique_triplets.push_back(Triplet<double>(it.first, it.first, it.second));
                    }
                    tripletList = unique_triplets;

                    SparseMatrix<double> M(IDDLC[0], IDDLC[0]);
                    M.setFromTriplets(tripletList.begin(), tripletList.end(), [](double a, double b)
                                      { return b; });
                    // checkM(M);

                    VectorXd evalues(10);

                    calcEigens1(M, 10, evalues, IDDLC[0]);

                    bool con = true;
                    for (int dd = 0; dd < 10; dd++)
                    {
                        if (eigenVq1[i][dd] <= 0)
                            break;
                        if (evalues[dd] < eigenVq1[i][dd])
                            if (eigenVq1[i][dd] - evalues[dd] > 0.01)
                            {
                                con = false;
                                // break;
                            }
                    }
                    if (!con)
                    {
                        prunC++;
                        removeVertexAndEgjesFK(FCS, cvertex, i, j);
                        returnhere = true;
                    }

                    //}
                    else
                    {
                        cvertex.change = false;
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
    }
    return returnhere;
}

int GQL(Graph *data_graph, Graph *query_graph)
{
    ui **candidates = NULL;
    ui *candidates_count = NULL;
    ui *tso_order = NULL;
    TreeNode *tso_tree = NULL;
    ui *cfl_order = NULL;
    TreeNode *cfl_tree = NULL;
    ui *dpiso_order = NULL;
    TreeNode *dpiso_tree = NULL;
    TreeNode *ceci_tree = NULL;
    ui *ceci_order = NULL;
    std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> TE_Candidates;
    std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> NTE_Candidates;
    GQLFilter(data_graph, query_graph, candidates, candidates_count);
    sortCandidates(candidates, candidates_count, query_graph->getVerticesCount());

    /*for (int i=0;i<query_graph->getVerticesCount();i++){
        for (ui j = 0; j < candidates_count[i]; ++j)
                cout<<candidates[i][j]<<",";
                cout<<" "<<endl;
}*/
    int totalCand = 0;
    for (int i = 0; i < query_graph->getVerticesCount(); i++)
    {

        // cout<<"C(i) "<<candidates_count[i]<<",";
        totalCand = candidates_count[i] + totalCand;
    }
    // cout<<"All Cand "<<totalCand<<",";
    Edges ***edge_matrix = NULL;

    edge_matrix = new Edges **[query_graph->getVerticesCount()];
    for (ui i = 0; i < query_graph->getVerticesCount(); ++i)
    {
        edge_matrix[i] = new Edges *[query_graph->getVerticesCount()];
    }

    buildTables(data_graph, query_graph, candidates, candidates_count, edge_matrix);

    ui *matching_order = NULL;
    ui *pivots = NULL;
    ui **weight_array = NULL;
    size_t embedding_count = 0;
    size_t order_num = 0;
    size_t call_count = 0;
    // sscanf(input_order_num.c_str(), "%zu", &order_num);
    size_t output_limit = numeric_limits<size_t>::max();
    output_limit = 1;
    return totalCand;
    generateGQLQueryPlan(data_graph, query_graph, candidates_count, matching_order, pivots);
    checkQueryPlanCorrectness(query_graph, matching_order, pivots);
    printSimplifiedQueryPlan(query_graph, matching_order);

    embedding_count = exploreGraphQLStyle(data_graph, query_graph, candidates, candidates_count,
                                          matching_order, output_limit, call_count);

    // cout<<embedding_count<<endl;
}

int main(int argc, char **argv)
{
    Graph *query_graph = new Graph(true);
    Graph *data_graph = new Graph(true);
    Graph *data_graph1 = new Graph(true);
    string input_data_graph_file = "dataset\\dblp\\data_graph\\dblp.graph";
    // input_data_graph_file="graphB.txt";
    string input_query_graph_file = "dataset\\wordnet\\query_graph\\query_dense_8_2.graph";

    data_graph->loadGraphFromFile(input_data_graph_file);
    data_graph1->loadGraphFromFileG(input_data_graph_file);
    ui u_nbrs_count;
    data_graph->printGraphMetaData();
    query_graph->printGraphMetaData();
    int sizd = data_graph->getVerticesCount();
    int sizq = query_graph->getVerticesCount();
    int eigenV = 30;
    MatrixXd eigenVD1(sizd, eigenV);
    MatrixXd eigenVDC(sizd, eigenV);
    MatrixXd eigenVq1(sizd, eigenV);

    std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> TE_Candidates;
    std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> NTE_Candidates;
    int degree = 4 * data_graph->getGraphMaxDegree();

    VectorXcd evalues;
    VectorXd e1values;

    auto start1 = high_resolution_clock::now();
    //MTcalc12(data_graph, degree, eigenVD1, true, eigenV); // need small improvements.memory wise like it is implemented later on

    auto stop1 = high_resolution_clock::now();
    auto duration1 = duration_cast<milliseconds>(stop1 - start1);

    // MTcalc12(query_graph, degree,eigenVq1,true);
    // calcres(sizd,eigenVD1,data_graph);
    //saveData("dblp.csv", eigenVD1);

    float **eigenData = NULL;
    eigenData = new float *[sizd];

    for (ui i = 0; i < sizd; ++i)
    {
        eigenData[i] = new float[30];
    }
    eigenVDC = (openData("dblp.csv"));
    openData1("dblp.csv", eigenData);
    for (int i = 0; i < sizd; i++)
        for (int j = 0; j < 30; j++)
            if (eigenVDC(i, j) != eigenData[i][j])
                if (abs(eigenVDC(i, j) - eigenData[i][j]) > 0.001)
                    cout << "ooo lal" << eigenVDC(i, j) << "jj" << eigenData[i][j] << endl;
    cout << "IndexConstruction" << duration1.count() << endl;
    // string inputgraph="dataset\\wordnet\\query_graph\\query_dense_24_1.graph";
    int GQLCand = 0;
    int SpecCand = 0;
    float GQLTime = 0;
    float SpecTime = 0;

    string input_query_graph_file1 = "dataset\\dblp\\query_graph\\query_dense_32_";

    int cc = 0;
    for (int i = 1; i <= 200; i++)
    {
        input_query_graph_file = input_query_graph_file1 + to_string(i) + ".graph";
        // input_query_graph_file="graphB.txt";
        cout << input_query_graph_file << endl;
        auto GQL1 = high_resolution_clock::now();
        Graph *query_graph = new Graph(true);
        query_graph->loadGraphFromFileG(input_query_graph_file);
        GQLCand = GQLCand + GQL(data_graph1, query_graph);
        auto GQL2 = high_resolution_clock::now();
        auto GQLdur = duration_cast<milliseconds>(GQL2 - GQL1);
        cout << "GQL" << GQLdur.count() << endl;
        GQLTime = GQLTime + GQLdur.count();
        auto start2 = high_resolution_clock::now();
        cc = SpectralMatching(sizd, eigenData, data_graph, input_query_graph_file); // make the loop outside
        SpecCand = SpecCand + cc;
        auto stop2 = high_resolution_clock::now();
        if (cc < query_graph->getVerticesCount())
            cout << "PRoblem " << cc << endl;
        auto duration2 = duration_cast<milliseconds>(stop2 - start2);
        cout << "Sprectral Matching" << duration2.count() << endl;
        SpecTime = SpecTime + duration2.count();
        // delete []query_graph;
    }
    cout << "Spectral Time: " << SpecTime / 200 << endl;
    cout << "Spectral Cand: " << SpecCand / 200 << endl;
    cout << "GQL Time: " << GQLTime / 200 << endl;
    cout << "GQLCand : " << GQLCand / 200 << endl;
    string qr[5] = {"4_", "8_", "16_", "24_", "32_"};
    // inputgraph="dataset\\yeast\\query_graph\\query_dense_";
    auto start3 = high_resolution_clock::now();
    // SpectralMatchingIter(sizd,eigenVD1,data_graph,inputgraph+qr[3],1,200);//make the loop outside
    auto stop3 = high_resolution_clock::now();
    auto duration3 = duration_cast<milliseconds>(stop3 - start3);
    cout << "CS_Construction200avg" << (duration3.count()) / 200 << endl;

    return 0;
}