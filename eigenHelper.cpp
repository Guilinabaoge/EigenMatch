#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <iostream>
#include <cmath>
#include <set>
#include <unordered_map>
#include <unordered_set>
// #include "GrM.h"
#include "eigenhelper.h"
#include "thread"
typedef Eigen::Triplet<double> T;
using namespace Eigen;
using namespace std;
using namespace Spectra;
void printSM(SparseMatrix<double> &M);
// bool checkM2(SparseMatrix<double> mat);
bool checkM2(SparseMatrix<double> mat, int rowO, int ColO);
bool checkM(SparseMatrix<double> mat)
{
    int count = 0;
    int row = 0;
    int col = 0;
    for (int k = 0; k < mat.outerSize(); ++k)
    {

        count = 0;
        for (SparseMatrix<double>::InnerIterator it(mat, k); it; ++it)
        {
            count = count + it.value();
            if (it.value() < -1)
                cout << "prob" << endl;
            row = it.row(); // row index
            col = it.col(); // col index (here it is equal to k)
            // it.index(); // inner index, here it is equal to it.row()
            if (checkM2(mat, row, col) == false)
            {
                cout << "Problem" << endl;
            }
        }
        if (count != 0)
        {
            "Problem";
            return true;
        }
    }
    return false;
}

bool checkM2(SparseMatrix<double> mat, int rowO, int ColO)
{
    int row = 0;
    int col = 0;
    for (int k = 0; k < mat.outerSize(); ++k)
    {
        // count=0;
        for (SparseMatrix<double>::InnerIterator it(mat, k); it; ++it)
        {
            // count=count+it.value();

            row = it.row(); // row index
            col = it.col(); // col index (here it is equal to k)
            // it.row();   // row index
            // it.col();   // col index (here it is equal to k)
            // it.index(); // inner index, here it is equal to it.row()
            if (row == ColO && col == rowO)
                return true;
        }
    }
    return false;
}
VertexID checkANX(vector<VertexID> ID, VertexID CID)
{
    for (int i = 0; i < ID.size(); i++)
        if (ID[i] == CID)
            return i;
    return 1000000;
}

VertexID checkA(VertexID *ID, VertexID vertex, int count)
{
    for (int i = 0; i <= count; i++)
    {
        if (ID[i] == vertex)
            return i;
    }
    return 1000000;
}

void printSM(SparseMatrix<double> &M)
{
    int count = 0;
    for (int k = 0; k < M.outerSize(); ++k)
    {
        for (SparseMatrix<double>::InnerIterator it(M, k); it; ++it)
        { // cout<<"query ID"<<i<<endl;
          // cout<<"test"<<endl;
            count = count + it.value();
            if (it.value() == -1)
            {

                // cout << "";
                // cout << " (" << it.row() << ",";
                // cout << it.col() << ")  ";
                // cout << "e " << it.col() << " " << it.row() << " "
                //     << "0" << endl;
                cout << "v" << it.row() << " " << it.col();
                // cout<<"v "<< it.col()<<" "<<"0";
                cout << " " << it.value() << endl;
            }
            // row index
            // col index (here it is equal to k)
            // it.index(); // inner index, here it is equal to it.row()
        } // cout<<" "<<endl;
        if (count != 0)
        {
            cout << "problem" << endl;
        }
        count = 0;
        // cout<<endl;
    }
    cout << " " << endl;
}

void ExtractAdL(SparseMatrix<double> &M, Graph *data_graph, int degree, int depth, VertexID vertex)
{
    VertexID *neighbors_;
    VertexID *ID; // add size then resize
    VertexID *IDL;
    ID = new VertexID[data_graph->getVerticesCount()];
    IDL = new VertexID[data_graph->getVerticesCount()];
    memset(ID, 0, sizeof(VertexID) * (data_graph->getVerticesCount()));
    memset(IDL, 0, sizeof(VertexID) * (data_graph->getVerticesCount()));
    VertexID vx1;
    VertexID vx2;
    VertexID vertexpair;
    int count = 0;

    ui u_nbrs_count;
    ui u_nbrs_count1;
    const VertexID *u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);

    ID[0] = vertex;
    queue<VertexID> q_curr;
    queue<VertexID> q_next;
    queue<VertexID> q_temp;

    M.insert(0, 0) = (float)u_nbrs_count;
    for (int j = 0; j < u_nbrs_count; ++j)
    {
        count++;
        M.insert(0, count) = -1.0;
        M.insert(count, 0) = -1.0;
        ID[count] = u_nbrs[j];
        q_curr.push(u_nbrs[j]);
    }

    for (int i = 1; i < depth; i++)
    {
        while (!q_curr.empty())
        {
            vertex = q_curr.front();
            q_curr.pop();
            u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);

            vx1 = checkA(ID, vertex, count);
            M.insert(vx1, vx1) = (float)u_nbrs_count;
            for (ui j = 0; j < u_nbrs_count; ++j)
            {
                vertexpair = u_nbrs[j];

                vx2 = checkA(ID, vertexpair, count);
                if (vx2 == 1000000)
                {
                    count++;
                    vx2 = count;
                    q_next.push(vertexpair);
                    ID[count] = vertexpair;
                    M.insert(vx1, vx2) = -1.0;
                    M.insert(vx2, vx1) = -1.0;
                    if (i == depth - 1)
                    {
                        IDL[vx2] = 1;
                    }
                }
                else
                {
                    M.coeffRef(vx1, vx2) = -1.0;
                    M.coeffRef(vx2, vx1) = -1.0;
                    if (i == depth - 1)
                    {
                        IDL[vx2]++;
                    }
                }
            }
        }
        if (!q_next.empty())
            q_curr.swap(q_next);
        else
        {
            i = depth;
        }
    }
    while (!q_curr.empty())
    {
        vertex = q_curr.front();
        q_curr.pop();
        vx1 = checkA(ID, vertex, count);

        M.insert(vx1, vx1) = IDL[vx1];
    }
}

void calcEigens1(SparseMatrix<double> M, int k, VectorXd &evalues, int count)
{
    int sizek = k * 2;
    int dek = k;

    if (count == 0 || count == 1)
    {
        evalues.resize(dek);

        evalues(0) = 1;
        for (int i = 1; i < dek; i++)
            evalues(i) = 0;

        return;
    }

    if (k >= count)
        k = count - 1;

    SparseGenMatProd<double> op(M);
    SymEigsSolver<SparseGenMatProd<double>> eigs(op, k, count);
    eigs.init();
    int nconv = eigs.compute(SortRule::LargestAlge);

    if (eigs.info() == CompInfo::Successful)
        evalues = eigs.eigenvalues();
    if (eigs.info() != CompInfo::Successful)
    {
        int info = static_cast<int>(eigs.info());

        evalues.resize(k);
        for (int ss = 0; ss < k; ss++)
            evalues(ss) = k;
    }

    if (evalues.size() < dek)
    {
        int sz = evalues.size();
        evalues.conservativeResize(dek);
        evalues(sz) = 0;
        sz++;
        for (int i = sz; i < dek; i++)
            evalues(i) = -1;
    }
}

void CompactADLEIGSet(Graph *data_graph, int degree, VectorXd &evalues, VertexID vertex, int depth, int Eprun)
{
    vector<T> tripletList;
    vector<VertexID> ID; // add size then resize
    vector<VertexID> IDD;

    VertexID *neighbors_;
    unordered_map<ui, ui> SID;
    VertexID *IDL;

    IDL = new VertexID[data_graph->getVerticesCount()];
   ;
    VertexID vx1;
    VertexID vx2 = 0;
    VertexID vertexpair;
    VertexID vertexprint;
    vertexprint = vertex;
    ui count = 1;

    ui u_nbrs_count;
    ui u_nbrs_count1;
    const VertexID *u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);
    int k = data_graph->getVerticesCount() - 1;
    if (k > 10)
        k = Eprun;
    if (k < 10)
        k = 10;
    if (u_nbrs_count > 150)
    {

        evalues.resize(k);
        for (int i = 0; i < k; i++)
            evalues(i) = 500;
        return;
    }
    SID.insert({vertex, 0});
    IDL[0] = 1;
    queue<VertexID> q_curr;
    queue<VertexID> q_next;
    queue<VertexID> q_temp;
    tripletList.push_back(T(0, 0, (float)u_nbrs_count));
    for (int j = 0; j < u_nbrs_count; ++j)
    {
        if (u_nbrs[j] == vertex)
            cout << "error :" << vertex << endl;

        tripletList.push_back(T(0, count, -1));
        tripletList.push_back(T(count, 0, -1));
        SID.insert({u_nbrs[j], count});


        IDL[count] = 1;
        count++;
        q_curr.push(u_nbrs[j]);
    }

    for (int i = 1; i < depth; i++)
    {
        while (!q_curr.empty())
        {
            vertex = q_curr.front();
            q_curr.pop();
            u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);
            if (u_nbrs_count > 150)
            {

                evalues.resize(k);
                for (int tt = 0; tt < k; tt++)
                    evalues(tt) = 500;
                
                return;
            }

            auto it = SID.find(vertex);
            vx1 = it->second;

            tripletList.push_back(T(vx1, vx1, (float)u_nbrs_count));


            for (ui j = 0; j < u_nbrs_count; ++j)
            {

                vertexpair = u_nbrs[j];
                if (vertexpair == vertex)
                    cout << "problem again!!" << endl;
                auto [it1, success] = SID.try_emplace(vertexpair, count);

                if (success)
                {
                    vx2 = count;
                    q_next.push(vertexpair);
                    IDL[count] = 1;
                    count++;
                }
                else
                {
                    vx2 = SID[vertexpair];
                    if (i == depth - 1)
                        IDL[vx2]++;
                }
                tripletList.push_back(T(vx1, vx2, -1));
                tripletList.push_back(T(vx2, vx1, -1));
            }
        }
        if (!q_next.empty())
            q_curr.swap(q_next);
        else
        {
            i = depth;
        }
    }
    while (!q_curr.empty())
    {
        vertex = q_curr.front();
        q_curr.pop();
        vx1 = SID[vertex];
        tripletList.push_back(T(vx1, vx1, IDL[vx1]));
    }

    SparseMatrix<double> M(count, count);
    M.setFromTriplets(tripletList.begin(), tripletList.end(), [](double a, double b)
                      { return b; });


    calcEigens1(M, k, evalues, count);
    tripletList.clear();
    ID.clear();
    SID.clear();

}
void CompactADLEIG(Graph *data_graph, int degree, VectorXd &evalues, VertexID vertex, int depth, int Eprun)
{
    vector<T> tripletList;
    vector<VertexID> ID; // add size then resize
    vector<VertexID> IDD;
    vector<VertexID> IDL;
    VertexID *neighbors_;
    VertexID vx1;
    VertexID vx2;
    VertexID vertexpair;
    VertexID vertexprint;
    vertexprint = vertex;
    int count = 0;

    ui u_nbrs_count;
    ui u_nbrs_count1;
    const VertexID *u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);
    int k = data_graph->getVerticesCount() - 1;
    if (k > 16)
        k = Eprun;
    if (k <= 15)
        k = 10;

    k = Eprun;
    if (u_nbrs_count > 100)
    {

        evalues.resize(k);
        for (int i = 0; i < k; i++)
            evalues(i) = 500;
        return;
    }
    ID.push_back(vertex);
    IDL.push_back(1);
    queue<VertexID> q_curr;
    queue<VertexID> q_next;
    queue<VertexID> q_temp;

    tripletList.push_back(T(0, 0, (float)u_nbrs_count));
    for (int j = 0; j < u_nbrs_count; ++j)
    {
        if (u_nbrs[j] == vertex)
            cout << "error :" << vertex << endl;
        count++;
        tripletList.push_back(T(0, count, -1));
        tripletList.push_back(T(count, 0, -1));
        ID.push_back(u_nbrs[j]);
        IDL.push_back(1);
        q_curr.push(u_nbrs[j]);
    }

    for (int i = 1; i < depth; i++)
    {
        while (!q_curr.empty())
        {
            vertex = q_curr.front();
            q_curr.pop();
            u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);
            if (u_nbrs_count > 100)
            {

                evalues.resize(k);
                for (int tt = 0; tt < k; tt++)
                    evalues(tt) = 500;
                return;
            }
            vx1 = checkANX(ID, vertex);
            tripletList.push_back(T(vx1, vx1, (float)u_nbrs_count));

            for (ui j = 0; j < u_nbrs_count; ++j)
            {

                vertexpair = u_nbrs[j];
                if (vertexpair == vertex)
                    cout << "problem again!!" << endl;

                vx2 = checkANX(ID, vertexpair);
                if (vx2 == 1000000)
                {
                    count++;
                    vx2 = count;
                    q_next.push(vertexpair);
                    ID.push_back(vertexpair);
                    IDL.push_back(1);
                    tripletList.push_back(T(vx1, vx2, -1));
                    tripletList.push_back(T(vx2, vx1, -1));
                    if (i == depth - 1)
                    {
                        IDL[vx2] = 1;
                    }
                }
                else
                {
                    tripletList.push_back(T(vx1, vx2, -1));
                    tripletList.push_back(T(vx2, vx1, -1));
                    if (i == depth - 1)
                    {
                        IDL[vx2]++;
                    }
                }
            }
        }
        if (!q_next.empty())
            q_curr.swap(q_next);
        else
        {
            i = depth;
        }
    }
    while (!q_curr.empty())
    {
        vertex = q_curr.front();
        q_curr.pop();
        vx1 = checkANX(ID, vertex);
        tripletList.push_back(T(vx1, vx1, IDL[vx1]));
    }
    count++; /*
                     map<int, int> count_uniques;
                     set<pair<int, int>> seen;
                     vector<Triplet<double>> unique_triplets;
                     for (auto t : tripletList)
                     {
                         if (seen.count({t.row(), t.col()}) == 0)
                         {
                             unique_triplets.push_back(Triplet<double>(t.row(), t.col(), t.value()));
                             seen.insert({t.row(), t.col()});
                             count_uniques[t.row()]++;
                         }
                     }
                     tripletList=unique_triplets;*/
    SparseMatrix<double> M(count, count);
    M.setFromTriplets(tripletList.begin(), tripletList.end(), [](double a, double b)
                      { return b; });
    // checkM(M);

    calcEigens1(M, k, evalues, count);
    tripletList.clear();
    ID.clear();
    IDL.clear();
}

void CompactADLEIGImpro(Graph *data_graph, int degree, VectorXd &evalues, VertexID vertex, int depth)
{
    vector<T> tripletList;
    vector<VertexID> ID; // add size then resize
    vector<VertexID> IDD;
    vector<VertexID> IDL;
    VertexID *neighbors_;
    unordered_map<ui, ui> SID;

    VertexID vx1;
    VertexID vx2;
    VertexID vertexpair;
    VertexID vertexprint;
    vertexprint = vertex;
    int count = 1;

    ui u_nbrs_count;
    ui u_nbrs_count1;
    const VertexID *u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);
    int k = 16;
    if (u_nbrs_count > 150)
    {
        evalues.resize(k);
        for (int i = 0; i < k; i++)
            evalues(i) = 500;
        return;
    }
    SID.insert({vertex, count});
    IDL.push_back(1);
    queue<VertexID> q_curr;
    queue<VertexID> q_next;
    queue<VertexID> q_temp;

    tripletList.push_back(T(0, 0, u_nbrs_count));
    for (int j = 0; j < u_nbrs_count; ++j)
    {
        if (u_nbrs[j] == vertex)
            continue;

        tripletList.push_back(T(0, count, -1));
        tripletList.push_back(T(count, 0, -1));
        SID.insert({u_nbrs[j], count});
        IDL.push_back(1);
        count++;
        q_curr.push(u_nbrs[j]);
    }

    for (int i = 1; i < depth; i++)
    {
        while (!q_curr.empty())
        {
            vertex = q_curr.front();
            q_curr.pop();
            u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);

            if (u_nbrs_count > 450 || count > 4500)
            {
                evalues.resize(k);
                for (int tt = 0; tt < k; tt++)
                    evalues(tt) = 500;
                return;
            }

            vx1 = SID[vertex];
            tripletList.push_back(T(vx1, vx1, u_nbrs_count));

            for (ui j = 0; j < u_nbrs_count; ++j)
            {
                vertexpair = u_nbrs[j];

                auto result = SID.insert({vertexpair, count});
                if (result.second)
                {

                    vx2 = count;
                    q_next.push(vertexpair);
                    IDL.push_back(1);
                    tripletList.push_back(T(vx1, vx2, -1));
                    tripletList.push_back(T(vx2, vx1, -1));
                    count++;
                    if (i == depth - 1)
                    {
                        IDL[vx2] = 1;
                    }
                }
                else
                {
                    vx2 = result.first->second;
                    if (vx2 == vx1)
                    {
                        continue;
                    }
                    tripletList.push_back(T(vx1, vx2, -1));
                    tripletList.push_back(T(vx2, vx1, -1));
                    if (i == depth - 1)
                    {
                        IDL[vx2]++;
                    }
                }
            }
        }

        if (!q_next.empty())
            q_curr.swap(q_next);
        else
        {
            i = depth;
        }
    }
    while (!q_curr.empty())
    {
        vertex = q_curr.front();
        q_curr.pop();
        vx1 = SID[vertex];
        tripletList.push_back(T(vx1, vx1, IDL[vx1]));
    }

    SparseMatrix<double> M(count, count);
    M.setFromTriplets(tripletList.begin(), tripletList.end(), [](double a, double b)
                      { return b; });

    calcEigens1(M, k, evalues, count);
    tripletList.clear();
    ID.clear();
    IDL.clear();
}

void MTEigCal(Graph *data_graph, int degree, MatrixXd &eigenVD, bool LE, int Eprun)
{

    auto AdJAdl1 = [](Graph *data_graph, int degree, VertexID svertex, VertexID evertex, MatrixXd &eigenVD, bool LE, int Eprun)
    {
        VectorXd evalues;
        int depth = 2;
        for (int i = svertex; i < evertex; i++)
        {

            if (LE)
                CompactADLEIG(data_graph, degree, evalues, i, depth, Eprun);

            else
                CompactADJEIG(data_graph, degree, evalues, i, depth);

            eigenVD.row(i) = evalues;
            evalues.setZero();
        }
    };

    int Tnum = 3;
    int siz = data_graph->getVerticesCount();

    int div = (int)(siz / Tnum);
    thread th[Tnum];
    for (int d = 0; d < Tnum - 1; d++)
    {
        th[d] = thread(AdJAdl1, data_graph, degree, div * d, div * (d + 1), ref(eigenVD), LE, Eprun);
    }
    th[Tnum - 1] = thread(AdJAdl1, data_graph, degree, div * (Tnum - 1), siz, ref(eigenVD), LE, Eprun);
    for (int d = 0; d < Tnum; d++)
        th[d].join();
}

void CompactADJEIG(Graph *data_graph, int degree, VectorXd &evalues, VertexID vertex, int depth)
{
    vector<T> tripletList;

    VertexID *neighbors_;
    VertexID *ID; // add size then resize
    VertexID *IDL;
    ID = new VertexID[data_graph->getVerticesCount()];
    IDL = new VertexID[data_graph->getVerticesCount()];
    memset(ID, 0, sizeof(VertexID) * (data_graph->getVerticesCount()));
    memset(IDL, 0, sizeof(VertexID) * (data_graph->getVerticesCount()));
    VertexID vx1;
    VertexID vx2;
    VertexID vertexpair;
    int count = 0;
    int k = 30;
    ui u_nbrs_count;
    ui u_nbrs_count1;
    const VertexID *u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);

    ID[0] = vertex;
    queue<VertexID> q_curr;
    queue<VertexID> q_next;
    queue<VertexID> q_temp;
    for (int j = 0; j < u_nbrs_count; ++j)
    {
        count++;
        tripletList.push_back(T(0, count, 1));
        tripletList.push_back(T(count, 0, 1));
        ID[count] = u_nbrs[j];
        q_curr.push(u_nbrs[j]);
    }

    for (int i = 1; i < depth; i++)
    {
        while (!q_curr.empty())
        {
            vertex = q_curr.front();
            q_curr.pop();
            u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);

            vx1 = checkA(ID, vertex, count);

            for (ui j = 0; j < u_nbrs_count; ++j)
            {
                vertexpair = u_nbrs[j];
                vx2 = checkA(ID, vertexpair, count);
                if (vx2 == 1000000)
                {
                    count++;
                    vx2 = count;
                    q_next.push(vertexpair);
                    ID[count] = vertexpair;
                    tripletList.push_back(T(vx1, vx2, 1));
                    tripletList.push_back(T(vx2, vx1, 1));
                }
                else
                {
                    tripletList.push_back(T(vx1, vx2, 1));
                    tripletList.push_back(T(vx2, vx1, 1));
                }
            }
        }

        if (!q_next.empty())
            q_curr.swap(q_next);
        else
        {
            i = depth;
        }
    }

    SparseMatrix<double> M(count, count);

    M.setFromTriplets(tripletList.begin(), tripletList.end(), [](double a, double b)
                      { return b; });
    if (checkM(M))
    {
        cout << "wrong Laplacian Matrix-Malakas" << endl;
    }

    calcEigens1(M, k, evalues, count);
    tripletList.clear();
    delete[] neighbors_;
    delete[] ID;
    delete[] IDL;
}
