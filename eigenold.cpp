
#include "eigenhelper.h"
#include "thread"
void MTcalc1(Graph *data_graph, int degree, MatrixXd &eigenVD)
{

    auto AdJAdl = [](Graph *data_graph, int degree, VertexID svertex, VertexID evertex, MatrixXd &eigenVD)
    {
        SparseMatrix<double> M(4000, 4000);
        SparseMatrix<double> Kt(5, 5);
        VectorXd evalues;
        int depth = 2;
        // cout<<"hi"<<endl;
        for (int i = svertex; i < evertex; i++)
        {
            // ExtractAdL(M,data_graph,degree,depth,i);
            ExtractAdj(M, data_graph, degree, depth, i);
            ///*

            if (i == 247)
            {
                // SparseMatrix<double> Kt(5,5);
                Kt.reserve(M.nonZeros());
                for (int k = 0; k < M.outerSize(); ++k)
                {
                    for (SparseMatrix<double>::InnerIterator it(M, k); it; ++it)
                    { // cout<<"query ID"<<i<<endl;
                        // cout<<"test"<<endl;
                        if (it.value() == 1)
                            Kt.insert(it.row(), it.col()) = it.value();
                        cout << " (" << it.row() << ",";
                        cout << it.col() << ")  ";
                        cout << it.value();
                        // row index
                        // col index (here it is equal to k)
                        // it.index(); // inner index, here it is equal to it.row()
                    }
                    // cout<<endl;
                }
                cout << endl;
                cout << endl;
                for (int tt = 0; tt < Kt.outerSize(); ++tt)
                {
                    for (SparseMatrix<double>::InnerIterator its(Kt, tt); its; ++its)
                    { // cout<<"query ID"<<i<<endl;
                        // cout<<"test"<<endl;

                        cout << " (" << its.row() << ",";
                        cout << its.col() << ")  ";
                        cout << its.value();
                        // row index
                        // col index (here it is equal to k)
                        // it.index(); // inner index, here it is equal to it.row()
                    }
                    // cout<<endl;
                }

            } //*/
            // evalues.setZero();
            // M.conservativeResize(M.nonZeros());
            // M.makeCompressed();
            calcEigens1(M, 10, evalues,5);

            if (i == 247)
            {
                // cout<<M.size()<<endl;
                cout << Kt.size() << endl;
                cout << Kt.nonZeros() << endl;

                cout << "compareEigen : " << endl;
                for (int jo = 0; jo < 10; jo++)
                    cout << evalues(jo) << ",";

                cout << " stop" << endl;
                calcEigens1(Kt, 10, evalues,5);
                cout << " stop" << endl;
                for (int jo = 0; jo < 10; jo++)
                    cout << evalues(jo) << ",";
            }

            // cout<<"hi4"<<endl;
            // std::cout << "Matrix:\n" << M.row(1) << std::endl;
            // calcEigens(M,4,eigenVD[i]);
            M.setZero();
            eigenVD.row(i) = evalues;
            evalues.setZero();
            // cout<<"hi3"<<endl;
            // std::cout << "Eigenvalues found:\n" <<  M << std::endl;
        }
    };
    // VectorXcd evalues;
    int Tnum = 1;
    int siz = data_graph->getVerticesCount();
    // siz=10;
    int div = (int)(siz / Tnum);
    thread th[Tnum];
    for (int d = 0; d < Tnum - 1; d++)
    {
        th[d] = thread(AdJAdl, data_graph, degree, div * d, div * (d + 1), ref(eigenVD));
    }
    th[Tnum - 1] = thread(AdJAdl, data_graph, degree, div * (Tnum - 1), siz, ref(eigenVD));
    // th[Tnum-1]=thread (AdJAdl,data_graph,degree,siz-2,siz,ref(eigenVD));
    for (int d = 0; d < Tnum; d++)
        th[d].join();
    // cout<<"hi"<<endl;
    // std::cout << "Eigenvalues found:\n" <<  eigenVD.row(0) << std::endl;
}

void ExtractAdj(SparseMatrix<double> &M, Graph *data_graph, int degree, int depth, VertexID vertex)
{
    VertexID *neighbors_;
    VertexID *ID; // add size then resize
    // ID = new VertexID[degree*depth];
    ID = new VertexID[data_graph->getVerticesCount()];
    // memset(ID, 0, sizeof(ui) * (degree*depth));
    memset(ID, 0, sizeof(VertexID) * (data_graph->getVerticesCount()));
    VertexID vx1;
    VertexID vx2;
    VertexID vertexpair;
    int count = 0;

    ui u_nbrs_count;
    const VertexID *u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);

    ID[0] = vertex;
    queue<VertexID> q_curr;
    queue<VertexID> q_next;
    queue<VertexID> q_temp;
    // cout<<u_nbrs_count<<" NC "<<endl;
    // cout<<vertex<<" vertex "<<endl;
    for (int j = 0; j < u_nbrs_count; ++j)
    {
        count++;
        M.insert(0, count) = 1.0;
        M.insert(count, 0) = 1.0;
        ID[count] = u_nbrs[j];
        q_curr.push(u_nbrs[j]);
        // cout<<u_nbrs[j]<<" edge pair "<<endl;
        if (u_nbrs[j] == vertex)
            cout << "hola" << endl;
    }

    for (int i = 1; i < depth; i++)
    {
        while (!q_curr.empty())
        {
            vertex = q_curr.front();
            q_curr.pop();
            u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);
            vx1 = checkA(ID, vertex, count);
            /// cout<<vertex<<" vertex "<<endl;
            // cout<<u_nbrs_count<<" NC "<<endl;
            for (ui j = 0; j < u_nbrs_count; ++j)
            {
                vertexpair = u_nbrs[j];
                // cout<<vertex<<" vertex " <<vertexpair<<" edge pair "<<endl;
                vx2 = checkA(ID, vertexpair, count);
                // cout<<vx1<<" vx1 " <<vx2<<" vx2 "<<endl;
                if (vx2 == 1000000)
                {
                    count++;
                    vx2 = count;
                    // cout<<"Indide "<<vx1<<" vx1 " <<vx2<<" vx2 "<<endl;
                    q_next.push(vertexpair);
                    ID[count] = u_nbrs[j];
                    M.insert(vx1, vx2) = 1.0;
                    M.insert(vx2, vx1) = 1.0;
                }
                else
                {
                    if (vx1 == vx2)
                        cout << "hola" << endl;
                    M.coeffRef(vx1, vx2) = 1.0;
                    M.coeffRef(vx2, vx1) = 1.0;
                }
            }
            // count=checkA(ID,vertex,count);
            // ui data_graph->offsets_[id + 1] - data_graph->offsets_[id];
            // for (ui j = 0; j < u_nbrs_count; ++j) {
            //   const VertexID* u_nbrs1 = data_graph->getVertexNeighbors(u_nbrs[j], u_nbrs_count1);
            // checkA(ID,vertex);
            //}
        }
        // q_temp=q_curr;
        // q_curr=q_next;
        // q_next=q_temp;
        // cout<<" I surrender "<<endl;
        if (!q_next.empty())
            q_curr.swap(q_next);
        else
        {
            i = depth;
        }
    }
    // cout<<" can we surrender "<<endl;
}
void LaplSparce(SparseMatrix<double> M, SparseMatrix<double> &L)
{
    int count = 0;
    for (int k = 0; k < M.outerSize(); ++k)
    {
        for (SparseMatrix<double>::InnerIterator it(M, k); it; ++it)
        {
            L.insert(it.row(), it.col()) = -1.0;
            count++;
            // row index
            // col index (here it is equal to k)
            // it.index(); // inner index, here it is equal to it.row()
        }
        L.insert(k, k) = count;
        count = 0;
    }
}
void calcEigens(SparseMatrix<double> M, int k, VectorXcd &evalues)
{
    // cout<<"hi"<<endl;

    int sizeM = M.nonZeros();
    int nev = 2 * k;
    if (sizeM == 0 || sizeM == 1)
    {
        evalues.resize(4);

        evalues << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
        return;
    }

    if (k > sizeM)
    {
        k = sizeM - 2;
        nev = k + 1;
    }

    if (nev > sizeM)
        nev = sizeM;
    // cout<<"hi"<<endl;
    SparseGenMatProd<double> op(M);
    GenEigsSolver<SparseGenMatProd<double>> eigs(op, k, k + 2);
    eigs.init();
    int nconv = eigs.compute(SortRule::LargestReal);

    // Retrieve results
    // cout<<"hi"<<endl;
    if (eigs.info() == CompInfo::Successful)
        evalues = eigs.eigenvalues();
    if (evalues.size() < 10)
    {
        int sz = evalues.size();
        // cout<<sz<<endl;
        // cout <<evalues<<endl;
        // cout <<" ev"<<endl;
        evalues.conservativeResize(10);
        // cout <<evalues<<endl;
        // cout <<" ev1"<<endl;
        for (int i = sz; i < 10; i++)
            evalues(i) = 0.0;
        // cout <<" ev2"<<endl;
        // cout <<evalues<<endl;
    }
    // std::cout << "Eigenvalues found:\n" << evalues << std::endl;
}
void MTcalc(Graph *data_graph, int degree, MatrixXcd &eigenVD)
{

    auto AdJAdl = [](Graph *data_graph, int degree, VertexID svertex, VertexID evertex, MatrixXcd &eigenVD)
    {
        SparseMatrix<double> M(4000, 4000);
        VectorXcd evalues;
        int depth = 2;
        // cout<<"hi"<<endl;
        for (int i = svertex; i < evertex; i++)
        {
            ExtractAdj(M, data_graph, degree, depth, i);
            // cout<<"hi"<<endl;
            if (i == 9458 || i == 9459)
                cout << "stop" << M.nonZeros() << endl;
            calcEigens(M, 10, evalues);
            // cout<<"hi"<<endl;
            // std::cout << "Matrix:\n" << M.row(1) << std::endl;
            // calcEigens(M,4,eigenVD[i]);
            M.setZero();
            eigenVD.row(i) = evalues;
            evalues.setZero();
            // std::cout << "Eigenvalues found:\n" <<  M << std::endl;
        }
    };
    // VectorXcd evalues;
    int Tnum = 1;
    int siz = data_graph->getVerticesCount();
    // siz=10;
    int div = (int)(siz / Tnum);
    thread th[Tnum];
    for (int d = 0; d < Tnum - 1; d++)
    {
        th[d] = thread(AdJAdl, data_graph, degree, div * d, div * (d + 1), ref(eigenVD));
    }
    th[Tnum - 1] = thread(AdJAdl, data_graph, degree, div * (Tnum - 1), siz, ref(eigenVD));
    // th[Tnum-1]=thread (AdJAdl,data_graph,degree,siz-2,siz,ref(eigenVD));
    for (int d = 0; d < Tnum; d++)
        th[d].join();
    // cout<<"hi"<<endl;
    // std::cout << "Eigenvalues found:\n" <<  eigenVD.row(0) << std::endl;
}
void calcEigensD(MatrixXd M, int k, VectorXd &evalues)
{
    int sizeM1 = M.nonZeros();
    int sizeM = M.size();
    int Tsize = (int)sqrt(sizeM);
    int nev = 2 * k;
    if (sizeM1 < 10000)
    {
        if (sizeM1 == 0 || sizeM1 == 1)
        {
            evalues.resize(k);
            evalues << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

            return;
        }

        if (k >= Tsize)
        {
            k = Tsize - 1;
            nev = Tsize;
            // nev=k+1;
        }
        if (k < 1)
            k = 1;

        if (nev > (Tsize))
            nev = (Tsize);
    }

    DenseGenMatProd<double> op(M);

    SymEigsSolver<DenseGenMatProd<double>> eigs(op, k, nev);
    eigs.init();
    int nconv = eigs.compute(SortRule::LargestAlge);

    if (eigs.info() == CompInfo::Successful)
        evalues = eigs.eigenvalues();
    if (evalues.size() < 10)
    {
        int sz = evalues.size();

        evalues.conservativeResize(10);

        evalues(sz) = 0;
        sz++;
        for (int i = sz; i < 10; i++)
            evalues(i) = -1;
    }
}