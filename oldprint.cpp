#include <Eigen/StdVector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "cs.h"
void printCSSize(vector<vector<CSV>> &FCS, int qsiz);
void printADJSparse(SparseMatrix<double> L);
void printCSSizeO(vector<vector<CSV>> &FCS, int qsiz);
float VarianceCalc(float val[], int valS)
void printCSSizeReal(vector<vector<CSV>> &FCS, int qsiz);
void printCSSizeReal(vector<vector<CSV>> &FCS, int qsiz);
void printEigenValues(MatrixXd &eigenVD1, int sizq, int ec);

{
    int count = 0;
    int Tcount = 0;
    for (int kk = 0; kk < qsiz; kk++)
    {
        for (int zz = 0; zz < FCS[kk].size(); zz++)
            // if(FCS[kk][zz].ID!=1000000)
            if (FCS[kk][zz].deleted == false)
                count++;
        Tcount = count + Tcount;
        cout << "i= " << kk << " size " << count << endl;
        count = 0;
    }
    cout << "Total candidate size " << Tcount << endl;
}
{
    float sum = 0.0, mean, variance = 0.0, stdDeviation;
    int i;
    for (i = 0; i < valS; ++i)
        sum += val[i];
    mean = sum / valS;
    for (i = 0; i < valS; ++i)
        variance += pow(val[i] - mean, 2);
    variance = variance / valS;
    stdDeviation = sqrt(variance);
    return variance;
    // return stdDeviation;
}
void printCSSizeO(vector<vector<CSV>> &FCS, int qsiz)
{
    int count = 0;
    for (int kk = 0; kk < qsiz; kk++)
    {
        // cout<<"size "<<FCS[kk].size()<<endl;
        count = count + FCS[kk].size();
    }

    cout << "Total candidate size " << count << endl;
}
void printADJSparse(SparseMatrix<double> L)
{
    for (int k = 0; k < L.outerSize(); ++k)
    {
        for (SparseMatrix<double>::InnerIterator it(L, k); it; ++it)
        {
            if (it.value() > 0)
                cout << "(" << it.row() << "," << it.col() << ") " << it.index() << endl;
        }
        // cout<<"Laplaciann" <<endl;
    }
}void printCSSize(vector<vector<CSV>> &FCS, int qsiz)
{
    int count = 0;

    for (int kk = 0; kk < qsiz; kk++)
    {
        cout << "i= " << kk << " size " << FCS[kk].size() << endl;
        count = count + FCS[kk].size();
    }

    cout << "Total candidate size " << count << endl;
}