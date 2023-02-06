//#include <Eigen/Core>
#include <Eigen/SparseCore>
//#include <Spectra/GenEigsSolver.h>
//#include <Spectra/MatOp/SparseGenMatProd.h>
#include <iostream>
#include "GrM.h"
//#include "thread"
using  namespace Eigen;
VertexID checkANX(vector<VertexID> ID,VertexID CID);
bool checkM(SparseMatrix<double> mat);
void LaplSparce(SparseMatrix<double> M,SparseMatrix<double> &L);
VertexID checkA(VertexID *ID,VertexID vertex,int count);
void calcEigens(SparseMatrix<double> M,int k,VectorXcd &evalues);
void ExtractAdj(SparseMatrix<double> &M,Graph *data_graph,int degree,int depth,VertexID vertex);
void ExtractAdL(SparseMatrix<double> &M,Graph *data_graph,int degree,int depth,VertexID vertex);
void MTcalc(Graph *data_graph,int degree, MatrixXcd &eigenVD);
void MTcalc1(Graph *data_graph,int degree, MatrixXd &eigenVD);
//void calcEigens1( SparseMatrix<double> M,int k,VectorXd &evalues);
void calcEigens1(SparseMatrix<double> M, int k, VectorXd &evalues,int count);
 void MTEigCal(Graph *data_graph,int degree,MatrixXd &eigenVD,bool LE,int Eprun);
 void CompactADLEIG(Graph *data_graph,int degree,VectorXd &evalues,VertexID vertex,int depth,int Eprun);
 void CompactADJEIG(Graph *data_graph,int degree,VectorXd& evalues,VertexID vertex,int depth);
