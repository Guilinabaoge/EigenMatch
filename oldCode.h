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
inline bool SecHopEigenPNS(vector<pair<VertexID, VertexID>> &q_curr, VertexID *ID, 
VertexID *IDD, int *IDDLC, vector<vector<CSV>> &FCS, vector<T> &tripletList, vector<pair<ui, int>> &EvalNeigb, Graph *query_graph, int countD);
void getAllEdges(CSV &cvertex, vector<VertexID> &temp2, vector<pair<VertexID, VertexID>> &q_curr, vector<vector<CSV>> &FCS);
inline VertexID checkAN(vector<VertexID> ID, VertexID CID);
void removeVertexAndEgjes(vector<vector<CSV>> &FCS, CSV cvertex, int i, int deli);
inline bool checkNeigb(Graph *query_graph, Graph *data_graph, vector<pair<ui, int>> QueryNlabel, VertexID vid);

inline VertexID checkANS(VertexID *ID, VertexID CID, int sizA);
inline VertexID checkANSBS(VertexID *ID, VertexID CID, int sizA);
bool allcomp(Graph *data_graph, Graph *query_graph, int di, int qi);
bool Labelcomp(Graph *data_graph, Graph *query_graph, int di, int qi);
inline bool eigcomp(VectorXd &qevalues, VectorXd &devalues, int comp);
VertexID findInd(vector<vector<CSV>> &FCS, VertexID IDC, VertexID IDQ);
bool degcomp(Graph *data_graph, Graph *query_graph, int di, int qi);

bool Refinement(vector<vector<pair<ui, int>>> &QueryNlabel, vector<vector<pair<ui, int>>> &QueryNlabel2, vector<int> &Qnew, vector<vector<CSV>> &FCS, int qsiz, Graph *query_graph, MatrixXd &eigenVq1);
void printPrunRes(MatrixXd &eigenVD1, MatrixXd &eigenVq1, Graph *data_graph, Graph *query_graph, int sizq, int sizd);
bool OneHopEigen(CSV &cvertex, vector<VertexID> &ID, vector<VertexID> &IDL, vector<VertexID> &IDD, vector<T> &tripletList, vector<pair<ui, int>> &EvalNeigb, Graph *query_graph, int countD, vector<pair<VertexID, VertexID>> &q_curr);
void InitPrunCS(vector<vector<CSV>> &FCS, int qsiz, Graph *query_graph);
void ExtractNI2H(vector<vector<pair<ui, int>>> &QueryNlabel, vector<vector<pair<ui, int>>> &QueryNlabel2, vector<int> &Dnew, Graph *query_graph, int qsiz);
bool OneHopEigenPN(CSV &cvertex, vector<VertexID> &ID, vector<VertexID> &IDL, vector<VertexID> &IDD, vector<T> &tripletList, vector<pair<ui, int>> &EvalNeigb, vector<pair<ui, int>> &EvalNeigb2, Graph *query_graph, int countD, vector<pair<VertexID, VertexID>> &q_curr);
bool ReverseRefinementDS(vector<vector<pair<ui, int>>> &QueryNlabel, vector<vector<pair<ui, int>>> &QueryNlabel2, vector<int> &Qnew, vector<vector<CSV>> &FCS, int qsiz, Graph *query_graph, MatrixXd &eigenVq1);
inline VertexID checkANSBSADD(VertexID *&ID, VertexID CID, int sizA);

inline bool OneHopEigenFPN(CSV &cvertex, vector<VertexID> &ID, vector<VertexID> &IDL, vector<VertexID> &IDD, vector<T> &tripletList, vector<pair<ui, int>> &EvalNeigb, vector<pair<ui, int>> &EvalNeigb2, Graph *query_graph, int countD, vector<pair<VertexID, VertexID>> &q_curr);
inline bool OneHopEigenPN(CSV &cvertex, vector<VertexID> &ID, vector<VertexID> &IDL, vector<VertexID> &IDD, vector<T> &tripletList, vector<pair<ui, int>> &EvalNeigb, vector<pair<ui, int>> &EvalNeigb2, Graph *query_graph, int countD, vector<pair<VertexID, VertexID>> &q_curr);
bool ReverseRefinementOC(vector<vector<pair<ui, int>>> &QueryNlabel, vector<vector<pair<ui, int>>> &QueryNlabel2, vector<int> &Qnew, vector<vector<CSV>> &FCS, int qsiz, Graph *query_graph, MatrixXd &eigenVq1);
bool SecHopEigenPN(vector<pair<VertexID, VertexID>> &q_curr, vector<VertexID> &ID, vector<VertexID> &IDD, vector<vector<CSV>> &FCS, vector<T> &tripletList, vector<pair<ui, int>> &EvalNeigb, Graph *query_graph, int countD);
bool ReverseRefinement(vector<vector<pair<ui, int>>> &QueryNlabel, vector<vector<pair<ui, int>>> &QueryNlabel2, vector<int> &Qnew, vector<vector<CSV>> &FCS, int qsiz, Graph *query_graph, MatrixXd &eigenVq1);
void SecHopEigen(vector<pair<VertexID, VertexID>> &q_curr, vector<VertexID> &ID, vector<vector<CSV>> &FCS, vector<T> &tripletList);
void PrunRes(Graph *query_graph, Graph *data_graph, MatrixXd &eigenVD1, MatrixXd &eigenVq1, int sizd, int sizq, int EN, float avgPruned[], float avgPrunedC[],
             float avgDprun[], float avgcheck[], float avgcheck1[], float avgVL[], float avgmagkas[],
             float avgDL[], float avgEL[]);
bool ReverseRefinementS(VertexID *ID, VertexID *IDD, VertexID *IDL, vector<vector<pair<ui, int>>> &QueryNlabel, vector<vector<pair<ui, int>>> &QueryNlabel2, vector<int> &Qnew, vector<vector<CSV>> &FCS, int qsiz, Graph *query_graph, MatrixXd &eigenVq1);
bool ReverseRefinementOCS(VertexID *ID, VertexID *IDD, vector<vector<pair<ui, int>>> &QueryNlabel,
                          vector<vector<pair<ui, int>>> &QueryNlabel2, vector<int> &Qnew, vector<vector<CSV>> &FCS,
                          int qsiz, Graph *query_graph, MatrixXd &eigenVq1);
void MTEdgesCS(vector<vector<CSV>> &FCS,int qsiz,int dsiz,Graph *data_graph,Graph *query_graph);
void edgesHelperCSPD(vector<vector<CSV>> &FCS,const VertexID u_nbrs[],const VertexID u_nbrsD[],
ui u_nbrs_count,ui u_nbrs_countD,int a,int b,VertexID VID);
void edgesHelperCS(vector<vector<CSV>> &FCS,const VertexID u_nbrs[],const VertexID u_nbrsD[],
ui u_nbrs_count,ui u_nbrs_countD,int a,int b);
void EdgesCS(vector<vector<CSV>> &FCS,int qsiz,int dsiz,Graph *data_graph,Graph *query_graph);
void VerticesCS(vector<vector<CSV>> &FCS,int qsiz,int dsiz,Graph *data_graph,Graph *query_graph,MatrixXd &eigenVD1,MatrixXd &eigenVq1);

void MTVerticesCS(vector<vector<CSV>> &FCS,int qsiz,int dsiz,Graph *data_graph,
Graph *query_graph,MatrixXd &eigenVD1,MatrixXd &eigenVq1);

void VerticesCSMN(ui **&candidates, ui *&candidates_count,int qsiz,
int dsiz,Graph *data_graph,Graph *query_graph,MatrixXd &eigenVD1,MatrixXd &eigenVq1);

void VerticesCSM(ui **&candidates, ui *&candidates_count,int qsiz,int dsiz,Graph *data_graph,
Graph *query_graph,MatrixXd &eigenVD1,MatrixXd &eigenVq1);
void calcres(int sizd,MatrixXd eigenVD1,Graph *data_graph);
