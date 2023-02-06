#include <iostream>
#include<vector>
#include <Eigen/Core>
#include "Grm.h"
using namespace std;
using  namespace Eigen;
 class CSV{
    public:  
//VectorXd eigenvalues;
VertexID ID;
vector<pair<VertexID,VertexID>> edges;
VertexID *edgesQ;
VertexID *edgesV;
bool change=true;
bool Ichange=true;
bool IPchange=true;
bool deleted=false;
public:  
CSV(int eigens,VertexID IDV){
//uniquelabels=new int(LBLSize);
//eigenvalues(eigens);
//eigenvalues.resize(eigens);
ID=IDV;
}
CSV(int eigens,VertexID IDV,ui maxDeg){

//eigenvalues.resize(eigens);
ID=IDV;
//edges.reserve(maxDeg);
deleted=false;

}
CSV(int eigens,ui IDV,ui maxDeg,ui MaxQDeg){
//eigenvalues.resize(eigens);
ID=IDV;
//edgesQ =new VertexID[maxDeg*MaxQDeg];
//edgesV =new VertexID[maxDeg*MaxQDeg];
//memset(edgesQ, 0, sizeof(ui) * maxDeg*MaxQDeg);
//memset(edgesV, 0, sizeof(ui) * maxDeg*MaxQDeg);
//change=true;
//Ichange=true;
}
CSV(ui IDV,ui totalD){
ID=IDV;
//edgesQ =new VertexID[totalD];
//edgesV =new VertexID[totalD];
//memset(edgesQ, 0, sizeof(ui) * totalD);
//memset(edgesV, 0, sizeof(ui) * totalD);

}
CSV(){
    ID=0;
    deleted=true;
   // eigenvalues.resize(10);
    //edges.reserve(500);
}
};

