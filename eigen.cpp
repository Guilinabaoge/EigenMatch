#include <Eigen/Dense>
#include <iostream>
#include <Eigen/Core>
#include <Spectra/SymEigsSolver.h>
using  namespace Eigen;
using namespace std;
int main(){
    //dynamic
    Eigen::MatrixXd d;
    //fix size Matrix
    Matrix4d f;
    f=Matrix4d::Random();
    Eigen::MatrixXd m1(2,3);
    cout<<"hi "<<m1.rows()<<m1.cols()<<endl;
    //resize
    //d.resize(4,4);
}