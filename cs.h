#include <iostream>
#include <vector>
#include <Eigen/Core>
#include "Grm.h"
using namespace std;
using namespace Eigen;
class CSV
{
public:
    // VectorXd eigenvalues;
    VertexID ID;
    vector<pair<VertexID, VertexID>> edges;
    VertexID *edgesQ;
    VertexID *edgesV;
    bool change = true;
    bool Ichange = true;
    bool IPchange = true;
    bool deleted = false;

public:
    CSV(int eigens, VertexID IDV)
    {

        ID = IDV;
    }
    CSV(int eigens, VertexID IDV, ui maxDeg)
    {

        ID = IDV;

        deleted = false;
    }
    CSV(int eigens, ui IDV, ui maxDeg, ui MaxQDeg)
    {

        ID = IDV;
    }
    CSV(ui IDV, ui totalD)
    {
        ID = IDV;
    }
    CSV()
    {
        ID = 0;
        deleted = true;
    }
};
