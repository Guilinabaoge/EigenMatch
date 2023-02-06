#include <Eigen/SparseCore>
bool tripletCompare(const Eigen::Triplet<double> &t1, const Eigen::Triplet<double> &t2)
{
    if (t1.row() != t2.row())
        return t1.row() < t2.row();
    return t1.col() < t2.col();
}
bool tripletUnique(const Eigen::Triplet<double> &t1, const Eigen::Triplet<double> &t2)
{
    if (t1.row() == t2.row())
        if (t1.col() == t2.col())
            return true;

    return false;
}
