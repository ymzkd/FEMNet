#include "BarElement.h"

std::ostream &operator<<(std::ostream &os, const BeamStressData &bsd)
{
    /*std::string force = std::format("Nx: {}, Qy: {}, Qz: {}, Mx: {}, My: {}, Mz: {}",
        bsd.Nx, bsd.Qy, bsd.Qz, bsd.Mx, bsd.My, bsd.Mz);
    os << force << std::endl;*/

    os << "Nx: " << bsd.Nx << ", Qy: " << bsd.Qy << ", Qz: " << bsd.Qz
       << ", Mx: " << bsd.Mx << ", My: " << bsd.My << ", Mz: " << bsd.Mz << std::endl;
    return os;
}

std::ostream &operator<<(std::ostream &os, const BeamStress &bsd)
{
    os << "Beam Stress" << std::endl;
    os << "Start  " << bsd.S0;
    os << "End    " << bsd.S1;
    return os;
}

double BarElementBase::length()
{
    return Nodes[0]->Location.distance_to(Nodes[1]->Location);
}

void BarElementBase::AssembleMassMatrix(Eigen::SparseMatrix<double> &mat)
{
    int indices[6];
    Eigen::VectorXd mass = NodeLumpedMass();
    for (size_t i = 0; i < 3; i++)
    {
        int idx = Nodes[0]->id * 6 + i;
        mat.coeffRef(idx, idx) += mass[0];
    }
    for (size_t i = 0; i < 3; i++)
    {
        int idx = Nodes[1]->id * 6 + i;
        mat.coeffRef(idx, idx) += mass[1];
    }
}
