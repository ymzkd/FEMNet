#ifndef _COMPLEX_BEAM_ELEMENT_
#define _COMPLEX_BEAM_ELEMENT_

#include "BarElement.h"
#include "BeamElement.h"

class ComplexBeamElement : public BeamElement
{
private:
    // static const int total_dof = 12;
    Eigen::MatrixXd stiffness_matrix_local() override;
    // double element_length();
    // Eigen::MatrixXd trans_matrix();
public:
    double Lambda_bz, Lambda_bz_, Lambda_sy, Lambda_sy_;
    double Lambda_by, Lambda_by_, Lambda_sz, Lambda_sz_;
    double lzi, lzj, lyi, lyj;

    ComplexBeamElement();
    ComplexBeamElement(Node *n0, Node *n1, Section *sec, Material mat, double beta = 0);
    ComplexBeamElement(int _id, Node *n0, Node *n1, Section *sec, Material mat, double beta = 0);

    Eigen::MatrixXd StiffnessMatrix() override;

    // Eigen::MatrixXd NodeConsistentMass();
};

#endif