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
    double Lambda_bzi, Lambda_bzj, Lambda_syi, Lambda_syj;
    double Lambda_byi, Lambda_byj, Lambda_szi, Lambda_szj;
    double lzi, lzj, lyi, lyj;

	// z軸周り, y方向たわみ, i端せん断バネ剛性の取得
    double get_ksyi();
    
    // z軸周り, y方向たわみ, i端回転バネ剛性の取得
    double get_kbzi();

	// z軸周り, y方向たわみ, j端せん断バネ剛性の取得
    double get_ksyj();

	// z軸周り, y方向たわみ, j端回転バネ剛性の取得
    double get_kbzj();

	// y軸周り, z方向たわみ, i端せん断バネ剛性の取得
    double get_kszi();

	// y軸周り, z方向たわみ, i端回転バネ剛性の取得
    double get_kbyi();

	// y軸周り, z方向たわみ, j端せん断バネ剛性の取得
    double get_kszj();

	// y軸周り, z方向たわみ, j端回転バネ剛性の取得
    double get_kbyj();


	// z軸周り, y方向たわみ, i端せん断バネ剛性の設定
    void set_ksyi(double ksyi);
    
	// z軸周り, y方向たわみ, i端回転バネ剛性の設定
    void set_kbzi(double kbzi);

	// z軸周り, y方向たわみ, j端せん断バネ剛性の設定
    void set_ksyj(double ksyj);

	// z軸周り, y方向たわみ, j端回転バネ剛性の設定
    void set_kbzj(double kbzj);

	// y軸周り, z方向たわみ, i端せん断バネ剛性の設定
    void set_kszi(double kszi);

	// y軸周り, z方向たわみ, i端回転バネ剛性の設定
    void set_kbyi(double kbyi);

	// y軸周り, z方向たわみ, j端せん断バネ剛性の設定
    void set_kszj(double kszj);

	// y軸周り, z方向たわみ, j端回転バネ剛性の設定
    void set_kbyj(double kbyj);

    ComplexBeamElement();
    ComplexBeamElement(Node *n0, Node *n1, Section *sec, Material mat, double beta = 0);
    ComplexBeamElement(int _id, Node *n0, Node *n1, Section *sec, Material mat, double beta = 0);

    Eigen::MatrixXd StiffnessMatrix() override;

    // Eigen::MatrixXd NodeConsistentMass();
    Displacement DisplaceAt(Displacement d0, Displacement d1, double p) override;
};

#endif