#ifndef _TENSION_TRUSS_ELEMENT_
#define _TENSION_TRUSS_ELEMENT_

#include "TrussElement.h"

// 状態依存要素のミックスイン抽象クラス。
// 反復解析で内部状態(引張のみ負担、圧縮のみ負担、スラック等)を更新し、
// 現在状態に応じた接線剛性・接線応力を返す Bar系要素の共通インターフェース。
class IStateDependentElement
{
public:
    // 規定剛性で機能している(true)か、縮退状態(false)か。
    // 反復収束後の最終状態を保持し、解析開始時は外部から true を代入してリセットする。
    bool IsActive = true;

    virtual ~IStateDependentElement() = default;

    // 現在の内部状態に応じた接線剛性を Triplet 形式で組立てる
    virtual void GetTangentStiffnessTriplets(
        std::vector<Eigen::Triplet<double>> &triplets) = 0;

    // 変位ベクトルに基づき内部状態(IsActive)を更新する。
    // 状態変化があれば true を返し、上位の反復ループは剛性を再構築する。
    virtual bool update(const std::vector<Displacement> &disp) = 0;

    // 現在の内部状態を反映した接線応力を返す
    virtual BeamStress tangent_stress(Displacement d0, Displacement d1) = 0;
};

class TensionTrussElement : public TrussElement, public IStateDependentElement
{
private:
    // ローカル剛性行列に係数を掛ける(TangentStiffnessMatrix と tangent_stress の両方に反映される)
    Eigen::MatrixXd tangentstiffness_matrix_local();

public:
    double CutoffTension = 0;
    double ReductionFactor = 0.0001;

    TensionTrussElement() {}
    TensionTrussElement(Node *n0, Node *n1, Section *sec, Material mat)
        : TrussElement(n0, n1, sec, mat) {}
    TensionTrussElement(int _id, Node *n0, Node *n1, Section *sec, Material mat)
        : TrussElement(_id, n0, n1, sec, mat) {}

    // 非抗圧性を反映したトラス要素の剛性行列(6x6)
    Eigen::MatrixXd TangentStiffnessMatrix();

    // === IStateDependentElement の実装 ===

    // 非抗圧性を反映したTriplet方式での行列組立
    void GetTangentStiffnessTriplets(
        std::vector<Eigen::Triplet<double>> &triplets) override;

    // tangent_stress を見て IsActive 状態を更新する。状態変化があれば true を返す。
    bool update(const std::vector<Displacement> &disp) override;

    // 非抗圧性を反映した応力計算
    BeamStress tangent_stress(Displacement d0, Displacement d1) override;
};

#endif
