#ifndef _TENSION_TRUSS_ELEMENT_
#define _TENSION_TRUSS_ELEMENT_

#include "TrussElement.h"

// 状態依存要素のミックスイン抽象クラス。
// 反復解析で内部状態(引張のみ負担、圧縮のみ負担、スラック等)に応じた
// 接線剛性・接線応力を返す Bar系要素の共通インターフェース。
// 状態(有効フラグ)は要素自身では保持せず、解析Operator側が所有して
// 各メソッドへ明示的に渡す。要素を状態レスに保つことで、同一モデルを
// 複数の荷重ケースで安全に共有できる。
class IStateDependentElement
{
public:
    virtual ~IStateDependentElement() = default;

    // 指定した状態(active)での接線剛性を Triplet 形式で組立てる
    virtual void GetTangentStiffnessTriplets(
        std::vector<Eigen::Triplet<double>> &triplets, bool active) = 0;

    // 変位ベクトルと現在状態から次の状態を判定して返す(要素は変更しない)。
    // 戻り値が current と異なる場合、上位の反復ループは剛性を再構築する。
    virtual bool NextState(const std::vector<Displacement> &disp, bool current) = 0;

    // 指定した状態(active)での接線応力を返す
    virtual BeamStress tangent_stress(Displacement d0, Displacement d1, bool active) = 0;
};

class TensionTrussElement : public TrussElement, public IStateDependentElement
{
public:
    double CutoffTension = 0;
    double ReductionFactor = 0.0001;

    TensionTrussElement() {}
    TensionTrussElement(Node *n0, Node *n1, Section *sec, Material mat)
        : TrussElement(n0, n1, sec, mat) {}
    TensionTrussElement(int _id, Node *n0, Node *n1, Section *sec, Material mat)
        : TrussElement(_id, n0, n1, sec, mat) {}

    // 指定した状態での非抗圧性を反映したトラス要素の剛性行列(6x6)
    Eigen::MatrixXd TangentStiffnessMatrix(bool active);

    // === IStateDependentElement の実装 ===

    // 非抗圧性を反映したTriplet方式での行列組立
    void GetTangentStiffnessTriplets(
        std::vector<Eigen::Triplet<double>> &triplets, bool active) override;

    // 接線応力の軸力を CutoffTension と比較して次状態を返す
    bool NextState(const std::vector<Displacement> &disp, bool current) override;

    // 指定した状態(active)での応力計算
    BeamStress tangent_stress(Displacement d0, Displacement d1, bool active) override;
};

#endif
