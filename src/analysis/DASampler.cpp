#include "DASampler.h"
#include "FEDynamic.h"

#include <cmath>

#include "ResponseNaming.h"

// 記録時点の応答をまとめて保存する
void DASampler::CaptureState(DynamicAnalysis &da)
{
    velocity = da.GetVelocities();
    displacement = da.GetDisplacements();
    acceleration = da.GetAccelerations();
    react_force = da.GetReactForces();
}

DASampler_MaxResponse::DASampler_MaxResponse(ResponseValueType value_type)
    : ValueType(value_type)
{
    Name = DefaultName();
}

DASampler_MaxResponse::DASampler_MaxResponse(ResponseValueType value_type, Vector direction)
    : ValueType(value_type), Direction(direction)
{
    Name = DefaultName();
}

std::string DASampler_MaxResponse::DefaultName() const
{
    return "Max" + ResponseValueTag(ValueType) + "." + ResponseDirectionTag(Direction);
}

void DASampler_MaxResponse::Sampling(DynamicAnalysis &da)
{
    // 既定コンストラクタ + プロパティ設定で生成された場合に備え、初回に名称を補う
    if (Name.empty())
        Name = DefaultName();

    // 評価対象の応答量を取得する
    std::vector<Displacement> response;
    switch (ValueType)
    {
    case ResponseValueType::Velocity:     response = da.GetVelocities();    break;
    case ResponseValueType::Acceleration: response = da.GetAccelerations(); break;
    default:                              response = da.GetDisplacements(); break;
    }

    // 方向が未設定(零ベクトル)なら大きさ、設定されていればその方向成分(絶対値)で評価する
    Vector dir = Direction;
    double dir_norm = dir.norm();
    bool directional = (dir_norm > 0.0);
    Vector unit = directional ? Vector::multiply(dir, 1.0 / dir_norm) : Vector();

    bool updated = false;

    // 任意節点の評価値が最大となるステップを記録
    for (size_t i = 0; i < da.model->Nodes.size(); i++)
    {
        Vector t = response[i].Translation();
        double value = directional
            ? std::abs(Vector::multiply(t, unit))
            : t.norm();
        if (value > MaxValue)
        {
            MaxValue = value;
            step = da.current_step;
            updated = true;
        }
    }

    if (updated)
    {
        CaptureState(da);
    }
}

// 全支点の並進反力の合力を求める(支点が無い場合は零ベクトル)
static Vector SumReactionForce(std::vector<NodeLoad> &reactions)
{
    Vector sum;
    for (auto &r : reactions)
        sum = sum + Vector(r.Px(), r.Py(), r.Pz());
    return sum;
}

void DASampler_MaxBaseShear::Sampling(DynamicAnalysis &da)
{
    std::vector<NodeLoad> reactions = da.GetReactForces();
    if (reactions.empty())
        return;

    // 合力(ベースシア)の大きさが最大となるステップを記録
    Vector sum = SumReactionForce(reactions);
    double length = sum.norm();
    if (length <= max_base_shear)
        return;

    max_base_shear = length;
    max_force = sum;
    step = da.current_step;

    CaptureState(da);
}

void DASampler_MaxBaseShearDirection::Sampling(DynamicAnalysis &da)
{
    // 方向が未設定(零ベクトル)の場合は評価できないためサンプリングしない
    Vector dir = direction;
    double dir_norm = dir.norm();
    if (dir_norm <= 0.0)
        return;
    Vector unit = Vector::multiply(dir, 1.0 / dir_norm);

    std::vector<NodeLoad> reactions = da.GetReactForces();
    if (reactions.empty())
        return;

    // 合力の指定方向成分(絶対値)が最大となるステップを記録
    Vector sum = SumReactionForce(reactions);
    double f_dir = std::abs(Vector::multiply(sum, unit));
    if (f_dir <= max_base_shear)
        return;

    max_base_shear = f_dir;
    max_force = sum;
    step = da.current_step;

    CaptureState(da);
}
