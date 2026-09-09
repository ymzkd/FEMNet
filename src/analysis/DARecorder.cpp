#include "DARecorder.h"
#include "FEDynamic.h"

#include <iostream>
#include <sstream>

#include "ResponseNaming.h"

void DARecorder_KineticEnergy::Record(DynamicAnalysis &da)
{
    double energy = 0.5 * da.current_vel.dot(da.matM_aa.selfadjointView<Eigen::Upper>() * da.current_vel);
    Values.push_back(energy);
}

void DARecorder_PotentialEnergy::Record(DynamicAnalysis &da)
{
    double energy = 0.5 * da.current_disp.dot(da.matK_aa.selfadjointView<Eigen::Upper>() * da.current_disp);
    Values.push_back(energy);
}

void DARecorder_DampingEnergy::Record(DynamicAnalysis &da)
{
    double energy = da.current_vel.dot(da.matC_aa.selfadjointView<Eigen::Upper>() * da.current_vel);
    Values.push_back(energy);
}

void DARecorder_InputEnergy::Record(DynamicAnalysis &da)
{
    // Record は current_step 更新後に呼ばれる。既存実装に合わせ 1つ前のステップの
    // 外力を参照する。入力エネルギーは現行の符号慣行に合わせ -F_ext·v を積算する
    // (地震では F_ext=-M·ι·a_g なので (M·ι·a_g)·v となり従来と一致)。
    Eigen::VectorXd f_reduced, f_fix;
    da.ReducedLoadVector(da.TimeAt(da.current_step - 1), f_reduced, f_fix);
    Values.push_back(-f_reduced.dot(da.current_vel));
}

DARecorder_NodeResponse::DARecorder_NodeResponse(int node_id, Vector direction, ResponseValueType value_type)
    : NodeId(node_id), Direction(direction), ValueType(value_type)
{
    Name = DefaultName();
    Description = DefaultDescription();
}

std::string DARecorder_NodeResponse::DefaultName() const
{
    std::ostringstream oss;
    oss << "Node" << ResponseValueTag(ValueType) << ".N" << NodeId
        << "." << ResponseDirectionTag(Direction);
    return oss.str();
}

std::string DARecorder_NodeResponse::DefaultDescription() const
{
    std::ostringstream oss;
    oss << "Nodal " << ResponseValueLabel(ValueType) << " of node " << NodeId
        << " along (" << Direction.x << ", " << Direction.y << ", " << Direction.z << ")";
    return oss.str();
}

void DARecorder_NodeResponse::Initialize(DynamicAnalysis &da)
{
    DARecorder::Initialize(da);

    // 節点番号と方向の妥当性を確認する。無効な場合でも記録数を他のレコーダと
    // 揃えるため、Record では 0 を積み続ける。
    int node_num = static_cast<int>(da.model ? da.model->Nodes.size() : 0);
    double dir_norm = Direction.norm();
    valid_ = (NodeId >= 0 && NodeId < node_num && dir_norm > 0.0);
    unit_ = valid_ ? Vector::multiply(Direction, 1.0 / dir_norm) : Vector();

    if (!valid_)
    {
        std::cerr << "DARecorder_NodeResponse: invalid node id (" << NodeId
                  << ") or zero direction. Zero values will be recorded." << std::endl;
    }

    if (Name.empty())
        Name = DefaultName();
    if (Description.empty())
        Description = DefaultDescription();

    // ステップ0の値を記録する(初期加速度は Initialize 内で計算済みのため 0 とは限らない)
    Values.push_back(CurrentValue(da));
}

void DARecorder_NodeResponse::Record(DynamicAnalysis &da)
{
    Values.push_back(CurrentValue(da));
}

double DARecorder_NodeResponse::CurrentValue(DynamicAnalysis &da)
{
    if (!valid_)
        return 0.0;

    std::vector<Displacement> response;
    switch (ValueType)
    {
    case ResponseValueType::Velocity:     response = da.GetVelocities(); break;
    case ResponseValueType::Acceleration: response = da.GetAccelerations(); break;
    default:                              response = da.GetDisplacements(); break;
    }

    if (NodeId >= static_cast<int>(response.size()))
        return 0.0;

    // 並進成分の指定方向成分(符号付き)
    return Vector::multiply(response[NodeId].Translation(), unit_);
}
