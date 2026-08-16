#include "DASampler.h"
#include "FEDynamic.h"

#include <cmath>

void DASampler_MaxDisplacement::Sampling(DynamicAnalysis &da)
{
    bool updated = false;
    std::vector<Displacement> disp = da.GetDisplacements();

    // 任意節点の変位が最大となるステップを記録
    for (int i = 0; i < da.model->Nodes.size(); i++)
    {
        double d_length = disp[i].Translation().norm();
        if (d_length > max_displacement)
        {
            max_displacement = d_length;
            step = da.current_step;
            updated = true;
        }
    }

    if (updated)
    {
        velocity = da.GetVelocities();
        displacement = disp;
        acceleration = da.GetAccelerations();
    }
}

void DASampler_MaxVelocity::Sampling(DynamicAnalysis &da)
{
    bool updated = false;
    std::vector<Displacement> vel = da.GetVelocities();

    // 任意節点の速度が最大となるステップを記録
    for (size_t i = 0; i < da.model->Nodes.size(); i++)
    {
        double v_length = vel[i].Translation().norm();
        if (v_length > max_velocity)
        {
            max_velocity = v_length;
            step = da.current_step;
            updated = true;
        }
    }

    if (updated)
    {
        velocity = vel;
        displacement = da.GetDisplacements();
        acceleration = da.GetAccelerations();
    }
}

void DASampler_MaxAcceleration::Sampling(DynamicAnalysis &da)
{
    bool updated = false;
    std::vector<Displacement> acc = da.GetAccelerations();

    // 任意節点の加速度が最大となるステップを記録
    for (size_t i = 0; i < da.model->Nodes.size(); i++)
    {
        double a_length = acc[i].Translation().norm();
        if (a_length > max_acceleration)
        {
            max_acceleration = a_length;
            step = da.current_step;
            updated = true;
        }
    }

    if (updated)
    {
        velocity = da.GetVelocities();
        displacement = da.GetDisplacements();
        acceleration = acc;
    }
}

void DASampler_MaxDispDirection::Sampling(DynamicAnalysis &da)
{
    // 方向が未設定(零ベクトル)の場合は評価できないためサンプリングしない
    Vector dir = direction;
    double dir_norm = dir.norm();
    if (dir_norm <= 0.0)
        return;
    Vector unit = Vector::multiply(dir, 1.0 / dir_norm);

    bool updated = false;
    std::vector<Displacement> disp = da.GetDisplacements();

    // 任意節点の指定方向並進変位成分(絶対値)が最大となるステップを記録
    for (size_t i = 0; i < da.model->Nodes.size(); i++)
    {
        double d_dir = std::abs(Vector::multiply(disp[i].Translation(), unit));
        if (d_dir > max_displacement)
        {
            max_displacement = d_dir;
            step = da.current_step;
            updated = true;
        }
    }

    if (updated)
    {
        velocity = da.GetVelocities();
        displacement = disp;
        acceleration = da.GetAccelerations();
    }
}

void DASampler_MaxVelocityDirection::Sampling(DynamicAnalysis &da)
{
    // 方向が未設定(零ベクトル)の場合は評価できないためサンプリングしない
    Vector dir = direction;
    double dir_norm = dir.norm();
    if (dir_norm <= 0.0)
        return;
    Vector unit = Vector::multiply(dir, 1.0 / dir_norm);

    bool updated = false;
    std::vector<Displacement> vel = da.GetVelocities();

    // 任意節点の指定方向並進速度成分(絶対値)が最大となるステップを記録
    for (size_t i = 0; i < da.model->Nodes.size(); i++)
    {
        double v_dir = std::abs(Vector::multiply(vel[i].Translation(), unit));
        if (v_dir > max_velocity)
        {
            max_velocity = v_dir;
            step = da.current_step;
            updated = true;
        }
    }

    if (updated)
    {
        velocity = vel;
        displacement = da.GetDisplacements();
        acceleration = da.GetAccelerations();
    }
}

void DASampler_MaxAccelDirection::Sampling(DynamicAnalysis &da)
{
    // 方向が未設定(零ベクトル)の場合は評価できないためサンプリングしない
    Vector dir = direction;
    double dir_norm = dir.norm();
    if (dir_norm <= 0.0)
        return;
    Vector unit = Vector::multiply(dir, 1.0 / dir_norm);

    bool updated = false;
    std::vector<Displacement> acc = da.GetAccelerations();

    // 任意節点の指定方向並進加速度成分(絶対値)が最大となるステップを記録
    for (size_t i = 0; i < da.model->Nodes.size(); i++)
    {
        double a_dir = std::abs(Vector::multiply(acc[i].Translation(), unit));
        if (a_dir > max_accel)
        {
            max_accel = a_dir;
            step = da.current_step;
            updated = true;
        }
    }

    if (updated)
    {
        velocity = da.GetVelocities();
        displacement = da.GetDisplacements();
        acceleration = acc;
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

    velocity = da.GetVelocities();
    displacement = da.GetDisplacements();
    acceleration = da.GetAccelerations();
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

    velocity = da.GetVelocities();
    displacement = da.GetDisplacements();
    acceleration = da.GetAccelerations();
}
