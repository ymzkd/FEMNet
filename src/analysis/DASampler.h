#ifndef _DASAMPLER_
#define _DASAMPLER_

#ifndef SWIG
#include <vector>
#include <string>
#endif

#include "Components.h"
#include "LoadComponent.h"

// 前方宣言
class DynamicAnalysis;

/// <summary>
/// 時刻歴応答解析の各ステップで応答を評価し、着目量が最大となるステップの
/// 状態量(変位・速度・加速度)を保持するサンプラーの基底クラス。
/// DynamicAnalysis へ複数登録でき、C#/Python 側でも派生できる。
/// </summary>
class DASampler
{

public:
    virtual ~DASampler() = default;

    int step;
    std::string Name;
    std::vector<Displacement> velocity, displacement, acceleration;
    std::vector<NodeLoad> react_force;

    DASampler() : step(0), Name("") {}
    DASampler(std::string name) : step(0), Name(name) {}

    virtual void Sampling(DynamicAnalysis &analysis) = 0;

    /// <summary>
    /// 記録時点の応答(変位・速度・加速度・反力)を保存する
    /// </summary>
    void CaptureState(DynamicAnalysis &analysis);
};

class DASampler_MaxDisplacement : public DASampler
{
public:
    double max_displacement = 0.0;

    DASampler_MaxDisplacement() : DASampler("MaxDisplacement") {}

    void Sampling(DynamicAnalysis &da) override;
};

class DASampler_MaxVelocity : public DASampler
{
public:
    double max_velocity = 0.0;

    DASampler_MaxVelocity() : DASampler("MaxVelocity") {}

    void Sampling(DynamicAnalysis &da) override;
};

class DASampler_MaxAcceleration : public DASampler
{
public:
    double max_acceleration = 0.0;

    DASampler_MaxAcceleration() : DASampler("MaxAcceleration") {}

    void Sampling(DynamicAnalysis &da) override;
};

class DASampler_MaxDispDirection : public DASampler
{
public:

    double max_displacement = 0.0;
    Vector direction;

    DASampler_MaxDispDirection() : DASampler("MaxDispDirection") {}
    DASampler_MaxDispDirection(Vector direction) : DASampler("MaxDispDirection")
    {
        this->direction = direction;
    }

    void Sampling(DynamicAnalysis &da) override;
};

class DASampler_MaxVelocityDirection : public DASampler
{
public:

    double max_velocity = 0.0;
    Vector direction;

    DASampler_MaxVelocityDirection() : DASampler("MaxVelocityDirection"){}
    DASampler_MaxVelocityDirection(Vector direction) : DASampler("MaxVelocityDirection")
    {
        this->direction = direction;
    }

    void Sampling(DynamicAnalysis &da) override;
};

class DASampler_MaxAccelDirection : public DASampler
{
public:
    double max_accel = 0.0;
    Vector direction;

    DASampler_MaxAccelDirection() : DASampler("MaxAccelDirection"){}
    DASampler_MaxAccelDirection(Vector direction) : DASampler("MaxAccelDirection")
    {
        this->direction = direction;
    }

    void Sampling(DynamicAnalysis &da) override;
};

/// <summary>
/// 全支点の並進反力の合力(ベースシア)の大きさが最大となるステップを記録する。
/// max_force には最大時の合力ベクトルが入る(成分の内訳を見る用)。
/// </summary>
class DASampler_MaxBaseShear : public DASampler
{
public:
    double max_base_shear = 0.0;
    Vector max_force; // 最大時の合力ベクトル

    DASampler_MaxBaseShear() : DASampler("MaxBaseShear") {}

    void Sampling(DynamicAnalysis &da) override;
};

/// <summary>
/// 全支点の並進反力の合力(ベースシア)の指定方向成分が最大となるステップを記録する。
/// </summary>
class DASampler_MaxBaseShearDirection : public DASampler
{
public:
    double max_base_shear = 0.0;
    Vector direction;
    Vector max_force; // 最大時の合力ベクトル

    DASampler_MaxBaseShearDirection() : DASampler("MaxBaseShearDirection") {}
    DASampler_MaxBaseShearDirection(Vector direction) : DASampler("MaxBaseShearDirection")
    {
        this->direction = direction;
    }

    void Sampling(DynamicAnalysis &da) override;
};

#endif
