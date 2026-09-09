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

/// <summary>
/// 全節点の並進応答(変位・速度・加速度)を評価し、着目量が最大となるステップを記録する。
/// Direction が零ベクトルの場合は応答の大きさ |v| を、非零の場合はその方向成分の
/// 絶対値 |v·n| を評価する(方向ベクトルは内部で正規化して用いる)。
/// </summary>
class DASampler_MaxResponse : public DASampler
{
public:
    /// 評価する応答量の種別
    ResponseValueType ValueType = ResponseValueType::Displacement;
    /// 評価方向(零ベクトル = 大きさで評価)
    Vector Direction;
    /// 記録時点の評価値(最大値)
    double MaxValue = 0.0;

    DASampler_MaxResponse() = default;
    explicit DASampler_MaxResponse(ResponseValueType value_type);
    DASampler_MaxResponse(ResponseValueType value_type, Vector direction);

    /// 現在の設定から生成した既定の名称(例: MaxDisp.Abs / MaxVel.X / MaxAccel.Dir)
    std::string DefaultName() const;

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
