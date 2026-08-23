#ifndef _RESPONSE_NAMING_
#define _RESPONSE_NAMING_

// 応答レコーダ/サンプラーの既定名を組み立てる共通ヘルパ。
// 実装ファイルからのみ使用するため、SWIG のラップ対象には含めない。

#include <cmath>
#include <string>

#include "Components.h"

/// 応答種別に対応する名称の断片(Disp / Vel / Accel)
inline std::string ResponseValueTag(ResponseValueType vt)
{
    switch (vt)
    {
    case ResponseValueType::Velocity:     return "Vel";
    case ResponseValueType::Acceleration: return "Accel";
    default:                              return "Disp";
    }
}

/// 応答種別に対応する説明文の断片(displacement / velocity / acceleration)
inline std::string ResponseValueLabel(ResponseValueType vt)
{
    switch (vt)
    {
    case ResponseValueType::Velocity:     return "velocity";
    case ResponseValueType::Acceleration: return "acceleration";
    default:                              return "displacement";
    }
}

/// 評価方向に対応する名称の断片
/// (零ベクトル = Abs、座標軸に一致すれば X/Y/Z、それ以外は Dir)
inline std::string ResponseDirectionTag(const Vector &direction)
{
    Vector dir = direction;
    double dir_norm = dir.norm();
    if (dir_norm <= 0.0)
        return "Abs";

    const double tol = 1e-9;
    Vector unit = Vector::multiply(dir, 1.0 / dir_norm);
    if (std::abs(std::abs(unit.x) - 1.0) <= tol) return "X";
    if (std::abs(std::abs(unit.y) - 1.0) <= tol) return "Y";
    if (std::abs(std::abs(unit.z) - 1.0) <= tol) return "Z";
    return "Dir";
}

#endif
