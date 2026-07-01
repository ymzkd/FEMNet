#include "SeismicModule.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


void ResponseMax(const std::vector<double> &accels, 
    const double dt, const double damping_ratio, const double T,
    double &max_disp, double &max_accel, double &max_vel) {
    
    const double omega = 2.0 * M_PI / T;  // 固有角振動数
    
    // Newmarkのβ法のパラメータ
    const double beta = 0.25;
    
    // 減衰係数
    const double cm = 2.0 * damping_ratio * omega;
    const double km = omega * omega;

    // 初期化
    double u = 0.0;      // 変位
    double v = 0.0;      // 速度
    double a = 0.0;      // 加速度
    
    max_disp = 0.0;
    max_accel = 0.0;
    max_vel = 0.0;
    
    // 時刻歴応答解析
    for (size_t i = 0; i < accels.size(); ++i) {
        double an1_up = accels[i] + cm * (v + a * 0.5 * dt) + km * (u + v * dt + (0.5 - beta) * a * dt * dt);
        double an1_dw = 1.0 + 0.5 * cm * dt + beta * km * dt * dt;
        double an1 = -an1_up / an1_dw;

        const double vn1 = v + (a + an1) * 0.5 * dt;
        const double un1 = u + v * dt + (0.5 - beta) * a * dt * dt + beta * an1 * dt * dt;

        // 最大値の更新
        max_disp = std::max(max_disp, std::abs(un1));
        max_accel = std::max(max_accel, std::abs(an1));
        max_vel = std::max(max_vel, std::abs(vn1));
        
        // 次ステップに更新
        u = un1;
        v = vn1;
        a = an1;
    }
}

void ResponseSpectrum(const std::vector<double> &accels, const std::vector<double> &periods,
                      const double dt, const double damping_ratio,
                      std::vector<double> &max_disp, std::vector<double> &max_accel, std::vector<double> &max_vel) {
    
    // 出力ベクトルのサイズを周期数に合わせて初期化
    const size_t n_periods = periods.size();
    max_disp.resize(n_periods);
    max_accel.resize(n_periods);
    max_vel.resize(n_periods);
    
    // 各周期に対して応答解析を実行
    for (size_t i = 0; i < n_periods; ++i) {
        const double T = periods[i];
        
        // 各周期に対してResponseMax関数を呼び出し
        ResponseMax(accels, dt, damping_ratio, T,
                   max_disp[i], max_accel[i], max_vel[i]);
    }
}

void ResponseSpectrum(const std::vector<double> &accels, const double Tmin, const double Tmax, const int n,
                      const double dt, const double damping_ratio,
                      std::vector<double> &max_disp, std::vector<double> &max_accel, std::vector<double> &max_vel) {
    
    // 周期のベクトルを生成
    std::vector<double> periods(n + 1);
    
    // 線形分割で周期を生成
    const double dT = (Tmax - Tmin) / n;
    for (int i = 0; i <= n; ++i) {
        periods[i] = Tmin + i * dT;
    }

    // 既存のResponseSpectrum関数を呼び出し
    ResponseSpectrum(accels, periods, dt, damping_ratio, max_disp, max_accel, max_vel);
}

