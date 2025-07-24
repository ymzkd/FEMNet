#ifndef SEISMIC_MODULE_H
#define SEISMIC_MODULE_H

#ifndef SWIGCSHARP
#include <vector>
#endif

void ResponseMax(const std::vector<double> &accels,
                 const double dt, const double damping_ratio, const double T,
                 double &max_disp, double &max_accel, double &max_vel);


void ResponseSpectrum(const std::vector<double> &accels, const std::vector<double> &periods,
                      const double dt, const double damping_ratio,
                      std::vector<double> &max_disp, std::vector<double> &max_accel, std::vector<double> &max_vel);


void ResponseSpectrum(const std::vector<double> &accels, const double Tmin, const double Tmax, const int n,
                      const double dt, const double damping_ratio,
                      std::vector<double> &max_disp, std::vector<double> &max_accel, std::vector<double> &max_vel);

#endif // SEISMIC_MODULE_H
