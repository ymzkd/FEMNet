#include "DARecorder.h"
#include "FEDynamic.h"

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
