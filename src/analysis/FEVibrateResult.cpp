#include "FEVibrateResult.h"

std::vector<double> FEVibrateResult::ParticipationFactors()
{
    Eigen::SparseMatrix<double> mass_mat = model->AssembleMassMatrix();

    std::vector<double> participation_factors;

    for (size_t i = 0; i < ModeNum(); i++)
    {
        Eigen::VectorXd v = Eigen::VectorXd::Zero(model->NodeNum() * 6);
        for (size_t j = 0; j < model->NodeNum(); j++)
        {
            int p = j * 6;
            v(p) = mode_vectors[i][j].Dx();
            v(p + 1) = mode_vectors[i][j].Dy();
            v(p + 2) = mode_vectors[i][j].Dz();
            v(p + 3) = mode_vectors[i][j].Rx();
            v(p + 4) = mode_vectors[i][j].Ry();
            v(p + 5) = mode_vectors[i][j].Rz();
        }

        Eigen::VectorXd vM = mass_mat.selfadjointView<Eigen::Upper>() * v;
        double beta_i = Eigen::VectorXd::Ones(v.size()).dot(vM) / v.dot(vM);

        participation_factors.push_back(beta_i);
    }

    return participation_factors;
}

std::vector<double> FEVibrateResult::ParticipationFactors(Vector direction)
{
    // MassMatrixの組み立て
    Eigen::SparseMatrix<double> mass_mat = model->AssembleMassMatrix();

    std::vector<double> participation_factors;

    for (size_t i = 0; i < ModeNum(); i++)
    {
        Eigen::VectorXd v = Eigen::VectorXd::Zero(model->NodeNum() * 6);
        Eigen::VectorXd f = Eigen::VectorXd::Zero(model->NodeNum() * 6);
        for (size_t j = 0; j < model->NodeNum(); j++)
        {
            int p = j * 6;
            v(p) = mode_vectors[i][j].Dx();
            v(p + 1) = mode_vectors[i][j].Dy();
            v(p + 2) = mode_vectors[i][j].Dz();
            v(p + 3) = mode_vectors[i][j].Rx();
            v(p + 4) = mode_vectors[i][j].Ry();
            v(p + 5) = mode_vectors[i][j].Rz();

            f(p) = direction.x;
            f(p + 1) = direction.y;
            f(p + 2) = direction.z;
        }

        Eigen::VectorXd vM = mass_mat.selfadjointView<Eigen::Upper>() * v;
        double beta_i = f.dot(vM) / v.dot(vM);

        participation_factors.push_back(beta_i);
    }

    return participation_factors;
}

std::vector<Displacement> FEVibrateResult::ParticipationDirectedFactors()
{
    // MassMatrixの組み立て
    Eigen::SparseMatrix<double> mass_mat = model->AssembleMassMatrix();
    std::vector<Displacement> participation_factors;

    for (size_t i = 0; i < ModeNum(); i++)
    {
        Eigen::VectorXd v = Eigen::VectorXd::Zero(model->NodeNum() * 6);
        for (size_t j = 0; j < model->NodeNum(); j++)
        {
            int p = j * 6;
            v(p) = mode_vectors[i][j].Dx();
            v(p + 1) = mode_vectors[i][j].Dy();
            v(p + 2) = mode_vectors[i][j].Dz();
            v(p + 3) = mode_vectors[i][j].Rx();
            v(p + 4) = mode_vectors[i][j].Ry();
            v(p + 5) = mode_vectors[i][j].Rz();
        }

        Eigen::VectorXd vM = mass_mat.selfadjointView<Eigen::Upper>() * v;
        vM /= v.dot(vM);
        Displacement facs;
        for (size_t j = 0; j < model->NodeNum(); j++)
        {
            int p = j * 6;
            // 方向別足し合わせ
            Displacement dj(
                vM(p), vM(p + 1), vM(p + 2),
                vM(p + 3), vM(p + 4), vM(p + 5));
            facs += dj;
        }
        participation_factors.push_back(facs);
    }

    return participation_factors;
}

std::vector<double> FEVibrateResult::EffectiveMassRates()
{
    Eigen::SparseMatrix<double> mass_mat = model->AssembleMassMatrix();

    double total_mass = mass_mat.diagonal().sum();

    std::vector<double> mass_rates;

    for (size_t i = 0; i < ModeNum(); i++)
    {
        Eigen::VectorXd v = Eigen::VectorXd::Zero(model->NodeNum() * 6);
        for (size_t j = 0; j < model->NodeNum(); j++)
        {
            int p = j * 6;
            v(p) = mode_vectors[i][j].Dx();
            v(p + 1) = mode_vectors[i][j].Dy();
            v(p + 2) = mode_vectors[i][j].Dz();
            v(p + 3) = mode_vectors[i][j].Rx();
            v(p + 4) = mode_vectors[i][j].Ry();
            v(p + 5) = mode_vectors[i][j].Rz();
        }

        Eigen::VectorXd vM = mass_mat.selfadjointView<Eigen::Upper>() * v;
        double mi = v.dot(vM);
        double beta_i = Eigen::VectorXd::Ones(v.size()).dot(vM) / mi;

        mass_rates.push_back(beta_i * beta_i * mi / total_mass);
    }

    return mass_rates;
}

std::vector<Displacement> FEVibrateResult::EffectiveDirectedMassRates()
{
    // MassMatrixの組み立て
    Eigen::SparseMatrix<double> mass_mat = model->AssembleMassMatrix();
    double total_mass = mass_mat.diagonal().sum();
    std::vector<Displacement> mass_rates;

    for (size_t i = 0; i < ModeNum(); i++)
    {
        Eigen::VectorXd v = Eigen::VectorXd::Zero(model->NodeNum() * 6);
        for (size_t j = 0; j < model->NodeNum(); j++)
        {
            int p = j * 6;
            v(p) = mode_vectors[i][j].Dx();
            v(p + 1) = mode_vectors[i][j].Dy();
            v(p + 2) = mode_vectors[i][j].Dz();
            v(p + 3) = mode_vectors[i][j].Rx();
            v(p + 4) = mode_vectors[i][j].Ry();
            v(p + 5) = mode_vectors[i][j].Rz();
        }

        Eigen::VectorXd vM = mass_mat.selfadjointView<Eigen::Upper>() * v;
        double vMv = v.dot(vM);
        Displacement facs;
        for (size_t j = 0; j < model->NodeNum(); j++)
        {
            int p = j * 6;
            // 方向別足し合わせ
            Displacement dj(
                vM(p), vM(p + 1), vM(p + 2),
                vM(p + 3), vM(p + 4), vM(p + 5));
            facs += dj;
        }

        facs = Displacement(
            facs.Dx() * facs.Dx() / vMv,
            facs.Dy() * facs.Dy() / vMv,
            facs.Dz() * facs.Dz() / vMv,
            facs.Rx() * facs.Rx() / vMv,
            facs.Ry() * facs.Ry() / vMv,
            facs.Rz() * facs.Rz() / vMv);

        facs = facs / (total_mass / 3.0);
        mass_rates.push_back(facs);
    }

    return mass_rates;
}
