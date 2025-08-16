#include "ResponseSpectrumMethod.h"

std::vector<Displacement> ResponseSpectrumMethod::calculate_responseCQC(ResponseValueType vt)
{
    std::vector<double> part_facs = VibrateResult.ParticipationFactors(Direction);
    std::vector<double> periods = VibrateResult.NaturalPeriods();
    std::vector<double> spectrums(periods.size());
    for (size_t i = 0; i < periods.size(); i++)
    {
        if (vt == ResponseValueType::Displacement)
            spectrums[i] = SpectrumFunction->Displacement(periods[i]);
        else if (vt == ResponseValueType::Velocity)
            spectrums[i] = SpectrumFunction->Velocity(periods[i]);
        else // (vt == ResponseValueType::Acceleration)
            spectrums[i] = SpectrumFunction->Acceleration(periods[i]);
    }

    std::vector<Displacement> responses(model->NodeNum());
    for (size_t j = 0; j < part_facs.size(); j++)
    {
        std::vector<Displacement> uj = VibrateResult.ModeVectors()[j];
        for (size_t k = 0; k < part_facs.size(); k++)
        {
            std::vector<Displacement> uk = VibrateResult.ModeVectors()[k];

            double rjk = periods[k] / periods[j];
            double correlation = 8.0 * damping_rate * damping_rate * (1.0 + rjk) * pow(rjk, 1.5) /
                                 (pow((1.0 - rjk * rjk), 2.0) + 4.0 * damping_rate * damping_rate * rjk * pow(1.0 + rjk, 2.0));
            double fac = spectrums[j] * spectrums[k] * part_facs[j] * part_facs[k] * correlation;

            for (size_t i = 0; i < model->NodeNum(); i++)
            {
                Displacement idj = uj[i];
                Displacement idk = uk[i];

                Displacement d = Displacement(
                    idj.Dx() * idk.Dx(), idj.Dy() * idk.Dy(), idj.Dz() * idk.Dz(),
                    idj.Rx() * idk.Rx(), idj.Ry() * idk.Ry(), idj.Rz() * idk.Rz());

                responses[i] += fac * d;
            }
        }
    }

    // 平方根を取る
    for (size_t j = 0; j < responses.size(); j++)
    {
        responses[j] = Displacement(
            sqrt(responses[j].Dx()), sqrt(responses[j].Dy()), sqrt(responses[j].Dz()),
            sqrt(responses[j].Rx()), sqrt(responses[j].Ry()), sqrt(responses[j].Rz()));
    }

    return responses;
}

std::vector<Displacement> ResponseSpectrumMethod::calculate_responseSRSS(ResponseValueType vt)
{
    std::vector<double> part_facs = VibrateResult.ParticipationFactors(Direction);
    std::vector<double> periods = VibrateResult.NaturalPeriods();
    std::vector<double> spectrums(periods.size());
    for (size_t i = 0; i < periods.size(); i++)
    {
        if (vt == ResponseValueType::Displacement)
            spectrums[i] = SpectrumFunction->Displacement(periods[i]);
        else if (vt == ResponseValueType::Velocity)
            spectrums[i] = SpectrumFunction->Velocity(periods[i]);
        else // (vt == ResponseValueType::Acceleration)
            spectrums[i] = SpectrumFunction->Acceleration(periods[i]);
    }

    std::vector<Displacement> responses(model->NodeNum());
    for (size_t i = 0; i < part_facs.size(); i++)
    {
        std::vector<Displacement> mode_vector = VibrateResult.ModeVectors()[i];

        // SRSS Method
        for (size_t j = 0; j < mode_vector.size(); j++)
        {
            Displacement d2 = spectrums[i] * part_facs[i] * mode_vector[j];
            // 変位の二乗和を計算
            d2 = Displacement(d2.Dx() * d2.Dx(), d2.Dy() * d2.Dy(), d2.Dz() * d2.Dz(),
                              d2.Rx() * d2.Rx(), d2.Ry() * d2.Ry(), d2.Rz() * d2.Rz());
            responses[j] += d2;
        }
    }

    // SRSS Methodの結果を平方根で正規化
    for (size_t j = 0; j < responses.size(); j++)
    {
        responses[j] = Displacement(
            sqrt(responses[j].Dx()), sqrt(responses[j].Dy()), sqrt(responses[j].Dz()),
            sqrt(responses[j].Rx()), sqrt(responses[j].Ry()), sqrt(responses[j].Rz()));
    }

    return responses;
}

std::vector<Displacement> ResponseSpectrumMethod::calculate_responseABS(ResponseValueType vt)
{
    std::vector<double> part_facs = VibrateResult.ParticipationFactors(Direction);
    std::vector<double> periods = VibrateResult.NaturalPeriods();
    std::vector<double> spectrums(periods.size());
    for (size_t i = 0; i < periods.size(); i++)
    {
        if (vt == ResponseValueType::Displacement)
            spectrums[i] = SpectrumFunction->Displacement(periods[i]);
        else if (vt == ResponseValueType::Velocity)
            spectrums[i] = SpectrumFunction->Velocity(periods[i]);
        else // (vt == ResponseValueType::Acceleration)
            spectrums[i] = SpectrumFunction->Acceleration(periods[i]);
    }

    std::vector<Displacement> responses(model->NodeNum());
    for (size_t i = 0; i < part_facs.size(); i++)
    {
        std::vector<Displacement> mode_vector = VibrateResult.ModeVectors()[i];
        // | sd x vector x beta_i |

        // ABS Method
        for (size_t j = 0; j < mode_vector.size(); j++)
        {
            // 変位の絶対値を計算
            Displacement d = spectrums[i] * part_facs[i] * mode_vector[j];
            responses[j] += Displacement(abs(d.Dx()), abs(d.Dy()), abs(d.Dz()),
                                         abs(d.Rx()), abs(d.Ry()), abs(d.Rz()));
        }
    }

    return responses;
}

std::vector<Displacement> ResponseSpectrumMethod::calculate_response(ResponseValueType vt)
{
    std::vector<Displacement> responses(model->NodeNum());
    if (MethodType == ResponseSpectrumMethodType::CQC)
        responses = calculate_responseCQC(vt);
    else if (MethodType == ResponseSpectrumMethodType::SRSS)
        responses = calculate_responseSRSS(vt);
    else // MethodType == ResponseSpectrumMethodType::ABS
        responses = calculate_responseABS(vt);

    return responses;
}

ResponseSpectrumMethod::ResponseSpectrumMethod(std::shared_ptr<FEModel> model,
                                               FEVibrateResult vibrate_result, Vector direction, IResponseSpectrum *spectrum_function, ResponseSpectrumMethodType type)
    : FEDeformOperator(model), VibrateResult(vibrate_result), SpectrumFunction(spectrum_function), Direction(direction), MethodType(type)
{

    Compute();
}

void ResponseSpectrumMethod::Compute()
{
    displacements = calculate_response(ResponseValueType::Displacement);
    velocities = calculate_response(ResponseValueType::Velocity);
    accelerations = calculate_response(ResponseValueType::Acceleration);
    m_computed = true;
}

std::vector<Displacement> ResponseSpectrumMethod::GetDisplacements()
{
    if (m_computed)
        return displacements;
    else
        return calculate_response(ResponseValueType::Displacement);
}

std::vector<Displacement> ResponseSpectrumMethod::GetVelocities()
{
    if (m_computed)
        return velocities;
    else
        return calculate_response(ResponseValueType::Velocity);
}

std::vector<Displacement> ResponseSpectrumMethod::GetAccelerations()
{
    if (m_computed)
        return accelerations;
    else
        return calculate_response(ResponseValueType::Acceleration);
}

BeamStressData ResponseSpectrumMethod::GetBeamStress(int eid, double p)
{
    if (!m_computed)
        throw std::runtime_error("ResponseSpectrumMethod: need to call Compute()");

    BarElementBase *be = dynamic_cast<BarElementBase *>(model->Elements[eid].get());

    BeamStress b_strs = be->stress(this->GetDisplacements()[be->Nodes[0]->id], this->GetDisplacements()[be->Nodes[1]->id]);
    BeamStressData strs = b_strs.Interpolate(p);
    return strs;
}

PlateStressData ResponseSpectrumMethod::GetPlateStressData(int eid, double xi, double eta)
{
    if (!m_computed)
        throw std::runtime_error("ResponseSpectrumMethod: need to call Compute()");

    PlateStressData data;
    std::vector<Displacement> displace = this->GetDisplacements();
    if (model->Elements[eid]->Type() == ElementType::DKT)
    {
        TriPlateElement *el = dynamic_cast<TriPlateElement *>(model->Elements[eid].get());
        data = el->stress(displace[el->Nodes[0]->id], displace[el->Nodes[1]->id],
                          displace[el->Nodes[2]->id], xi, eta);
    }
    else if (model->Elements[eid]->Type() == ElementType::DKQ)
    {
        QuadPlateElement *el = dynamic_cast<QuadPlateElement *>(model->Elements[eid].get());
        data = el->stress(displace[el->Nodes[0]->id], displace[el->Nodes[1]->id],
                          displace[el->Nodes[2]->id], displace[el->Nodes[3]->id], xi, eta);
    }
    return data;
}

Displacement ResponseSpectrumMethod::GetBeamDisplace(int eid, double p)
{
    if (!m_computed)
        throw std::runtime_error("ResponseSpectrumMethod: need to call Compute()");

    BeamElement *elm = model->GetBeamElement(eid);
    std::vector<Displacement> displace = this->GetDisplacements();
    Displacement disp = elm->DisplaceAt(
        displace[elm->Nodes[0]->id], displace[elm->Nodes[1]->id], p);

    return disp;
}

FELinearStaticOp ResponseSpectrumMethod::GetLinearStaticCase()
{
    std::vector<std::shared_ptr<LoadBase>> loads;
    std::vector<Displacement> accels = this->GetAccelerations();
    for (size_t i = 0; i < model->NodeNum(); i++)
    {
        // とりあえず並進だけ
        double x = accels[i].Dx() / model->GraityAccel;
        double y = accels[i].Dy() / model->GraityAccel;
        double z = accels[i].Dz() / model->GraityAccel;
        loads.push_back(std::make_shared<NodeBodyForce>(NodeBodyForce(&model->Nodes[i], x, y, z)));
    }
    return FELinearStaticOp(model, loads);
}
