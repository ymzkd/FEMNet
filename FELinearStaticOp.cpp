#include "FELinearStaticOp.h"

void FELinearStaticOp::Compute()
{
    this->model->SolveLinearStatic(this->loads, this->displace, this->react_force);
    m_computed = true;
}

BeamStressData FELinearStaticOp::GetBeamStress(int eid, double p)
{
    if (!m_computed)
        throw std::runtime_error("FELinearStaticOP: need to call compute()");

    BarElementBase *be = dynamic_cast<BarElementBase *>(model->Elements[eid].get());

    // BeamElement* elm = GetBeamElement(eid);
    BeamStress b_strs = be->stress(displace[be->Nodes[0]->id], displace[be->Nodes[1]->id]);
    BeamStressData strs = b_strs.Interpolate(p);

    // Add beam stresses due to beam loads
    BeamElement *beamElement = dynamic_cast<BeamElement *>(be);
    if (beamElement != nullptr)
    {
        for (const auto &l : loads)
        {
            // BeamLoadBase* bpl = dynamic_cast<BeamLoadBase*>(l);
            std::shared_ptr<BeamLoadBase> bpl = std::dynamic_pointer_cast<BeamLoadBase>(l);
            if (bpl == NULL)
                continue;
            if (bpl->element->id != beamElement->id)
                continue;
            BeamStressData bsd = bpl->GetBeamStress(p);
            strs.Nx += bsd.Nx;
            strs.My += bsd.My;
            strs.Mz += bsd.Mz;
            strs.Mx += bsd.Mx;
            strs.Qy += bsd.Qy;
            strs.Qz += bsd.Qz;
        }
    }
    return strs;
}

PlateStressData FELinearStaticOp::GetPlateStressData(int eid, double xi, double eta)
{
    if (!m_computed)
        throw std::runtime_error("FELinearStaticOP: need to call compute()");

    PlateStressData data;
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

Displacement FELinearStaticOp::GetBeamDisplace(int eid, double p)
{
    if (!m_computed)
        throw std::runtime_error("FELinearStaticOP: need to call compute()");

    BeamElement *elm = model->GetBeamElement(eid);
    Displacement disp = elm->DisplaceAt(
        displace[elm->Nodes[0]->id], displace[elm->Nodes[1]->id], p);

    // BeamStressData strs = b_strs.Interpolate(p);
    for (const auto &l : loads)
    {
        std::shared_ptr<BeamLoadBase> bpl = std::dynamic_pointer_cast<BeamLoadBase>(l);
        if (bpl == NULL)
            continue;
        if (bpl->element->id != elm->id)
            continue;

        // Not Implemented
        if (bpl->axis == BeamLoadAxis::XAxis)
            continue;

        Displacement disp_i = bpl->GetDisplacement(p);
        disp = Displacement(
            disp.Dx() + disp_i.Dx(),
            disp.Dy() + disp_i.Dy(),
            disp.Dz() + disp_i.Dz(),
            disp.Rx() + disp_i.Rx(),
            disp.Ry() + disp_i.Ry(),
            disp.Rz() + disp_i.Rz());
    }
    return disp;
}

std::vector<Displacement> FELinearStaticOp::GetDisplacements()
{
    if (!m_computed)
        throw std::runtime_error("FELinearStaticOP: need to call compute()");

    return this->displace;
}

std::vector<NodeLoad> FELinearStaticOp::GetReactForces()
{
    if (!m_computed)
        throw std::runtime_error("FELinearStaticOP: need to call compute()");

    return this->react_force;
}

BeamStressData LinearStaticCombinationOperator::GetBeamStress(int eid, double p)
{
    BeamStressData data;
    for (LinearStaticDeformFactor op : this->cases)
    {
        if (!op.op->Computed())
            throw std::runtime_error("FELinearStaticOP: need to call compute()");

        data += op.factor * op.op->GetBeamStress(eid, p);
    }
    return data;
}

PlateStressData LinearStaticCombinationOperator::GetPlateStressData(int eid, double xi, double eta)
{
    PlateStressData data;
    for (LinearStaticDeformFactor op : this->cases)
    {
        if (!op.op->Computed())
            throw std::runtime_error("FELinearStaticOP: need to call compute()");

        data += op.factor * op.op->GetPlateStressData(eid, xi, eta);
    }

    return data;
}

Displacement LinearStaticCombinationOperator::GetBeamDisplace(int eid, double p)
{
    Displacement disp;
    for (LinearStaticDeformFactor op : this->cases)
    {
        if (!op.op->Computed())
            throw std::runtime_error("FELinearStaticOP: need to call compute()");

        disp += op.factor * op.op->GetBeamDisplace(eid, p);
    }
    return disp;
}

std::vector<Displacement> LinearStaticCombinationOperator::GetDisplacements()
{
    std::vector<Displacement> disp(model->NodeNum());
    for (LinearStaticDeformFactor op : this->cases)
    {
        if (!op.op->Computed())
            throw std::runtime_error("FELinearStaticOP: need to call compute()");

        std::vector<Displacement> d = op.op->GetDisplacements();
        for (size_t i = 0; i < d.size(); i++)
            disp[i] += d[i] * op.factor;
    }
    return disp;
}

std::vector<NodeLoad> LinearStaticCombinationOperator::GetReactForces()
{
    std::vector<NodeLoad> combined_react;
    for (const auto &case_factor : cases)
    {
        if (!case_factor.op->Computed())
            throw std::runtime_error("FELinearStaticOP: need to call compute()");

        auto react = case_factor.op->GetReactForces();
        for (size_t i = 0; i < react.size(); i++)
        {
            bool need_fallback = false;
            // compine_reactがi番目要素を含む場合はi番目要素をチェック
            if (i < combined_react.size())
            {
                // i番目要素が対象のノードか？
                if (combined_react[i].id == react[i].id)
                {
                    combined_react[i].data += (case_factor.factor * react[i].data);
                    continue;
                }
            }

            for (auto &cr : combined_react)
            {
                if (cr.id == react[i].id)
                {
                    cr.data += (case_factor.factor * react[i].data);
                    break;
                }
                need_fallback = true;
            }

            if (need_fallback)
                combined_react.push_back(react[i]);
        }
    }
    return combined_react;
}