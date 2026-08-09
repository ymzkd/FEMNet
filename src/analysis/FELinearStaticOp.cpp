#include "FELinearStaticOp.h"
#include "ReducedSystem.h"

// 状態依存要素(引張専用トラス等)を考慮した反復解法による線形静的解析。
// 要素状態はOperator(m_states)が所有し、剛性組立時にFEModelへ明示的に渡す。
// 収束時は Iterations に反復回数が入り Converged=true、
// 最大反復数に達した場合は Converged=false となる。
void FELinearStaticOp::Compute()
{
    const int max_iter = 100;
    FEModel &m = *model;
    const int full_size = (int)m.Nodes.size() * 6;

    Eigen::VectorXd force_vec = m.AssembleLoadVector(loads);
    Eigen::VectorXd residual_vec = force_vec;
    Eigen::VectorXd disp_vec = Eigen::VectorXd::Zero(full_size);

    ReducedSystem rs(m);

    // 状態依存要素の状態を初期化 (規定剛性で機能する状態へ)
    m_states.Clear();

    Eigen::SparseMatrix<double> full_stiffmat = m.AssembleStiffnessMatrix(&m_states);

    int iter_result = -max_iter;
    for (int iter = 0; iter < max_iter; iter++)
    {
        // 縮約系の構築
        Eigen::SparseMatrix<double> mii, mij;
        rs.Reduce(full_stiffmat, mii, &mij);

        Eigen::VectorXd f_input, f_fix;
        rs.ReduceVector(residual_vec, f_input, f_fix);

        // Solve
        auto solver_static = createSolver();
        solver_static->compute(mii);
        Eigen::VectorXd d_result = solver_static->solve(f_input);
        Eigen::VectorXd r_fix = mij.transpose() * d_result - f_fix;

        // 反力データ整理
        Eigen::VectorXd r = Eigen::VectorXd::Zero(full_size);
        for (size_t i = 0; i < rs.fixed_indices.size(); i++)
            r(rs.fixed_indices[i]) = r_fix(i);
        react_force.clear();
        for (size_t i = 0; i < m.Nodes.size(); i++)
        {
            if (!m.Nodes[i].Fix.IsAnyFix())
                continue;
            int pos = i * 6;
            react_force.push_back(NodeLoad(i, r[pos], r[pos + 1], r[pos + 2], r[pos + 3], r[pos + 4], r[pos + 5]));
        }

        // 変形データ整理
        disp_vec += rs.ExpandVector(d_result, full_size);
        displace.clear();
        for (size_t i = 0; i < m.Nodes.size(); i++)
        {
            int pos = i * 6;
            displace.push_back(Displacement(disp_vec[pos], disp_vec[pos + 1], disp_vec[pos + 2], disp_vec[pos + 3], disp_vec[pos + 4], disp_vec[pos + 5]));
        }

        // 判定と更新
        // 状態依存要素の次状態を判定し、状態変化があれば剛性を再構築する
        bool any_change = false;
        for (size_t i = 0; i < m.Elements.size(); i++)
        {
            if (auto sde = std::dynamic_pointer_cast<IStateDependentElement>(m.Elements[i]))
            {
                bool current = m_states.Get((int)i);
                bool next = sde->NextState(displace, current);
                if (next != current)
                {
                    m_states.Set((int)i, next);
                    any_change = true;
                }
            }
        }

        if (!any_change)
        {
            iter_result = iter + 1; // 収束判定: 状態変化なしなら終了
            break;
        }
        full_stiffmat = m.AssembleStiffnessMatrix(&m_states); // 状態変化あり: 剛性行列を再構築
        residual_vec = force_vec - full_stiffmat.selfadjointView<Eigen::Upper>() * disp_vec; // 内力を再計算
    }

    Converged = iter_result > 0;
    Iterations = iter_result > 0 ? iter_result : -iter_result;
    m_computed = true;
}

std::shared_ptr<FELinearStaticOp> FELinearStaticOp::FromCombination(
    std::shared_ptr<FEModel> model,
    const std::vector<LinearStaticDeformFactor>& cases)
{
    std::vector<std::shared_ptr<LoadBase>> merged;
    for (const auto& c : cases)
    {
        if (c.op == nullptr)
            throw std::invalid_argument("FromCombination: case operator is null");
        for (const auto& load : c.op->loads)
            merged.push_back(load->scaled(c.factor));
    }
    return std::make_shared<FELinearStaticOp>(model, merged);
}

BeamStressData FELinearStaticOp::GetBeamStress(int eid, double p)
{
    if (!m_computed)
        throw std::runtime_error("FELinearStaticOP: need to call compute()");

    BarElementBase *be = dynamic_cast<BarElementBase *>(model->Elements[eid].get());

    if (auto sde = std::dynamic_pointer_cast<IStateDependentElement>(model->Elements[eid]))
    {
        // Operatorが所有するこのケースの収束時状態で応力を復元する
        BeamStress b_strs = sde->tangent_stress(
            displace[be->Nodes[0]->id], displace[be->Nodes[1]->id], m_states.Get(eid));
        BeamStressData strs = b_strs.Interpolate(p);
        return strs;
    }

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

std::vector<NodeLoadData> FELinearStaticOp::GetPlateNodalForces(int eid, bool local)
{
    if (!m_computed)
        throw std::runtime_error("FELinearStaticOP: need to call compute()");

    if (model->Elements[eid]->Type() == ElementType::DKT)
    {
        TriPlateElement *el = dynamic_cast<TriPlateElement *>(model->Elements[eid].get());
        return el->NodalForces(displace[el->Nodes[0]->id], displace[el->Nodes[1]->id],
                               displace[el->Nodes[2]->id], local);
    }
    else if (model->Elements[eid]->Type() == ElementType::DKQ)
    {
        QuadPlateElement *el = dynamic_cast<QuadPlateElement *>(model->Elements[eid].get());
        return el->NodalForces(displace[el->Nodes[0]->id], displace[el->Nodes[1]->id],
                               displace[el->Nodes[2]->id], displace[el->Nodes[3]->id], local);
    }
    return std::vector<NodeLoadData>();
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

std::vector<BeamStressData> LinearStaticCombinationOperator::GetBeamStressComponents(int eid, double p)
{
    std::vector<BeamStressData> results;
    for (auto &c : cases)
    {
        BeamStressData s = c.op->GetBeamStress(eid, p);
        s *= c.factor; // 係数を適用
        results.push_back(s);
    }
    return results;
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

std::vector<PlateStressData> LinearStaticCombinationOperator::GetPlateStressDataComponents(int eid, double xi, double eta)
{
    std::vector<PlateStressData> results;
    for (auto &c : cases)
    {
        PlateStressData s = c.op->GetPlateStressData(eid, xi, eta);
        s *= c.factor; // 演算子オーバーロード
        results.push_back(s);
    }
    return results;
}

std::vector<NodeLoadData> LinearStaticCombinationOperator::GetPlateNodalForces(int eid, bool local)
{
    // 合成節点変位 u_combined から f = K_e * u_combined を計算（f は u に線形なので厳密）
    std::vector<Displacement> disp = GetDisplacements();

    if (model->Elements[eid]->Type() == ElementType::DKT)
    {
        TriPlateElement *el = dynamic_cast<TriPlateElement *>(model->Elements[eid].get());
        return el->NodalForces(disp[el->Nodes[0]->id], disp[el->Nodes[1]->id],
                               disp[el->Nodes[2]->id], local);
    }
    else if (model->Elements[eid]->Type() == ElementType::DKQ)
    {
        QuadPlateElement *el = dynamic_cast<QuadPlateElement *>(model->Elements[eid].get());
        return el->NodalForces(disp[el->Nodes[0]->id], disp[el->Nodes[1]->id],
                               disp[el->Nodes[2]->id], disp[el->Nodes[3]->id], local);
    }
    return std::vector<NodeLoadData>();
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
            // インデックスが一致する場合の高速パス
            if (i < combined_react.size() && combined_react[i].id == react[i].id)
            {
                combined_react[i].data += (case_factor.factor * react[i].data);
                continue;
            }

            // 全検索で既存ノードを探す
            bool found = false;
            for (auto &cr : combined_react)
            {
                if (cr.id == react[i].id)
                {
                    cr.data += (case_factor.factor * react[i].data);
                    found = true;
                    break;
                }
            }

            // 見つからなければ新規追加（係数を適用）
            if (!found)
            {
                NodeLoad nl = react[i];
                nl.data = case_factor.factor * react[i].data;
                combined_react.push_back(nl);
            }
        }
    }
    return combined_react;
}