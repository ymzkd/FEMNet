#include <unordered_map>
#include <unordered_set>
#include <algorithm>

#define MODEL_PI 3.141592653589793238462643

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include "Elements/Elements.h"
#include "Components.h"
#include "LoadComponent.h"

#include "Model.h"
#include "SparseMatrixUtils.h"

#include "analysis/FEDynamic.h"

FEModel::FEModel()
    : RigidLinkData(std::make_shared<RigidLinks>())
{
}

std::vector<int> FEModel::FreeIndices(bool rigid_link)
{
    std::vector<int> slave_indices = RigidLinkData->SlaveDOFIndices();
    std::unordered_set<int> slaveid_set(slave_indices.begin(), slave_indices.end());
    std::vector<int> indices;
    int idx = 0;
    for (Node n : Nodes)
        for (const bool fixed : n.Fix.isdof_fixed()) {
            
            bool is_slave = (slaveid_set.find(idx) != slaveid_set.end());
            if (rigid_link){
                if (!fixed && !is_slave)
                    indices.push_back(idx);
            }
            else if (!fixed)
                indices.push_back(idx);
            
            idx++;
        }
    return indices;
}

std::vector<int> FEModel::SlaveIndices()
{
    return RigidLinkData->SlaveDOFIndices();
}

/// <summary>
/// 固定された自由度について全体自由度におけるインデックスを格納
/// </summary>
/// <returns>長さが固定自由度数で全体自由度インデックスが格納されたint型vector</returns>
std::vector<int> FEModel::FixIndices()
{
    std::vector<int> indices;
    int idx = 0;
    for (Node n : Nodes) {
        auto fixed = n.Fix.isdof_fixed();
        for (bool f : fixed) {
            if (f) indices.push_back(idx);
            idx++;
        }
    }
    return indices;
}

void FEModel::add_element(BeamElement data)
{
    std::shared_ptr<BeamElement> ptr = std::make_shared<BeamElement>(data);
    Elements.push_back(ptr);
    
    Nodes[data.Nodes[0]->id].Fix.UnlockAllRot();
    Nodes[data.Nodes[1]->id].Fix.UnlockAllRot();
}

void FEModel::add_element(ComplexBeamElement data)
{
    std::shared_ptr<ComplexBeamElement> ptr = std::make_shared<ComplexBeamElement>(data);
    Elements.push_back(ptr);

    Nodes[data.Nodes[0]->id].Fix.UnlockAllRot();
    Nodes[data.Nodes[1]->id].Fix.UnlockAllRot();
}

void FEModel::add_element(TrussElement data)
{
    std::shared_ptr<TrussElement> ptr = std::make_shared<TrussElement>(data);
    Elements.push_back(ptr);
}

void FEModel::add_element(TensionTrussElement data)
{
    std::shared_ptr<TensionTrussElement> ptr = std::make_shared<TensionTrussElement>(data);
    Elements.push_back(ptr);
}

void FEModel::add_element(TriPlaneElement data)
{
    std::shared_ptr<TriPlaneElement> ptr = std::make_shared<TriPlaneElement>(data);
    Elements.push_back(ptr);
}

void FEModel::add_element(TriPlateElement data)
{
    std::shared_ptr<TriPlateElement> ptr = std::make_shared<TriPlateElement>(data);
    Elements.push_back(ptr);

    Nodes[ptr->Nodes[0]->id].Fix.UnlockAllRot();
    Nodes[ptr->Nodes[1]->id].Fix.UnlockAllRot();
    Nodes[ptr->Nodes[2]->id].Fix.UnlockAllRot();
}

void FEModel::add_element(QuadPlaneElement data)
{
    std::shared_ptr<QuadPlaneElement> ptr = std::make_shared<QuadPlaneElement>(data);
    Elements.push_back(ptr);
}

void FEModel::add_element(QuadPlateElement data)
{
    std::shared_ptr<QuadPlateElement> ptr = std::make_shared<QuadPlateElement>(data);
    Elements.push_back(ptr);

    Nodes[ptr->Nodes[0]->id].Fix.UnlockAllRot();
    Nodes[ptr->Nodes[1]->id].Fix.UnlockAllRot();
    Nodes[ptr->Nodes[2]->id].Fix.UnlockAllRot();
    Nodes[ptr->Nodes[3]->id].Fix.UnlockAllRot();
}

void FEModel::add_truss_element(int id, int n1_id, int n2_id, int sec_id, int mat_id)
{
    Node* n1 = &Nodes[n1_id];
    Node* n2 = &Nodes[n2_id];
    Section* sec = &Sections[sec_id];
    Material mat = Materials[mat_id];

    std::shared_ptr<TrussElement> ptr = std::make_shared<TrussElement>(id, n1, n2, sec, mat);
    Elements.push_back(ptr);
}

void FEModel::add_beam_element(int id, int n1_id, int n2_id, int sec_id, int mat_id, double beta)
{
    Node* n1 = &Nodes[n1_id];
    Node* n2 = &Nodes[n2_id];
    Section* sec = &Sections[sec_id];
    Material mat = Materials[mat_id];

    std::shared_ptr<BeamElement> ptr = std::make_shared<BeamElement>(id, n1, n2, sec, mat, beta);
    Elements.push_back(ptr);

    Nodes[ptr->Nodes[0]->id].Fix.UnlockAllRot();
    Nodes[ptr->Nodes[1]->id].Fix.UnlockAllRot();
}

void FEModel::add_tri_plate_element(int id, int n1_id, int n2_id, int n3_id, double thickness, int mat_id)
{
    Node* n1 = &Nodes[n1_id];
    Node* n2 = &Nodes[n2_id];
    Node* n3 = &Nodes[n3_id];
    Material mat = Materials[mat_id];

    std::shared_ptr<TriPlateElement> ptr = std::make_shared<TriPlateElement>(id, n1, n2, n3, thickness, mat);
    Elements.push_back(ptr);

    Nodes[ptr->Nodes[0]->id].Fix.UnlockAllRot();
    Nodes[ptr->Nodes[1]->id].Fix.UnlockAllRot();
    Nodes[ptr->Nodes[2]->id].Fix.UnlockAllRot();

}

void FEModel::add_quad_plate_element(int id, int n1_id, int n2_id, int n3_id, int n4_id, double thickness, int mat_id)
{
    Node* n1 = &Nodes[n1_id];
    Node* n2 = &Nodes[n2_id];
    Node* n3 = &Nodes[n3_id];
    Node* n4 = &Nodes[n4_id];
    Material mat = Materials[mat_id];

    std::shared_ptr<QuadPlateElement> ptr = 
        std::make_shared<QuadPlateElement>(id, n1, n2, n3, n4, thickness, mat);
    Elements.push_back(ptr);

    Nodes[ptr->Nodes[0]->id].Fix.UnlockAllRot();
    Nodes[ptr->Nodes[1]->id].Fix.UnlockAllRot();
    Nodes[ptr->Nodes[2]->id].Fix.UnlockAllRot();
    Nodes[ptr->Nodes[3]->id].Fix.UnlockAllRot();
}

BarElementBase* FEModel::GetBarElement(int id)
{
    BarElementBase* be = dynamic_cast<BarElementBase*>(Elements[id].get());
    return be;
}

BeamElement* FEModel::GetBeamElement(int id)
{
    BeamElement* be = dynamic_cast<BeamElement*>(Elements[id].get());
    return be;
}

TrussElement* FEModel::GetTrussElement(int id)
{
    TrussElement* be = dynamic_cast<TrussElement*>(Elements[id].get());
    return be;
}

QuadPlateElement* FEModel::GetQuadPlateElement(int id)
{
    QuadPlateElement* be = dynamic_cast<QuadPlateElement*>(Elements[id].get());
    return be;
}

TriPlateElement* FEModel::GetTriPlateElement(int id)
{
    TriPlateElement* be = dynamic_cast<TriPlateElement*>(Elements[id].get());
    return be;
}

std::vector<NodeLoadData> FEModel::InnertialForceToNodeLoads(const InertialForce inertial_force)
{
    std::vector<NodeLoadData> node_loads;
    //if (std::shared_ptr<InertialForce> inertial = std::dynamic_pointer_cast<InertialForce>(load)) {
    for (auto& e : Elements)
    {
        std::vector<NodeLoadData> elem_loads =
            e->InertialForceToNodeLoadData(
                Eigen::Vector3d(
                    inertial_force.accels.x, 
                    inertial_force.accels.y, 
                    inertial_force.accels.z));
        node_loads.insert(node_loads.end(), elem_loads.begin(), elem_loads.end());
    }
    return node_loads;
}

double FEModel::SumNodeMass()
{
    Vector mass;
    double sum = 0;
    for (Node n : Nodes)
        sum += n.MassData.SumMass();

    return sum;
}

Eigen::SparseMatrix<double> FEModel::AssembleStiffnessMatrix(const ElementStates *states)
{
    int mat_size = Nodes.size() * 6;

    // 非零要素数を推定
    size_t estimated_nnz = Elements.size() * 100;
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(estimated_nnz);

    // 各要素からTripletを収集
    for (size_t i = 0; i < Elements.size(); i++) {
        const std::shared_ptr<ElementBase>& eh = Elements[i];
        if (states) {
            // 状態依存要素は指定状態に応じた接線剛性を組立てる
            if (auto sde = std::dynamic_pointer_cast<IStateDependentElement>(eh)) {
                sde->GetTangentStiffnessTriplets(tripletList, states->Get((int)i));
                continue;
            }
        }
        eh->GetStiffnessTriplets(tripletList);
    }

    // Tripletから疎行列を一括構築
    Eigen::SparseMatrix<double> mat(mat_size, mat_size);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());

    return mat;
}

Eigen::SparseMatrix<double> FEModel::AssembleMassMatrix()
{
    std::vector<Eigen::Triplet<double>> tripletList;
    for (int i = 0; i < NodeNum(); ++i) {
        tripletList.push_back(Eigen::Triplet<double>(i * 6, i * 6, 0.0));
        tripletList.push_back(Eigen::Triplet<double>(i * 6 + 1, i * 6 + 1, 0.0));
        tripletList.push_back(Eigen::Triplet<double>(i * 6 + 2, i * 6 + 2, 0.0));
    }
    Eigen::SparseMatrix<double> mass_mat(NodeNum() * 6, NodeNum() * 6);
    mass_mat.setFromTriplets(tripletList.begin(), tripletList.end());

    // Construct MassMatrix
    //for each(std::shared_ptr<ElementBase> eh in Elements)
    //    eh->AssembleMassMatrix(mass_mat);

    for (const Node& n : Nodes)
    {
        int i = n.id * 6;
        mass_mat.coeffRef(i, i) += n.MassData.SumMass();
        mass_mat.coeffRef(i + 1, i + 1) += n.MassData.SumMass();
        mass_mat.coeffRef(i + 2, i + 2) += n.MassData.SumMass();
    }

    mass_mat *= (1.0 / GraityAccel);
    return mass_mat;
}

Eigen::SparseMatrix<double> FEModel::AssembleGeometricStiffnessMatrix(
    const std::vector<Displacement>& displacements)
{
    int mat_size = Nodes.size() * 6;

    size_t estimated_nnz = Elements.size() * 100;
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(estimated_nnz);

    // 各要素の節点変位を準備して Triplet を収集
    for (const std::shared_ptr<ElementBase>& eh : Elements) {
        std::vector<Displacement> disp_vec(eh->NodeNum());
        std::vector<Node*> nodes = eh->NodesList();
        for (size_t i = 0; i < eh->NodeNum(); i++)
            disp_vec[i] = displacements[nodes[i]->id];

        eh->GetGeometricStiffnessTriplets(disp_vec, tripletList);
    }

    Eigen::SparseMatrix<double> mat(mat_size, mat_size);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());

    return mat;
}

void FEModel::ComputeElementNodeMass()
{
    for (size_t i = 0; i < Nodes.size(); i++)
        Nodes[i].MassData.ElementMass = 0.0;

    for (std::shared_ptr<ElementBase> eh : Elements) {
        Eigen::VectorXd masses = eh->NodeLumpedMass();
        for (size_t i = 0; i < eh->NodeNum(); i++) {
            int node_id = eh->NodesList()[i]->id;
            Nodes[node_id].MassData.ElementMass += masses(i);
        }
    }
}

Eigen::VectorXd FEModel::AssembleLoadVector(const std::vector<std::shared_ptr<LoadBase>> &loads)
{
    Eigen::VectorXd f(Nodes.size() * 6);
    f.setZero();
    for (auto& load : loads)
    {

        std::vector<NodeLoadData> node_loads;
        if (std::shared_ptr<InertialForce> inertial = std::dynamic_pointer_cast<InertialForce>(load)) {
            for (auto& e : Elements)
            {
                std::vector<NodeLoadData> elem_loads =
                    e->InertialForceToNodeLoadData(Eigen::Vector3d(inertial->accels.x, inertial->accels.y, inertial->accels.z));
                node_loads.insert(node_loads.end(), elem_loads.begin(), elem_loads.end());
            }
        }
        else {
            node_loads = load->NodeLoads();
        }

        for (auto &nl : node_loads)
        {
            if (nl.id < 0) continue;
            int pos = nl.id * 6;
            f[pos] += nl.Px();
            f[pos + 1] += nl.Py();
            f[pos + 2] += nl.Pz();
            f[pos + 3] += nl.Mx();
            f[pos + 4] += nl.My();
            f[pos + 5] += nl.Mz();
        }
    }
    return f;
}

double IResponseSpectrum::DampingCorrectionFactor(const double t)
{
    if (DampInitializer == nullptr) return 1.0;

    double h = DampInitializer->DampRateAtPeriod(t);
    // 算定不能(t<=0のZPAや初期化前など、h<0)は補正なし。
    // これにより零周期応答(t=0)は減衰非依存となり、Fhが負に化けるのも防ぐ。
    if (h < 0.0) return 1.0;
    return 1.5 / (1.0 + 10.0 * h);
}

double IResponseSpectrum::acceleration_factored(double t)
{
    if (enable_damp_factor)
        return DampingCorrectionFactor(t) * Acceleration(t);
    else
        return Acceleration(t);
}

double IResponseSpectrum::velocity_factored(double t)
{
    if (enable_damp_factor)
        return DampingCorrectionFactor(t) * Velocity(t);
    else
        return Velocity(t);
}

double IResponseSpectrum::displacement_factored(double t)
{
    if (enable_damp_factor)
        return DampingCorrectionFactor(t) * Displacement(t);
    else
        return Displacement(t);
}
