// 荷重(LoadBase 派生)テキスト形式 SaveLoads/LoadLoads の往復テスト
//
// 全6種(NodeLoad/InertialForce/NodeBodyForce/PlateLoad/BeamPolyLoad/
// AxialPolyLoad)を含む荷重リストを構築し、
//   SaveLoads -> LoadLoads(model) -> SaveLoads の2ファイルがバイト一致すること、
//   主要フィールドが一致すること、
//   復元した荷重で FELinearStaticOp が解けること
// を検証する。

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <memory>
#include <vector>

#include "Model.h"
#include "Elements/Elements.h"
#include "LoadComponent.h"
#include "LoadIO.h"
#include "FELinearStaticOp.h"

namespace {

int g_failures = 0;

void Check(bool cond, const std::string &msg)
{
    if (!cond)
    {
        std::cout << "  [FAIL] " << msg << "\n";
        g_failures++;
    }
}

std::string ReadFile(const std::string &path)
{
    std::ifstream ifs(path, std::ios::binary);
    std::ostringstream ss;
    ss << ifs.rdbuf();
    return ss.str();
}

// 梁2本(0-1, 1-2)の片持ち + 四角形板(0,1,4,3) を持つ可解モデル
std::shared_ptr<FEModel> BuildModel()
{
    auto m = std::make_shared<FEModel>();

    // 全節点が要素で使用される(未使用の自由節点は特異行列の原因になる)
    double coords[5][3] = {
        {0, 0, 0}, {1000, 0, 0}, {2000, 0, 0},
        {0, 1000, 0}, {1000, 1000, 0}};
    for (int i = 0; i < 5; i++)
        m->Nodes.push_back(Node(i, coords[i][0], coords[i][1], coords[i][2]));

    // 支点: 板側・梁根元を固定
    for (int id : {0, 3, 4})
        for (int i = 0; i < 6; i++)
            m->Nodes[id].Fix.flags[i] = true;
    // 節点2に質量(NodeBodyForce が非ゼロになるよう)
    m->Nodes[2].MassData.Mass = 100.0;

    Material mat(205000.0, 0.3, 7.85e-9);
    mat.id = 0;
    m->Materials.push_back(mat);

    Section s;
    s.id = 0;
    s.A = 1000.0; s.Iy = 2.0e6; s.Iz = 3.0e6; s.Iyz = 0.0; s.K = 4.0e6;
    m->Sections.push_back(s);

    // 要素 (Nodes/Sections は以後 resize しないのでポインタ安定)
    m->Elements.push_back(std::make_shared<BeamElement>(
        0, &m->Nodes[0], &m->Nodes[1], &m->Sections[0], m->Materials[0], 0.0));
    m->Elements.push_back(std::make_shared<BeamElement>(
        1, &m->Nodes[1], &m->Nodes[2], &m->Sections[0], m->Materials[0], 0.0));
    m->Elements.push_back(std::make_shared<QuadPlateElement>(
        2, &m->Nodes[0], &m->Nodes[1], &m->Nodes[4], &m->Nodes[3], 12.0, m->Materials[0], 0.0));

    return m;
}

std::vector<std::shared_ptr<LoadBase>> BuildLoads(FEModel &m)
{
    std::vector<std::shared_ptr<LoadBase>> loads;

    // NodeLoad (節点2)
    loads.push_back(std::make_shared<NodeLoad>(2, 10.0, -1000.0, 0.0, 0.0, 0.0, 5.0));

    // InertialForce
    loads.push_back(std::make_shared<InertialForce>(0.0, -9806.65, 0.0));

    // NodeBodyForce (節点2, 非既定 selector)
    auto nbf = std::make_shared<NodeBodyForce>(&m.Nodes[2], 0.0, -1.0, 0.0);
    nbf->selector = static_cast<BodyForaceSelector>(
        static_cast<unsigned int>(BodyForaceSelector::NodeMass) |
        static_cast<unsigned int>(BodyForaceSelector::LoadMass)); // = 3
    loads.push_back(nbf);

    // PlateLoad (要素2 = 四角形板)
    PlaneElementBase *plate = dynamic_cast<PlaneElementBase *>(m.Elements[2].get());
    loads.push_back(std::make_shared<PlateLoad>(plate, 0.0, 0.0, -0.01));

    // BeamPolyLoad (要素0, YAxis, 等分布)
    loads.push_back(std::make_shared<BeamPolyLoad>(
        std::vector<double>{-5.0, -5.0}, std::vector<double>{0.0, 1.0},
        dynamic_cast<BeamElement *>(m.Elements[0].get()), BeamLoadAxis::YAxis));

    // AxialPolyLoad (要素1)
    loads.push_back(std::make_shared<AxialPolyLoad>(
        std::vector<double>{2.0, 2.0}, std::vector<double>{0.0, 1.0},
        dynamic_cast<BeamElement *>(m.Elements[1].get())));

    return loads;
}

} // namespace

int main()
{
    std::cout << "LoadIO round-trip test\n";

    const std::string path_a = "loadio_a.txt";
    const std::string path_b = "loadio_b.txt";

    std::shared_ptr<FEModel> model = BuildModel();
    std::vector<std::shared_ptr<LoadBase>> loads = BuildLoads(*model);

    SaveLoads(path_a, loads);
    std::vector<std::shared_ptr<LoadBase>> loads2 = LoadLoads(path_a, *model);
    SaveLoads(path_b, loads2);

    // 1) バイト一致
    std::string a = ReadFile(path_a);
    std::string b = ReadFile(path_b);
    Check(!a.empty(), "出力ファイルが空でない");
    Check(a == b, "SaveLoads->LoadLoads->SaveLoads がバイト一致する");

    // 2) 件数・種別
    Check(loads2.size() == 6, "荷重数 = 6");

    Check(std::dynamic_pointer_cast<NodeLoad>(loads2[0]) != nullptr, "load0 = NodeLoad");
    Check(std::dynamic_pointer_cast<InertialForce>(loads2[1]) != nullptr, "load1 = InertialForce");
    Check(std::dynamic_pointer_cast<NodeBodyForce>(loads2[2]) != nullptr, "load2 = NodeBodyForce");
    Check(std::dynamic_pointer_cast<PlateLoad>(loads2[3]) != nullptr, "load3 = PlateLoad");
    Check(std::dynamic_pointer_cast<BeamPolyLoad>(loads2[4]) != nullptr, "load4 = BeamPolyLoad");
    Check(std::dynamic_pointer_cast<AxialPolyLoad>(loads2[5]) != nullptr, "load5 = AxialPolyLoad");

    // 3) 主要フィールド
    if (auto nl = std::dynamic_pointer_cast<NodeLoad>(loads2[0]))
    {
        Check(nl->data.id == 2, "NodeLoad nodeId の往復");
        Check(std::abs(nl->data.Py() + 1000.0) < 1e-9 && std::abs(nl->data.Mz() - 5.0) < 1e-12,
              "NodeLoad 成分の往復");
    }
    if (auto inf = std::dynamic_pointer_cast<InertialForce>(loads2[1]))
        Check(std::abs(inf->accels.y + 9806.65) < 1e-9, "InertialForce 加速度の往復");
    if (auto nbf = std::dynamic_pointer_cast<NodeBodyForce>(loads2[2]))
    {
        Check(nbf->node != nullptr && nbf->node->id == 2, "NodeBodyForce 節点の往復");
        Check(static_cast<unsigned int>(nbf->selector) == 3u, "NodeBodyForce selector の往復");
    }
    if (auto pl = std::dynamic_pointer_cast<PlateLoad>(loads2[3]))
    {
        Check(pl->element != nullptr && pl->element->id == 2, "PlateLoad 要素の往復");
        Check(pl->load_vecs.size() == 4, "PlateLoad load_vecs 数の往復");
        if (pl->load_vecs.size() == 4)
            Check(std::abs(pl->load_vecs[0].z + 0.01) < 1e-12, "PlateLoad ベクトルの往復");
    }
    if (auto bpl = std::dynamic_pointer_cast<BeamPolyLoad>(loads2[4]))
    {
        Check(bpl->element != nullptr && bpl->element->id == 0, "BeamPolyLoad 要素の往復");
        Check(bpl->axis == BeamLoadAxis::YAxis, "BeamPolyLoad axis の往復");
        Check(bpl->w.size() == 2 && std::abs(bpl->w[0] + 5.0) < 1e-12, "BeamPolyLoad w の往復");
        Check(bpl->params.size() == 2 && std::abs(bpl->params[1] - 1.0) < 1e-12, "BeamPolyLoad params の往復");
    }
    if (auto apl = std::dynamic_pointer_cast<AxialPolyLoad>(loads2[5]))
    {
        Check(apl->element != nullptr && apl->element->id == 1, "AxialPolyLoad 要素の往復");
        Check(apl->w.size() == 2 && std::abs(apl->w[0] - 2.0) < 1e-12, "AxialPolyLoad w の往復");
    }

    // 4) 復元した荷重で解析が走ること
    try
    {
        FELinearStaticOp op(model, loads2);
        op.Compute();
        std::vector<Displacement> disp = op.GetDisplacements();
        Check(disp.size() == model->Nodes.size(), "FELinearStaticOp が解けて変位を返す");
    }
    catch (const std::exception &e)
    {
        Check(false, std::string("FELinearStaticOp 実行で例外: ") + e.what());
    }

    std::cout << "\n";
    if (g_failures == 0)
    {
        std::cout << "RESULT: PASS (all checks)\n";
        return 0;
    }
    std::cout << "RESULT: FAIL (" << g_failures << " checks failed)\n";
    return 1;
}
