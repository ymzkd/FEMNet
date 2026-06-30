// FEModel テキスト形式 Save/Load の往復テスト
//
// 全要素種別(Truss/Beam/ComplexBeam/TriMembrane/QuadMembrane/DKT/DKQ) +
// 支点条件 + 剛体連結 + 非既定の重力加速度を含むモデルを構築し、
//   save -> load -> save の2ファイルがバイト一致すること、および
//   主要フィールドが一致すること
// を検証する。

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <memory>

#include "Model.h"
#include "Elements/Elements.h"
#include "RigidLink.h"

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

FEModel BuildSampleModel()
{
    FEModel m;
    m.GraityAccel = 9800.0; // 非既定値

    // --- Nodes (0..11) ---
    double coords[12][3] = {
        {0, 0, 0}, {1000, 0, 0}, {2000, 0, 0}, {3000, 0, 0},
        {0, 0, 1000}, {1000, 0, 1000}, {1000, 1000, 1000}, {0, 1000, 1000},
        {0, 0, 2000}, {1000, 0, 2000}, {1000, 1000, 2000}, {0, 1000, 2000}};
    for (int i = 0; i < 12; i++)
        m.Nodes.push_back(Node(i, coords[i][0], coords[i][1], coords[i][2]));

    // 支点: node0 完全固定
    for (int i = 0; i < 6; i++)
        m.Nodes[0].Fix.flags[i] = true;
    // node1 に lockflags のパターンを設定(往復確認用)
    m.Nodes[1].Fix.lockflags.flags[3] = true;
    m.Nodes[1].Fix.lockflags.flags[5] = true;

    // --- Materials ---
    Material mat0(205000.0, 0.3, 7.85e-9);
    mat0.id = 0;
    m.Materials.push_back(mat0);
    Material mat1(70000.0, 0.33, 2.7e-9);
    mat1.id = 1;
    m.Materials.push_back(mat1);

    // --- Sections ---
    Section s0;
    s0.id = 0;
    s0.A = 1000.0; s0.Iy = 2.0e6; s0.Iz = 3.0e6; s0.Iyz = 1.5e5; s0.K = 4.0e6;
    m.Sections.push_back(s0);
    Section s1;
    s1.id = 1;
    s1.A = 500.0; s1.Iy = 1.0e6; s1.Iz = 1.2e6; s1.Iyz = 0.0; s1.K = 2.0e6;
    m.Sections.push_back(s1);

    // --- Elements (Nodes/Sections は以後 resize しないのでポインタ安定) ---
    m.Elements.push_back(std::make_shared<TrussElement>(
        0, &m.Nodes[0], &m.Nodes[1], &m.Sections[0], m.Materials[0]));
    m.Elements.push_back(std::make_shared<BeamElement>(
        1, &m.Nodes[1], &m.Nodes[2], &m.Sections[1], m.Materials[0], 0.5));

    auto cb = std::make_shared<ComplexBeamElement>(
        2, &m.Nodes[2], &m.Nodes[3], &m.Sections[0], m.Materials[1], 0.25);
    cb->Lambda_bzi = 0.7; cb->Lambda_bzj = 0.8;
    cb->Lambda_syi = 0.6; cb->Lambda_syj = 0.9;
    cb->lzi = 1.0; cb->lzj = 2.0; cb->lyi = 3.0; cb->lyj = 4.0;
    m.Elements.push_back(cb);

    Thickness th(10.0, 12.0, 8.0); // plane/plate/weight を別値に
    m.Elements.push_back(std::make_shared<TriPlaneElement>(
        3, &m.Nodes[4], &m.Nodes[5], &m.Nodes[6], th, m.Materials[0], 0.1));
    m.Elements.push_back(std::make_shared<QuadPlaneElement>(
        4, &m.Nodes[4], &m.Nodes[5], &m.Nodes[6], &m.Nodes[7], th, m.Materials[1], 0.2));
    m.Elements.push_back(std::make_shared<TriPlateElement>(
        5, &m.Nodes[8], &m.Nodes[9], &m.Nodes[10], th, m.Materials[0], 0.3));
    m.Elements.push_back(std::make_shared<QuadPlateElement>(
        6, &m.Nodes[8], &m.Nodes[9], &m.Nodes[10], &m.Nodes[11], th, m.Materials[1], 0.4));

    // --- RigidLink ---
    RigidLink link(true, false, true, false, true, false);
    link.Master = &m.Nodes[4];
    link.Slaves.push_back(m.Nodes[5]);
    link.Slaves.push_back(m.Nodes[6]);
    m.RigidLinkData = std::make_shared<RigidLinks>(std::vector<RigidLink>{link});

    return m;
}

} // namespace

int main()
{
    std::cout << "ModelIO round-trip test\n";

    const std::string path_a = "modelio_a.txt";
    const std::string path_b = "modelio_b.txt";

    FEModel m1 = BuildSampleModel();
    m1.Save(path_a);

    FEModel m2;
    m2.Load(path_a);
    m2.Save(path_b);

    // 1) save -> load -> save のバイト一致
    std::string a = ReadFile(path_a);
    std::string b = ReadFile(path_b);
    Check(!a.empty(), "出力ファイルが空でない");
    Check(a == b, "save->load->save がバイト一致する");

    // 2) 主要フィールドの一致
    Check(m2.Nodes.size() == 12, "Node 数 = 12");
    Check(m2.Materials.size() == 2, "Material 数 = 2");
    Check(m2.Sections.size() == 2, "Section 数 = 2");
    Check(m2.Elements.size() == 7, "Element 数 = 7");
    Check(std::abs(m2.GraityAccel - 9800.0) < 1e-12, "重力加速度の往復");

    // 支点
    bool node0_fixed = true;
    for (int i = 0; i < 6; i++)
        node0_fixed = node0_fixed && m2.Nodes[0].Fix.flags[i];
    Check(node0_fixed, "node0 完全固定の往復");
    Check(m2.Nodes[1].Fix.lockflags.flags[3] && m2.Nodes[1].Fix.lockflags.flags[5] &&
              !m2.Nodes[1].Fix.lockflags.flags[0],
          "node1 lockflags の往復");

    // 座標
    Check(std::abs(m2.Nodes[10].Location.x - 1000.0) < 1e-9 &&
              std::abs(m2.Nodes[10].Location.y - 1000.0) < 1e-9 &&
              std::abs(m2.Nodes[10].Location.z - 2000.0) < 1e-9,
          "node10 座標の往復");

    // 要素種別
    Check(m2.Elements[0]->Type() == ElementType::Truss, "elem0 = Truss");
    Check(m2.Elements[1]->Type() == ElementType::Beam, "elem1 = Beam");
    Check(m2.Elements[2]->Type() == ElementType::ComplexBeam, "elem2 = ComplexBeam");
    Check(m2.Elements[3]->Type() == ElementType::TriMembrane, "elem3 = TriMembrane");
    Check(m2.Elements[4]->Type() == ElementType::QuadMembrane, "elem4 = QuadMembrane");
    Check(m2.Elements[5]->Type() == ElementType::DKT, "elem5 = DKT");
    Check(m2.Elements[6]->Type() == ElementType::DKQ, "elem6 = DKQ");

    // Beam beta / 接続 / 材料・断面 id
    if (auto *beam = dynamic_cast<BeamElement *>(m2.Elements[1].get()))
    {
        Check(std::abs(beam->Beta - 0.5) < 1e-12, "beam beta の往復");
        Check(beam->Nodes[0]->id == 1 && beam->Nodes[1]->id == 2, "beam 接続の往復");
        Check(beam->Sec->id == 1, "beam section id の往復");
        // 材料は要素ごとにインライン保存されるため定数値で検証(Mat.id は復元対象外)
        Check(std::abs(beam->Mat.Young - 205000.0) < 1e-9 &&
                  std::abs(beam->Mat.Poisson - 0.3) < 1e-12,
              "beam material 定数の往復");
    }
    else
        Check(false, "elem1 を BeamElement にキャスト");

    // ComplexBeam バネパラメータ
    if (auto *c = dynamic_cast<ComplexBeamElement *>(m2.Elements[2].get()))
    {
        Check(std::abs(c->Lambda_bzi - 0.7) < 1e-12 &&
                  std::abs(c->Lambda_syj - 0.9) < 1e-12 &&
                  std::abs(c->lyj - 4.0) < 1e-12,
              "ComplexBeam バネパラメータの往復");
        Check(std::abs(c->Beta - 0.25) < 1e-12, "ComplexBeam beta の往復");
    }
    else
        Check(false, "elem2 を ComplexBeamElement にキャスト");

    // 板厚(3値)の往復
    if (auto *pe = dynamic_cast<PlaneElementBase *>(m2.Elements[6].get()))
    {
        Check(std::abs(pe->thickness.plane_thick - 10.0) < 1e-12 &&
                  std::abs(pe->thickness.plate_thick - 12.0) < 1e-12 &&
                  std::abs(pe->thickness.weight_thick - 8.0) < 1e-12,
              "DKQ 板厚(3値)の往復");
        Check(std::abs(pe->Beta - 0.4) < 1e-12, "DKQ beta の往復");
    }
    else
        Check(false, "elem6 を PlaneElementBase にキャスト");

    // DKQ 接続
    Check(m2.Elements[6]->NodesList().size() == 4 &&
              m2.Elements[6]->NodesList()[3]->id == 11,
          "DKQ 接続の往復");

    // RigidLink
    Check(m2.RigidLinkData->links.size() == 1, "RigidLink 数 = 1");
    if (m2.RigidLinkData->links.size() == 1)
    {
        RigidLink &l = m2.RigidLinkData->links[0];
        Check(l.flags[0] && !l.flags[1] && l.flags[2] && !l.flags[3] && l.flags[4] && !l.flags[5],
              "RigidLink flags の往復");
        Check(l.Master && l.Master->id == 4, "RigidLink master の往復");
        Check(l.Slaves.size() == 2 && l.Slaves[0].id == 5 && l.Slaves[1].id == 6,
              "RigidLink slaves の往復");
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
