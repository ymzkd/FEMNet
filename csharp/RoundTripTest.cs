// FEModel テキスト形式 Save/Load の C# 往復テスト
//
// Truss/Beam/DKT/DKQ 要素 + 支点 + 剛体連結 + 非既定の重力加速度を含むモデルを
// 構築し、save -> load -> save の2ファイルがバイト一致すること、および主要
// フィールドが一致することを検証する。
//
// 注: ComplexBeam/TriMembrane/QuadMembrane は専用の add_*_element ヘルパーが
// 無いため本テストでは扱わない(C++側 modelio_test が全7種を網羅)。
// save->load->save のバイト一致はモデル全体の忠実性を担保する。

using System;
using System.IO;
using FEMNet;

namespace FemNetTest
{
    static class RoundTripTest
    {
        static int failures = 0;

        static void Check(bool cond, string msg)
        {
            if (!cond)
            {
                Console.WriteLine("  [FAIL] " + msg);
                failures++;
            }
        }

        static FEModel BuildModel()
        {
            var m = new FEModel();
            m.GraityAccel = 9800.0; // 非既定値

            double[,] c = {
                {0, 0, 0}, {1000, 0, 0}, {2000, 0, 0},
                {0, 0, 1000}, {1000, 0, 1000}, {1000, 1000, 1000}, {0, 1000, 1000},
            };
            for (int i = 0; i < 7; i++)
                m.AddNode(i, c[i, 0], c[i, 1], c[i, 2]);

            // 支点: node0 完全固定
            m.GetNode(0).Fix.FixAll();

            // 材料・断面
            m.AddMaterialWithDensity(205000.0, 0.3, 7.85e-9);
            m.AddSection(1000.0, 2.0e6, 3.0e6, 4.0e6);

            // 要素
            m.add_truss_element(0, 0, 1, 0, 0);
            m.add_beam_element(1, 1, 2, 0, 0, 0.5);
            m.add_tri_plate_element(2, 3, 4, 5, 12.0, 0);
            m.add_quad_plate_element(3, 3, 4, 5, 6, 12.0, 0);

            // 剛体連結
            var link = new RigidLink(true, false, true, false, true, false);
            link.Master = m.GetNode(3);
            link.Slaves.Add(m.GetNode(4));
            link.Slaves.Add(m.GetNode(5));
            var links = new VectorRigidLink();
            links.Add(link);
            m.RigidLinkData = new RigidLinks(links);

            return m;
        }

        static int Main()
        {
            Console.WriteLine("ModelIO C# round-trip test");

            const string pathA = "csharp_a.txt";
            const string pathB = "csharp_b.txt";

            FEModel m1 = BuildModel();
            m1.Save(pathA);

            var m2 = new FEModel();
            m2.Load(pathA);
            m2.Save(pathB);

            // 1) save -> load -> save のバイト一致
            string a = File.ReadAllText(pathA);
            string b = File.ReadAllText(pathB);
            Check(a.Length > 0, "出力ファイルが空でない");
            Check(a == b, "save->load->save がバイト一致する");

            // 2) 主要フィールド
            Check(m2.NodeNum() == 7, "Node 数 = 7");
            Check(m2.Materials.Count == 1, "Material 数 = 1");
            Check(m2.Sections.Count == 1, "Section 数 = 1");
            Check(m2.Elements.Count == 4, "Element 数 = 4");
            Check(Math.Abs(m2.GraityAccel - 9800.0) < 1e-12, "重力加速度の往復");

            // 支点
            bool node0Fixed = true;
            for (int i = 0; i < 6; i++)
                node0Fixed = node0Fixed && m2.GetNode(0).Fix.Get(i);
            Check(node0Fixed, "node0 完全固定の往復");

            // 座標
            var p = m2.GetNode(5).Location;
            Check(Math.Abs(p.x - 1000.0) < 1e-9 && Math.Abs(p.y - 1000.0) < 1e-9 &&
                      Math.Abs(p.z - 1000.0) < 1e-9,
                  "node5 座標の往復");

            // 要素種別
            Check(m2.Elements[0].Type() == ElementType.Truss, "elem0 = Truss");
            Check(m2.Elements[1].Type() == ElementType.Beam, "elem1 = Beam");
            Check(m2.Elements[2].Type() == ElementType.DKT, "elem2 = DKT");
            Check(m2.Elements[3].Type() == ElementType.DKQ, "elem3 = DKQ");

            // Beam の詳細(GetBeamElement で型付き取得)
            BeamElement beam = m2.GetBeamElement(1);
            Check(beam != null, "GetBeamElement(1)");
            if (beam != null)
            {
                Check(Math.Abs(beam.Beta - 0.5) < 1e-12, "beam beta の往復");
                Check(beam.Sec.id == 0, "beam section id の往復");
                Check(Math.Abs(beam.Mat.Young - 205000.0) < 1e-9, "beam material 定数の往復");
            }

            // 剛体連結
            Check(m2.RigidLinkData.links.Count == 1, "RigidLink 数 = 1");
            if (m2.RigidLinkData.links.Count == 1)
            {
                RigidLink l = m2.RigidLinkData.links[0];
                Check(l.Master != null && l.Master.id == 3, "RigidLink master の往復");
                Check(l.Slaves.Count == 2 && l.Slaves[0].id == 4 && l.Slaves[1].id == 5,
                      "RigidLink slaves の往復");
            }

            // --- 荷重の往復 (LoadIO) ---
            // m2 の要素: 0=Truss, 1=Beam, 2=DKT, 3=DKQ
            const string loadA = "csharp_loads_a.txt";
            const string loadB = "csharp_loads_b.txt";

            var loads = new VectorLoad();
            loads.Add(new NodeLoad(2, 10.0, -1000.0, 0.0, 0.0, 0.0, 5.0));
            loads.Add(new InertialForce(0.0, -9806.65, 0.0));
            loads.Add(new BeamPolyLoad(
                new VectorDouble(new double[] { -5.0, -5.0 }),
                new VectorDouble(new double[] { 0.0, 1.0 }),
                m2.GetBeamElement(1), BeamLoadAxis.YAxis));

            global::FEMNet.FEMNet.SaveLoads(loadA, loads);
            VectorLoad loads2 = global::FEMNet.FEMNet.LoadLoads(loadA, m2);
            global::FEMNet.FEMNet.SaveLoads(loadB, loads2);

            Check(File.ReadAllText(loadA) == File.ReadAllText(loadB),
                  "loads save->load->save がバイト一致する");
            Check(loads2.Count == 3, "荷重数 = 3");

            Console.WriteLine();
            if (failures == 0)
            {
                Console.WriteLine("RESULT: PASS (all checks)");
                return 0;
            }
            Console.WriteLine($"RESULT: FAIL ({failures} checks failed)");
            return 1;
        }
    }
}
