// 構造モデル(FEModel)のテキスト形式ファイル入出力
//
// 設計方針:
//   コードベースは「id == ベクタ添字」を前提に動作するが、AddMaterial/
//   AddSection/AddNodeXYZ 等は .id を設定しない(-1 のまま)ことがある。
//   そのため本IOでは .id フィールドに依存せず、ベクタの位置(インデックス)
//   を正準IDとして用いる。
//   - 節点・断面は要素からポインタ参照されるため、要素側はその位置インデックス
//     (ポインタ演算で算出)で参照する。
//   - 材料は要素が値コピーで保持し位置を復元できないため、要素ごとに
//     材料定数(Young/Poisson/dense)をインライン保存する。
//
// フォーマット(v3, 空白・改行非依存のトークン列):
//   FEMNET_MODEL_TEXT_V3
//   GRAVITY <g>
//   NODES <count>      ... <idx> <x> <y> <z> <btype0..5>
//     btype は ConstraintType (0=Free, 1=Fix)。
//     V2 は <fix0..5> <lock0..5> (0/1) で、fix をそのまま btype として読み、lock は無視する
//     (剛性の付かない回転自由度の拘束は解析側が剛性行列から判定するため)。
//   MATERIALS <count>  ... <idx> <Young> <Poisson> <dense>
//   SECTIONS <count>   ... <idx> <A> <Iy> <Iz> <Iyz> <K>
//   ELEMENTS <count>
//     棒系:   <eid> <Kind> <niIdx> <njIdx> <secIdx> <Young> <Poisson> <dense>
//             [<beta>] [<バネ12値>]
//     面・板: <eid> <Kind> <nodeIdx...> <Young> <Poisson> <dense>
//             <tplane> <tplate> <tweight> <beta>
//     支点ばね: <eid> SupportSpring <nodeIdx> <k0..5>
//   RIGIDLINKS <count>
//     <flag0..5> <mx> <my> <mz> <mbtype0..5> <slaveCount> <slaveIdx...>
//     マスタは仮想節点(剛床の重心など)を含むため、座標＋支持条件(btype)を
//     実体で常にインライン保持。節点と同じ並び(V2 は fix/lock)。
//   END
//
// 要素種別タグは ElementType に対応する名前(Truss/Beam/ComplexBeam/
// TriMembrane/QuadMembrane/DKT/DKQ/SupportSpring)。

#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <string>

#include "Model.h"
#include "Elements/Elements.h"
#include "RigidLink.h"

namespace {

const char *kMagic = "FEMNET_MODEL_TEXT_V3";
const char *kMagicV2 = "FEMNET_MODEL_TEXT_V2"; // 旧形式(読み込みのみ対応)

// 支持条件の書き出し/読み込み(節点・剛体リンクのマスタ節点で共通)
void WriteSupport(std::ostream &os, const Support &sup)
{
    for (int i = 0; i < 6; i++)
        os << " " << static_cast<int>(sup.BoundaryTypes[i]);
}

// ElementType -> シリアライズ用タグ名
std::string KindName(ElementType t)
{
    switch (t)
    {
    case ElementType::Truss:        return "Truss";
    case ElementType::Beam:         return "Beam";
    case ElementType::ComplexBeam:  return "ComplexBeam";
    case ElementType::TriMembrane:  return "TriMembrane";
    case ElementType::QuadMembrane: return "QuadMembrane";
    case ElementType::DKT:          return "DKT";
    case ElementType::DKQ:          return "DKQ";
    case ElementType::SupportSpring: return "SupportSpring";
    default:
        throw std::runtime_error("ModelIO: 未対応の要素種別です (Type=" +
                                 std::to_string(static_cast<int>(t)) + ")");
    }
}

// タグ名 -> ElementType
ElementType ParseKind(const std::string &name)
{
    if (name == "Truss")        return ElementType::Truss;
    if (name == "Beam")         return ElementType::Beam;
    if (name == "ComplexBeam")  return ElementType::ComplexBeam;
    if (name == "TriMembrane")  return ElementType::TriMembrane;
    if (name == "QuadMembrane") return ElementType::QuadMembrane;
    if (name == "DKT")          return ElementType::DKT;
    if (name == "DKQ")          return ElementType::DKQ;
    if (name == "SupportSpring") return ElementType::SupportSpring;
    throw std::runtime_error("ModelIO: 不明な要素種別タグ: " + name);
}

// 次の必須トークンを読む。読めなければ例外。
template <typename T>
T ReadToken(std::istream &is, const char *what)
{
    T v;
    if (!(is >> v))
        throw std::runtime_error(std::string("ModelIO: ") + what + " の読み取りに失敗しました");
    return v;
}

// legacy_v2: V2 形式(<fix0..5> <lock0..5>)として読む。lock は読み捨てる。
void ReadSupport(std::istream &is, Support &sup, bool legacy_v2)
{
    for (int i = 0; i < 6; i++)
    {
        int t = ReadToken<int>(is, "Support boundary type");
        if (t < static_cast<int>(ConstraintType::Free) || t > static_cast<int>(ConstraintType::Fix))
            throw std::runtime_error("ModelIO: 支持条件の種別が範囲外です (0=Free, 1=Fix): " +
                                     std::to_string(t));
        sup.BoundaryTypes[i] = static_cast<ConstraintType>(t);
    }
    if (legacy_v2)
    {
        for (int i = 0; i < 6; i++)
            ReadToken<int>(is, "Support lock (V2, 無視)");
    }
}

// 期待するキーワードを読み、一致しなければ例外。
void ExpectKeyword(std::istream &is, const char *keyword)
{
    std::string tok = ReadToken<std::string>(is, keyword);
    if (tok != keyword)
        throw std::runtime_error(std::string("ModelIO: '") + keyword +
                                 "' を期待しましたが '" + tok + "' でした");
}

} // namespace

void FEModel::Save(const std::string &path)
{
    std::ofstream ofs(path, std::ios::out | std::ios::trunc);
    if (!ofs)
        throw std::runtime_error("ModelIO: ファイルを開けません: " + path);

    // double を完全に往復させるため有効桁17桁
    ofs << std::setprecision(17);

    // 節点・断面のポインタ→位置インデックス算出用の基点
    const Node *node_base = Nodes.empty() ? nullptr : &Nodes[0];
    const Section *sec_base = Sections.empty() ? nullptr : &Sections[0];

    ofs << kMagic << "\n";
    ofs << "GRAVITY " << GraityAccel << "\n";

    // --- Nodes (位置インデックスを ID として書き出す) ---
    ofs << "NODES " << Nodes.size() << "\n";
    for (size_t i = 0; i < Nodes.size(); i++)
    {
        const Node &n = Nodes[i];
        ofs << i << " "
            << n.Location.x << " " << n.Location.y << " " << n.Location.z;
        WriteSupport(ofs, n.Fix);
        ofs << "\n";
    }

    // --- Materials ---
    ofs << "MATERIALS " << Materials.size() << "\n";
    for (size_t i = 0; i < Materials.size(); i++)
    {
        const Material &m = Materials[i];
        ofs << i << " " << m.Young << " " << m.Poisson << " " << m.dense << "\n";
    }

    // --- Sections ---
    ofs << "SECTIONS " << Sections.size() << "\n";
    for (size_t i = 0; i < Sections.size(); i++)
    {
        const Section &s = Sections[i];
        ofs << i << " " << s.A << " " << s.Iy << " " << s.Iz << " "
            << s.Iyz << " " << s.K << "\n";
    }

    // --- Elements ---
    ofs << "ELEMENTS " << Elements.size() << "\n";
    for (const std::shared_ptr<ElementBase> &el : Elements)
    {
        ElementType t = el->Type();
        ofs << el->id << " " << KindName(t);

        if (IsBarType(t))
        {
            BarElementBase *bar = dynamic_cast<BarElementBase *>(el.get());
            long ni = bar->Nodes[0] - node_base;
            long nj = bar->Nodes[1] - node_base;
            long sec = bar->Sec - sec_base;
            ofs << " " << ni << " " << nj << " " << sec
                << " " << el->Mat.Young << " " << el->Mat.Poisson << " " << el->Mat.dense;

            if (t == ElementType::Beam || t == ElementType::ComplexBeam)
            {
                BeamElement *beam = dynamic_cast<BeamElement *>(el.get());
                ofs << " " << beam->Beta;
            }
            if (t == ElementType::ComplexBeam)
            {
                ComplexBeamElement *cb = dynamic_cast<ComplexBeamElement *>(el.get());
                ofs << " " << cb->Lambda_bzi << " " << cb->Lambda_bzj
                    << " " << cb->Lambda_syi << " " << cb->Lambda_syj
                    << " " << cb->Lambda_byi << " " << cb->Lambda_byj
                    << " " << cb->Lambda_szi << " " << cb->Lambda_szj
                    << " " << cb->lzi << " " << cb->lzj
                    << " " << cb->lyi << " " << cb->lyj;
            }
        }
        else if (t == ElementType::SupportSpring)
        {
            SupportSpringElement *sp = dynamic_cast<SupportSpringElement *>(el.get());
            ofs << " " << (sp->Nodes[0] - node_base);
            for (int i = 0; i < 6; i++)
                ofs << " " << sp->K[i];
        }
        else // 平面・板要素
        {
            PlaneElementBase *pe = dynamic_cast<PlaneElementBase *>(el.get());
            std::vector<Node *> nodes = el->NodesList();
            for (Node *nd : nodes)
                ofs << " " << (nd - node_base);
            ofs << " " << el->Mat.Young << " " << el->Mat.Poisson << " " << el->Mat.dense
                << " " << pe->thickness.plane_thick
                << " " << pe->thickness.plate_thick
                << " " << pe->thickness.weight_thick
                << " " << pe->Beta;
        }
        ofs << "\n";
    }

    // --- RigidLinks ---
    ofs << "RIGIDLINKS " << RigidLinkData->links.size() << "\n";
    for (const RigidLink &link : RigidLinkData->links)
    {
        for (int i = 0; i < 6; i++)
            ofs << (link.flags[i] ? 1 : 0) << " ";

        // マスタ節点: 実体(値)で保持しているため、座標＋拘束(fix/lock)を常にインライン出力。
        const Node &m = link.Master;
        ofs << m.Location.x << " " << m.Location.y << " " << m.Location.z;
        WriteSupport(ofs, m.Fix);

        ofs << " " << link.Slaves.size();
        for (const Node &slave : link.Slaves)
            ofs << " " << slave.id;
        ofs << "\n";
    }

    ofs << "END\n";

    if (!ofs)
        throw std::runtime_error("ModelIO: ファイル書き込み中にエラーが発生しました: " + path);
}

void FEModel::Load(const std::string &path)
{
    std::ifstream ifs(path, std::ios::in);
    if (!ifs)
        throw std::runtime_error("ModelIO: ファイルを開けません: " + path);

    std::string magic = ReadToken<std::string>(ifs, "magic");
    const bool legacy_v2 = (magic == kMagicV2);
    if (magic != kMagic && !legacy_v2)
        throw std::runtime_error("ModelIO: フォーマット識別子が一致しません: " + magic);

    // 既存内容をクリア
    Nodes.clear();
    Materials.clear();
    Sections.clear();
    Elements.clear();
    RigidLinkData = std::make_shared<RigidLinks>();

    // --- Gravity ---
    ExpectKeyword(ifs, "GRAVITY");
    GraityAccel = ReadToken<double>(ifs, "GRAVITY 値");

    // --- Nodes ---
    ExpectKeyword(ifs, "NODES");
    int node_count = ReadToken<int>(ifs, "NODES 個数");
    Nodes.resize(node_count);
    for (int k = 0; k < node_count; k++)
    {
        int idx = ReadToken<int>(ifs, "Node idx");
        if (idx < 0 || idx >= node_count)
            throw std::runtime_error("ModelIO: Node idx が範囲外です: " + std::to_string(idx));
        double x = ReadToken<double>(ifs, "Node x");
        double y = ReadToken<double>(ifs, "Node y");
        double z = ReadToken<double>(ifs, "Node z");
        Node n(idx, x, y, z); // id = 添字に正規化
        ReadSupport(ifs, n.Fix, legacy_v2);
        Nodes[idx] = n;
    }

    // --- Materials ---
    ExpectKeyword(ifs, "MATERIALS");
    int mat_count = ReadToken<int>(ifs, "MATERIALS 個数");
    Materials.resize(mat_count);
    for (int k = 0; k < mat_count; k++)
    {
        int idx = ReadToken<int>(ifs, "Material idx");
        if (idx < 0 || idx >= mat_count)
            throw std::runtime_error("ModelIO: Material idx が範囲外です: " + std::to_string(idx));
        double young = ReadToken<double>(ifs, "Material Young");
        double poisson = ReadToken<double>(ifs, "Material Poisson");
        double dense = ReadToken<double>(ifs, "Material dense");
        Material m(young, poisson, dense);
        m.id = idx;
        Materials[idx] = m;
    }

    // --- Sections ---
    ExpectKeyword(ifs, "SECTIONS");
    int sec_count = ReadToken<int>(ifs, "SECTIONS 個数");
    Sections.resize(sec_count);
    for (int k = 0; k < sec_count; k++)
    {
        int idx = ReadToken<int>(ifs, "Section idx");
        if (idx < 0 || idx >= sec_count)
            throw std::runtime_error("ModelIO: Section idx が範囲外です: " + std::to_string(idx));
        Section s;
        s.id = idx;
        s.A = ReadToken<double>(ifs, "Section A");
        s.Iy = ReadToken<double>(ifs, "Section Iy");
        s.Iz = ReadToken<double>(ifs, "Section Iz");
        s.Iyz = ReadToken<double>(ifs, "Section Iyz");
        s.K = ReadToken<double>(ifs, "Section K");
        Sections[idx] = s;
    }

    // --- Elements ---
    // この時点で Nodes / Sections は確定済み(以後 resize しない)なので
    // &Nodes[i] / &Sections[i] のポインタは安定。
    ExpectKeyword(ifs, "ELEMENTS");
    int elem_count = ReadToken<int>(ifs, "ELEMENTS 個数");
    Elements.reserve(elem_count);
    for (int k = 0; k < elem_count; k++)
    {
        int id = ReadToken<int>(ifs, "Element id");
        std::string kind = ReadToken<std::string>(ifs, "Element kind");
        ElementType t = ParseKind(kind);

        auto read_material = [&]() -> Material {
            double young = ReadToken<double>(ifs, "Element Young");
            double poisson = ReadToken<double>(ifs, "Element Poisson");
            double dense = ReadToken<double>(ifs, "Element dense");
            return Material(young, poisson, dense);
        };

        if (IsBarType(t))
        {
            int ni = ReadToken<int>(ifs, "Element ni");
            int nj = ReadToken<int>(ifs, "Element nj");
            int sec = ReadToken<int>(ifs, "Element section");
            Material mat = read_material();

            if (t == ElementType::Truss)
            {
                Elements.push_back(std::make_shared<TrussElement>(
                    id, &Nodes[ni], &Nodes[nj], &Sections[sec], mat));
            }
            else
            {
                double beta = ReadToken<double>(ifs, "Element beta");
                if (t == ElementType::Beam)
                {
                    Elements.push_back(std::make_shared<BeamElement>(
                        id, &Nodes[ni], &Nodes[nj], &Sections[sec], mat, beta));
                }
                else // ComplexBeam
                {
                    auto cb = std::make_shared<ComplexBeamElement>(
                        id, &Nodes[ni], &Nodes[nj], &Sections[sec], mat, beta);
                    cb->Lambda_bzi = ReadToken<double>(ifs, "Lambda_bzi");
                    cb->Lambda_bzj = ReadToken<double>(ifs, "Lambda_bzj");
                    cb->Lambda_syi = ReadToken<double>(ifs, "Lambda_syi");
                    cb->Lambda_syj = ReadToken<double>(ifs, "Lambda_syj");
                    cb->Lambda_byi = ReadToken<double>(ifs, "Lambda_byi");
                    cb->Lambda_byj = ReadToken<double>(ifs, "Lambda_byj");
                    cb->Lambda_szi = ReadToken<double>(ifs, "Lambda_szi");
                    cb->Lambda_szj = ReadToken<double>(ifs, "Lambda_szj");
                    cb->lzi = ReadToken<double>(ifs, "lzi");
                    cb->lzj = ReadToken<double>(ifs, "lzj");
                    cb->lyi = ReadToken<double>(ifs, "lyi");
                    cb->lyj = ReadToken<double>(ifs, "lyj");
                    Elements.push_back(cb);
                }
            }
        }
        else if (t == ElementType::SupportSpring)
        {
            int ni = ReadToken<int>(ifs, "SupportSpring node");
            if (ni < 0 || ni >= node_count)
                throw std::runtime_error("ModelIO: SupportSpring の節点が範囲外です: " + std::to_string(ni));
            double k[6];
            for (int i = 0; i < 6; i++)
                k[i] = ReadToken<double>(ifs, "SupportSpring k");
            Elements.push_back(std::make_shared<SupportSpringElement>(
                id, &Nodes[ni], k[0], k[1], k[2], k[3], k[4], k[5]));
        }
        else // 平面・板要素
        {
            int nnum = (t == ElementType::TriMembrane || t == ElementType::DKT) ? 3 : 4;
            int n[4] = {0, 0, 0, 0};
            for (int i = 0; i < nnum; i++)
                n[i] = ReadToken<int>(ifs, "Element node");
            Material mat = read_material();
            double tp = ReadToken<double>(ifs, "thickness plane");
            double tt = ReadToken<double>(ifs, "thickness plate");
            double tw = ReadToken<double>(ifs, "thickness weight");
            double beta = ReadToken<double>(ifs, "Element beta");
            Thickness th(tp, tt, tw);

            switch (t)
            {
            case ElementType::TriMembrane:
                Elements.push_back(std::make_shared<TriPlaneElement>(
                    id, &Nodes[n[0]], &Nodes[n[1]], &Nodes[n[2]], th, mat, beta));
                break;
            case ElementType::QuadMembrane:
                Elements.push_back(std::make_shared<QuadPlaneElement>(
                    id, &Nodes[n[0]], &Nodes[n[1]], &Nodes[n[2]], &Nodes[n[3]], th, mat, beta));
                break;
            case ElementType::DKT:
                Elements.push_back(std::make_shared<TriPlateElement>(
                    id, &Nodes[n[0]], &Nodes[n[1]], &Nodes[n[2]], th, mat, beta));
                break;
            case ElementType::DKQ:
                Elements.push_back(std::make_shared<QuadPlateElement>(
                    id, &Nodes[n[0]], &Nodes[n[1]], &Nodes[n[2]], &Nodes[n[3]], th, mat, beta));
                break;
            default:
                throw std::runtime_error("ModelIO: 未対応の面要素種別です");
            }
        }
    }

    // --- RigidLinks ---
    ExpectKeyword(ifs, "RIGIDLINKS");
    int link_count = ReadToken<int>(ifs, "RIGIDLINKS 個数");
    std::vector<RigidLink> links;
    links.reserve(link_count);
    for (int k = 0; k < link_count; k++)
    {
        RigidLink link;
        for (int i = 0; i < 6; i++)
            link.flags[i] = (ReadToken<int>(ifs, "RigidLink flag") != 0);

        // マスタ節点: 座標＋拘束(fix/lock)を読み、実体(値)として保持する。
        double mx = ReadToken<double>(ifs, "RigidLink master x");
        double my = ReadToken<double>(ifs, "RigidLink master y");
        double mz = ReadToken<double>(ifs, "RigidLink master z");
        Node m(-1, mx, my, mz);
        ReadSupport(ifs, m.Fix, legacy_v2);
        link.Master = m;

        int slave_count = ReadToken<int>(ifs, "RigidLink slave 数");
        for (int s = 0; s < slave_count; s++)
        {
            int slave_idx = ReadToken<int>(ifs, "RigidLink slave idx");
            link.Slaves.push_back(Nodes[slave_idx]);
        }
        links.push_back(link);
    }
    RigidLinkData = std::make_shared<RigidLinks>(links);

    ExpectKeyword(ifs, "END");
}
