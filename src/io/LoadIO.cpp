// 荷重(LoadBase 派生群)のテキスト形式ファイル入出力
//
// フォーマット(v1, 空白・改行非依存のトークン列):
//   FEMNET_LOADS_TEXT_V1
//   LOADS <count>
//     NodeLoad      <nodeId> <Px> <Py> <Pz> <Mx> <My> <Mz>
//     InertialForce <ax> <ay> <az>
//     NodeBodyForce <nodeIdx> <ax> <ay> <az> <selector>
//     PlateLoad     <elemIdx> <nNodes> <vx vy vz>*nNodes
//     BeamPolyLoad  <elemIdx> <axis> <nw> <w...> <np> <params...>
//     AxialPolyLoad <elemIdx> <nw> <w...> <np> <params...>
//     LinePolyLoad  <nodeIdx_i> <nodeIdx_j> <dx> <dy> <dz> <nw> <w...> <np> <params...>
//   END
//
// 設計方針: 要素・節点はポインタ参照だが、保存時はモデルを持たないため
//   ->id (= ベクタ添字, モデルの不変条件) を識別子として書き出す。
//   読み込み時はモデルを与え &model.Nodes[i] / model.Elements[i] で復元する。

#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <cmath>

#include "LoadIO.h"
#include "Elements/Elements.h"

namespace {

const char *kMagic = "FEMNET_LOADS_TEXT_V1";

// BeamLoadAxis <-> 名前
std::string AxisName(BeamLoadAxis a)
{
    switch (a)
    {
    case BeamLoadAxis::YAxis: return "YAxis";
    case BeamLoadAxis::ZAxis: return "ZAxis";
    case BeamLoadAxis::XAxis: return "XAxis";
    default:
        throw std::runtime_error("LoadIO: 未対応の BeamLoadAxis です");
    }
}

BeamLoadAxis ParseAxis(const std::string &name)
{
    if (name == "YAxis") return BeamLoadAxis::YAxis;
    if (name == "ZAxis") return BeamLoadAxis::ZAxis;
    if (name == "XAxis") return BeamLoadAxis::XAxis;
    throw std::runtime_error("LoadIO: 不明な BeamLoadAxis タグ: " + name);
}

template <typename T>
T ReadToken(std::istream &is, const char *what)
{
    T v;
    if (!(is >> v))
        throw std::runtime_error(std::string("LoadIO: ") + what + " の読み取りに失敗しました");
    return v;
}

void ExpectKeyword(std::istream &is, const char *keyword)
{
    std::string tok = ReadToken<std::string>(is, keyword);
    if (tok != keyword)
        throw std::runtime_error(std::string("LoadIO: '") + keyword +
                                 "' を期待しましたが '" + tok + "' でした");
}

// double 列を書く
void WriteDoubleVec(std::ostream &os, const std::vector<double> &v)
{
    os << " " << v.size();
    for (double d : v)
        os << " " << d;
}

// double 列を読む
std::vector<double> ReadDoubleVec(std::istream &is, const char *what)
{
    int n = ReadToken<int>(is, what);
    if (n < 0)
        throw std::runtime_error(std::string("LoadIO: ") + what + " の個数が不正です");
    std::vector<double> v;
    v.reserve(n);
    for (int i = 0; i < n; i++)
        v.push_back(ReadToken<double>(is, what));
    return v;
}

// 外部入力の検証は入出力側で行い、荷重クラスは計算に専念させる。
void CheckLinePolyLoad(const LinePolyLoad &load)
{
    if (load.node_i == nullptr || load.node_j == nullptr)
        throw std::runtime_error("LoadIO: LinePolyLoad の節点が未設定です");

    const double length = (load.node_j->Location - load.node_i->Location).norm();
    if (!std::isfinite(length) || length <= 0.0)
        throw std::runtime_error("LoadIO: LinePolyLoad の線分長さが不正です");

    Vector direction = load.direction;
    const double direction_length = direction.norm();
    if (!std::isfinite(load.direction.x) || !std::isfinite(load.direction.y) ||
        !std::isfinite(load.direction.z) || !std::isfinite(direction_length) || direction_length <= 0.0)
        throw std::runtime_error("LoadIO: LinePolyLoad の方向が不正です");

    if (load.w.size() != load.params.size() || load.w.size() < 2)
        throw std::runtime_error("LoadIO: LinePolyLoad の荷重値と位置は同数で2点以上が必要です");

    for (size_t k = 0; k < load.w.size(); ++k)
    {
        if (!std::isfinite(load.w[k]))
            throw std::runtime_error("LoadIO: LinePolyLoad の荷重値が不正です");
        if (!std::isfinite(load.params[k]) || load.params[k] < 0.0 || load.params[k] > 1.0)
            throw std::runtime_error("LoadIO: LinePolyLoad の位置が範囲外です");
        if (k > 0 && load.params[k] < load.params[k - 1])
            throw std::runtime_error("LoadIO: LinePolyLoad の位置が逆順です");
    }
}

} // namespace

void SaveLoads(const std::string &path,
               const std::vector<std::shared_ptr<LoadBase>> &loads)
{
    // 荷重の検証で例外になっても既存のファイルを壊さないよう、全体を書き出せてからファイルを開く
    std::ostringstream os;
    os << std::setprecision(17);
    os << kMagic << "\n";
    os << "LOADS " << loads.size() << "\n";

    for (const std::shared_ptr<LoadBase> &load : loads)
    {
        // 派生の深い順に判定する
        if (auto axial = std::dynamic_pointer_cast<AxialPolyLoad>(load))
        {
            if (axial->element == nullptr)
                throw std::runtime_error("LoadIO: AxialPolyLoad の要素が未設定です");
            os << "AxialPolyLoad " << axial->element->id;
            WriteDoubleVec(os, axial->w);
            WriteDoubleVec(os, axial->params);
            os << "\n";
        }
        else if (auto beam = std::dynamic_pointer_cast<BeamPolyLoad>(load))
        {
            if (beam->element == nullptr)
                throw std::runtime_error("LoadIO: BeamPolyLoad の要素が未設定です");
            os << "BeamPolyLoad " << beam->element->id << " " << AxisName(beam->axis);
            WriteDoubleVec(os, beam->w);
            WriteDoubleVec(os, beam->params);
            os << "\n";
        }
        else if (auto line = std::dynamic_pointer_cast<LinePolyLoad>(load))
        {
            CheckLinePolyLoad(*line);
            os << "LinePolyLoad " << line->node_i->id << " " << line->node_j->id
                << " " << line->direction.x << " " << line->direction.y << " " << line->direction.z;
            WriteDoubleVec(os, line->w);
            WriteDoubleVec(os, line->params);
            os << "\n";
        }
        else if (auto nbf = std::dynamic_pointer_cast<NodeBodyForce>(load))
        {
            if (nbf->node == nullptr)
                throw std::runtime_error("LoadIO: NodeBodyForce の節点が未設定です");
            os << "NodeBodyForce " << nbf->node->id
                << " " << nbf->Accels.x << " " << nbf->Accels.y << " " << nbf->Accels.z
                << " " << static_cast<unsigned int>(nbf->selector) << "\n";
        }
        else if (auto inertial = std::dynamic_pointer_cast<InertialForce>(load))
        {
            os << "InertialForce "
                << inertial->accels.x << " " << inertial->accels.y << " " << inertial->accels.z << "\n";
        }
        else if (auto nl = std::dynamic_pointer_cast<NodeLoad>(load))
        {
            os << "NodeLoad " << nl->data.id
                << " " << nl->data.Px() << " " << nl->data.Py() << " " << nl->data.Pz()
                << " " << nl->data.Mx() << " " << nl->data.My() << " " << nl->data.Mz() << "\n";
        }
        else if (auto pl = std::dynamic_pointer_cast<PlateLoad>(load))
        {
            if (pl->element == nullptr)
                throw std::runtime_error("LoadIO: PlateLoad の要素が未設定です");
            os << "PlateLoad " << pl->element->id << " " << pl->load_vecs.size();
            for (const Vector &v : pl->load_vecs)
                os << " " << v.x << " " << v.y << " " << v.z;
            os << "\n";
        }
        else
        {
            throw std::runtime_error("LoadIO: 未対応の荷重種別です");
        }
    }

    os << "END\n";

    std::ofstream ofs(path, std::ios::out | std::ios::trunc);
    if (!ofs)
        throw std::runtime_error("LoadIO: ファイルを開けません: " + path);

    ofs << os.str();

    if (!ofs)
        throw std::runtime_error("LoadIO: ファイル書き込み中にエラーが発生しました: " + path);
}

std::vector<std::shared_ptr<LoadBase>> LoadLoads(const std::string &path, FEModel &model)
{
    std::ifstream ifs(path, std::ios::in);
    if (!ifs)
        throw std::runtime_error("LoadIO: ファイルを開けません: " + path);

    std::string magic = ReadToken<std::string>(ifs, "magic");
    if (magic != kMagic)
        throw std::runtime_error("LoadIO: フォーマット識別子が一致しません: " + magic);

    const int node_count = static_cast<int>(model.Nodes.size());
    const int elem_count = static_cast<int>(model.Elements.size());

    auto node_ptr = [&](int idx) -> Node * {
        if (idx < 0 || idx >= node_count)
            throw std::runtime_error("LoadIO: 節点インデックスが範囲外です: " + std::to_string(idx));
        return &model.Nodes[idx];
    };
    auto beam_ptr = [&](int idx) -> BeamElement * {
        if (idx < 0 || idx >= elem_count)
            throw std::runtime_error("LoadIO: 要素インデックスが範囲外です: " + std::to_string(idx));
        BeamElement *b = dynamic_cast<BeamElement *>(model.Elements[idx].get());
        if (b == nullptr)
            throw std::runtime_error("LoadIO: 要素 " + std::to_string(idx) + " は梁要素ではありません");
        return b;
    };
    auto plane_ptr = [&](int idx) -> PlaneElementBase * {
        if (idx < 0 || idx >= elem_count)
            throw std::runtime_error("LoadIO: 要素インデックスが範囲外です: " + std::to_string(idx));
        PlaneElementBase *p = dynamic_cast<PlaneElementBase *>(model.Elements[idx].get());
        if (p == nullptr)
            throw std::runtime_error("LoadIO: 要素 " + std::to_string(idx) + " は平面・板要素ではありません");
        return p;
    };

    ExpectKeyword(ifs, "LOADS");
    int count = ReadToken<int>(ifs, "LOADS 個数");

    std::vector<std::shared_ptr<LoadBase>> loads;
    loads.reserve(count);

    for (int k = 0; k < count; k++)
    {
        std::string kind = ReadToken<std::string>(ifs, "荷重種別");

        if (kind == "NodeLoad")
        {
            int id = ReadToken<int>(ifs, "NodeLoad nodeId");
            double px = ReadToken<double>(ifs, "Px");
            double py = ReadToken<double>(ifs, "Py");
            double pz = ReadToken<double>(ifs, "Pz");
            double mx = ReadToken<double>(ifs, "Mx");
            double my = ReadToken<double>(ifs, "My");
            double mz = ReadToken<double>(ifs, "Mz");
            loads.push_back(std::make_shared<NodeLoad>(id, px, py, pz, mx, my, mz));
        }
        else if (kind == "InertialForce")
        {
            double ax = ReadToken<double>(ifs, "ax");
            double ay = ReadToken<double>(ifs, "ay");
            double az = ReadToken<double>(ifs, "az");
            loads.push_back(std::make_shared<InertialForce>(ax, ay, az));
        }
        else if (kind == "NodeBodyForce")
        {
            int idx = ReadToken<int>(ifs, "NodeBodyForce nodeIdx");
            double ax = ReadToken<double>(ifs, "ax");
            double ay = ReadToken<double>(ifs, "ay");
            double az = ReadToken<double>(ifs, "az");
            unsigned int sel = ReadToken<unsigned int>(ifs, "selector");
            auto nbf = std::make_shared<NodeBodyForce>(node_ptr(idx), ax, ay, az);
            nbf->selector = static_cast<BodyForaceSelector>(sel);
            loads.push_back(nbf);
        }
        else if (kind == "PlateLoad")
        {
            int idx = ReadToken<int>(ifs, "PlateLoad elemIdx");
            int nn = ReadToken<int>(ifs, "PlateLoad nNodes");
            if (nn < 0)
                throw std::runtime_error("LoadIO: PlateLoad の節点数が不正です");
            auto pl = std::make_shared<PlateLoad>();
            pl->element = plane_ptr(idx);
            for (int i = 0; i < nn; i++)
            {
                double vx = ReadToken<double>(ifs, "vx");
                double vy = ReadToken<double>(ifs, "vy");
                double vz = ReadToken<double>(ifs, "vz");
                pl->load_vecs.push_back(Vector(vx, vy, vz));
            }
            loads.push_back(pl);
        }
        else if (kind == "BeamPolyLoad")
        {
            int idx = ReadToken<int>(ifs, "BeamPolyLoad elemIdx");
            std::string axis = ReadToken<std::string>(ifs, "axis");
            std::vector<double> w = ReadDoubleVec(ifs, "BeamPolyLoad w");
            std::vector<double> params = ReadDoubleVec(ifs, "BeamPolyLoad params");
            loads.push_back(std::make_shared<BeamPolyLoad>(w, params, beam_ptr(idx), ParseAxis(axis)));
        }
        else if (kind == "LinePolyLoad")
        {
            int idx_i = ReadToken<int>(ifs, "LinePolyLoad nodeIdx_i");
            int idx_j = ReadToken<int>(ifs, "LinePolyLoad nodeIdx_j");
            double dx = ReadToken<double>(ifs, "LinePolyLoad dx");
            double dy = ReadToken<double>(ifs, "LinePolyLoad dy");
            double dz = ReadToken<double>(ifs, "LinePolyLoad dz");
            std::vector<double> w = ReadDoubleVec(ifs, "LinePolyLoad w");
            std::vector<double> params = ReadDoubleVec(ifs, "LinePolyLoad params");
            auto line = std::make_shared<LinePolyLoad>(
                w, params, node_ptr(idx_i), node_ptr(idx_j), Vector(dx, dy, dz));
            CheckLinePolyLoad(*line);
            loads.push_back(line);
        }
        else if (kind == "AxialPolyLoad")
        {
            int idx = ReadToken<int>(ifs, "AxialPolyLoad elemIdx");
            std::vector<double> w = ReadDoubleVec(ifs, "AxialPolyLoad w");
            std::vector<double> params = ReadDoubleVec(ifs, "AxialPolyLoad params");
            loads.push_back(std::make_shared<AxialPolyLoad>(w, params, beam_ptr(idx)));
        }
        else
        {
            throw std::runtime_error("LoadIO: 不明な荷重種別タグ: " + kind);
        }
    }

    ExpectKeyword(ifs, "END");

    return loads;
}
