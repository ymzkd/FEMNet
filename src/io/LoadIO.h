#ifndef _LOAD_IO_
#define _LOAD_IO_

#ifndef SWIG
#include <string>
#include <vector>
#include <memory>

#include "Model.h"
#include "LoadComponent.h"
#endif

// 荷重(LoadBase 派生群)のテキスト形式ファイル入出力。
//
// 荷重は FEModel に属さず、解析(FELinearStaticOp)へ
// std::vector<std::shared_ptr<LoadBase>> として渡される。荷重は節点・要素を
// ポインタ参照するため、読み込み時はモデルを与えてポインタを復元する。
//
// 対象: NodeLoad / InertialForce / NodeBodyForce / PlateLoad /
//       BeamPolyLoad / AxialPolyLoad
// (DynamicAccelLoad は LoadBase 派生でなく対象外)

// 荷重リストを path に書き出す。
void SaveLoads(const std::string &path,
               const std::vector<std::shared_ptr<LoadBase>> &loads);

// path から荷重リストを読み込む。要素・節点参照は model から解決する。
std::vector<std::shared_ptr<LoadBase>> LoadLoads(const std::string &path, FEModel &model);

#endif
