#ifndef _ELEMENTS_H_
#define _ELEMENTS_H_

//
// FEMNet Elements Library
// 統合ヘッダーファイル - 全ての要素クラスを一括include
//

// 基底要素クラス
#include "ElementBase.h"

// 棒要素系（1次元要素）
#include "BarElement.h"        // 棒要素基底クラス
#include "TrussElement.h"      // トラス要素
#include "BeamElement.h"       // 梁要素
#include "ComplexBeamElement.h" // 複合梁要素

// 平面要素系（2次元要素）
#include "PlaneElement.h"      // 平面要素基底クラス

// 膜要素（平面応力要素）
#include "TriPlaneElement.h"   // 三角形膜要素
#include "QuadPlaneElement.h"  // 四角形膜要素

// 板要素（曲げ要素）
#include "TriPlateElement.h"   // 三角形板要素
#include "QuadPlateElement.h"  // 四角形板要素

#endif // _ELEMENTS_H_
