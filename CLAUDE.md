# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## プロジェクト概要

FEMNetは有限要素法(FEM)による構造解析ライブラリで、C++17で実装されています。SWIGを使用してC#(.NET)とPython向けのバインディングを生成し、StructureMM Grasshopperプラグイン（C#）とスタンドアロンのPythonパッケージの両方で使用されます。

**重要**: このディレクトリは独立したgitリポジトリです。変更は親プロジェクト(StructureMM)とは別にコミットされます。

## ビルドコマンド

**ビルド環境**: 想定する主開発環境は Visual Studio 2026 (VS18) です。CMake で新規に構成する場合は `-G "Visual Studio 18 2026"` を指定してください。既存の `build/` ディレクトリがあれば、その CMakeCache.txt のジェネレータ設定をそのまま使い回してください。

**ビルドディレクトリの方針**: ビルドは原則 `build/` ディレクトリ単一で行います。既存の `build/` が現在の環境に合わない、CMake変数を変えて再構成したい等の理由で作り直しが必要な場合も、新しいビルドディレクトリを増やさず、**既存の `build/` をクリーン（削除して作り直す）してからリビルド**してください。ビルドディレクトリが増えるのを避けるためです（絶対ではありませんが、できる限りこの方針に従ってください）。

### C#バインディング（StructureMM用）

```bash
# 初回のみ
mkdir build-csharp
cd build-csharp
cmake .. -DBUILD_CSHARP=ON -DBUILD_PYTHON=OFF

# ビルド（変更後は毎回実行）
cmake --build .
```

生成物: `fem.dll`（ネイティブライブラリ）と`*.cs`ファイルが`../../FEMNet/`に出力されます。

### Pythonバインディング

```bash
# 初回のみ
mkdir build-python
cd build-python
cmake .. -DBUILD_CSHARP=OFF -DBUILD_PYTHON=ON

# ビルド
cmake --build .

# Pythonパッケージのインストール（開発モード）
cd ../python
pip install -e .
```

生成物: `_femnet.pyd` (Windows) / `_femnet.so` (Unix)と`femnet.py`が`python/femnet/`に出力されます。

### C++サンプルアプリ

```bash
mkdir build
cd build
cmake ..
cmake --build .
```

生成物: `sampleapp.exe` - C++ライブラリの動作確認用

## アーキテクチャ

### コア構造

**ソースツリー構成** (`src/` 配下にドメイン別で配置):
- `src/core/` - `FEModel`([Model.h](src/core/Model.h)), `Components`, `RigidLink`, `LoadComponent`
- `src/Elements/` - 全要素クラス
- `src/analysis/` - 解析オペレーター群と`SeismicModule`
- `src/io/` - モデル・荷重のファイル入出力（`ModelIO`, `LoadIO`）
- `src/math/` - 疎行列ユーティリティ・ソルバー（`SparseMatrixUtils`, `SparseSolver`）
- `src/cuda/` - GPUソルバー（オプショナル）
- `apps/` - サンプルアプリ（[apps/app.cpp](apps/app.cpp)）

**モデルクラス階層**:
- `FEModel` ([src/core/Model.h](src/core/Model.h)) - メインモデルコンテナ（Nodes, Materials, Sections, Elements保持）
- `ElementBase` ([src/Elements/ElementBase.h](src/Elements/ElementBase.h)) - 全要素の基底クラス
  - `TrussElement`, `BeamElement`, `ComplexBeamElement` - 1次元要素
  - `TriPlaneElement`, `QuadPlaneElement` - 平面応力要素
  - `TriPlateElement`, `QuadPlateElement` - 板要素

**解析オペレーター**:
- `FELinearStaticOp` ([src/analysis/FELinearStaticOp.h](src/analysis/FELinearStaticOp.h)) - 線形静的解析
- `FEVibrationAnalysis` ([src/analysis/FEVibrationAnalysis.h](src/analysis/FEVibrationAnalysis.h)) - 固有値解析
- `FEBucklingAnalysis` ([src/analysis/FEBucklingAnalysis.h](src/analysis/FEBucklingAnalysis.h)) - 座屈解析
- `FEDynamic` ([src/analysis/FEDynamic.h](src/analysis/FEDynamic.h)) - 動的解析
- `ResponseSpectrumMethod` ([src/analysis/ResponseSpectrumMethod.h](src/analysis/ResponseSpectrumMethod.h)) - 応答スペクトル法

**基本コンポーネント** ([src/core/Components.h](src/core/Components.h)):
- `Node`, `Material`, `Section`, `Displacement`, `Point`, `Vector`など

### SWIGインターフェース構造

SWIGバインディングは3層構造:

1. **共通定義** (`swig/common/*.i`) - 言語非依存の定義
   - [fem_common.i](swig/common/fem_common.i) - メインエントリーポイント、STLテンプレート定義、クラス拡張
   - [Elements_common.i](swig/common/Elements_common.i) - 要素クラスの共通定義
   - [Operator_common.i](swig/common/Operator_common.i) - 解析オペレーターの共通定義
   - [SeismicModule_common.i](swig/common/SeismicModule_common.i) - 地震解析モジュール

2. **C#固有** (`swig/csharp/fem_csharp.i`)
   - C#の型マップ（typemap）、`ToString()`オーバーライド、GC対策

3. **Python固有** (`swig/python/fem_python.i`)
   - Python例外処理、`__str__`/`__repr__`メソッド、プロパティ定義

**重要**: 共通の変更は`swig/common/`で、言語固有の動作は各言語のディレクトリで行います。

### Eigenとの相互運用

- 内部計算にはEigen3を使用（疎行列、密行列、ベクトル演算）
- Eigenの型（`MatrixXd`, `VectorXd`等）はSWIGでラップ不可のため`%ignore`
- `FEModel`の拡張メソッド（`AddNode`, `AddMaterial`等）で簡易インターフェースを提供

## 開発ワークフロー

### C++コードの変更

1. `.h`/`.cpp`ファイルを修正
2. 新しいクラス/メソッドを公開する場合:
   - 適切な`swig/common/*.i`ファイルに`%include`または`%extend`を追加
   - STLコンテナが必要なら`fem_common.i`に`%template`を追加
3. CMakeで再ビルド
4. 生成された`.cs`または`.py`ファイルを確認

### SWIGインターフェースの変更

**新しいクラスを公開する場合**:
1. 対応する`.h`ファイルを`fem_common.i`または個別モジュールの`.i`ファイルに`%include`
2. 必要に応じて`%shared_ptr`宣言を追加
3. Eigenの型を使うメンバーは`%ignore`で除外
4. STLコンテナが必要なら`%template(VectorXxx) std::vector<Xxx>;`を追加

**言語固有の動作を追加する場合**:
- C#: `fem_csharp.i`で`%typemap(cscode)`を使用
- Python: `fem_python.i`で`%extend`内の`%pythoncode`を使用

### テスト

**C++**: [apps/app.cpp](apps/app.cpp)を編集してサンプルアプリで動作確認

**Python**: [python/test_femnet.py](python/test_femnet.py)でユニットテストまたは手動テスト

**C#**: 親プロジェクト(StructureMM)でGrasshopperコンポーネントを通じてテスト

## 依存関係

- **CMake** 3.25+ - ビルドシステム
- **C++17対応コンパイラ** - MSVC 2019+, GCC 9+, Clang 10+
- **Eigen3** 3.4+ - 線形代数演算（vcpkgで自動取得）
- **Spectra** - 固有値ソルバー（vcpkgで自動取得）
- **SWIG** 4.0+ - バインディング生成
- **Intel MKL** (オプション) - 高速化用、PardisoソルバーとBLAS/LAPACK

MKLが見つかった場合は自動的に有効化され、`EIGEN_USE_MKL_ALL`が定義されます。

## プロジェクト間の関係

このFEMNetライブラリは以下で使用されます:

1. **StructureMM** (親プロジェクト): C#バインディング(`fem.dll` + `*.cs`)を`../../FEMNet/`プロジェクトで参照
2. **Pythonパッケージ**: スタンドアロンの`femnet` Pythonパッケージ（`pip install`可能）

C#バインディングをビルドした後は、親プロジェクトのVisual Studioソリューションで`FEMNet`プロジェクトが生成されたファイルを認識します。
