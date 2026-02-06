# FEMNet

有限要素法(FEM)による構造解析ライブラリ。C++で実装され、SWIGを通じてPythonとC#から利用可能。

## 機能

### 対応要素
- **Truss要素** - 軸力のみを伝達する2節点要素
- **Beam要素** - 曲げ・せん断・ねじりを考慮した2節点梁要素
- **ComplexBeam要素** - 端部剛性を考慮した梁要素
- **TriPlane要素** - 三角形平面応力要素
- **QuadPlane要素** - 四角形平面応力要素
- **TriPlate要素** - 三角形板要素（曲げ）
- **QuadPlate要素** - 四角形板要素（曲げ）

### 解析機能
- **線形静的解析** (`FELinearStaticOp`)
- **固有値解析** (`FEVibrateResult`)
- **座屈解析** (`FEBucklingAnalysis`)
- **動的解析** (`FEDynamic`)
- **応答スペクトル解析** (`ResponseSpectrumMethod`)

## 依存関係

- CMake 3.25+
- C++17対応コンパイラ
- Eigen3 3.4+
- Spectra
- SWIG 4.0+
- Intel MKL (オプション、高速化)

## ビルド

### C++ライブラリ（サンプルアプリ）

```bash
mkdir build && cd build
cmake .. 
cmake --build .
```

### C#バインディング

```bash
mkdir build-csharp && cd build-csharp
cmake .. -DBUILD_CSHARP=ON
cmake --build . {--config Release}
```

生成されるファイル:
- `fem.dll` - ネイティブライブラリ
- `*.cs` - C#ラッパークラス

### Pythonバインディング

```bash
mkdir build-python && cd build-python
cmake .. -DBUILD_PYTHON=ON -DBUILD_CSHARP=OFF
cmake --build .
```

ビルド後、`python/femnet/`ディレクトリに以下が生成されます:
- `_femnet.pyd` (Windows) / `_femnet.so` (Linux/Mac) - ネイティブ拡張
- `_femnet.py` - SWIGが生成したPythonラッパー

## Pythonパッケージのインストール

### 開発モード（editable install）

```bash
cd python
pip install -e .
```

ソースコードを編集しながら開発する場合に便利です。

### 通常インストール

```bash
cd python
pip install .
```

## 使い方

### Python

```python
from femnet import *

# モデル作成
model = FEModel()

# 節点追加
model.AddNode(0, 0.0, 0.0, 0.0)
model.AddNode(1, 1000.0, 0.0, 0.0)
model.AddNode(2, 2000.0, 0.0, 0.0)

# 境界条件（固定）
model.GetNode(0).Fix.FixAll()

# 材料追加（E=205GPa, ν=0.3）
model.AddMaterial(205e3, 0.3)

# 断面追加（A, Iy, Iz, K）
model.AddSection(1000.0, 8333.0, 833333.0, 100.0)

# 梁要素追加
model.add_beam_element(0, 0, 1, 0, 0, 0.0)  # id, n1, n2, sec_id, mat_id, beta
model.add_beam_element(1, 1, 2, 0, 0, 0.0)

# 荷重設定
loads = VectorLoad()
loads.append(NodeLoad(2, 0, -1000, 0))  # node_id, Px, Py, Pz

# 線形静的解析
solver = FELinearStaticOp(model, loads)
solver.Compute()

if solver.Computed():
    # 変位取得
    displacements = solver.GetDisplacements()
    for i, disp in enumerate(displacements):
        print(f"Node {i}: Dx={disp.Dx():.6f}, Dy={disp.Dy():.6f}")

    # 反力取得
    reactions = solver.GetReactionData()
    for react in reactions:
        print(f"Reaction at node {react.id}: Py={react.Py():.2f} N")

    # 要素応力取得
    stress = solver.GetBeamStress(0, 0.5)  # 要素ID, 位置(0-1)
    print(f"Bending moment Mz: {stress.Mz:.2f} N-mm")
```

### C#

```csharp
using FEMNet;

// モデル作成
var model = new FEModel();

// 節点追加
model.AddNode(0, 0.0, 0.0, 0.0);
model.AddNode(1, 1000.0, 0.0, 0.0);
model.AddNode(2, 2000.0, 0.0, 0.0);

// 境界条件
model.GetNode(0).Fix.FixAll();

// 材料・断面追加
model.AddMaterial(205e3, 0.3);
model.AddSection(1000.0, 8333.0, 833333.0, 100.0);

// 梁要素追加
model.add_beam_element(0, 0, 1, 0, 0, 0.0);
model.add_beam_element(1, 1, 2, 0, 0, 0.0);

// 荷重設定
var loads = new VectorLoad();
loads.Add(new NodeLoad(2, 0, -1000, 0));

// 解析実行
var solver = new FELinearStaticOp(model, loads);
solver.Compute();

if (solver.Computed())
{
    var displacements = solver.GetDisplacements();
    foreach (var (disp, i) in displacements.Select((d, i) => (d, i)))
    {
        Console.WriteLine($"Node {i}: Dy = {disp.Dy():F6} mm");
    }
}
```

