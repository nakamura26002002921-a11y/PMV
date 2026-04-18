# pmv_calc

C++（CGAL）と pybind11 を用いた体積計算モジュール。

---

## Overview

`pmv_calc` は、C++で実装された体積計算関数を Python から呼び出せるようにした拡張モジュールです。
高精度な幾何計算のために CGAL を使用しています。

主に、分子動力学（MD）シミュレーションデータからの体積評価に利用することを想定しています。

---

## Environment Setup

このリポジトリには `environment.yml` が含まれています。
以下のコマンドで環境を構築してください：

```bash
conda env create -f environment.yml
conda activate pmv_calc_env
```

---

## Requirements

主な依存関係：

* CMake >= 3.14
* Python (3.x)
* pybind11
* CGAL
* Eigen3
* Boost
* GMP / MPFR

### Python packages

* numpy
* scipy
* MDAnalysis

---

## Build

```bash
mkdir build
cd build
cmake .. \
  -DCMAKE_PREFIX_PATH=$CONDA_PREFIX \
  -DPython3_EXECUTABLE=$CONDA_PREFIX/bin/python
make -j
```

ビルドが成功すると、`pmv_calc*.so` が生成されます。

## Build

```bash
mkdir build
cd build

cmake .. \
  -DCMAKE_PREFIX_PATH=$CONDA_PREFIX \
  -DPython3_EXECUTABLE=$CONDA_PREFIX/bin/python

make -j
```

ビルドが成功すると、`pmv_calc*.so` が `build/` ディレクトリに生成されます。

Pythonから利用するには、以下のいずれかを行ってください：

* `build/` ディレクトリで実行する
* または `PYTHONPATH` に追加する

例：

```bash
export PYTHONPATH=$PWD/build:$PYTHONPATH
```

---

## Usage

### Pythonからの基本的な利用

```python
import pmv_calc

result = pmv_calc.compute_volume(
    REGIONS,   # Voronoi領域インデックス
    VERTICES,  # Voronoi頂点座標
    POINTS,    # 原子座標
    RADIUS     # 半径
)
```

---

### MD trajectory を使った例

```python
import MDAnalysis as mda
import pmv_calc

u = mda.Universe("md.pdb", "centered.xtc")
protein = u.select_atoms("protein")

for ts in u.trajectory:
    coords = protein.positions
    # Voronoi分割などを行い pmv_calc.compute_volume に渡す
```

---

### 実行例

```bash
python python/main.py
```


## Preprocessing (Important)

正しい体積計算のためには、周期境界条件（PBC）の補正が必要です。
以下のように GROMACS の `trjconv` を用いて前処理してください：

```bash
echo 0 | gmx trjconv -f md.xtc -s md.tpr -o md_w.xtc -pbc whole
printf "1\n0\n" | gmx trjconv -f md_w.xtc -s md.tpr -o md_w_cluster.xtc -pbc cluster
printf "1\n0\n" | gmx trjconv -f md_w_cluster.xtc -s md.tpr -o md_w_cluster_center.xtc -pbc mol -ur compact -center
```

その後、`md_w_cluster_center.xtc` を入力として使用してください。

---

## Running Example

```bash
python python/main.py
```

---

## Project Structure

```
pmv_calc/
├── CMakeLists.txt
├── environment.yml
├── src/
│   └── module.cpp
├── python/
│   └── main.py
├── README.md
└── LICENSE
```

---

## Notes

* CGAL は GMP / MPFR に依存しています
* 環境によっては `CMAKE_PREFIX_PATH` の指定が必要です

例：

```bash
cmake -DCMAKE_PREFIX_PATH="/path/to/CGAL;/path/to/Eigen" ..
```

* `build/` ディレクトリはリポジトリに含めないでください（`.gitignore`で除外）

---

## License

This project is licensed under the terms of the LICENSE file.
