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
  -DPython3_EXECUTABLE=$CONDA_PREFIX/bin/python \
  -DCMAKE_BUILD_TYPE=Release

make -j
```

ビルドが成功すると、`pmv_calc*.so` が `build/` ディレクトリに生成されます。

---

## 実行方法

### 入力データ

以下の2つが必要：

* `md.pdb`
* `centered.xtc`


### 1. Monte Carlo 法（低速・近似）

```bash
export PYTHONPATH=$PWD/build:$PYTHONPATH
python python/monte.py
```


### 2. Voronoi 法（高速・高精度）

```bash
export PYTHONPATH=$PWD/build:$PYTHONPATH
python python/voronoi.py
```


## 出力ファイル

どちらも `result.txt` を生成

---

### Monte Carlo の出力

```
frame   time   volume   variance   water_count
```

#### 各列の意味

* **frame**: フレーム番号
* **time**: 計算時間（秒）
* **volume**: 推定体積(ファンデルワールス体積)
* **variance**: Monte Carlo のばらつき
* **water_count**: タンパク質近傍の水分子数

---

### Voronoi の出力

```
frame   time   volume   water_count
```

#### 各列の意味

* **frame**: フレーム番号
* **time**: 計算時間（秒）
* **volume**: 厳密体積
* **water_count**: タンパク質近傍の水分子数



## Project Structure

```
pmv_calc/
├── src/
│   └── module.cpp
├── python/
│   ├── monte.py
│   └── voronoi.py
├── example/
│   ├── md.pdb
│   └── centered.xtc
├── CMakeLists.txt
├── environment.yml
└── README.md
```

## 補足

### PBC補正

正しい体積計算のためには、周期境界条件（PBC）の補正が必要です。
以下のように GROMACS の `trjconv` を用いて前処理してください：

```bash
echo 0 | gmx trjconv -f md.xtc -s md.tpr -o md_w.xtc -pbc whole
printf "1\n0\n" | gmx trjconv -f md_w.xtc -s md.tpr -o md_w_cluster.xtc -pbc cluster
printf "1\n0\n" | gmx trjconv -f md_w_cluster.xtc -s md.tpr -o centered.xtc -pbc mol -ur compact -center
```

## License

This project is licensed under the terms of the LICENSE file.
