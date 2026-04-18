# pmv_calc

C++（CGAL）と pybind11 を用いた体積計算モジュール。

## Overview

`pmv_calc` は、C++で実装された体積計算関数を Python から呼び出せるようにした拡張モジュールです。
高精度な幾何計算のために CGAL を使用しています。

---

## Requirements

以下の環境が必要です：

* CMake >= 3.14
* Python (3.x)
* pybind11
* CGAL
* Eigen3
* Boost
* GMP / MPFR

---

## Build

```bash
mkdir build
cd build
cmake ..
make -j
```

ビルドが成功すると、`pmv_calc*.so` が生成されます。

---

## Usage

Python から以下のように利用できます：

```python
import pmv_calc

result = pmv_calc.compute_volume(
    REGIONS,
    VERTICES,
    POINTS,
    RADIUS
)
```

---

## Project Structure

```
pmv_calc/
├── CMakeLists.txt
├── src/
│   └── module.cpp
└── README.md
```

---

## Notes

* CGAL は GMP / MPFR に依存しています
* 環境によっては `CMAKE_PREFIX_PATH` の指定が必要です

例：

```bash
cmake -DCMAKE_PREFIX_PATH="/path/to/CGAL;/path/to/Eigen" ..
```

---

## License

This project is licensed under the terms of the LICENSE file.
