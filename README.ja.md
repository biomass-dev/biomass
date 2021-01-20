# BioMASS

[![Actions Status](https://github.com/okadalabipr/biomass/workflows/Tests/badge.svg)](https://github.com/okadalabipr/biomass/actions)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/okadalabipr/biomass.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/okadalabipr/biomass/context:python)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Downloads](https://pepy.tech/badge/biomass)](https://pepy.tech/project/biomass)
[![PyPI version](https://img.shields.io/pypi/v/biomass.svg?logo=PyPI&color=blue)](https://pypi.python.org/pypi/biomass/)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/biomass.svg)](https://pypi.python.org/pypi/biomass/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

![logo2](https://github.com/okadalabipr/biomass_docs.jl/blob/master/docs/src/assets/logo2.png?raw=true)

当ソフトウェアは，細胞内シグナル伝達機構の数理モデリングを目的としたツールです．実験から得られた実測データに基づいてモデル内のパラメータを最適化し，その結果に基づいた解析により，細胞応答に対する反応ネットワーク内の重要な要素を同定することができます．

ここでは，早期転写の反応ネットワーク ([Nakakuki _et al._, **_Cell_**, 2010](https://doi.org/10.1016/j.cell.2010.03.054)) を例に使用しています.

## インストール

```
$ pip3 install biomass
```

## 使い方

#### 実行可能なモデルの生成

```python
from biomass.models import Nakakuki_Cell_2010

Nakakuki_Cell_2010.show_info()
```

```
Nakakuki_Cell_2010 information
------------------------------
36 species
115 parameters, of which 75 to be estimated
```

```python
model = Nakakuki_Cell_2010.create()
```

#### ODE モデルのパラメータ推定 (n = 1, 2, 3, · · ·)

遺伝的アルゴリズムを使用して，目的関数（ここではモデルのシミュレーション値と実験データの誤差）を最小化します．

```python
from biomass import optimize

optimize(
    model=model, start=1, options={
        "popsize": 3,
        "max_generation": 1000,
        "allowable_error": 0.5,
        "local_search_method": "DE",
    }
)
```

世代交代毎に，最良のパラメータセットが`out/n/`に保存されます．

進捗は`out/n/optimization.log`で見ることができます．

```
Generation1: Best Fitness = 1.726069e+00
Generation2: Best Fitness = 1.726069e+00
Generation3: Best Fitness = 1.726069e+00
Generation4: Best Fitness = 1.645414e+00
Generation5: Best Fitness = 1.645414e+00
Generation6: Best Fitness = 1.645414e+00
Generation7: Best Fitness = 1.645414e+00
Generation8: Best Fitness = 1.645414e+00
Generation9: Best Fitness = 1.645414e+00
Generation10: Best Fitness = 1.645414e+00
Generation11: Best Fitness = 1.645414e+00
Generation12: Best Fitness = 1.645414e+00
Generation13: Best Fitness = 1.645414e+00
Generation14: Best Fitness = 1.645414e+00
Generation15: Best Fitness = 1.645414e+00
Generation16: Best Fitness = 1.249036e+00
Generation17: Best Fitness = 1.171606e+00
Generation18: Best Fitness = 1.171606e+00
Generation19: Best Fitness = 1.171606e+00
Generation20: Best Fitness = 1.171606e+00
```

- 途中で中断したところから再開したい場合，

```python
from biomass import optimize_continue

optimize_continue(
    model=model, start=1, options={
        "popsize": 3,
        "max_generation": 1000,
        "allowable_error": 0.5,
        "local_search_method": "DE",
    }
)
```

- 複数のパラメータセット（例：1 から 10）を同時に探索したい場合,

```python
from biomass import optimize

optimize(
    model=model, start=1, end=10, options={
        "popsize": 5,
        "max_generation": 2000,
        "allowable_error": 0.5,
        "local_search_method": "mutation",
        "n_children": 50
    }
)
```

- 最適化したパラメータの取得

```python
from biomass.result import OptimizationResults

res = OptimizationResults(model)
res.to_csv()
```

---

#### シミュレーション結果の可視化

パラメータ推定で得た複数のパラメータセットでのシミュレーション結果を出力します．結果は`figure/`に保存されます．

```python
from biomass import run_simulation

run_simulation(model, viz_type='average', show_all=False, stdev=True)
```

関数`run_simulation`の引数を設定することで，出力されるグラフの表示法を変更することができます．

**viz_type** : _str_

- `'average'` : `out/`にある複数のパラメータセットでのシミュレーション結果の平均を表示します．

- `'best'` : `out/`にある複数のパラメータセットでのシミュレーション結果のうち，最良のものを表示します．

- `'original'` : `set_model.py`に記述されているパラメータ，初期値を使用したシミュレーション結果を表示します．

- `'n(=1,2,...)'` : `out/n(=1,2,...)` における最新のパラメータセットでのシミュレーション結果を表示します．

- `'experiment'` : `observable.py` に記述されている実験値の結果のみを表示します．

**show_all** : _bool_ (default: False)

- `out/n(=1,2,...)`に格納されたパラメータセットでの全てのシミュレーション結果を淡色で表示します．

**stdev** : _bool_ (default: False)

- `viz_type == 'average'`の際，標準偏差も含めてシミュレーション結果を表示します．

**save_format** : _str_ (default: "pdf")

- 保存する図の拡張子（`"pdf"`または`"png"`）を指定します．

![simulation_best](https://github.com/okadalabipr/biomass_docs.jl/blob/master/docs/src/assets/simulation_best.png?raw=true)

点（青, EGF; 赤, HRG）は実験データ，線はシミュレーション結果を表す

---

#### 感度解析

```python
from biomass import run_analysis

run_analysis(model, target='reaction', metric='integral', style='barplot')
```

感度係数は以下の式で記述されます．

_s<sub>i</sub>_(_q_(**v**),_v<sub>i</sub>_) = _∂_ ln(_q_(**v**)) / _∂_ ln(_v<sub>i</sub>_) = _∂_ _q_(**v**) / _∂_ _v<sub>i</sub>_ · _v<sub>i</sub>_ / _q_(**v**)

ここで _v<sub>i</sub>_ は*i*番目の反応速度を表し, **v** は反応速度のベクトル **v** = (_v<sub>1</sub>_, _v<sub>2</sub>_, ...)，_q_(**v**) は出力を定量する関数です（例：応答の積分値，最大値，持続時間など）． 感度係数は微分を 1%の反応速度の変化で有限差分近似して計算されます．

各反応における感度係数を求めるためには，[`model/set_model.py`](biomass/model/set_model.py)中で，反応速度を 'v' で表す場合，全ての反応式を記述した直後に，以下を書いておく必要があります．

```python
if self.perturbation:
    for i, dv in self.perturbation.items():
        v[i] = v[i] * dv
```

**target** : _str_

何に対する感度解析かを選択します（反応速度・初期値・パラメータ）.

- `'reaction'`
- `'initial_condition'`
- `'parameter'`

**metric** : _str_ (default: 'integral')

出力に用いる基準を設定します．

- `'maximum'`
  : 最大値．

- `'minimum'`
  : 最小値．

- `'argmax'`
  : 最大値に到達するまでの時間．

- `'argmin'`
  : 最小値に到達するまでの時間．

- `'timepoint'`
  : options['timepoint']で指定した時間点におけるシミュレーション値.

- `'duration'`
  : options['duration']で指定した閾値まで減少するまでにかかる時間．

- `'integral'`
  : シミュレーション時間内における濃度の積分値．

**style** : _str_ (default: 'barplot')

グラフを選択します．

- `'barplot'`
- `'heatmap'`

**options** : _dict, optional_

詳細な metric を設定します．

- **timepoint** : _int_ (default: model.sim.t[-1])

  - (metric == `'timepoint'`) どの時間点を使用するかを指定します．

- **duration**: _float_ (default: 0.5)

  - (metric == `'duration'`) 0 から 1 の間の実数を指定します．例えば，最大値の 10%まで減少する時間を metric に使用する場合には 0.1 に設定します．

![sensitivity_PcFos](https://github.com/okadalabipr/biomass_docs.jl/blob/master/docs/src/assets/sensitivity_PcFos.png?raw=true)

pc-Fos の積分値に対する感度係数（青, EGF; 赤, HRG）．棒の上下の数字は反応番号を，エラーバーは標準偏差を表す．

## フィードバックのお願い

ソフトウェアの継続的改善のため，使用感の報告やバグレポート，アドバイスをいただきたいと考えております．また，お手持ちの実験データの背後に潜む制御メカニズムを数理モデルを用いて明らかにしたいというご要望も承りますので，どうぞお気軽にご相談ください．

## 免責事項

当ソフトウェアの出力結果の正確性や妥当性につきまして，一切の保障はいたしません．また当ソフトウェアを用いることで生じたあらゆる損害について，一切の責任を負いません．
