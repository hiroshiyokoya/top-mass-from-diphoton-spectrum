# コードレビューと修正提案（2026-10）

対象：フィット部分（2018 年版の Python。`fit/TMDP.py`, `fit/Scan*.py`）と、信号を計算する Fortran（2016 年、arXiv:1607.00990。`fortran/lib/`, `fortran/legacy/`）。
**このドキュメントは提案である。** 修正は承認されたものから別イシュー・別 PR で行う（イシュー #5）。

> リポジトリは 2026-10 に `pyTMDP` から `top-mass-from-diphoton-spectrum` に改名し、Fortran を主、フィットを従とする構成に再編した（イシュー #12）。Docker イメージ名は `tmdp`。本文中のパスは再編後のもの。2018 年当時の版を指すときは「2018 年版」と書く。

環境：`docker/Dockerfile`（ROOT 6.34.00 / Python 3.12 / gfortran 13 / LHAPDF 6.5.5 / CHAPLIN 1.2）。
各項目の確度：**verified**＝コンテナで実行して確かめた／**read**＝コードを読んで確認したが実行していない／**inference**＝推論。

## 優先度：高

### P1. 現行 ROOT では疑似データが空になる（`TMDP.genEvents`）— verified → **修正済み（#14）**
- `hGenBG.FillRandom(self.hBG, N)` と `hGenSig.FillRandom(self.hSig, N)` は、元ヒスト（200 ビン／1001 ビン）と先ヒスト（`hbin`=100 ビン）のビニングが違う。
- ROOT 6.34 では、ビニングが異なり `N` がおよそ `10 × nbins` を超えると **エラーも警告も出さずに 1 事象も詰めない**（最小例：200→100 ビン、N=1000 は詰まり、N=1999 以上は 0）。結果、空ヒストに対してフィットが走り、χ² は数値として返る。
- 2018 年当時の ROOT では動いていたと推定（inference）。どの版で変わったかは未確認。
- **修正（#14）**：`FillRandom` をやめ、背景と信号の形をフィット用ビンに積分した期待値 $\mu_i$（`bin_fractions`）から `numpy.random.Generator.poisson(μ_i)` でビンごとに生成する。乱数の種は fit 入力の `seed` で固定できる。P4 も同時に解決。
- 確認（verified）：`tests/fit/test_smoke.py::test_pseudo_data` で、総数・ビンごとの Poisson ゆらぎ（χ²）・信号の形を検査する。合成背景（2000 万事象）での擬似実験 30 回では、χ²/ndf = 1.01、$k_{gg}$ = 0.099（真値 0.1）だった。

### P2. 最良質量がテンプレート格子に量子化される（`fit/ScanMass.py`, `fit/TMDP.py` の `__main__`）— read
- 各疑似実験の最良値が「χ² 最小のテンプレート質量」なので、結果は格子間隔（0.5–1 GeV）の離散値になり、`np.std(list_best)` は格子より細かい精度を表せない。
- **提案**：χ²(m_t) の最小点近傍 3–5 点に放物線を当てて連続的な最小値と Δχ²=1 幅を出す。疑似実験ごとの統計誤差も同時に得られる。

### P3. Fortran の積分関数がカット外で戻り値を設定しない — verified → **修正済み（#15）**
- `fortran/legacy/` の積分関数 7 本（`DSig_gg2aa`, `DSig_gg2aa0`, `DSig_terms`, `MKD_gg2aa`, `Scl_gg2aa`, `Sig_gg2aa`, `Sig_qqb2aa`）の `INT2`/`INT1` は、`IF (cut) RETURN` で **関数値を代入せずに戻る**（未定義値）。
- 同じソースから `INT2 = 0D0` の 1 行だけを除いたものを、最適化レベルだけ変えて実行した結果（$\sqrt s$=100 TeV, CT14lo, $m_t$=171.6, $\Gamma_t$=1.5, $m_{\gamma\gamma}$=300 GeV, 同じ乱数列）：

  | ビルド | $d\sigma/dm_{\gamma\gamma}$ [fb/GeV] |
  |---|---|
  | `-O0`（初期化なし） | 11.9744 |
  | `-O2`（初期化なし） | 5.80641 |
  | `-O2`（`INT2 = 0D0` あり＝`fortran/src/mktemplate.f`） | 5.80641 |
  | 同梱 `fortran/legacy/MKD_gg2aa.f`（修正前）を `-O2` でビルド | 11.9744 |

  **結果がコンパイラの最適化次第で約 2.06 倍変わる。** カットで捨てるべき点に直前の値が返り、実質的にカット（特に $p_T>0.4\,m_{\gamma\gamma}$）が効いていない状態になる。
  裏づけ：同じ条件で $p_T>0.4\,m_{\gamma\gamma}$ カットだけを外す（`ptratio: 0`）と 11.8767 になり、`-O0` の 11.9744 に近い（差 0.8%）。
- 未定義動作の中身は「捨てた点で直前に受理した積分値を返す」こと。`mktemplate.f` の `LEGACYCUT=1` でこれを明示的に再現すると、上の `-O0` の結果、および修正前の `MKD_gg2aa.exe` と全桁一致する（verified）。
- **1607.00990 の図は影響を受けていない**（verified）。論文の図から取り出した曲線（`reference/1607.00990/`）は、正しいカット処理（`LEGACYCUT=0`）で 0.5% 以内に再現される。FCC（$p_T>0.4\,m_{\gamma\gamma}$）では、旧挙動だと 2.05 倍ずれる。詳しくは `config/repro/1607.00990.yml`、`scripts/reproduce_1607_00990.py`、`docs/repro-1607.00990/`。
- 2018 年版のテンプレート（当時の `Template/`）がどちらで作られたかは、ファイルが無いので確認できない。
- **修正（#15）**：`fortran/legacy/` の 7 本の積分関数の先頭に `INT2 = 0D0` を 1 行ずつ足した（それ以外は無変更）。修正後の `MKD_gg2aa.exe` は、`-O0` でも `-O2` でも 5.80641 を返し、`mktemplate`（`LEGACYCUT=0`）と一致する（verified）。元の挙動は git 履歴と `mktemplate` の `LEGACYCUT=1` で再現できる。

## 優先度：中

### P4. 疑似実験の総事象数が固定（`genEvents`）— read → **修正済み（#14、P1 と同時）**
- `FillRandom(h, int(N))` は総数を固定するので、総数と信号比の Poisson 揺らぎが入らない。精度をやや過小評価する。P1 の修正で、ビンごとの Poisson 生成になった。

### P5. テンプレート関数が格子決め打ちの「最近傍（左側）参照」（`GG2AA.interpolate_func`）— verified（格子）／read（その他）
- `i = int(10*(x-299.95))` は 300–400 GeV・0.1 GeV 刻み・1001 点を前提にしている。別の格子だと `None` が返り、`TF1` 評価で `TypeError` になる（スモークテストで確認）。
- `x < 299.95` では `int()` が 0 方向に丸めるため負のインデックスになり、**Python の負インデックスで配列の末尾（≈399 GeV）を黙って返す**。現行の `fmin=330` では起きない。
- 補間ではなく階段関数で、`TF1` はビン中心の値×ビン幅で期待値を近似している（`I` オプションなし）。dip–bump 構造が 1 GeV ビン内で変わるので、ビン積分した期待値を使うべき。
- **提案**：テンプレートを `numpy` 配列で持ち、`np.interp`／ビン境界での累積積分から期待値を作る。格子はファイルから読む。

### P6. ファイル名からの $m_t,\Gamma_t$ 取得がパス全体を `_` で割る（`get_mass_width_from_filename`）— verified → **修正済み（#16）**
- `fname.split('_')` にディレクトリを含むフルパスを渡しているので、ディレクトリ名に `_` があると壊れる（pytest の一時ディレクトリで `ValueError` を確認）。
- **修正（#16）**：`os.path.basename` を取ってから正規表現 `Tab_<mt>_<Gt>.dat` で読む。形式が違えば `ValueError`。テストは `_` を含む pytest の一時ディレクトリをそのまま使う。

### P7. フィットの収束状態を見ていない — read
- `hGen.Fit(...)` の戻り値（`TFitResultPtr`／status）を確認せずに χ² を比較している。失敗したフィットも最良候補になりうる。
- **提案**：`S` オプションで結果を受け、`IsValid()` でなければ除外・記録する。

### P8. `GG2AAG` の Green 関数キャッシュのキーがエネルギーだけ（`fortran/lib/amp/gg2aaG.f`）— read
- `E = RS - 2MT` が前回と同じなら、`MT`・`GT`・`MU`・`J` が変わっても前回の `GRN` を返す。現在の使い方（1 本の積分中は全部固定）では問題にならないが、幅を変えて同じ $m_{\gamma\gamma}$ を続けて呼ぶ使い方では古い値を返す。
- **提案**：キャッシュのキーを (E, MT, GT, MU, J) にする。

### P9. `GRNNLOMSB` の数値的な不備（`fortran/lib/green/GrnMSBNLO.f`）— read
- 収束判定の比較値 `BP` が 1 回目のループで未初期化。
- 収束しないと `NLOOP` と `X1`（DATA 文／COMMON）を増やしてやり直すが、増えた値が以後の呼び出しにも残る（呼び出し順で結果と計算時間が変わりうる）。
- `PI=3.14159254D0`（正しくは …265）。相対 $3\times10^{-8}$ で実害はない。
- **提案**：`BP` を初期化し、`NLOOP`・`X1` は局所変数にコピーしてから増やす。Green 関数単体の検証（$\alpha_s\to0$ で $G_0$ に一致、LO ポテンシャルで 1S 極が $E=-m_t(C_F\alpha_s)^2/4$）をテストにする。

### P10. VEGAS の乱数が点ごとに続きから使われる — read／inference
- `fortran/extern/intvegas.f` の `RANF` は種が固定（234612947）で、プロセス内では状態が続く。旧 `fortran/legacy/MKD_gg2aa.f` は 1 プロセスで 8 質量を回すので、質量ごとに別の乱数列になる。
- 新しい `scripts/make_templates.py` は 1 質量 1 プロセスなので、全テンプレートが同じ乱数列を使う（質量間の差では統計揺らぎが相関して打ち消し合う方向。inference）。
- 2 次元の滑らかな積分なので、**決定論的な求積（Gauss–Legendre）に替えればテンプレートのジッターが消える**。提案。

## 優先度：低

- **P11.** `fortran/legacy/MKD_gg2aa.f` の主プログラムが `INT2` を `DOUBLE COMPLEX` と宣言しているが、関数は `DOUBLE PRECISION`（read。主プログラムからは直接呼ばないので実害なし。`-std=legacy` でコンパイル可）。
- **P12.** 単位換算 `389429.57D6`（GeV⁻² → fb）。$(\hbar c)^2 = 0.3893794$ GeV² mb なら `389379.4D6` で、相対 $+1.3\times10^{-4}$（inference。出典を確認）。
- **P13.** ROOT オブジェクト名の重複（`'dir'`, `'sig'`, 空名の `TF1`）。複数の `TMDP` を作ると上書き・警告になる（read）。
- **P14.** `set_init(self, param, dict)` が組み込みの `dict` を隠す。yml の `Nloop` は読まれていない。docstring のファイル名（`Scan_mass.py` など）が実ファイルと違う。`i % (Nloop/10)` が float 演算。`plt.show()` はヘッドレス環境で無意味（read）。
- **P15.** `calc_simpson_integral` は点数が奇数・等間隔を前提にし、偶数だと最後の区間を黙って落とす（read）。
- **P16.** `fortran/legacy/Sig_tot.f` は BASES で積分しており、BASES を同梱していないのでビルドしていない。
- **P17.** ハード係数の $\mathcal{O}(\alpha_s)$ 定数 `AT10` は 0 固定（"not known yet"）。Green 関数の $\alpha_s$ スケールは 40 GeV 固定、PDF は 2018 年版の設定（`config/templates/LHC13T.yml` など）では CT14lo、論文の設定（`paper1607_*`）では CT14nlo（物理の仕様。変えるなら別途議論）。
- **P18.** 依存が PyROOT に集中している。P1・P5 の修正後は numpy + iminuit（または scipy）だけでフィットまで書け、ROOT なしで動かせる（大きめの変更。任意）。
- **P19.** $\alpha_s$ を走らせる次数は libQCD の `NQCD`（`fortran/legacy/MKD_gg2aa.f` は 0 = LO、`Sig_gg2aa.f` は 1 = NLO）で、ドライバごとに違う。PDF は CT14nlo でも、1607.00990 の図は LO の走り（`NQCD=0`、$\alpha_s(M_Z)=0.1185$）で計算されている。NLO の走りにすると、全 $m_{\gamma\gamma}$ で 1.4% 小さくなる（verified）。PDF と $\alpha_s$ の次数をそろえるなら、物理の仕様として別途議論する。
- **P20.** 擬似データの背景の形は、背景イベントのヒストグラム（`hBG`、200 ビン）そのものなので、背景サンプルの統計ゆらぎが毎回同じ形で擬似データに入る。1 ビンあたり 1000 事象（合成サンプル 20 万事象）では、形のゆらぎ ±3% が擬似データの Poisson ゆらぎ（1 GeV ビンで約 1.5%）を上回り、χ²/ndf ≈ 2.5 になった。2000 万事象では χ²/ndf = 1.01（verified）。背景サンプルは、擬似データ 1 回分よりも十分大きい統計で用意するか、背景の形を滑らかな関数に当ててから使う。

## 範囲外（文書化のみ）

- 背景イベント（`Direct.dat`, `OneF.dat`, `TwoF.dat`）は Fortran では作れず、リポにも無い。形式は「$m_{\gamma\gamma}$ の値を空白区切りで並べたテキスト」（`fit/TMDP.py` の `read_events_from_file`）。出所（DIPHOX など）は未確認。
