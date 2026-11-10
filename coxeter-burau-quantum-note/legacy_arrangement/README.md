# 超平面配置の代数的構成とモノドロミー計算

このディレクトリには研究ノートの一般配置に関する付録を再現する
Python 3 / SymPyコードを収録しています。
本文はリポジトリ直下の `coxeter_burau_quantum.pdf` と `.tex` です。

- `arrangement_burau.py`：階数1入力のconvolution、商、Burau行列。
- `general_convolution.py`：任意階数入力と最大降下商。
- `fox_burau.py`：Fox微分によるBurau行列。
- `monodromy_matrices.py`：純組紐と6直線例の基点付きモノドロミー。
- `ALGEBRAIC_CONSTRUCTION.md` / `MONODROMY.md`：詳しい規約と実行例。

リポジトリ直下で `make verify-arrangements` を実行すると、各スクリプトの
既存検証が実行され、同じディレクトリのJSONレポートが更新されます。

## 構成の範囲

Adachi–Hiroe の代数的 middle convolution を階数1に適用すると、射影方向を
選んだ超平面配置に対して留数行列を直接書ける。等重みの場合には通常の
Burau 表現との対応を配置空間への引き戻しとして説明できる。

ただし出力の定義域は一般には元の配置 A でなく、その方向の閉包 B。
元の配置の基本群表現を求めるなら、追加した超平面の周りのモノドロミーが
自明になることが必要。ここを省いて「任意の元の配置で構成済み」とは言えない。

## 入力と一般式

座標を (x,y) とし、射影方向を y 軸とする。非平行な超平面を

    H_i: y - f_i(x) = 0,  i=1,...,n

と書く。f_i はアフィン一次式。平行な超平面は D: g_D(x)=0。
階数1の入力接続を

\[
\Omega=\sum_i a_i\,d\log(y-f_i(x))+\sum_D b_D\,d\log g_D(x)
\]

とする。重み a_i,b_D は任意のスカラーなので入力は平坦。
相異なる f_i,f_j の差が非定数なら、その零点超平面を追加して B を作る。
重複する追加超平面は同一視する。差が非零定数なら追加は不要。

畳み込みパラメータを lambda、a=(a_1,...,a_n)^T とすると、W=C^n 上で

\[
C_{H_i}=e_i a^T+\lambda E_{ii}.
\]

衝突対 i<j の寄与を

\[
P_{ij}=a_j(E_{ii}-E_{ij})+a_i(E_{jj}-E_{ji})
\]

とおくと、平行超平面 D の留数は

\[
C_D=b_DI_n+\sum_{\{i,j\}:\,\{f_i=f_j\}=D}P_{ij}.
\]

追加した D では入力 b_D=0。ここで「=D」は base の超平面を total space に
引き戻して比較する意味。三重点など複数対が同じ D に衝突する場合も、全対を足す。
以上が [1] の式 (3.2)–(3.4) の階数1への特殊化。

middle quotient は次の部分空間を割る：

\[
K=\operatorname{span}\{e_i:a_i=0\},\qquad
L_\lambda=\ker(\mathbf1a^T+\lambda I_n),\qquad
W_\lambda=W/(K+L_\lambda).
\]

lambda≠0 では、L_lambda は sum(a_i)+lambda=0 のとき対角線 C*1、
そうでなければ零。lambda=0 のとき L_0=ker(a^T)。
係数を特殊化してから商を計算すること。一般パラメータでの商を特殊値へ代入する
だけでは、特殊値で生じる階数低下を取り逃す。

## Burau との対応を明示する

ここでは全ての非平行重みを a_i=alpha、平行重みを b_D=0、
畳み込みパラメータも lambda=alpha とする。alpha は一般の値とする。
このとき K=L=0 で出力階数は n。

閉包の補空間上には具体的な写像

\[
F:M(B)\longrightarrow\operatorname{Conf}_{n+1}(\mathbb C),\quad
(x,y)\longmapsto(f_1(x),\ldots,f_n(x),y)
\]

がある。衝突超平面を除いたことで、右辺の点は必ず相異なる。

**留数レベルでの直接検証。** V=C^(n+1)/C*1 を
u_i=w_i-w_(n+1) で表示する。通常の A_n 反射接続の留数
T_ij=(e_i-e_j)(e_i-e_j)^T のこの商上での表示に対し、
S=I_n+J_n（J_n は全成分1）とおくと

\[
S(C_{H_i}/\alpha)S^{-1}=T_{i,n+1},\qquad
S(P_{ij}/\alpha)S^{-1}=T_{ij}.
\]

右辺の T はいずれも商座標で表示した作用。従って今回の接続は、定数基底変換後に
\(\alpha\sum T_{ij}d\log(z_i-z_j)\) の F による引き戻しとなる。
この行列恒等式は任意の n で上式から証明でき、コードでも n=2,3,4 を検算。

一般パラメータでは、同じ内容を n+1 個の穿孔を持つ直線の階数1ねじれホモロジー
で読める。全ての有限穿孔のモノドロミーを q=exp(2*pi*i*alpha) とした
H_1 の次元は n で、その組紐群作用は reduced Burau。
したがって今回のモノドロミーは reduced Burau of B_(n+1) の
P_(n+1) への制限を F_* で引き戻したもの（同型を除く）。
ねじれホモロジーと KZ 型接続の対応は [3] §3 を参照。
ここでは q≠1、q^(n+1)≠1 などを満たす一般パラメータで述べており、
共鳴値の同定まで一括して主張してはいない。

これは **B_(n+1) の reduced Burau なので階数 n**。
下の lambda=-sum(a_i) による別の階数低下とは区別する。
平行重み b_D を非零にすると対応する階数1の捩りが加わる。
重みを不均等にした場合は、対応するねじれホモロジーの Gassner 型の側になる。

また y を含む全空間でなく射影の base の群を考える場合には、
x -> (f_1(x),...,f_n(x)) から P_n への写像が得られる。
等重みの平行留数 sum P_ij は C^n の対角線を零にし、商 C^n/C*1 上では
階数 n-1 の Burau 型接続となる。どちらの群の表現かを最初に決める必要がある。

## 計算例：三角形をなす3直線

\[
\mathcal A=\{y=0,\ y=x,\ y=1-x\},\qquad (a_1,a_2,a_3)=(a,b,c).
\]

追加される直線は x=0, x=1, x=1/2。出力は一般に階数3で、元の3直線の留数は

\[
C_1=\begin{pmatrix}a+\lambda&b&c\\0&0&0\\0&0&0\end{pmatrix},\quad
C_2=\begin{pmatrix}0&0&0\\a&b+\lambda&c\\0&0&0\end{pmatrix},\quad
C_3=\begin{pmatrix}0&0&0\\0&0&0\\a&b&c+\lambda\end{pmatrix}.
\]

追加直線の留数はそれぞれ P_12, P_13, P_23。

**階数2の特殊化。** a,b,c と a+b+c が非零のもとで lambda=-a-b-c とし、
u=v_1-v_3, v=v_2-v_3 を商座標とする。このとき元の3直線の留数は

\[
B_1=\begin{pmatrix}-b-c&b\\0&0\end{pmatrix},\quad
B_2=\begin{pmatrix}0&0\\a&-a-c\end{pmatrix},\quad
B_3=\begin{pmatrix}-a&-b\\-a&-b\end{pmatrix}.
\]

追加直線の留数は

\[
B_{x=0}=\begin{pmatrix}b&-b\\-a&a\end{pmatrix},\quad
B_{x=1}=\begin{pmatrix}a+c&0\\a&0\end{pmatrix},\quad
B_{x=1/2}=\begin{pmatrix}0&b\\0&b+c\end{pmatrix}.
\]

それぞれの固有値は (0,a+b), (0,a+c), (0,b+c)。いずれかの和が非整数なら、
その追加直線の局所モノドロミーは非自明で、元の3直線配置には降りない。
特に等重みでは、2*alpha が非整数なら降下できない。
整数の場合に降下できると逆向きに結論することはしていない。

## 最小の Burau 確認例

元から閉じている配置 {y=0,y=x,x=0} に、入力重み (alpha,alpha,0)、
lambda=alpha を使う。出力は階数2で

\[
C_{y=0}=\alpha\begin{pmatrix}2&1\\0&0\end{pmatrix},\quad
C_{y=x}=\alpha\begin{pmatrix}0&0\\1&2\end{pmatrix},\quad
C_{x=0}=\alpha\begin{pmatrix}1&-1\\-1&1\end{pmatrix}.
\]

(0,x,y) への写像で P_3 と対応し、一般パラメータで B_3 の reduced Burau に戻る。
比較用にコードが使う B_3 の具体的な braid generator matrices は

\[
\beta_q(\sigma_1)=\begin{pmatrix}1-q&q\\1&0\end{pmatrix},\qquad
\beta_q(\sigma_2)=\begin{pmatrix}1&-1\\0&-q\end{pmatrix}.
\]

これらは留数とは別の基底。純組紐生成元なら例えば
A_12=sigma_1^2, A_23=sigma_2^2, A_13=sigma_2*sigma_1^2*sigma_2^(-1)
を評価する。一般配置の実際の基点付きループに対する行列が必要なら、まず F_* に
よるそのループの braid word を決めてから `burau_word` に渡す。
本コードは任意の配置からその braid word を自動復元する部分は含まない。

**留数行列 C_H を全て同時に exp(2*pi*i*C_H) にしても、同一基点での
生成元行列には一般にならない。** 基点までの平行移動が必要で、共鳴の場合には
局所指数化の扱いにも注意が要る。今回の主出力は接続／holonomy Lie algebra の行列。

## 実行方法

```bash
python3 -m pip install sympy
python3 arrangement_burau.py
```

同じディレクトリに `results.json` を書き出す。
別の配置の例：

```python
import sympy as s
from arrangement_burau import convolve, middle_quotient, kohno_check
x,y=s.symbols('x y')
a,b,c,lam=s.symbols('a b c lambda')
m=convolve([y, y-x, y+x-1], [a,b,c], (x,y), y, lam)
for f,C in zip(m.equations,m.residues):
    print('Hyperplane:',f,'=0'); print(C)
print(kohno_check(m.equations,m.residues,m.variables))
Q,R,B,N=middle_quotient(m,{lam:-a-b-c})
```

幾何の係数は厳密な有理数を入力する（1/2 は `s.Rational(1,2)`）。
重みと lambda は記号でもよい。複素係数や記号係数の幾何はこの実装の対象外。
一般方向には、あらかじめその方向を y 軸にする座標変換を行う。
大きな配置は手元で実行できる。留数生成は衝突対について O(n^2) 個の寄与を作るが、
全 codimension-two 関係の検算はそれより重いので、必要な場合だけ呼ぶ。

## 実行済みの検証

- 三角形例：一般パラメータの階数3、および特殊化後の階数2で、それぞれ
  全6個の codimension-two flat に付随する15個の可積分性関係を記号検算。
- 三重点と既存の平行超平面、差が定数の平行な直線、3次元の配置を検算。
- lambda=0、a=0 を含む商の階数変化、非平行超平面がない場合を検算。
- middle quotient の不変性 `Q*C=B*Q` を全留数について検算。
- A_n 反射接続との定数共役を n=2,3,4 で検算。
- B_3,B_4,B_5 の braid relations、離れた生成元の可換性、full twist の
  像 q^n*I を検算。全て PASS。

## 一次資料

1. S. Adachi–K. Hiroe, *Middle convolution for affine hyperplane arrangements
   and Aomoto realizability*, arXiv:2601.00544v4 (2026-09-29),
   §3、特に式 (3.2)–(3.6)、平坦性は Proposition 4.6。
   https://arxiv.org/abs/2601.00544v4
2. K. Hiroe, *Middle convolution for Lie algebra representations*,
   arXiv:2605.09828v3, §§3–4。
   https://arxiv.org/abs/2605.09828
3. T. Kohno, *Hyperplane arrangements, local system homology and iterated
   integrals*, Advanced Studies in Pure Mathematics 62 (2012), 157–174, §3。
   https://www.ms.u-tokyo.ac.jp/~kohno/papers/kohno_aspm.pdf

一般の留数公式・平坦性・middle quotient は既存構成の使用。
このパッケージは階数1への特殊化、Burau との対応の行列表現、具体例と検算をまとめたもの。
