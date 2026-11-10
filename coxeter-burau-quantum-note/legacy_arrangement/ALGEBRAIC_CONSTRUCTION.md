# 任意の超平面配置を入力とする代数的構成

前の README の階数1計算を、任意の階数へ拡張したもの。
配置は任意の有限アフィン超平面配置。理論は複素数上、実装の幾何は有理数係数。
結論は「任意の配置から、その方向の閉包上の表現を構成できる」。
元の配置上に非自明な Burau 型表現を常に作れるという主張ではない。

## 1. Lie 代数の表現としての完全に代数的な構成

任意の配置 A の holonomy Lie algebra は

\[
\mathfrak h(A)=\operatorname{Lie}\langle t_H:H\in A\rangle/
\left\langle[t_H,\sum_{H'\supset X}t_{H'}]:
X\in L_2(A),\ H\supset X\right\rangle.
\]

従って入力は有限次元 E と行列 A_H であり、上の括弧関係を満たすことだけを要求する。
方向 Y を選んで B=A^Y を作り、B\A の生成元には零を割り当てる。
これは全射 p:h(B)->h(A) による表現の引き戻しである。
廣惠の Theorem 4.7 は、この引き戻し後の Lie 代数版 MC と
Haraoka の行列構成の一致を与える [1]。

非平行超平面を H_i={y=f_i(x)}, i=1,...,n とする。
W=E^n、r=dim E とし、ブロックの順序を i、次に E の座標とする。

\[
C_{H_i}=\sum_{j=1}^n E_{ij}\otimes A_{H_j}
             +\lambda E_{ii}\otimes I_E.
\]

平行超平面 D では

\[
C_D=I_n\otimes A_D+
\sum_{\{i,j\}:\,\{f_i=f_j\}=D}
\left((E_{ii}-E_{ij})\otimes A_{H_j}
 +(E_{jj}-E_{ji})\otimes A_{H_i}\right).
\]

これらは h(B) の関係を満たす。次に

\[
K=\bigoplus_i e_i\otimes\ker A_{H_i},\qquad
L_\lambda=\ker\left(\sum_i C_{H_i}\right),\qquad
\operatorname{mc}^Y_\lambda(E)=W/(K+L_\lambda)
\]

とする。K,L は各 C_H で不変なので、行列を商へ落とせる。
lambda≠0 では

\[
L_\lambda=\{(v,\ldots,v):(\sum_i A_{H_i}+\lambda I_E)v=0\},
\qquad K\cap L_\lambda=0.
\]

よって階数は

\[
nr-\sum_i\dim\ker A_{H_i}
  -\dim\ker(\sum_i A_{H_i}+\lambda I_E).
\]

ここまで行列の加算、核、商だけで完了し、微分方程式を解く必要はない。
射は I_n tensor phi で誘導され、

\[
\operatorname{mc}^Y_\lambda:
\operatorname{Rep}_{fd}\mathfrak h(A)
\longrightarrow\operatorname{Rep}_{fd}\mathfrak h(A^Y)
\]

という関手になる。A=A^Y なら元の Lie 代数の表現圏の自己関手。
行列式と middle quotient は [2] §3 に対応する。

## 2. 元の配置へ降りる最大商：追加の代数的整理

これは上の MC そのものではなく、出力に任意に適用する普遍的な商操作。
V を出力の h(B)-加群、B_H をその作用行列とする。

\[
N_0=\sum_{D\in B\setminus A}\operatorname{im}B_D,\qquad
N_{k+1}=N_k+\sum_{H\in B}B_HN_k.
\]

次元が安定したところを N とすると、V/N は h(A)-加群であり、
V から得られる h(A)-加群の最大商である。
厳密な次元増加は高々 dim V 回なので、この操作も有限回の線形代数で終わる。

証明：h(B)->h(A) の核は追加 t_D が生成する Lie ideal。
N は追加作用の像を含む最小の h(B)-部分加群なので、V/N 上でこの ideal は零。
また、その ideal が零になる任意の商写像の核は N を含む。これが普遍性を与える。

**この最大商は零になり得る。** 特に階数1の等重み a_i=lambda=alpha≠0、
平行重み0の Burau 型について、実際に追加超平面 D が一つでもあれば、
この Lie 代数としての最大商は零になる。

証明：C_i=alpha*e_i*(1^T+e_i^T)。追加 C_D は衝突クラスごとの完全グラフの
Laplacian の alpha 倍であり、非零かつ像は ker(1^T) に含まれる。
その像の非零元 v を取り、v_i≠0 とすると C_i v=alpha*v_i*e_i。
従って N は e_i を含む。さらに j≠i で C_j e_i=alpha*e_j なので N=W。
この場合 lambda+sum(a_i)=(n+1)*alpha≠0 で middle quotient 自体は W のまま。

この観察は「追加作用を商で消せば、一般配置上の非自明な Burau が必ずできる」
という期待への障害である。他の入力、方向、特殊パラメータ、別構成の不存在までは言わない。

## 3. 群表現を積分せず作る場合：Fox 微分

等重み Burau 型なら、B=A^Y 上の写像

\[
F:M(B)\to\operatorname{Conf}_{n+1}(\mathbb C),\quad
(x,y)\mapsto(f_1(x),\ldots,f_n(x),y)
\]

による braid monodromy と Artin action を使える。
m=n+1、K=Q(q)、F_m=<x_1,...,x_m> とし、全 x_i を q に送る指標を取る。
ねじれ鎖複体と H_1 は

\[
K^m\xrightarrow{(q-1,\ldots,q-1)}K,
\qquad H_1=\ker(q-1,\ldots,q-1),\quad\dim_K H_1=m-1.
\]

ループが与える自由群自己同型 phi に対し、Fox Jacobian

\[
J(\phi)_{ij}
=\left.\frac{\partial\phi(x_j)}{\partial x_i}\right|_{x_1=\cdots=x_m=q}
\]

を H_1 に制限する。この作用が Burau のホモロジー模型。
自由群の語と Laurent 多項式の演算だけで群表現行列が出る。
`fox_burau.py` はこの計算を実装し、以前の quotient 模型との基底変換を
B_3,B_4,B_5 で厳密検算する。

ただし、任意の配置と任意の基点付きループから braid monodromy の語を自動計算する
実装は含めていない。Fox 計算への入力にはその語／Artin 自己同型が必要。
また q=1 ではねじれ H_1 の次元が変わるので、この Q(q) 上の模型をそのまま
「q=1 のねじれ H_1」と同一視しない。Burau 行列の特殊化は別途定義できる。

元の群 G_A=pi_1(M(A)) については G_B -> G_A の核は追加 meridian の正規閉包。
既知の群表現 rho:G_B->GL(V) から降りる最大商は

\[
V/\operatorname{span}\{\rho(g)(\rho(m_D)-I)v:
g\in G_B,D\in B\setminus A,v\in V\}.
\]

有限生成元の行列で像の和を不変閉包にする操作で計算できる。
`maximal_descent(...,mode='group')` がこれを行う。
群では rho(m_D)=I、Lie 代数では B_D=0 を要求するので、両判定は区別する。
特に共鳴値で「Lie 代数として降りない」だけから「群表現としても降りない」とは言えない。

## 4. 使用と検証

```bash
python3 general_convolution.py
python3 fox_burau.py
```

前者は `general_results.json` を生成する。

```python
import sympy as s
from general_convolution import matrix_convolve, matrix_middle, maximal_descent
x,y=s.symbols('x y')
A1=s.Matrix([[2,1],[0,0]])/7
A2=s.Matrix([[0,0],[1,2]])/7
A3=s.Matrix([[1,-1],[-1,1]])/7
m=matrix_convolve([y,y-x,x],[A1,A2,A3],(x,y),y,s.Rational(1,11))
Q,R,B,N=matrix_middle(m)
answer=maximal_descent(B,m.added)
print('output rank:',Q.rows,'descended rank:',answer['quotient_dimension'])
```

入力行列が holonomy 関係を満たさなければ、計算開始時の検査で止まる。
この例は非可換な階数2の入力で、閉じた方向では middle 後も階数2。
同じ3直線配置で全ての直線に横断的な方向を選ぶと、直線が1本追加され、
middle 後の階数は3、Lie 代数として元の配置に降りる最大商は零だった。
三角形・等重み例では N の次元が 2 -> 3 と増えて零商になる。
lambda=0 の別の特殊化では階数1が残り、追加留数は零となる。
以上に加え、階数1の旧コードとの一致、各例の全可積分性関係と商の不変性を検算済み。

## 一次資料

[1] K. Hiroe, *Middle convolution for Lie algebra representations*,
arXiv:2605.09828v1, Theorem 3.2, §3.3, Theorem 4.7.
https://arxiv.org/html/2605.09828v1

[2] S. Adachi–K. Hiroe, *Middle convolution for affine hyperplane arrangements
and Aomoto realizability*, arXiv:2601.00544v4, §3 and Proposition 4.6.
https://arxiv.org/abs/2601.00544v4

§1 は既存構成の具体化。§2 の最大商と等重みの零商の議論は、この計算を整理する
ためにここで示した基本的な代数的帰結であり、新規性を主張するものではない。
