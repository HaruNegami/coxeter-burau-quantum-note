# 基点を揃えたモノドロミー行列

2026-11-06。対象は a_i=lambda=alpha、入力の平行留数0の等重み Burau 型。
q=exp(2*pi*i*alpha) とする。接続／MC とこの Burau 模型の同定は一般パラメータで述べる。
Burau 行列自体は全ての q≠0 で定義できるが、共鳴値での MC の階数変化等は別途扱う。
以下は留数を個別に指数化した行列ではなく、同一基点・同一基底での群表現行列。

## 1. 一般の純組紐生成元の表示

m=n+1 本の組紐、reduced Burau の階数 n とする。
以前から使う quotient 模型 C^m/C*1 に座標 w_i-w_m を入れる。
unreduced の sigma_i は 2×2 ブロック [[1-q,q],[1,0]] を持つ規約。
列ベクトルに作用し、群の積では右の経路を先にたどる。
正の半捻りは2点の相対位置が反時計回りに半回転するものとする。

標準純組紐を

\[
A_{ij}=\sigma_{j-1}\cdots\sigma_{i+1}\sigma_i^2
\sigma_{i+1}^{-1}\cdots\sigma_{j-1}^{-1}
\]

とすると、この基底で次の閉形式が成り立つ：

\[
\boxed{\beta_q(A_{ij})=I_n+(q-1)u_{ij}v_{ij}^{\mathsf T}.}
\]

j≤n の場合、

\[
u_{ij}=q e_i+(q-1)\sum_{i<k<j}e_k-e_j,
\qquad v_{ij}=e_i-e_j.
\]

j=m=n+1 の場合、

\[
u_{i,m}=\sum_{k<i}e_k+(q+1)e_i+q\sum_{k>i}e_k,
\qquad v_{i,m}=e_i.
\]

証明：beta(sigma_i)^2-I を階数1に分解し、
beta(sigma_(j-1)...sigma_(i+1)) で共役する。
j≤n では -e_(i+1) が右へ移るたび、中間成分に q-1 が現れる。
j=m の最後の sigma_n を作用させると後者のベクトルになる。
この演算は各 n に共通であり、コードでも m=3,4,5,6 の全34生成元を検算した。

全ての i<j について

\[
v_{ij}^{\mathsf T}u_{ij}=q+1,\quad
\det\beta_q(A_{ij})=q^2,\quad
(\beta_q(A_{ij})-I)(\beta_q(A_{ij})-q^2I)=0.
\]

q≠±1 なら固有値1（n-1重）とq^2（1重）の複素擬鏡映。
q=-1 ではこの式だけから単位行列とは結論できず、非自明なunipotentの場合がある。

一般配置の方向閉包 B=A^Y に対して

\[
F(x,y)=(f_1(x),\ldots,f_n(x),y),\qquad
\rho(\gamma)=\beta_q(F_*\gamma).
\]

したがって、基点付きループの braid word を決めれば行列積で計算できる。
固定したxのファイバーでyがf_i(x)を回る標準meridianなら、上のA_(i,m)になる。
平行超平面の局所衝突が複数点なら、そのクラスのfull twistを評価する。
共通基点への接続路による共役を含める必要があり、一般に各衝突対のA_ijを
接続路なしに機械的に割り当てることはできない。

## 2. 三角形の3直線から作った6直線配置

\[
\mathcal B=\{y=0,y=x,y=1-x,x=0,x=1,x=\tfrac12\}.
\]

基点は (x,y)=(1/4,2)、順序付き点配置は (0,1/4,3/4,2)。
以下の各小円は反時計回り。

- y=0,y=x,y=1-x のmeridian：x=1/4を固定。y=2から下半平面を通って
  それぞれの穿孔のすぐ右へ進み、小円を一周して同じ接続路を戻る。
- x=0 のmeridian：y=2を固定、実軸に沿ってx=0のすぐ右へ進む。
- x=1/2 のmeridian：y=2を固定、実軸に沿ってx=1/2のすぐ左へ進む。
- x=1 のmeridian：y=2を固定、x=1/2を上半平面で迂回してx=1のすぐ左へ進む。

最後の3つも、それぞれ小円を一周して接続路を逆に戻る。
この規約で braid words は次の通り。

| 超平面 | 純組紐 | 群の語 |
|---|---|---|
| y=0 | A_14 | sigma_3 sigma_2 sigma_1^2 sigma_2^-1 sigma_3^-1 |
| y=x | A_24 | sigma_3 sigma_2^2 sigma_3^-1 |
| y=1-x | A_34 | sigma_3^2 |
| x=0 | A_12 | sigma_1^2 |
| x=1 | A_13 | sigma_2 sigma_1^2 sigma_2^-1 |
| x=1/2 | A_23 | sigma_2^2 |

小円半径1/20、迂回半径と下半平面の深さ1/10とする具体的な経路もコードに含む。
上の語は、接続路での点の交換と小円での正の全捻りから直接読める。
独立チェックとして点の軌跡をサンプリングし、実部の交差順序と虚部から
6個すべての語を再現した。この経路チェックは浮動小数点であり、区間演算の証明ではない。
行列の計算と恒等式の検証は全て有理関数による厳密計算。

d=q-1 と書けば、共通基底での行列は

\[
M_{y=0}=\begin{pmatrix}q^2&0&0\\qd&1&0\\qd&0&1\end{pmatrix},\qquad
M_{y=x}=\begin{pmatrix}1&d&0\\0&q^2&0\\0&qd&1\end{pmatrix},\qquad
M_{y=1-x}=\begin{pmatrix}1&0&d\\0&1&d\\0&0&q^2\end{pmatrix}.
\]

追加直線については

\[
M_{x=0}=\begin{pmatrix}1+qd&-qd&0\\-d&q&0\\0&0&1\end{pmatrix},\qquad
M_{x=1}=\begin{pmatrix}1+qd&0&-qd\\d^2&1&-d^2\\-d&0&q\end{pmatrix},\qquad
M_{x=1/2}=\begin{pmatrix}1&0&0\\0&1+qd&-qd\\0&-d&q\end{pmatrix}.
\]

同じ行列の階数1表示も列挙すると：

| 超平面H | u_H^T | v_H^T |
|---|---|---|
| y=0 | (q+1,q,q) | (1,0,0) |
| y=x | (1,q+1,q) | (0,1,0) |
| y=1-x | (1,1,q+1) | (0,0,1) |
| x=0 | (q,-1,0) | (1,-1,0) |
| x=1 | (q,q-1,-1) | (1,0,-1) |
| x=1/2 | (0,q,-1) | (0,1,-1) |

前に計算した留数のalpha倍を除いたものをC_Hとし、S=I_3+J_3とすると

\[
\left.\frac{dM_H}{dq}\right|_{q=1}=S C_H S^{-1}
\]

が6本全てで成立する。q=exp(2*pi*i*alpha) なので、alphaの一次項は
2*pi*i*S*C_H*S^-1。これはパラメータと符号を検証する一次の比較であり、
有限のalphaでM_H=exp(2*pi*i*alpha*S*C_H*S^-1)を主張するものではない。

さらに標準full twistの行列積は

\[
M_{12}M_{13}M_{23}M_{14}M_{24}M_{34}=q^4I_3.
\]

追加直線の行列はq≠1では非自明で、q=-1でも単位行列にはならない。
この明示的なBurau表現はq≠1では元の3直線配置にそのまま降りない。
q=1で全行列はIとなるが、その特殊化と共鳴したMCとの同定は別問題。

## 3. 実行

```bash
python3 monodromy_matrices.py
```

`monodromy_results.json` に6個の行列、LaTeX、純組紐の語、検算結果を出力する。
任意の純組紐生成元は `pure_matrix(m,i,j,q)`、一般の語は既存の `burau_word` を使う。
このファイルは等重みのBurau型を具体化したもので、任意階数の任意入力に対する
モノドロミーの一律の閉形式を与えるものではない。

標準純組紐生成元の規約の参照：H. A. Haidar–M. N. Abdulrahim,
*On the irreducibility of the representation of the pure braid group on three strands*,
International Journal of Applied Mathematics 32 (2019), 249–257, Definition 4。
https://www.diogenes.bg/ijam/contents/2019-32-2/7/7.pdf
この文献のBurauの基底は本ファイルとは異なる。
