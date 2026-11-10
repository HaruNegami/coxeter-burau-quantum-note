# 一般化Burau表現と量子ゲート — v2

根上春（Haru Negami）による研究ノートと再現用コードです。

[本文PDF](coxeter_burau_quantum.pdf) · [TeX](coxeter_burau_quantum.tex) ·
[English README](README.md) · [変更履歴](CHANGELOG.md)

初版の本文・付録・計算コードを残し、第7節に
**一般の超平面配置 → 1次元seedへのMC → 不変形式と符号判定 → 量子ゲート**
の構成を追加しました。不変形式と符号条件には著者の先行研究を用いています。

## 追加した内容

- 実際のMC行列と、異なるseed重みを許す色付き組紐の移送。
- 不変部分空間 `K+L` によるmiddle quotientと商上の共通ユニタリ化。
- 不変エルミート形式の符号・正定値条件と、共鳴時の階数低下。
- 無理数位相に対する完全なユニタリ閉包・非Clifford性の十分条件。
- 有理位相かつ有限像のbase characterに対するGalois符号による有限像判定。
- 3準位ゲート例、共鳴する2次元例、正定値だが無限像となる有理例。
- Dynkin型をMCの入力としたときの方向、穿孔数、配置の拡大。

第1–6節のSquierのCoxeter–Burau表現と、第7節のMC構成は区別します。
前者のE8型8次元ゲートを、後者のMC出力と同一視していません。
初版のPDF・TeXは [`versions/v1/`](versions/v1/) に保存しています。

## 実行

Python 3.11以降、NumPy、SymPyを使用します。

```bash
python -m pip install -r requirements.txt
python mc_arrangement_gates.py --verify --verify-symbolic
```

MC行列・不変形式・商・ユニタリ行列と検証結果を `mc_gate_example.json` に出力します。
標準例は3直線 `y=0, y=x, y=1-x` に衝突直線を加えた配置上の3準位ゲートです。
seedの角度は `(1/20,1/15,1/12)`、MCパラメータの角度は `sqrt(2)/20` です。

```bash
make verify-mc          # MCの164例・記号計算・3つの出力例
make verify-coxeter     # 初版のE6/E7/E8計算
make verify-arrangements # 初版の付録の4スクリプト
make verify            # 上記すべて
make pdf               # v2本文PDFを再生成
```

Makeを使わない場合の実行方法は英語READMEにあります。
初版のE8ゲートは `python verify_coxeter_gates.py` で生成できます。

## MCの条件と定義域

`q_i=exp(2*pi*i*a_i)`、`lambda=exp(2*pi*i*ell)`、
`0<=a_i<1`、`0<ell<1` とします。非零の `a_i` の数を `p`、
`T=ell+sum(a_i)`、`m=floor(T)`、`delta=1`（Tが整数）、それ以外は0とすると、
MC商の次元は `p-delta`、形式の符号は `(p-m,m-delta)` です。

非零の商に対し、正定値条件は `0<T<=1`、負定値条件は `p<=T<p+1` です。
負定値の場合も符号を反転すればユニタリ化できます。2次元以上で、形式が定値かつ
すべてのactiveな `ell+a_i` が無理数なら、閉包は全ユニタリ群になります。
これは十分条件であり、全パラメータの閉包分類を主張するものではありません。

一般にはMC出力は元の配置Aの方向閉包 `A^Y` 上にあります。
追加した超平面の周りのモノドロミーが自明な場合に元の配置へ降ります。
コードは与えたbase loopの純組紐語を評価できますが、任意の配置方程式から
全メリディアンの語を自動決定するものではありません。

数値行列は近似値です。密度と有限像は本文の証明に基づきます。
本文・コードの作成にはOpenAI ChatGPTを大幅に使用しています。
検算の成功は証明の査読や新規性の確認を意味しません。
