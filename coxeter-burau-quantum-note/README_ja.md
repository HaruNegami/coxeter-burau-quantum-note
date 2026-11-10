# 一般化Burau表現と量子ゲート

根上春（Haru Negami）による研究ノートと再現用コードです。

[本文PDF](coxeter_burau_quantum.pdf) · [TeX](coxeter_burau_quantum.tex) ·
[English README](README.md)

一般化Burau表現の行列、不変エルミート形式、正定値条件、量子ゲートへの解釈を
まとめています。E6・E7・E8型の閉包については、既存のFLW二固有値分類からの
導出とrank-twoの直接証明を収録しています。一般配置のconvolutionと基点付き
モノドロミーは付録で扱っています。

## 計算できること

- E6・E7・E8の行列と不変形式を生成し、関係式を厳密検証する。
- 共通ユニタリ化を計算し、E8の場合は3量子ビットのゲート行列を出力する。
- 一般配置の留数、中間商、最大降下商を計算する。
- Fox微分、純組紐生成元、6直線例のモノドロミー行列を検算する。

## 実行

Python 3.10以降とNumPy、SymPyを使用します。仮想環境の作成後、
リポジトリ直下で実行してください。

```bash
python -m pip install -r requirements.txt
python verify_coxeter_gates.py
```

標準では `theta/pi = 1/60` のE8ゲートを出力します。厳密検算と数値検算を
区別して記録しており、密度や有限像は数値実験によって判定していません。

一般配置のコードは `legacy_arrangement/` にあります。
`make verify` で全5スクリプトの既存検証を実行できます。
詳しいコマンドと制約は英語READMEと各ディレクトリの説明を参照してください。

本文を再生成するには `make pdf`、または `pdflatex` を2回実行します。

## 文献上の位置づけ

既知の構成とその帰結を、明示計算・証明・実行可能なコードとともにまとめた
ノートです。一般化Burau表現と不変形式はSquierの既存構成を使用し、
E型の閉包についてはFLWの分類定理への依存を明記しています。
任意の配置にDynkin型が付くことや、全Coxeter型のBurau表現が今回の
middle convolutionから得られることは主張していません。
