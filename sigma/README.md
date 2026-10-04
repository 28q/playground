# 和・積の計算機

有限和・積を厳密に評価し、式の構造と値から漸化式・母関数を調べ、式専用のC++17コードを生成します。高速版の既定値は **mod 998244353**。Boostによる厳密有理数モードも残しています。

## 起動とGitHub Pages

Node.js 18以上で、このフォルダーから `npm start` を実行し、http://127.0.0.1:4173/ を開きます。追加パッケージは不要です。

GitHub Pagesには、次の10ファイルを**同じ公開フォルダー**へアップロードしてください。ビルド不要です。

```text
index.html
style.css
app.js
engine.js
optimize.js
cpp.js
fast-cpp.js
worker.js
queries.js
query-cpp.js
```

今回の版に `queries.js` と `query-cpp.js` が追加されています。古いファイルも含めて上記一式を更新してください。server.cjs、テスト、Node.js、Boostは公開サイトの動作には不要です。

## 入力例

```text
motzkin(n)*sum(binomial(n,k)^2,k,0,n)
sum(k^3,k,1,n)
sum(binomial(n,k)*bell(k),k,0,n)
sum(bell(k)*fibonacci(n-k),k,0,n)
sum(binomial(n,k)^2*binomial(n+k,k)^2,k,0,n)
Sum[Binomial[n,k]^2,{k,0,n}]
```

自由変数は n。38種類の関数（うち5種類は配列クエリ用）をボタンと辞典に全件表示します。掛け算の * は省略できません。自然言語・LaTeX・無限和は未対応。空和=0、空積=1、0^0=1です。

## 式変形と推定

- 多項式の展開、Bernoulli多項式による有限和の消去、入れ子の多項式和、二項定理、Vandermonde、hockey-stick、望遠鏡和など。空区間を誤って多項式へ延長しないよう適用条件を確認します。
- 多項式・既知の数列の和や項別積から、定数係数漸化式またはP-recursiveの階数上界を計算。合成時の係数次数は別途探索します。
- Berlekamp–Masseyで定数係数候補を得て、有理数の連立方程式で係数を求め直し、厳密値と照合します。
- P-recursiveは最大階数・係数次数を画面で指定。既定は8階・4次、設定可能範囲は16階・8次。未知係数・探索候補は各36まで。既知のC-finite階数上界がある場合は探索階数をその上界まで拡張します。
- 最後の8項を係数の決定に使わず検算。既知の漸化式から有理関数の恒等式として検証できた候補は「記号的検証済み」、有限値の一致だけなら「推測・未証明」です。
- 有理型母関数と既知のOGF/EGFを扱います。未検出なら、Gの2次・xの2次までの代数的母関数も推測します。

参照記事の構文木に対する多項式和の消去を取り入れています。記事にある一般の条件式・合同式の完全な場合分けまでは実装していません。

## 生成C++

高速版は入力式に必要な演算のみを出力し、関数名を判定する汎用ディスパッチャーを含みません。外部ライブラリー不要です。

| 式の形 | 第N項 | 0〜Nの列 |
| --- | --- | --- |
| 多項式 | Horner法 O(d) | O(Nd) |
| 有理型母関数・定数係数漸化式 | Bostan–Mori O(M(r) log N) | FPS inverse |
| 4階以下のP-recursive | 多項式行列のBSGS・多点評価。小さいNは逐次 | 専用漸化式・分母の一括反転 |
| 5階以上のP-recursive | 専用漸化式・分母の一括反転 | 同左 |
| ベル数・ベルヌーイ数・Euler数・分割数 | FPSでNまで展開して抽出 | NTT・FPS exp/inverse |
| 通常/二項畳み込み・二項変換 | Nまで展開して抽出 | NTTによる母関数演算 |
| 非特異な代数的母関数候補 | Nまで展開して抽出 | FPS Newton法（推測） |

M(t)=O(t log t)。低階のP-recursiveは概ね Õ(r³√(dN)) の体演算です。ベル数等の第N項を対数時間で求めるものではありません。代数的母関数の枝が特異な場合は例外にします。

`coefficient(N)` は a(N)、`coefficients(N)` は初項から添字Nまで（通常N+1個）を返します。途中の開始添字から推定した式はコードに記載したSTARTからです。標準入力にNを渡すmainも含み、SIGMA_NO_MAINで省略できます。

```text
g++ -std=c++17 -O2 sigma.cpp -o sigma
echo 1000000 | ./sigma
```

FPSと列出力は N < 2²¹、P-recursive単項は N < 998244353。多項式・有理型の単項はこの上限を受けません（付属mainは符号付き64bitの非負整数入力）。法上で分母が0となる、初期項の範囲外の特異点では例外にします。

厳密値の直接計算版は従来の汎用処理を含みます。厳密値の専用漸化式版は初期値と漸化式のみです。いずれもBoost.Multiprecisionが必要です。

## 範囲と検証

ブラウザーの厳密評価は24〜128項、開始添字0〜128。通常の関数の引数は0〜256、二項係数と和・積の端点は−512〜512、和・積は最大513項。演算回数・桁数・20秒の時間制限もあります。

任意の名前付き特殊関数を網羅する万能CASではありません。一般的な創造的望遠鏡法、多変数母関数、一般のP-recursiveからの微分方程式解法は未実装です。対応する高速計算法を見つけられない式には、その旨を表示します。有限サンプルの一致は証明ではありません。

`npm test`、`npm run test:optimizer`、`npm run test:fast` で検証できます。最後のコマンドはg++を使用します（CXXで変更可）。詳細は [VALIDATION.md](VALIDATION.md)。

## 構成と参考

engine.jsは厳密演算と構文・関数、optimize.jsは式変形・推定・記号的検証・計算法選択、fast-cpp.jsは高速コード生成、cpp.jsは厳密値コード生成、worker.jsは解析の分離実行です。

- [元の計算機](https://sigma.yosupo.jp/)
- [wata氏の記事：多項式和の式変形](https://wata-orz.hatenadiary.org/entry/20091223/1261582436)
- [Bostan–Mori論文](https://arxiv.org/abs/2008.08822)
- [Bostan–Gaudry–Schost：多項式行列積の高速計算](https://cs.uwaterloo.ca/~eschost/publications/pollard.pdf)
- [関数の規約：DLMF](https://dlmf.nist.gov/26)

外部サービスへの数式送信・外部スクリプトの読み込みはありません。入力式・モード・固定パラメータ・検算用クエリと配列をブラウザーに保存します。

更新時は10ファイルを一式置き換えてください。2.0.1以降は依存ファイルのURLに版を付け、Workerの版も照合します。読み込み不整合の回帰テストは `npm run test:loading` です。

数式プレビューと簡約式はブラウザー標準のMathMLで組版します。Cambria Math / STIX Two Mathなどの端末内の数式フォントを使用し、外部フォントのダウンロードは不要です。

## 複数変数とクエリ（3.0）

「計算対象」で次を選べます。

- **数列 a(n)**：第N項、または0〜Nを出力。mなどの追加変数は「固定パラメータ」に `m=2, a=3` のように指定します。
- **複数の (n,m) クエリ**：両変数を自由に動かして一括処理。a,bなど他の変数は固定値を指定します。検算用の入力は1行に `n m`。自由に動く座標は現在2変数までです。
- クエリ結果の「出力」で **0..N × 0..Mの全表** を選ぶと、N Mを入力として各点の値を出すmainを生成します。全表は100万セルまでで、各式の定義域を満たす必要があります。

```text
sum(binomial(n,k),k,0,m)                        # 二項部分和のMo
stirling2(n,m)                                 # 行NTT / 行共有DP
sum(fibonacci(k),k,n,m)                        # 数列生成＋累積和
sum(sum(gcd(i,j),j,0,m),i,0,n)                 # 二次元累積和
sum(sum(bell(i)*fibonacci(j),j,0,m),i,0,n)      # 変数分離
sum(floor((3*k+2)/m),k,0,n-1)                 # Euclidean floor_sum
sum(mobius(k)*floor(n/k)*floor(m/k),k,1,min(n,m)) # 互いに素な組数
sum(array(k),k,n,m)                            # 累積和 / Fenwick / Segment
prod(array(k),k,n,m)                           # ゼロ対応累積積 / Segment
range_min(n,m)                                # Sparse table / Segment
range_max(n,m)
range_gcd(n,m)
distinct(n,m)                                 # 最終出現＋Fenwick / Mo
bell(n)+fibonacci(m)                           # 独立な数列をFPS等で前計算
```

上記の `#` 以降は説明です。入力欄には式だけを貼り付けてください。計算法と実装状況40項目は画面の「計算法・対応状況」と [ALGORITHMS.md](ALGORITHMS.md) に掲載しています。

複数の計算法が適用できる場合は選択できます。DPとNTTのどちらが速いかはクエリ数・異なるnの個数・最大添字によって変わります。常に最適な手法が先頭になることは保証しません。クエリ用のプランは式の構造・恒等式から選び、有限値だけの推測を採用しません。

### C++の入力とAPI

通常は `Q` の後にQ組の `n m`。配列を使う式は `L`、L個の配列値、`Q`、Q組の `n m` の順です。区間は0始まり・閉区間 [n,m] です。関数 `solve_queries(vector<Query>, data)` は入力順の答えを返します。配列を使わない式はdataを省略できます。

配列のmin/max/gcdは元の整数で比較・計算してからmodに写します。distinctも元の整数の等値性で判定します。たとえば0と998244353はdistinctでは異なる値ですが、積ではどちらもmod上のゼロです。

Fenwick生成コードは `add(i, mint)` と `prefix(end)`、Segment生成コードは `set(i, value)` と `query(l,r)` を持ちます。これらのクラスの区間は **[l,r)** です。付属mainと画面は静的クエリのみで、更新操作列のUIは未実装です。

クエリ用C++はGCC/Clangの `__int128` を使用。法は998244353。全般にQ≤100万、各座標≤10^12。配列値・算術関数の整数引数は絶対値10^12以下。配列・階乗表・FPSは添字2²¹未満、二次元累積和は境界込み400万セル、行DPは(U+1)(V+1)≤3000万。手法ごとの制限・ゼロ除算は実行時に検査します。区間クエリ生成コードは n≤m を必要とし、数学的な空区間の延長はしません。

画面での直接検算は各座標256以下、最大128クエリ、配列1万要素までです。大きい座標は検算範囲外と表示しますが、生成C++は各手法の上限まで使用できます。

### 実装と検証

追加の公開用ファイルは **queries.js と query-cpp.js** です。公開時には既存8ファイルと合わせて10ファイルを同時更新してください。

```text
npm run test:queries
npm run test:query-edges
```

生成C++の値比較、行NTTと独立なDPの比較、大きな床関数和、配列更新API、入力mainを検証します。

データ構造の定義・計算量の参考：[ACL Fenwick tree](https://atcoder.github.io/ac-library/production/document_en/fenwicktree.html)、[ACL Segment tree](https://atcoder.github.io/ac-library/production/document_en/segtree.html)、[ACL floor_sum](https://atcoder.github.io/ac-library/production/document_en/math.html)。生成コードは自己完結した実装で、ACLのインストールは不要です。
