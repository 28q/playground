(function(root){
'use strict';
const E=root.SigmaEngine||(typeof require!=='undefined'?require('./engine.js'):null),O=root.SigmaOptimize||(typeof require!=='undefined'?require('./optimize.js'):null);
const num=x=>({type:'num',value:E.q(x)}),v=name=>({type:'var',name}),bin=(op,left,right)=>({type:'bin',op,left,right});
const isVar=(a,name)=>a&&a.type==='var'&&a.name===name,isNum=(a,x)=>a&&a.type==='num'&&a.value.eq(x),same=(a,b)=>E.format(a)===E.format(b);
function parameters(text=''){const out={};for(const item of text.split(/[,;\n]/).map(x=>x.trim()).filter(Boolean)){const match=/^([a-zA-Z_][\w]*)\s*=\s*(-?\d+)(?:\s*\/\s*(\d+))?$/.exec(item);if(!match)throw Error('固定パラメータは a=2, b=3 の形式で入力してください');const name=match[1].toLowerCase();if(name==='n'||name==='__proto__'||name==='constructor'||name==='prototype')throw Error('このパラメータ名は使えません: '+name);if(Object.hasOwn(out,name))throw Error('パラメータが重複しています: '+name);out[name]=new E.Q(match[2],match[3]||1);}return out;}
function bind(source,params={},query=false){if(query&&Object.hasOwn(params,'m'))throw Error('クエリモードでは m は入力座標です。固定パラメータから m を外してください');let ast=E.parse(source,[...(query?['n','m']:['n']),...Object.keys(params)]);for(const [name,value] of Object.entries(params))ast=O.substitute(ast,name,num(value));return ast;}
function points(text){let lines=text.trim().split(/\n/).filter(x=>x.trim());if(!lines.length||lines.length>128)throw Error('画面で検算するクエリは1〜128行です');return lines.map((line,i)=>{let a=line.trim().split(/[\s,]+/).map(Number);if(a.length!==2||a.some(x=>!Number.isSafeInteger(x)||x<0||x>1000000000000))throw Error((i+1)+'行目は非負整数 n m（各10^12以下）を入力してください');return a;});}
function dataset(text){if(!text.trim())return [];const a=text.trim().split(/[\s,]+/).map(x=>{if(!/^-?\d+$/.test(x)||BigInt(x)<-1000000000000n||BigInt(x)>1000000000000n)throw Error('配列は絶対値10^12以下の整数を入力してください');return x;});if(a.length>10000)throw Error('画面の配列は10,000要素までです');return a;}
function arrayFunctions(data){const index=q=>{if(q.b!==1n||q.a<0||q.a>=data.length)throw Error('配列の添字が範囲外です');return Number(q.a)};const range=(l,r)=>{let a=index(l),b=index(r);if(a>b)throw Error('配列の区間は n ≤ m が必要です');return data.slice(a,b+1).map(E.q)};return {array:x=>E.q(data[index(x)]),distinct:(l,r)=>E.q(new Set(range(l,r).map(String)).size),range_min:(l,r)=>range(l,r).reduce((a,b)=>a.sub(b).a<0n?a:b),range_max:(l,r)=>range(l,r).reduce((a,b)=>a.sub(b).a>0n?a:b),range_gcd:(l,r)=>range(l,r).reduce((a,b)=>E.fn('gcd',[a,b]),E.q(0))};}
function evaluate(ast,n,m,data=[]){return E.evaluate(ast,n,{steps:0,functions:arrayFunctions(data)},{n:E.q(n),m:E.q(m)});}
function flatten(a){return a.type==='bin'&&a.op==='*'?[...flatten(a.left),...flatten(a.right)]:[a]}
const product=a=>a.reduce((x,y)=>x?bin('*',x,y):y,null)||num(1);
function scalar(a){return a.type==='num'||a.type==='var'&&['n','m','i','j'].includes(a.name)||a.type==='neg'&&scalar(a.x)||a.type==='bin'&&['+','-','*'].includes(a.op)&&scalar(a.left)&&scalar(a.right);}
function direct(a){if(a.type==='num'||a.type==='var')return true;if(a.type==='neg')return direct(a.x);if(a.type==='bin')return direct(a.left)&&direct(a.right)&&(a.op!=='^'||scalar(a.right));return a.type==='call'&&['binomial','factorial','rising','falling','gcd','lcm','min','max','abs'].includes(a.name)&&a.args.every(scalar);}
function sequence(ast,index){let a=O.substitute(ast,index,v('n'));try{return O.seriesPlan(O.simplify(a))}catch{return null}}
const choice=(id,name,cost,condition,extra={})=>({id,name,cost,condition,...extra});
function independentSeries(a){
 if(a.type==='num')return {kind:'constant',value:a.value};
 for(const name of ['n','m'])if(!O.free(a,name==='n'?'m':'n')){const p=sequence(a,name);if(p)return {kind:'lookup',variable:name,sequence:p};}
 if(a.type==='neg'){let child=independentSeries(a.x);return child&&{kind:'neg',child};}
 if(a.type==='bin'&&['+','-','*','/'].includes(a.op)){let left=independentSeries(a.left),right=independentSeries(a.right);if(left&&right)return {kind:'combine',op:a.op,left,right};}
 return null;
}
function plan(ast){let choices=[],terms=O.multivariate(ast);if(terms&&terms.every(t=>Object.keys(t.powers).every(x=>['n','m'].includes(x))))choices.push(choice('polynomial','多変数多項式の和の消去','1クエリ O(Td)、追加メモリ O(T)','多項式・有限和。Tは単項式数、dは次数。',{terms}));
 if(ast.type==='call'&&['distinct','range_min','range_max','range_gcd'].includes(ast.name)&&isVar(ast.args[0],'n')&&isVar(ast.args[1],'m')){
  if(ast.name==='distinct'){
   choices.push(choice('distinct-fenwick','右端でソート＋最終出現位置＋Fenwick tree','O((L+Q) log L)、メモリ O(L+Q)','静的配列、0 ≤ n ≤ m < L。値を座標圧縮。',{array:true}));
   choices.push(choice('distinct-mo','Mo’s algorithm＋座標圧縮','O((L+Q)√L + L log L + Q log Q)、メモリ O(L+Q)','静的配列、端点を動かして頻度を更新。',{array:true}));
  }else {const op=ast.name.slice(6);choices.push(choice('sparse','Sparse table','前処理 O(L log L)、クエリ O(1)','結合的かつ冪等な min / max / gcd。',{array:true,op}));choices.push(choice('segment','Segment tree','構築 O(L)、クエリ・点代入 O(log L)','生成クラスの set() で更新可能。画面とmainは静的クエリ。',{array:true,op}));}
 }
 if(ast.type==='call'&&['stirling1','stirling1s','stirling2','eulerian'].includes(ast.name)&&isVar(ast.args[0],'n')&&isVar(ast.args[1],'m')){
  if(ast.name==='stirling2')choices.push(choice('stirling-fps','スターリング第二種の行をNTTで一括生成','O(D(V log V + V log U) + Q)','同じnをまとめ、各行 S(n,0..V) を畳み込み。V < 2²¹。',{function:ast.name}));
  choices.push(choice('row-dp','行を共有するDP＋オフライン整列','O(UV+Q log Q)、メモリ O(V+Q)','U=max n、V=min(max m,U)。UV ≤ 3×10^7。',{function:ast.name}));
 }
 if(ast.type==='call'&&['sum','prod'].includes(ast.name)){
  const [body,k,lo,hi]=ast.args,kn=k.name;
  if(ast.name==='sum'&&isNum(lo,0)&&isVar(hi,'m')&&body.type==='call'&&body.name==='binomial'&&isVar(body.args[0],'n')&&same(body.args[1],k))choices.push(choice('binomial-mo','二項係数の部分和をMo’s algorithmで処理','O(U + Q log Q + (U+Q)√U)、メモリ O(U+Q)','S(n,m)=Σₖ₌₀ᵐ C(n,k)。n方向とm方向にO(1)で遷移。U < 2²¹。'));
  if(isVar(lo,'n')&&isVar(hi,'m')&&!O.free(body,'n')&&!O.free(body,'m')){
   const array=body.type==='call'&&body.name==='array'&&same(body.args[0],k),sp=array?null:sequence(body,kn);
   if(array||sp){let op=ast.name==='sum'?'sum':'product';choices.push(choice(op==='sum'?'prefix':'product-prefix',op==='sum'?'累積和':'累積積＋逆元＋ゼロ個数','前処理は数列生成＋O(U)、クエリ O(1)','0 ≤ n ≤ m の静的区間。積は法上のゼロを別管理。',{array,sequence:sp,op}));if(array){choices.push(choice('segment','Segment tree','構築 O(L)、クエリ・点代入 O(log L)','set() で更新可能。画面とmainは静的クエリ。',{array:true,op}));if(op==='sum')choices.push(choice('fenwick','Fenwick tree','構築 O(L)、クエリ・点加算 O(log L)','add() で更新可能。画面とmainは静的クエリ。',{array:true,op}));}}
  }
  if(ast.name==='sum'&&isNum(lo,0)&&same(hi,bin('-',v('n'),num(1)))&&body.type==='call'&&body.name==='floor'&&body.args[0].type==='bin'&&body.args[0].op==='/'&&isVar(body.args[0].right,'m')){
   let p=O.polynomial(O.substitute(body.args[0].left,kn,v('n')));if(p&&p.length<=2&&p.every(x=>x.b===1n&&x.a>=0n&&x.a<=1000000000000n))choices.push(choice('floor-sum','Euclidean floor_sum','1クエリ O(log(m+a+1))、メモリ O(1)','Σₖ₌₀ⁿ⁻¹ floor((ak+b)/m)。m ≥ 1、固定 a,b ≥ 0。',{a:p[1]||E.q(0),b:p[0]}));
  }
  if(ast.name==='sum'&&isNum(lo,1)&&hi.type==='call'&&hi.name==='min'&&isVar(hi.args[0],'n')&&isVar(hi.args[1],'m')){
   let factors=flatten(body),weight='one',seen=new Set(),ok=true;for(let f of factors){if(f.type==='call'&&f.name==='floor'&&f.args[0].type==='bin'&&f.args[0].op==='/'&&same(f.args[0].right,k)&&f.args[0].left.type==='var'&&['n','m'].includes(f.args[0].left.name)&&!seen.has(f.args[0].left.name))seen.add(f.args[0].left.name);else if(weight==='one'&&f.type==='call'&&['totient','mobius'].includes(f.name)&&same(f.args[0],k))weight=f.name;else ok=false;}if(ok&&seen.size===2)choices.push(choice('quotient','商が一定の区間でまとめる'+(weight==='one'?'':'＋線形篩'),'1クエリ O(√n+√m)'+(weight==='one'?'':'、篩 O(U)'),'Σ w(k) floor(n/k) floor(m/k)。w=1, φ, μ。',{weight}));
  }
  if(ast.name==='sum'&&isNum(lo,0)&&isVar(hi,'n')&&body.type==='call'&&body.name==='sum'&&isNum(body.args[2],0)&&isVar(body.args[3],'m')&&kn!==body.args[1].name){
   let [cell,j]=body.args;const jn=j.name;if(!O.free(cell,'n')&&!O.free(cell,'m')){
    let left=[],right=[],separable=true;for(let f of flatten(cell)){if(!O.free(f,jn))left.push(f);else if(!O.free(f,kn))right.push(f);else separable=false;}
    if(separable){let a=sequence(product(left),kn),b=sequence(product(right),jn);if(a&&b)choices.push(choice('separable','変数分離＋一次元累積和','数列生成＋O(U+V+Q)、メモリ O(U+V)','被加数が f(i)g(j) に分離できる矩形和。',{left:a,right:b}));}
    cell=O.substitute(O.substitute(cell,kn,v('_row')),jn,v('_col'));cell=O.substitute(O.substitute(cell,'_row',v('i')),'_col',v('j'));if(direct(cell))choices.push(choice('prefix2d','二次元累積和','前処理 O(UV × セル評価)、クエリ O(1)','矩形 [0,n]×[0,m]。最大4,000,000セル。',{cell}));
   }
  }
 }
 const independent=independentSeries(ast);if(independent&&!choices.length&&!direct(ast))choices.push(choice('sequence-expression','独立な数列の一括生成＋参照','各数列の生成コスト＋O(Q × 式の大きさ)','nだけ / mだけに依存する部分式をFPS等で前計算。各最大添字 < 2²¹。',{expression:independent}));
 if(direct(ast)&&!choices.some(x=>x.id==='polynomial'))choices.push(choice('direct','式専用の算術・階乗テーブル','前処理 O(F)、1クエリは式の演算数（累乗はlog）','二項係数等の上側引数は非負、階乗表 F < 2²¹。',{expression:ast}));
 return choices;
}
function analyze(source,options={}){const params=parameters(options.parameters||''),ast=bind(source,params,true),queries=points(options.queries||'0 0\n3 2\n5 3'),data=dataset(options.array||''),choices=plan(ast);let values=queries.map(([n,m])=>{try{if(n>256||m>256)throw Error('画面の直接検算は各座標256以下。生成C++では各手法の上限まで対応');return {n,m,value:evaluate(ast,n,m,data),error:null}}catch(e){return {n,m,value:null,error:e.message}}});return {mode:'queries',source,ast,params,choices,values,queries,dataLength:data.length};}
const catalog=[
 ['多変数多項式・Faulhaber','実装済み','複数の有限和を厳密に消去。固定パラメータも代入可能。'],
 ['Horner法','実装済み','多項式の単項・全項評価。'],['NTT畳み込み','実装済み','通常畳み込み、二項畳み込み、二項変換。'],['FPS inverse / log / exp','実装済み','ベル数・ベルヌーイ数・Euler数・分割数。'],['FPS Newton法','実装済み','非特異な2次代数的母関数。推測は未証明と表示。'],['Berlekamp–Massey','実装済み','有限列から定数係数候補を推定し厳密検算。'],['Bostan–Mori','実装済み','有理型母関数の第N項。'],['P-recursive行列BSGS','実装済み','4階以下の単項計算。'],['多点評価・部分積木','実装済み','多項式行列BSGSの内部で使用。'],['分母の一括反転','実装済み','漸化式による列の生成。'],['累積和・累積積','実装済み','静的区間クエリ。積のゼロも処理。'],['二次元累積和','実装済み','二変数の矩形和。'],['変数分離','実装済み','二重和を独立した一次元列へ変換。'],['Mo’s algorithm','実装済み','二項係数の部分和・配列のdistinct。'],['オフライン整列＋Fenwick','実装済み','右端順に処理するdistinct。'],['座標圧縮','実装済み','負数・大きい値を含む配列のdistinct。'],['Fenwick tree','実装済み','区間和。生成クラスは点加算APIを持つ。'],['Segment tree','実装済み','和・積・min・max・gcd。生成クラスは点代入APIを持つ。'],['Sparse table','実装済み','静的なmin・max・gcd。'],['商の区間分割','実装済み','floor(n/k)floor(m/k)とφ・μの重み。'],['線形篩・Möbius反転の形','実装済み','φとμの重み付き商和、互いに素な組数など。'],['Euclidean floor_sum','実装済み','線形分子の床関数和。'],['DPの行共有・ローリング配列','実装済み','スターリング数・Eulerian数の複数クエリ。'],['スターリング第二種の行NTT','実装済み','同じnをまとめて全kを畳み込みで取得。'],
 ['遅延Segment tree・双対Segment tree','未実装','区間更新の操作列・作用素の定義が必要。'],['永続Segment tree・Wavelet matrix','未実装','区間k番目・頻度など。'],['Disjoint sparse table・SWAG','未実装','一般の結合的演算・移動窓。'],['更新付きMo・木上Mo','未実装','時刻軸や木の入力が必要。'],['CDQ分割統治・オフライン二次元BIT','未実装','順序制約と座標を持つクエリ。'],['Convex hull trick・Li Chao tree','未実装','線形遷移・最小値DPの証明条件が必要。'],['分割統治DP・Knuth最適化・SMAWK','未実装','Monge性や最適位置の単調性を証明する必要。'],['subset / superset zeta・Möbius変換','未実装','集合添字の入力形式が必要。'],['FWHT・部分集合畳み込み','未実装','XOR / AND / ORや部分集合添字。'],['Relaxed convolution・CDQ畳み込み','未実装','オンラインの相互依存漸化式。'],['Kitamasa・汎用行列累乗','未実装','定数係数単項にはBostan–Moriを使用。'],['Lucas・素数冪の二項係数・CRT','未実装','大きい階乗添字や別の法への拡張。'],['Min_25・杜教篩・商集合DP','未実装','巨大な算術関数の和。'],['多変数母関数・creative telescoping','未実装','任意の多変数式の自動漸化式証明。'],['木DP・rerooting・HLD・binary lifting','未実装','木・グラフの入力形式が必要。'],['FFT任意mod・多項式補間・Taylor shift','未実装','専用の式認識と生成APIが必要。']
];
function run(source,options={}){if(options.mode==='queries')return analyze(source,options);const params=parameters(options.parameters||''),ast=bind(source,params,false);const r=E.analyze(E.format(ast),options);r.originalSource=source;r.parameters=params;return r;}
root.SigmaQueries={run,parameters,bind,points,dataset,evaluate,plan,analyze,catalog};if(typeof module!=='undefined')module.exports=root.SigmaQueries;
})(typeof globalThis!=='undefined'?globalThis:this);
