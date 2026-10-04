(function(root){
'use strict';
const abs=x=>x<0n?-x:x, gcd=(a,b)=>{a=abs(a);b=abs(b);while(b){[a,b]=[b,a%b]}return a};
class Q{
 constructor(a=0,b=1){if(a instanceof Q){this.a=a.a;this.b=a.b;return}a=BigInt(a);b=BigInt(b);if(!b)throw Error('ゼロ除算です');if(b<0n){a=-a;b=-b}let g=gcd(a,b);this.a=a/g;this.b=b/g;if(abs(this.a).toString(2).length>40000)throw Error('結果が大きすぎます（40,000 bit 制限）')}
 add(y){y=q(y);return new Q(this.a*y.b+y.a*this.b,this.b*y.b)} sub(y){return this.add(q(y).neg())} neg(){return new Q(-this.a,this.b)} mul(y){y=q(y);return new Q(this.a*y.a,this.b*y.b)} div(y){y=q(y);return new Q(this.a*y.b,this.b*y.a)} eq(y){y=q(y);return this.a===y.a&&this.b===y.b} zero(){return this.a===0n} toString(){return this.b===1n?''+this.a:this.a+'/'+this.b} toJSON(){return this.toString()}
}
const q=x=>x instanceof Q?x:new Q(x), zero=()=>q(0), one=()=>q(1);
function integer(x,lo=-10000,hi=10000){x=q(x);if(x.b!==1n||x.a<BigInt(lo)||x.a>BigInt(hi))throw Error(`整数 ${lo}〜${hi} が必要です`);return Number(x.a)}
function pow(x,k){k=integer(k,-1000,1000);if(k<0)return one().div(pow(x,-k));let r=one();for(;k;k=Math.floor(k/2),x=x.mul(x))if(k%2)r=r.mul(x);return r}
const memo=new Map();
function fact(n){let r=one();for(let i=2;i<=n;i++)r=r.mul(i);return r}
function choose(n,k){if(k<0)return zero();if(n>=0&&k>n)return zero();let r=one();for(let i=1;i<=k;i++)r=r.mul(n-i+1).div(i);return r}
const CATALOG=[
 ['array','配列要素','array(n)','クエリモードで入力した配列の0始まりの要素。','配列','array(i) = A[i]'],
 ['distinct','区間の異なる値の個数','distinct(n,m)','配列の閉区間 [n,m] に現れる値の種類数。','配列','distinct(l,r) = |{A[i] : l ≤ i ≤ r}|'],
 ['range_min','区間最小値','range_min(n,m)','配列の閉区間 [n,m] の最小値。','配列','min A[l..r]'],
 ['range_max','区間最大値','range_max(n,m)','配列の閉区間 [n,m] の最大値。','配列','max A[l..r]'],
 ['range_gcd','区間GCD','range_gcd(n,m)','配列の閉区間 [n,m] の最大公約数。','配列','gcd A[l..r]'],
 ['binomial','二項係数','binomial(n,k)','n 個から k 個を選ぶ方法の数。整数 n の一般化にも対応。','組合せ','binomial(n,k) = n(n−1)…(n−k+1)/k!'],
 ['factorial','階乗','factorial(n)','n 個のものを並べる方法の数。0! = 1。','組合せ','n! = n·(n−1)!'],
 ['doublefactorial','二重階乗','doublefactorial(n)','2 ずつ減らして掛ける階乗。(-1)!! = 0!! = 1。','組合せ','n!! = n·(n−2)!!'],
 ['rising','上昇階乗','rising(n,3)','連続する k 個の数を増加方向に掛ける。','組合せ','rising(x,k) = x(x+1)…(x+k−1)'],
 ['falling','下降階乗','falling(n,3)','連続する k 個の数を減少方向に掛ける。','組合せ','falling(x,k) = x(x−1)…(x−k+1)'],
 ['stirling1','第一種スターリング数（符号なし）','stirling1(n,k)','n 個の要素から k 個の巡回を持つ順列を作る方法の数。','組合せ','c(n,k) = c(n−1,k−1)+(n−1)c(n−1,k)'],
 ['stirling1s','第一種スターリング数（符号付き）','stirling1s(n,k)','符号規約 s(n,k) = (−1)^(n−k)c(n,k)。','組合せ','s(n,k) = s(n−1,k−1)−(n−1)s(n−1,k)'],
 ['stirling2','第二種スターリング数','stirling2(n,k)','n 個の要素を k 個の空でない集合に分割する方法の数。','組合せ','S(n,k) = S(n−1,k−1)+kS(n−1,k)'],
 ['bell','ベル数','bell(n)','n 個の要素からなる集合の分割総数。B₀ = 1。','数列','EGF: exp(exp(x)−1)'],
 ['catalan','カタラン数','catalan(n)','正しい括弧列や根付き二分木を数える数列。','数列','Cₙ = binomial(2n,n)/(n+1)'],
 ['fibonacci','フィボナッチ数','fibonacci(n)','F₀ = 0, F₁ = 1 で、前の二つの項の和。','数列','OGF: x/(1−x−x²)'],
 ['lucas','リュカ数','lucas(n)','L₀ = 2, L₁ = 1 のフィボナッチ型数列。','数列','OGF: (2−x)/(1−x−x²)'],
 ['derangement','完全順列数','derangement(n)','元の位置に残る要素が一つもない順列の数。','数列','Dₙ = (n−1)(Dₙ₋₁+Dₙ₋₂)'],
 ['partition','整数分割数','partition(n)','n を正の整数の和に分割する方法の数（順序を区別しない）。','数列','OGF: ∏ⱼ≥₁ 1/(1−xʲ)'],
 ['motzkin','モツキン数','motzkin(n)','水平・上昇・下降の歩幅で高さが負にならない道の数。','数列','M₀=M₁=1; (n+2)Mₙ=(2n+1)Mₙ₋₁+(3n−3)Mₙ₋₂'],
 ['schroeder','大シュレーダー数','schroeder(n)','大シュレーダー数。初項は 1, 2, 6, 22。','数列','Rₙ = Rₙ₋₁ + Σₖ₌₀ⁿ⁻¹ RₖRₙ₋₁₋ₖ'],
 ['eulerian','オイラー数（Eulerian）','eulerian(n,k)','n 要素の順列で降下が k 個あるものの数。Euler 数とは別。','組合せ','A(n,k)=(k+1)A(n−1,k)+(n−k)A(n−1,k−1)'],
 ['lah','ラー数（符号なし）','lah(n,k)','n 要素を k 個の空でない順序付きリストに分割する数。','組合せ','L(n,k)=binomial(n−1,k−1)n!/k!'],
 ['bernoulli','ベルヌーイ数','bernoulli(n)','有理数の数列。ここでは B₁ = −1/2 の規約。','有理数','EGF: x/(exp(x)−1)'],
 ['euler','オイラー数（Euler）','euler(n)','sech(x) の指数型母関数係数。E₂ = −1, E₄ = 5。','有理数','EGF: 1/cosh(x)'],
 ['harmonic','調和数','harmonic(n)','Hₙ = 1 + 1/2 + … + 1/n。H₀ = 0。','有理数','OGF: −log(1−x)/(1−x)'],
 ['harmonic2','一般化調和数','harmonic2(n,2)','Hₙ⁽ʳ⁾ = Σₖ₌₁ⁿ 1/kʳ。整数 r に対応。','有理数','harmonic2(n,r) = Σₖ₌₁ⁿ k^(−r)'],
 ['totient','オイラーのトーシェント','totient(n)','1 から n までの整数で n と互いに素なものの数。n ≥ 1。','整数論','φ(n) = n ∏ₚ∣ₙ (1−1/p)'],
 ['mobius','メビウス関数','mobius(n)','平方因子があれば 0、それ以外は素因子数の偶奇による ±1。n ≥ 1。','整数論','μ(1)=1'],
 ['divisors','約数個数','divisors(n)','正の約数の個数。n ≥ 1。','整数論','τ(n) = Σd|n 1'],
 ['divisorsum','約数和','divisorsum(n)','正の約数の和。n ≥ 1。','整数論','σ(n) = Σd|n d'],
 ['gcd','最大公約数','gcd(n,12)','二つの整数の非負の最大公約数。','整数論','gcd(0,0)=0'],
 ['lcm','最小公倍数','lcm(n,12)','二つの整数の非負の最小公倍数。','整数論','lcm(a,b)=|ab|/gcd(a,b)'],
 ['abs','絶対値','abs(n)','整数または有理数の絶対値。','基本','abs(x) = |x|'],
 ['floor','床関数','floor(n/2)','与えられた値以下の最大の整数。','基本','floor(−1/2) = −1'],
 ['ceil','天井関数','ceil(n/2)','与えられた値以上の最小の整数。','基本','ceil(1/2) = 1'],
 ['min','最小値','min(n,3)','二つの有理数の小さい方。','基本','min(a,b)'],
 ['max','最大値','max(n,3)','二つの有理数の大きい方。','基本','max(a,b)']
].map(([id,name,example,description,category,formula])=>({id,name,example,description,category,formula}));
const two=new Set(['distinct','range_min','range_max','range_gcd','binomial','rising','falling','stirling1','stirling1s','stirling2','eulerian','lah','harmonic2','gcd','lcm','min','max']);
const aliases={choose:'binomial',binom:'binomial',c:'binomial',s1:'stirling1',s2:'stirling2',stirlingstirling:'stirling2',stirlings1:'stirling1s',stirlings2:'stirling2',bellb:'bell',catalannumber:'catalan',fibonacci:'fibonacci',fib:'fibonacci',bernoullib:'bernoulli',eulere:'euler',harmonicnumber:'harmonic',partitionsp:'partition',eulerphi:'totient',moebiusmu:'mobius',pochhammer:'rising',product:'prod'};
function fn(name,args){let key=name+':'+args.join(',');if(memo.has(key))return memo.get(key);if(args.length!==(two.has(name)?2:1))throw Error(name+' の引数の数が違います');let a=args[0],b=args[1],r;if(['stirling1','stirling1s','stirling2','eulerian'].includes(name)){let n=integer(a,0,256),k=integer(b,-256,256),row=memo.get(name+':row:'+n);if(row)return k<0||k>n?zero():row[k]}
 if(name==='abs')r=new Q(abs(a.a),a.b);
 else if(name==='floor')r=q(a.a/a.b-(a.a<0n&&a.a%a.b?1n:0n));
 else if(name==='ceil')r=q(a.a/a.b+(a.a>0n&&a.a%a.b?1n:0n));
 else if(name==='min'||name==='max')r=(a.sub(b).a<0n)===(name==='min')?a:b;
 else if(name==='gcd'||name==='lcm'){if(a.b!==1n||b.b!==1n)throw Error('整数が必要です');let g=gcd(a.a,b.a);r=q(name==='gcd'?g:g?abs(a.a*b.a)/g:0n)}
 else if(name==='rising'||name==='falling'){let k=integer(b,0,256);r=one();for(let i=0;i<k;i++)r=r.mul(a.add(name==='rising'?i:-i))}
 else if(name==='binomial'){r=choose(integer(a,-512,512),integer(b,-512,512))}
 else {let n=integer(a,name==='doublefactorial'?-1:0,256),k=b===undefined?0:integer(b,-256,256);
  const f=(id,i,j)=>fn(id,j===undefined?[q(i)]:[q(i),q(j)]);
  if(name==='factorial')r=fact(n);
  else if(name==='doublefactorial'){r=one();for(let i=n;i>0;i-=2)r=r.mul(i)}
  else if(['stirling1','stirling1s','stirling2','eulerian'].includes(name)){
   let row=[one(),...Array.from({length:n+1},zero)];for(let i=1;i<=n;i++){let next=Array.from({length:n+2},zero);for(let j=0;j<=i;j++){let left=j?row[j-1]:zero();next[j]= name==='eulerian'?row[j].mul(j+1).add(left.mul(i-j)):left.add(row[j].mul(name==='stirling2'?j:name==='stirling1s'?1-i:i-1))}row=next}r=k<0||k>n?zero():row[k];memo.set(name+':row:'+n,row);
  }else if(name==='bell'){r=zero();for(let j=0;j<=n;j++)r=r.add(f('stirling2',n,j))}
  else if(name==='catalan')r=choose(2*n,n).div(n+1);
  else if(name==='fibonacci'||name==='lucas'){let x=q(name==='fibonacci'?0:2),y=one();for(let i=0;i<n;i++)[x,y]=[y,x.add(y)];r=x}
  else if(name==='derangement'){let x=one(),y=zero();if(n===0)r=x;else {for(let i=2;i<=n;i++)[x,y]=[y,x.add(y).mul(i-1)];r=y}}
  else if(name==='partition'){let d=Array.from({length:n+1},zero);d[0]=one();for(let i=1;i<=n;i++)for(let j=i;j<=n;j++)d[j]=d[j].add(d[j-i]);r=d[n]}
  else if(name==='motzkin'){let x=one(),y=one();for(let i=2;i<=n;i++)[x,y]=[y,y.mul(2*i+1).add(x.mul(3*i-3)).div(i+2)];r=n===0?x:y}
  else if(name==='schroeder'){r=n===0?one():f(name,n-1);for(let j=0;j<n;j++)r=r.add(f(name,j).mul(f(name,n-1-j)))}
  else if(name==='lah')r=n===0?q(k===0?1:0):k<=0||k>n?zero():choose(n-1,k-1).mul(fact(n)).div(fact(k));
  else if(name==='bernoulli'){r=n===0?one():zero();if(n>0){for(let j=0;j<n;j++)r=r.sub(choose(n+1,j).mul(f(name,j)));r=r.div(n+1)}}
  else if(name==='euler'){r=n===0?one():zero();if(n%2===0&&n>0)for(let j=0;j<n;j+=2)r=r.sub(choose(n,j).mul(f(name,j)))}
  else if(name==='harmonic'||name==='harmonic2'){r=zero();for(let j=1;j<=n;j++)r=r.add(pow(q(j),q(name==='harmonic'?-1:-k)))}
  else if(['totient','mobius','divisors','divisorsum'].includes(name)){if(n<1)throw Error(name+' は n ≥ 1 のみ定義されます');r=zero();if(name==='totient'){for(let i=1;i<=n;i++)if(gcd(BigInt(i),BigInt(n))===1n)r=r.add(1)}else if(name==='mobius'){r=one();let v=n;for(let p=2;p<=v;p++)if(v%p===0){v/=p;r=r.neg();if(v%p===0){r=zero();break}}}else for(let i=1;i<=n;i++)if(n%i===0)r=r.add(name==='divisors'?1:i)}
  else throw Error('未対応の関数: '+name);
 }if(memo.size>30000)memo.clear();memo.set(key,r);return r;
}
function parse(source,variables=['n']){
 if(source.length>3000)throw Error('入力は 3,000 文字以内にしてください');
 source=source.replace(/\[/g,'(').replace(/\]/g,')').replace(/\{/g,'').replace(/\}/g,'').replace(/Σ|∑/g,'sum').replace(/Π|∏/g,'prod').replace(/×|·/g,'*').replace(/÷/g,'/').replace(/−/g,'-').replace(/\*\*/g,'^');
 let ts=[],re=/\s*(?:(\d+(?:\.\d+)?)|([a-zA-Z_][a-zA-Z_0-9]*)|([+\-*/^(),!]))/gy,m,pos=0;
 while(pos<source.length){re.lastIndex=pos;m=re.exec(source);if(!m){if(!source.slice(pos).trim())break;throw Error((pos+1)+' 文字目: 読み取れない文字です')}ts.push(m[1]||m[2]||m[3]);pos=re.lastIndex}
 let at=0,depth=0;const peek=()=>ts[at],eat=t=>{if(ts[at++]!==t)throw Error('「'+t+'」が必要です')};
 function expr(min=0){if(++depth>80)throw Error('式の入れ子が深すぎます');let t=ts[at++],x;
  if(t==='-'||t==='+'){x=expr(25);if(t==='-')x={type:'neg',x}}
  else if(t==='('){x=expr();eat(')')}
  else if(t&&/^\d/.test(t)){let [a,b='']=t.split('.');x={type:'num',value:new Q(a+b,10n**BigInt(b.length))}}
  else if(t&&/^[a-zA-Z_]/.test(t)){let token=t.toLowerCase(),name=peek()==='('?(aliases[token]||token):token;if(peek()==='('){at++;let args=[];if(peek()!==')'){do {args.push(expr());if(peek()!==',')break;at++}while(true)}eat(')');x={type:'call',name,args}}else x={type:'var',name}}
  else throw Error('数・変数・関数を入力してください');
  while(true){let op=peek();if(op==='!'){at++;x={type:'call',name:'factorial',args:[x]};continue}let p={'+':10,'-':10,'*':20,'/':20,'^':30}[op];if(p===undefined||p<min)break;at++;x={type:'bin',op,left:x,right:expr(p+(op==='^'?0:1))}}
  depth--;return x;
 }let ast=expr();if(at!==ts.length)throw Error('式の途中に余分な文字があります。掛け算は * を使ってください');validate(ast,new Set(variables===null?freeVariables(ast):variables));return ast;
}
function freeVariables(ast){const names=new Set();function visit(a,bound){if(a.type==='var'){if(!bound.has(a.name))names.add(a.name);}else if(a.type==='neg')visit(a.x,bound);else if(a.type==='bin'){visit(a.left,bound);visit(a.right,bound);}else if(a.type==='call'){if(a.name==='sum'||a.name==='prod'){if(a.args.length!==4||a.args[1].type!=='var')throw Error('和・積は4引数で添字変数を指定してください');visit(a.args[2],bound);visit(a.args[3],bound);visit(a.args[0],new Set([...bound,a.args[1].name]));}else for(const x of a.args)visit(x,bound);}}visit(ast,new Set());return [...names].sort((a,b)=>{const rank=x=>x==='n'?0:x==='m'?1:2;return rank(a)-rank(b)||a.localeCompare(b);});}
function validate(a,scope){if(a.type==='var'){if(!scope.has(a.name))throw Error('未定義の変数: '+a.name+'（計算対象と固定パラメータを確認してください）')}else if(a.type==='neg')validate(a.x,scope);else if(a.type==='bin'){validate(a.left,scope);validate(a.right,scope)}else if(a.type==='call'){if(a.name==='sum'||a.name==='prod'){if(a.args.length!==4||a.args[1].type!=='var')throw Error('sum(式, k, 下限, 上限) の形式で入力してください');validate(a.args[2],scope);validate(a.args[3],scope);let s=new Set(scope);s.add(a.args[1].name);validate(a.args[0],s)}else {if(!CATALOG.some(x=>x.id===a.name))throw Error('未対応の関数: '+a.name);if(a.args.length!==(two.has(a.name)?2:1))throw Error(a.name+' の引数の数が違います');a.args.forEach(x=>validate(x,scope))}}}
function evaluate(ast,n,budget={steps:0},env={n:q(n)}){if(++budget.steps>600000)throw Error('計算量の上限に達しました。項数または式を小さくしてください');let E=x=>evaluate(x,n,budget,env);
 if(ast.type==='num')return ast.value;if(ast.type==='var'){if(env[ast.name]===undefined)throw Error('変数 '+ast.name+' の値が必要です');return env[ast.name];}if(ast.type==='neg')return E(ast.x).neg();if(ast.type==='bin'){let a=E(ast.left),b=E(ast.right);return ast.op==='+'?a.add(b):ast.op==='-'?a.sub(b):ast.op==='*'?a.mul(b):ast.op==='/'?a.div(b):pow(a,b)}
 if(ast.name==='sum'||ast.name==='prod'){let [body,v,lo,hi]=ast.args,l=integer(E(lo),-512,512),h=integer(E(hi),-512,512);if(h-l>512)throw Error('和・積は 513 項までです');let r=q(ast.name==='sum'?0:1);for(let i=l;i<=h;i++){let x=evaluate(body,n,budget,{...env,[v.name]:q(i)});r=ast.name==='sum'?r.add(x):r.mul(x)}return r}return budget.functions&&budget.functions[ast.name]?budget.functions[ast.name](...ast.args.map(E)):fn(ast.name,ast.args.map(E));
}
function format(a){if(a.type==='num')return a.value.toString();if(a.type==='var')return a.name;if(a.type==='neg')return '-('+format(a.x)+')';if(a.type==='bin')return '('+format(a.left)+a.op+format(a.right)+')';return a.name+'('+a.args.map(format).join(',')+')'}
const num=x=>({type:'num',value:q(x)}),call=(name,...args)=>({type:'call',name,args}),v=name=>({type:'var',name}),bin=(op,left,right)=>({type:'bin',op,left,right});
function simplify(a,steps=[],allowIdentities=true){let old=format(a);if(a.type==='neg'){a={...a,x:simplify(a.x,steps,allowIdentities)};if(a.x.type==='num')a=num(a.x.value.neg())}
 if(a.type==='bin'){a={...a,left:simplify(a.left,steps,allowIdentities),right:simplify(a.right,steps,allowIdentities)};let l=a.left,r=a.right;if(l.type==='num'&&r.type==='num')a=num(evaluate(a,0));else if((a.op==='+'||a.op==='-')&&r.type==='num'&&r.value.zero())a=l;else if(a.op==='*'&&r.type==='num'&&r.value.eq(1))a=l;else if(a.op==='*'&&l.type==='num'&&l.value.eq(1))a=r;else if(a.op==='/'&&r.type==='num'&&r.value.eq(1))a=l}
 if(a.type==='call'){a={...a,args:a.args.map((x,i)=>simplify(x,steps,allowIdentities&&!(i===0&&(a.name==='sum'||a.name==='prod')&&a.args[1].name==='n')))};if(a.name==='sum'||a.name==='prod'){let [body,k,lo,hi]=a.args,N=v('n');let bounds=allowIdentities&&k.name!=='n'&&lo.type==='num'&&hi.type==='var'&&hi.name==='n';if(bounds&&lo.value.eq(0)&&a.name==='sum'&&body.type==='call'&&body.args.length===2&&format(body.args[0])==='n'&&format(body.args[1])===format(k)){if(body.name==='binomial')a=bin('^',num(2),N);if(body.name==='stirling2')a=call('bell',N);if(body.name==='stirling1')a=call('factorial',N)}
 if(bounds&&lo.value.eq(0)&&a.name==='sum'&&body.type==='bin'&&body.op==='^'&&format(body.right)==='2'&&body.left.type==='call'&&body.left.name==='binomial'&&format(body.left.args[0])==='n'&&format(body.left.args[1])===format(k))a=call('binomial',bin('*',num(2),N),N);
 if(bounds&&lo.value.eq(1)&&format(body)===format(k)&&a.name==='prod')a=call('factorial',N);
 if(bounds&&lo.value.eq(1)&&format(body)===format(k)&&a.name==='sum')a=bin('/',bin('*',N,bin('+',N,num(1))),num(2));
 }}if(format(a)!==old)steps.push(old+' → '+format(a));return a;
}
// Exact rational nullspace. At least eight equations beyond the parameter count.
function nullVector(rows,cols){let a=rows.map(r=>r.map(q)),piv=[],row=0;for(let c=0;c<cols&&row<a.length;c++){let k=a.findIndex((r,i)=>i>=row&&!r[c].zero());if(k<0)continue;[a[row],a[k]]=[a[k],a[row]];let d=a[row][c];for(let j=c;j<cols;j++)a[row][j]=a[row][j].div(d);for(let i=0;i<a.length;i++)if(i!==row&&!a[i][c].zero()){let f=a[i][c];for(let j=c;j<cols;j++)a[i][j]=a[i][j].sub(f.mul(a[row][j]))}piv.push(c);row++}if(piv.length!==cols-1)return null;let free=Array.from({length:cols},(_,i)=>i).find(i=>!piv.includes(i)),x=Array.from({length:cols},zero);x[free]=one();for(let i=piv.length-1;i>=0;i--)x[piv[i]]=a[i][free].neg();let den=1n;for(let y of x)den=den/gcd(den,y.b)*y.b;let z=x.map(y=>y.a*(den/y.b)),g=z.reduce(gcd,0n);z=z.map(y=>y/g);if(z[z.length-1]<0n)z=z.map(y=>-y);return z.map(q)}
function polyEval(p,n){return p.reduceRight((s,c)=>s.mul(n).add(c),zero())}
function guess(seq,start=0){if(seq.every(x=>x.zero()))return {order:1,degree:0,p:[[q(0)],[q(1)]],train:seq.length-8,holdout:8};let train=seq.length-8;
 for(let cost=2;cost<=12;cost++)for(let d=0;d<=3;d++)for(let r=1;r<=4;r++){let cols=(r+1)*(d+1);if(cols!==cost||train-r<cols+2)continue;let rows=[];for(let i=0;i<train-r;i++){let row=[];for(let j=0;j<=r;j++)for(let k=0;k<=d;k++)row.push(seq[i+j].mul(pow(q(start+i),q(k))));rows.push(row)}let x=nullVector(rows,cols);if(!x)continue;let p=Array.from({length:r+1},(_,j)=>x.slice(j*(d+1),(j+1)*(d+1)));if(p[r].every(c=>c.zero()))continue;let good=true;for(let i=train-r;i<seq.length-r;i++){let s=zero();for(let j=0;j<=r;j++)s=s.add(polyEval(p[j],start+i).mul(seq[i+j]));if(!s.zero()){good=false;break}}if(good)return {order:r,degree:d,p,train,holdout:8}}return null;
}
function polyText(p,variable='n'){let s='';p.forEach((c,k)=>{if(c.zero())return;let positive=c.a>0n,magnitude=c.a<0n?c.neg():c,term=(k&&magnitude.eq(1)?'':magnitude.toString())+(k?variable+(k>1?['','','²','³'][k]||'^'+k:''):'');s+=(s?(positive?' + ':' − '):(positive?'':'−'))+term});return s||'0'}
function recurrenceText(g){return g.p.map((p,j)=>'['+polyText(p)+'] a(n'+(j?'+'+j:'')+')').join(' + ')+' = 0'}
function algebraicGF(seq){
 const length=seq.length,train=length-8,powers=[Array.from({length},(_,i)=>q(i===0?1:0)),seq];let squared=Array.from({length},zero);for(let n=0;n<length;n++)for(let k=0;k<=n;k++)squared[n]=squared[n].add(seq[k].mul(seq[n-k]));powers.push(squared);
 for(let dx=1;dx<=2;dx++){let cols=3*(dx+1);if(train<cols+4)continue;let row=n=>powers.flatMap(p=>Array.from({length:dx+1},(_,k)=>n<k?zero():p[n-k]));let c=nullVector(Array.from({length:train},(_,n)=>row(n)),cols);if(!c||c.slice(2*(dx+1)).every(x=>x.zero()))continue;let valid=true;for(let n=train;n<length;n++){let r=row(n).reduce((s,x,i)=>s.add(x.mul(c[i])),zero());if(!r.zero()){valid=false;break}}if(valid){let formula=[];for(let j=0;j<3;j++){let p=c.slice(j*(dx+1),(j+1)*(dx+1));if(p.every(x=>x.zero()))continue;formula.push('('+polyText(p,'x')+')'+(j?' G(x)'+(j===2?'²':''):''))}return {kind:'代数的な通常母関数（推測）',formula:formula.join(' + ')+' = 0',proven:false,algebraic:true,coefficients:Array.from({length:3},(_,j)=>c.slice(j*(dx+1),(j+1)*(dx+1)))}}}return null;
}
function generating(ast,g,seq){let s=format(ast),known={ 'bell(n)':['指数型母関数','exp(exp(x) − 1)'],'factorial(n)':['指数型母関数','1 / (1 − x)'],'fibonacci(n)':['通常母関数','x / (1 − x − x²)'],'lucas(n)':['通常母関数','(2 − x) / (1 − x − x²)'],'catalan(n)':['通常母関数','(1 − √(1 − 4x)) / (2x)'],'partition(n)':['通常母関数','∏ⱼ≥₁ (1 − xʲ)⁻¹'],'bernoulli(n)':['指数型母関数','x / (exp(x) − 1)'],'euler(n)':['指数型母関数','1 / cosh(x)'],'harmonic(n)':['通常母関数','−log(1 − x) / (1 − x)'],'derangement(n)':['指数型母関数','exp(−x) / (1 − x)']};if(known[s])return {kind:known[s][0],formula:known[s][1],proven:true};if(g&&g.degree===0){let den=g.p.map(p=>p[0]).reverse(),numer=[];for(let i=0;i<g.order;i++){let x=zero();for(let j=0;j<=i;j++)x=x.add(den[j].mul(seq[i-j]));numer.push(x)}return {kind:'通常母関数（推測）',formula:'('+polyText(numer,'x')+') / ('+polyText(den,'x')+')',proven:false}}return null}
function analyze(source,{count=40,start=0,skipGuess=false}={}){if(!Number.isInteger(count)||count<24||count>128)throw Error('項数は 24〜128 の整数です');if(!Number.isInteger(start)||start<0||start>128)throw Error('開始 n は 0〜128 の整数です');let ast=parse(source),steps=[],simple=simplify(ast,steps),seq=[];let budget={steps:0};for(let n=start;n<start+count;n++){try{seq.push(evaluate(ast,n,budget))}catch(e){throw Error('n = '+n+': '+e.message)}}let g=skipGuess?null:guess(seq,start),gf=start===0?generating(simple,g,seq):null;if(start===0){let extra={'binomial((2*n),n)':'1 / √(1 − 4x)','(2^n)':'1 / (1 − 2x)'};if(extra[format(simple)])gf={kind:'通常母関数',formula:extra[format(simple)],proven:true}}if(start===0&&!gf&&!skipGuess)gf=algebraicGF(seq);return {source,ast,simple,steps,seq,g,gf,start,count}}
root.SigmaEngine={Q,q,freeVariables,parse,evaluate,format,fn,CATALOG,analyze,recurrenceText,polyText,polyEval,simplify,nullVector,algebraicGF};if(typeof module!=='undefined')module.exports=root.SigmaEngine;
})(typeof globalThis!=='undefined'?globalThis:this);
