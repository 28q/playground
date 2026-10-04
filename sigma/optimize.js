(function(root){
'use strict';
const E=root.SigmaEngine||(typeof require!=='undefined'?require('./engine.js'):null),{q,format,evaluate,parse,polyEval}=E;
const MOD=998244353n, norm=x=>(x%MOD+MOD)%MOD;
function mpow(a,n){let r=1n;for(a=norm(a);n;n>>=1n,a=a*a%MOD)if(n&1n)r=r*a%MOD;return r}
function residue(x){x=q(x);let b=norm(x.b);if(!b)throw Error('分母が 998244353 の倍数です');return norm(x.a)*mpow(b,MOD-2n)%MOD}
const num=x=>({type:'num',value:q(x)}),v=name=>({type:'var',name}),bin=(op,left,right)=>({type:'bin',op,left,right}),call=(name,...args)=>({type:'call',name,args});
const z=()=>q(0),one=()=>q(1);
const trim=a=>{while(a.length>1&&a[a.length-1].zero())a.pop();return a.length?a:[z()]};
const padd=(a,b)=>trim(Array.from({length:Math.max(a.length,b.length)},(_,i)=>(a[i]||z()).add(b[i]||z())));
const pscale=(a,b)=>trim(a.map(x=>x.mul(b)));
function pmul(a,b){let c=Array.from({length:a.length+b.length-1},z);for(let i=0;i<a.length;i++)for(let j=0;j<b.length;j++)c[i+j]=c[i+j].add(a[i].mul(b[j]));return trim(c)}
function pshift(a,k){let r=[z()],p=[one()];for(let c of a){r=padd(r,pscale(p,c));p=pmul(p,[q(k),one()])}return r}
function pdivmod(a,b){a=trim(a.slice());b=trim(b.slice());if(b.every(x=>x.zero()))throw Error('polynomial zero divisor');let out=Array.from({length:Math.max(1,a.length-b.length+1)},z);while(a.length>=b.length&&!a.every(x=>x.zero())){let d=a.length-b.length,c=a[a.length-1].div(b[b.length-1]);out[d]=c;for(let j=0;j<b.length;j++)a[d+j]=a[d+j].sub(c.mul(b[j]));trim(a)}return [trim(out),a]}
function pgcd(a,b){while(!b.every(x=>x.zero()))[a,b]=[b,pdivmod(a,b)[1]];return pscale(a,one().div(a[a.length-1]))}
function rat(a,b=[one()]){a=trim(a);b=trim(b);if(a.every(x=>x.zero()))return {a:[z()],b:[one()]};let g=pgcd(a,b);a=pdivmod(a,g)[0];b=pdivmod(b,g)[0];let c=b[b.length-1];return {a:pscale(a,one().div(c)),b:pscale(b,one().div(c))}}
const radd=(x,y)=>rat(padd(pmul(x.a,y.b),pmul(y.a,x.b)),pmul(x.b,y.b));
const rmul=(x,y)=>rat(pmul(x.a,y.a),pmul(x.b,y.b));
const rneg=x=>({a:pscale(x.a,q(-1)),b:x.b});
const rdiv=(x,y)=>rat(pmul(x.a,y.b),pmul(x.b,y.a));

// Sparse multivariate polynomial arithmetic; all coefficients remain rational.
function monomial(exps){return Object.entries(exps).filter(([,k])=>k).sort(([a],[b])=>a.localeCompare(b)).map(([x,k])=>x+':'+k).join(',')}
function exponents(key){return Object.fromEntries(key?key.split(',').map(x=>{let [v,k]=x.split(':');return [v,+k]}):[])}
const constant=x=>new Map([['',q(x)]]);
function madd(a,b){let c=new Map(a);for(let [key,x] of b)c.set(key,(c.get(key)||z()).add(x));for(let [k,c0] of c)if(c0.zero())c.delete(k);return c}
function mscale(a,b){return new Map([...a].map(([k,c])=>[k,c.mul(b)]).filter(([,c])=>!c.zero()))}
function mmul(a,b){let c=new Map();for(let [ka,ca] of a)for(let [kb,cb] of b){let e=exponents(ka);for(let [x,k] of Object.entries(exponents(kb)))e[x]=(e[x]||0)+k;if(Object.values(e).reduce((s,k)=>s+k,0)>24)throw Error('polynomial degree budget');let key=monomial(e);c.set(key,(c.get(key)||z()).add(ca.mul(cb)));if(c.size>400)throw Error('polynomial term budget')}return c}
function mpower(a,n){let r=constant(1);for(let i=0;i<n;i++)r=mmul(r,a);return r}
function msubst(a,x,b){let r=new Map();for(let [key,c] of a){let e=exponents(key),k=e[x]||0;delete e[x];r=madd(r,mmul(new Map([[monomial(e),c]]),mpower(b,k)))}return r}
function nonnegative(a,domain){return [...a].every(([key,c])=>c.a>=0n&&Object.keys(exponents(key)).every(x=>domain.has(x)))}
function mpolynomial(ast,domain=new Set(['n'])){
 if(ast.type==='num')return constant(ast.value);if(ast.type==='var')return new Map([[ast.name+':1',one()]]);
 if(ast.type==='neg'){let p=mpolynomial(ast.x,domain);return p&&mscale(p,-1)}
 if(ast.type==='bin'){let l=mpolynomial(ast.left,domain),r=mpolynomial(ast.right,domain);if(!l||!r)return null;if(ast.op==='+')return madd(l,r);if(ast.op==='-')return madd(l,mscale(r,-1));if(ast.op==='*')return mmul(l,r);if(ast.op==='/'&&r.size===1&&r.has('')&&!r.get('').zero())return mscale(l,one().div(r.get('')));if(ast.op==='^'&&r.size===1&&r.has('')&&r.get('').b===1n&&r.get('').a>=0n&&r.get('').a<=12n)return mpower(l,Number(r.get('').a));return null}
 if(ast.name==='binomial'&&ast.args[1].type==='num'){let k=ast.args[1].value;if(k.b!==1n||k.a<0n||k.a>12n)return null;let p=mpolynomial(ast.args[0],domain);if(!p)return null;let r=constant(1);for(let j=0;j<Number(k.a);j++)r=mscale(mmul(r,madd(p,constant(-j))),one().div(j+1));return r}
 if(ast.name==='sum'){
  let [body,x,l,h]=ast.args,lo=mpolynomial(l,domain),hi=mpolynomial(h,domain);if(!lo||!hi||[...lo.values(),...hi.values()].some(c=>c.b!==1n)||!nonnegative(madd(madd(hi,mscale(lo,-1)),constant(1)),domain))return null;
  let innerDomain=new Set(domain);innerDomain.delete(x.name);if(nonnegative(lo,domain))innerDomain.add(x.name);
  let p=mpolynomial(body,innerDomain);if(!p)return null;let r=new Map();
  for(let [key,c] of p){let e=exponents(key),d=e[x.name]||0;delete e[x.name];let anti=new Map();for(let j=0;j<=d;j++){let coefficient=E.fn('binomial',[q(d+1),q(j)]).mul(E.fn('bernoulli',[q(j)])).div(d+1);anti=madd(anti,new Map([[x.name+':'+(d+1-j),coefficient]]))}let definite=madd(msubst(anti,x.name,madd(hi,constant(1))),mscale(msubst(anti,x.name,lo),-1));r=madd(r,mmul(new Map([[monomial(e),c]]),definite));}
  return r;
 }return null;
}
function polynomial(ast){try{let p=mpolynomial(ast);if(!p||[...p.keys()].some(k=>Object.keys(exponents(k)).some(x=>x!=='n')))return null;let a=[z()];for(let [key,c] of p){let k=exponents(key).n||0;while(a.length<=k)a.push(z());a[k]=a[k].add(c)}return trim(a)}catch{return null}}
function polynomialAST(p){let a=null;p.forEach((c,k)=>{if(c.zero())return;let term=k===0?num(c):k===1?v('n'):bin('^',v('n'),num(k));if(k&&!c.eq(1))term=bin('*',num(c),term);a=a?bin('+',a,term):term});return a||num(0)}
function substitute(a,name,b){
 if(a.type==='var')return a.name===name?b:a;
 if(a.type==='neg')return {...a,x:substitute(a.x,name,b)};
 if(a.type==='bin')return {...a,left:substitute(a.left,name,b),right:substitute(a.right,name,b)};
 if(a.type==='call'){
  if(a.name==='sum'||a.name==='prod'){
   let [body,index,lo,hi]=a.args;
   // The binder scopes only the body. Rename it before introducing a free variable.
   if(index.name!==name&&free(body,name)&&free(b,index.name)){
    let i=0,fresh;do{fresh='_s'+i++}while(fresh===name||free(body,fresh)||free(b,fresh));
    body=substitute(body,index.name,v(fresh));index=v(fresh);
   }
   return {...a,args:[index.name===name?body:substitute(body,name,b),index,substitute(lo,name,b),substitute(hi,name,b)]};
  }
  return {...a,args:a.args.map(x=>substitute(x,name,b))};
 }
 return a;
}
function free(a,x){if(a.type==='var')return a.name===x;if(a.type==='num')return false;if(a.type==='neg')return free(a.x,x);if(a.type==='bin')return free(a.left,x)||free(a.right,x);if(a.type==='call'&&(a.name==='sum'||a.name==='prod'))return free(a.args[2],x)||free(a.args[3],x)||(a.args[1].name!==x&&free(a.args[0],x));return a.args.some(y=>free(y,x))}
function constantValue(a){if(free(a,'n'))return null;try{return evaluate(a,0)}catch{return null}}
const same=(a,b)=>format(a)===format(b),isnum=(a,n)=>a.type==='num'&&a.value.eq(n),isn=a=>a.type==='var'&&a.name==='n';
function flatten(a){return a.type==='bin'&&a.op==='*'?[...flatten(a.left),...flatten(a.right)]:[a]}
function product(a){return a.reduce((x,y)=>x?bin('*',x,y):y,null)||num(1)}
function simplify(ast,steps=[]){
 let a=E.simplify(ast,steps),p=polynomial(a);if(p){let b=polynomialAST(p);if(!same(a,b)){steps.push(format(a)+' → '+format(b)+' （多項式の展開・有限和の厳密消去）');return b}return a}
 if(a.type==='bin')a={...a,left:simplify(a.left,steps),right:simplify(a.right,steps)};
 if(a.type==='neg')a={...a,x:simplify(a.x,steps)};
 // Only whole sums with n free are transformed here; bound-variable domains stay explicit.
 if(a.type!=='call'||a.name!=='sum')return a;
 let [body,k,lo,hi]=a.args;if(k.name==='n'||!isnum(lo,0)||!isn(hi))return a;let b=null;
 if(body.type==='call'&&body.name==='binomial'&&same(body.args[0],k)&&body.args[1].type==='num'&&body.args[1].value.b===1n&&body.args[1].value.a>=0n)b=call('binomial',bin('+',v('n'),num(1)),num(body.args[1].value.add(1)));
 if(body.type==='bin'&&body.op==='^'&&same(body.right,k)){let c=constantValue(body.left);if(c)b=c.eq(1)?bin('+',v('n'),num(1)):bin('/',bin('-',bin('^',num(c),bin('+',v('n'),num(1))),num(1)),num(c.sub(1)))}
 let factors=flatten(body),choose=factors.find(f=>f.type==='call'&&f.name==='binomial'&&isn(f.args[0])&&same(f.args[1],k));
 if(choose){let rest=factors.filter(f=>f!==choose),left=q(1),right=q(1),valid=true;for(let f of rest){if(f.type!=='bin'||f.op!=='^'){valid=false;break}let c=constantValue(f.left);if(!c){valid=false;break}if(same(f.right,k))left=left.mul(c);else if(f.right.type==='bin'&&f.right.op==='-'&&isn(f.right.left)&&same(f.right.right,k))right=right.mul(c);else {valid=false;break}}if(valid)b=bin('^',num(left.add(right)),v('n'));}
 if(factors.length===2&&factors.every(f=>f.type==='call'&&f.name==='binomial')){for(let [f,g] of [[factors[0],factors[1]],[factors[1],factors[0]]]){let nk=bin('-',v('n'),k);if(same(f.args[1],k)&&same(g.args[1],nk)&&!free(f.args[0],k.name)&&!free(g.args[0],k.name)){let x=polynomial(f.args[0]),y=polynomial(g.args[0]);if(x&&y&&[...x,...y].every(c=>c.b===1n&&c.a>=0n))b=call('binomial',bin('+',f.args[0],g.args[0]),v('n'))}}}
 if(body.type==='bin'&&body.op==='-'&&same(body.left,substitute(body.right,k.name,bin('+',k,num(1)))))b=bin('-',substitute(body.right,k.name,bin('+',v('n'),num(1))),substitute(body.right,k.name,num(0)));
 if(b){b=E.simplify(b);steps.push(format(a)+' → '+format(b)+' （恒等式）');return b}return a;
}

// Berlekamp–Massey in F_998244353, followed by rational lifting and exact checks.
function berlekampMassey(sequence){let s=sequence.map(residue),C=[1n],B=[1n],L=0,m=1,b=1n;for(let n=0;n<s.length;n++){let d=s[n];for(let j=1;j<=L;j++)d=(d+(C[j]||0n)*s[n-j])%MOD;if(!d){m++;continue}let T=C.slice(),coef=d*mpow(b,MOD-2n)%MOD;while(C.length<B.length+m)C.push(0n);for(let j=0;j<B.length;j++)C[j+m]=norm(C[j+m]-coef*B[j]);if(2*L<=n){L=n+1-L;B=T;b=d;m=1}else m++}return {order:L,coefficients:C.slice(0,L+1).map(Number)}}
function verifyRec(g,seq,start,from=0){for(let i=from;i<seq.length-g.order;i++){let s=z();for(let j=0;j<=g.order;j++)s=s.add(polyEval(g.p[j],start+i).mul(seq[i+j]));if(!s.zero())return false}return true}
function fit(seq,start,order,degree){let train=seq.length-8,cols=(order+1)*(degree+1);if(train-order<cols+2)return null;let rows=[];for(let i=0;i<Math.min(train-order,cols+3);i++){let row=[];for(let j=0;j<=order;j++)for(let k=0;k<=degree;k++)row.push(seq[i+j].mul(q(start+i).a**BigInt(k)));rows.push(row)}let x=E.nullVector(rows,cols);if(!x)return null;let p=Array.from({length:order+1},(_,j)=>trim(x.slice(j*(degree+1),(j+1)*(degree+1))));if(p[order].every(c=>c.zero()))return null;let g={order,degree,p,train,holdout:8,method:'有理数係数の線形方程式'};return verifyRec(g,seq,start,0)?g:null}
const KNOWN={
 fibonacci:{p:[[-1],[-1],[1]],init:[0,1]},lucas:{p:[[-1],[-1],[1]],init:[2,1]},factorial:{p:[[-1,-1],[1]],init:[1]},
 catalan:{p:[[-2,-4],[2,1]],init:[1]},central:{p:[[-2,-4],[1,1]],init:[1]},
 motzkin:{p:[[-3,-3],[-5,-2],[4,1]],init:[1,1]},schroeder:{p:[[0,1],[-9,-6],[3,1]],init:[1,2]},
 derangement:{p:[[-1,-1],[-1,-1],[1]],init:[1,0]},harmonic:{p:[[1,1],[-3,-2],[2,1]],init:[0,1]}
};
function leaf(a){if(a.type==='call'&&a.args.length===1&&isn(a.args[0])&&KNOWN[a.name])return {name:a.name,...KNOWN[a.name]};if(a.type==='call'&&a.name==='binomial'&&isn(a.args[1])){let p=polynomial(a.args[0]);if(p&&p.length===2&&p[0].zero()&&p[1].eq(2))return {name:'central',...KNOWN.central}}return null}
function structure(a){let p=polynomial(a);if(p)return {type:'C-finite',order:p.length,degree:0,reason:'n の '+(p.length-1)+' 次多項式。差分の '+p.length+' 回反復で 0。'};let l=leaf(a);if(l)return {type:l.p.every(p=>p.length===1)?'C-finite':'P-recursive',order:l.p.length-1,degree:Math.max(...l.p.map(p=>p.length-1)),reason:l.name+' の既知の漸化式。'};
 if(a.type==='bin'&&a.op==='^'&&isn(a.right)&&constantValue(a.left))return {type:'C-finite',order:1,degree:0,reason:'定数の n 乗。1階の定数係数漸化式。'};
 if(a.type==='neg')return structure(a.x);
 if(a.type==='bin'&&['+','-','*'].includes(a.op)){let x=structure(a.left),y=structure(a.right);if(x.order&&y.order){let order=a.op==='*'?x.order*y.order:x.order+y.order;return {type:x.type==='C-finite'&&y.type==='C-finite'?'C-finite':'P-recursive',order:Math.min(order,256),degree:null,reason:(a.op==='*'?'項ごとの積':'和・差')+'の閉性から階数 ≤ '+order+'。多項式係数の次数は別途探索。'}}}
 if(a.type==='call'&&['bell','bernoulli','euler','partition'].includes(a.name)&&isn(a.args[0]))return {type:'母関数',reason:'既知の母関数を FPS 演算で展開可能。有限階の P-recursive 性は仮定しない。'};
 return {type:'未分類',reason:'有限サンプルから探索。見つからないことは不存在の証明ではない。'};
}
function shiftVectors(ast,max){let p=polynomial(ast);if(p)return Array.from({length:max+1},(_,j)=>[rat(pshift(p,j))]);let l=leaf(ast);if(l){let p=l.p.map(p=>p.map(q)),r=p.length-1,vs=[];for(let j=0;j<=max;j++){if(j<r){vs.push(Array.from({length:r},(_,i)=>rat([q(i===j?1:0)])));continue}let s=j-r,out=Array.from({length:r},()=>rat([z()]));for(let i=0;i<r;i++)for(let k=0;k<r;k++)out[k]=radd(out[k],rmul(rat(pscale(pshift(p[i],s),-1),pshift(p[r],s)),vs[s+i][k]));vs.push(out)}return vs}
 if(ast.type==='bin'&&ast.op==='^'&&isn(ast.right)){let c=constantValue(ast.left);if(c&&!c.zero())return Array.from({length:max+1},(_,j)=>[rat([new E.Q(c.a**BigInt(j),c.b**BigInt(j))])])}
 if(ast.type==='neg'){let x=shiftVectors(ast.x,max);return x&&x.map(row=>row.map(rneg))}
 if(ast.type==='bin'&&['+','-','*'].includes(ast.op)){let x=shiftVectors(ast.left,max),y=shiftVectors(ast.right,max);if(!x||!y||x[0].length*y[0].length>32)return null;return x.map((row,j)=>ast.op==='*'?row.flatMap(a=>y[j].map(b=>rmul(a,b))):[...row,...y[j].map(b=>ast.op==='-'?rneg(b):b)])}return null;
}
function certify(ast,g){try{let vectors=shiftVectors(ast,g.order);if(!vectors)return false;for(let k=0;k<vectors[0].length;k++){let s=rat([z()]);for(let j=0;j<=g.order;j++)s=radd(s,rmul(rat(g.p[j]),vectors[j][k]));if(!s.a.every(c=>c.zero()))return false}return true}catch{return false}}

function rationalGF(ast){let p=polynomial(ast);if(p){let d=p.length,seq=Array.from({length:d},(_,n)=>polyEval(p,n)),den=[one()];for(let k=0;k<d;k++)den=pmul(den,[one(),q(-1)]);return rat(pmul(seq,den).slice(0,d),den)}let l=leaf(ast);if(l&&l.p.every(p=>p.length===1)){let den=l.p.map(p=>q(p[0])).reverse(),seq=l.init.map(q);return rat(pmul(seq,den).slice(0,seq.length),den)}
 if(ast.type==='bin'&&ast.op==='^'&&isn(ast.right)){let c=constantValue(ast.left);if(c)return rat([one()],[one(),c.neg()])}
 if(ast.type==='neg'){let x=rationalGF(ast.x);return x&&rneg(x)}
 if(ast.type==='bin'&&['+','-'].includes(ast.op)){let x=rationalGF(ast.left),y=rationalGF(ast.right);if(x&&y)return radd(x,ast.op==='-'?rneg(y):y)}
 if(ast.type==='bin'&&['*','/'].includes(ast.op)){let c=constantValue(ast.right),x=rationalGF(ast.left);if(c&&x&&(ast.op!=='/'||!c.zero()))return rmul(x,rat([ast.op==='*'?c:one().div(c)]));if(ast.op==='*'){c=constantValue(ast.left);x=rationalGF(ast.right);if(c&&x)return rmul(x,rat([c]))}}
 let s=seriesShape(ast);if(s&&s.kind==='prefix'){let x=rationalGF(s.child);if(x)return rat(x.a,pmul(x.b,[one(),q(-1)]))}if(s&&s.kind==='convolution'){let x=rationalGF(s.left),y=rationalGF(s.right);if(x&&y)return rmul(x,y)}return null;
}
function seriesShape(a){if(a.type!=='call'||a.name!=='sum')return null;let [body,k,lo,hi]=a.args;if(k.name==='n'||!isnum(lo,0)||!isn(hi))return null;let factors=flatten(body),bi=factors.find(x=>x.type==='call'&&x.name==='binomial'&&isn(x.args[0])&&same(x.args[1],k));if(bi)factors=factors.filter(x=>x!==bi);
 let left=[],right=[],nk=bin('-',v('n'),k);for(let f of factors){if(!free(f,'n'))left.push(substitute(f,k.name,v('n')));else{let code=format(f),target=format(nk);if(!code.includes(target))return null;let temp=replaceExact(f,nk,v('_i'));if(free(temp,'n')||free(temp,k.name))return null;right.push(substitute(temp,'_i',v('n')))}}
 if(!right.length)return {kind:bi?'binomial-transform':'prefix',child:product(left)};
 return {kind:bi?'binomial-convolution':'convolution',left:product(left),right:product(right)};
}
function replaceExact(a,from,to){if(same(a,from))return to;if(a.type==='bin')return {...a,left:replaceExact(a.left,from,to),right:replaceExact(a.right,from,to)};if(a.type==='neg')return {...a,x:replaceExact(a.x,from,to)};if(a.type==='call'&&!['sum','prod'].includes(a.name))return {...a,args:a.args.map(x=>replaceExact(x,from,to))};return a}
function seriesPlan(a,depth=0){if(depth>12)return null;let p=polynomial(a);if(p)return {kind:'polynomial',p,proven:true};let l=leaf(a);if(l)return {kind:'known-recurrence',...l,proven:true};let gf=rationalGF(a);if(gf)return {kind:'rational',P:gf.a,Q:gf.b,proven:true};
 if(a.type==='call'&&a.args.length===1&&isn(a.args[0])&&['bell','bernoulli','euler','partition'].includes(a.name))return {kind:a.name,proven:true};
 if(a.type==='neg'){let x=seriesPlan(a.x,depth+1);return x&&{kind:'scale',child:x,scale:q(-1),proven:true}}
 if(a.type==='bin'&&['+','-','*'].includes(a.op)){let left=seriesPlan(a.left,depth+1),right=seriesPlan(a.right,depth+1);if(left&&right)return {kind:a.op==='*'?'hadamard':a.op==='+'?'add':'subtract',left,right,proven:true}}
 if(a.type==='bin'&&a.op==='/'){let c=constantValue(a.right),x=seriesPlan(a.left,depth+1);if(c&&!c.zero()&&x)return {kind:'scale',child:x,scale:one().div(c),proven:true}}
 let s=seriesShape(a);if(s){if(s.child){let child=seriesPlan(s.child,depth+1);return child&&{kind:s.kind,child,proven:true}}let left=seriesPlan(s.left,depth+1),right=seriesPlan(s.right,depth+1);if(left&&right)return {kind:s.kind,left,right,proven:true}}return null;
}
function complexity(plan,target){if(plan.kind==='polynomial')return target==='nth'?'O(d)':'O(Nd)';if(plan.kind==='rational')return target==='nth'?'O(M(r) log N)':'O(M(N) + M(r))';if(plan.kind==='p-recursive')return target==='nth'&&plan.order<=4?'Õ(r³ √(dN))（多項式行列BSGS）':'O(Nr(d+1) + (N/4096) log p)';if(plan.kind==='known-recurrence')return 'O(Nr(d+1) + (N/4096) log p)';return target==='nth'?'O(M(N) + N log N)（Nまで展開して抽出）':'O(M(N) + N log N)';}
function planFor(result){let a=result.simple,poly=polynomial(a);if(poly)return {kind:'polynomial',p:poly,start:0,proven:true,method:'多項式の厳密消去 → Horner 法'};let rational=rationalGF(a);if(rational)return {kind:'rational',P:rational.a,Q:rational.b,start:0,proven:true,method:'有理型母関数 → Bostan–Mori / FPS inverse'};
 let sp=seriesPlan(a),g=result.g;
 if(g&&g.degree===0){let den=g.p.map(p=>p[0]).reverse(),seq=result.seq.slice(0,g.order);return {kind:'rational',P:trim(pmul(seq,den).slice(0,g.order)),Q:den,start:result.start,proven:!!g.proven,method:'Berlekamp–Massey / 漸化式 → Bostan–Mori / FPS inverse'}}
 if(sp&&['bell','bernoulli','euler','partition','convolution','binomial-convolution','binomial-transform','prefix'].includes(sp.kind))return {...sp,start:0,method:'母関数の構造 → NTT / FPS 演算'};
 if(g)return {kind:'p-recursive',p:g.p,order:g.order,degree:g.degree,seeds:result.seq,start:result.start,proven:!!g.proven,method:g.order<=4?'P-recursive → 多項式行列BSGS / 列は専用漸化式':'P-recursive → 専用漸化式（分母を一括反転）'};
 if(sp)return {...sp,start:0,method:'既知の数列と閉性 → 専用の係数列計算'};
 if(result.gf&&result.gf.algebraic&&result.gf.coefficients)return {kind:'algebraic',polynomials:result.gf.coefficients,initial:result.seq[0],start:0,proven:false,method:'代数的母関数の推測 → FPS Newton 法'};
 return null;
}
const originalAnalyze=E.analyze;
function analyze(source,options={}){
 let res=originalAnalyze(source,{...options,skipGuess:true});res.simple=simplify(res.simple,res.steps);res.structure=structure(res.simple);let orderCap=Math.max(1,Math.min(16,Number(options.maxOrder)||8)),degreeCap=Math.max(0,Math.min(8,options.maxDegree===undefined?4:Number(options.maxDegree)));let train=res.seq.length-8;
 // A proven C-finite upper bound takes priority over the generic search cap when samples suffice.
 if(res.structure.type==='C-finite'&&res.structure.order<=16)orderCap=Math.max(orderCap,res.structure.order);
 let g=null,bm=null;let known=leaf(res.simple),knownFPS=res.structure.type==='母関数';if(known){g={order:known.p.length-1,degree:Math.max(...known.p.map(p=>p.length-1)),p:known.p.map(p=>p.map(q)),train,holdout:8,method:'既知の漸化式'};if(g.order>orderCap||g.degree>degreeCap)g=null}try{if(!g&&!knownFPS){bm=berlekampMassey(res.seq.slice(0,train));if(bm.order>0&&bm.order<=orderCap&&train>=2*bm.order+2){g=fit(res.seq,res.start,bm.order,0);if(g)g.method='Berlekamp–Massey (mod 998244353) → 有理数への持ち上げ → 厳密検算'}}}catch{}
 if(!g&&res.g&&res.g.order<=orderCap&&res.g.degree<=degreeCap)g=res.g;
 if(!g&&!knownFPS){let attempts=0;outer:for(let cost=4;cost<=Math.min(36,train-3);cost++)for(let d=1;d<=degreeCap;d++)for(let r=1;r<=Math.min(orderCap,res.structure.order||orderCap);r++){if((r+1)*(d+1)!==cost||train-r<cost+2)continue;if(++attempts>36)break outer;g=fit(res.seq,res.start,r,d);if(g)break outer}}
 if(g){g.proven=certify(res.simple,g);if(!g.proven&&g.degree===0&&res.structure.type==='C-finite'&&res.seq.length-g.order>=res.structure.order){g.proven=true;g.certificate='既知の定数係数漸化式の階数上界と、十分な連続初期値の厳密一致'}else if(g.proven)g.certificate='既知の漸化式をシフト・加算・項別乗算し、有理関数の恒等式として検証';}
 res.g=g;if(!res.gf&&res.start===0)res.gf=E.algebraicGF(res.seq);res.search={maxOrder:orderCap,maxDegree:degreeCap,bmOrder:bm&&bm.order,holdout:8,train};res.fast=planFor(res);
 if(res.fast&&res.fast.kind==='rational'&&res.fast.start===0)res.gf={kind:'通常母関数'+(res.fast.proven?'':'（推測）'),formula:'('+E.polyText(res.fast.P,'x')+') / ('+E.polyText(res.fast.Q,'x')+')',proven:res.fast.proven};
 if(res.fast){res.fast.nthComplexity=complexity(res.fast,'nth');res.fast.prefixComplexity=complexity(res.fast,'prefix');res.fast.modulus=Number(MOD)}
 return res;
}
function multivariate(ast,variables=['n','m']){try{let p=mpolynomial(ast,new Set(variables));if(!p)return null;return [...p].map(([key,c])=>({powers:exponents(key),coefficient:c}));}catch{return null}}
const API={multivariate,substitute,free,analyze,simplify,polynomial,structure,berlekampMassey,certify,seriesPlan,rationalGF,residue,planFor,MOD:Number(MOD)};
root.SigmaOptimize=API;E.analyze=analyze;if(typeof module!=='undefined')module.exports=API;
})(typeof globalThis!=='undefined'?globalThis:this);
