const E=require('./engine'),O=require('./optimize'),R=require('./rewrite'),M=require('./modular'),Q=require('./queries'),C=require('./query-cpp'),F=require('./fast-cpp'),assert=require('assert/strict'),fs=require('fs'),cp=require('child_process'),path=require('path'),os=require('os');let checks=0;
const expressions=[
 'sum(C(n,k)*C(m,r-k),k,0,n)', 'sum(C(m,r-k)*C(n,k),k,0,min(n,r))',
 'sum(C(n,k+a)*C(m,k+b),k,-a,n-a)', 'sum(C(n,n-k)*C(m,r-k),k,0,n)',
 'sum(3*C(n,k)*C(m,r-k),k,0,n)', 'sum(k^5*C(n,k),k,0,n)',
 'sum((k^3+a*k+2)*C(n,k),k,0,n)', 'sum(C(n,k),k,0,min(m,n))',
 'sum(k*factorial(k),k,0,n)', 'sum((k^2+3*k+1)*2^k,k,0,n)',
 'sum((k^2+3*k+1)*factorial(k+1),k,0,n)', 'sum(k*factorial(k),k,a,a+n)',
 'sum(sum(C(n,i)*C(m,j)*min(i,j),j,0,m),i,0,n)',
 'sum(sum(i+2*j,j,0,i),i,0,n)', 'sum(sum(fibonacci(i),i,j,n),j,0,n)',
 'sum(sum(factorial(j+1)-factorial(j),j,0,i),i,0,n)',
 'sum(k^2,k,3,n)', 'sum(k*factorial(k),k,3,n)',
 'sum(C(-2,k),k,0,n)', 'sum(sum(C(i,k),k,0,i),i,0,n)',
 'sum(sum(n+j,j,0,n),n,0,m)', 'sum(C(n,k)*factorial(n-k),k,0,min(n,m))',
 'sum(sum(gcd(i,j),j,0,m),i,0,n)', 'sum(sum(gcd(i,j),j,1,m),i,1,n)',
 'sum(fibonacci(k),k,a+1,b+2)', 'prod(k-2,k,a,b+c)',
 'sum(factorial(-1)*k,k,1,0)', 'sum(1/k-1/(k-1),k,1,0)'
];
for(const source of expressions){const ast=E.parse(source,null),names=E.freeVariables(ast),variants=R.variants(ast);for(let seed=0;seed<32;seed++){let env=Object.fromEntries(names.map((x,i)=>[x,(seed*(2*i+1)+(seed>>i))%6])),want=Q.evaluateAt(ast,env);for(const candidate of variants){assert.equal(String(Q.evaluateAt(candidate.ast,env)),String(want),source+' '+JSON.stringify(env)+' => '+E.format(candidate.ast));checks++;}}}
assert(R.normalize(E.parse(expressions[0],null)).steps.some(x=>x.reason.includes('Vandermonde')));
assert(R.normalize(E.parse(expressions[8],null)).steps.some(x=>x.reason.includes('多項式証明書')));
assert(!R.normalize(E.parse('sum(C(n,k)*factorial(n-k),k,0,min(n,m))',null)).steps.some(x=>x.reason.includes('上限を拡大')));
assert(!R.normalize(E.parse('sum(k*factorial(k),k,3,n)')).steps.some(x=>x.reason.includes('証明書')),'Cannot telescope across potentially empty interval');
for(const d of [3,23,63]){const seq=Array.from({length:2*(d+1)+14},(_,n)=>E.q(BigInt(n+1)**BigInt(d))),g=M.fit(seq,0,1,d);assert(g);assert.equal(g.unknowns,2*(d+1));assert(g.primes.length>=3);assert.equal(g.proven,false);const extra=Array.from({length:seq.length+15},(_,n)=>E.q(BigInt(n+1)**BigInt(d)));assert(M.verify(g,extra,0));const changed=seq.slice();changed[changed.length-1]=changed.at(-1).add(1);assert.equal(M.fit(changed,0,1,d),null);checks+=6;}
const rationalSeq=Array.from({length:24},(_,n)=>new E.Q(1,2*n+1)),rationalGuess=M.fit(rationalSeq,0,1,1);assert(rationalGuess);assert.equal(rationalGuess.reconstruction,'有理数復元');assert(M.verify(rationalGuess,Array.from({length:40},(_,n)=>new E.Q(1,2*n+1)),0));checks+=3;
const high=E.analyze('(n+1)^249',{count:512,maxOrder:1,maxDegree:249,maxUnknowns:512,exactSearch:true});assert(high.g);assert.equal(high.g.unknowns,500);assert(high.g.primes.length>=3);assert.equal(high.g.proven,false);assert(F.generate(high).includes('NOT MATHEMATICALLY PROVEN'));checks+=5;
const univariate=['sum(k^3*C(n,k),k,0,n)','sum((-1)^k*C(n,k)*k^4,k,0,n)','sum((-1)^(n-k)*C(n,k)*(k^4+2),k,0,n)','sum(k^2*3^k,k,0,n)','sum(k*factorial(k),k,0,n)'];let cases=[];
for(const source of univariate){const r=E.analyze(source);assert(r.fast?.proven,source);for(let n=0;n<28;n++){assert.equal(String(E.evaluate(r.simple,n)),String(E.evaluate(r.ast,n)),source);checks++;}if(r.differential){const df=r.differential;for(let n=0;n<35;n++){let lhs=E.q(0);for(const t of df.terms)if(n>=t.xPower)lhs=lhs.add(E.polyEval(t.theta,n-t.xPower).mul(E.evaluate(r.ast,n-t.xPower)));assert.equal(String(lhs),String(df.rhs[n]||E.q(0)));checks++;}}cases.push({source,r,kind:'sequence',points:Array.from({length:30},(_,i)=>i)});}
for(const source of [expressions[0],expressions[2],expressions[4],expressions[6],expressions[7],expressions[12],...expressions.slice(22,26)]){const r=Q.analyze(source);assert(r.choices.length,source);for(const plan of r.choices){const points=Array.from({length:35},(_,i)=>r.variables.map((_,j)=>(i*(j*2+1)+(i>>j))%8));cases.push({source,r,plan,kind:'query',points});}}
const dir=fs.mkdtempSync(path.join(os.tmpdir(),'sigma-math-'));let cpp='#include <bits/stdc++.h>\n#define SIGMA_NO_MAIN\n';cases.forEach((t,i)=>{const code=t.kind==='sequence'?F.generate(t.r):C.generate(t.r,t.plan.id);cpp+='namespace test'+i+' {\n'+code.replace(/^#include .*$/gm,'')+'\n}\n'});cpp+='int main(){\n';cases.forEach((t,i)=>{if(t.kind==='sequence')cpp+='for(int n=0;n<30;n++)std::cout<<"'+i+' "<<n<<" "<<test'+i+'::coefficient(n).v<<"\\n";\n';else cpp+='{std::vector<test'+i+'::Query> qs={'+t.points.map((row,j)=>'{'+[...row,j].join(',')+'}').join(',')+'};auto a=test'+i+'::solve_queries(qs);for(size_t j=0;j<a.size();j++)std::cout<<"'+i+' "<<j<<" "<<a[j].v<<"\\n";}\n'});cpp+='}\n';fs.writeFileSync(path.join(dir,'test.cpp'),cpp);cp.execFileSync('g++',['-std=c++17','-O2',path.join(dir,'test.cpp'),'-o',path.join(dir,'test.exe')],{timeout:120000});const output=cp.execFileSync(path.join(dir,'test.exe'),[],{encoding:'utf8',timeout:120000});for(const line of output.trim().split(/\r?\n/)){const [i,j,value]=line.split(' '),t=cases[i],want=t.kind==='sequence'?E.evaluate(t.r.ast,Number(j)):Q.evaluateAt(t.r.ast,Object.fromEntries(t.r.variables.map((x,k)=>[x,t.points[j][k]])));assert.equal(value,String(O.residue(want)),t.source+' case '+j);checks++;}
fs.writeFileSync(path.join(dir,'high.cpp'),F.generate(high));cp.execFileSync('g++',['-std=c++17','-O2',path.join(dir,'high.cpp'),'-o',path.join(dir,'high.exe')],{timeout:120000});assert.equal(cp.execFileSync(path.join(dir,'high.exe'),[],{input:'600\n',encoding:'utf8'}).trim(),String(601n**249n%998244353n));checks++;
console.log('PASS mathematical rewrites / modular recovery / differential equations / generated C++:',checks,'C++ plans:',cases.length,'500 unknowns primes:',high.g.primes.length);console.log('Artifacts:',dir);
