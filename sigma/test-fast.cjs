const E=require('./engine.js'),O=require('./optimize.js'),C=require('./fast-cpp.js'),fs=require('fs'),cp=require('child_process'),assert=require('assert/strict'),path=require('path');
const testDir=fs.mkdtempSync(path.join(require('os').tmpdir(),'sigma-fast-test-'));
const artifact=name=>path.join(testDir,name);
const expressions=['n^5+2*n+1','sum(k^3,k,1,n)','sum(sum(j^2,j,0,k),k,0,n)','fibonacci(n)','lucas(n)','2^n+3^n','sum(fibonacci(k),k,0,n)','fibonacci(n)*lucas(n)','motzkin(n)*sum(binomial(n,k)^2,k,0,n)','catalan(n)','derangement(n)','schroeder(n)','harmonic(n)','bell(n)','bernoulli(n)','euler(n)','partition(n)','sum(bell(k)*fibonacci(n-k),k,0,n)','sum(binomial(n,k)*bell(k),k,0,n)','sum(binomial(n,k)*bell(k)*bell(n-k),k,0,n)','sum(binomial(n,k)^2*binomial(n+k,k)^2,k,0,n)','floor(n/3)','sum(catalan(k)*catalan(n-k),k,0,n)','binomial(2*n,n)/(n+1)'];
let results=expressions.map(source=>{let r=E.analyze(source,source==='binomial(2*n,n)/(n+1)'?{maxOrder:1,maxDegree:0}:{});console.log(r.fast&&r.fast.kind,source);assert(r.fast);return r});
let shifted=E.analyze(expressions[8],{start:8,count:64});results.push(shifted);
let harness='#include <algorithm>\n#include <cmath>\n#include <cstdint>\n#include <iostream>\n#include <stdexcept>\n#include <vector>\n#define SIGMA_NO_MAIN\n';
for(let i=0;i<results.length;i++)for(let target of ['nth','prefix']){let source=C.generate(results[i],target);assert(!source.includes('Q F('));assert(!source.includes('name=='));harness+='namespace test_'+i+'_'+target+' {\n'+source.replace(/^#include .*$/gm,'')+'\n}\n';}
harness+='int main(){\n';for(let i=0;i<results.length;i++){let start=results[i].fast.start||0;harness+='for(int n='+start+';n<85;n++)std::cout<<"'+i+' nth "<<n<<" "<<test_'+i+'_nth::coefficient(n).v<<"\\n";\n';harness+='{auto a=test_'+i+'_prefix::coefficients(84);for(int i=0;i<(int)a.size();i++)std::cout<<"'+i+' prefix "<<(i+'+start+')<<" "<<a[i].v<<"\\n";}\n'}
// Trigger BSGS, NTT, multipoint remainders and Newton iteration beyond naive cutoffs.
for(let i of [8,9,10,11,12,20,24])harness+='std::cout<<"'+i+' large 5000 "<<test_'+i+'_nth::coefficient(5000).v<<"\\n";\n';
harness+='std::cout<<"3 huge 1000000000000000000 "<<test_3_nth::coefficient(1000000000000000000ULL).v<<"\\n";\n';
for(let i of [13,14,15,16,17,18,19,22,23])harness+='{auto a=test_'+i+'_prefix::coefficients(1024);std::cout<<"'+i+' prefix-large 1024 "<<a.back().v<<"\\n";}\n';
harness+='}\n';fs.writeFileSync(artifact('fast-test.cpp'),harness);cp.execFileSync(process.env.CXX||'g++',['-std=c++17','-O2',artifact('fast-test.cpp'),'-o',artifact('fast-test.exe')],{stdio:'pipe'});let output=cp.execFileSync(path.resolve(artifact('fast-test.exe')),[],{encoding:'utf8',timeout:120000});fs.writeFileSync(artifact('fast-test-output.txt'),output);
let checks=0,P=998244353n,mod=x=>(x%P+P)%P,pow=(a,k)=>{let r=1n;for(;k;k>>=1n,a=a*a%P)if(k&1n)r=r*a%P;return r};
function reference(g,source,start,N){let r=g.order,a=Array.from({length:r},(_,i)=>O.residue(E.evaluate(E.parse(source),start+i)));for(let n=start;n<=N-r;n++){let s=0n;for(let j=0;j<r;j++)s=mod(s+O.residue(E.polyEval(g.p[j],n))*a[j]);let den=O.residue(E.polyEval(g.p[r],n)),next=mod(-s*pow(den,P-2n));a.shift();a.push(next)}return a[r-1]}
function fib(n){if(!n)return [0n,1n];let [a,b]=fib(n>>1n),c=mod(a*(2n*b-a)),d=mod(a*a+b*b);return n&1n?[d,mod(c+d)]:[c,d]}

let largeCache=null;
function largeReference(index,N){
 if(!largeCache){
  const fac=Array(N+2).fill(1n),invfac=Array(N+2).fill(1n);for(let i=1;i<=N+1;i++)fac[i]=fac[i-1]*BigInt(i)%P;invfac[N+1]=pow(fac[N+1],P-2n);for(let i=N+1;i>0;i--)invfac[i-1]=invfac[i]*BigInt(i)%P;
  const bell=[1n];let row=[1n];for(let n=1;n<=N;n++){let next=[row[row.length-1]];for(let k=1;k<=n;k++)next.push(mod(next[k-1]+row[k-1]));bell.push(next[0]);row=next;}
  const invert=f=>{let out=[pow(f[0],P-2n)];for(let n=1;n<=N;n++){let sum=0n;for(let k=1;k<=n;k++)sum=mod(sum+f[k]*out[n-k]);out.push(mod(-sum*out[0]))}return out;};
  let bern=invert(Array.from({length:N+1},(_,k)=>invfac[k+1])).map((x,k)=>x*fac[k]%P),euler=invert(Array.from({length:N+1},(_,k)=>k%2?0n:invfac[k])).map((x,k)=>x*fac[k]%P);
  let part=Array(N+1).fill(0n);part[0]=1n;for(let k=1;k<=N;k++)for(let n=k;n<=N;n++)part[n]=mod(part[n]+part[n-k]);
  let f=[0n,1n];for(let k=2;k<=N;k++)f[k]=mod(f[k-1]+f[k-2]);
  let conv=0n,trans=0n,bconv=0n;for(let k=0;k<=N;k++){let choose=fac[N]*invfac[k]%P*invfac[N-k]%P;conv=mod(conv+bell[k]*f[N-k]);trans=mod(trans+choose*bell[k]);bconv=mod(bconv+choose*bell[k]%P*bell[N-k]);}
  let cat=1n;for(let k=0;k<=N;k++)cat=cat*BigInt(4*k+2)%P*pow(BigInt(k+2),P-2n)%P;let catN=cat*BigInt(N+2)%P*pow(BigInt(4*N+2),P-2n)%P;
  largeCache={13:bell[N],14:bern[N],15:euler[N],16:part[N],17:conv,18:trans,19:bconv,22:cat,23:catN};
 }return largeCache[index];
}
for(let line of output.trim().split(/\r?\n/)){let [index,kind,n,value]=line.split(' '),r=results[+index],want;if(kind==='large')want=reference(r.g,r.source,r.start,+n);else if(kind==='huge')want=fib(BigInt(n))[0];else if(kind==='prefix-large')want=largeReference(+index,+n);else want=O.residue(E.evaluate(r.ast,+n));assert.equal(value,String(want),line+' expected '+want);checks++}
console.log('PASS fast C++ checks:',checks);fs.writeFileSync(artifact('fast-results.json'),JSON.stringify(results.map(r=>({source:r.source,start:r.start,g:r.g,fast:r.fast}))));

console.log("Compiler artifacts:",testDir);
