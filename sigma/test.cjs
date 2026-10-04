const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),os=require('node:os'),cp=require('node:child_process');
const E=require('./engine.js'),C=require('./cpp.js');let checks=0;
function eq(a,b){assert.equal(String(a),String(b));checks++}
function value(s,n=0){return E.evaluate(E.parse(s),n).toString()}
const cases=[['binomial(10,3)','120'],['binomial(-3,2)','6'],['factorial(10)','3628800'],['doublefactorial(-1)','1'],['doublefactorial(8)','384'],['rising(1/2,3)','15/8'],['falling(5,3)','60'],['stirling1(5,2)','50'],['stirling1s(5,2)','-50'],['stirling2(5,2)','15'],['bell(8)','4140'],['catalan(7)','429'],['fibonacci(15)','610'],['lucas(10)','123'],['derangement(6)','265'],['partition(10)','42'],['motzkin(6)','51'],['schroeder(5)','394'],['eulerian(5,2)','66'],['lah(5,2)','240'],['bernoulli(1)','-1/2'],['bernoulli(12)','-691/2730'],['euler(6)','-61'],['harmonic(5)','137/60'],['harmonic2(3,2)','49/36'],['totient(36)','12'],['mobius(30)','-1'],['mobius(12)','0'],['divisors(36)','9'],['divisorsum(12)','28'],['gcd(-12,18)','6'],['lcm(-12,18)','36'],['abs(-2/3)','2/3'],['floor(-3/2)','-2'],['ceil(-3/2)','-1'],['min(2/3,3/4)','2/3'],['max(2/3,3/4)','3/4'],['0.1+0.2','3/10'],['2^3^2','512'],['-2^2','-4'],['3!+0!','7'],['sum(k,k,3,2)','0'],['prod(k,k,3,2)','1'],['sum(sum(j,j,0,k),k,0,4)','20'],['sum(k+sum(k,k,1,2),k,1,2)','9'],['Sum[Binomial[n,k]^2,{k,0,n}]','252',5]];
for(let [s,expected,n] of cases)eq(value(s,n),expected);
for(let s of ['unknown(n)','n+k','sum(k,k,0)','sum(k,n+1,0,n)','factorial(-1)','1/0','n;alert(1)','2n','sum(k,k,0,1/2)','factorial(300)','bell(n,n)']){assert.throws(()=>value(s));checks++}
for(const s of ['sum(sum(k,k,1,n),n,-2,2)','sum(binomial(n,n),n,0,n)']){const ast=E.parse(s),simple=E.simplify(ast);for(let n=0;n<8;n++)eq(E.evaluate(ast,n),E.evaluate(simple,n));}
const expressions=['sum(binomial(n,k),k,0,n)','sum(binomial(n,k)^2,k,0,n)','sum(stirling2(n,k),k,0,n)','sum(stirling1(n,k),k,0,n)','prod(k,k,1,n)','sum(k,k,1,n)','sum(k^3,k,1,n)','fibonacci(n)','catalan(n)','harmonic(n)','sum(binomial(n,k)^2*binomial(n+k,k)^2,k,0,n)','prod((k+1)/k,k,1,n)'];
const results=expressions.map(s=>E.analyze(s));
for(let a of results){for(let n=0;n<15;n++)eq(E.evaluate(a.ast,n),E.evaluate(a.simple,n));if(a.g){for(let n=40;n<50;n++){let total=E.q(0);a.g.p.forEach((p,j)=>{total=total.add(E.polyEval(p,n).mul(E.evaluate(a.ast,n+j)))});eq(total,0)}}}
eq(results[0].gf.formula,'1 / (1 − 2x)');eq(results[1].gf.formula,'1 / √(1 − 4x)');eq(results[2].gf.formula,'exp(exp(x) − 1)');eq(results[2].g,null);eq(results[10].g.order,2);eq(results[10].g.degree,3);
const shifted=E.analyze('fibonacci(n)',{start:5});eq(shifted.seq[0],5);eq(shifted.gf,null);
const zero=E.analyze('0');eq(zero.g.order,1);
const convolution=E.analyze('sum(catalan(k)*catalan(n-k),k,0,n)');
eq(convolution.gf.algebraic,true);eq(convolution.gf.proven,false);eq(convolution.gf.formula,'(1) + (−1 + 2x) G(x) + (x²) G(x)² = 0');
// Verify the inferred algebraic relation beyond all coefficients used by discovery.
const conv=Array.from({length:52},(_,n)=>E.evaluate(convolution.ast,n));
for(let n=40;n<52;n++){let t=conv[n].neg().add(conv[n-1].mul(2));for(let k=0;k<=n-2;k++)t=t.add(conv[k].mul(conv[n-2-k]));eq(t,0)}
console.log('PASS: '+checks+' mathematical/parser/holdout checks');
// Optional integration suite: BOOST_INCLUDE=/path/to/boost_root node test.cjs --cpp
if(process.argv.includes('--cpp')){
 const temp=fs.mkdtempSync(path.join(os.tmpdir(),'sigma-lab-test-')),include=process.env.BOOST_INCLUDE;
 function compileRun(code,input){let cpp=path.join(temp,'test.cpp'),exe=path.join(temp,process.platform==='win32'?'test.exe':'test');fs.writeFileSync(cpp,code);cp.execFileSync(process.env.CXX||'g++',['-std=c++17','-O0',...(include?['-I',include]:[]),cpp,'-o',exe],{stdio:'pipe'});return cp.execFileSync(exe,[],{input:String(input)+'\n',encoding:'utf8'}).trim()}
 // One compilation exercises every function and the parser-to-code path.
 let a=E.analyze('n'),base=C.generate(a),cut=base.indexOf('int main()');let harness=base.slice(0,cut)+'int main(){\n';
 // Generate direct() for each test, with distinct names, to compare exact text.
 for(let i=0;i<cases.length;i++){let [s,expected,n=0]=cases[i],ast=E.parse(s),snippet=C.generate({...a,source:s,simple:ast});let direct=snippet.slice(snippet.indexOf('Q direct('),snippet.indexOf('int main()'));harness=harness.replace('int main(){\n',direct.replace('Q direct(','Q f'+i+'(')+'int main(){\n');harness+='cout<<f'+i+'(Q('+n+')).str()<<"\\n";\n'}harness+='}\n';let lines=compileRun(harness,0).split(/\r?\n/);cases.forEach((c,i)=>eq(lines[i],c[1]));console.log('PASS: generated C++ function parity, '+cases.length+' cases');
 for(let i of [1,2,6,10]){let a=results[i];for(let mode of a.g?['direct','recurrence']:['direct']){let lines=compileRun(C.generate(a,mode),45).split(/\r?\n/);for(let n=0;n<=45;n++)eq(lines[n],n+' '+E.evaluate(a.ast,n));console.log('PASS: C++ '+mode+' '+a.source)}}
 let shiftedLines=compileRun(C.generate(shifted,'recurrence'),45).split(/\r?\n/);for(let i=0;i<shiftedLines.length;i++)eq(shiftedLines[i],(i+5)+' '+E.evaluate(shifted.ast,i+5));console.log('PASS: C++ nonzero start; total assertions '+checks);console.log('Temporary compiler artifacts: '+temp);
}
