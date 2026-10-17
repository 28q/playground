const E=require('./engine.js'),O=require('./optimize.js'),C=require('./fast-cpp.js'),assert=require('assert/strict');
let checks=0;
// Include empty ranges, negative endpoints, shadowed variables and varying bounds.
for(const source of ['sum(sum(fibonacci(k+1),n,0,k+1)-sum(fibonacci(k),n,0,k),k,0,n)','sum(k^5,k,0,n)','sum(k^2,k,-n,n)','sum(k,k,n,0)','sum(k,k,1,0)','sum(k,k,3,n)','sum(sum(j^2,j,0,k),k,0,n)','sum(sum(n,n,0,k),k,0,n)','sum(sum(j,j,k,n),k,0,n)','sum(binomial(k,3),k,0,n)','sum(binomial(n,k)*2^k*3^(n-k),k,0,n)','sum(binomial(n,k)*binomial(n,n-k),k,0,n)','sum(3^k,k,0,n)','sum(factorial(k+1)-factorial(k),k,0,n)']){
 const ast=E.parse(source),simple=O.simplify(ast);
 for(let n=0;n<18;n++){assert.equal(String(E.evaluate(simple,n)),String(E.evaluate(ast,n)),source+' n='+n);checks++;}
}
const certified=E.analyze('motzkin(n)*sum(binomial(n,k)^2,k,0,n)');
assert(certified.g.proven);assert.equal(certified.g.order,2);assert.equal(certified.g.degree,2);checks+=3;
const guess=E.analyze('floor(n/3)');assert.equal(guess.g.proven,false);checks++;
const fps=E.analyze('bell(n)');assert.equal(fps.fast.kind,'bell');assert.equal(fps.g,null);checks+=2;
const algebraic=E.analyze('binomial(2*n,n)/(n+1)',{maxOrder:1,maxDegree:0});assert.equal(algebraic.fast.kind,'algebraic');assert.equal(algebraic.fast.proven,false);checks+=2;
assert.throws(()=>O.residue(E.q('1/998244353')));checks++;
for(const [r,name] of [[fps,'bell'],[certified,'P-recursive']]){const code=C.generate(r);assert(!code.includes('name=='));assert(!code.includes('Q F('));assert(code.includes('coefficient(u64 N)'));checks+=3;}
console.log('PASS optimizer checks:',checks);
