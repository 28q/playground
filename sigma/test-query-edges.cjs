const fs=require('fs'),cp=require('child_process'),path=require('path'),os=require('os'),assert=require('assert/strict'),Q=require('./queries.js'),C=require('./query-cpp.js');
const dir=fs.mkdtempSync(path.join(os.tmpdir(),'sigma-query-edge-'));
function compile(code,name){const file=path.join(dir,name+'.cpp'),exe=path.join(dir,name+'.exe');fs.writeFileSync(file,code);try{cp.execFileSync(process.env.CXX||'g++',['-std=c++17','-O2',file,'-o',exe],{stdio:'pipe',timeout:120000})}catch(e){console.error(String(e.stderr));throw e}return exe;}
const specs=[['mo','sum(binomial(n,k),k,0,m)','binomial-mo'],['fps','stirling2(n,m)','stirling-fps'],['dp','stirling2(n,m)','row-dp'],['floor_test','sum(floor((3*k+2)/m),k,0,n-1)','floor-sum'],['fenwick','sum(array(k),k,n,m)','fenwick'],['seg','range_min(n,m)','segment'],['product','prod(array(k),k,n,m)','product-prefix']];
let code='#include <bits/stdc++.h>\n#define SIGMA_NO_MAIN\n';for(let [name,source,id] of specs)code+='namespace '+name+' {\n'+C.generate(Q.analyze(source),id).replace(/^#include .*$/gm,'')+'\n}\n';
code+=String.raw`int main(){
 auto moAnswers=mo::solve_queries({{10000,20,0},{9999,9999,1},{10000,9998,2},{2,10000,3},{0,0,4}});
 mo::Comb comb;mo::mint ref=0;for(int k=0;k<=20;k++)ref+=comb.choose(10000,k);assert(moAnswers[0].v==ref.v);assert(moAnswers[1].v==mo::mint(2).pow(9999).v);assert(moAnswers[2].v==(mo::mint(2).pow(10000)-comb.choose(10000,9999)-1).v);assert(moAnswers[3].v==4);assert(moAnswers[4].v==1);
 std::vector<fps::Query> fq;std::vector<dp::Query> dq;for(int k=0;k<=100;k++){fq.push_back({1000,k,k});dq.push_back({1000,k,k});}auto a=fps::solve_queries(fq);auto b=dp::solve_queries(dq);for(int i=0;i<=100;i++)assert(a[i].v==b[i].v);
 auto fl=floor_test::solve_queries({{1000000000000LL,1,0},{0,1,1}});using W=unsigned __int128;W N=1000000000000ULL;assert(fl[0].v==int((3*N*(N-1)/2+2*N)%998244353));assert(fl[1].v==0);
 fenwick::Fenwick bit({2,-5,7,3});bit.add(1,10);assert((bit.prefix(4)-bit.prefix(1)).v==15);
 seg::SegmentTree tree({2,-5,7,3});tree.set(1,9);assert(tree.query(0,4)==2);tree.set(2,-8);assert(tree.query(1,4)==-8);
 auto prod=product::solve_queries({{0,0,0},{0,2,1},{2,3,2}},{998244353,7,-2,5});assert(prod[0].v==0&&prod[1].v==0&&prod[2].v==998244343);
 bool caught=false;try{mo::solve_queries({{-1,2,0}});}catch(...){caught=true;}assert(caught);caught=false;try{product::solve_queries({{2,1,0}},{1,2,3});}catch(...){caught=true;}assert(caught);caught=false;try{floor_test::solve_queries({{10,0,0}});}catch(...){caught=true;}assert(caught);caught=false;try{seg::solve_queries({{0,2,0}},{1,2});}catch(...){caught=true;}assert(caught);
 std::cout<<"PASS large Mo / NTT-vs-DP / 128-bit floor sum / update APIs / zero product / domain guards\n";
}
`;
console.log(cp.execFileSync(compile(code,'edge'),[],{encoding:'utf8'}).trim());
const grid=compile(C.generate(Q.analyze('binomial(n+m,n)'),'direct','grid'),'grid');let lines=cp.execFileSync(grid,[],{input:'2 3\n',encoding:'utf8'}).trim().split(/\r?\n/);assert.equal(lines.length,12);assert.equal(lines[11],'2 3 10');
const cli=compile(C.generate(Q.analyze('distinct(n,m)'),'distinct-mo'),'cli');assert.equal(cp.execFileSync(cli,[],{input:'5\n2 1 2 3 1\n3\n0 4\n1 1\n1 3\n',encoding:'utf8'}).trim().replace(/\r/g,''),'3\n1\n3');assert.notEqual(cp.spawnSync(cli,[],{input:'1\n5\n1\n0 2\n'}).status,0);
const fixed=Q.run('sum(binomial(n,k)^m,k,0,n)',{parameters:'m=2'});assert.equal(fixed.fast.proven,true);assert.equal(String(fixed.seq[5]),'252');
console.log('PASS grid main / array main / fixed parameters');console.log('Artifacts:',dir);


