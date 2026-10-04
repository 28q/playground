(function(root){
'use strict';
const O=root.SigmaOptimize||(typeof require!=='undefined'?require('./optimize.js'):null);
const units={
mod:{deps:[],code:String.raw`#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;
using u64 = uint64_t;
constexpr int MOD = 998244353;
struct mint {
 int v;
 mint(long long x=0):v(int(x%MOD)){if(v<0)v+=MOD;}
 mint operator+(mint b)const{int x=v+b.v;return mint(x>=MOD?x-MOD:x);}
 mint operator-(mint b)const{return mint(v-b.v);}
 mint operator-()const{return mint(v?MOD-v:0);}
 mint operator*(mint b)const{return mint(1LL*v*b.v);}
 mint& operator+=(mint b){return *this=*this+b;}
 mint& operator-=(mint b){return *this=*this-b;}
 mint& operator*=(mint b){return *this=*this*b;}
 mint pow(u64 n)const{mint a=*this,r=1;for(;n;n>>=1,a*=a)if(n&1)r*=a;return r;}
 mint inv()const{if(!v)throw runtime_error("noninvertible denominator modulo 998244353");return pow(MOD-2);}
 mint operator/(mint b)const{return *this*b.inv();}
};
using Poly=vector<mint>;
void trim(Poly& a){while(a.size()>1&&!a.back().v)a.pop_back();if(a.empty())a.push_back(0);}
mint eval(const Poly& p,mint x){mint r=0;for(int i=int(p.size())-1;i>=0;--i)r=r*x+p[i];return r;}
// This is an allocation limit for prefix/FPS methods, not a 256-term sample limit.
void prefix_limit(int N){if(N<0||N>(1<<21)-1)throw runtime_error("prefix requires 0 <= N < 2^21");}
`},
ntt:{deps:['mod'],code:String.raw`void ntt(Poly& a,bool inverse){
 int n=int(a.size());if(n==0||(n&(n-1))||n>(1<<23))throw runtime_error("NTT length exceeds 2^23");
 for(int i=1,j=0;i<n;i++){int bit=n>>1;for(;j&bit;bit>>=1)j^=bit;j^=bit;if(i<j)swap(a[i],a[j]);}
 for(int len=2;len<=n;len<<=1){mint root=mint(3).pow((MOD-1)/len);if(inverse)root=root.inv();for(int i=0;i<n;i+=len){mint w=1;for(int j=0;j<len/2;j++){mint u=a[i+j],v=a[i+j+len/2]*w;a[i+j]=u+v;a[i+j+len/2]=u-v;w*=root;}}}
 if(inverse){mint inv=mint(n).inv();for(auto& x:a)x*=inv;}
}
Poly mul(Poly a,Poly b){
 if(a.empty()||b.empty())return {};int need=int(a.size()+b.size()-1);
 if(min(a.size(),b.size())<=32){Poly c(need);for(size_t i=0;i<a.size();i++)for(size_t j=0;j<b.size();j++)c[i+j]+=a[i]*b[j];return c;}
 int n=1;while(n<need)n<<=1;if(n>(1<<23))throw runtime_error("convolution exceeds NTT capacity");a.resize(n);b.resize(n);ntt(a,false);ntt(b,false);for(int i=0;i<n;i++)a[i]*=b[i];ntt(a,true);a.resize(need);return a;
}
`},
inv:{deps:['ntt'],code:String.raw`Poly fps_inv(const Poly& f,int n){
 if(n==0)return {};if(f.empty()||!f[0].v)throw runtime_error("FPS inverse needs nonzero constant term");Poly g{f[0].inv()};
 while(int(g.size())<n){int m=min(n,int(g.size())*2);Poly a(f.begin(),f.begin()+min(m,int(f.size())));Poly t=mul(a,g);t.resize(m);for(auto& x:t)x=-x;t[0]+=2;g=mul(g,t);g.resize(m);}return g;
}
`},
log:{deps:['inv'],code:String.raw`Poly fps_log(const Poly& f,int n){
 if(f.empty()||f[0].v!=1)throw runtime_error("FPS log needs constant term 1");if(n==0)return {};
 Poly d(max(0,min(n,int(f.size()))-1));for(int i=1;i<=int(d.size());i++)d[i-1]=f[i]*i;Poly h=mul(d,fps_inv(f,n));h.resize(max(0,n-1));Poly out(n),iv(max(2,n));iv[1]=1;for(int i=2;i<n;i++)iv[i]=-mint(MOD/i)*iv[MOD%i];for(int i=1;i<n;i++)out[i]=h[i-1]*iv[i];return out;
}
`},
exp:{deps:['log'],code:String.raw`Poly fps_exp(const Poly& f,int n){
 if(!f.empty()&&f[0].v)throw runtime_error("FPS exp needs constant term 0");if(n==0)return {};Poly g{1};
 while(int(g.size())<n){int m=min(n,int(g.size())*2);Poly h=fps_log(g,m);for(int i=0;i<m;i++)h[i]=(i<int(f.size())?f[i]:mint(0))-h[i];h[0]+=1;g=mul(g,h);g.resize(m);}return g;
}
`},
factorials:{deps:['mod'],code:String.raw`void factorial_tables(int n,Poly& fac,Poly& ifac){
 if(n<0||n>=MOD)throw runtime_error("factorial table requires n < MOD");fac.assign(n+1,1);ifac.assign(n+1,1);for(int i=1;i<=n;i++)fac[i]=fac[i-1]*i;ifac[n]=fac[n].inv();for(int i=n;i>0;i--)ifac[i-1]=ifac[i]*i;
}
`},
bostan:{deps:['ntt'],code:String.raw`mint bostan_mori(Poly P,Poly Q,u64 N){
 if(Q.empty()||!Q[0].v)throw runtime_error("rational series needs Q(0) != 0");trim(P);trim(Q);
 while(N){Poly minus=Q;for(size_t i=1;i<minus.size();i+=2)minus[i]=-minus[i];Poly s=mul(P,minus),t=mul(Q,minus);P.clear();Q.clear();for(size_t i=N&1;i<s.size();i+=2)P.push_back(s[i]);for(size_t i=0;i<t.size();i+=2)Q.push_back(t[i]);trim(P);trim(Q);N>>=1;}return P[0]/Q[0];
}
`},
batch:{deps:['mod'],code:String.raw`Poly batch_inverse(const Poly& a){
 Poly pref(a.size()+1,1),out(a.size());for(size_t i=0;i<a.size();i++){if(!a[i].v)throw runtime_error("recurrence has a singular denominator beyond the embedded initial terms");pref[i+1]=pref[i]*a[i];}mint x=pref.back().inv();for(size_t i=a.size();i-->0;){out[i]=x*pref[i];x*=a[i];}return out;
}
`},
rec_prefix:{deps:['batch'],code:String.raw`Poly recurrence_prefix(int N,int start,const vector<Poly>& p,const Poly& seeds){
 prefix_limit(N);if(N<start)throw runtime_error("N is before START");int r=int(p.size())-1;if(int(seeds.size())<r)throw runtime_error("insufficient initial terms");Poly a=seeds;a.resize(min(int(a.size()),N-start+1));
 for(int t=start+int(a.size());t<=N;){int end=min(N+1,t+4096);Poly den(end-t);for(int u=t;u<end;u++)den[u-t]=eval(p[r],u-r);Poly inv=batch_inverse(den);for(int u=t;u<end;u++){mint s=0;for(int j=0;j<r;j++)s+=eval(p[j],u-r)*a[u-r+j-start];a.push_back(-s*inv[u-t]);}t=end;}return a;
}
`},
rec_nth:{deps:['batch'],code:String.raw`mint recurrence_nth(u64 N,int start,const vector<Poly>& p,const Poly& seeds){
 if(N<u64(start)||N>=MOD)throw runtime_error("P-recursive evaluation requires START <= N < MOD");int r=int(p.size())-1;if(int(seeds.size())<r)throw runtime_error("insufficient initial terms");if(N-start<seeds.size())return seeds[N-start];
 Poly ring(r);for(int i=0;i<r;i++)ring[(seeds.size()-r+i)%r]=seeds[seeds.size()-r+i];
 for(int t=start+int(seeds.size());t<=int(N);){int end=min(int(N)+1,t+4096);Poly den(end-t);for(int u=t;u<end;u++)den[u-t]=eval(p[r],u-r);Poly inv=batch_inverse(den);for(int u=t;u<end;u++){mint s=0;for(int j=0;j<r;j++)s+=eval(p[j],u-r)*ring[(u-r+j-start)%r];ring[(u-start)%r]=-s*inv[u-t];}t=end;}return ring[(N-start)%r];
}
`},
multipoint:{deps:['inv'],code:String.raw`Poly poly_rem(Poly a,const Poly& b){
 trim(a);if(b.size()==1)return {0};if(a.size()<b.size())return a;int k=int(a.size()-b.size()+1);Poly ar=a,br=b;reverse(ar.begin(),ar.end());reverse(br.begin(),br.end());ar.resize(k);Poly qr=mul(ar,fps_inv(br,k));qr.resize(k);reverse(qr.begin(),qr.end());Poly prod=mul(qr,b);a.resize(b.size()-1);for(size_t i=0;i<a.size();i++)a[i]-=prod[i];trim(a);return a;
}
struct Multipoint {
 int n,size;vector<Poly> tree;
 explicit Multipoint(const Poly& xs):n(int(xs.size())),size(1){while(size<n)size*=2;tree.resize(2*size);for(int i=0;i<size;i++)tree[size+i]=i<n?Poly{-xs[i],1}:Poly{1};for(int i=size-1;i;i--)tree[i]=mul(tree[2*i],tree[2*i+1]);}
 void dfs(int u,const Poly& a,Poly& out)const{if(u>=size){if(u-size<n)out[u-size]=a.empty()?mint(0):a[0];return;}dfs(u*2,poly_rem(a,tree[u*2]),out);dfs(u*2+1,poly_rem(a,tree[u*2+1]),out);}
 Poly evaluate(const Poly& p)const{Poly out(n);dfs(1,poly_rem(p,tree[1]),out);return out;}
};
Poly poly_shift_small(const Poly& p,mint c){Poly out(p.size()),power{1};for(size_t j=0;j<p.size();j++){for(size_t k=0;k<power.size();k++)out[k]+=p[j]*power[k];Poly next(power.size()+1);for(size_t k=0;k<power.size();k++){next[k]+=c*power[k];next[k+1]+=power[k];}power=next;}trim(out);return out;}
`},
bsgs:{deps:['multipoint','rec_nth'],code:String.raw`using PMatrix=vector<vector<Poly>>;
struct Block { PMatrix M;Poly den; };
PMatrix matrix_product(const PMatrix& A,const PMatrix& B){
 int r=int(A.size());PMatrix C(r,vector<Poly>(r,Poly{0}));for(int i=0;i<r;i++)for(int k=0;k<r;k++)for(int j=0;j<r;j++){Poly t=mul(A[i][k],B[k][j]);if(C[i][j].size()<t.size())C[i][j].resize(t.size());for(size_t x=0;x<t.size();x++)C[i][j][x]+=t[x];}for(auto& row:C)for(auto& p:row)trim(p);return C;
}
Block matrix_block(const vector<Poly>& p,int lo,int hi,int origin){
 int r=int(p.size())-1;if(hi-lo==1){Poly den=poly_shift_small(p[r],origin+lo);PMatrix M(r,vector<Poly>(r,Poly{0}));for(int j=0;j<r-1;j++)M[j][j+1]=den;for(int j=0;j<r;j++){M[r-1][j]=poly_shift_small(p[j],origin+lo);for(auto& x:M[r-1][j])x=-x;}return {M,den};}
 int mid=(lo+hi)/2;Block a=matrix_block(p,lo,mid,origin),b=matrix_block(p,mid,hi,origin);return {matrix_product(b.M,a.M),mul(a.den,b.den)};
}
mint recurrence_bsgs(u64 N,int start,const vector<Poly>& p,const Poly& seeds){
 if(N<u64(start)||N>=MOD)throw runtime_error("P-recursive BSGS requires START <= N < MOD");int r=int(p.size())-1;if(r<1||r>4)throw runtime_error("BSGS supports order 1..4");if(N-start<seeds.size())return seeds[N-start];u64 steps=N-(start+seeds.size()-1);if(steps<4096)return recurrence_nth(N,start,p,seeds);
 int degree=0;for(auto& a:p)degree=max(degree,int(a.size())-1);int B=max(1,int(sqrt(double(steps)/max(1,degree)))),blocks=int(steps/B);int origin=start+int(seeds.size())-r;
 Block block=matrix_block(p,0,B,origin);Poly xs(blocks);for(int i=0;i<blocks;i++)xs[i]=1LL*i*B;Multipoint points(xs);Poly denInv=batch_inverse(points.evaluate(block.den));vector<vector<Poly>> entries(r,vector<Poly>(r));for(int i=0;i<r;i++)for(int j=0;j<r;j++)entries[i][j]=points.evaluate(block.M[i][j]);
 Poly state(seeds.end()-r,seeds.end());for(int k=0;k<blocks;k++){Poly next(r);for(int i=0;i<r;i++)for(int j=0;j<r;j++)next[i]+=entries[i][j][k]*state[j];for(auto& x:next)x*=denInv[k];state=next;}
 int tail=int(steps-u64(blocks)*B);Poly den(tail);for(int t=0;t<tail;t++)den[t]=eval(p[r],origin+blocks*B+t);Poly inv=batch_inverse(den);for(int t=0;t<tail;t++){mint next=0;for(int j=0;j<r;j++)next-=eval(p[j],origin+blocks*B+t)*state[j];for(int j=0;j<r-1;j++)state[j]=state[j+1];state[r-1]=next*inv[t];}return state.back();
}
`},
algebraic:{deps:['inv'],code:String.raw`Poly algebraic_series(const vector<Poly>& p,mint initial,int n){
 if(p.size()!=3)throw runtime_error("quadratic algebraic equation required");mint residual=p[0][0]+p[1][0]*initial+p[2][0]*initial*initial;if(residual.v)throw runtime_error("invalid initial branch");if(!(p[1][0]+mint(2)*p[2][0]*initial).v)throw runtime_error("singular algebraic branch at x=0");Poly y{initial};
 while(int(y.size())<n){int m=min(n,int(y.size())*2);Poly yy=mul(y,y),value=mul(p[2],yy),linear=mul(p[1],y);value.resize(m);linear.resize(m);for(int i=0;i<m;i++)value[i]+=linear[i]+(i<int(p[0].size())?p[0][i]:mint(0));Poly derivative=mul(p[2],y);derivative.resize(m);for(int i=0;i<m;i++)derivative[i]=mint(2)*derivative[i]+(i<int(p[1].size())?p[1][i]:mint(0));Poly correction=mul(value,fps_inv(derivative,m));y.resize(m);for(int i=0;i<m;i++)y[i]-=correction[i];}return y;
}
`}
};
const lit=x=>String(O.residue(x)),poly=p=>'{'+p.map(lit).join(',')+'}',matrix=p=>'{'+p.map(poly).join(',')+'}';
function generate(result,target='nth'){
 const plan=result.fast;if(!plan)throw Error('この式に対応する高速計算方法はまだ検出できません。探索上限・項数を増やすか、厳密計算モードを使ってください。');
 let used=new Set(),definitions=[],serial=0;function need(name){if(used.has(name))return;for(let x of units[name].deps)need(x);used.add(name)}need('mod');
 function series(p){let id='series_'+serial++,body='';
  if(p.kind==='polynomial')body='const Poly p='+poly(p.p)+'; Poly a(N+1);for(int n=0;n<=N;n++)a[n]=eval(p,n);return a;';
  else if(p.kind==='rational'){need('inv');body='const Poly P='+poly(p.P)+',Q='+poly(p.Q)+';Poly a=mul(P,fps_inv(Q,N+1));a.resize(N+1);return a;'}
  else if(p.kind==='known-recurrence'){need('rec_prefix');body='return recurrence_prefix(N,0,'+matrix(p.p)+','+poly(p.init)+');'}
  else if(p.kind==='bell'){need('exp');need('factorials');body='Poly fac,ifac;factorial_tables(N,fac,ifac);Poly f=ifac;f[0]=0;Poly a=fps_exp(f,N+1);for(int n=0;n<=N;n++)a[n]*=fac[n];return a;'}
  else if(p.kind==='bernoulli'){need('inv');need('factorials');body='Poly fac,ifac;factorial_tables(N+1,fac,ifac);Poly f(N+1);for(int n=0;n<=N;n++)f[n]=ifac[n+1];Poly a=fps_inv(f,N+1);for(int n=0;n<=N;n++)a[n]*=fac[n];return a;'}
  else if(p.kind==='euler'){need('inv');need('factorials');body='Poly fac,ifac;factorial_tables(N,fac,ifac);for(int n=1;n<=N;n+=2)ifac[n]=0;Poly a=fps_inv(ifac,N+1);for(int n=0;n<=N;n++)a[n]*=fac[n];return a;'}
  else if(p.kind==='partition'){need('exp');body='Poly sigma(N+1),iv(max(2,N+1));iv[1]=1;for(int d=1;d<=N;d++)for(int n=d;n<=N;n+=d)sigma[n]+=d;for(int n=2;n<=N;n++)iv[n]=-mint(MOD/n)*iv[MOD%n];for(int n=1;n<=N;n++)sigma[n]*=iv[n];return fps_exp(sigma,N+1);'}
  else if(p.kind==='algebraic'){need('algebraic');body='return algebraic_series('+matrix(p.polynomials)+','+lit(p.initial)+',N+1);'}
  else if(p.kind==='scale'){let child=series(p.child);body='Poly a='+child+'(N);for(auto& x:a)x*=mint('+lit(p.scale)+');return a;'}
  else if(p.kind==='prefix'){let child=series(p.child);body='Poly a='+child+'(N);for(int n=1;n<=N;n++)a[n]+=a[n-1];return a;'}
  else if(p.kind==='binomial-transform'){let child=series(p.child);need('ntt');need('factorials');body='Poly a='+child+'(N),fac,ifac;factorial_tables(N,fac,ifac);for(int n=0;n<=N;n++)a[n]*=ifac[n];a=mul(a,ifac);a.resize(N+1);for(int n=0;n<=N;n++)a[n]*=fac[n];return a;'}
  else if(['add','subtract','hadamard','convolution','binomial-convolution'].includes(p.kind)){let left=series(p.left),right=series(p.right);body='Poly a='+left+'(N),b='+right+'(N);';if(['add','subtract','hadamard'].includes(p.kind))body+='for(int n=0;n<=N;n++)a[n]'+({add:'+=',subtract:'-=',hadamard:'*='}[p.kind])+'b[n];return a;';else{need('ntt');if(p.kind==='binomial-convolution'){need('factorials');body+='Poly fac,ifac;factorial_tables(N,fac,ifac);for(int n=0;n<=N;n++){a[n]*=ifac[n];b[n]*=ifac[n];}'}body+='a=mul(a,b);a.resize(N+1);';if(p.kind==='binomial-convolution')body+='for(int n=0;n<=N;n++)a[n]*=fac[n];';body+='return a;'}}else throw Error('未対応の生成プラン: '+p.kind);
  definitions.push('Poly '+id+'(int N){ '+body+' }');return id;
 }
 let wrapper='',start=plan.start||0;
 if(target==='nth'){
  if(plan.kind==='polynomial')wrapper='mint coefficient(u64 N){return eval(Poly'+poly(plan.p)+',mint(N%MOD));}';
  else if(plan.kind==='rational'){need('bostan');wrapper='mint coefficient(u64 N){if(N<'+start+'ULL)throw runtime_error("N is before START");return bostan_mori(Poly'+poly(plan.P)+',Poly'+poly(plan.Q)+',N-'+start+');}';}
  else if(plan.kind==='p-recursive'){let bsgs=plan.order<=4;need(bsgs?'bsgs':'rec_nth');wrapper='mint coefficient(u64 N){const vector<Poly> p='+matrix(plan.p)+';const Poly initial='+poly(plan.seeds)+';return '+(bsgs?'recurrence_bsgs':'recurrence_nth')+'(N,'+start+',p,initial);}';}
  else {let s=series(plan);wrapper='mint coefficient(u64 N){if(N>=(1ULL<<21))throw runtime_error("this FPS method requires N < 2^21");return '+s+'(int(N))[N];}';}
 }else{
  if(plan.kind==='p-recursive'){need('rec_prefix');wrapper='Poly coefficients(int N){return recurrence_prefix(N,'+start+','+matrix(plan.p)+','+poly(plan.seeds)+');}';}
  else {let s=series(plan);wrapper='Poly coefficients(int N){prefix_limit(N);if(N<'+start+')throw runtime_error("N is before START");return '+s+'(N-'+start+');}';}
 }
 let comment='// C++17, self-contained, modulo 998244353. No Boost / ACL required.\n// Input: '+result.source.replace(/[\r\n]/g,' ')+'\n// '+(plan.proven?'Certified from identities / known recurrence relations.':'CONJECTURE: fitted to finite samples and checked on 8 held-out terms; not a proof.')+'\n// Method: '+plan.method+'\n// '+(target==='nth'?'coefficient(N) returns a(N).':'coefficients(N)[i] returns a(START+i), START='+start+'.')+'\n// Cost: '+(target==='nth'?plan.nthComplexity:plan.prefixComplexity)+' field operations.\n';
 if(plan.kind==='p-recursive')comment+='// Requires START <= N < MOD. Singular denominators beyond embedded seeds throw.\n';
 if(!['rational','polynomial','p-recursive'].includes(plan.kind))comment+='// Computes the prefix through N. It is not a logarithmic-time N-th coefficient method.\n';
 let main=target==='nth'?'cout<<coefficient(u64(N)).v<<"\\n";':'if(N>=(1LL<<21))throw runtime_error("prefix requires N < 2^21");auto a=coefficients(int(N));for(size_t i=0;i<a.size();i++)cout<<('+start+'+i)<<" "<<a[i].v<<"\\n";';
 return comment+[...used].map(k=>units[k].code).join('\n')+'\n'+definitions.join('\n')+'\n'+wrapper+'\n\n#ifndef SIGMA_NO_MAIN\nint main(){try{long long N;if(!(cin>>N)||N<0)throw runtime_error("N must be nonnegative");'+main+'}catch(const exception& e){cerr<<e.what()<<"\\n";return 1;}}\n#endif\n';
}
function library(names){let used=new Set();function add(k){if(used.has(k))return;for(let dep of units[k].deps)add(dep);used.add(k)}for(let name of names)add(name);return [...used].map(k=>units[k].code).join('\n')}
root.SigmaFastCPP={generate,library};if(typeof module!=='undefined')module.exports=root.SigmaFastCPP;
})(typeof globalThis!=='undefined'?globalThis:this);
