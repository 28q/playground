(function(root){
'use strict';
const runtime=String.raw`#include <boost/multiprecision/cpp_int.hpp>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include <algorithm>
using namespace std;
using Z = boost::multiprecision::cpp_int;
Z az(Z a){return a<0?-a:a;}
Z gz(Z a,Z b){a=az(a);b=az(b);while(b!=0){Z t=a%b;a=b;b=t;}return a;}
struct Q {
 Z a,b;
 Q(long long n=0):a(n),b(1){}
 Q(Z n,Z d):a(n),b(d){if(b==0)throw runtime_error("division by zero");if(b<0){a=-a;b=-b;}Z g=gz(a,b);a/=g;b/=g;}
 Q operator+(Q y)const{return Q(a*y.b+y.a*b,b*y.b);}
 Q operator-()const{return Q(-a,b);}
 Q operator-(Q y)const{return *this+(-y);}
 Q operator*(Q y)const{return Q(a*y.a,b*y.b);}
 Q operator/(Q y)const{return Q(a*y.b,b*y.a);}
 bool iszero()const{return a==0;}
 string str()const{return a.convert_to<string>()+(b==1?"":"/"+b.convert_to<string>());}
};
long long ix(Q x,long long lo=-10000,long long hi=10000){if(x.b!=1||x.a<lo||x.a>hi)throw runtime_error("integer argument out of range");return x.a.convert_to<long long>();}
Q pw(Q a,Q e){long long k=ix(e,-1000,1000);if(k<0)return Q(1)/pw(a,Q(-k));Q r(1);for(;k;k/=2,a=a*a)if(k%2)r=r*a;return r;}
Q fac(int n){Q r(1);for(int i=2;i<=n;i++)r=r*Q(i);return r;}
Q choose(int n,int k){if(k<0||(n>=0&&k>n))return Q(0);Q r(1);for(int i=1;i<=k;i++)r=r*Q(n-i+1)/Q(i);return r;}
Q F(string name,vector<Q> a){
 static map<string,Q> memo;
 string key=name;for(Q x:a)key+=":"+x.str();auto it=memo.find(key);if(it!=memo.end())return it->second;
 Q x=a[0],y=a.size()>1?a[1]:Q(0),r;
 if(name=="abs")r=Q(az(x.a),x.b);
 else if(name=="floor")r=Q(Z(x.a/x.b-(x.a<0&&x.a%x.b!=0?1:0)),Z(1));
 else if(name=="ceil")r=Q(Z(x.a/x.b+(x.a>0&&x.a%x.b!=0?1:0)),Z(1));
 else if(name=="min"||name=="max")r=((x-y).a<0)==(name=="min")?x:y;
 else if(name=="gcd"||name=="lcm"){if(x.b!=1||y.b!=1)throw runtime_error("integer required");Z g=gz(x.a,y.a);r=Q(name=="gcd"?g:g==0?Z(0):Z(az(x.a*y.a)/g),Z(1));}
 else if(name=="rising"||name=="falling"){int k=ix(y,0,256);r=Q(1);for(int i=0;i<k;i++)r=r*(x+Q(name=="rising"?i:-i));}
 else if(name=="binomial")r=choose(ix(x,-512,512),ix(y,-512,512));
 else {
 int n=ix(x,name=="doublefactorial"?-1:0,256),k=ix(y,-256,256);
 auto f=[&](string id,int i){return F(id,{Q(i)});};
 if(name=="factorial")r=fac(n);
 else if(name=="doublefactorial"){r=Q(1);for(int i=n;i>0;i-=2)r=r*Q(i);}
 else if(name=="stirling1"||name=="stirling1s"||name=="stirling2"||name=="eulerian"){
 vector<Q> row(n+2);row[0]=Q(1);for(int i=1;i<=n;i++){vector<Q> next(n+2);for(int j=0;j<=i;j++){Q left=j?row[j-1]:Q(0);next[j]=name=="eulerian"?row[j]*Q(j+1)+left*Q(i-j):left+row[j]*Q(name=="stirling2"?j:name=="stirling1s"?1-i:i-1);}row=next;}r=k<0||k>n?Q(0):row[k];
 }else if(name=="bell"){for(int j=0;j<=n;j++)r=r+F("stirling2",{Q(n),Q(j)});}
 else if(name=="catalan")r=choose(2*n,n)/Q(n+1);
 else if(name=="fibonacci"||name=="lucas"){Q u(name=="fibonacci"?0:2),v(1);for(int i=0;i<n;i++){Q t=u+v;u=v;v=t;}r=u;}
 else if(name=="derangement"){Q u(1),v(0);for(int i=2;i<=n;i++){Q t=(u+v)*Q(i-1);u=v;v=t;}r=n==0?Q(1):v;}
 else if(name=="partition"){vector<Q>d(n+1);d[0]=Q(1);for(int i=1;i<=n;i++)for(int j=i;j<=n;j++)d[j]=d[j]+d[j-i];r=d[n];}
 else if(name=="motzkin"){Q u(1),v(1);for(int i=2;i<=n;i++){Q t=(v*Q(2*i+1)+u*Q(3*i-3))/Q(i+2);u=v;v=t;}r=v;}
 else if(name=="schroeder"){r=n==0?Q(1):f(name,n-1);for(int j=0;j<n;j++)r=r+f(name,j)*f(name,n-1-j);}
 else if(name=="lah")r=n==0?Q(k==0):k<=0||k>n?Q(0):choose(n-1,k-1)*fac(n)/fac(k);
 else if(name=="bernoulli"){if(n==0)r=Q(1);else {for(int j=0;j<n;j++)r=r-choose(n+1,j)*f(name,j);r=r/Q(n+1);}}
 else if(name=="euler"){if(n==0)r=Q(1);else if(n%2==0)for(int j=0;j<n;j+=2)r=r-choose(n,j)*f(name,j);}
 else if(name=="harmonic"||name=="harmonic2"){for(int j=1;j<=n;j++)r=r+pw(Q(j),Q(name=="harmonic"?-1:-k));}
 else if(name=="totient"||name=="mobius"||name=="divisors"||name=="divisorsum"){
 if(n<1)throw runtime_error("n must be positive");if(name=="totient"){for(int i=1;i<=n;i++)if(gz(Z(i),Z(n))==1)r=r+Q(1);}
 else if(name=="mobius"){r=Q(1);int v=n;for(int p=2;p<=v;p++)if(v%p==0){v/=p;r=-r;if(v%p==0){r=Q(0);break;}}}
 else for(int i=1;i<=n;i++)if(n%i==0)r=r+Q(name=="divisors"?1:i);
 }else throw runtime_error("unknown function");
 }
 memo[key]=r;return r;
}
template<class T> Q range(bool product,Q lo,Q hi,T f){int l=ix(lo,-512,512),h=ix(hi,-512,512);if(h-l>512)throw runtime_error("range too large");Q r(product?1:0);for(int i=l;i<=h;i++){Q v=f(Q(i));r=product?r*v:r+v;}return r;}
`;
function literal(q){return 'Q(Z("'+q.a+'"),Z("'+q.b+'"))'}
function emit(a,env={n:'n'},depth=0){if(a.type==='num')return literal(a.value);if(a.type==='var')return env[a.name];if(a.type==='neg')return '(-'+emit(a.x,env,depth)+')';if(a.type==='bin'){let l=emit(a.left,env,depth),r=emit(a.right,env,depth);return a.op==='^'?'pw('+l+','+r+')':'('+l+a.op+r+')'}if(a.name==='sum'||a.name==='prod'){let [body,v,lo,hi]=a.args,b='bound'+depth;return 'range('+String(a.name==='prod')+','+emit(lo,env,depth)+','+emit(hi,env,depth)+',[&](Q '+b+'){ return '+emit(body,{...env,[v.name]:b},depth+1)+'; })'}return 'F("'+a.name+'", {'+a.args.map(x=>emit(x,env,depth)).join(',')+'})'}
function generateRecurrence(result){
 if(!result.g)throw Error('漸化式がありません');let g=result.g,start=result.start;
 let source='// C++17, exact rational recurrence; requires Boost headers.\n// '+(g.proven?'Certified from known recurrences.':'CONJECTURE, not a proof.')+'\n// Initial terms are embedded; no function-name dispatcher.\n'+runtime.slice(0,runtime.indexOf('long long ix('));
 source+='\nQ coefficient(int N){if(N<'+start+'||N>10000)throw runtime_error("N out of range");static const vector<Q> initial={'+result.seq.map(literal).join(',')+'};static const vector<vector<Q>> p={'+g.p.map(p=>'{'+p.map(literal).join(',')+'}').join(',')+'};if(N-'+start+'<(int)initial.size())return initial[N-'+start+'];static vector<Q> a=initial;for(int t='+start+'+(int)a.size();t<=N;t++){vector<Q> c;for(auto& poly:p){Q v;for(int j=(int)poly.size()-1;j>=0;j--)v=v*Q(t-'+g.order+')+poly[j];c.push_back(v);}if(c.back().iszero())throw runtime_error("unresolved recurrence singularity");Q s;for(int j=0;j<'+g.order+';j++)s=s+c[j]*a[t-'+g.order+'+j-'+start+'];a.push_back(-s/c.back());}return a[N-'+start+'];}\n';
 source+='int main(){try{int N;if(!(cin>>N)||N<'+start+'||N>10000)throw runtime_error("N out of range");for(int n='+start+';n<=N;n++)cout<<n<<" "<<coefficient(n).str()<<"\\n";}catch(const exception& e){cerr<<e.what()<<"\\n";return 1;}}\n';return source;
}
function generate(result,mode='direct'){
 if(mode==='recurrence')return generateRecurrence(result);
 let recurrence=mode==='recurrence'&&result.g,g=result.g;
 let source='// Sigma Lab | C++17 | exact rational arithmetic\n// Requires Boost headers. Build: g++ -std=c++17 -O2 main.cpp -o main\n// stdin: N (last index); stdout: n a(n), starting at '+result.start+'\n// Original: '+result.source.replace(/[\r\n]/g,' ')+'\n// Domain: function arguments <= 256; sum/product bounds [-512,512].\n';
 source+=recurrence?'// WARNING: empirically guessed recurrence, NOT a proof.\n// Validated on '+result.count+' exact terms; '+g.holdout+' held out from fitting.\n':'// Uses the input definition with only proven identity rewrites.\n';
 source+=runtime+'\nQ direct(Q n){ return '+emit(result.simple)+'; }\n';
 if(recurrence){source+='Q coeff(int j,int n){ static const vector<vector<Q>> p={'+g.p.map(p=>'{'+p.map(literal).join(',')+'}').join(',')+'}; Q r; for(int k='+g.degree+';k>=0;k--)r=r*Q(n)+p[j][k];return r;}\n';}
 source+='int main(){try{int N;if(!(cin>>N)||N<'+result.start+'||N>256)throw runtime_error("N out of range ['+result.start+',256]");\n';
 if(recurrence){source+='vector<Q> a; for(int t='+result.start+';t<=N;t++){ Q value; int n=t-'+g.order+'; if(t<'+(result.start+g.order)+') value=direct(Q(t)); else { Q den=coeff('+g.order+',n); if(den.iszero())value=direct(Q(t)); else {Q s;for(int j=0;j<'+g.order+';j++)s=s+coeff(j,n)*a[n+j-'+result.start+'];value=-s/den;}} a.push_back(value);cout<<t<<" "<<value.str()<<"\\n";}\n';}
 else source+='for(int n='+result.start+';n<=N;n++)cout<<n<<" "<<direct(Q(n)).str()<<"\\n";\n';
 return source+'}catch(const exception& e){cerr<<e.what()<<"\\n";return 1;} }\n';
}
root.SigmaCPP={generate};if(typeof module!=='undefined')module.exports=root.SigmaCPP;
})(typeof globalThis!=='undefined'?globalThis:this);
