(function(root){
'use strict';
const E=root.SigmaEngine||(typeof require!=='undefined'?require('./engine.js'):null),q=E.q;
const abs=x=>x<0n?-x:x,gcd=(a,b)=>{a=abs(a);b=abs(b);while(b)[a,b]=[b,a%b];return a},mod=(a,m)=>(a%m+m)%m;
function power(a,n,p){let r=1;for(;n;n=Math.floor(n/2),a=a*a%p)if(n%2)r=r*a%p;return r}
function primes(count){let out=[];for(let p=1000003;out.length<count;p+=2){let ok=true;for(let d=3;d*d<=p;d+=2)if(p%d===0){ok=false;break}if(ok)out.push(p)}return out}
const PRIMES=primes(48);
function sqrt(n){if(n<2n)return n;let x=1n<<BigInt(Math.ceil(n.toString(2).length/2));for(;;){let y=(x+n/x)/2n;if(y>=x)return x;x=y}}
function reconstruct(a,m){const bound=sqrt(m/2n);let r0=m,r1=mod(a,m),s0=0n,s1=1n;while(abs(r1)>bound){if(!r1)return null;let c=r0/r1;[r0,r1]=[r1,r0-c*r1];[s0,s1]=[s1,s0-c*s1]}if(!s1||abs(s1)>bound||gcd(r1,s1)!==1n||mod(r1-a*s1,m)!==0n)return null;return new E.Q(r1,s1)}
function nullMod(rows,cols,p,deadline){let pivots=[],r=0;for(let c=0;c<cols&&r<rows.length;c++){if(Date.now()>deadline)throw Error('漸化式探索の時間上限です。次数・階数・候補数を下げてください');let k=r;while(k<rows.length&&!rows[k][c])k++;if(k===rows.length)continue;[rows[r],rows[k]]=[rows[k],rows[r]];let row=rows[r],inv=power(row[c],p-2,p);for(let j=c;j<cols;j++)row[j]=row[j]*inv%p;for(let i=r+1;i<rows.length;i++){let f=rows[i][c];if(!f)continue;rows[i][c]=0;for(let j=c+1;j<cols;j++){let x=(rows[i][j]-f*row[j])%p;rows[i][j]=x<0?x+p:x;}}pivots.push(c);r++;}
 if(pivots.length!==cols-1)return null;let free=0;while(pivots.includes(free))free++;const x=new Float64Array(cols);x[free]=1;for(let i=pivots.length-1;i>=0;i--){let c=pivots[i],s=0;for(let j=c+1;j<cols;j++)s=(s+rows[i][j]*x[j])%p;x[c]=s?p-s:0;}return {x,profile:pivots.join(',')};
}
function verify(g,seq,start){for(let i=0;i<seq.length-g.order;i++){let s=q(0);for(let j=0;j<=g.order;j++)s=s.add(E.polyEval(g.p[j],start+i).mul(seq[i+j]));if(!s.zero())return false}return true}
function fit(seq,start,order,degree,options={}){
 const cols=(order+1)*(degree+1),holdout=options.holdout||8,train=seq.length-holdout;if(cols>512||order<1||degree<0||train-order<cols+2)return null;
 function candidate(ints,used){let common=ints.reduce(gcd,0n);if(!common)return null;ints=ints.map(x=>x/common);if(ints.findLast(x=>x!==0n)<0n)ints=ints.map(x=>-x);const p=Array.from({length:order+1},(_,j)=>{let row=ints.slice(j*(degree+1),(j+1)*(degree+1)).map(q);while(row.length>1&&row.at(-1).zero())row.pop();return row});if(p[order].every(x=>x.zero()))return null;const g={order,degree:Math.max(...p.map(a=>a.length-1)),p,train,holdout,method:'複数素数の線形方程式 → CRT → 整数・有理数復元 → 全項の厳密検算',primes:used.slice(),unknowns:cols,proven:false};return verify(g,seq,start)?g:null;}
 const deadline=options.deadline||Date.now()+90000;let residues=null,modulus=1n,profile=null,used=[],bad=0;
 for(const prime of PRIMES){if(Date.now()>deadline)throw Error('漸化式探索の時間上限です。指定範囲を狭めてください');let P=BigInt(prime),s=[];for(let v of seq.slice(0,train)){let b=Number(mod(v.b,P));if(!b){s=null;break}s.push(Number(mod(v.a,P))*power(b,prime-2,prime)%prime)}if(!s)continue;
  let rows=[];for(let i=0;i<Math.min(train-order,cols+3);i++){let row=new Float64Array(cols),powers=[1];for(let k=1;k<=degree;k++)powers.push(powers[k-1]*(start+i)%prime);for(let j=0;j<=order;j++)for(let k=0;k<=degree;k++)row[j*(degree+1)+k]=s[i+j]*powers[k]%prime;rows.push(row)}
  const ns=nullMod(rows,cols,prime,deadline);if(!ns){if(++bad>=2)return null;continue}if(profile&&profile!==ns.profile)continue;profile=ns.profile;
  if(!residues)residues=Array.from(ns.x,BigInt);else{let inv=BigInt(power(Number(modulus%P),prime-2,prime));residues=residues.map((a,i)=>a+mod((BigInt(ns.x[i])-a)*inv,P)*modulus)}modulus*=P;used.push(prime);if(used.length<3)continue;
  const centered=candidate(residues.map(x=>x>modulus/2n?x-modulus:x),used);if(centered){centered.reconstruction='整数への対称復元';return centered;}
  const coeff=residues.map(a=>reconstruct(a,modulus));if(coeff.some(x=>!x))continue;let den=1n;for(const x of coeff)den=den/gcd(den,x.b)*x.b;const recovered=candidate(coeff.map(x=>x.a*(den/x.b)),used);if(recovered){recovered.reconstruction='有理数復元';return recovered;}
 }
 return null;
}
root.SigmaModular={fit,verify,reconstruct,primes:PRIMES};if(typeof module!=='undefined')module.exports=root.SigmaModular;
})(typeof globalThis!=='undefined'?globalThis:this);
