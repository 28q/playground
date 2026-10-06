(function(root){
'use strict';
const F=root.SigmaFastCPP||(typeof require!=='undefined'?require('./fast-cpp.js'):null),C=root.SigmaCPP||(typeof require!=='undefined'?require('./cpp.js'):null),QC=root.SigmaQueryCPP||(typeof require!=='undefined'?require('./query-cpp.js'):null);
// Only whitespace changes: comments, literals, operators and preprocessor lines
// are single tokens. The same formatter is used for display, copy and download.
function tokens(source){return source.match(/^[ \t]*#[^\r\n]*|\/\/[^\r\n]*|\/\*[\s\S]*?\*\/|"(?:\\[\s\S]|[^"\\])*"|'(?:\\[\s\S]|[^'\\])*'|[A-Za-z_][A-Za-z_0-9]*|\d+(?:\.\d*)?(?:[eE][+-]?\d+)?[A-Za-z_0-9]*|::|->\*|->|\.\*|<<=|>>=|<<|>>|\+\+|--|&&|\|\||[+\-*\/%&|^!=<>]=|[^\s]/gm)||[];}
function format(source){
 const ts=tokens(source),out=[];let line='',indent=0,paren=0,braces=[],angles=0;
 function flush(){if(line.trim())out.push('    '.repeat(Math.max(0,indent))+line.trim());line='';}
 function add(t){const prev=line.at(-1)||'',word=line.match(/[A-Za-z_0-9]+$/)?.[0];const tightParen=t==='('&&!['if','for','while','switch','catch','sizeof'].includes(word);if(line&&![';',',',')',']','[','::','.','->'].includes(t)&&!tightParen&&!['(','[','.'].includes(prev)&&!(prev==='<'&&angles>0)&&!line.endsWith('::')&&!line.endsWith('->'))line+=' ';line+=t;}
 for(let i=0;i<ts.length;i++){let t=ts[i],next=ts[i+1];
  if(t.trimStart().startsWith('#')){flush();out.push(t.trim());continue;}
  if(t==='<'&&/^(vector|array|map|pair|tuple|set|function|static_cast|dynamic_cast|reinterpret_cast|const_cast)$/.test(ts[i-1]||'')){line+='<';angles++;continue;}
  if((t==='>'||t==='>>')&&angles>0){line=line.trimEnd()+t;angles=Math.max(0,angles-t.length);continue;}
  if(t.startsWith('//')||t.startsWith('/*')){flush();add(t);flush();continue;}
  if(t==='{'){add(t);flush();braces.push(paren);indent++;continue;}
  if(t==='}'){flush();indent--;braces.pop();add(t);if(![';',',',')',']','(','else','catch'].includes(next))flush();continue;}
  if(t==='('){add(t);paren++;continue;}
  if(t===')'){paren--;add(t);continue;}
  add(t);
  if(t===';'&&paren<=(braces.at(-1)||0))flush();
  else if(t===','&&line.length>100)flush();
 }
 flush();return out.join('\n')+'\n';
}
function parts(raw){let cut=raw.indexOf('#ifndef SIGMA_NO_MAIN');if(cut<0)cut=raw.lastIndexOf('int main()');if(cut>=0)raw=raw.slice(0,cut);const includes=[...new Set(raw.match(/^#include[^\r\n]*/gm)||[])];return {includes,body:raw.replace(/^#include[^\r\n]*\r?\n/gm,'').replace(/^\/\/[^\r\n]*\r?\n/,'')};}
function assemble(raw,api,types,bridge,body,main,withMain,forward=''){
 const p=parts(raw);p.includes.push('#include <array>','#include <string>');
 let source=api.comment+'\n'+[...new Set(p.includes)].join('\n')+'\n\n'+types+'\n'+forward+'\n';
 if(body){source+='namespace sigma_detail { '+api.returnType+' solve_impl('+api.args+'); }\n\n';source+='// Public API: '+api.description+'\n'+api.signature+' { '+body+' }\n\n';}
 else source+='// Public API: '+api.description+'\n'+api.signature+';\n\n';
 source+='// Algorithm implementation and required helpers.\nnamespace sigma_detail {\n'+p.body+'\n'+bridge+'\n}\n';
 if(!body)source+=api.definition+'\n';
 if(withMain)source+='\n#ifndef SIGMA_NO_MAIN\nint main() { try { '+main+' } catch (const std::exception& error) { std::cerr << error.what() << "\\n"; return 1; } }\n#endif\n';
 return {...api,code:format(source),withMain};
}
function sequence(result,mode='fast-nth',withMain=false){
 const fast=mode.startsWith('fast-'),prefix=mode==='fast-prefix',start=fast?(result.fast?.start||0):result.start,ret=prefix?'std::vector<int>':fast?'int':'Rational',args=prefix||!fast?'int n':'long long n',signature=ret+' solve('+args+')';
 const proven=fast?result.fast?.proven:mode==='direct'||result.g?.proven;
 const description=prefix?'Returns a(START)..a(n), inclusive. START = '+start+'.':fast?'Returns a(n) modulo 998244353.':'Returns the exact rational a(n). Use .str() for output.';
 const comment='// C++17 | '+(fast?'modulo 998244353':'exact rational arithmetic; Boost headers required')+'\n// '+(proven?'PROVEN: identities / definitions / symbolic verification.':'GUESSED: NOT MATHEMATICALLY PROVEN. Finite checks do not establish an identity.')+'\n// Expression: '+result.source.replace(/[\r\n]/g,' ')+'\n// '+description;
 const api={signature,returnType:ret,args,description,comment,example:prefix?'auto values = solve(100); // values[i] = a('+start+' + i)':fast?'int answer = solve(100);':'Rational answer = solve(100);\nstd::cout << answer.str();'};
 let raw=fast?F.generate(result,prefix?'prefix':'nth'):C.generate(result,mode),bridge='',body='',types='',forward='',main='';
 if(fast){bridge=ret+' solve_impl('+args+') { if(n < '+start+') throw std::runtime_error("n is before START");';if(prefix)bridge+='auto values = coefficients(n); std::vector<int> answer; answer.reserve(values.size()); for(auto value : values) answer.push_back(value.v); return answer; }';else bridge+='return coefficient(static_cast<u64>(n)).v; }';body='return sigma_detail::solve_impl(n);';}
 else{forward='namespace sigma_detail { struct Q; }\nusing Rational = sigma_detail::Q;';api.definition=signature+' { if(n < '+start+' || n > '+(mode==='recurrence'?10000:256)+') throw std::runtime_error("n outside supported range"); return sigma_detail::'+(mode==='recurrence'?'coefficient(n)':'direct(sigma_detail::Q(n))')+'; }';}
 main=args+'; if(!(std::cin >> n)) throw std::runtime_error("expected n");';
 main+=prefix?'auto values = solve(n); for(size_t i=0;i<values.size();i++) std::cout << ('+start+'+i) << " " << values[i] << "\\n";':fast?'std::cout << solve(n) << "\\n";':'std::cout << solve(n).str() << "\\n";';
 return assemble(raw,api,types,bridge,body,main,withMain,forward);
}
const reserved=new Set(('alignas alignof and and_eq asm atomic_cancel atomic_commit atomic_noexcept auto bitand bitor bool break case catch char char16_t char32_t class compl concept const constexpr const_cast continue co_await co_return co_yield decltype default delete do double dynamic_cast else enum explicit export extern false float for friend goto if inline int long mutable namespace new noexcept not not_eq nullptr operator or or_eq private protected public register reinterpret_cast requires return short signed sizeof static static_assert static_cast struct switch synchronized template this thread_local throw true try typedef typeid typename union unsigned using virtual void volatile wchar_t while xor xor_eq Query GridValue').split(' '));
function fields(names){let used=new Set(names.filter(x=>!reserved.has(x)&&!x.startsWith('_')&&!x.includes('__')));return names.map((x,i)=>{if(!reserved.has(x)&&!x.startsWith('_')&&!x.includes('__'))return x;let f='argument_'+i;while(used.has(f))f+='_';used.add(f);return f;});}
function queries(result,id,target='queries',withMain=false){
 const plan=result.choices.find(p=>p.id===(id||result.choices[0]?.id));if(!plan)throw Error('適用できる高速計算法がありません');
 const raw=QC.generate(result,plan.id,target),names=result.variables,D=names.length,fs=fields(names),grid=target==='grid';
 const ret=grid?'std::vector<GridValue>':'std::vector<int>',args=grid?'const std::array<long long, '+D+'>& upper_bounds':(plan.array?'const std::vector<long long>& data, ':'')+'const std::vector<Query>& queries';
 const description=grid?'Returns every coordinate and value in [0, upper_bounds], last coordinate varying fastest.':'Returns answers modulo 998244353 in the same order as queries.';
 const api={signature:ret+' solve('+args+')',returnType:ret,args,description,comment:'// C++17 | modulo 998244353 | GCC / Clang (__int128)\n// PROVEN under the algorithm conditions below.\n// Expression: '+result.source.replace(/[\r\n]/g,' ')+'\n// Coordinate order: '+names.join(', '),example:grid?'auto values = solve(std::array<long long, '+D+'>{'+Array(D).fill('10').join(', ')+'});':(plan.array?'std::vector<long long> data = {2, 1, 2, 3, 1, 0, 4, 2};\n':'')+'std::vector<Query> queries = {{'+Array(D).fill('2').join(', ')+'}, {'+Array(D).fill('5').join(', ')+'}};\nauto answers = solve('+(plan.array?'data, ':'')+'queries);'};
 let types=grid?'struct GridValue { std::array<long long, '+D+'> coordinates; int value; };':'struct Query {\n'+fs.map((f,i)=>(f===names[i]?'':'// Expression variable: '+names[i]+'\n')+'long long '+f+';').join('\n')+'\n};';
 let bridge=ret+' solve_impl('+args+') { std::vector<sigma_detail::Query> work;';
 if(grid){bridge+='std::array<long long, '+D+'> stride; long long total=1; for(int d='+D+'-1;d>=0;--d){if(upper_bounds[d]<0||upper_bounds[d]>=1000000||upper_bounds[d]+1>1000000/total)throw std::runtime_error("table limit: 1000000 cells");stride[d]=total;total*=upper_bounds[d]+1;} work.reserve(total);for(long long flat=0;flat<total;++flat)work.push_back({'+fs.map((_,i)=>'(flat/stride['+i+'])%(upper_bounds['+i+']+1)').join(',')+',int(flat)});auto values=solve_queries(work);std::vector<GridValue> answer;answer.reserve(values.size());for(size_t flat=0;flat<values.size();++flat)answer.push_back({{'+fs.map((_,i)=>'static_cast<long long>((flat/stride['+i+'])%(upper_bounds['+i+']+1))').join(',')+'},values[flat].v});return answer;}';}
 else{bridge+='if(queries.size()>1000000)throw std::runtime_error("Q exceeds 1000000");work.reserve(queries.size());for(const auto& query:queries)work.push_back({'+fs.map(f=>'query.'+f).join(',')+',0});auto values=solve_queries(work'+(plan.array?',data':'')+');std::vector<int> answer;answer.reserve(values.size());for(auto value:values)answer.push_back(value.v);return answer;}';}
 // Qualify the public Query; the private implementation has its own indexed Query.
 const publicArgs=args.replaceAll('std::vector<Query>','std::vector<::Query>');bridge=bridge.replace(ret+' solve_impl('+args+')',ret+' solve_impl('+publicArgs+')');api.args=publicArgs;
 let main='';
 if(grid)main='std::array<long long, '+D+'> upper_bounds;for(auto& bound:upper_bounds)if(!(std::cin>>bound))throw std::runtime_error("expected upper bounds");for(const auto& row:solve(upper_bounds)){for(auto coordinate:row.coordinates)std::cout<<coordinate<<" ";std::cout<<row.value<<"\\n";}';
 else{main=plan.array?'int length;if(!(std::cin>>length)||length<=0||length>=(1<<21))throw std::runtime_error("invalid array length");std::vector<long long> data(length);for(auto& value:data)if(!(std::cin>>value))throw std::runtime_error("expected array value");':'';main+='int count;if(!(std::cin>>count)||count<0||count>1000000)throw std::runtime_error("invalid query count");std::vector<Query> queries(count);for(auto& query:queries)if(!(std::cin>>'+fs.map(f=>'query.'+f).join('>>')+'))throw std::runtime_error("expected query coordinates");for(auto answer:solve('+(plan.array?'data, ':'')+'queries))std::cout<<answer<<"\\n";';}
 return assemble(raw,api,types,bridge,'return sigma_detail::solve_impl('+(grid?'upper_bounds':(plan.array?'data, ':'')+'queries')+');',main,withMain);
}
function highlight(source){const escape=s=>s.replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));return source.replace(/\/\/[^\n]*|"(?:\\.|[^"\\])*"|'(?:\\.|[^'\\])*'|\b(?:int|long|const|auto|return|if|else|for|while|struct|namespace|using|void|bool|throw|try|catch|static|constexpr|template|typename|class|sizeof|false|true)\b|\b\d+(?:ULL|LL)?\b|[&<>]/g,t=>{let kind=t.startsWith('//')?'comment':/^["']/.test(t)?'string':/^\d/.test(t)?'number':/^[a-z]/.test(t)?'keyword':'';return kind?'<span class="syntax-'+kind+'">'+escape(t)+'</span>':escape(t);});}
root.SigmaCodeOutput={sequence,queries,format,tokens,highlight,fields};if(typeof module!=='undefined')module.exports=root.SigmaCodeOutput;
})(typeof globalThis!=='undefined'?globalThis:this);
