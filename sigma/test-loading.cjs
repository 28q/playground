const fs=require('fs'),vm=require('vm'),path=require('path'),assert=require('assert/strict');
const dir=__dirname,source='sum(2*k*k+3*k+4,k,0,n)',E=require('./engine.js'),legacy=E.analyze(source);
let reply;const worker=vm.createContext({self:{postMessage:x=>reply=x}});worker.importScripts=(...urls)=>{for(const url of urls)vm.runInContext(fs.readFileSync(path.join(dir,url.split('?')[0]),'utf8'),worker)};
vm.runInContext(fs.readFileSync(path.join(dir,'worker.js'),'utf8'),worker);worker.self.onmessage({data:{source,options:{}}});assert(reply.result.structure);assert.equal(reply.result.fast.kind,'polynomial');assert.equal(String(reply.result.seq[3]),'62');
const script=fs.readFileSync(path.join(dir,'app.js'),'utf8');
async function app(response){const nodes=new Map();function get(id){if(!nodes.has(id))nodes.set(id,{value:({expression:source,start:'0',count:'40','max-order':'8','max-degree':'4','code-mode':'fast-nth',category:'すべて'})[id]||'',hidden:false,textContent:'',innerHTML:'',addEventListener(){},classList:{toggle(){}},setAttribute(){}});return nodes.get(id)}
 const context=vm.createContext({SigmaEngine:E,SigmaFastCPP:{generate:()=> '// specialized'},SigmaCPP:{generate:()=>''},document:{getElementById:get,querySelectorAll:()=>[],querySelector:()=>get('option'),addEventListener(){}},location:{protocol:'http:'},localStorage:{getItem:()=>null,setItem(){}},performance:{now:()=>0},setTimeout,clearTimeout,console:{error(){}},Worker:class{constructor(url){assert(url.includes('?v='))}terminate(){}postMessage(){queueMicrotask(()=>this.onmessage({data:structuredClone(response)}))}}});
 vm.runInContext(script,context);await new Promise(r=>setTimeout(r,10));return {context,get};}
(async()=>{
 const old=await app({result:legacy});assert.equal(old.get('error').hidden,false);assert.match(old.get('error').textContent,/版が一致/);assert.equal(old.get('results').hidden,true);
 // Reproduce the reported exception by rendering a legacy payload without the new guard.
 old.context.legacyPayload=legacy;
 assert.throws(()=>vm.runInContext('result=legacyPayload;render()',old.context),/reading 'type'/);
 const modern=await app(reply);assert.equal(modern.get('error').hidden,true);assert.equal(modern.get('results').hidden,false);assert.match(modern.get('identity-note').textContent,/変形規則/);assert.match(modern.get('structure-note').textContent,/C-finite/);assert.match(modern.get('cpp-code').textContent,/specialized/);
 const incomplete=await app({build:reply.build,result:legacy});assert.match(incomplete.get('error').textContent,/必要な情報/);assert.equal(incomplete.get('results').hidden,true);
 console.log('PASS: reproduced legacy payload type error; version rejection, incomplete payload rejection, current worker polynomial rendering');
})().catch(e=>{console.error(e);process.exitCode=1});

