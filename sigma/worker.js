importScripts('engine.js','optimize.js');
self.onmessage=event=>{try{self.postMessage({result:SigmaEngine.analyze(event.data.source,event.data.options)})}catch(e){self.postMessage({error:e.message})}};
