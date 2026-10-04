const BUILD='20261004-3';
importScripts('engine.js?v='+BUILD,'optimize.js?v='+BUILD);
self.onmessage=event=>{try{self.postMessage({build:BUILD,result:SigmaEngine.analyze(event.data.source,event.data.options)})}catch(e){self.postMessage({build:BUILD,error:e.message})}};
