const BUILD='20261006-2';
importScripts('engine.js?v='+BUILD,'optimize.js?v='+BUILD,'modular.js?v='+BUILD,'rewrite.js?v='+BUILD,'queries.js?v='+BUILD);
self.onmessage=event=>{try{self.postMessage({build:BUILD,result:SigmaQueries.run(event.data.source,event.data.options)})}catch(e){self.postMessage({build:BUILD,error:e.message})}};
