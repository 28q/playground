const BUILD='20261004-7';
importScripts('engine.js?v='+BUILD,'optimize.js?v='+BUILD,'queries.js?v='+BUILD);
self.onmessage=event=>{try{self.postMessage({build:BUILD,result:SigmaQueries.run(event.data.source,event.data.options)})}catch(e){self.postMessage({build:BUILD,error:e.message})}};
