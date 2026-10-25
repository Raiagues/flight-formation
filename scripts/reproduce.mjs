import {readFileSync,writeFileSync} from 'node:fs';
import {EKF, initialState} from '../web/physics.js';
const root=new URL('../original/NOVO/CODIGOS/NAO_LINEAR/tests/CODIGO_PROFESSOR/',import.meta.url);
const read=name=>readFileSync(new URL(name,root),'utf8').trim().split(/\r?\n/).map(l=>l.trim().split(/\s+/).map(Number));
const data=read('data.dat'), truth=read('data_verdade.dat');
const rms=a=>Math.sqrt(a.reduce((s,v)=>s+v*v,0)/a.length);
const std=a=>{const m=a.reduce((s,v)=>s+v,0)/a.length;return rms(a.map(v=>v-m));};
function run(legacy) {
  const filter=new EKF(initialState()), rows=[];
  for(let i=0;i<data.length;i++) {
    if(!legacy) filter.predict(data[i][0]-(data[i-1]?.[0]??0));
    const m=data[i].slice(1), r=filter.update(m);
    rows.push({innovation:r.innovation,correction:r.correction,
      measurementError:m.map((v,j)=>v-truth[i][j+1]),
      estimateError:filter.state.slice(0,3).map((v,j)=>v-truth[i][j+1])});
    if(legacy) filter.predict(5);
  }
  // MATLAB uses 500:end (one-based), i.e. index 499 here.
  const tail=rows.slice(499), stats={};
  for(const key of ['innovation','correction','measurementError','estimateError'])
    stats[key]={rms:[0,1,2].map(j=>rms(tail.map(r=>r[key][j]))),std:[0,1,2].map(j=>std(tail.map(r=>r[key][j])))};
  return stats;
}
const historical=run(true), aligned=run(false);
const output={samples:data.length,firstTime:data[0][0],lastTime:data.at(-1)[0],
  retainedFromSample:500,report:{measured:[15.2502,16.2529,14.9726],estimated:[3.1836,4.7768,2.9203]},
  historical,aligned,
  notes:['Historical: update at t=0 with the first measurement stamped t=5; innovation and correction are the metrics labeled RMSE in the original.',
    'Aligned: predict to measurement timestamp; measure errors against data_verdade.dat.',
    'RK4 dt<=5 s replaces MATLAB ode45. Input ASCII files round state values. Original files remain unchanged.']};
writeFileSync(new URL('../research/reproduction.json',import.meta.url),JSON.stringify(output,null,2)+'\n');
console.log(JSON.stringify(output,null,2));
