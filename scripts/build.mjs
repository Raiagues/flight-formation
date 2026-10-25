import {cpSync,mkdirSync,rmSync,writeFileSync,readFileSync,existsSync} from 'node:fs';
const root=new URL('../',import.meta.url),dist=new URL('dist/',root);
rmSync(dist,{recursive:true,force:true});mkdirSync(dist,{recursive:true});
cpSync(new URL('web/',root),dist,{recursive:true});
cpSync(new URL('research/relatorio.pdf',root),new URL('relatorio.pdf',dist));
writeFileSync(new URL('.nojekyll',dist),'');
for(const p of ['index.html','app.js','physics.js','style.css','assets/earth-day.jpg','assets/earth-night.jpg','assets/cubesat-mirata.glb','vendor/draco/draco_decoder.wasm']) {
  if(!existsSync(new URL(p,dist))||!readFileSync(new URL(p,dist)).length)throw new Error(`Missing asset: ${p}`);
}
console.log('Static site ready in dist/; all assets served locally.');
