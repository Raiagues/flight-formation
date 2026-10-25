import {readFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
const root=new URL('../',import.meta.url);
const manifest=JSON.parse(readFileSync(new URL('research/original-sha256.json',root)));
for(const [name,expected] of Object.entries(manifest)) {
  const actual=createHash('sha256').update(readFileSync(new URL(`original/${name}`,root))).digest('hex');
  if(actual!==expected)throw new Error(`Original file changed: ${name}`);
}
console.log(`${Object.keys(manifest).length} original files verified byte-for-byte.`);
