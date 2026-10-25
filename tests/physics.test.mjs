import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {MU,A,EARTH_R,PERIOD,initialState,dynamics,jacobian,propagate,toECI,astronomy,Simulation} from '../web/physics.js';
test('nonlinear LVLH dynamics agree with independent inertial two-body integration',()=>{
  let s=initialState(); const [x,y,z,vx,vy,vz,r,vr,,w]=s;
  let inertial=[r+x,y,z,vr+vx-w*y,vy+w*(r+x),vz];
  function f(q) {const g=-MU/Math.hypot(...q.slice(0,3))**3; return [...q.slice(3),...q.slice(0,3).map(v=>v*g)];}
  function step(q,h) {
    const a=f(q),b=f(q.map((v,i)=>v+h*a[i]/2)),c=f(q.map((v,i)=>v+h*b[i]/2)),d=f(q.map((v,i)=>v+h*c[i]));
    return q.map((v,i)=>v+h*(a[i]+2*b[i]+2*c[i]+d[i])/6);
  }
  for(let t=0;t<23000;t+=5) {
    s=propagate(s,5); for(let k=0;k<5;k++)inertial=step(inertial,1);
    if(t%1000===0) {
      const c=Math.cos(s[8]),sn=Math.sin(s[8]);
      const p=[(s[6]+s[0])*c-s[1]*sn,(s[6]+s[0])*sn+s[1]*c,s[2]];
      assert.ok(Math.hypot(...p.map((v,i)=>v-inertial[i]))<.3);
    }
  }
});
test('leader conserves energy and angular momentum; original perigee is below surface',()=>{
  const s=initialState(),end=propagate(s,PERIOD),energy=q=>.5*(q[7]**2+(q[6]*q[9])**2)-MU/q[6];
  assert.ok(Math.abs((energy(end)-energy(s))/energy(s))<1e-9);
  assert.ok(Math.abs(end[6]-s[6])<.03);
  assert.ok(Math.abs(end[6]**2*end[9]/(s[6]**2*s[9])-1)<1e-9);
  assert.ok(Math.abs(A*.9-EARTH_R+97814)<1e-6);
});
test('analytic Jacobian matches central finite differences',()=>{
  const s=propagate(initialState(),800),F=jacobian(s);
  for(let j=0;j<10;j++) {
    const h=j===9?1e-8:j===8?1e-5:Math.max(Math.abs(s[j])*1e-6,1e-3);
    const lo=s.slice(),hi=s.slice();lo[j]-=h;hi[j]+=h;
    const a=dynamics(lo),b=dynamics(hi);
    for(let i=0;i<10;i++) {
      const numeric=(b[i]-a[i])/(2*h);
      assert.ok(Math.abs(F[i*10+j]-numeric)<1e-6*Math.max(1,Math.abs(numeric)),`F[${i},${j}]`);
    }
  }
});
test('trajectory agrees with archived MATLAB truth (rounded ASCII)',()=>{
  const rows=readFileSync(new URL('../original/NOVO/CODIGOS/NAO_LINEAR/tests/CODIGO_PROFESSOR/data_verdade.dat',import.meta.url),'utf8').trim().split(/\r?\n/);
  let s=initialState();let last=0,maxError=0;
  for(const line of rows) {
    const v=line.trim().split(/\s+/).map(Number);s=propagate(s,v[0]-last);last=v[0];
    maxError=Math.max(maxError,...s.slice(0,3).map((x,i)=>Math.abs(x-v[i+1])));
  }
  assert.ok(maxError<.02,`max error ${maxError} m`);
});
test('EKF reduces actual position error, covariance remains symmetric, noise zero is finite',()=>{
  for(const sigma of [0,14.1421356237,50]) {
    const sim=new Simulation({sigma});let m=0,e=0;
    for(let i=0;i<700;i++) {
      const row=sim.step();if(i>100)for(let j=0;j<3;j++){m+=(row.measurement[j]-row.truth[j])**2;e+=(row.estimate[j]-row.truth[j])**2;}
      assert.ok(sim.filter.state.every(Number.isFinite));
    }
    if(sigma>0)assert.ok(e<m*.25,`sigma ${sigma}: ${e/m}`);
    for(let i=0;i<10;i++)for(let j=0;j<10;j++) {
      assert.ok(Math.abs(sim.filter.P[i*10+j]-sim.filter.P[j*10+i])<1e-10);
      if(i===j)assert.ok(sim.filter.P[i*10+i]>=0);
    }
  }
});
test('second follower is independent and deterministic; ECI distance preserves relative norm',()=>{
  const a=new Simulation(),b=new Simulation({second:true}),c=new Simulation();
  a.step();b.step();c.step();assert.deepEqual(a.history,c.history);
  assert.equal(a.truth[2],-b.truth[2]);
  const p=toECI(a.truth),leader=toECI(a.truth,[0,0,0]);
  assert.ok(Math.abs(Math.hypot(...p.map((v,i)=>v-leader[i]))-Math.hypot(...a.truth.slice(0,3)))<1e-8);
});
test('astronomy: J2000 sidereal angle and solar solstices',()=>{
  const j=astronomy(new Date('2000-01-01T12:00:00Z'));
  assert.ok(Math.abs(j.gmst*180/Math.PI-280.46061837)<1e-8);
  const summer=astronomy(new Date('2026-06-21T12:00:00Z')),winter=astronomy(new Date('2026-12-21T12:00:00Z'));
  assert.ok(summer.sun[2]>.397 && winter.sun[2]<-.397);
  assert.ok(Math.abs(Math.hypot(...summer.sun)-1)<1e-12);
});
