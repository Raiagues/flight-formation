// SI units; state = [x,y,z,vx,vy,vz,r,rdot,theta,thetadot].
// Nonlinear relative dynamics, report equations (1)–(5).
export const MU = 3.986004418e14;
export const EARTH_R = 6378140;
export const A = EARTH_R + 600e3;
export const PERIOD = 2 * Math.PI * Math.sqrt(A ** 3 / MU);
export const SAMPLE = 5;
export const DEFAULT_SIGMA = 10 * Math.sqrt(2);
export function initialState(e = 0.1, second = false) {
  const n = Math.sqrt(MU / A ** 3);
  return [-0.01127*A, 0, (second ? -0.1 : 0.1)*A,
    0.02*n*A, 0.02*n*A, 0, A*(1-e), 0, 0,
    Math.sqrt(MU/(A*(1-e*e))**3)*(1+e)**2];
}
export function dynamics(s) {
  const [x,y,z,vx,vy,vz,r,vr,,w] = s;
  const g = MU / ((r+x)**2+y*y+z*z)**1.5;
  const wd = -2*vr*w/r;
  return [vx,vy,vz, 2*w*vy+wd*y+w*w*x-g*(r+x)+MU/(r*r),
    -2*w*vx-wd*x+w*w*y-g*y, -g*z, vr, r*w*w-MU/(r*r), w, wd];
}
export function jacobian(s) {
  const [x,y,z,vx,vy,,r,vr,,w] = s;
  const q = [r+x,y,z], d2=q.reduce((a,v)=>a+v*v,0);
  const g=MU/d2**1.5, b=3*MU/d2**2.5, F=new Float64Array(100);
  F[3]=F[14]=F[25]=F[67]=F[89]=1;
  for(let i=0;i<3;i++) for(let j=0;j<3;j++)
    F[(i+3)*10+j]=b*q[i]*q[j]-(i===j?g:0);
  F[30]+=w*w; F[31]-=2*vr*w/r; F[34]=2*w;
  F[40]+=2*vr*w/r; F[41]+=w*w; F[43]=-2*w;
  F[36]=b*q[0]*q[0]-g-2*MU/r**3+2*vr*w*y/r**2;
  F[37]=-2*w*y/r; F[39]=2*vy+2*w*x-2*vr*y/r;
  F[46]=b*y*q[0]-2*vr*w*x/r**2;
  F[47]=2*w*x/r; F[49]=-2*vx+2*w*y+2*vr*x/r;
  F[56]=b*z*q[0];
  F[76]=w*w+2*MU/r**3; F[79]=2*r*w;
  F[96]=2*vr*w/r**2; F[97]=-2*w/r; F[99]=-2*vr/r;
  return F;
}
export function rk4(s, dt, f = dynamics) {
  const a=f(s), b=f(s.map((v,i)=>v+a[i]*dt/2));
  const c=f(s.map((v,i)=>v+b[i]*dt/2)), d=f(s.map((v,i)=>v+c[i]*dt));
  return s.map((v,i)=>v+dt*(a[i]+2*b[i]+2*c[i]+d[i])/6);
}
export function propagate(s, dt) {
  let out = Array.from(s);
  while(dt>1e-9) { const h=Math.min(dt,5); out=rk4(out,h); dt-=h; }
  return out;
}
export const PROCESS_Q = [1,1,1,1e-4,1e-4,1e-4,1e-10,1e-40,(1e-20*Math.PI/180)**2,(1e-20*Math.PI/180)**2];
export function initialCovariance() {
  const p=new Float64Array(100), sig=[20,20,20,.1,.1,.1,1e-20,1e-20,1e-20*Math.PI/180,1e-20*Math.PI/180];
  sig.forEach((v,i)=>p[i*10+i]=v*v); return p;
}
function augmentedDynamics(s) {
  const f=dynamics(s), F=jacobian(s), out=new Array(110).fill(0);
  for(let i=0;i<10;i++) {
    out[i]=f[i];
    for(let j=0;j<10;j++) {
      let v=i===j?PROCESS_Q[i]:0;
      for(let k=0;k<10;k++) v+=F[i*10+k]*s[10+k*10+j]+s[10+i*10+k]*F[j*10+k];
      out[10+i*10+j]=v;
    }
  }
  return out;
}
export class EKF {
  constructor(state=initialState(), sigma=DEFAULT_SIGMA) {
    this.state=Array.from(state); this.P=initialCovariance(); this.sigma=sigma;
  }
  predict(dt=SAMPLE) {
    let s=[...this.state,...this.P];
    while(dt>1e-9) { const h=Math.min(5,dt); s=rk4(s,h,augmentedDynamics); dt-=h; }
    this.state=s.slice(0,10); this.P=Float64Array.from(s.slice(10));
  }
  update(measurement) {
    // Scalar sequential updates are equivalent to a vector update for diagonal R.
    // Joseph form protects covariance symmetry and positive semidefiniteness.
    const innovation=measurement.map((v,i)=>v-this.state[i]);
    const before=this.state.slice();
    const R=Math.max(this.sigma**2,1e-12);
    for(let axis=0;axis<3;axis++) {
      const P=this.P, denom=P[axis*10+axis]+R;
      const K=Array.from({length:10},(_,i)=>P[i*10+axis]/denom);
      const res=measurement[axis]-this.state[axis];
      this.state=this.state.map((v,i)=>v+K[i]*res);
      const next=new Float64Array(100);
      for(let i=0;i<10;i++) for(let j=0;j<10;j++)
        next[i*10+j]=P[i*10+j]-K[i]*P[axis*10+j]-P[i*10+axis]*K[j]+K[i]*denom*K[j];
      this.P=next;
    }
    return {innovation,correction:this.state.slice(0,3).map((v,i)=>v-before[i])};
  }
}
export function normalRandom(seed=2023) {
  let n=seed>>>0;
  const uniform=()=>{n=(Math.imul(1664525,n)+1013904223)>>>0; return (n+.5)/4294967296;};
  return ()=>Math.sqrt(-2*Math.log(uniform()))*Math.cos(2*Math.PI*uniform());
}
export class Simulation {
  constructor({e=.1,sigma=DEFAULT_SIGMA,second=false}={}) {
    this.e=e; this.sigma=sigma; this.second=second; this.t=0;
    this.truth=initialState(e,second); this.filter=new EKF(this.truth,sigma);
    this.random=normalRandom(second?2024:2023); this.history=[];
  }
  step() {
    this.truth=propagate(this.truth,SAMPLE); this.t+=SAMPLE;
    this.filter.predict();
    const measurement=this.truth.slice(0,3).map(v=>v+this.sigma*this.random());
    this.filter.update(measurement);
    const row={t:this.t,truth:this.truth.slice(0,3),measurement,estimate:this.filter.state.slice(0,3)};
    this.history.push(row); if(this.history.length>240) this.history.shift();
    return row;
  }
}
// Orientation is only a display choice: the original has no inclination/RAAN.
export const INCLINATION = 35*Math.PI/180;
export function toECI(s, relative=s.slice(0,3)) {
  const c=Math.cos(s[8]), sn=Math.sin(s[8]);
  const X=(s[6]+relative[0])*c-relative[1]*sn;
  const Y=(s[6]+relative[0])*sn+relative[1]*c;
  return [X,Y*Math.cos(INCLINATION)-relative[2]*Math.sin(INCLINATION),
    Y*Math.sin(INCLINATION)+relative[2]*Math.cos(INCLINATION)];
}
export function astronomy(date=new Date()) {
  const jd=date.getTime()/86400000+2440587.5, d=jd-2451545, rad=Math.PI/180;
  const wrap=x=>((x%360)+360)%360;
  const g=wrap(357.529+.98560028*d)*rad;
  const lambda=(wrap(280.459+.98564736*d)+1.915*Math.sin(g)+.020*Math.sin(2*g))*rad;
  const eps=(23.439-.00000036*d)*rad;
  return {gmst:wrap(280.46061837+360.98564736629*d)*rad,
    sun:[Math.cos(lambda),Math.cos(eps)*Math.sin(lambda),Math.sin(eps)*Math.sin(lambda)]};
}
