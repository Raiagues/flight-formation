import * as THREE from './vendor/three.module.js';
import {OrbitControls} from './vendor/OrbitControls.js';
import {GLTFLoader} from './vendor/GLTFLoader.js';
import {DRACOLoader} from './vendor/DRACOLoader.js';
import {RoomEnvironment} from './vendor/RoomEnvironment.js';
import {A,EARTH_R,PERIOD,SAMPLE,DEFAULT_SIGMA,Simulation,initialState,propagate,toECI,astronomy,normalRandom} from './physics.js';

const $=id=>document.getElementById(id);
const format=(x,d=1)=>x.toLocaleString('pt-BR',{minimumFractionDigits:d,maximumFractionDigits:d});
let sim=new Simulation(),sim2=new Simulation({second:true}),view='earth',playing=!matchMedia('(prefers-reduced-motion: reduce)').matches;
let time=0,speed=120,last=performance.now(),lastUI=0,renderer,scene,camera,controls;
let earth,atmosphere,sunLight,sunSprite,stars,grid,leader,follower,second,relativeTrail,secondTrail,orbitTrail;
let sunDirection=new THREE.Vector3(),relativePath=[],path2=[],modelPrototype;
let measurementMarker,estimateMarker;
const eciVector=p=>new THREE.Vector3(p[0],p[2],-p[1]);
const scaleECI=s=>eciVector(toECI(s)).multiplyScalar(1/EARTH_R);
let toastTimer;
function toast(text) { $('message').textContent=text;$('message').hidden=false;clearTimeout(toastTimer);toastTimer=setTimeout(()=>$('message').hidden=true,4300); }
function makeLine(color,opacity=.5) {return new THREE.Line(new THREE.BufferGeometry(),new THREE.LineBasicMaterial({color,transparent:true,opacity}));}
function setPoints(line,points) {line.geometry.dispose();line.geometry=new THREE.BufferGeometry().setFromPoints(points);}
function basis(s) {
  const c=Math.cos(s[8]),sn=Math.sin(s[8]),inc=35*Math.PI/180;
  return new THREE.Matrix4().makeBasis(eciVector([c,sn*Math.cos(inc),sn*Math.sin(inc)]),
    eciVector([-sn,c*Math.cos(inc),c*Math.sin(inc)]),eciVector([0,-Math.sin(inc),Math.cos(inc)]));
}
function relativeVector(p) {return new THREE.Vector3(p[0],p[2],-p[1]).multiplyScalar(4/A);}
async function init() {
  try {renderer=new THREE.WebGLRenderer({antialias:true,alpha:true,powerPreference:'high-performance'});}
  catch {throw new Error('A visualização precisa de WebGL 2. Ative a aceleração gráfica ou abra em outro navegador.');}
  renderer.setPixelRatio(Math.min(devicePixelRatio,2));renderer.setClearColor(0x080d15,0);
  renderer.toneMapping=THREE.ACESFilmicToneMapping;renderer.toneMappingExposure=1.15;
  $('viewport').appendChild(renderer.domElement);
  renderer.domElement.addEventListener('webglcontextlost',event=>{event.preventDefault();playing=false;updatePlay();toast('A conexão gráfica foi interrompida. Recarregue a página para retomar.');});
  scene=new THREE.Scene();camera=new THREE.PerspectiveCamera(39,1,.001,300);
  controls=new OrbitControls(camera,renderer.domElement);controls.enableDamping=true;controls.dampingFactor=.075;controls.enablePan=true;controls.maxDistance=8;controls.minDistance=1.15;
  const pmrem=new THREE.PMREMGenerator(renderer);const env=new RoomEnvironment();scene.environment=pmrem.fromScene(env,.04).texture;env.dispose();pmrem.dispose();
  scene.add(new THREE.AmbientLight(0xb0c4e6,.15));
  sunLight=new THREE.DirectionalLight(0xfff3dc,3.3);scene.add(sunLight);
  const textureLoader=new THREE.TextureLoader();
  const draco=new DRACOLoader().setDecoderPath('./vendor/draco/');
  const modelLoader=new GLTFLoader().setDRACOLoader(draco);
  const [day,night,gltf]=await Promise.all([
    textureLoader.loadAsync('./assets/earth-day.jpg'),textureLoader.loadAsync('./assets/earth-night.jpg'),
    modelLoader.loadAsync('./assets/cubesat-mirata.glb')]);
  draco.dispose();
  for(const t of [day,night]) {t.colorSpace=THREE.SRGBColorSpace;t.anisotropy=Math.min(8,renderer.capabilities.getMaxAnisotropy());}
  const material=new THREE.ShaderMaterial({uniforms:{dayMap:{value:day},nightMap:{value:night},sun:{value:new THREE.Vector3()}},
    vertexShader:`varying vec2 vUv;varying vec3 vN;void main(){vUv=uv;vN=normalize(mat3(modelMatrix)*normal);gl_Position=projectionMatrix*modelViewMatrix*vec4(position,1.);}`,
    fragmentShader:`uniform sampler2D dayMap;uniform sampler2D nightMap;uniform vec3 sun;varying vec2 vUv;varying vec3 vN;
    void main(){float d=dot(normalize(vN),sun);float light=smoothstep(-.07,.13,d);vec3 day=texture2D(dayMap,vUv).rgb;vec3 night=texture2D(nightMap,vUv).rgb;
      vec3 color=day*(.022+1.25*max(d,0.))*light+night*(1.-light)*1.6+day*.007;
      gl_FragColor=vec4(color,1.);
      #include <tonemapping_fragment>
      #include <colorspace_fragment>
    }`});
  earth=new THREE.Mesh(new THREE.SphereGeometry(1,128,80),material);scene.add(earth);
  atmosphere=new THREE.Mesh(new THREE.SphereGeometry(1.025,96,64),new THREE.ShaderMaterial({transparent:true,depthWrite:false,side:THREE.BackSide,blending:THREE.AdditiveBlending,
    uniforms:{sun:{value:new THREE.Vector3()}},vertexShader:`varying vec3 vN;varying vec3 vP;void main(){vN=normalize(mat3(modelMatrix)*normal);vP=(modelMatrix*vec4(position,1.)).xyz;gl_Position=projectionMatrix*modelViewMatrix*vec4(position,1.);}`,
    fragmentShader:`uniform vec3 sun;varying vec3 vN;varying vec3 vP;void main(){float rim=pow(1.-abs(dot(normalize(vN),normalize(cameraPosition-vP))),3.);float light=smoothstep(-.3,.7,dot(normalize(vN),sun));gl_FragColor=vec4(.12,.4,.85,rim*(.12+.48*light));}` }));scene.add(atmosphere);
  // A procedural glow is a light effect; the planet maps and spacecraft are real assets.
  const glow=document.createElement('canvas');glow.width=glow.height=128;const ctx=glow.getContext('2d');
  const gradient=ctx.createRadialGradient(64,64,0,64,64,64);gradient.addColorStop(0,'#fffef2');gradient.addColorStop(.12,'#fff6d7');gradient.addColorStop(.2,'#ffd699aa');gradient.addColorStop(.55,'#ffb85b18');gradient.addColorStop(1,'#ffb85b00');ctx.fillStyle=gradient;ctx.fillRect(0,0,128,128);
  sunSprite=new THREE.Sprite(new THREE.SpriteMaterial({map:new THREE.CanvasTexture(glow),blending:THREE.AdditiveBlending,depthWrite:false}));sunSprite.scale.setScalar(6);scene.add(sunSprite);
  const random=normalRandom(42),positions=[];
  for(let i=0;i<1400;i++){const p=new THREE.Vector3(random(),random(),random()).normalize().multiplyScalar(95);positions.push(p.x,p.y,p.z);}
  const starGeo=new THREE.BufferGeometry();starGeo.setAttribute('position',new THREE.Float32BufferAttribute(positions,3));
  stars=new THREE.Points(starGeo,new THREE.PointsMaterial({color:0x879db8,size:.08,sizeAttenuation:true,transparent:true,opacity:.5}));scene.add(stars);
  grid=new THREE.GridHelper(2,12,0x253d51,0x152433);grid.position.y=-.55;grid.material.transparent=true;grid.material.opacity=.35;grid.visible=false;scene.add(grid);
  modelPrototype=gltf.scene;const box=new THREE.Box3().setFromObject(modelPrototype),size=box.getSize(new THREE.Vector3()),center=box.getCenter(new THREE.Vector3());
  modelPrototype.position.sub(center);const normalized=new THREE.Group();normalized.add(modelPrototype);normalized.scale.setScalar(1/Math.max(size.x,size.y,size.z));modelPrototype=normalized;
  leader=spacecraft(0x90baff);follower=spacecraft(0x8ce4c1);second=spacecraft(0xedc487);
  orbitTrail=makeLine(0x90baff,.32);relativeTrail=makeLine(0x8ce4c1,.65);secondTrail=makeLine(0xedc487,.5);scene.add(orbitTrail,relativeTrail,secondTrail);
  measurementMarker=new THREE.Mesh(new THREE.SphereGeometry(.004,8,8),new THREE.MeshBasicMaterial({color:0x829bc7,transparent:true,opacity:.8}));
  estimateMarker=new THREE.Mesh(new THREE.SphereGeometry(.003,8,8),new THREE.MeshBasicMaterial({color:0xa2fbd4}));scene.add(measurementMarker,estimateMarker);
  rebuildPaths();resetCamera();resize();new ResizeObserver(resize).observe($('viewport'));
  $('loading').hidden=true;bind();updatePlay();requestAnimationFrame(frame);
  // Read-only diagnostics for the browser smoke test, no remote telemetry.
  window.flightFormation={snapshot:()=>({view,time,second:$('second').checked,axes:$('axes').checked,
    satellites:$('second').checked?3:2,sigma:sim.sigma,state:sim.truth.slice(),samples:sim.history.length,
    modelLoaded:!!modelPrototype,earthRotation:earth.rotation.y,sun:sunDirection.toArray(),playing}),
    render:()=>renderer.render(scene,camera)};
}
function spacecraft(color) {
  const group=new THREE.Group(),model=modelPrototype.clone(true);model.scale.multiplyScalar(.075);group.add(model);
  model.traverse(o=>{if(o.isMesh){o.material=o.material.clone();o.material.envMapIntensity=.8;}});
  const axes=new THREE.AxesHelper(.16);axes.visible=false;group.add(axes);
  const dot=new THREE.Mesh(new THREE.SphereGeometry(.005,10,8),new THREE.MeshBasicMaterial({color}));group.add(dot);
  group.userData={axes,model,dot};scene.add(group);return group;
}
function rebuildPaths() {
  relativePath=[];path2=[];let a=initialState(sim.e),b=initialState(sim.e,true);
  const leaderPoints=[],followerPoints=[],secondPoints=[];
  for(let i=0;i<=360;i++) {
    relativePath.push(relativeVector(a));path2.push(relativeVector(b));
    leaderPoints.push(eciVector(toECI(a,[0,0,0])).multiplyScalar(1/EARTH_R));followerPoints.push(scaleECI(a));secondPoints.push(scaleECI(b));
    a=propagate(a,PERIOD/360);b=propagate(b,PERIOD/360);
  }
  setPoints(orbitTrail,leaderPoints);setPoints(relativeTrail,view==='earth'?followerPoints:relativePath);setPoints(secondTrail,view==='earth'?secondPoints:path2);
}
function resetCamera() {
  controls.target.set(0,0,0);
  if(view==='earth'){
    const sun=eciVector(astronomy(new Date()).sun);
    const side=new THREE.Vector3().crossVectors(sun,new THREE.Vector3(0,1,0)).normalize();
    camera.position.copy(sun).multiplyScalar(3.6).addScaledVector(side,2.7).add(new THREE.Vector3(0,.8,0));
    controls.minDistance=1.15;controls.maxDistance=8;
  }
  else {camera.position.set(1.15,.9,1.4);controls.target.set(0,0,-.05);controls.minDistance=.12;controls.maxDistance=4;}
  controls.update();
}
function resize() {
  if(!renderer)return;const r=$('viewport').getBoundingClientRect();camera.aspect=r.width/r.height;camera.fov=camera.aspect<1?55:39;camera.updateProjectionMatrix();renderer.setSize(r.width,r.height);drawChart();
}
function updatePlay(){ $('play').innerHTML=playing?'Ⅱ <span>Pausar</span>':'▷ <span>Continuar</span>';$('play').setAttribute('aria-label',playing?'Pausar simulação':'Continuar simulação');$('running-status').textContent=playing?'Em execução':'Pausada'; }
function restart(notify=true) {
  const options={e:$('scenario').value==='original'?.1:0,sigma:Number($('noise').value)};
  sim=new Simulation(options);sim2=new Simulation({...options,second:true});time=0;last=performance.now();
  rebuildPaths();drawChart();if(notify)toast('Simulação reiniciada com os parâmetros selecionados.');
}
function bind() {
  $('play').onclick=()=>{playing=!playing;updatePlay();};$('reset').onclick=()=>restart();
  $('speed').onchange=()=>{speed=Number($('speed').value);};$('camera-reset').onclick=resetCamera;
  for(const mode of ['earth','relative'])$( `view-${mode}`).onclick=()=>{
    view=mode;for(const id of ['earth','relative']){$(`view-${id}`).classList.toggle('selected',id===mode);$(`view-${id}`).setAttribute('aria-pressed',String(id===mode));}
    $('frame-label').textContent=mode==='earth'?'Referencial geocêntrico':'LVLH · referencial fixo ao líder';
    $('scene-note').textContent=mode==='earth'?'Terra em escala · CubeSats ampliados para visualização':'X radial · Y tangencial · Z normal · CubeSats ampliados';
    rebuildPaths();resetCamera();
  };
  $('noise').oninput=()=>{$('noise-value').textContent=`${format(Number($('noise').value),2)} m`;};
  $('noise').onchange=()=>restart();$('scenario').onchange=()=>{restart();resetCamera();};
  $('axes').onchange=()=>{$('axes-key').hidden=!$('axes').checked;};
  $('second').onchange=()=>{$('satellite-count').textContent=$('second').checked?'03 SATÉLITES':'02 SATÉLITES';toast($('second').checked?'Seguidor 02 ativado com condição inicial espelhada em Z.':'Seguidor 02 desativado.');};
  $('chart-axis').onchange=drawChart;
  document.addEventListener('visibilitychange',()=>{last=performance.now();});
  $('viewport').onkeydown=e=>{if(e.code==='Space'){e.preventDefault();playing=!playing;updatePlay();}if(e.key.toLowerCase()==='r')resetCamera();};
}
function place(group,s,relative,isLeader=false) {
  group.position.copy(view==='earth'?eciVector(toECI(s,isLeader?[0,0,0]:relative)).multiplyScalar(1/EARTH_R):isLeader?new THREE.Vector3():relativeVector(relative));
  if(view==='earth')group.quaternion.setFromRotationMatrix(basis(s));
  else group.quaternion.setFromRotationMatrix(new THREE.Matrix4().makeBasis(new THREE.Vector3(1,0,0),new THREE.Vector3(0,0,-1),new THREE.Vector3(0,1,0)));
  group.userData.axes.visible=$('axes').checked;
  // Only the dot is an X-ray locator when the mathematical orbit is inside Earth.
  // Spacecraft geometry remains depth-tested and is never moved to a fake orbit.
  const underground=view==='earth'&&group.position.length()<1;
  group.userData.dot.material.depthTest=!underground;group.userData.dot.renderOrder=underground?5:0;
  group.userData.dot.scale.setScalar(underground?1.6:1);
}
function label(id,p,underground=false) {
  const el=$(id),v=p.clone().project(camera),size=$('viewport').getBoundingClientRect();
  const to=p.clone().sub(camera.position),distance=to.length(),direction=to.normalize();
  const b=camera.position.dot(direction),disc=b*b-camera.position.lengthSq()+1;
  const occluded=view==='earth'&&!underground&&disc>0&&(-b-Math.sqrt(disc)>0)&&(-b-Math.sqrt(disc)<distance-.02);
  el.hidden=v.z>1||v.z<-1||occluded||Math.abs(v.x)>1||Math.abs(v.y)>1;
  el.style.transform=`translate(${(v.x*.5+.5)*size.width+13}px,${(-v.y*.5+.5)*size.height-13}px)`;
}
function frame(now) {
  requestAnimationFrame(frame);const dt=Math.min((now-last)/1000,.1);last=now;
  if(playing&&!document.hidden) {
    time=Math.min(23005,time+dt*speed);
    // Run both independent followers so enabling #02 preserves phase continuity.
    while(sim.t+SAMPLE<=time){sim.step();sim2.step();}
    if(time>=23005){time=23005;playing=false;updatePlay();toast('Fim da janela de 23.005 s. Reinicie para executar novamente.');}
  }
  const s=propagate(sim.truth,Math.max(0,time-sim.t)),s2=propagate(sim2.truth,Math.max(0,time-sim2.t));
  place(leader,s,[0,0,0],true);place(follower,s,s);place(second,s2,s2);
  second.visible=$('second').checked;earth.visible=atmosphere.visible=stars.visible=sunSprite.visible=view==='earth';grid.visible=view==='relative';
  const sky=astronomy(new Date());sunDirection.copy(eciVector(sky.sun));
  earth.rotation.y=sky.gmst;earth.material.uniforms.sun.value.copy(sunDirection);atmosphere.material.uniforms.sun.value.copy(sunDirection);
  sunLight.position.copy(sunDirection).multiplyScalar(40);sunSprite.position.copy(sunDirection).multiplyScalar(70);
  orbitTrail.visible=$('trails').checked&&view==='earth';relativeTrail.visible=$('trails').checked;secondTrail.visible=$('trails').checked&&second.visible;
  const row=sim.history.at(-1),show=$('estimation').checked&&!!row;
  measurementMarker.visible=estimateMarker.visible=show;
  if(show) {
    // Measurement markers use the actual 5 s sample epoch, never noise inflation.
    measurementMarker.position.copy(view==='earth'?eciVector(toECI(sim.truth,row.measurement)).multiplyScalar(1/EARTH_R):relativeVector(row.measurement));
    estimateMarker.position.copy(view==='earth'?eciVector(toECI(sim.truth,row.estimate)).multiplyScalar(1/EARTH_R):relativeVector(row.estimate));
  }
  controls.update();renderer.render(scene,camera);
  label('label-leader',leader.position,leader.position.length()<1);label('label-follower',follower.position,follower.position.length()<1);
  label('label-second',second.position,second.position.length()<1);if(!second.visible)$('label-second').hidden=true;
  label('label-sun',sunSprite.position);if(view!=='earth')$('label-sun').hidden=true;
  if(now-lastUI>100) {
    lastUI=now;$('distance').innerHTML=`${format(Math.hypot(...s.slice(0,3))/1000)}<small> km</small>`;
    $('altitude').innerHTML=`${format((s[6]-EARTH_R)/1000)}<small> km</small>`;
    $('elapsed').textContent=new Date(time*1000).toISOString().slice(11,19);
    $('utc').textContent=new Date().toISOString().replace('T',' · ').slice(0,21)+' UTC';
    $('orbit-note').hidden=sim.e===0;
    $('orbit-note').textContent=view==='earth'&&s[6]<EARTH_R?'Cenário original: líder sob a superfície. Ponto luminoso indica sua posição matemática.':'Cenário original: perigeu de −98 km. Modelo matemático com passagem pelo interior da Terra.';
    drawChart();
  }
}
function drawChart() {
  const canvas=$('error-chart'),rect=canvas.getBoundingClientRect();if(!rect.width)return;
  const ratio=Math.min(devicePixelRatio,2),w=rect.width,h=rect.height;canvas.width=w*ratio;canvas.height=h*ratio;
  const ctx=canvas.getContext('2d');ctx.scale(ratio,ratio);ctx.clearRect(0,0,w,h);
  const rows=sim.history,axis=Number($('chart-axis').value),m=rows.map(r=>r.measurement[axis]-r.truth[axis]),e=rows.map(r=>r.estimate[axis]-r.truth[axis]);
  $('chart-empty').hidden=rows.length>0;
  const max=Math.max(10,...m.map(Math.abs),...e.map(Math.abs))*1.12,left=36,right=w-4,top=17,bottom=h-20;
  ctx.font='8px Consolas,monospace';ctx.strokeStyle='#ffffff0b';ctx.fillStyle='#75859c';
  for(let j=-1;j<=1;j++) {const y=top+(1-j)*(bottom-top)/2;ctx.beginPath();ctx.moveTo(left,y);ctx.lineTo(right,y);ctx.stroke();ctx.fillText(`${Math.round(j*max)}`,0,y+3);}
  ctx.fillText('m',0,7);if(rows.length){ctx.fillText(`${Math.round(rows[0].t)} s`,left,h-3);ctx.textAlign='right';ctx.fillText(`${Math.round(rows.at(-1).t)} s`,right,h-3);ctx.textAlign='left';}
  for(const [values,color] of [[m,'#7892bc'],[e,'#8ce4c1']]) {
    ctx.beginPath();ctx.strokeStyle=color;ctx.lineWidth=1.2;
    values.forEach((v,i)=>{const x=left+i/Math.max(1,values.length-1)*(right-left),y=top+(1-v/max)*(bottom-top)/2;i?ctx.lineTo(x,y):ctx.moveTo(x,y);});ctx.stroke();
  }
  const rmse=values=>values.length?`${format(Math.sqrt(values.reduce((a,v)=>a+v*v,0)/values.length),2)} m`:'—';
  $('rmse-measured').textContent=rmse(m);$('rmse-estimated').textContent=rmse(e);
  $('window-size').textContent=`${rows.length} amostras · eixo ${'XYZ'[axis]}`;
}
init().catch(error=>{
  console.error(error);$('loading').textContent=`${error.message} Se o carregamento falhou, recarregue a página.`;$('loading').classList.add('error');
  for(const id of ['play','reset','speed','scenario','second','axes','trails','estimation','noise','view-earth','view-relative','camera-reset'])$(id).disabled=true;
  $('running-status').textContent='Indisponível';
});
