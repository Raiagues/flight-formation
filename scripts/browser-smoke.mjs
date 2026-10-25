import {mkdirSync} from 'node:fs';
import assert from 'node:assert/strict';
const {chromium}=await import(process.env.PLAYWRIGHT_MODULE||'playwright');
const browser=await chromium.launch({headless:true,args:['--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
const errors=[],page=await browser.newPage({viewport:{width:1440,height:1080},deviceScaleFactor:1});
page.on('pageerror',e=>errors.push(e.message));
page.on('console',m=>{if(m.type()==='error')errors.push(m.text());});
page.on('response',r=>{if(r.status()>=400)errors.push(`${r.status()} ${r.url()}`);});
mkdirSync('test-results',{recursive:true});
try {
  await page.goto(process.env.PREVIEW_URL||'http://127.0.0.1:4173',{waitUntil:'networkidle'});
  await page.waitForFunction(()=>!!window.flightFormation,{timeout:45000});
  await page.waitForFunction(()=>window.flightFormation.snapshot().samples>1);
  assert.equal(await page.locator('#second').isChecked(),false);
  assert.equal(await page.locator('#axes').isChecked(),false);
  await page.locator('#play').click();
  const snapshot=()=>page.evaluate(()=>window.flightFormation.snapshot());
  const paused=await snapshot();await page.waitForTimeout(250);assert.equal((await snapshot()).time,paused.time);
  assert.ok(paused.modelLoaded);
  await page.screenshot({path:'test-results/earth-desktop.png',fullPage:true});
  await page.locator('#axes').check();await page.locator('#second').check();
  await page.locator('#view-relative').click();await page.waitForTimeout(250);
  const relative=await snapshot();assert.equal(relative.view,'relative');assert.equal(relative.satellites,3);assert.equal(relative.axes,true);
  await page.screenshot({path:'test-results/formation-desktop.png',fullPage:true});
  await page.locator('#noise').evaluate(el=>{
    el.value='50';el.dispatchEvent(new Event('input',{bubbles:true}));el.dispatchEvent(new Event('change',{bubbles:true}));
  });
  assert.equal((await snapshot()).sigma,50);assert.equal((await snapshot()).time,0);
  await page.locator('#scenario').selectOption('circular');assert.ok(Math.abs((await snapshot()).state[6]-6978140)<1);
  await page.locator('#play').click();await page.waitForFunction(()=>window.flightFormation.snapshot().time>5);
  await page.locator('#chart-axis').selectOption('2');await page.locator('#trails').uncheck();await page.locator('#estimation').uncheck();
  await page.locator('#view-earth').click();await page.locator('#camera-reset').click();
  await page.locator('#second').uncheck();assert.equal((await snapshot()).satellites,2);
  await page.setViewportSize({width:390,height:844});await page.waitForTimeout(350);
  assert.equal(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth),true);
  await page.screenshot({path:'test-results/mobile.png',fullPage:true});
  assert.deepEqual(errors,[]);
  console.log('Browser smoke passed: WebGL, NASA model, pause, views, axes, follower, noise, scenario, chart, mobile.');
} finally {await browser.close();}
