(function(){
'use strict';
const root=document.getElementById('fluidic-cfd'), $=id=>root.querySelector('#cfd-'+id), D=window.FLUIDIC_CFD;
const pct=n=>(n*100).toFixed(n<1e-5?6:3)+'%', fmt=n=>Math.abs(n)<.00001?'0':n.toFixed(3);
const order=['mix-s1','mix-s2','signed-s1','e0-s1','e1-s1','e2-s1','e3-s1'];
order.filter(n=>D.cases[n]).forEach(n=>{const o=document.createElement('option');o.value=n;o.textContent=D.cases[n].label;$('case').append(o);});
const rows=(a,id,digits)=>{$(id).replaceChildren();a.forEach((row,j)=>{const tr=document.createElement('tr');const th=document.createElement('th');th.scope='row';th.textContent='Q'+(j+1);tr.append(th);row.forEach(v=>{const td=document.createElement('td');td.className='text-end';td.textContent=digits===0?Math.round(v):v.toFixed(digits);tr.append(td);});$(id).append(tr);});};
rows(D.matrix,'matrix',3);rows(D.widths_um,'widths',0);
const V=D.validation;
const checks=[['Mixed-input superposition',pct(V.checks['mix-s1'].relative_L2_superposition_error)],['Signed-input superposition',pct(V.checks['signed-s1'].relative_L2_superposition_error)],['10 → 5 µm grid change',pct(V.grid_refinement.relative_L2_change)],['Square-duct error · 5 µm',pct(V.square_duct[1].relative_error)],['Peak mass imbalance',pct(V.maximum_device_mass_imbalance)],['Largest lattice Mach number',V.maximum_device_mach.toFixed(5)],['Water viscosity','1 mPa·s'],['Water density','1,000 kg/m³'],['Manifold cross-section','120 × 120 µm'],['Device envelope','880 × 880 × 480 µm']];
checks.forEach(([label,value])=>{const r=document.createElement('div');r.className='cfd-check';const a=document.createElement('span'),b=document.createElement('span');a.textContent=label;b.textContent=value;r.append(a,b);$('checks').append(r);});
let renderer,scene,camera,cloud,arrows,outline,labels=[],data=[],picked=-1,theta=-.72,phi=1.01,distance=1620,colorLo,colorHi,colorNeg,minValue,maxValue;
let themeProbe=document.createElement('span');themeProbe.hidden=true;root.append(themeProbe);
const colorCanvas=document.createElement('canvas');colorCanvas.width=colorCanvas.height=1;const colorContext=colorCanvas.getContext('2d');
function color(token){colorContext.clearRect(0,0,1,1);themeProbe.style.color='var(--background)';colorContext.fillStyle=getComputedStyle(themeProbe).color;colorContext.fillRect(0,0,1,1);themeProbe.style.color='var('+token+')';colorContext.fillStyle=getComputedStyle(themeProbe).color;colorContext.fillRect(0,0,1,1);const a=colorContext.getImageData(0,0,1,1).data;return new THREE.Color('rgb('+a[0]+','+a[1]+','+a[2]+')');}
function coords(a){return new THREE.Vector3(a[0]-440,a[2]-240,a[1]-440);}
function cameraUpdate(){if(!renderer)return;const d=distance*Math.max(1,1.35/camera.aspect);camera.position.set(d*Math.sin(phi)*Math.cos(theta),d*Math.cos(phi),d*Math.sin(phi)*Math.sin(theta));camera.lookAt(0,0,0);draw();}
function draw(){if(!renderer)return;renderer.render(scene,camera);const w=$('stage').clientWidth,h=$('stage').clientHeight;labels.forEach(l=>{const p=l.p.clone().project(camera),half=l.el.offsetWidth/2;l.el.style.left=Math.max(half+4,Math.min(w-half-4,(p.x+1)*w/2))+'px';l.el.style.top=Math.max(12,Math.min(h-12,(1-p.y)*h/2))+'px';l.el.hidden=p.z>1||l.height>Number($('slice').value);});}
function dispose(obj){if(!obj)return;scene.remove(obj);if(obj.geometry)obj.geometry.dispose();if(obj.material)obj.material.dispose();}
function recolor(){if(!renderer)return;colorLo=color('--muted');colorHi=color('--viz-series-1');colorNeg=color('--viz-series-2');updateField();outline.material.color.copy(color('--foreground'));draw();}
function value(a){return $('field').value==='pressure'?a[3]:Math.hypot(a[4],a[5],a[6]);}
function scalarColor(v){if(minValue<0&&$('field').value==='pressure'){return colorLo.clone().lerp(v<0?colorNeg:colorHi,Math.min(1,Math.abs(v)/Math.max(Math.abs(minValue),Math.abs(maxValue))));}return colorLo.clone().lerp(colorHi,.13+.87*(v-minValue)/(maxValue-minValue||1));}
function updateField(){
  if(!renderer)return;
  const slice=Number($('slice').value);$('slice-value').textContent=slice===480?'All layers':slice+' µm';
  const active=D.cases[$('case').value];data=active.samples.filter(a=>a[2]<=slice);const vals=active.samples.map(value);minValue=Math.min(0,...vals);maxValue=Math.max(...vals);
  $('low').textContent=fmt(minValue);$('high').textContent=fmt(maxValue);$('field-label').textContent=$('field').value==='pressure'?'Pressure (Pa)':'Speed (mm/s)';
  $('gradient').style.background=minValue<0&&$('field').value==='pressure'?'linear-gradient(to right,var(--viz-series-2),var(--muted),var(--viz-series-1))':'linear-gradient(to right,var(--muted),var(--viz-series-1))';
  const positions=[],colors=[];data.forEach(a=>{const p=coords(a),c=scalarColor(value(a));positions.push(p.x,p.y,p.z);colors.push(c.r,c.g,c.b);});
  dispose(cloud);const g=new THREE.BufferGeometry();g.setAttribute('position',new THREE.Float32BufferAttribute(positions,3));g.setAttribute('color',new THREE.Float32BufferAttribute(colors,3));cloud=new THREE.Points(g,new THREE.PointsMaterial({size:15,vertexColors:true,sizeAttenuation:true,transparent:true,opacity:.93}));scene.add(cloud);
  dispose(arrows);const vectors=[],vcolors=[];const peak=Math.max(...active.samples.map(a=>Math.hypot(a[4],a[5],a[6])));
  if($('vectors').checked)data.forEach((a,i)=>{if(i%9)return;const mag=Math.hypot(a[4],a[5],a[6]);if(mag<peak*.04)return;const p=coords(a),dir=new THREE.Vector3(a[4],a[6],a[5]).normalize(),len=12+36*mag/peak,end=p.clone().addScaledVector(dir,len);let side=new THREE.Vector3(0,1,0).cross(dir);if(side.length()<.1)side=new THREE.Vector3(1,0,0).cross(dir);side.normalize();const b=end.clone().addScaledVector(dir,-7);[p,end,end,b.clone().addScaledVector(side,3),end,b.clone().addScaledVector(side,-3)].forEach(v=>{vectors.push(v.x,v.y,v.z);const c=color('--foreground');vcolors.push(c.r,c.g,c.b);});});
  const ag=new THREE.BufferGeometry();ag.setAttribute('position',new THREE.Float32BufferAttribute(vectors,3));ag.setAttribute('color',new THREE.Float32BufferAttribute(vcolors,3));arrows=new THREE.LineSegments(ag,new THREE.LineBasicMaterial({vertexColors:true,transparent:true,opacity:.65}));scene.add(arrows);draw();
}
function inspect(index){if(index<0||index>=data.length)return;const a=data[index];picked=index;$('detail').textContent='('+a[0]+', '+a[1]+', '+a[2]+') µm · '+fmt(a[3])+' Pa · '+fmt(Math.hypot(a[4],a[5],a[6]))+' mm/s';}
function updateCase(){
  const active=D.cases[$('case').value],m=active.meta;
  $('inputs').replaceChildren();m.inputs.forEach((x,i)=>{const s=document.createElement('span');s.textContent='P'+(i+1)+' = '+(x*m.pressure_Pa).toFixed(1)+' Pa';$('inputs').append(s);});
  $('output').replaceChildren();D.matrix.forEach((row,j)=>{const prediction=row.reduce((s,a,i)=>s+a*m.inputs[i]*m.pressure_Pa,0),actual=m.output_mass_equivalent_nLs[j],err=100*(actual-prediction)/(Math.abs(actual)||1);const tr=document.createElement('tr');tr.innerHTML='<th scope="row">Q'+(j+1)+'</th><td class="text-end">'+fmt(prediction)+'</td><td class="text-end">'+fmt(actual)+'</td><td class="text-end">'+(err>=0?'+':'')+err.toFixed(2)+'</td>';$('output').append(tr);});
  $('output-note').textContent=m.scale===2?'nL/s · comparison to the 10 µm matrix':'nL/s · mass-equivalent flow';
  $('stats').textContent=m.fluid_cells.toLocaleString()+' fluid voxels · mass imbalance '+pct(m.mass_balance_relative)+' · Mach '+m.max_mach.toFixed(5);
  $('status').textContent='D3Q19 LBM · no-slip rigid walls · low-Mach, isothermal flow · '+active.samples.length.toLocaleString()+' displayed samples';
  $('detail').textContent='Select a fluid sample to inspect its pressure and velocity.';picked=-1;updateField();
}
try{
  scene=new THREE.Scene();camera=new THREE.PerspectiveCamera(40,1,1,10000);renderer=new THREE.WebGLRenderer({alpha:true,antialias:true});renderer.setPixelRatio(Math.min(window.devicePixelRatio||1,2));$('stage').prepend(renderer.domElement);
  const lines=[],centers=[160,340,520,700];
  function box(center,size){const g=new THREE.EdgesGeometry(new THREE.BoxGeometry(size[0],size[2],size[1]));const p=g.attributes.position.array;const cp=coords(center);for(let i=0;i<p.length;i+=3)lines.push(p[i]+cp.x,p[i+1]+cp.y,p[i+2]+cp.z);g.dispose();}
  centers.forEach(x=>box([x,400,360],[120,800,120]));centers.forEach(y=>box([480,y,100],[800,120,120]));D.widths_um.forEach((row,j)=>row.forEach((w,i)=>{if(w)box([centers[i],centers[j],230],[w,w,260]);}));
  const geom=new THREE.BufferGeometry();geom.setAttribute('position',new THREE.Float32BufferAttribute(lines,3));outline=new THREE.LineSegments(geom,new THREE.LineBasicMaterial({color:color('--foreground'),transparent:true,opacity:.2}));scene.add(outline);
  centers.forEach((x,i)=>{const el=document.createElement('span');el.className='cfd-label text-small';el.textContent='P'+(i+1);$('stage').append(el);labels.push({el,p:coords([x,-30,360]),height:360});});
  centers.forEach((y,i)=>{const el=document.createElement('span');el.className='cfd-label text-small';el.textContent='Q'+(i+1);$('stage').append(el);labels.push({el,p:coords([920,y,100]),height:100});});
  const resize=()=>{const w=$('stage').clientWidth,h=$('stage').clientHeight;renderer.setSize(w,h,false);camera.aspect=w/h;camera.updateProjectionMatrix();cameraUpdate();};new ResizeObserver(resize).observe($('stage'));
  const ray=new THREE.Raycaster();ray.params.Points.threshold=12;let drag=null;
  renderer.domElement.addEventListener('pointerdown',e=>{drag={x:e.clientX,y:e.clientY,lastX:e.clientX,lastY:e.clientY,moved:false};renderer.domElement.setPointerCapture(e.pointerId);});
  renderer.domElement.addEventListener('pointermove',e=>{if(!drag)return;const dx=e.clientX-drag.lastX,dy=e.clientY-drag.lastY;if(Math.hypot(e.clientX-drag.x,e.clientY-drag.y)>4)drag.moved=true;theta-=dx*.008;phi=Math.max(.04,Math.min(Math.PI-.04,phi+dy*.008));drag.lastX=e.clientX;drag.lastY=e.clientY;cameraUpdate();});
  renderer.domElement.addEventListener('pointerup',e=>{if(drag&&!drag.moved&&cloud){const rect=renderer.domElement.getBoundingClientRect();ray.setFromCamera(new THREE.Vector2(2*(e.clientX-rect.left)/rect.width-1,1-2*(e.clientY-rect.top)/rect.height),camera);const hit=ray.intersectObject(cloud);if(hit.length)inspect(hit[0].index);}drag=null;});
  renderer.domElement.addEventListener('pointercancel',()=>drag=null);
  renderer.domElement.addEventListener('wheel',e=>{e.preventDefault();distance=Math.max(850,Math.min(3000,distance*Math.exp(e.deltaY*.001)));cameraUpdate();},{passive:false});
  $('view').addEventListener('change',()=>{const v=$('view').value;if(v==='top'){phi=.04;theta=-Math.PI/2;}else if(v==='side'){phi=Math.PI/2;theta=-Math.PI/2;}else{theta=-.72;phi=1.01;}cameraUpdate();});
  window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change',recolor);new MutationObserver(recolor).observe(document.documentElement,{attributes:true,attributeFilter:['class','style','data-theme']});
  recolor();resize();cameraUpdate();
}catch(err){$('webgl-error').hidden=false;console.error(err);}
$('case').addEventListener('change',updateCase);$('field').addEventListener('change',updateField);$('slice').addEventListener('input',updateField);$('vectors').addEventListener('change',updateField);updateCase();
root.dataset.ready='true';
})();
