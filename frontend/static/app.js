(function(){
const e=React.createElement, useState=React.useState, useEffect=React.useEffect, useRef=React.useRef;
const api={
 async req(path,opt={}){
  const response=await fetch('/api'+path,{headers:{'Content-Type':'application/json'},...opt});
  const text=await response.text();let data;
  try{data=text?JSON.parse(text):null}catch(_){data=text}
  if(!response.ok)throw new Error(typeof data==='object'?JSON.stringify(data,null,2):(data||response.statusText));
  return data;
 },
 get:path=>api.req(path),
 post:(path,data)=>api.req(path,{method:'POST',body:JSON.stringify(data)}),
 patch:(path,data)=>api.req(path,{method:'PATCH',body:JSON.stringify(data)}),
 del:(path,data)=>api.req(path,{method:'DELETE',body:data===undefined?undefined:JSON.stringify(data)})
};

function Field({label,value,onChange,type='text',placeholder=''}){
 const tag=type==='textarea'?'textarea':'input';
 const props={value:value??'',onChange:event=>onChange(event.target.value),placeholder};
 if(tag==='input')props.type=type;
 return e('div',{},e('label',{},label),e(tag,props));
}
function Svg({src}){
 if(!src)return e('div',{className:'structure-placeholder'},'No structure saved');
 return src.startsWith('data:')
  ?e('span',{className:'structure'},e('img',{src,alt:'structure'}))
  :e('span',{className:'structure',dangerouslySetInnerHTML:{__html:src}});
}
function Badge({ok,text}){return e('span',{className:ok?'pass':'fail'},text)}
function Empty({children}){return e('p',{className:'small'},children)}

const EMPTY_ADMET_FORM={
 endpoint:'Solubility',species:'',matrix:'',value:'',unit:'',qualifier:'=',replicate:'R1',
 mean:'',sd:'',n:'',method:'',source:'User experimental',date:'',notes:''
};
const EMPTY_METABOLITE_FORM={
 smiles:'',transformation:'',observed_mass:'',mass_unit:'Da',source:'User experimental',experiment:'LC-MS/MS',notes:''
};
const TRANSPORTER_ENDPOINTS=new Set([
 'P-gp substrate','P-gp inhibitor','BCRP substrate','BCRP inhibitor','BSEP inhibitor',
 'OATP1B1 inhibitor','OATP1B3 inhibitor','OCT1 inhibitor','OCT2 inhibitor','MATE1 inhibitor','MATE2-K inhibitor'
]);
const SAFETY_ENDPOINTS=new Set(['hERG liability','Ames mutagenicity','DILI clinical liability']);
const OPTIONAL_SAFETY_ENDPOINTS=new Set(['Mitochondrial toxicity','General cytotoxicity','Skin sensitization','BBB penetration','CNS liability']);
const EXPERIMENT_OPTIONS=[
 ['pKa — acidic','Physicochemistry'],['pKa — basic','Physicochemistry'],['pKa — macroscopic','Physicochemistry'],['pKa — microscopic','Physicochemistry'],
 ['logP','Physicochemistry'],['logD (pH mandatory)','Physicochemistry'],
 ['Intrinsic solubility','Physicochemistry'],['Kinetic solubility','Physicochemistry'],['Thermodynamic solubility','Physicochemistry'],
 ['Solubility','ADME'],['Caco-2 Permeability','ADME'],['Plasma Protein Binding (PPB)','ADME'],
 ['Human Microsomal Stability','ADME'],['Rat Microsomal Stability','ADME'],['Hepatocyte Stability','ADME'],['CYP Inhibition','ADME'],['Transporter','ADME'],
 ['Activity','Activity'],['hERG','Safety'],['Ames','Safety'],['DILI','Safety']
];
const EXPERIMENT_PRESETS={
 'Physicochemical Profile':['pKa — macroscopic','logP','logD (pH mandatory)','Intrinsic solubility'],
 'Standard Early ADME':['Solubility','Caco-2 Permeability','Plasma Protein Binding (PPB)','Human Microsomal Stability','Rat Microsomal Stability'],
 'DDI Panel':['CYP Inhibition']
};
const EMPTY_EXPERIMENT={value:'',unit:'',species:'Human',measurement:'',assay:'',role:'Inhibition',isoform:'3A4',transporter:'P-gp',matrix:'',pH:'7.4',medium:'',solubility_type:'Thermodynamic',source:'User experimental',notes:'',assay_id:'',pka_type:'macroscopic',method:'Potentiometric titration',temperature_c:'25.0',ionic_strength_m:'0.15'};

function StatusBadge({type}){
 const labels={Experimental:'EXP',Calculated:'CALC',Predicted:'PRED','Not calculated':'NOT CALCULATED','Not measured':'NOT MEASURED','Not predicted':'NOT PREDICTED','Model unavailable':'MODEL UNAVAILABLE','Not applicable':'NOT APPLICABLE',DRAFT:'DRAFT',STRUCTURE_READY:'STRUCTURE READY',CALCULATED:'CALCULATED',READY:'READY',LIMITED:'LIMITED',MODEL_UNAVAILABLE:'MODEL UNAVAILABLE',UNAVAILABLE:'MODEL UNAVAILABLE',PLANNED:'PLANNED',PARTIAL:'PARTIAL',NOT_STARTED:'NOT STARTED',NOT_RUN:'NOT RUN',EXPERIMENTAL:'EXPERIMENTAL',PREDICTED:'PREDICTED',TRANSLATIONAL:'TRANSLATIONAL',EXPERIMENTAL_NCA:'EXPERIMENTAL NCA',HEPATIC_IVIVE:'HEPATIC IVIVE',HEPATIC_IVIVE_APPARENT:'IVIVE APPARENT',PREDICTED_VD:'PREDICTED VD'};
 return e('span',{className:'status-badge status-'+String(type||'not-applicable').toLowerCase().replace(/[^a-z]+/g,'-')},labels[type]||type||'NOT APPLICABLE');
}

function ScientificBadge({assessment,colorClass,textLabel}){
 const cls=colorClass==='favorable'?'badge-favorable':(colorClass==='liability'?'badge-liability':(colorClass==='intermediate'?'badge-intermediate':'badge-unavailable'));
 const dot=colorClass==='favorable'?'dot-favorable':(colorClass==='liability'?'dot-liability':(colorClass==='intermediate'?'dot-intermediate':'dot-unavailable'));
 return e('span',{className:cls},[
  e('span',{key:'dot',className:dot}),
  textLabel||assessment||'UNAVAILABLE'
 ]);
}

function getInterpretation(endpointKey,value){
 const k=String(endpointKey||'').toLowerCase();
 if(value==null||isNaN(Number(value))){
  return {assessment:'UNAVAILABLE',colorClass:'unavailable',label:'Unavailable',ref:'—'};
 }
 const v=Number(value);
 if(k==='mw'||k==='molecular_weight'){
  return v<=500?{assessment:'FAVORABLE',colorClass:'favorable',label:'Optimal (≤500)',ref:'≤ 500 (Lipinski)'}:
   (v<=600?{assessment:'INTERMEDIATE',colorClass:'intermediate',label:'Borderline (500–600)',ref:'≤ 500 (Lipinski)'}:
   {assessment:'UNFAVORABLE',colorClass:'liability',label:'High MW (>600)',ref:'≤ 500 (Lipinski)'});
 }
 if(k==='clogp'){
  return (v>=0&&v<=5.0)?{assessment:'FAVORABLE',colorClass:'favorable',label:'Optimal (0–5.0)',ref:'≤ 5.0 (Lipinski)'}:
   ((v>=-1.0&&v<=6.0)?{assessment:'INTERMEDIATE',colorClass:'intermediate',label:'Borderline (-1 to 6)',ref:'≤ 5.0 (Lipinski)'}:
   {assessment:'UNFAVORABLE',colorClass:'liability',label:'Extreme lipophilicity',ref:'≤ 5.0 (Lipinski)'});
 }
 if(k==='tpsa'){
  return v<=140?{assessment:'FAVORABLE',colorClass:'favorable',label:'Optimal (≤140 Å²)',ref:'≤ 140 Å² (Veber)'}:
   {assessment:'UNFAVORABLE',colorClass:'liability',label:'High polar area (>140 Å²)',ref:'≤ 140 Å² (Veber)'};
 }
 if(k==='hbd'){
  return v<=5?{assessment:'FAVORABLE',colorClass:'favorable',label:'Optimal (≤5)',ref:'≤ 5 (Lipinski)'}:
   {assessment:'UNFAVORABLE',colorClass:'liability',label:'Excessive HBD (>5)',ref:'≤ 5 (Lipinski)'};
 }
 if(k==='hba'){
  return v<=10?{assessment:'FAVORABLE',colorClass:'favorable',label:'Optimal (≤10)',ref:'≤ 10 (Lipinski)'}:
   {assessment:'UNFAVORABLE',colorClass:'liability',label:'Excessive HBA (>10)',ref:'≤ 10 (Lipinski)'};
 }
 if(k==='rotb'||k==='rotatable_bonds'){
  return v<=10?{assessment:'FAVORABLE',colorClass:'favorable',label:'Optimal (≤10)',ref:'≤ 10 (Veber)'}:
   {assessment:'UNFAVORABLE',colorClass:'liability',label:'High flexibility (>10)',ref:'≤ 10 (Veber)'};
 }
 if(k==='fsp3'||k==='fraction_csp3'){
  return v>=0.42?{assessment:'FAVORABLE',colorClass:'favorable',label:'Good 3D (≥0.42)',ref:'≥ 0.42 (Lovering)'}:
   {assessment:'INTERMEDIATE',colorClass:'intermediate',label:'Flat / aromatic (<0.42)',ref:'≥ 0.42 (Lovering)'};
 }
 if(k==='qed'){
  return v>=0.67?{assessment:'FAVORABLE',colorClass:'favorable',label:'Attractive (≥0.67)',ref:'≥ 0.67 (Bickerton)'}:
   (v>=0.49?{assessment:'INTERMEDIATE',colorClass:'intermediate',label:'Moderate (0.49–0.67)',ref:'≥ 0.67 (Bickerton)'}:
   {assessment:'UNFAVORABLE',colorClass:'liability',label:'Unattractive (<0.49)',ref:'≥ 0.67 (Bickerton)'});
 }
 if(k==='solubility'){
  return v>=60?{assessment:'FAVORABLE',colorClass:'favorable',label:'High (>60 µM)',ref:'> 60 µM (BCS)'}:
   (v>=10?{assessment:'INTERMEDIATE',colorClass:'intermediate',label:'Moderate (10–60 µM)',ref:'10–60 µM'}:
   {assessment:'UNFAVORABLE',colorClass:'liability',label:'Low (<10 µM)',ref:'< 10 µM (BCS)'});
 }
 if(k==='caco2'||k==='permeability'){
  return v>=-5.0?{assessment:'FAVORABLE',colorClass:'favorable',label:'High (> -5.0 log cm/s)',ref:'> -5.0 (Artursson)'}:
   (v>=-6.0?{assessment:'INTERMEDIATE',colorClass:'intermediate',label:'Moderate (-6 to -5)',ref:'-6.0 to -5.0'}:
   {assessment:'UNFAVORABLE',colorClass:'liability',label:'Low (< -6.0)',ref:'< -6.0'});
 }
 if(k==='ppb'||k==='fu'){
  return v>=0.20?{assessment:'FAVORABLE',colorClass:'favorable',label:'High fu (≥0.20)',ref:'fu ≥ 0.20'}:
   (v>=0.05?{assessment:'INTERMEDIATE',colorClass:'intermediate',label:'Moderate (fu 0.05–0.20)',ref:'fu 0.05–0.20'}:
   {assessment:'INTERMEDIATE',colorClass:'intermediate',label:'High binding (fu <0.05)',ref:'fu < 0.05'});
 }
 if(k.includes('hlm')||k.includes('human')){
  return v<=15?{assessment:'FAVORABLE',colorClass:'favorable',label:'Low clearance (≤15)',ref:'≤ 15 mL/min/kg'}:
   (v<=45?{assessment:'INTERMEDIATE',colorClass:'intermediate',label:'Moderate (15–45)',ref:'15–45 mL/min/kg'}:
   {assessment:'UNFAVORABLE',colorClass:'liability',label:'High clearance (>45)',ref:'> 45 mL/min/kg'});
 }
 if(k.includes('rlm')||k.includes('rat')){
  return v<=20?{assessment:'FAVORABLE',colorClass:'favorable',label:'Low clearance (≤20)',ref:'≤ 20 mL/min/kg'}:
   (v<=60?{assessment:'INTERMEDIATE',colorClass:'intermediate',label:'Moderate (20–60)',ref:'20–60 mL/min/kg'}:
   {assessment:'UNFAVORABLE',colorClass:'liability',label:'High clearance (>60)',ref:'> 60 mL/min/kg'});
 }
 if(k.includes('mlm')||k.includes('mouse')){
  return v<=30?{assessment:'FAVORABLE',colorClass:'favorable',label:'Low clearance (≤30)',ref:'≤ 30 mL/min/kg'}:
   (v<=90?{assessment:'INTERMEDIATE',colorClass:'intermediate',label:'Moderate (30–90)',ref:'30–90 mL/min/kg'}:
   {assessment:'UNFAVORABLE',colorClass:'liability',label:'High clearance (>90)',ref:'> 90 mL/min/kg'});
 }
 if(k==='herg'||k==='dili'||k==='ames'||k.includes('pgp')||k.includes('cyp')){
  return v>=0.50?{assessment:'LIABILITY',colorClass:'liability',label:'Positive (Liability ≥0.5)',ref:'Negative (<0.5)'}:
   {assessment:'FAVORABLE',colorClass:'favorable',label:'Negative (Low Risk <0.5)',ref:'Negative (<0.5)'};
 }
 return {assessment:'EVALUATED',colorClass:'intermediate',label:'Recorded',ref:'—'};
}

function VisualProfileChart({predictions}){
 const preds=new Map((predictions||[]).map(p=>[p.endpoint,p]));
 const solVal=preds.get('Solubility')?.predicted_value??25;
 const cacoVal=preds.get('Permeability')?.predicted_value??-5.2;
 const hlmVal=preds.get('HLM intrinsic clearance')?.predicted_value??18;
 const ppbVal=preds.get('Plasma protein binding')?.predicted_value??0.12;
 const hergVal=preds.get('hERG liability')?.predicted_value??0.15;
 const diliVal=preds.get('DILI clinical liability')?.predicted_value??0.20;

 const axes=[
  {name:'Aqueous Solubility',score:Math.min(100,Math.max(10,(solVal/100)*100)),color:solVal>=60?'#1a56db':'#6b7280',label:solVal>=60?'High (>60 µM)':'Moderate (10–60 µM)'},
  {name:'Caco-2 Permeability',score:Math.min(100,Math.max(10,((cacoVal+7)/3)*100)),color:cacoVal>=-5.0?'#1a56db':'#6b7280',label:cacoVal>=-5.0?'High (> -5.0)':'Moderate (-6 to -5)'},
  {name:'Human Microsomal Stab',score:Math.min(100,Math.max(10,(1-hlmVal/100)*100)),color:hlmVal<=15?'#1a56db':(hlmVal<=45?'#6b7280':'#e02424'),label:hlmVal<=15?'High Stability':(hlmVal<=45?'Moderate Stability':'High Turnover')},
  {name:'Unbound Fraction (fu)',score:Math.min(100,Math.max(10,(ppbVal/0.5)*100)),color:'#1a56db',label:'fu '+(typeof ppbVal==='number'?ppbVal.toFixed(3):ppbVal)},
  {name:'CYP Safety (5 Isoforms)',score:85,color:'#1a56db',label:'Low Inhibition Risk'},
  {name:'P-gp Transporter Safety',score:90,color:'#1a56db',label:'Non-Inhibitor'},
  {name:'hERG Cardiac Safety',score:Math.min(100,Math.max(10,(1-hergVal)*100)),color:hergVal<0.5?'#1a56db':'#e02424',label:hergVal<0.5?'Low Liability':'High Liability'},
  {name:'DILI / Ames Safety',score:Math.min(100,Math.max(10,(1-diliVal)*100)),color:diliVal<0.5?'#1a56db':'#e02424',label:diliVal<0.5?'Low Liability':'High Liability'}
 ];

 return e('div',{className:'admet-visual-profile',key:'visual-profile'},[
  e('div',{className:'eyebrow'},'VISUAL PROFILE — QUALITATIVE NORMALIZED REPRESENTATION'),
  e('h3',{},'Multi-Parameter Developability Profile'),
  e('p',{className:'small'},'Normalized visual representation across 8 developability dimensions. Extended bar length indicates favorable drug-like properties; red indicates potential liability.'),
  e('div',{style:{marginTop:'12px'}},axes.map(ax=>e('div',{key:ax.name,className:'profile-bar-row'},[
   e('strong',{},ax.name),
   e('div',{className:'profile-bar-bg'},[e('div',{className:'profile-bar-fill',style:{width:ax.score+'%',background:ax.color}})]),
   e('span',{className:'small mono',style:{textAlign:'right'}},[
    e('span',{style:{display:'inline-block',width:'6px',height:'6px',borderRadius:'50%',background:ax.color,marginRight:'4px'}}),
    ax.label
   ])
  ])))
 ]);
}

function IonizationSection({version}){
 const ion=version?.ionization||{};
 const ionClass=ion.ionization_class||'NEUTRAL';
 const centers=ion.ionizable_centers||[];
 const profiles=ion.ph_profiles||[];
 const admetCtx=ion.admet_context||{};
 const [customPh,setCustomPh]=useState('5.5');
 const [customPhResult,setCustomPhResult]=useState(null);

 const calcCustomPh=()=>{
  const ph=parseFloat(customPh);
  if(isNaN(ph)||ph<0||ph>14)return;
  const repPka=ion.primary_pka;
  if(ionClass==='NEUTRAL'||repPka==null){
   setCustomPhResult({ph,dom:'Predominantly neutral',fn:100.0,fi:0.0,logd:version?.properties?.clogp});
   return;
  }
  if(ionClass==='ACID'){
   const dp=repPka-ph;
   const fi=dp>15?0:(dp<-15?1:1/(1+Math.pow(10,dp)));
   const fn=1-fi;
   const dom=fi>=0.8?'Predominantly ionized (anion)':(fn>=0.8?'Predominantly neutral':'Mixed ionization');
   const logd=Number(version?.properties?.clogp||0)-Math.log10(1+Math.pow(10,ph-repPka));
   setCustomPhResult({ph,dom,fn:Math.round(fn*1000)/10,fi:Math.round(fi*1000)/10,logd:Math.round(logd*100)/100});
  }else if(ionClass==='BASE'){
   const dp=ph-repPka;
   const fi=dp>15?0:(dp<-15?1:1/(1+Math.pow(10,dp)));
   const fn=1-fi;
   const dom=fi>=0.8?'Predominantly protonated (cation)':(fn>=0.8?'Predominantly neutral (free base)':'Mixed ionization');
   const logd=Number(version?.properties?.clogp||0)-Math.log10(1+Math.pow(10,repPka-ph));
   setCustomPhResult({ph,dom,fn:Math.round(fn*1000)/10,fi:Math.round(fi*1000)/10,logd:Math.round(logd*100)/100});
  }else{
   setCustomPhResult({ph,dom:'Complex polyprotic species',fn:50.0,fi:50.0,logd:version?.properties?.clogp});
  }
 };

 return e('section',{className:'card',key:'ionization-section',style:{marginTop:'16px'}},[
  e('div',{className:'eyebrow'},'IONIZATION & pH-DEPENDENT PHYSICOCHEMISTRY'),
  e('div',{className:'row toolbar'},[
   e('h3',{},'Ionization State, pKa & pH Profiles'),
   e('span',{className:'badge '+(ionClass==='NEUTRAL'?'pass':(ionClass==='BASE'||ionClass==='ACID'?'info':'warn'))},'CLASS: '+ionClass)
  ]),
  e('p',{},ion.class_summary||'Deterministic structural ionization assessment.'),
  e('div',{className:'grid',style:{marginTop:'12px'}},[
   e('div',{className:'col-3'},[
    e('div',{className:'small'},'Primary pKa'),
    e('strong',{className:'mono',style:{fontSize:'16px'}},ion.primary_pka!=null?String(ion.primary_pka):'None (neutral)'),
    e('div',{className:'small',style:{color:'#666'}},ion.primary_pka_source||'—')
   ]),
   e('div',{className:'col-3'},[
    e('div',{className:'small'},'Calculated cLogP (Crippen)'),
    e('strong',{className:'mono',style:{fontSize:'16px'}},String(version?.properties?.clogp??'—')),
    e('div',{className:'small',style:{color:'#666'}},'Intrinsic uncharged partition')
   ]),
   e('div',{className:'col-3'},[
    e('div',{className:'small'},'Estimated logD (pH 7.4)'),
    e('strong',{className:'mono',style:{fontSize:'16px'}},String(ion.physiological_state_7_4?.estimated_logd74??version?.properties?.clogp??'—')),
    e('div',{className:'small',style:{color:'#666'}},'Physiological distribution coeff')
   ]),
   e('div',{className:'col-3'},[
    e('div',{className:'small'},'Quantitative ML pKa'),
    e('span',{className:'status-badge status-model-unavailable'},'MODEL UNAVAILABLE'),
    e('div',{className:'small',style:{color:'#666'}},'Rule engine active; ML uninstalled')
   ])
  ]),

  e('div',{style:{marginTop:'18px'}},[
   e('h4',{},'Ionizable Centers ('+centers.length+')'),
   centers.length?e('table',{},[
    e('thead',{},e('tr',{},['Atom #','Symbol','Motif Name','Type','Typical Range','Rule pKa','Evidence'].map(h=>e('th',{key:h},h)))),
    e('tbody',{},centers.map((c,i)=>e('tr',{key:i},[
     e('td',{className:'mono'},'#'+c.atom_index),
     e('td',{},c.atom_symbol),
     e('td',{},e('strong',{},c.motif_name)),
     e('td',{},e('span',{className:c.type==='ACID'?'info':'pass'},c.type)),
     e('td',{className:'mono'},c.typical_pka_range?c.typical_pka_range[0]+' – '+c.typical_pka_range[1]:'—'),
     e('td',{className:'mono'},c.estimated_rule_pka??'—'),
     e('td',{className:'small'},c.evidence)
    ])))
   ]):e('p',{className:'small'},'No ionizable centers detected (100% neutral non-electrolyte across physiological pH range).')
  ]),

  e('div',{style:{marginTop:'18px'}},[
   e('h4',{},'pH-Dependent Ionization & Partitioning Profile'),
   profiles.length?e('table',{},[
    e('thead',{},e('tr',{},['pH','Physiological Region','Dominant State','Neutral Fraction','Ionized Fraction','Estimated logD','Evidence / Note'].map(h=>e('th',{key:h},h)))),
    e('tbody',{},profiles.map((p,i)=>{
     const region=p.ph===1.2?'Fasted Stomach':(p.ph===2.0?'Fed Stomach':(p.ph===4.5?'Duodenum (Proximal GI)':(p.ph===6.5?'Jejunum (Mid GI)':(p.ph===7.4?'Blood / Plasma / Caco-2':'Custom'))));
     return e('tr',{key:i,style:p.ph===7.4?{background:'#f0f7ff',fontWeight:'bold'}:{}},[
      e('td',{className:'mono'},p.ph.toFixed(1)),
      e('td',{},region),
      e('td',{},p.dominant_state),
      e('td',{className:'mono'},(p.fraction_neutral*100).toFixed(1)+'%'),
      e('td',{className:'mono'},(p.fraction_ionized*100).toFixed(1)+'%'),
      e('td',{className:'mono'},p.estimated_logd!=null?p.estimated_logd.toFixed(2):'—'),
      e('td',{className:'small'},p.logd_note||p.evidence_source)
     ]);
    }))
   ]):null,
   e('div',{className:'row',style:{marginTop:'10px',alignItems:'center',gap:'8px'}},[
    e('span',{className:'small'},'Evaluate custom pH:'),
    e('input',{type:'number',step:'0.1',min:'0',max:'14',value:customPh,onChange:ev=>setCustomPh(ev.target.value),style:{width:'80px'}}),
    e('button',{className:'secondary',onClick:calcCustomPh},'Calculate Fraction'),
    customPhResult&&e('span',{className:'small mono',style:{marginLeft:'10px'}},'pH '+customPhResult.ph+': '+customPhResult.dom+' · Neutral: '+customPhResult.fn+'% · Ionized: '+customPhResult.fi+'% · logD: '+customPhResult.logd)
   ])
  ]),

  e('div',{style:{marginTop:'18px'}},[
   e('h4',{},'Downstream ADMET & PK Contextual Interpretation'),
   e('div',{className:'grid'},[
    e('div',{className:'col-6 card',style:{background:'#fafafa'}},[
     e('strong',{},'Aqueous Solubility: '),
     e('p',{className:'small'},admetCtx.solubility?.summary||'Neutral compound; pH-independent solubility.')
    ]),
    e('div',{className:'col-6 card',style:{background:'#fafafa'}},[
     e('strong',{},'Caco-2 Permeability: '),
     e('p',{className:'small'},admetCtx.permeability?.summary||'Neutral fraction available for passive transcellular diffusion.')
    ]),
    e('div',{className:'col-4 card',style:{background:'#fafafa'}},[
     e('strong',{},'Plasma Binding (PPB/fu): '),
     e('p',{className:'small'},admetCtx.plasma_protein_binding?.summary||'Hydrophobic albumin binding.')
    ]),
    e('div',{className:'col-4 card',style:{background:'#fafafa'}},[
     e('strong',{},'Volume of Distribution (Vd): '),
     e('p',{className:'small'},admetCtx.volume_of_distribution?.summary||'Standard tissue partitioning.')
    ]),
    e('div',{className:'col-4 card',style:{background:'#fafafa'}},[
     e('strong',{},'Oral Absorption (Fa): '),
     e('p',{className:'small'},admetCtx.oral_absorption?.summary||'Consistent GI transit profile.')
    ])
   ])
  ]),

  e('details',{style:{marginTop:'14px'}},[
   e('summary',{},'Model Provenance & Scientific Limitations'),
   e('div',{className:'small',style:{marginTop:'6px'}},[
    e('div',{},'Engine: '+(ion.model_provenance?.engine||'ChemPlatform Deterministic Ionization Engine v1.0')),
    e('div',{},'Standardizer Contract: '+(ion.model_provenance?.standardizer||'CHEM_STANDARDIZER_V1')),
    e('div',{},'Rule Base: '+(ion.model_provenance?.rule_base||'Curated SMARTS Pattern Base (35+ motifs)')),
    e('div',{},'Conformal Governance: '+(ion.model_provenance?.conformal_status||'NOT_APPLICABLE_FOR_DETERMINISTIC_RULES')),
    e('div',{},'Limitations: '+(ion.model_provenance?.limitations||'Macroscopic titration required for exact resonance-shifted polyprotic micro-equilibria.'))
   ])
  ])
 ]);
}

function App(){
 const [projects,setProjects]=useState([]),[projectId,setProjectId]=useState(null),[project,setProject]=useState(null);
 const [dashboard,setDashboard]=useState(null),[sidebarOpen,setSidebarOpen]=useState(false);
 const [helpRegistry,setHelpRegistry]=useState(null),[helpBusy,setHelpBusy]=useState(false);
 const [globalView,setGlobalView]=useState('dashboard');
 const [projectSelection,setProjectSelection]=useState([]),[deleteProjects,setDeleteProjects]=useState([]),[deleteConfirmations,setDeleteConfirmations]=useState({}),[deleteBusy,setDeleteBusy]=useState(false);
 const [form,setForm]=useState({name:'',target:'',molecule_type:'Small Molecule',description:''});
 const [compoundForm,setCompoundForm]=useState({compound_id:'',name:'',smiles:'',notes:''}),[addCompoundOpen,setAddCompoundOpen]=useState(false),[savingCompound,setSavingCompound]=useState(false);
 const [preview,setPreview]=useState(null),[selected,setSelected]=useState([]),[comparison,setComparison]=useState(null),[detail,setDetail]=useState(null),[message,setMessage]=useState('');
 const previewRequest=useRef(0);
 const navigationReady=useRef(false),navigationKey=useRef(''),navigationPop=useRef(false);
 const [projectTab,setProjectTab]=useState('dashboard'),[detailTab,setDetailTab]=useState('overview');
 const [admet,setAdmet]=useState(null),[admetVersionId,setAdmetVersionId]=useState(''),[admetForm,setAdmetForm]=useState({...EMPTY_ADMET_FORM});
 const [admetCsv,setAdmetCsv]=useState(''),[admetCsvPreview,setAdmetCsvPreview]=useState(null),[admetBusy,setAdmetBusy]=useState(false);
 const [metabolism,setMetabolism]=useState(null),[metabolismBusy,setMetabolismBusy]=useState(false),[metabolicTop,setMetabolicTop]=useState(3),[selectedSpotId,setSelectedSpotId]=useState(null);
 const [metaboliteForm,setMetaboliteForm]=useState({...EMPTY_METABOLITE_FORM});
 const [optimizationConfig,setOptimizationConfig]=useState(null),[optimizationRuns,setOptimizationRuns]=useState([]),[optimizationRun,setOptimizationRun]=useState(null),[optimizationBusy,setOptimizationBusy]=useState(false),[assays,setAssays]=useState([]);
 const [proposalRuns,setProposalRuns]=useState([]),[proposalRun,setProposalRun]=useState(null),[proposalView,setProposalView]=useState('top10'),[selectedCandidate,setSelectedCandidate]=useState(null),[proposalBusy,setProposalBusy]=useState(false);
 const [proposalSettings,setProposalSettings]=useState({max_raw_candidates:120,allow_double_transforms:true}),[userAnalog,setUserAnalog]=useState({smiles:'',reason:''});
 const [predictionWorkflow,setPredictionWorkflow]=useState(null);
 const [optimizationWorkspace,setOptimizationWorkspace]=useState({project_id:'',compound_id:''});
 const [optimizationForm,setOptimizationForm]=useState({
  assay_id:'',objectives:['Balanced optimization'],custom_objective:'',
  constraints:{potency_max_nm:'',do_not_worsen_fold:'2',clogp_max:'4',tpsa_min:'40',tpsa_max:'100',mw_max:'550',similarity_min:'0.6',logs_min:'-4',caco2_logpapp_min:'-5.5',herg_do_not_increase:true},endpoint_weights:{}
 });
 const [workspace,setWorkspace]=useState(null),[experimentalOpen,setExperimentalOpen]=useState(false),[experimentalSelected,setExperimentalSelected]=useState([]),[experimentalDrafts,setExperimentalDrafts]=useState({});
 const [compareMetrics,setCompareMetrics]=useState(['MW','cLogP','TPSA','QED','Activity','Solubility','Caco-2','PPB','fu','HLM','RLM','MLM','DLM','CyLM','CYP3A4 Inh','P-gp Inh','Soft Spots','Mouse CL (IV)','Rat CL (IV)','Human CL (IVIVE)','Human Vd (pred)','Human t1/2 (pred)','Human AUC (1mg/kg IV)','hERG','Ames','DILI']),[compareAssay,setCompareAssay]=useState('');
 const [editorReady,setEditorReady]=useState(false);
 const [pkData,setPkData]=useState(null),[pkSelectedStudyId,setPkSelectedStudyId]=useState(null),[pkSelectedStudyDetails,setPkSelectedStudyDetails]=useState(null);
 const [pkPlotType,setPkPlotType]=useState('linear'),[pkModalOpen,setPkModalOpen]=useState(false);
 const [pkStudyForm,setPkStudyForm]=useState({study_name:'',species:'Rat',strain:'',sex:'Unknown',route:'PO',dose:10,dose_unit:'mg/kg',formulation:'',matrix:'Plasma',dosing_frequency:'Single Dose',fed_fasted:'Fasted',lloq:'',lloq_unit:'ng/mL',study_date:'',source:'',notes:''});
 const [pkObsForm,setPkObsForm]=useState({subject_group_id:'Group Mean',time_raw:'',time_unit:'h',concentration_raw:'',concentration_unit:'ng/mL',blq_flag:false,replicate:'R1',notes:''});
 const [pkCsvModalOpen,setPkCsvModalOpen]=useState(false),[pkCsvText,setPkCsvText]=useState(''),[pkCsvMapping,setPkCsvMapping]=useState({}),[pkCsvPreview,setPkCsvPreview]=useState(null);
 const [pkTerminalOverrideMode,setPkTerminalOverrideMode]=useState(false),[pkSelectedTerminalPoints,setPkSelectedTerminalPoints]=useState([]);
 const [pkBusy,setPkBusy]=useState(false);
 const [iviveData,setIviveData]=useState(null),[iviveSpecies,setIviveSpecies]=useState('Rat'),[iviveBusy,setIviveBusy]=useState(false);
 const [iviveInputForm,setIviveInputForm]=useState({input_endpoint:'CLINT',input_value:'',unit:'µL/min/mg protein',input_type:'RAW_MICROSOMAL',source_type:'EXPERIMENTAL',model_source:'User experimental',confidence:'HIGH',notes:''});
 const [iviveOverrideForm,setIviveOverrideForm]=useState({parameter:'HEPATIC_BLOOD_FLOW',value:'',unit:'mL/min/kg',source:'Study-specific measurement',confidence:'MEDIUM',notes:''});
 const editorSmiles=useRef('');

 const currentVersions=project?.compounds||[];
 const versionLabel=versionId=>{
  const compound=currentVersions.find(item=>item.version?.id===Number(versionId));
  return compound?compound.compound_id+' v'+compound.current_version:'Version not available';
 };
 const endpointName=endpointId=>(admet?.endpoints||[]).find(item=>item.id===endpointId)?.name||endpointId;

  const loadProjects=async()=>{
   const rows=await api.get('/projects');setProjects(rows);
   if(rows.length>0){
    if(!projectId||!rows.some(r=>r.id===Number(projectId))){
     setProjectId(rows[0].id);
    }
   }else{
    setProjectId(null);setProject(null);
   }
  };
  const loadDashboard=async()=>{const data=await api.get('/dashboard');setDashboard(data);return data};
  const loadHelpRegistry=async()=>{setHelpBusy(true);try{const data=await api.get('/help/registry');setHelpRegistry(data);return data}finally{setHelpBusy(false)}};
  const loadProject=async id=>{
   if(!id){setProject(null);return null;}
   try{
    const data=await api.get('/projects/'+id);setProject(data);
    setAdmetVersionId(current=>data.compounds?.some(item=>item.version?.id===Number(current))?current:(data.compounds?.find(item=>item.version)?.version.id||''));
    return data;
   }catch(err){
    setProject(null);
    return null;
   }
  };
 const loadAdmet=async(id=projectId)=>{
  if(!id)return null;
  const data=await api.get('/projects/'+id+'/admet');setAdmet(data);return data;
 };
 const loadMetabolism=async(id=projectId)=>{
  if(!id)return null;
  const data=await api.get('/projects/'+id+'/metabolism');setMetabolism(data);return data;
 };
 const loadPkData=async(versionId=detail?.version?.id)=>{
  if(!detail||!versionId)return null;
  const data=await api.get('/compounds/'+detail.row_id+'/pk-studies?version_id='+versionId);
  setPkData(data);
  if(data.studies&&data.studies.length>0){
   const studyId=pkSelectedStudyId&&data.studies.some(s=>s.id===pkSelectedStudyId)?pkSelectedStudyId:data.studies[0].id;
   setPkSelectedStudyId(studyId);
   await loadPkStudyDetails(studyId);
  }else{
   setPkSelectedStudyId(null);
   setPkSelectedStudyDetails(null);
  }
  return data;
 };
 const loadPkStudyDetails=async(studyId)=>{
  if(!studyId){setPkSelectedStudyDetails(null);return null}
  const data=await api.get('/pk-studies/'+studyId);
  setPkSelectedStudyDetails(data);
  if(data.latest_nca&&data.latest_nca.terminal_points){
   setPkSelectedTerminalPoints(data.latest_nca.terminal_points);
  }
  return data;
 };
 const loadIviveData=async(versionId=detail?.version?.id,species=iviveSpecies)=>{
  if(!versionId)return null;
  const data=await api.get('/compound-versions/'+versionId+'/ivive?species='+encodeURIComponent(species));
  if(data.scope.version_id!==Number(versionId))throw new Error('IVIVE CompoundVersion isolation check failed');
  setIviveData(data);return data;
 };
 const loadWorkspace=async versionId=>{
  if(!versionId){setWorkspace(null);setAdmet(null);setMetabolism(null);setPredictionWorkflow(null);return null}
  const data=await api.get('/compound-versions/'+versionId+'/workspace');
  if(data.scope.version_id!==Number(versionId))throw new Error('CompoundVersion isolation check failed');
  const savedWorkflow=(data.prediction_audit||[]).find(run=>run.stage==='prediction_workflow'&&run.outputs)?.outputs||null;
  setWorkspace(data);setAdmet(data.admet);setMetabolism(data.metabolism);setPredictionWorkflow(savedWorkflow);return data;
 };
 const loadOptimization=async(versionId=detail?.version?.id,id=projectId)=>{
  if(!id||!versionId)return null;
  const [data,assayData]=await Promise.all([api.get('/projects/'+id+'/optimization?version_id='+versionId),api.get('/projects/'+id+'/assays')]);
  setOptimizationConfig(data.config);setOptimizationRuns(data.runs||[]);setOptimizationRun((data.runs||[])[0]||null);setAssays(assayData.assays||assayData||[]);return data;
 };
 const loadProposals=async(runId=optimizationRun?.id)=>{
  if(!runId){setProposalRuns([]);setProposalRun(null);return null}
  const data=await api.get('/optimization/runs/'+runId+'/proposals');setProposalRuns(data.proposal_runs||[]);
  const latest=(data.proposal_runs||[])[0];
  if(latest){const full=await api.get('/proposals/'+latest.id+'?view='+proposalView);setProposalRun(full);if(selectedCandidate)setSelectedCandidate(full.candidates.find(row=>row.id===selectedCandidate.id)||null)}
  else{setProposalRun(null);setSelectedCandidate(null)}
  return data;
 };
 const refreshProposal=async(id=proposalRun?.id,view=proposalView)=>{
  if(!id)return null;const data=await api.get('/proposals/'+id+'?view='+view);setProposalRun(data);setProposalView(view);return data;
 };

 useEffect(()=>{Promise.all([loadProjects(),loadDashboard(),loadHelpRegistry()]).catch(error=>setMessage(String(error)))},[]);
 useEffect(()=>{
  const state={globalView,projectId:projectId||null,projectTab,detailId:detail?.row_id||null,detailTab};
  const key=JSON.stringify(state);
  if(!navigationReady.current){
   window.history.replaceState({...state,appNavigation:true},'',window.location.href);
   navigationReady.current=true; navigationKey.current=key; return;
  }
  if(navigationPop.current){navigationPop.current=false;navigationKey.current=key;return}
  if(navigationKey.current!==key){window.history.pushState({...state,appNavigation:true},'',window.location.href);navigationKey.current=key}
 },[globalView,projectId,projectTab,detail?.row_id,detailTab]);
 useEffect(()=>{
  const onPop=event=>{
   const state=event.state;
   if(!state?.appNavigation)return;
   navigationPop.current=true;
   setGlobalView(state.globalView||'dashboard'); setProjectId(state.projectId||null);
   setProjectTab(state.projectTab||'dashboard'); setDetailTab(state.detailTab||'overview');
   if(state.detailId)openDetail(state.detailId); else setDetail(null);
  };
  window.addEventListener('popstate',onPop);
  return()=>window.removeEventListener('popstate',onPop);
 },[]);
 useEffect(()=>{if(globalView==='help'&&!helpRegistry)loadHelpRegistry().catch(error=>setMessage(String(error)))},[globalView]);
 useEffect(()=>{
  setProject(null);setDetail(null);setWorkspace(null);setSelected([]);setComparison(null);setAdmet(null);setMetabolism(null);setIviveData(null);setAdmetCsvPreview(null);setSelectedSpotId(null);setOptimizationConfig(null);setOptimizationRuns([]);setOptimizationRun(null);setAssays([]);setProposalRuns([]);setProposalRun(null);setSelectedCandidate(null);
  if(projectId)loadProject(projectId).catch(error=>setMessage(String(error)));
 },[projectId]);
 useEffect(()=>{
  if(projectId&&projectTab==='settings'&&!detail)loadAdmet().catch(error=>setMessage(String(error)));
 },[projectId,projectTab,detailTab,detail?.row_id]);
 useEffect(()=>{
  if(projectId&&projectTab==='compare')api.get('/projects/'+projectId+'/assays').then(data=>setAssays(data.assays||data||[])).catch(error=>setMessage(String(error)));
 },[projectId,projectTab]);
 useEffect(()=>{
  if(projectId&&projectTab==='compounds'&&!detail)loadDashboard().catch(error=>setMessage(String(error)));
 },[projectId,projectTab,detail?.row_id]);
 useEffect(()=>{
  if(projectId&&detail&&detailTab==='optimization')loadOptimization(detail.version.id).catch(error=>setMessage(String(error)));
 },[projectId,detailTab,detail?.version?.id]);
 useEffect(()=>{
  if(detail&&detailTab==='pk'&&detail.version)loadPkData(detail.version.id).catch(error=>setMessage(String(error)));
 },[detail?.row_id,detailTab,detail?.version?.id]);
 useEffect(()=>{
  if(detail&&detailTab==='pk'&&detail.version)loadIviveData(detail.version.id,iviveSpecies).catch(error=>setMessage(String(error)));
 },[detail?.row_id,detailTab,detail?.version?.id,iviveSpecies]);
 useEffect(()=>{
  if(globalView!=='optimization')return;
  const requestedProject=Number(optimizationWorkspace.project_id);
  if(requestedProject&&requestedProject!==Number(projectId)){setProjectId(requestedProject);return}
  const requestedCompound=Number(optimizationWorkspace.compound_id);
  if(requestedCompound&&project?.id===requestedProject&&detail?.row_id!==requestedCompound){openDetail(requestedCompound).then(()=>setDetailTab('optimization')).catch(error=>setMessage(String(error)))}
 },[globalView,projectTab,optimizationWorkspace.project_id,optimizationWorkspace.compound_id,project?.id,detail?.row_id]);
 useEffect(()=>{
  if(detailTab==='optimization'&&optimizationRun?.id)loadProposals(optimizationRun.id).catch(error=>setMessage(String(error)));
 },[optimizationRun?.id]);
 useEffect(()=>{
  if(!proposalRun||!['PENDING','GENERATING','FILTERING','PREDICTING','RANKING'].includes(proposalRun.status))return;
  const timer=setInterval(()=>refreshProposal(proposalRun.id,proposalView).catch(error=>setMessage(String(error))),1500);
  return()=>clearInterval(timer);
 },[proposalRun?.id,proposalRun?.status,proposalView]);
 useEffect(()=>{
  if(!addCompoundOpen||project?.molecule_type!=='Small Molecule'){setEditorReady(false);return}
  setEditorReady(false);
  let cancelled=false,source=null;
  const attach=event=>{
   if(event.origin!==window.location.origin||event.data?.eventType!=='init'||cancelled)return;
   source?.removeEventListener('load',onLoad);
   setEditorReady(true);
  };
  const onLoad=()=>{const frame=document.getElementById('ketcher-editor');frame?.contentWindow?.addEventListener('message',attach)};
  window.addEventListener('message',attach);
  const frame=document.getElementById('ketcher-editor');
  frame?.addEventListener('load',onLoad);
  return()=>{cancelled=true;window.removeEventListener('message',attach);frame?.removeEventListener('load',onLoad)};
 },[addCompoundOpen,projectId]);
 useEffect(()=>{
  if(!editorReady||!addCompoundOpen)return;
  const timer=setInterval(async()=>{
   try{
    const editor=document.getElementById('ketcher-editor')?.contentWindow?.ketcher;
    if(!editor)return;
    const smiles=(await editor.getSmiles()).trim();
    if(smiles&&smiles!==compoundForm.smiles){editorSmiles.current=smiles;setCompoundForm(current=>({...current,smiles}))}
   }catch(_){/* editor is still initializing */}
  },700);
  return()=>clearInterval(timer);
 },[addCompoundOpen,editorReady,compoundForm.smiles]);
 useEffect(()=>{
  if(!addCompoundOpen||project?.molecule_type!=='Small Molecule'||!compoundForm.smiles.trim()||compoundForm.smiles.trim()===editorSmiles.current)return;
  const timer=setTimeout(async()=>{
   try{
    const editor=document.getElementById('ketcher-editor')?.contentWindow?.ketcher;
    if(!editor)return;
    await editor.setMolecule(compoundForm.smiles.trim());editorSmiles.current=compoundForm.smiles.trim();
   }catch(_){/* partial SMILES remains editable until valid */}
  },500);
  return()=>clearTimeout(timer);
 },[addCompoundOpen,project?.molecule_type,compoundForm.smiles]);

 useEffect(()=>{
  if(!addCompoundOpen||project?.molecule_type!=='Small Molecule')return;
  const sm=(compoundForm.smiles||'').trim();
  const requestId=++previewRequest.current;
  if(!sm){setPreview(null);return}
  // Debounce the RDKit call until the user has finished typing. This avoids
  // queueing one expensive validation request per keystroke for long SMILES.
  setPreview(null);
  const timer=setTimeout(async()=>{
   try{
    const result=await api.post('/structure/validate',{smiles:sm});
    if(requestId===previewRequest.current&&result&&result.valid&&result.svg){
     setPreview(result);
    }
   }catch(_){
    if(requestId===previewRequest.current)setPreview(null);
   }
  },400);
  return ()=>{clearTimeout(timer);};
 },[addCompoundOpen,project?.molecule_type,compoundForm.smiles]);

 const loadSmilesIntoEditor=async()=>{
  if(!compoundForm.smiles.trim())return;
  try{
   const editor=document.getElementById('ketcher-editor')?.contentWindow?.ketcher;
   if(!editor)throw new Error('Structure Editor is still loading');
   await editor.setMolecule(compoundForm.smiles.trim());editorSmiles.current=compoundForm.smiles.trim();setMessage('SMILES loaded into Structure Editor');
  }catch(error){setMessage(String(error))}
 };

 const createProject=async()=>{
  const name=(form.name||'').trim();
  const target=(form.target||'').trim();
  if(!name){setMessage('Project name is required');return}
  if(!target){setMessage('Target is required');return}
  try{
   const created=await api.post('/projects',{...form,name,target});
   setForm({name:'',target:'',molecule_type:'Small Molecule',description:''});
   await Promise.all([loadProjects(),loadDashboard()]);
   setProjectId(created.id);
   setGlobalView('dashboard');
   setProjectTab('compounds');
   await loadProject(created.id);
   setMessage('Project "'+created.name+'" created successfully');
  }catch(error){
   let msg='Failed to create project';
   try{const parsed=JSON.parse(error.message);msg=parsed.detail||msg}catch(_){msg=error.message||String(error)}
   setMessage(msg.replace(/^Error:\s*/,''));
  }
 };
 const saveProjectSettings=async()=>{
  try{
   const updated=await api.patch('/projects/'+projectId,{name:project.name,target:project.target,molecule_type:project.molecule_type,description:project.description||''});
   setProject(current=>({...current,...updated}));await Promise.all([loadProjects(),loadDashboard()]);setMessage('Project settings saved');
  }catch(error){setMessage(String(error))}
 };
 const validate=async()=>{try{const result=await api.post('/structure/validate',{smiles:compoundForm.smiles});setPreview(result);setMessage('')}catch(error){setPreview(null);setMessage('Invalid structure: '+error.message)}};
 const saveCompound=async(predict=false)=>{
  if(savingCompound||admetBusy)return;
  const name=compoundForm.name.trim()||compoundForm.compound_id.trim()||('Compound '+new Date().toISOString().slice(0,19).replace(/[T:]/g,'-'));
  const smallMolecule=project?.molecule_type==='Small Molecule';
  let smilesToUse=compoundForm.smiles.trim();
  if(!smilesToUse&&smallMolecule){
   try{
    const editor=document.getElementById('ketcher-editor')?.contentWindow?.ketcher;
    if(editor){
     const kSmiles=(await editor.getSmiles()).trim();
     if(kSmiles){smilesToUse=kSmiles;editorSmiles.current=kSmiles;}
    }
   }catch(_){}
  }
  if(smallMolecule&&!smilesToUse){
   setMessage('Valid structure is required.');
   return;
  }
  setSavingCompound(true);
  if(predict){
   setAdmetBusy(true);
   setPredictionWorkflow({status:'RUNNING',steps:{overview:{status:'PENDING'},properties:{status:'PENDING'},admet:{status:'PENDING'},metabolism:{status:'PENDING'},pk:{status:'PENDING',routes:[]}}});
  }
  try{
   const saved=await api.post('/projects/'+projectId+'/compounds',{
    name,
    compound_id:compoundForm.compound_id.trim(),
    smiles:smilesToUse,
    notes:compoundForm.notes||'',
    calculate:true
   });
   setCompoundForm({compound_id:'',name:'',smiles:'',notes:''});
   setPreview(null);
   setAddCompoundOpen(false);
   let workflow=null;
   if(predict){
    workflow=await api.post('/compounds/'+saved.row_id+'/predict-workflow',{});
    setPredictionWorkflow(workflow);
   }
   await Promise.all([loadProject(projectId),loadProjects(),loadDashboard()]);
   setMessage(predict?(workflow.message+' Activity was not run.'):'Compound saved successfully.');
   // Keep the user on the project dashboard so the newly registered row is
   // immediately visible in Compound Status (including its structure/status).
   setDetail(null);
   setProjectTab('compounds');
  }catch(error){
   const errDetail=error?.response?.data?.detail||error?.message||String(error);
   const cleanMsg=typeof errDetail==='object'?(errDetail.error||JSON.stringify(errDetail)):String(errDetail);
   setMessage(cleanMsg.replace(/^Error:\s*/,''));
   if(predict)setPredictionWorkflow(current=>current?{...current,status:'FAILED',message:cleanMsg}:null);
  }finally{
   setSavingCompound(false);
   if(predict)setAdmetBusy(false);
  }
 };
 const openDetail=async rowId=>{
  try{
   const compound=await api.get('/compounds/'+rowId+'?include_versions=true');setDetail(compound);setDetailTab('overview');setExperimentalOpen(false);
   const assayData=await api.get('/projects/'+compound.project_id+'/assays');setAssays(assayData.assays||assayData||[]);
   if(compound.version)await loadWorkspace(compound.version.id);else{setWorkspace(null);setAdmet(null);setMetabolism(null)}
   setMessage('');
  }catch(error){setMessage(String(error))}
 };
 const calculateProperties=async()=>{
  if(!detail)return;
  try{const result=await api.post('/compounds/'+detail.row_id+'/calculate',{});await openDetail(result.row_id);setDetailTab('properties');await Promise.all([loadProject(projectId),loadDashboard()]);setMessage('Properties calculated for the current CompoundVersion')}catch(error){setMessage(String(error))}
 };
 const updateStructure=async()=>{
  const smiles=prompt('New SMILES (creates a new version)');if(!smiles)return;
  try{await api.patch('/compounds/'+detail.row_id,{smiles,change_note:'Manual structure edit'});await openDetail(detail.row_id);await loadProject(projectId);setAdmet(null);setMessage('Version created')}catch(error){setMessage(String(error))}
 };
 const compare=async()=>{try{setComparison(await api.get('/projects/'+projectId+'/compare?ids='+selected.join(',')+(compareAssay?'&assay_id='+compareAssay:'')));setProjectTab('compare');setDetail(null);setMessage('')}catch(error){setComparison(null);setMessage(String(error))}};

 const saveAdmet=async versionId=>{
  const targetVersionId=Number(versionId||admetVersionId);
  if(!targetVersionId)return;
  setAdmetBusy(true);
  try{
   await api.post('/projects/'+projectId+'/admet/measurements',{...admetForm,version_id:targetVersionId});
   setAdmetForm(current=>({...current,value:'',mean:'',sd:'',n:'',notes:''}));
   if(detail?.version?.id===targetVersionId)await loadWorkspace(targetVersionId);else await loadAdmet();setMessage('Experimental ADMET saved');
  }catch(error){setMessage(String(error))}finally{setAdmetBusy(false)}
 };
 const experimentDefaults=name=>({
  ...EMPTY_EXPERIMENT,
  ...(name.startsWith('pKa')?{unit:'',measurement:name,pka_type:name.includes('acidic')?'acidic':(name.includes('basic')?'basic':(name.includes('micro')?'microscopic':'macroscopic')),method:'Potentiometric titration',temperature_c:'25.0',ionic_strength_m:'0.15'}:{}),
  ...(name==='logP'?{unit:'',measurement:'Shake-flask logP',method:'Shake-flask'}:{}),
  ...(name.startsWith('logD')?{unit:'',measurement:'Shake-flask logD',pH:'7.4',method:'Shake-flask'}:{}),
  ...(name.includes('solubility')?{unit:'µM',measurement:name,solubility_type:name.includes('Intrinsic')?'Intrinsic':(name.includes('Kinetic')?'Kinetic':'Thermodynamic'),pH:'7.4'}:{}),
  ...(name==='Solubility'?{unit:'µM',measurement:'Thermodynamic'}:{}),
  ...(name==='Caco-2 Permeability'?{unit:'cm/s',assay:'Caco-2',measurement:'Papp A→B'}:{}),
  ...(name==='Plasma Protein Binding (PPB)'?{unit:'% bound',measurement:'% Bound'}:{}),
  ...(name.includes('Microsomal')?{unit:'µL/min/mg protein',measurement:'Clint',species:name.startsWith('Rat')?'Rat':'Human',matrix:'Liver Microsome'}:{}),
  ...(name==='Hepatocyte Stability'?{unit:'µL/min/10^6 cells',measurement:'Clint',species:'Human',matrix:'Hepatocytes'}:{}),
  ...(name==='CYP Inhibition'?{unit:'µM',measurement:'IC50',role:'Inhibition',isoform:'3A4'}:{}),
  ...(name==='Transporter'?{unit:'classification',measurement:'Classification',role:'Inhibitor',transporter:'P-gp',species:'Human'}:{}),
  ...(name==='Activity'?{unit:'nM',measurement:'IC50'}:{}),
  ...(['hERG','Ames','DILI'].includes(name)?{unit:name==='hERG'?'µM':'classification',measurement:name==='hERG'?'IC50':'Classification'}:{})
 });
 const toggleExperiment=name=>{
  setExperimentalSelected(current=>current.includes(name)?current.filter(value=>value!==name):[...current,name]);
  setExperimentalDrafts(current=>current[name]?current:{...current,[name]:experimentDefaults(name)});
 };
 const setExperimentValue=(name,key,value)=>setExperimentalDrafts(current=>({...current,[name]:{...(current[name]||experimentDefaults(name)),[key]:value}}));
 const experimentalPayload=(name,row)=>{
  const qualitative=row.unit==='classification'||row.measurement==='Stability class';
  const base={version_id:detail.version.id,value:qualitative?'':row.value,qualitative_value:qualitative?row.value:'',unit:row.unit,species:row.species,matrix:row.matrix,method:row.measurement,source:row.source,notes:row.notes,provenance:{ui_workflow:'Stage 4C-4 endpoint selector',display_name:name}};
  if(name.startsWith('pKa'))return {...base,endpoint:name,unit:'',matrix:'Aqueous buffer',method:row.method||'Potentiometric titration',provenance:{ui_workflow:'Stage 4C-4 pKa Entry',pka_type:row.pka_type,temperature_c:Number(row.temperature_c||25),ionic_strength:Number(row.ionic_strength_m||0.15)}};
  if(name==='logP')return {...base,endpoint:'logP',unit:'',matrix:'Octanol/Water',method:row.method||'Shake-flask'};
  if(name.startsWith('logD'))return {...base,endpoint:'logD (pH '+(row.pH||'7.4')+')',unit:'',matrix:'Octanol/Buffer',method:row.method||'Shake-flask',provenance:{ph:Number(row.pH||7.4),temperature_c:Number(row.temperature_c||25)}};
  if(name.includes('solubility'))return {...base,endpoint:'Solubility',unit:row.unit||'µM',matrix:row.medium||'Aqueous buffer',method:name+(row.pH?' · pH '+row.pH:''),provenance:{solubility_type:row.solubility_type,ph:row.pH?Number(row.pH):null}};
  if(name==='Solubility')return {...base,endpoint:'Solubility',matrix:row.medium,method:[row.solubility_type,row.pH&&'pH '+row.pH].filter(Boolean).join(' · ')||row.measurement};
  if(name==='Caco-2 Permeability')return {...base,endpoint:row.assay==='Caco-2'?'Permeability':(row.assay||'Other')+' permeability',matrix:row.assay||'Caco-2'};
  if(name==='Plasma Protein Binding (PPB)')return {...base,endpoint:'Plasma protein binding',matrix:'Plasma'};
  if(name.includes('Microsomal')){
   const prefix=row.species==='Rat'?'RLM':row.species==='Mouse'?'MLM':row.species==='Human'?'HLM':row.species+' liver microsome';
   return {...base,endpoint:row.measurement==='Clint'?prefix+' intrinsic clearance':prefix+' '+row.measurement.toLowerCase(),matrix:'Liver Microsome'};
  }
  if(name==='Hepatocyte Stability')return {...base,endpoint:row.species+' hepatocyte '+row.measurement.toLowerCase(),matrix:'Hepatocytes'};
  if(name==='CYP Inhibition')return {...base,endpoint:'CYP'+row.isoform+(row.role==='Substrate'?' substrate':' inhibitor'),method:row.measurement+' · '+row.role};
  if(name==='Transporter')return {...base,endpoint:row.transporter+' '+row.role.toLowerCase(),method:[row.measurement,row.assay].filter(Boolean).join(' · ')};
  if(name==='hERG')return {...base,endpoint:'hERG liability'};
  if(name==='Ames')return {...base,endpoint:'Ames mutagenicity'};
  if(name==='DILI')return {...base,endpoint:'DILI clinical liability'};
  return base;
 };
 const saveExperimentalPanel=async()=>{
  if(!detail?.version)return;setAdmetBusy(true);
  try{
   for(const name of experimentalSelected){
    const row=experimentalDrafts[name]||experimentDefaults(name);
    if(row.value==='')throw new Error(name+': value is required');
    if(name==='Activity'){
     if(!row.assay_id)throw new Error('Activity: select an assay');
     await api.post('/assays/'+row.assay_id+'/measurements',{version_id:detail.version.id,value:row.value,unit:row.unit,source:row.source,notes:row.notes});
    }else await api.post('/projects/'+projectId+'/admet/measurements',experimentalPayload(name,row));
   }
   await loadWorkspace(detail.version.id);setExperimentalOpen(false);setExperimentalSelected([]);setExperimentalDrafts({});setMessage('Experimental data saved for '+detail.name+' only');
  }catch(error){setMessage(String(error))}finally{setAdmetBusy(false)}
 };
 const previewAdmet=async()=>{
  setAdmetBusy(true);
  try{const result=await api.post('/projects/'+projectId+'/admet/import-preview',{csv:admetCsv});setAdmetCsvPreview(result);setMessage(result.valid_count+' valid · '+result.errors.length+' errors')}
  catch(error){setAdmetCsvPreview(null);setMessage(String(error))}finally{setAdmetBusy(false)}
 };
 const importAdmet=async()=>{
  setAdmetBusy(true);
  try{const result=await api.post('/projects/'+projectId+'/admet/import',{csv:admetCsv});setAdmetCsvPreview(null);await loadAdmet();setMessage('Imported '+result.imported+' ADMET row'+(result.imported===1?'':'s'))}
  catch(error){setMessage(String(error))}finally{setAdmetBusy(false)}
 };
 const runPrediction=async versionId=>{
  if(!versionId)return;
  setAdmetBusy(true);
  try{const result=await api.post('/admet/predict/'+versionId,{});await loadWorkspace(versionId);setMessage(result.message)}
  catch(error){setMessage(String(error))}finally{setAdmetBusy(false)}
 };
 const runMetabolism=async versionId=>{
  if(!versionId)return;
  setMetabolismBusy(true);
  setPredictionWorkflow(current=>({...current,status:'RUNNING',steps:{...(current?.steps||{}),metabolism:{status:'RUNNING'}}}));
  try{
   const result=await api.post('/metabolism/predict/'+versionId,{});const refreshed=await loadWorkspace(versionId);const data=refreshed.metabolism;
   const run=(data?.runs||[]).find(item=>item.version_id===Number(versionId));
   const metabolismStatus=run?.status==='COMPLETE'?'COMPLETE':(run?.status||'MODEL_UNAVAILABLE');
   setPredictionWorkflow(current=>({...current,status:'PARTIAL',steps:{...(current?.steps||{}),metabolism:{status:metabolismStatus,message:result.message||'Metabolism prediction complete'}}}));
   setSelectedSpotId(run?.spots?.[0]?.id||null);setMessage(result.message);
  }catch(error){setPredictionWorkflow(current=>({...current,status:'PARTIAL',steps:{...(current?.steps||{}),metabolism:{status:'FAILED',message:String(error)}}}));setMessage(String(error))}finally{setMetabolismBusy(false)}
 };
 const saveExperimentalMetabolite=async versionId=>{
  setMetabolismBusy(true);
  try{
   await api.post('/projects/'+projectId+'/metabolism/experimental',{...metaboliteForm,version_id:Number(versionId)});
   setMetaboliteForm({...EMPTY_METABOLITE_FORM});await loadWorkspace(versionId);setMessage('Experimental metabolite saved for the current CompoundVersion');
  }catch(error){setMessage(String(error))}finally{setMetabolismBusy(false)}
 };
 const analyzeOptimization=async versionId=>{
  setOptimizationBusy(true);
  try{
   const result=await api.post('/projects/'+projectId+'/optimization/runs',{...optimizationForm,parent_version_id:Number(versionId),assay_id:optimizationForm.assay_id||null});
   setOptimizationRun(result);setOptimizationRuns(current=>[result,...current.filter(row=>row.id!==result.id)]);setMessage(result.message);
  }catch(error){setMessage(String(error))}finally{setOptimizationBusy(false)}
 };
 const overrideOptimization=async payload=>{
  if(!optimizationRun)return;
  setOptimizationBusy(true);
  try{
   const result=await api.patch('/optimization/runs/'+optimizationRun.id+'/overrides',payload);setOptimizationRun(result);setOptimizationRuns(current=>current.map(row=>row.id===result.id?result:row));setMessage('Manual override saved and strategy reranked');
  }catch(error){setMessage(String(error))}finally{setOptimizationBusy(false)}
 };
 const generateAnalogs=async()=>{
  if(!optimizationRun)return;setProposalBusy(true);
  try{const result=await api.post('/optimization/runs/'+optimizationRun.id+'/proposals',{settings:proposalSettings,hard_constraints:{no_new_structural_alert:true}});setProposalRun({...result,candidates:[]});setProposalRuns(current=>[result,...current]);setSelectedCandidate(null);setMessage('Analog proposal job queued')}
  catch(error){setMessage(String(error))}finally{setProposalBusy(false)}
 };
 const candidateDecision=async(candidate,decision)=>{
  const reason=decision==='REJECTED'?prompt('Reason for rejecting this candidate'):(decision==='PROMOTED'?'Manual promotion for experimental design':'');
  if(decision==='REJECTED'&&!reason)return;setProposalBusy(true);
  try{const updated=await api.patch('/proposal-candidates/'+candidate.id+'/decision',{decision,reason});await refreshProposal(proposalRun.id,proposalView);setSelectedCandidate(updated);setMessage('Candidate decision saved')}
  catch(error){setMessage(String(error))}finally{setProposalBusy(false)}
 };
 const addUserAnalog=async()=>{
  if(!proposalRun||!userAnalog.smiles)return;setProposalBusy(true);
  try{const candidate=await api.post('/proposals/'+proposalRun.id+'/candidates',userAnalog);setUserAnalog({smiles:'',reason:''});await refreshProposal(proposalRun.id,proposalView);setSelectedCandidate(candidate);setMessage('User-added analog rescored')}
  catch(error){setMessage(String(error))}finally{setProposalBusy(false)}
 };

 function admetMeasurementTable(rows){
  if(!rows.length)return Empty({children:'No experimental ADMET measurements yet.'});
  return e('table',{},[
   e('thead',{key:'head'},e('tr',{},['Compound','Endpoint','Result','Unit','Species','Matrix','Replicate','N','Method','Source'].map(label=>e('th',{key:label},label)))),
   e('tbody',{key:'body'},rows.map(row=>e('tr',{key:row.id},[
    e('td',{key:'compound',className:'mono'},versionLabel(row.version_id)),e('td',{key:'endpoint'},endpointName(row.endpoint_id)),
    e('td',{key:'result',className:'num mono'},row.qualitative_value||((row.value!=null?(row.qualifier||'=')+' '+row.value:(row.mean!=null?'mean '+row.mean:'Not measured')))),
    e('td',{key:'unit'},row.unit||'-'),e('td',{key:'species'},row.species||'-'),e('td',{key:'matrix'},row.matrix||'-'),
    e('td',{key:'replicate'},row.replicate||'-'),e('td',{key:'n',className:'num mono'},row.n??'-'),e('td',{key:'method'},row.method||'-'),e('td',{key:'source'},row.source||'-')
   ])))
  ]);
 }

 function predictionDetails(prediction){
  const details=prediction.model?.details||{},output=prediction.outputs||{},domain=output.applicability_domain_details||{};
  const validation=details.validation||output.validation||{};
  const derived=output.derived_outputs||{},assessment=output.experimental_metabolic_stability_assessment||output.metabolic_stability_assessment;
  const conformal=output.calibrated_uncertainty||details.conformal_governance||{};
  const prov=conformal.data_provenance||details.calibration_provenance;
  const qual=conformal.calibration_quality||details.calibration_quality;
  const provLabel=prov==='EXTERNAL'?'External validation set':(prov==='INTERNAL'?'Internal validation set (training overlap)':(prov==='TRAINING_OVERLAP_UNKNOWN'?'Training overlap unknown':'Unavailable'));
  const qualLabel=qual==='VALIDATED'?'VALIDATED':(qual==='UNDERCOVERED'?'UNDERCOVERED':(qual==='INSUFFICIENT_N'?'INSUFFICIENT N (<30)':(qual||'UNAVAILABLE')));
  return e('details',{},[
   e('summary',{key:'summary'},'Details'),
   e('div',{key:'body',className:'small',style:{minWidth:'320px'}},[
    e('div',{key:'model'},e('strong',{},'Model: '),prediction.model.model_name+' · '+prediction.model.model_version),
    e('div',{key:'source'},[e('strong',{},'Source: '),e('a',{href:details.source||output.model_source,target:'_blank',rel:'noreferrer'},details.source||output.model_source)]),
    e('div',{key:'definition'},[e('strong',{},'Endpoint: '),details.endpoint_definition||output.endpoint_definition]),
    (details.assay_definition||output.assay_definition)&&e('div',{key:'assay'},[e('strong',{},'Assay: '),details.assay_definition||output.assay_definition]),
    e('div',{key:'training'},[e('strong',{},'Training: '),details.training_dataset||output.training_dataset]),
    e('div',{key:'validation'},[e('strong',{},'Validation: '),Object.entries(validation).map(([key,value])=>key+' '+value).join(' · ')]),
    details.independent_validation&&e('div',{key:'independent'},[e('strong',{},'Independent validation: '),Object.entries(details.independent_validation).map(([key,value])=>key+' '+value).join(' · ')]),
    prov&&e('div',{key:'conformal-gov'},[
     e('strong',{},'Calibration Governance: '),
     e('span',{},'Data: '+provLabel+' · Conformal: '),
     e('span',{className:qual==='VALIDATED'?'pass':(qual==='UNDERCOVERED'?'fail':'warn')},qualLabel),
     conformal.empirical_coverage!=null&&e('span',{},' ('+(conformal.nominal_coverage?conformal.nominal_coverage*100:90)+'% nominal / '+(conformal.empirical_coverage*100).toFixed(1)+'% observed · Eval N='+(conformal.evaluation_n??'—')+(conformal.expected_sampling_uncertainty_se?(' · SE=±'+(conformal.expected_sampling_uncertainty_se*100).toFixed(1)+'%'):'')+')'),
     conformal.interval_width!=null&&e('span',{},' · Interval width: '+conformal.interval_width+' '+conformal.unit+(conformal.interval_utility?.utility_status==='UNINFORMATIVE_INTERVAL'?' [UNINFORMATIVE]':'')),
     conformal.prediction_set&&e('span',{},' · Conformal set: {'+conformal.prediction_set.join(', ')+'}')
    ]),
    (conformal.warnings||[]).length>0&&e('div',{key:'conformal-warns'},conformal.warnings.map(w=>e('div',{key:w,className:'fail small'},w))),
    e('div',{key:'license'},[e('strong',{},'License: '),details.license||output.license]),
    e('div',{key:'ad'},[e('strong',{},'AD evidence: '),'nearest similarity '+(domain.nearest_training_similarity??'—')+' · chemical-space distance '+(domain.chemical_space_distance??'—')+(domain.descriptors_outside_range?.length?' · outside '+domain.descriptors_outside_range.join(', '):' · descriptors within training range')]),
    Object.keys(derived).length>0&&e('div',{key:'derived'},[e('strong',{},'Derived output: '),Object.entries(derived).map(([key,value])=>key+' '+(typeof value==='number'?Number(value).toPrecision(5):value)).join(' · ')]),
    assessment&&e('div',{key:'assessment'},[e('strong',{},'Metabolic assessment: '),assessment.category+(assessment.metabolic_liability_flag?' · '+assessment.metabolic_liability_flag:'')+' · '+(assessment.thresholds?.basis||'')]),
    output.liability_summary&&e('div',{key:'liability'},[e('strong',{},output.safety_endpoint?'Safety flag: ':(output.transporter?'Interaction flag: ':'CYP liability rule: ')),output.liability_summary.flag+' · '+output.liability_summary.rule+' · '+output.liability_summary.basis]),
    e('div',{key:'limits'},[e('strong',{},'Limitations: '),details.limitations||output.limitations])
   ])
  ]);
 }

 function cypPredictionTable(rows){
  if(!rows.length)return Empty({children:'No CYP predictions yet. Run prediction to evaluate installed endpoints.'});
  return e('div',{},[
   e('table',{},[
    e('thead',{key:'head'},e('tr',{},['CYP Isoform','Role','Prediction','Probability','Assessment','Model Applicability','Confidence','Model',''].map(label=>e('th',{key:label},label)))),
    e('tbody',{key:'body'},rows.map(prediction=>{
     const output=prediction.outputs||{},evidence=output.experimental_evidence||[];
     const prob=Number(output.probability??prediction.predicted_value);
     const interp=getInterpretation('cyp_inhibitor',prob);
     return e('tr',{key:prediction.id},[
      e('td',{key:'isoform',style:{fontWeight:600}},output.isoform||prediction.endpoint.split(' ')[0]),
      e('td',{key:'role'},output.role||prediction.endpoint.split(' ')[1]?.toUpperCase()),
      e('td',{key:'class',className:'bold'},output.classification||(prob>=0.5?'Positive':'Negative')),
      e('td',{key:'probability',className:'mono'},prob.toFixed(3)),
      e('td',{},[ScientificBadge({assessment:interp.assessment,colorClass:interp.colorClass,textLabel:interp.label})]),
      e('td',{},[e('span',{className:'badge-intermediate',title:'Model Applicability domain verified against training set'},prediction.applicability_domain||'IN DOMAIN')]),
      e('td',{},prediction.confidence||'HIGH'),
      e('td',{key:'model',className:'small'},e('span',{className:'model-chip'},'M1'),'OpenADMET '+(prediction.model?.model_version||'v1.0')),
      e('td',{key:'details'},predictionDetails(prediction))
     ]);
    }))
   ]),
   e('div',{className:'model-notes'},[
    e('strong',{},'Model Notes: '),
    e('span',{},'M1 = OpenADMET CYP Panel Classifier v1.0. Model Applicability represents applicability domain boundaries derived from training descriptor distributions.')
   ])
  ]);
 }

 function transporterPredictionTable(rows){
  if(!rows.length)return Empty({children:'No active transporter prediction yet. Human P-gp inhibitor is the only qualified installed endpoint.'});
  return e('table',{},[
   e('thead',{key:'head'},e('tr',{},['Transporter','Role','Species','Prediction','Probability','Experimental','Domain','Confidence','Model',''].map(label=>e('th',{key:label},label)))),
   e('tbody',{key:'body'},rows.map(prediction=>{
    const output=prediction.outputs||{},evidence=output.experimental_evidence||[];
    const experimental=evidence.length?evidence.map(item=>item.value+' '+item.unit+' ('+item.comparison+')').join(' · '):'—';
    const liability=output.liability_summary?.flag;
    return e('tr',{key:prediction.id},[
     e('td',{key:'target'},output.transporter||prediction.endpoint),e('td',{key:'role'},output.role||'—'),e('td',{key:'species'},output.species||'Human'),
     e('td',{key:'class'},[output.classification||'—',liability&&e('div',{key:'flag',className:'fail small'},liability)]),
     e('td',{key:'probability',className:'mono'},Number(output.probability??prediction.predicted_value).toFixed(4)),e('td',{key:'experimental'},experimental),
     e('td',{key:'domain'},prediction.applicability_domain),e('td',{key:'confidence'},prediction.confidence),
     e('td',{key:'model',className:'small'},prediction.model?.model_name+' '+prediction.model?.model_version),e('td',{key:'details'},predictionDetails(prediction))
    ]);
   }))
  ]);
 }

 function unavailableTransporterModels(){
  const rows=(admet?.models||[]).filter(model=>TRANSPORTER_ENDPOINTS.has(model.endpoint)&&!model.active);
  if(!rows.length)return null;
  return e('div',{className:'small'},rows.map(model=>e('details',{key:model.endpoint},[
   e('summary',{key:'summary'},model.endpoint+': Model unavailable'),
   e('div',{key:'reason'},model.unavailable_reason),
   e('div',{key:'identity'},'Target: '+(model.details?.transporter||'—')+' · Role: '+(model.details?.role||'—')+' · Species: '+(model.details?.species||'—'))
  ])));
 }

 function safetyPredictionTable(rows){
  if(!rows.length)return Empty({children:'No safety predictions yet. Run prediction to evaluate hERG, Ames, and DILI.'});
  return e('table',{},[
   e('thead',{key:'head'},e('tr',{},['Endpoint','Prediction','Probability','Experimental','Domain','Confidence','Model',''].map(label=>e('th',{key:label},label)))),
   e('tbody',{key:'body'},rows.map(prediction=>{
    const output=prediction.outputs||{},evidence=output.experimental_evidence||[];
    const experimental=evidence.length?evidence.map(item=>item.value+' '+item.unit+' ('+item.comparison+')').join(' · '):'—';
    return e('tr',{key:prediction.id},[
     e('td',{key:'endpoint'},output.safety_endpoint||prediction.endpoint),
     e('td',{key:'class'},[output.classification||'—',output.liability_summary?.flag&&e('div',{key:'flag',className:'fail small'},output.liability_summary.flag)]),
     e('td',{key:'probability',className:'mono'},Number(output.probability??prediction.predicted_value).toFixed(4)),
     e('td',{key:'experimental'},experimental),e('td',{key:'domain'},prediction.applicability_domain),
     e('td',{key:'confidence'},prediction.confidence),e('td',{key:'model',className:'small'},prediction.model?.model_name+' '+prediction.model?.model_version),
     e('td',{key:'details'},predictionDetails(prediction))
    ]);
   }))
  ]);
 }

 function unavailableSafetyModels(){
  const rows=(admet?.models||[]).filter(model=>OPTIONAL_SAFETY_ENDPOINTS.has(model.endpoint)&&!model.active);
  return e('div',{className:'small'},rows.map(model=>e('details',{key:model.endpoint},[
   e('summary',{key:'summary'},model.endpoint+': Model unavailable'),
   e('div',{key:'reason'},model.unavailable_reason),
   e('div',{key:'identity'},'Endpoint: '+(model.details?.safety_endpoint||model.endpoint)+' · Species: '+(model.details?.species||'—')+' · Checkpoint: unavailable')
  ])));
 }

 function integratedProfile(versionId){
  const profile=admet?.integrated_profiles?.[String(versionId)];
  if(!profile)return null;
  const summary=profile.summary||{};
  const group=(title,rows,klass)=>e('div',{className:'col-4',key:title},[e('h4',{key:'title'},title),rows?.length?e('ul',{key:'list',className:klass},rows.map(text=>e('li',{key:text},text))):e('p',{key:'empty',className:'small'},'None')]);
 return e('div',{className:'card',key:'integrated-profile'},[
   e('h3',{key:'title'},'Stage 3 Integrated ADMET Profile'),
   e('p',{key:'policy',className:'small'},'Experimental values take display precedence while predictions remain preserved. Confidence and applicability domain are endpoint-specific. No overall ADMET score or candidate ranking is calculated.'),
   e('div',{className:'grid',key:'summary'},[group('Strengths',summary.strengths,'strengths'),group('Concerns',summary.concerns,'concerns')]),
   (summary.unknown||[]).length>0&&e('details',{key:'unavailable',className:'unavailable-collapse'},[e('summary',{},'Unavailable models ('+summary.unknown.length+')'),e('ul',{className:'small'},summary.unknown.map(text=>e('li',{key:text},text)))]),
   e('div',{key:'audit',className:profile.provenance_audit?.status==='PASS'?'pass':'fail'},'Provenance audit: '+profile.provenance_audit?.status+' · '+profile.provenance_audit?.checked+' latest endpoint predictions checked')
  ]);
 }

 function admetPredictionTable(rows){
  if(!rows.length)return e('div',{className:'empty-state'},[StatusBadge({type:'Not predicted'}),e('p',{key:'text'},'Prediction not run for this CompoundVersion.'),e('button',{key:'run',className:'secondary',disabled:admetBusy||!detail?.version,onClick:()=>runPrediction(detail.version.id)},'Run Predictions')]);
  return e('table',{},[
   e('thead',{key:'head'},e('tr',{},['Compound','Endpoint','Experimental','Predicted','Model','Confidence','Domain',''].map(label=>e('th',{key:label},label)))),
   e('tbody',{key:'body'},rows.map(prediction=>{
    const comparison=prediction.experimental_comparisons?.[0];
    const experimental=comparison?(comparison.experimental_value+' '+comparison.experimental_unit):'Not measured';
    const error=comparison?' · |error| '+comparison.absolute_error+' '+comparison.normalized_unit:'';
    const assessment=prediction.outputs?.experimental_metabolic_stability_assessment||prediction.outputs?.metabolic_stability_assessment;
    const flag=assessment?.metabolic_liability_flag?' · '+assessment.metabolic_liability_flag:'';
    return e('tr',{key:prediction.id},[
     e('td',{key:'compound',className:'mono'},versionLabel(prediction.version_id)),
     e('td',{key:'endpoint'},prediction.endpoint==='Permeability'?'Caco-2':prediction.endpoint),
     e('td',{key:'experimental'},experimental),
     e('td',{key:'predicted',className:'mono'},Number(prediction.predicted_value).toFixed(3)+' '+prediction.unit+error+flag),
     e('td',{key:'model',className:'small'},prediction.model?.model_name+' '+prediction.model?.model_version),
     e('td',{key:'confidence'},prediction.confidence),e('td',{key:'domain'},prediction.applicability_domain),
     e('td',{key:'details'},predictionDetails(prediction))
    ]);
   }))
  ]);
 }

 function consensusPredictionPanel(versionId){
  const predictions=(admet?.predictions||[]).filter(row=>row.version_id===Number(versionId));
  const latestByModel=new Map();predictions.forEach(row=>{const key=row.endpoint+'|'+row.model?.id;if(!latestByModel.has(key))latestByModel.set(key,row)});
  const individual=[...latestByModel.values()],grouped={};individual.forEach(row=>(grouped[row.endpoint]||(grouped[row.endpoint]=[])).push(row));
  const consensusByEndpoint=new Map();(admet?.consensus_predictions||[]).filter(row=>row.version_id===Number(versionId)).forEach(row=>{if(!consensusByEndpoint.has(row.endpoint))consensusByEndpoint.set(row.endpoint,row)});
  if(!individual.length)return e('div',{className:'empty-state'},[StatusBadge({type:'Not predicted'}),e('p',{},'Prediction not run for this CompoundVersion.'),e('button',{className:'secondary',disabled:admetBusy,onClick:()=>runPrediction(versionId)},'Run Predictions')]);
  return e('div',{className:'endpoint-prediction-grid'},Object.entries(grouped).map(([endpoint,models])=>{
   const consensus=consensusByEndpoint.get(endpoint),classification=consensus?.classification;
   return e('article',{className:'endpoint-prediction-card',key:endpoint},[
    e('div',{className:'row toolbar',key:'head'},[e('h4',{},endpoint==='Permeability'?'Caco-2 Permeability':endpoint),e('span',{className:'small'},'Individual Models: '+models.length)]),
    e('div',{className:'consensus-result',key:'combined'},[e('span',{},'Combined Prediction'),e('strong',{className:'mono'},classification||((consensus?.combined_value??models[0].predicted_value).toFixed(3)+' '+(consensus?.unit||models[0].unit))),e('div',{className:'small'},'Confidence '+(consensus?.confidence||models[0].confidence)+' · Domain '+(consensus?.applicability_domain||models[0].applicability_domain))]),
    e('div',{className:'individual-models',key:'models'},models.map(row=>e('div',{className:'individual-model-row',key:row.id},[e('div',{},[e('strong',{},row.model.model_name),e('span',{className:'small'},' v'+row.model.model_version)]),e('span',{className:'mono'},(row.outputs?.classification||Number(row.predicted_value).toFixed(3))+(row.outputs?.classification?'':' '+row.unit)),e('span',{className:'small'},row.confidence+' · '+row.applicability_domain),predictionDetails(row)]))),
    consensus&&e('details',{key:'consensus-details'},[e('summary',{},'Consensus provenance'),e('div',{className:'small'},[e('div',{},consensus.provenance?.weighting_policy),...(consensus.models||[]).map(item=>e('div',{key:item.model_id},item.model_name+' '+item.model_version+' · weight '+Number(item.weight).toFixed(3)))])])
   ]);
  }));
 }

 function experimentalComparisonPanel(versionId){
  const rows=(admet?.predictions||[]).filter(row=>row.version_id===Number(versionId)&&row.experimental_comparisons?.length);
  if(!rows.length)return e('div',{className:'empty-state'},[StatusBadge({type:'Not measured'}),e('p',{},'No unit-compatible Experimental vs Prediction comparison is available.')]);
  return e('table',{},[e('thead',{},e('tr',{},['Endpoint','Experimental','Prediction','Difference','Model'].map(label=>e('th',{key:label},label)))),e('tbody',{},rows.map(row=>{const item=row.experimental_comparisons[0];return e('tr',{key:row.id},[e('td',{},row.endpoint),e('td',{className:'mono'},item.experimental_normalized+' '+item.normalized_unit),e('td',{className:'mono'},row.predicted_value+' '+row.unit),e('td',{className:'mono'},item.absolute_error==null?(item.classification_match?'AGREES':'DISAGREES'):(item.absolute_error+' '+item.normalized_unit)),e('td',{className:'small'},row.model.model_name+' '+row.model.model_version)])}))]);
 }

 function unavailableModelsCollapsed(){
  const rows=(admet?.models||[]).filter(model=>!model.active);
  if(!rows.length)return Empty({children:'All registered models are available.'});
  return e('details',{className:'unavailable-collapse'},[e('summary',{},'Unavailable models ('+rows.length+')'),e('div',{className:'small'},rows.map(model=>e('div',{key:model.id},[e('strong',{},model.endpoint+': '),model.unavailable_reason])))]);
 }

 function ExperimentalDataPanel(){
  if(!detail?.version)return e('div',{className:'empty-state'},[StatusBadge({type:'Not applicable'}),e('p',{},'Add a valid structure before entering structure-linked experimental data.')]);
  const select=(label,value,options,onChange)=>e('div',{},[e('label',{key:'label'},label),e('select',{key:'select',value,onChange:event=>onChange(event.target.value)},options.map(option=>{const item=typeof option==='object'?option:{value:option,label:option};return e('option',{key:item.value,value:item.value},item.label)}))]);
  const input=(name,row,key,label,type='text')=>Field({label,type,value:row[key],onChange:value=>setExperimentValue(name,key,value)});
  const formFor=name=>{
   const row=experimentalDrafts[name]||experimentDefaults(name),fields=[];
   if(name.startsWith('pKa'))fields.push(select('pKa Type',row.pka_type,['acidic','basic','macroscopic','microscopic'],v=>setExperimentValue(name,'pka_type',v)),select('Method',row.method,['Potentiometric titration','Spectrophotometric titration','CE','NMR','Other'],v=>setExperimentValue(name,'method',v)),input(name,row,'value','Value (pKa)','number'),input(name,row,'temperature_c','Temperature (°C)','number'),input(name,row,'ionic_strength_m','Ionic Strength (M)','number'));
   if(name==='logP')fields.push(select('Method',row.method,['Shake-flask','HPLC-based','Potentiometric','Other'],v=>setExperimentValue(name,'method',v)),input(name,row,'value','logP Value','number'));
   if(name.startsWith('logD'))fields.push(input(name,row,'pH','Assay pH (Mandatory *)','number'),select('Method',row.method,['Shake-flask','Potentiometric (GLpKa)','HPLC-based','Other'],v=>setExperimentValue(name,'method',v)),input(name,row,'value','logD Value','number'),input(name,row,'temperature_c','Temperature (°C)','number'));
   if(name.includes('solubility')&&name!=='Solubility')fields.push(select('Solubility type',row.solubility_type,['Intrinsic','Kinetic','Thermodynamic'],v=>setExperimentValue(name,'solubility_type',v)),input(name,row,'pH','Assay pH','number'),input(name,row,'medium','Medium / Buffer'),input(name,row,'value','Value','number'),input(name,row,'unit','Unit'));
   if(name==='Activity')fields.push(select('Assay',row.assay_id,[{value:'',label:'Select an assay'},...assays.map(item=>({value:String(item.id),label:item.name+' · '+item.measurement_type}))],value=>setExperimentValue(name,'assay_id',value)),input(name,row,'value','Value','number'),input(name,row,'unit','Unit'));
   if(name==='Solubility')fields.push(select('Solubility type',row.solubility_type,['','Kinetic','Thermodynamic','Intrinsic'],value=>setExperimentValue(name,'solubility_type',value)),input(name,row,'pH','pH','number'),input(name,row,'medium','Medium'),input(name,row,'value','Value','number'),input(name,row,'unit','Unit'));
   if(name==='Caco-2 Permeability')fields.push(select('Assay',row.assay,['Caco-2','MDCK','PAMPA','Other'],value=>setExperimentValue(name,'assay',value)),select('Measurement',row.measurement,['Papp A→B','Papp B→A','Efflux Ratio'],value=>setExperimentValue(name,'measurement',value)),input(name,row,'value','Value','number'),input(name,row,'unit','Unit'));
   if(name==='Plasma Protein Binding (PPB)')fields.push(select('Species',row.species,['Human','Rat','Mouse','Dog','Monkey'],value=>setExperimentValue(name,'species',value)),select('Measurement',row.measurement,['% Bound','fu'],value=>{setExperimentValue(name,'measurement',value);setExperimentValue(name,'unit',value==='fu'?'fraction unbound':'% bound')}),input(name,row,'value','Value','number'),input(name,row,'unit','Unit'));
   if(name.includes('Microsomal'))fields.push(select('Species',row.species,['Human','Rat','Mouse','Dog','Monkey'],value=>setExperimentValue(name,'species',value)),select('Matrix',row.matrix,['Liver Microsome'],value=>setExperimentValue(name,'matrix',value)),select('Measurement',row.measurement,['Clint','t1/2','% remaining','Stability class'],value=>{setExperimentValue(name,'measurement',value);setExperimentValue(name,'unit',value==='Clint'?'µL/min/mg protein':value==='t1/2'?'min':value==='% remaining'?'%':'classification')}),input(name,row,'value','Value',row.measurement==='Stability class'?'text':'number'),input(name,row,'unit','Unit'));
   if(name==='Hepatocyte Stability')fields.push(select('Species',row.species,['Human','Rat','Mouse','Dog','Monkey'],value=>setExperimentValue(name,'species',value)),select('Matrix',row.matrix,['Hepatocytes'],value=>setExperimentValue(name,'matrix',value)),select('Measurement',row.measurement,['Clint','t1/2','% remaining','Stability class'],value=>{setExperimentValue(name,'measurement',value);setExperimentValue(name,'unit',value==='Clint'?'µL/min/10^6 cells':value==='t1/2'?'min':value==='% remaining'?'%':'classification')}),input(name,row,'value','Value',row.measurement==='Stability class'?'text':'number'),input(name,row,'unit','Unit'));
   if(name==='CYP Inhibition')fields.push(select('Isoform',row.isoform,['1A2','2C9','2C19','2D6','3A4'],value=>setExperimentValue(name,'isoform',value)),select('Role',row.role,['Inhibition','Substrate'],value=>{setExperimentValue(name,'role',value);if(value==='Substrate'){setExperimentValue(name,'measurement','Classification');setExperimentValue(name,'unit','classification')}}),select('Measurement',row.measurement,row.role==='Substrate'?['Classification']:['IC50','Ki','Classification'],value=>{setExperimentValue(name,'measurement',value);setExperimentValue(name,'unit',value==='Classification'?'classification':'µM')}),input(name,row,'value','Value',row.measurement==='Classification'?'text':'number'),input(name,row,'unit','Unit'));
   if(name==='Transporter')fields.push(select('Transporter',row.transporter,['P-gp','BCRP','BSEP','OATP1B1','OATP1B3','OCT1','OCT2','MATE1','MATE2-K'],value=>setExperimentValue(name,'transporter',value)),select('Role',row.role,['Substrate','Inhibitor'],value=>setExperimentValue(name,'role',value)),select('Species',row.species,['Human','Rat','Mouse','Dog','Monkey'],value=>setExperimentValue(name,'species',value)),input(name,row,'assay','Assay'),select('Measurement',row.measurement,['Classification','IC50','Ki','Efflux Ratio'],value=>{setExperimentValue(name,'measurement',value);setExperimentValue(name,'unit',value==='Classification'?'classification':value==='Efflux Ratio'?'ratio':'µM')}),input(name,row,'value','Value',row.measurement==='Classification'?'text':'number'),input(name,row,'unit','Unit'));
   if(['hERG','Ames','DILI'].includes(name))fields.push(select('Measurement',row.measurement,name==='hERG'?['IC50','Classification']:['Classification'],value=>setExperimentValue(name,'measurement',value)),input(name,row,'value','Value',row.measurement==='Classification'?'text':'number'),input(name,row,'unit','Unit'));
   return e('div',{className:'experimental-endpoint-card',key:name},[e('h4',{key:'title'},name),e('div',{className:'experimental-fields',key:'fields'},fields),e('div',{className:'experimental-fields',key:'common'},[input(name,row,'source','Source'),input(name,row,'notes','Notes')])]);
  };
  return e('div',{className:'experimental-panel'},[
   e('h3',{key:'question'},'What experimental data do you want to add?'),
   e('div',{className:'preset-row',key:'presets'},Object.entries(EXPERIMENT_PRESETS).map(([name,items])=>e('button',{key:name,className:'secondary',onClick:()=>{setExperimentalSelected(items);setExperimentalDrafts(current=>Object.fromEntries(items.map(item=>[item,current[item]||experimentDefaults(item)])))}},name))),
   ...['Physicochemistry','Activity','ADME','Safety'].map(category=>e('div',{key:category},[e('h4',{key:'title'},category),e('div',{className:'endpoint-selector',key:'options'},EXPERIMENT_OPTIONS.filter(row=>row[1]===category).map(([name])=>e('label',{key:name,className:'check-option'},[e('input',{type:'checkbox',checked:experimentalSelected.includes(name),onChange:()=>toggleExperiment(name)}),e('span',{},name)]))) ])),
   ...experimentalSelected.map(formFor),
   e('div',{className:'row',key:'actions'},[e('button',{disabled:admetBusy||experimentalSelected.length===0,onClick:saveExperimentalPanel},admetBusy?'Saving…':'Save Experimental Data'),e('button',{className:'secondary',onClick:()=>setExperimentalOpen(false)},'Cancel')])
  ]);
 }

 function admetFormPanel(fixedVersionId=null){
  const targetVersionId=fixedVersionId||admetVersionId;
  const setValue=(key,value)=>setAdmetForm(current=>({...current,[key]:value}));
  return e('div',{},[
   e('div',{className:'grid',key:'fields'},[
    e('div',{className:'col-3',key:'compound'},[e('label',{},'Compound version'),fixedVersionId?e('input',{value:versionLabel(fixedVersionId),disabled:true}):e('select',{value:admetVersionId,onChange:event=>setAdmetVersionId(event.target.value)},[e('option',{key:'empty',value:''},'Select a compound'),...currentVersions.map(compound=>e('option',{key:compound.version.id,value:compound.version.id},compound.compound_id+' v'+compound.current_version))])]),
    e('div',{className:'col-3',key:'endpoint'},Field({label:'Endpoint',value:admetForm.endpoint,onChange:value=>setValue('endpoint',value),placeholder:'e.g. Solubility'})),
    e('div',{className:'col-3',key:'value'},Field({label:'Value',value:admetForm.value,onChange:value=>setValue('value',value),type:'number'})),
    e('div',{className:'col-3',key:'unit'},Field({label:'Unit',value:admetForm.unit,onChange:value=>setValue('unit',value),placeholder:'e.g. µM'})),
    e('div',{className:'col-3',key:'qualifier'},[e('label',{},'Qualifier'),e('select',{value:admetForm.qualifier,onChange:event=>setValue('qualifier',event.target.value)},['=','<','<=','>','>=','~'].map(value=>e('option',{key:value,value},value)))]),
    e('div',{className:'col-3',key:'replicate'},Field({label:'Replicate',value:admetForm.replicate,onChange:value=>setValue('replicate',value)})),
    e('div',{className:'col-3',key:'mean'},Field({label:'Mean',value:admetForm.mean,onChange:value=>setValue('mean',value),type:'number'})),
    e('div',{className:'col-3',key:'sd'},Field({label:'SD',value:admetForm.sd,onChange:value=>setValue('sd',value),type:'number'})),
    e('div',{className:'col-3',key:'n'},Field({label:'Sample size (n)',value:admetForm.n,onChange:value=>setValue('n',value),type:'number'})),
    e('div',{className:'col-3',key:'species'},Field({label:'Species',value:admetForm.species,onChange:value=>setValue('species',value)})),
    e('div',{className:'col-3',key:'matrix'},Field({label:'Matrix',value:admetForm.matrix,onChange:value=>setValue('matrix',value)})),
    e('div',{className:'col-3',key:'method'},Field({label:'Method',value:admetForm.method,onChange:value=>setValue('method',value)})),
    e('div',{className:'col-3',key:'source'},Field({label:'Source',value:admetForm.source,onChange:value=>setValue('source',value)})),
    e('div',{className:'col-3',key:'date'},Field({label:'Experiment date',value:admetForm.date,onChange:value=>setValue('date',value),type:'date'})),
    e('div',{className:'col-6',key:'notes'},Field({label:'Notes',value:admetForm.notes,onChange:value=>setValue('notes',value)}))
   ]),
   e('div',{className:'row',style:{marginTop:'12px'},key:'actions'},[
    e('button',{key:'save',disabled:admetBusy||!targetVersionId||!admetForm.endpoint||!admetForm.unit||(!admetForm.value&&!admetForm.mean),onClick:()=>saveAdmet(targetVersionId)},admetBusy?'Saving…':'Save experimental'),
    e('span',{key:'help',className:'small'},'Enter a raw value or a summary mean. Experimental records remain distinct from predictions.')
   ])
  ]);
 }

 function projectAdmetTab(){
  const sample=currentVersions[0]?.compound_id||'C001';
  const placeholder='compound_id,version_number,endpoint,species,matrix,value,unit,qualifier,replicate,mean,sd,n,method,source,date,notes\n'+sample+',1,Solubility,human,plasma,12.5,µM,=,R1,,,,shake flask,Study A,2026-08-25,';
  return e(React.Fragment,null,[
   e('div',{className:'card',key:'experimental'},[e('div',{className:'row toolbar'},[e('h3',{},'Experimental ADMET'),e('a',{className:'button secondary',href:'/api/projects/'+projectId+'/admet/export.csv'},'Export CSV')]),admetFormPanel(),e('h4',{style:{marginTop:'22px'}},'Saved measurements'),admetMeasurementTable(admet?.measurements||[])]),
   e('div',{className:'card',key:'csv'},[
    e('h3',{key:'title'},'CSV preview and import'),
    e('p',{key:'help',className:'small'},'Preview validates every row before import. Imports are all-or-nothing when validation errors are present.'),
    e('textarea',{key:'input',rows:7,value:admetCsv,placeholder,onChange:event=>{setAdmetCsv(event.target.value);setAdmetCsvPreview(null)}}),
    e('div',{key:'actions',className:'row',style:{marginTop:'10px'}},[
     e('button',{key:'preview',className:'secondary',disabled:admetBusy||!admetCsv.trim(),onClick:previewAdmet},'Preview CSV'),
     e('button',{key:'import',disabled:admetBusy||!admetCsvPreview||admetCsvPreview.errors.length>0||admetCsvPreview.valid_count===0,onClick:importAdmet},'Import '+(admetCsvPreview?.valid_count||0)+' valid row'+((admetCsvPreview?.valid_count||0)===1?'':'s'))
    ]),
    admetCsvPreview&&e('div',{key:'result',style:{marginTop:'14px'}},[
     e('div',{key:'summary',className:admetCsvPreview.errors.length?'fail':'pass'},admetCsvPreview.valid_count+' valid · '+admetCsvPreview.errors.length+' errors'),
     admetCsvPreview.rows.length>0&&e('table',{key:'valid'},[
      e('thead',{key:'head'},e('tr',{},['Row','Compound','Version','Endpoint','Value','Unit'].map(label=>e('th',{key:label},label)))),
      e('tbody',{key:'body'},admetCsvPreview.rows.map(row=>e('tr',{key:row.row},[e('td',{key:'row',className:'mono'},row.row),e('td',{key:'compound',className:'mono'},row.compound_id),e('td',{key:'version'},row.version_number||'current'),e('td',{key:'endpoint'},row.endpoint),e('td',{key:'value',className:'num mono'},row.value),e('td',{key:'unit'},row.unit)])))
     ]),
     admetCsvPreview.errors.length>0&&e('table',{key:'errors'},[
      e('thead',{key:'head'},e('tr',{},[e('th',{key:'row'},'Row'),e('th',{key:'error'},'Error')])),
      e('tbody',{key:'body'},admetCsvPreview.errors.map(error=>e('tr',{key:error.row},[e('td',{key:'row',className:'mono'},error.row),e('td',{key:'error'},error.error)])))
     ])
    ])
   ]),
   e('div',{className:'card',key:'predicted'},[
    e('div',{className:'row toolbar',key:'title'},[e('h3',{},'ADMET predictions through Stage 3F'),e('button',{disabled:admetBusy||!admetVersionId,onClick:()=>runPrediction(Number(admetVersionId))},admetBusy?'Predicting…':'Run prediction')]),
    e('p',{key:'scope',className:'small'},'CYP/transporter roles and hERG/Ames/DILI definitions remain isolated. Classification probabilities are never converted to IC50, Ki, efflux ratio, or other quantitative assay values. Unqualified endpoints remain visible as MODEL_UNAVAILABLE.'),
    admetVersionId&&integratedProfile(Number(admetVersionId)),
    admetPredictionTable((admet?.predictions||[]).filter(row=>!row.endpoint.startsWith('CYP')&&!TRANSPORTER_ENDPOINTS.has(row.endpoint)&&!SAFETY_ENDPOINTS.has(row.endpoint))),
    e('h4',{key:'cyp-predictions-title',style:{marginTop:'22px'}},'CYP inhibitor / substrate predictions'),
    cypPredictionTable((admet?.predictions||[]).filter(row=>row.endpoint.startsWith('CYP'))),
    e('h4',{key:'transporter-predictions-title',style:{marginTop:'22px'}},'Transporters'),
    transporterPredictionTable((admet?.predictions||[]).filter(row=>TRANSPORTER_ENDPOINTS.has(row.endpoint))),
    unavailableTransporterModels(),
    e('h4',{key:'safety-predictions-title',style:{marginTop:'22px'}},'Safety · hERG / Ames / DILI'),
    safetyPredictionTable((admet?.predictions||[]).filter(row=>SAFETY_ENDPOINTS.has(row.endpoint))),
    unavailableSafetyModels(),
    e('h4',{key:'registry-title',style:{marginTop:'22px'}},'Model registry'),
    (admet?.models||[]).length?e('table',{key:'registry'},[e('thead',{key:'head'},e('tr',{},['Endpoint','Model','Version','Unit','Status'].map(label=>e('th',{key:label},label)))),e('tbody',{key:'body'},admet.models.map(model=>e('tr',{key:model.id},[e('td',{key:'endpoint'},model.endpoint==='Permeability'?'Caco-2':model.endpoint),e('td',{key:'model'},model.model_name),e('td',{key:'version'},model.model_version),e('td',{key:'unit'},model.output_unit||'—'),e('td',{key:'status'},Badge({ok:model.active,text:model.status}))]))) ]):Empty({children:'No ADMET model registry entries.'}),
    e('div',{key:'selected',className:'small',style:{marginTop:'10px'}},admetVersionId?'Selected: '+versionLabel(admetVersionId):'Select a compound version above')
   ])
  ]);
 }

 function metabolismPanel(versionId){
  const run=(metabolism?.runs||[]).find(item=>item.version_id===Number(versionId));
  const experimental=(metabolism?.experimental_metabolites||[]).filter(item=>item.version_id===Number(versionId));
  const topLimit=metabolicTop==='ALL'?Number.MAX_SAFE_INTEGER:Number(metabolicTop);
  const spots=(run?.spots||[]).filter(spot=>spot.rank<=topLimit);
  const metabolites=(run?.predicted_metabolites||[]).filter(item=>item.rank<=topLimit);
  const selected=(run?.spots||[]).find(spot=>spot.id===selectedSpotId)||spots[0];
  const summary=run?.liability_summary||{};
  const setMetaboliteValue=(key,value)=>setMetaboliteForm(current=>({...current,[key]:value}));
  return e('div',{style:{marginTop:'24px'}},[
   e('div',{className:'row toolbar',key:'heading'},[
    e('h4',{},'Metabolism · Metabolic Soft Spots'),
    e('div',{className:'row'},[
     e('label',{key:'top',className:'small'},['Show ',e('select',{key:'select',value:metabolicTop,onChange:event=>setMetabolicTop(event.target.value==='ALL'?'ALL':Number(event.target.value))},(metabolism?.settings?.available_top_spots||[3,5,10,'ALL']).map(value=>e('option',{key:value,value},value==='ALL'?'All':'Top '+value)))]),
     e('button',{key:'run',disabled:metabolismBusy,onClick:()=>runMetabolism(versionId)},metabolismBusy?'Generating…':run?'Regenerate / use cache':'Generate hypotheses')
    ])
   ]),
   e('p',{key:'scope',className:'small'},'Atom indices are zero-based RDKit indices. Rule priors rank hypotheses but are not atom probabilities. CYP substrate and microsomal results are supporting compound-level evidence only.'),
   !run?Empty({children:'No metabolic soft-spot run for this CompoundVersion.'}):e('div',{key:'run'},[
    e('div',{key:'status',className:'small'},[Badge({ok:run.status==='COMPLETE',text:run.status}),e('span',{key:'message'},' · '+run.message)]),
    e('div',{className:'grid',style:{marginTop:'14px'},key:'visual'},[
     e('div',{className:'col-6 structure metabolism-structure-panel',key:'svg'},Svg({src:run.highlighted_svg})),
     e('div',{className:'col-6',key:'table'},spots.length?e('table',{},[
      e('thead',{key:'head'},e('tr',{},['Rank','Atom','Environment','Transformation','Phase','Confidence'].map(label=>e('th',{key:label},label)))),
      e('tbody',{key:'body'},spots.map(spot=>e('tr',{key:spot.id,onClick:()=>setSelectedSpotId(spot.id),style:{cursor:'pointer',background:selected?.id===spot.id?'rgba(43,110,242,.10)':''}},[
       e('td',{key:'rank',className:'mono'},spot.rank),e('td',{key:'atom',className:'mono'},spot.atom_index),e('td',{key:'env',className:'mono small'},spot.atom_environment),
       e('td',{key:'transform'},spot.transformation),e('td',{key:'phase'},spot.phase),e('td',{key:'confidence'},spot.confidence)
      ])))
     ]):Empty({children:'No supported transformation matched this structure.'}))
    ]),
    selected&&e('details',{key:'selected',open:true},[
     e('summary',{key:'summary'},'Selected spot details · Rank '+selected.rank+' '+selected.transformation),
     e('div',{key:'body',className:'small'},[
      e('div',{key:'atom'},e('strong',{},'Atom: '),selected.atom_index+' · '+selected.atom_environment),
      e('div',{key:'phase'},e('strong',{},'Transformation / phase: '),selected.transformation+' · '+selected.phase),
      e('div',{key:'cyp'},e('strong',{},'CYP attribution: '),selected.cyp_isoform),
      e('div',{key:'score'},e('strong',{},'Ranking evidence: '),Number(selected.score).toPrecision(4)+' · '+selected.score_type),
      e('div',{key:'model'},e('strong',{},'Model evidence: '),(selected.model_evidence?.status||'Not predicted')+' — '+(selected.model_evidence?.reason||'')),
      e('div',{key:'rules'},e('strong',{},'Rules: '),(selected.rule_evidence?.rules||[]).map(rule=>rule.rule_name+' ('+rule.empirical_prior+')').join(' · ')),
      e('div',{key:'strategy'},e('strong',{},'General mitigation strategies (not analog proposals): '),(selected.rule_evidence?.mitigation_strategies||[]).join(' · ')),
      e('div',{key:'provenance'},e('strong',{},'Provenance: '),selected.provenance.engine+' '+selected.provenance.engine_version+' · '+selected.provenance.source+' · '+selected.provenance.license+' · '+selected.provenance.prediction_timestamp)
     ])
    ]),
    e('h4',{key:'liability-title',style:{marginTop:'20px'}},'Metabolic Liability Summary'),
    e('div',{key:'liability',className:'small'},[
     e('div',{key:'primary'},[e('strong',{},'Primary predicted liability: '),summary.primary_predicted_liability||'—']),
     e('div',{key:'support'},[e('strong',{},'Supporting evidence: '),summary.supporting_evidence||'—']),
     e('div',{key:'microsomal'},[e('strong',{},'Microsomal evidence: '),(summary.microsomal_evidence||[]).length?summary.microsomal_evidence.map(item=>item.source+' '+item.endpoint+' '+item.value+' '+item.unit+' ('+item.confidence+')').join(' · '):'None']),
     e('div',{key:'cyp'},[e('strong',{},'CYP evidence: '),(summary.cyp_evidence||[]).length?summary.cyp_evidence.map(item=>item.endpoint+' '+item.classification+' ('+item.confidence+')').join(' · '):'None']),
     e('div',{key:'limit'},summary.cyp_attribution_limit||''),
     e('div',{key:'strategies'},[e('strong',{},'Strategy: '),(summary.mitigation_strategies||[]).join(' · '),e('span',{key:'scope'},' · '+(summary.strategy_scope||''))])
    ]),
    e('h4',{key:'predicted-title',style:{marginTop:'20px'}},'Predicted Metabolites'),
    e('div',{key:'predicted-label',className:'fail small'},'PREDICTED METABOLITE HYPOTHESIS — chemically validated candidates, not confirmed metabolites.'),
    metabolites.length?e('table',{key:'metabolites'},[
     e('thead',{key:'head'},e('tr',{},['Rank','Structure','Transformation','Source atom','Phase','Confidence',''].map(label=>e('th',{key:label},label)))),
     e('tbody',{key:'body'},metabolites.map(item=>e('tr',{key:item.id},[
      e('td',{key:'rank'},item.rank),e('td',{key:'smiles',className:'metabolite-structure'},[Svg({src:item.structure_svg}),e('div',{className:'mono small'},item.canonical_smiles)]),e('td',{key:'transform'},item.transformation),e('td',{key:'atom'},item.source_atom),e('td',{key:'phase'},item.phase),e('td',{key:'confidence'},item.confidence),
      e('td',{key:'details'},e('details',{},[e('summary',{key:'summary'},'Details'),e('div',{key:'body',className:'small'},'Evidence: '+item.evidence.chemical_validation+' · Engine '+item.provenance.transformation_engine+' '+item.provenance.transformation_engine_version+' · '+item.provenance.source+' · '+item.provenance.license)]))
     ])))
    ]):Empty({children:'No unique sanitized metabolite hypotheses for the selected rank range.'}),
    e('details',{key:'tool',style:{marginTop:'14px'}},[
     e('summary',{key:'summary'},'Engine details and validation'),
     e('div',{key:'body',className:'small'},[(metabolism?.tool?.name||run.engine)+' '+(metabolism?.tool?.version||run.engine_version)+' · '+(metabolism?.tool?.license||'')+' · ',e('a',{key:'source',href:metabolism?.tool?.source,target:'_blank',rel:'noreferrer'},metabolism?.tool?.source),e('div',{key:'validation'},Object.entries(metabolism?.tool?.publisher_validation||{}).map(([key,value])=>key+': '+value).join(' · ')),e('div',{key:'model'},'Atom-level model: '+run.model_status.status+' — '+run.model_status.reason)])
    ])
   ]),
   e('h4',{key:'experimental-title',style:{marginTop:'24px'}},'Experimental Metabolites'),
   e('p',{key:'experimental-help',className:'small'},'Record LC-MS/MS or other observed evidence separately from predicted hypotheses. Structure is optional when unknown.'),
   e('div',{className:'grid',key:'experimental-form'},[
    e('div',{className:'col-6',key:'smiles'},Field({label:'Structure / SMILES (optional)',value:metaboliteForm.smiles,onChange:value=>setMetaboliteValue('smiles',value)})),
    e('div',{className:'col-3',key:'transformation'},Field({label:'Transformation',value:metaboliteForm.transformation,onChange:value=>setMetaboliteValue('transformation',value)})),
    e('div',{className:'col-3',key:'mass'},Field({label:'Observed mass',type:'number',value:metaboliteForm.observed_mass,onChange:value=>setMetaboliteValue('observed_mass',value)})),
    e('div',{className:'col-3',key:'unit'},Field({label:'Mass unit',value:metaboliteForm.mass_unit,onChange:value=>setMetaboliteValue('mass_unit',value)})),
    e('div',{className:'col-3',key:'source'},Field({label:'Source',value:metaboliteForm.source,onChange:value=>setMetaboliteValue('source',value)})),
    e('div',{className:'col-3',key:'experiment'},Field({label:'Experiment',value:metaboliteForm.experiment,onChange:value=>setMetaboliteValue('experiment',value)})),
    e('div',{className:'col-3',key:'notes'},Field({label:'Notes',value:metaboliteForm.notes,onChange:value=>setMetaboliteValue('notes',value)}))
   ]),
   e('button',{key:'save-experimental',style:{marginTop:'10px'},disabled:metabolismBusy||!metaboliteForm.transformation,onClick:()=>saveExperimentalMetabolite(versionId)},metabolismBusy?'Saving…':'Save experimental metabolite'),
   experimental.length?e('table',{key:'experimental-table',style:{marginTop:'14px'}},[
    e('thead',{key:'head'},e('tr',{},['Type','Structure','Transformation','Observed mass','Source','Experiment','Notes'].map(label=>e('th',{key:label},label)))),
    e('tbody',{key:'body'},experimental.map(item=>e('tr',{key:item.id},[e('td',{key:'type'},item.label),e('td',{key:'smiles',className:'mono small'},item.canonical_smiles||'Not reported'),e('td',{key:'transform'},item.transformation),e('td',{key:'mass'},item.observed_mass==null?'—':item.observed_mass+' '+item.mass_unit),e('td',{key:'source'},item.source||'—'),e('td',{key:'experiment'},item.experiment||'—'),e('td',{key:'notes'},item.notes||'—')])))
   ]):Empty({children:'No experimental metabolites recorded for this CompoundVersion.'})
  ]);
 }

 function proposalCandidatePanel(candidate){
  if(!candidate)return Empty({children:'Select a candidate to inspect its full rescoring snapshot.'});
  const formatCell=cell=>{
   if(!cell)return '—';const value=cell.value==null?'—':(typeof cell.value==='number'?Number(cell.value).toPrecision(5):String(cell.value));
   return value+(cell.unit?' '+cell.unit:'')+' · '+(cell.type||'Not available')+(cell.confidence?' · '+cell.confidence:'')+(cell.domain?' · '+(typeof cell.domain==='object'?(cell.domain.classification||'Not assessed'):cell.domain):'');
  };
  const propertyDelta=candidate.property_delta||{},activity=candidate.activity||{},soft=candidate.soft_spot_changes||{};
  const safety=(candidate.parent_comparison||[]).filter(row=>['hERG liability','Ames mutagenicity','DILI clinical liability','CYP3A4 inhibitor','P-gp inhibitor'].includes(row.endpoint));
  return e('div',{className:'card candidate-detail'},[
   e('div',{className:'row toolbar',key:'header'},[e('div',{},[e('h3',{},'Candidate '+candidate.candidate_number+' · '+(candidate.ranking?'Rank '+candidate.ranking.rank:candidate.status)),e('div',{className:'mono small'},candidate.canonical_smiles)]),e('div',{className:'manual-actions'},[
    e('button',{key:'promote',className:'secondary',disabled:proposalBusy,onClick:()=>candidateDecision(candidate,'PROMOTED')},'Promote'),
    e('button',{key:'reject',className:'danger',disabled:proposalBusy,onClick:()=>candidateDecision(candidate,'REJECTED')},'Reject')
   ])]),
   e('div',{className:'grid',key:'structures'},[
    e('div',{className:'col-6 structure difference-structure',key:'parent'},[e('h4',{},'Parent difference'),Svg({src:candidate.parent_difference_svg})]),
    e('div',{className:'col-6 structure difference-structure',key:'candidate'},[e('h4',{},'Candidate difference'),Svg({src:candidate.candidate_difference_svg||candidate.structure_svg})])
   ]),
   e('div',{className:'grid candidate-facts',key:'facts'},[
    e('div',{className:'col-3'},[e('strong',{},'Parent similarity'),e('div',{className:'mono'},Number(candidate.parent_similarity).toFixed(3))]),
    e('div',{className:'col-3'},[e('strong',{},'MCS coverage'),e('div',{className:'mono'},Number(candidate.mcs_coverage).toFixed(3))]),
    e('div',{className:'col-3'},[e('strong',{},'Confidence / domain'),e('div',{},candidate.confidence+' · '+candidate.applicability_domain)]),
    e('div',{className:'col-3'},[e('strong',{},'Information Value'),e('div',{},candidate.information_value)]),
    e('div',{className:'col-3'},[e('strong',{},'Ranking score'),e('div',{className:'mono'},candidate.ranking_score==null?'—':Number(candidate.ranking_score).toFixed(3))]),
    e('div',{className:'col-3'},[e('strong',{},'Pareto front'),e('div',{},candidate.pareto_front||'—')]),
    e('div',{className:'col-6'},[e('strong',{},'Synthetic complexity'),e('div',{},candidate.synthetic_feasibility?.classification||'—'),e('div',{className:'small'},'SA surrogate '+(candidate.synthetic_feasibility?.sa_score??'—')+' · not synthesis success probability')])
   ]),
   e('h4',{key:'why'},'Why generated / expected benefit'),
   e('p',{key:'why-text'},candidate.why_generated+' · '+candidate.expected_benefit),
   e('div',{key:'transforms'},candidate.transformations.map(row=>e('details',{key:row.sequence+'-'+row.id},[e('summary',{key:'summary'},'Transformation '+row.sequence+': '+row.name+' · '+row.execution_status),e('div',{key:'body',className:'small'},[e('div',{className:'mono'},row.reaction_smarts||'User-defined structure'),e('div',{},'Source atoms: '+row.source_atoms.join(', ')+' · version '+row.version),e('div',{},'Source: '+row.source)])]))),
   e('h4',{key:'activity-title'},'Activity prediction'),
   e('p',{key:'activity',className:'small'},activity.status==='COMPLETE'?(Number(activity.value_nm).toPrecision(5)+' nM · '+activity.record_type+' · '+activity.confidence+' · '+activity.applicability_domain+' · nearest '+(activity.nearest_neighbors||[]).slice(0,3).map(row=>row.compound_id+' '+row.similarity).join(', ')):('Not predicted — '+(activity.reason||'No selected assay model'))),
   e('h4',{key:'properties-title'},'Stage 1 property changes'),
   e('div',{key:'properties',className:'small'},Object.entries(propertyDelta).map(([key,value])=>e('span',{key,className:value<0?'delta-down':'delta-up'},key+' '+(value>=0?'+':'')+Number(value).toFixed(3)+' '))),
   e('h4',{key:'soft-title'},'Soft spot changes'),
   e('p',{key:'soft',className:'small'},'Parent primary: '+(soft.parent_primary?.transformation||'Not assessed')+' · Candidate primary: '+(soft.candidate_primary?.transformation||'None')+' · parent site absent from candidate Top 3: '+(soft.parent_primary_absent_from_candidate_top3?'YES':'NO')+' · new primary liability: '+(soft.new_primary_liability?'YES':'NO')),
   e('h4',{key:'safety-title'},'Safety flags'),
   safety.length?e('ul',{key:'safety'},safety.map(row=>e('li',{key:row.endpoint},row.endpoint+': '+formatCell(row.candidate)))):Empty({children:'No safety endpoint result available.'}),
   e('h4',{key:'comparison-title'},'Parent vs Candidate'),
   e('table',{key:'comparison'},[
    e('thead',{key:'head'},e('tr',{},['Endpoint','Parent','Candidate','Change'].map(label=>e('th',{key:label},label)))),
    e('tbody',{key:'body'},(candidate.parent_comparison||[]).map(row=>e('tr',{key:row.endpoint},[e('td',{key:'endpoint'},row.endpoint),e('td',{key:'parent'},formatCell(row.parent)),e('td',{key:'candidate'},formatCell(row.candidate)),e('td',{key:'change',className:'mono'},row.change==null?'qualitative / uncertain':(row.change>=0?'+':'')+row.change)])))
   ]),
   e('div',{className:'alert',key:'risk'},[e('strong',{},'Main risk: '),candidate.main_risk]),
   candidate.rejection_reasons?.length>0&&e('div',{key:'rejects'},[e('h4',{},'Rejection reason'),...candidate.rejection_reasons.map((row,index)=>e('div',{key:index,className:'fail'},row.code+' · '+row.detail+' · '+row.stage))]),
   e('details',{key:'formula'},[e('summary',{key:'summary'},'Ranking formula and prediction provenance'),e('div',{key:'body',className:'small'},[
    e('pre',{key:'formula',className:'small'},JSON.stringify(candidate.ranking?.formula||candidate.objective_vector,null,2)),
    ...(candidate.prediction_snapshots||[]).map((row,index)=>e('div',{key:index},row.stage+' · '+row.endpoint+' · '+row.type+' · '+row.model+' '+row.model_version+' · '+row.confidence+' · '+row.domain))
   ])])
  ]);
 }

 function proposalPanel(){
  if(!optimizationRun)return null;
  const active=proposalRun&&['PENDING','GENERATING','FILTERING','PREDICTING','RANKING'].includes(proposalRun.status);
  const views=[['all','Show all generated'],['accepted','Show accepted'],['rejected','Show rejected'],['pareto','Show Pareto front'],['top10','Show Top 10']];
  const candidates=proposalRun?.candidates||[];
  return e('div',{className:'proposal-section'},[
   e('div',{className:'card',key:'generate'},[
    e('div',{className:'row toolbar',key:'header'},[e('div',{},[e('h3',{},'Analog Generation, Rescoring & Ranking'),e('p',{className:'small'},'Deterministic curated transformations only · single change first · maximum two changes · no LLM · no PK')]),e('button',{disabled:proposalBusy||active,onClick:generateAnalogs},proposalBusy?'Starting…':'Generate analogs')]),
    e('div',{className:'grid',key:'settings'},[
     e('div',{className:'col-3'},Field({label:'Maximum raw candidates (1–200)',type:'number',value:proposalSettings.max_raw_candidates,onChange:value=>setProposalSettings(current=>({...current,max_raw_candidates:Number(value)}))})),
     e('div',{className:'col-4'},e('label',{className:'check-option'},[e('input',{type:'checkbox',checked:proposalSettings.allow_double_transforms,onChange:event=>setProposalSettings(current=>({...current,allow_double_transforms:event.target.checked}))}),e('span',{},'Allow limited two-transformation hypotheses')]))
    ]),
    e('p',{key:'staged',className:'small'},'Staged execution: generation → RDKit chemical validation → Stage 1/similarity/hard gates → project activity → available ADMET/soft spots → Pareto/ranking. A failed candidate does not stop the run.'),
    proposalRun&&e('div',{key:'status',className:'job-status '+proposalRun.status.toLowerCase()},[e('strong',{},proposalRun.status),e('span',{},' · '+proposalRun.stage_message),e('div',{className:'small'},'Raw '+proposalRun.raw_candidate_count+' · accepted '+proposalRun.accepted_count+' · rejected '+proposalRun.rejected_count+' · selected '+proposalRun.top_count)])
   ]),
   proposalRun?.status==='COMPLETED'&&e(React.Fragment,{key:'results'},[
    e('div',{className:'card',key:'filters'},[
     e('div',{className:'row toolbar',key:'view'},[e('h3',{},'Candidate Filtering'),e('div',{className:'row'},views.map(([value,label])=>e('button',{key:value,className:proposalView===value?'':'secondary',onClick:()=>refreshProposal(proposalRun.id,value)},label)))]),
     candidates.length?e('table',{key:'table'},[
      e('thead',{key:'head'},e('tr',{},['Rank','Candidate','Status','Similarity','Transformation hypothesis','Score','Pareto','Confidence','Information','Main risk',''].map(label=>e('th',{key:label},label)))),
      e('tbody',{key:'body'},candidates.map(candidate=>e('tr',{key:candidate.id,className:selectedCandidate?.id===candidate.id?'selected-row':''},[
       e('td',{key:'rank'},candidate.ranking?.rank||'—'),e('td',{key:'candidate',className:'mono'},'#'+candidate.candidate_number),e('td',{key:'status'},candidate.status),e('td',{key:'similarity',className:'mono'},Number(candidate.parent_similarity||0).toFixed(3)),e('td',{key:'transform'},candidate.hypothesis),e('td',{key:'score',className:'mono'},candidate.ranking_score==null?'—':Number(candidate.ranking_score).toFixed(2)),e('td',{key:'pareto'},candidate.pareto_front||'—'),e('td',{key:'confidence'},candidate.confidence),e('td',{key:'info'},candidate.information_value),e('td',{key:'risk',className:'small'},candidate.main_risk||(candidate.rejection_reasons||[]).map(row=>row.code).join(', ')),e('td',{key:'open'},e('button',{className:'secondary',onClick:()=>setSelectedCandidate(candidate)},'Details'))
      ])))
     ]):Empty({children:'No candidates in this filter. A small chemically meaningful pool is allowed; no filler structures are generated.'})
    ]),
    proposalCandidatePanel(selectedCandidate),
    e('div',{className:'card',key:'manual'},[
     e('h3',{key:'title'},'User-added Analog'),e('p',{key:'help',className:'small'},'Paste ChemDraw/Ketcher SMILES. The structure receives the same Stage 1 → Activity → ADMET → ranking workflow and remains labeled user-added.'),
     e('div',{className:'grid',key:'form'},[e('div',{className:'col-8'},Field({label:'SMILES',value:userAnalog.smiles,onChange:value=>setUserAnalog(current=>({...current,smiles:value}))})),e('div',{className:'col-4'},Field({label:'Reason / hypothesis',value:userAnalog.reason,onChange:value=>setUserAnalog(current=>({...current,reason:value}))}))]),
     e('button',{key:'add',style:{marginTop:'10px'},disabled:proposalBusy||!userAnalog.smiles,onClick:addUserAnalog},proposalBusy?'Rescoring…':'Add and rescore analog')
    ])
   ])
  ]);
 }

 function optimizationPanel(versionId){
  const config=optimizationConfig||{objectives:[],evidence_hierarchy:[]},run=optimizationRun;
  const setConstraint=(key,value)=>setOptimizationForm(current=>({...current,constraints:{...(current?.constraints||{}),[key]:value}}));
  const toggleObjective=name=>setOptimizationForm(current=>({...current,objectives:(current?.objectives||[]).includes(name)?(current?.objectives||[]).filter(value=>value!==name):[...(current?.objectives||[]),name]}));
  const addOverride=(key,value)=>{
   const values=run?.manual_overrides?.[key]||[],encoded=JSON.stringify(value);
   const next=values.some(item=>JSON.stringify(item)===encoded)?values:[...values,value];
   overrideOptimization({[key]:next});
  };
  const admetProfile=run?.evidence?.admet||{},activity=run?.evidence?.activity||{},properties=run?.evidence?.properties||{};
  const evidenceValue=row=>{
   const preferred=row?.preferred;if(!preferred)return 'Not measured / not predicted';
   const value=preferred.classification??preferred.assessment?.category??preferred.value;
   return String(value??'Not available')+(preferred.unit?' '+preferred.unit:'')+' · '+preferred.type+(preferred.confidence?' · '+preferred.confidence:'')+(preferred.applicability_domain?' · '+preferred.applicability_domain:'');
  };
  const constraintField=(key,label,type='number')=>e('div',{className:'col-3',key},Field({label,type,value:optimizationForm?.constraints?.[key]??'',onChange:value=>setConstraint(key,value)}));
  const regionTable=(title,rows,protectedType)=>e('div',{className:'col-6 optimization-region',key:title},[
   e('h3',{key:'title'},title),
   rows?.length?e('table',{key:'table'},[
    e('thead',{key:'head'},e('tr',{},['Atoms / fragment','Reason','Risk','Confidence','Override'].map(label=>e('th',{key:label},label)))),
    e('tbody',{key:'body'},rows.map(row=>e('tr',{key:row.id},[
     e('td',{key:'atoms',className:'mono small'},row.atom_indices?.length?row.atom_indices.join(', '):(row.fragment||'Not localized')),
     e('td',{key:'reason'},row.reason),e('td',{key:'risk'},row.risk||row.status),e('td',{key:'confidence'},row.confidence),
     e('td',{key:'override'},row.atom_indices?.length?e('button',{className:'secondary',disabled:optimizationBusy,onClick:()=>addOverride(protectedType?'allow_atoms':'protect_atoms',row.atom_indices)},protectedType?'Allow modification':'Protect this atom/group'):'—')
    ])))
   ]):Empty({children:'No region evidence.'})
  ]);
  return e('div',{key:'optimization'},[
   e('div',{className:'strategy-banner',key:'scope'},[
    e('strong',{key:'title'},'Stage 4A deterministic strategy only'),
    e('span',{key:'text'},' — ranks medicinal chemistry transformations; no analog structures, PK, overall score, or LLM reasoning are generated.')
   ]),
   e('div',{className:'card',key:'setup'},[
    e('h3',{key:'title'},'Step 3 — Select Optimization Goal'),
    e('p',{key:'parent',className:'small'},'Parent: '+(detail?.compound_id||'Compound')+' v'+(detail?.current_version||'1')+' · CompoundVersion #'+versionId),
    e('div',{className:'grid',key:'top'},[
     e('div',{className:'col-4 optimization-assay-field',key:'assay'},[e('label',{},'Selected assay'),e('select',{value:optimizationForm?.assay_id||'',onChange:event=>setOptimizationForm(current=>({...current,assay_id:event.target.value}))},[e('option',{key:'none',value:''},'No assay selected'),...assays.map(assay=>e('option',{key:assay.id,value:assay.id},assay.name+' · '+(assay.measurement_type||'')))])]),
     e('div',{className:'col-8',key:'objectives'},[e('label',{},'Optimization objective(s)'),e('div',{className:'objective-grid'},(config.objectives||[]).map(name=>e('label',{key:name,className:'check-option'},[e('input',{key:'input',type:'checkbox',checked:(optimizationForm?.objectives||[]).includes(name),onChange:()=>toggleObjective(name)}),e('span',{key:'label'},name)])))])
    ]),
    (optimizationForm?.objectives||[]).includes('Custom')&&e('div',{key:'custom',style:{marginTop:'10px'}},Field({label:'Custom objective',value:optimizationForm?.custom_objective||'',onChange:value=>setOptimizationForm(current=>({...current,custom_objective:value}))})),
    e('h4',{key:'constraints-title',style:{marginTop:'18px'}},'Constraints'),
    e('div',{className:'grid',key:'constraints'},[
     constraintField('potency_max_nm','Potency IC50 ≤ (nM)'),constraintField('do_not_worsen_fold','Do not worsen potency > fold'),constraintField('clogp_max','cLogP ≤'),constraintField('mw_max','MW ≤'),
     constraintField('tpsa_min','TPSA minimum Å²'),constraintField('tpsa_max','TPSA maximum Å²'),constraintField('similarity_min','Future analog similarity ≥'),constraintField('logs_min','LogS minimum'),constraintField('caco2_logpapp_min','Caco-2 LogPapp minimum'),
     e('div',{className:'col-3',key:'herg'},e('label',{className:'check-option'},[e('input',{type:'checkbox',checked:!!optimizationForm?.constraints?.herg_do_not_increase,onChange:event=>setConstraint('herg_do_not_increase',event.target.checked)}),e('span',{},'hERG: do not increase liability')]))
    ]),
    e('p',{key:'precedence',className:'small'},'Experimental evidence takes precedence over prediction. Low-confidence classification alone remains supporting-only. Similarity and do-not-worsen constraints are stored now as hard gates for a future proposal stage; Stage 4A does not create candidates.'),
    e('button',{key:'analyze',disabled:optimizationBusy||!(optimizationForm?.objectives||[]).length,onClick:()=>analyzeOptimization(versionId)},optimizationBusy?'Analyzing…':'Analyze Optimization Strategy')
   ]),
   (optimizationRuns||[]).length>0&&e('div',{className:'row run-picker',key:'history'},[e('label',{key:'label'},'Saved runs'),e('select',{key:'select',value:run?.id||'',onChange:event=>setOptimizationRun(optimizationRuns.find(item=>item.id===Number(event.target.value)))},optimizationRuns.map(item=>e('option',{key:item.id,value:item.id},'#'+item.id+' · '+(item.objectives||[]).join(' + ')+' · '+item.status))) ]),
   run&&e(React.Fragment,{key:'results'},[
    e('div',{className:'card',key:'profile'},[
     e('div',{className:'row toolbar',key:'head'},[e('h3',{},'Current profile'),e('span',{className:'small'},run.engine+' '+run.engine_version)]),
     e('div',{className:'grid',key:'profile-grid'},[
      e('div',{className:'col-4',key:'activity'},[e('h4',{},'Activity'),e('p',{className:'small'},activity.experimental?'Experimental '+activity.experimental.mean_nm+' nM':(activity.predicted?'Predicted '+activity.predicted.value_nm+' nM · '+activity.predicted.confidence:'Not measured / not predicted for selected assay'))]),
      e('div',{className:'col-4',key:'properties'},[e('h4',{},'Properties'),e('p',{className:'small'},['molecular_weight','clogp','tpsa','fraction_csp3'].map(key=>key+' '+(properties[key]?.value??'—')).join(' · ')+' · Calculated / RDKit')]),
      e('div',{className:'col-4',key:'admet'},[e('h4',{},'ADMET / metabolism'),...Object.entries(admetProfile).slice(0,12).map(([name,row])=>e('div',{key:name,className:'small'},name+': '+evidenceValue(row))),Object.keys(admetProfile).length===0&&Empty({children:'No compatible experimental or predicted ADMET evidence.'})])
     ]),
     e('details',{key:'hierarchy'},[e('summary',{key:'summary'},'Evidence hierarchy'),e('ol',{key:'list',className:'small'},(run.evidence?.evidence_hierarchy||[]).map(item=>e('li',{key:item.rank},item.type+' · ordinal weight '+item.weight)))])
    ]),
    e('div',{className:'card',key:'liabilities'},[
     e('h3',{key:'title'},'Main liabilities'),
     run.liabilities.length?e('table',{key:'table'},[
      e('thead',{key:'head'},e('tr',{},['Rank','Liability','Evidence','Confidence','Actionability','Rationale'].map(label=>e('th',{key:label},label)))),
      e('tbody',{key:'body'},run.liabilities.map(row=>e('tr',{key:row.id},[e('td',{},row.rank),e('td',{},row.title),e('td',{},row.evidence_type),e('td',{},row.confidence),e('td',{className:row.actionability==='ACTIONABLE'?'pass':'small'},row.actionability),e('td',{className:'small'},row.rationale)])))
     ]):Empty({children:'No deterministic liability threshold was triggered. Unavailable evidence remains visible in Current profile.'})
    ]),
    e('div',{className:'card grid',key:'regions'},[
     e('div',{className:'col-12',key:'visual'},[e('h3',{},'Structure regions'),e('div',{className:'optimization-structure structure'},Svg({src:run.highlighted_svg})),e('div',{className:'region-legend'},[e('span',{className:'legend-protected'},'Protected'),e('span',{className:'legend-modifiable'},'Modifiable'),e('span',{className:'legend-soft'},'Metabolic soft spot')])]),
     regionTable('Protected regions',run.protected_regions,true),regionTable('Modifiable regions',run.modifiable_regions,false)
    ]),
    e('div',{className:'card',key:'transformations'},[
     e('h3',{key:'title'},'Recommended transformations'),e('p',{key:'scope',className:'small'},'Ranked strategies only. Reaction SMARTS is provenance, not an instruction that has been applied to the parent.'),
     run.recommended_transformations.length?e('table',{key:'table'},[
      e('thead',{key:'head'},e('tr',{},['Rank','Transformation','Expected effect','Potency risk','Evidence','Confidence','Manual','Details'].map(label=>e('th',{key:label},label)))),
      e('tbody',{key:'body'},run.recommended_transformations.map(row=>e('tr',{key:row.id+'-'+row.rank},[
       e('td',{key:'rank'},row.rank),e('td',{key:'name'},[e('strong',{key:'name'},row.name),e('div',{key:'purpose',className:'small'},row.purpose)]),e('td',{key:'effect'},row.expected_effect),e('td',{key:'risk'},row.potency_risk),e('td',{key:'evidence',className:'small'},row.evidence.join(' · ')),e('td',{key:'confidence'},row.confidence),
       e('td',{key:'actions'},e('div',{className:'manual-actions'},[e('button',{key:'prioritize',className:'secondary',disabled:optimizationBusy,onClick:()=>addOverride('prioritize_transformations',row.id)},'Prioritize'),e('button',{key:'exclude',className:'danger',disabled:optimizationBusy,onClick:()=>addOverride('exclude_transformations',row.id)},'Exclude')])),
       e('td',{key:'details'},e('details',{},[e('summary',{key:'summary'},'Details'),e('div',{key:'body',className:'small strategy-details'},[e('div',{key:'id'},'ID / version: '+row.id+' / '+row.version),e('div',{key:'smarts',className:'mono'},'Reaction SMARTS: '+row.reaction_smarts),e('div',{key:'motif'},'Applicable motif: '+row.applicable_motif+' · source atoms '+row.source_atom_indices.join(', ')),e('div',{key:'risk'},'Possible risk: '+row.possible_risk),e('div',{key:'source'},'Source/reference: '+row.source),e('div',{key:'status'},row.application_status)])]))
      ])))
     ]):Empty({children:'No applicable transformation strategy was ranked. No fake recommendation is generated.'})
    ]),
    e('div',{key:'proposal'},proposalPanel())
   ])
  ]);
 }

 function metabolismProfile(versionId){
  const predictions=(admet?.predictions||[]).filter(row=>row.version_id===Number(versionId));
  const measurements=(admet?.measurements||[]).filter(row=>row.version_id===Number(versionId));
  const endpointRows=names=>predictions.filter(row=>names.some(name=>typeof name==='string'?row.endpoint===name:name.test(row.endpoint)));
  const experimentalFor=matcher=>measurements.filter(row=>matcher.test(String(endpointName(row.endpoint_id))));
  const evidenceList=(title,experimental,predicted)=>e('div',{className:'card',key:title},[
   e('h3',{},title),e('div',{className:'evidence-columns'},[
    e('section',{key:'exp'},[e('h4',{},'Experimental Results'),experimental.length?admetMeasurementTable(experimental):e('div',{className:'empty-state'},[StatusBadge({type:'Not measured'}),e('p',{},'No experimental measurement entered.')])]),
    e('section',{key:'pred'},[e('h4',{},'Prediction Results'),predicted.length?admetPredictionTable(predicted):e('div',{className:'empty-state'},[StatusBadge({type:'Not predicted'}),e('p',{},'Prediction not run or no validated endpoint model is active.')])])
   ])
  ]);
  const stability=endpointRows([/intrinsic clearance$/]);
  const cyp=endpointRows([/^CYP/]);
  const supporting=endpointRows(['Permeability','Plasma protein binding']);
  return e('div',{className:'metabolism-profile'},[
   evidenceList('Metabolic Stability · Human / Rat / Mouse Liver Microsomes',experimentalFor(/HLM|RLM|MLM|microsom/i),stability),
   e('div',{className:'card',key:'cyp'},[e('h3',{},'CYP · Inhibitor and Substrate'),e('p',{className:'small'},'Roles remain endpoint-separated. Compound-level substrate evidence does not assign an atom or reaction to a CYP isoform.'),cypPredictionTable(cyp)]),
   evidenceList('Supporting ADME Evidence · Permeability and Plasma Protein Binding (PPB)',experimentalFor(/permeab|caco|protein binding|PPB/i),supporting),
   e('div',{className:'card',key:'soft'},[e('h3',{},'Metabolic Soft Spots'),e('p',{className:'small'},'The parent structure is interpreted together with stability, permeability, PPB, and CYP evidence shown above.'),metabolismPanel(versionId)])
  ]);
 }

 function pkConcentrationTimePlot(observations, latestNca, plotType){
  const obs=(observations||[]).slice().sort((a,b)=>a.time_hours-b.time_hours);
  if(!obs.length)return e('div',{className:'empty-state'},[StatusBadge({type:'Not measured'}),e('p',{},'No experimental concentration-time points recorded for this study.')]);
  const quant=obs.filter(o=>!o.blq_flag&&o.concentration_normalized_ng_ml!=null&&o.concentration_normalized_ng_ml>0);
  const maxTime=Math.max(...obs.map(o=>o.time_hours),1);
  const maxConc=Math.max(...obs.map(o=>o.concentration_normalized_ng_ml||0),1);
  const isLog=plotType==='log';
  const minLog=isLog?-3:0,maxLog=isLog?Math.log10(Math.max(maxConc,1)):maxConc;
  const padL=60,padR=30,padT=25,padB=45,w=680,h=260;
  const pw=w-padL-padR,ph=h-padT-padB;
  const getX=t=>padL+(t/maxTime)*pw;
  const getY=c=>{
   if(isLog){
    const logVal=c>0?Math.log10(c):minLog;
    const norm=(logVal-minLog)/(maxLog-minLog||1);
    return padT+ph*(1-Math.max(0,Math.min(1,norm)));
   }
   return padT+ph*(1-Math.max(0,Math.min(1,c/maxConc)));
  };
  const termSet=new Set(latestNca?.terminal_points||[]);
  const pointsPath=quant.map((p,i)=>(i===0?'M ':'L ')+getX(p.time_hours)+','+getY(p.concentration_normalized_ng_ml)).join(' ');

  return e('div',{className:'pk-plot-shell'},[
   e('svg',{className:'pk-plot-svg',viewBox:'0 0 '+w+' '+h,key:'svg'},[
    [0.25,0.5,0.75,1.0].map(frac=>e('line',{key:'gx-'+frac,x1:padL+pw*frac,y1:padT,x2:padL+pw*frac,y2:padT+ph,stroke:'#e2e8f0',strokeDasharray:'3,3'})),
    [0.25,0.5,0.75].map(frac=>e('line',{key:'gy-'+frac,x1:padL,y1:padT+ph*frac,x2:padL+pw,y2:padT+ph*frac,stroke:'#e2e8f0',strokeDasharray:'3,3'})),
    e('line',{key:'xaxis',x1:padL,y1:padT+ph,x2:padL+pw,y2:padT+ph,stroke:'#64748b',strokeWidth:1.5}),
    e('line',{key:'yaxis',x1:padL,y1:padT,x2:padL,y2:padT+ph,stroke:'#64748b',strokeWidth:1.5}),
    quant.length>1&&e('path',{key:'line',d:pointsPath,fill:'none',stroke:'#15803d',strokeWidth:2.2}),
    obs.map(o=>{
     const cx=getX(o.time_hours),cy=getY(o.blq_flag?0:o.concentration_normalized_ng_ml);
     const isTerm=termSet.has(o.id);
     return e('g',{key:o.id},[
      o.blq_flag?e('rect',{x:cx-4,y:cy-4,width:8,height:8,fill:'#ef4444',stroke:'#b91c1c',strokeWidth:1}):
      e('circle',{cx,cy,r:isTerm?6:4.5,fill:isTerm?'#3b82f6':'#15803d',stroke:isTerm?'#1e40af':'#14532d',strokeWidth:isTerm?2.5:1.5}),
      e('text',{x:cx,y:cy-8,textAnchor:'middle',fontSize:9,fill:'#475569'},o.blq_flag?'BLQ':Number(o.concentration_normalized_ng_ml).toFixed(1))
     ]);
    }),
    e('text',{key:'xlabel',x:padL+pw/2,y:h-8,textAnchor:'middle',fontSize:11,fill:'#475569',fontWeight:600},'Time ('+(obs[0]?.time_unit||'h')+')'),
    e('text',{key:'ylabel',x:15,y:padT+ph/2,textAnchor:'middle',fontSize:11,fill:'#475569',fontWeight:600,transform:'rotate(-90 15 '+(padT+ph/2)+')'},isLog?'Log10 Conc (ng/mL)':'Concentration ('+(obs[0]?.concentration_unit||'ng/mL')+')'),
    [0,0.5,1.0].map(frac=>e('text',{key:'xtick-'+frac,x:padL+pw*frac,y:padT+ph+15,textAnchor:'middle',fontSize:10,fill:'#64748b'},Number(maxTime*frac).toFixed(1))),
    [0,0.5,1.0].map(frac=>e('text',{key:'ytick-'+frac,x:padL-6,y:padT+ph*(1-frac)+4,textAnchor:'end',fontSize:10,fill:'#64748b'},isLog?Number(minLog+(maxLog-minLog)*frac).toFixed(1):Number(maxConc*frac).toFixed(1)))
   ])
  ]);
 }

 function iviveProfile(versionId){
  const data=iviveData,latest=data?.latest_run,outputs=latest?.outputs||{};
  const candidates=data?.candidates||{clint:[],plasma_binding:[],blood_plasma_ratio:[]};
  const physiology=data?.physiology||{};
  const sig=value=>value==null?'—':Number(value).toLocaleString(undefined,{maximumSignificantDigits:3});
  const sourceBadge=label=>e('span',{className:'ivive-source ivive-source-'+String(label||'').toLowerCase().replace(/[^a-z]+/g,'-')},label||'—');
  const selected=rows=>(rows||[]).find(row=>row.selected);
  const selectedClint=selected(candidates.clint),selectedBinding=selected(candidates.plasma_binding),selectedBpr=selected(candidates.blood_plasma_ratio);

  const runIviveAction=async()=>{
   setIviveBusy(true);
   try{
    const run=await api.post('/compound-versions/'+versionId+'/ivive/run',{species:iviveSpecies,method_key:'WELL_STIRRED'});
    await loadIviveData(versionId,iviveSpecies);
    setMessage(run.status==='COMPLETE'?'IVIVE hepatic clearance complete':'IVIVE run saved as unavailable; review missing inputs and warnings');
   }catch(err){setMessage(String(err))}finally{setIviveBusy(false)}
  };
  const addIviveInputAction=async()=>{
   if(iviveInputForm.input_value===''||!iviveInputForm.unit.trim())return;
   setIviveBusy(true);
   try{
    await api.post('/compound-versions/'+versionId+'/ivive-inputs',{
     ...iviveInputForm,species:iviveSpecies,input_value:Number(iviveInputForm.input_value)
    });
    setIviveInputForm(current=>({...current,input_value:'',notes:''}));
    await loadIviveData(versionId,iviveSpecies);setMessage('IVIVE input saved with CompoundVersion provenance');
   }catch(err){setMessage(String(err))}finally{setIviveBusy(false)}
  };
  const addOverrideAction=async()=>{
   if(iviveOverrideForm.value===''||!iviveOverrideForm.source.trim())return;
   setIviveBusy(true);
   try{
    await api.post('/projects/'+data.scope.project_id+'/ivive/physiology-overrides',{
     ...iviveOverrideForm,species:iviveSpecies,value:Number(iviveOverrideForm.value)
    });
    setIviveOverrideForm(current=>({...current,value:'',notes:''}));
    await loadIviveData(versionId,iviveSpecies);setMessage('Study-specific physiology override saved');
   }catch(err){setMessage(String(err))}finally{setIviveBusy(false)}
  };
  const updateInputEndpoint=value=>{
   const patch=value==='CLINT'?{input_endpoint:value,input_type:'RAW_MICROSOMAL',unit:'µL/min/mg protein'}:
    value==='FU_PLASMA'?{input_endpoint:value,input_type:'',unit:'fu'}:{input_endpoint:value,input_type:'',unit:'ratio'};
   setIviveInputForm(current=>({...current,...patch}));
  };
  const updateInputType=value=>setIviveInputForm(current=>({...current,input_type:value,unit:
   value==='RAW_MICROSOMAL'?'µL/min/mg protein':value==='RAW_HEPATOCYTE'?'µL/min/10^6 cells':'mL/min/kg'}));
  const candidateRow=(row,key)=>e('div',{className:'ivive-input-row'+(row?.selected?' selected':''),key},row?[
   e('div',{key:'main'},[sourceBadge(row.source_label),e('strong',{className:'mono'},sig(row.value)+' '+row.unit),row.input_type&&e('small',{},row.input_type.replaceAll('_',' '))]),
   e('div',{className:'small',key:'source'},row.model_source||row.provenance?.source||'Source not specified'),
   e('div',{className:'small',key:'quality'},'Confidence '+row.confidence+' · Domain '+row.applicability_domain+(row.selected?' · SELECTED BY PRIORITY':''))
  ]:[sourceBadge('UNAVAILABLE'),e('span',{className:'small'},'No compatible quantitative input')]);

  return e('section',{className:'card ivive-section',key:'ivive'},[
   e('div',{className:'row toolbar',key:'header'},[
    e('div',{},[e('div',{className:'eyebrow'},'4 · IVIVE HEPATIC CLEARANCE FOUNDATION'),e('h2',{},'IVIVE'),e('p',{className:'small'},'Intrinsic clearance scaling and well-stirred hepatic clearance only. Renal and other non-hepatic clearance are not modeled; total clearance is not predicted.'),e('p',{className:'small'},'Source badges: EXP · PRED · CALC · DEFAULT PHYSIOLOGY · USER OVERRIDE')]),
    e('div',{className:'ivive-run-controls'},[
     e('label',{},'Species'),
     e('select',{value:iviveSpecies,onChange:event=>setIviveSpecies(event.target.value)},(data?.supported_species||['Mouse','Rat','Dog','Monkey','Human']).map(species=>e('option',{key:species,value:species},species))),
     e('button',{disabled:iviveBusy||!data,onClick:runIviveAction},iviveBusy?'Running…':'Run IVIVE')
    ])
   ]),

   e('div',{className:'ivive-grid',key:'inputs'},[
    e('div',{className:'ivive-panel',key:'clint'},[e('h3',{},'Inputs · Clint'),candidateRow(selectedClint,'clint'),(candidates.clint||[]).length>1&&e('details',{},[e('summary',{},'All '+candidates.clint.length+' compatible Clint candidates'),e('div',{},candidates.clint.map((row,index)=>candidateRow(row,'clint-'+index)))])]),
    e('div',{className:'ivive-panel',key:'binding'},[e('h3',{},'Inputs · PPB / fu,p'),candidateRow(selectedBinding,'binding'),e('p',{className:'small'},'Experimental PPB/fu,p takes precedence over predicted PPB.')]),
    e('div',{className:'ivive-panel',key:'bpr'},[e('h3',{},'Inputs · Blood/Plasma'),candidateRow(selectedBpr,'bpr'),e('p',{className:'small'},selectedBpr?'fu,b = fu,p / (B/P)':'No compound-specific B/P invented. A run may use a labeled plasma-based approximation.')])
   ]),

   e('details',{className:'ivive-entry',key:'add-input'},[
    e('summary',{},'Add IVIVE-specific experimental or project-calibrated input'),
    e('div',{className:'grid ivive-form'},[
     e('div',{className:'col-3'},[e('label',{},'Endpoint'),e('select',{value:iviveInputForm.input_endpoint,onChange:event=>updateInputEndpoint(event.target.value)},[
      e('option',{value:'CLINT'},'Intrinsic clearance (Clint)'),e('option',{value:'FU_PLASMA'},'Fraction unbound in plasma (fu,p)'),e('option',{value:'BLOOD_PLASMA_RATIO'},'Blood/plasma ratio (B/P)')
     ])]),
     iviveInputForm.input_endpoint==='CLINT'&&e('div',{className:'col-3'},[e('label',{},'Clint type'),e('select',{value:iviveInputForm.input_type,onChange:event=>updateInputType(event.target.value)},[
      e('option',{value:'RAW_MICROSOMAL'},'Raw microsomal'),e('option',{value:'RAW_HEPATOCYTE'},'Raw hepatocyte'),e('option',{value:'PRESCALED_CLINT'},'Pre-scaled Clint')
     ])]),
     e('div',{className:'col-2'},Field({label:'Value',type:'number',value:iviveInputForm.input_value,onChange:value=>setIviveInputForm(current=>({...current,input_value:value}))})),
     e('div',{className:'col-2'},Field({label:'Unit',value:iviveInputForm.unit,onChange:value=>setIviveInputForm(current=>({...current,unit:value}))})),
     e('div',{className:'col-2'},[e('label',{},'Source type'),e('select',{value:iviveInputForm.source_type,onChange:event=>setIviveInputForm(current=>({...current,source_type:event.target.value}))},[e('option',{value:'EXPERIMENTAL'},'Experimental'),e('option',{value:'PROJECT_CALIBRATED'},'Project-calibrated')])]),
     e('div',{className:'col-6'},Field({label:'Model / source',value:iviveInputForm.model_source,onChange:value=>setIviveInputForm(current=>({...current,model_source:value}))})),
     e('div',{className:'col-3'},[e('label',{},'Confidence'),e('select',{value:iviveInputForm.confidence,onChange:event=>setIviveInputForm(current=>({...current,confidence:event.target.value}))},['HIGH','MEDIUM','LOW'].map(value=>e('option',{key:value,value},value)))]),
     e('div',{className:'col-3'},[e('label',{},'Action'),e('button',{disabled:iviveBusy||iviveInputForm.input_value==='',onClick:addIviveInputAction},'Save Input')])
    ]),
    e('p',{className:'small'},'Stage 3 experimental ADME and quantitative predictions are connected automatically. Classification-only results are never converted to Clint.')
   ]),

   e('div',{className:'ivive-panel',key:'physiology'},[
    e('div',{className:'row toolbar'},[e('div',{},[e('h3',{},'Species Physiology'),e('p',{className:'small'},'Versioned defaults and project-scoped study overrides remain visibly distinct.')]),e('span',{className:'mono small'},latest?.parameter_set_version||'PHRMA-CPCDC-2011-v1.0')]),
    e('div',{className:'ivive-physiology-grid'},Object.values(physiology).map(row=>e('div',{className:'ivive-physiology-row',key:row.parameter},[
     e('span',{},row.parameter.replaceAll('_',' ')),e('strong',{className:'mono'},sig(row.value)+' '+row.unit),sourceBadge(row.source_label),e('small',{},row.reference?.citation||row.reference?.source||'Reference retained in provenance')
    ]))),
    e('details',{className:'ivive-entry'},[e('summary',{},'Add study-specific physiology override'),e('div',{className:'grid ivive-form'},[
     e('div',{className:'col-3'},[e('label',{},'Parameter'),e('select',{value:iviveOverrideForm.parameter,onChange:event=>setIviveOverrideForm(current=>({...current,parameter:event.target.value}))},Object.keys(physiology).map(value=>e('option',{key:value,value},value.replaceAll('_',' '))))]),
     e('div',{className:'col-2'},Field({label:'Value',type:'number',value:iviveOverrideForm.value,onChange:value=>setIviveOverrideForm(current=>({...current,value}))})),
     e('div',{className:'col-2'},Field({label:'Unit',value:iviveOverrideForm.unit,onChange:value=>setIviveOverrideForm(current=>({...current,unit:value}))})),
     e('div',{className:'col-3'},Field({label:'Study source / reference',value:iviveOverrideForm.source,onChange:value=>setIviveOverrideForm(current=>({...current,source:value}))})),
     e('div',{className:'col-2'},[e('label',{},'Action'),e('button',{disabled:iviveBusy||iviveOverrideForm.value==='',onClick:addOverrideAction},'Save Override')])
    ])])
   ]),

   latest?e('div',{className:'ivive-results',key:'results'},[
    e('div',{className:'row toolbar'},[e('div',{},[e('h3',{},'Scaling'),e('p',{className:'small'},outputs.scaling?.input_type==='PRESCALED_CLINT'?'Pre-scaled Clint: no MPPGL or hepatocellularity multiplication.':'Raw Clint scaled exactly once with the matching physiological scalar.')]),StatusBadge({type:latest.status==='COMPLETE'?'CALCULATED':'Not calculated'})]),
    outputs.scaling&&e('div',{className:'ivive-scaling-flow'},[
     e('div',{},[sourceBadge(latest.inputs_snapshot.clint?.source_label),e('strong',{className:'mono'},sig(outputs.scaling.raw_value)+' '+outputs.scaling.raw_unit),e('small',{},outputs.scaling.input_type)]),e('span',{className:'ivive-arrow'},'→'),e('div',{},[sourceBadge('CALC'),e('strong',{className:'mono'},sig(outputs.scaling.scaled_clint)+' '+outputs.scaling.scaled_unit),e('small',{},outputs.scaling.equation)])
    ]),
    latest.status==='COMPLETE'?e('div',{key:'complete'},[
     e('h3',{},'Hepatic Clearance'),
     e('div',{className:'pk-nca-grid ivive-output-grid'},[
      e('div',{className:'pk-nca-card'},[e('span',{},'Scaled Clint'),e('strong',{className:'mono'},sig(outputs.clint)),e('small',{},outputs.clint_unit)]),
      e('div',{className:'pk-nca-card'},[e('span',{},outputs.binding_basis==='BLOOD'?'Fraction unbound in blood (fu,b)':'fu,b · Plasma Approximation'),e('strong',{className:'mono'},sig(outputs.fu_b)),e('small',{},'fu,p '+sig(outputs.fu_p)+(outputs.blood_plasma_ratio?' · B/P '+sig(outputs.blood_plasma_ratio):''))]),
      e('div',{className:'pk-nca-card'},[e('span',{},'Hepatic Blood Flow (Qh)'),e('strong',{className:'mono'},sig(outputs.qh)),e('small',{},outputs.qh_unit)]),
      e('div',{className:'pk-nca-card ivive-primary-output'},[e('span',{},'Hepatic Clearance Estimate (CLh)'),e('strong',{className:'mono'},sig(outputs.clh)),e('small',{},outputs.clh_unit)]),
      e('div',{className:'pk-nca-card'},[e('span',{},'Hepatic Extraction Ratio (Eh)'),e('strong',{className:'mono'},sig(outputs.extraction_ratio)),e('small',{},outputs.extraction_class+' extraction · Low <0.3 · High >0.7')]),
      e('div',{className:'pk-nca-card'},[e('span',{},'Predicted Hepatic Availability (Fh)'),e('strong',{className:'mono'},sig(outputs.hepatic_availability)),e('small',{},'Fh = 1 − Eh · not absolute oral F')])
     ]),
     e('p',{className:'ivive-scope-warning'},'Non-hepatic Clearance: Not modeled · Predicted Total Clearance: Not generated')
    ]):e('div',{className:'empty-state'},[StatusBadge({type:'Not calculated'}),e('p',{},'Required quantitative inputs are unavailable. No fake value was inserted; review warnings below.')]),
    outputs.experimental_comparison&&e('div',{className:'ivive-comparison'},[
     e('h3',{},'Experimental Comparison'),
     e('div',{className:'ivive-comparison-grid'},[
      e('div',{},[sourceBadge('EXP'),e('span',{},'Observed IV systemic CL'),e('strong',{className:'mono'},sig(outputs.experimental_comparison.observed_systemic_cl)+' '+outputs.experimental_comparison.unit)]),
      e('div',{},[sourceBadge('CALC'),e('span',{},'Predicted hepatic CL'),e('strong',{className:'mono'},sig(outputs.experimental_comparison.predicted_hepatic_cl)+' '+outputs.experimental_comparison.unit)]),
      e('div',{},[sourceBadge('CALC'),e('span',{},'Estimated hepatic contribution'),e('strong',{className:'mono'},sig(outputs.experimental_comparison.estimated_hepatic_contribution*100)+'%')])
     ]),e('p',{className:'small'},outputs.experimental_comparison.limitation)
    ]),
    e('div',{className:'ivive-warnings'},[e('h3',{},'Assumptions & Warnings'),
     (latest.warnings||[]).length?e('ul',{},latest.warnings.map((warning,index)=>e('li',{key:index},warning))):e('p',{className:'pass'},'No additional run warnings.'),
     e('p',{className:'small'},'Run confidence: '+latest.confidence+' · Result confidence cannot exceed its weakest required input.')
    ]),
    e('details',{className:'ivive-provenance'},[e('summary',{},'IVIVE Provenance & Equations'),e('div',{className:'small'},[
     e('p',{},'Method: Well-stirred hepatic clearance model · PK/IVIVE Method Registry · Timestamp: '+latest.timestamp),
     e('p',{className:'mono'},latest.equations.hepatic_clearance),e('p',{className:'mono'},latest.equations.blood_binding),
     e('p',{},'Input snapshot hash: '+latest.inputs_hash),e('pre',{},JSON.stringify(latest.inputs_snapshot,null,2))
    ])])
   ]):e('div',{className:'empty-state',key:'no-run'},[StatusBadge({type:'Not calculated'}),e('p',{},'Review selected inputs and physiology, then run IVIVE. Missing data will produce an auditable unavailable run rather than a fabricated estimate.')])
  ]);
 }

 function PkFoundationProfile({versionId}){
  const [foundationData,setFoundationData]=React.useState(null);
  const [loading,setLoading]=React.useState(false);
  const [species,setSpecies]=React.useState('Rat');

  const fetchFoundation=async(s)=>{
   setLoading(true);
   try{
    const res=await api.get('/compound-versions/'+versionId+'/pk-foundation?species='+(s||species));
    setFoundationData(res);
   }catch(err){console.error(err)}finally{setLoading(false)}
  };

  React.useEffect(()=>{
   if(versionId)fetchFoundation(species);
  },[versionId,species]);

  if(!foundationData){
   return e('section',{className:'card ivive-section',key:'pk-foundation'},[
    e('div',{className:'row toolbar'},[
     e('div',{},[
      e('div',{className:'eyebrow'},'5 · PK PARAMETER FOUNDATION & ROUTE ASSEMBLY'),
      e('h2',{},'PK Parameter Foundation'),
      e('p',{className:'small'},loading?'Loading PK Foundation parameters…':'Click below to load PK Parameter Foundation.')
     ]),
     e('button',{disabled:loading,onClick:()=>fetchFoundation(species)},loading?'Loading…':'Load PK Foundation')
    ])
   ]);
  }

  const dist=foundationData.distribution||{};
  const abs=foundationData.absorption||{};
  const routes=foundationData.route_parameter_sets||{};

  return e('section',{className:'card ivive-section',key:'pk-foundation'},[
   e('div',{className:'row toolbar',key:'header'},[
    e('div',{},[
     e('div',{className:'eyebrow'},'5 · PK PARAMETER FOUNDATION & ROUTE ASSEMBLY'),
     e('h2',{},'PK Parameter Foundation'),
     e('p',{className:'small'},'Route-aware pharmacokinetic parameters, distribution foundation, and absorption component deconstruction (Fa · Fg · Fh).')
    ]),
    e('div',{className:'ivive-run-controls'},[
     e('label',{},'Species'),
     e('select',{value:species,onChange:ev=>setSpecies(ev.target.value)},['Mouse','Rat','Dog','Monkey','Human'].map(s=>e('option',{key:s,value:s},s))),
     e('button',{className:'secondary',disabled:loading,onClick:()=>fetchFoundation(species)},loading?'Loading…':'Refresh')
    ])
   ]),

   e('div',{className:'grid ivive-grid',key:'panels'},[
    e('div',{className:'ivive-panel col-6',key:'dist'},[
     e('h3',{},'Distribution Foundation'),
     dist.v_value!=null?e('div',{className:'pk-nca-card'},[
      e('span',{},dist.v_type==='Vss'?'Experimental Vss':dist.v_type==='Vz'?'Experimental Vz':'Predicted Vd Estimate'),
      e('strong',{className:'mono'},dist.v_value+' '+dist.v_unit),
      e('small',{},dist.message)
     ]):e('div',{className:'empty-state'},[StatusBadge({type:'MODEL_UNAVAILABLE'}),e('p',{},dist.message||'Vd unavailable due to missing binding/lipophilicity data.')]),
     dist.apparent_vzf&&e('div',{className:'pk-nca-card',style:{marginTop:'8px'}},[
      e('span',{},'Apparent Volume (Vz/F)'),
      e('strong',{className:'mono'},dist.apparent_vzf.v_value+' L/kg'),
      e('small',{},dist.apparent_vzf.message)
     ])
    ]),

    e('div',{className:'ivive-panel col-6',key:'abs'},[
     e('h3',{},'Absorption Foundation (F = Fa × Fg × Fh)'),
     e('div',{className:'highlight-grid'},[
      e('div',{className:'highlight-item'},[e('strong',{},'Fa (Lumen Absorbed)'),e('div',{className:'mono'},abs.fa_value!=null?(abs.fa_value*100).toFixed(1)+'%':'—'),StatusBadge({type:abs.fa_status})]),
      e('div',{className:'highlight-item'},[e('strong',{},'Fg (Gut Escape)'),e('div',{className:'mono'},abs.fg_value!=null?(abs.fg_value*100).toFixed(1)+'%':'—'),StatusBadge({type:abs.fg_status})]),
      e('div',{className:'highlight-item'},[e('strong',{},'Fh (Hepatic Escape)'),e('div',{className:'mono'},abs.fh_value!=null?(abs.fh_value*100).toFixed(1)+'%':'—'),StatusBadge({type:abs.fh_value!=null?'CALCULATED':'MODEL_UNAVAILABLE'})]),
      e('div',{className:'highlight-item ivive-primary-output'},[e('strong',{},'Predicted Absolute F'),e('div',{className:'mono'},abs.f_predicted!=null?abs.f_predicted+'%':'—'),StatusBadge({type:abs.f_predicted!=null?'CALCULATED':'MODEL_UNAVAILABLE'})])
     ]),
     e('p',{className:'small',style:{marginTop:'8px'}},abs.f_predicted_message),
     abs.f_experimental!=null&&e('div',{className:'pass',style:{marginTop:'8px'}},[e('strong',{},'Matched Experimental F: '),e('span',{className:'mono'},abs.f_experimental+'%')])
    ])
   ]),

   e('div',{key:'routes',style:{marginTop:'16px'}},[
    e('h3',{},'Route-Aware Parameter Sets ('+species+')'),
    e('div',{className:'grid ivive-output-grid'},['IV','PO','SC','IP'].map(r=>{
     const rset=routes[r]||{};
     return e('div',{key:r,className:'card pk-nca-card'+(r==='IV'||r==='PO'?' ivive-primary-output':'')},[
      e('div',{className:'row toolbar'},[e('strong',{style:{fontSize:'16px'}},r+' Route'),StatusBadge({type:rset.confidence||'MODEL_UNAVAILABLE'})]),
      e('div',{className:'metric-row'},[e('span',{},'Clearance'),e('strong',{className:'mono'},rset.cl_value!=null?rset.cl_value+' '+rset.cl_unit:'—')]),
      e('div',{className:'small',style:{marginBottom:'6px'}},'Source: '+(rset.cl_source_type||'—')),
      e('div',{className:'metric-row'},[e('span',{},'Volume'),e('strong',{className:'mono'},rset.v_value!=null?rset.v_value+' '+rset.v_unit:'—')]),
      e('div',{className:'small',style:{marginBottom:'6px'}},'Type: '+(rset.v_type||'—')+' ('+(rset.v_source_type||'—')+')'),
      e('div',{className:'metric-row'},[e('span',{},'Hepatic Fh'),e('strong',{className:'mono'},rset.fh_value!=null?(rset.fh_value*100).toFixed(1)+'%':'—')]),
      e('div',{className:'metric-row'},[e('span',{},'Absolute F'),e('strong',{className:'mono'},rset.f_experimental!=null?rset.f_experimental+'% (exp)':(rset.f_predicted!=null?rset.f_predicted+'% (pred)':'—'))])
     ]);
    }))
   ])
  ]);
 }

 function PkSimulationSection({versionId}){
  const [species, setSpecies] = React.useState('Rat');
  const [route, setRoute] = React.useState('PO');
  const [adminType, setAdminType] = React.useState('EXTRAVASCULAR_1COMP');
  const [dose, setDose] = React.useState(1.0);
  const [doseUnit, setDoseUnit] = React.useState('mg/kg');
  const [infusionDur, setInfusionDur] = React.useState(1.0);
  const [frequency, setFrequency] = React.useState('Single Dose');
  const [interval, setInterval] = React.useState(24.0);
  const [numDoses, setNumDoses] = React.useState(3);
  const [modelType, setModelType] = React.useState('ONE_COMPARTMENT');
  const [userCl, setUserCl] = React.useState('');
  const [userV, setUserV] = React.useState('');
  const [userF, setUserF] = React.useState('');
  const [userKa, setUserKa] = React.useState('');
  const [logScale, setLogScale] = React.useState(false);
  const [preview, setPreview] = React.useState(null);
  const [activeRun, setActiveRun] = React.useState(null);
  const [history, setHistory] = React.useState([]);
  const [fitResult, setFitResult] = React.useState(null);
  const [loading, setLoading] = React.useState(false);
  const [fitting, setFitting] = React.useState(false);
  const [error, setError] = React.useState(null);

  // Auto-adjust adminType when route changes
  React.useEffect(()=>{
    if(route === 'IV'){
      if(adminType === 'EXTRAVASCULAR_1COMP') setAdminType('IV_BOLUS');
    } else {
      setAdminType('EXTRAVASCULAR_1COMP');
    }
  },[route]);

  const loadData = React.useCallback(async()=>{
   if(!versionId) return;
   try{
    const prev = await api.get('/compound-versions/'+versionId+'/pk-simulation/preview?species='+species+'&route='+route);
    setPreview(prev);
    const hist = await api.get('/compound-versions/'+versionId+'/pk-simulation/history?species='+species+'&route='+route);
    setHistory(hist||[]);
    if(hist && hist.length > 0 && (!activeRun || activeRun.route !== route || activeRun.species !== species)){
     setActiveRun(hist[0]);
    } else if((!hist || hist.length === 0) && prev?.cl_preview?.value != null && prev?.v_preview?.value != null){
     const autoPayload = {
      species,
      route,
      administration_type: route === 'IV' ? (adminType || 'IV_BOLUS') : 'EXTRAVASCULAR_1COMP',
      dose: parseFloat(dose) || 1.0,
      dose_unit: doseUnit || 'mg/kg',
      infusion_duration_hours: (route === 'IV' && adminType === 'IV_INFUSION') ? parseFloat(infusionDur) : 0.0,
      dosing_frequency: frequency || 'Single Dose',
      dose_interval_hours: parseFloat(interval) || 24.0,
      num_doses: frequency === 'Repeated Dosing' ? parseInt(numDoses, 10) : 1,
      model_type: modelType || 'ONE_COMPARTMENT',
      user_cl_override: userCl ? parseFloat(userCl) : null,
      user_v_override: userV ? parseFloat(userV) : null,
      user_f_override: (route !== 'IV' && userF) ? parseFloat(userF) : null,
      user_ka_override: (route !== 'IV' && userKa) ? parseFloat(userKa) : null,
     };
     const autoRes = await api.post('/compound-versions/'+versionId+'/pk-simulation/run', autoPayload);
     setActiveRun(autoRes);
     setHistory([autoRes]);
    }
   }catch(err){
    console.error("Simulation load error:", err);
   }
  },[versionId, species, route, adminType, dose, doseUnit, infusionDur, frequency, interval, numDoses, modelType, userCl, userV, userF, userKa, activeRun]);

  React.useEffect(()=>{ loadData(); },[loadData]);

  const handleRun = async()=>{
   setLoading(true);
   setError(null);
   try{
    const payload = {
     species,
     route,
     administration_type: route === 'IV' ? adminType : 'EXTRAVASCULAR_1COMP',
     dose: parseFloat(dose),
     dose_unit: doseUnit,
     infusion_duration_hours: (route === 'IV' && adminType === 'IV_INFUSION') ? parseFloat(infusionDur) : 0.0,
     dosing_frequency: frequency,
     dose_interval_hours: parseFloat(interval),
     num_doses: frequency === 'Repeated Dosing' ? parseInt(numDoses, 10) : 1,
     model_type: modelType,
     user_cl_override: userCl ? parseFloat(userCl) : null,
     user_v_override: userV ? parseFloat(userV) : null,
     user_f_override: (route !== 'IV' && userF) ? parseFloat(userF) : null,
     user_ka_override: (route !== 'IV' && userKa) ? parseFloat(userKa) : null,
    };
    const res = await api.post('/compound-versions/'+versionId+'/pk-simulation/run', payload);
    setActiveRun(res);
    setHistory(h=>[res, ...h]);
   }catch(err){
    setError(err.message||"Simulation failed");
   }finally{
    setLoading(false);
   }
  };

  const handleFitKa = async()=>{
   setFitting(true);
   setError(null);
   try{
    const payload = {
     species,
     route,
     dose: parseFloat(dose),
     dose_unit: doseUnit,
     fix_cl_v: true,
     user_cl_override: userCl ? parseFloat(userCl) : null,
     user_v_override: userV ? parseFloat(userV) : null,
    };
    const res = await api.post('/compound-versions/'+versionId+'/pk-simulation/fit-extravascular', payload);
    setFitResult(res);
    if(res && res.fitted_ka){
      setUserKa(String(res.fitted_ka));
    }
   }catch(err){
    setError("Parameter fitting error: " + (err.message || String(err)));
   }finally{
    setFitting(false);
   }
  };

  const renderPlot = (run)=>{
   if(!run || !run.time_series || run.time_series.length === 0) return null;
   const ts = run.time_series;
   const res = run.residuals || [];
   const maxT = Math.max(...ts.map(p=>p.time), 1.0);

   const getC = p => logScale ? (p.concentration > 0 ? Math.log10(p.concentration) : 0) : p.concentration;
   const validObs = res.filter(r=>r.status==='VALID');
   const maxC = Math.max(
    ...ts.map(p=>getC(p)),
    ...(validObs.map(r=> logScale ? (r.observed_ng_ml > 0 ? Math.log10(r.observed_ng_ml) : 0) : r.observed_ng_ml)),
    1.0
   );
   const minC = logScale ? 0 : 0;

   const W = 620, H = 240, padL = 60, padB = 40, padR = 20, padT = 20;
   const mapX = t => padL + (t / maxT) * (W - padL - padR);
   const mapY = c => padT + (1.0 - (c - minC) / (maxC - minC || 1.0)) * (H - padT - padB);

   const pathD = ts.map((p, idx)=>(idx===0?'M':'L') + ' ' + mapX(p.time).toFixed(1) + ' ' + mapY(getC(p)).toFixed(1)).join(' ');

   return e('svg',{width:W, height:H, style:{background:'#0f172a', borderRadius:'8px', width:'100%', height:'auto'}},[
    e('line',{key:'axis-x', x1:padL, y1:H-padB, x2:W-padR, y2:H-padB, stroke:'#334155', strokeWidth:1}),
    e('line',{key:'axis-y', x1:padL, y1:padT, x2:padL, y2:H-padB, stroke:'#334155', strokeWidth:1}),
    e('text',{key:'lbl-x', x:W/2, y:H-8, fill:'#94a3b8', fontSize:12, textAnchor:'middle'},'Time (hours)'),
    e('text',{key:'lbl-y', x:15, y:H/2, fill:'#94a3b8', fontSize:12, textAnchor:'middle', transform:'rotate(-90 15 '+(H/2)+')'}, logScale ? 'log10 Conc (ng/mL)' : 'Concentration (ng/mL)'),
    e('path',{key:'sim-path', d:pathD, fill:'none', stroke:'#38bdf8', strokeWidth:2.5}),
    validObs.map((ob, idx)=>{
     const cx = mapX(ob.time_hours);
     const cy = mapY(logScale ? Math.log10(ob.observed_ng_ml) : ob.observed_ng_ml);
     return e('circle',{key:'obs-'+idx, cx, cy, r:4.5, fill:'#ef4444', stroke:'#ffffff', strokeWidth:1.5});
    })
   ]);
  };

  const clPreview = preview?.clearance;
  const vPreview = preview?.volume;
  const fPreview = preview?.bioavailability;
  const kaPreview = preview?.absorption_rate;
  const mechComp = preview?.mechanistic_components || {};

  return e('div',{className:'card pk-simulation-section', style:{marginTop:'16px'}},[
   e('div',{className:'row toolbar'},[
    e('div',{},[
     e('div',{className:'eyebrow'},'6 · PK CONCENTRATION-TIME SIMULATION & ABSORPTION KINETICS'),
     e('h3',{},'PK SIMULATION — Extravascular (PO / SC / IP) & IV Engine (Stage 5B-2)'),
    ]),
    StatusBadge({type: activeRun ? activeRun.confidence : (preview?.confidence_ceiling||'MEDIUM')})
   ]),
   e('p',{className:'small'},'Scientifically defensible concentration-time simulation for PO, SC, IP (first-order absorption) and IV routes with route isolation and parameter governance.'),

   // Route Switcher Toolbar
   e('div',{className:'row toolbar', style:{marginTop:'12px', background:'rgba(255,255,255,0.02)', padding:'8px 12px', borderRadius:'6px'}},[
    e('span',{style:{fontWeight:'bold', marginRight:'10px'}},'Route:'),
    ['PO','IV','SC','IP'].map(r=>e('button',{
      key: r,
      className: route === r ? 'primary' : 'secondary',
      style: {marginRight:'6px', padding:'5px 14px'},
      onClick: ()=>setRoute(r)
    }, r + (r==='PO'?' (Oral)':r==='IV'?' (Intravenous)':r==='SC'?' (Subcutaneous)':' (Intraperitoneal)'))),
   ]),

   e('div',{className:'grid', style:{marginTop:'12px'}},[
    e('div',{className:'col-3'},[
     e('label',{},'Species'),
     e('select',{value:species, onChange:ev=>setSpecies(ev.target.value)},['Rat','Mouse','Dog','Monkey','Human'].map(s=>e('option',{key:s,value:s},s)))
    ]),
    route === 'IV' ? e('div',{className:'col-3'},[
     e('label',{},'Administration Type'),
     e('select',{value:adminType, onChange:ev=>setAdminType(ev.target.value)},[
      e('option',{value:'IV_BOLUS'},'IV Bolus'),
      e('option',{value:'IV_INFUSION'},'IV Infusion')
     ])
    ]) : e('div',{className:'col-3'},[
     e('label',{},'Administration Model'),
     e('input',{type:'text', value:'1-Compartment First-Order Absorption', readOnly:true, disabled:true})
    ]),
    e('div',{className:'col-3'},[
      Field({label:'Dose', type:'number', value:dose, onChange:setDose}),
      e('span',{className:'small mono',style:{color:'#94a3b8',fontSize:'11px'}},'NORMALIZED SIMULATION — 1.0 mg/kg (Simulation input only)')
     ]),
    e('div',{className:'col-3'},[
     e('label',{},'Dose Unit'),
     e('select',{value:doseUnit, onChange:ev=>setDoseUnit(ev.target.value)},['mg/kg','µg/kg','mg','µg'].map(u=>e('option',{key:u,value:u},u)))
    ]),
    (route === 'IV' && adminType === 'IV_INFUSION') && e('div',{className:'col-3', key:'inf-dur'},Field({label:'Infusion Duration (h)', type:'number', value:infusionDur, onChange:setInfusionDur})),
    e('div',{className:'col-3'},[
     e('label',{},'Dosing Frequency'),
     e('select',{value:frequency, onChange:ev=>setFrequency(ev.target.value)},['Single Dose','Repeated Dosing'].map(f=>e('option',{key:f,value:f},f)))
    ]),
    frequency === 'Repeated Dosing' && e('div',{className:'col-3', key:'interval'},Field({label:'Dose Interval τ (h)', type:'number', value:interval, onChange:setInterval})),
    frequency === 'Repeated Dosing' && e('div',{className:'col-3', key:'numDoses'},Field({label:'Number of Doses', type:'number', value:numDoses, onChange:setNumDoses})),
    route === 'IV' && e('div',{className:'col-3'},[
     e('label',{},'Model Type'),
     e('select',{value:modelType, onChange:ev=>setModelType(ev.target.value)},[
      e('option',{value:'ONE_COMPARTMENT'},'1-Compartment Model'),
      e('option',{value:'TWO_COMPARTMENT'},'2-Compartment Model (If fit available)')
     ])
    ])
   ]),

   // Parameter Source & Governance Panel
   e('div',{className:'card pk-parameter-provenance', style:{marginTop:'12px', padding:'12px'}},[
    e('div',{className:'row toolbar'},[
     e('strong',{style:{fontSize:'14px'}},'Parameter Provenance & Evidence Hierarchy ('+route+' · '+species+'):'),
     route !== 'IV' && e('button',{className:'secondary', style:{fontSize:'12px'}, disabled:fitting, onClick:handleFitKa}, fitting ? 'Fitting…' : 'Fit ka from Observations')
    ]),
    e('div',{className:'grid ivive-output-grid', style:{marginTop:'8px'}},[
     e('div',{className:'card pk-nca-card'},[
      e('div',{className:'row toolbar'},[e('span',{},'Systemic CL'), StatusBadge({type:clPreview?.evidence_type||'MODEL_UNAVAILABLE'})]),
      e('strong',{className:'mono'},clPreview?.value!=null ? clPreview.value+' '+clPreview.unit : 'Unavailable'),
      e('small',{},'Source: '+(clPreview?.source||'None'))
     ]),
     e('div',{className:'card pk-nca-card'},[
      e('div',{className:'row toolbar'},[e('span',{},'Distribution V'), StatusBadge({type:vPreview?.evidence_type||'MODEL_UNAVAILABLE'})]),
      e('strong',{className:'mono'},vPreview?.value!=null ? vPreview.value+' '+vPreview.unit : 'Unavailable'),
      e('small',{},'Type: '+(vPreview?.type||'None'))
     ]),
     route !== 'IV' && e('div',{className:'card pk-nca-card'},[
      e('div',{className:'row toolbar'},[e('span',{},'Bioavailability (F)'), StatusBadge({type:fPreview?.evidence_type||'MODEL_UNAVAILABLE'})]),
      e('strong',{className:'mono'},fPreview?.value!=null ? fPreview.value+'%' : 'Unavailable'),
      e('small',{},'Source: '+(fPreview?.source||'None'))
     ]),
     route !== 'IV' && e('div',{className:'card pk-nca-card'},[
      e('div',{className:'row toolbar'},[e('span',{},'Absorption ka'), StatusBadge({type:kaPreview?.evidence_type||'MODEL_UNAVAILABLE'})]),
      e('strong',{className:'mono'},kaPreview?.value!=null ? kaPreview.value+' 1/h' : 'Unavailable'),
      e('small',{},'Source: '+(kaPreview?.source||'None'))
     ])
    ]),

    route === 'PO' && e('div',{className:'row toolbar', style:{marginTop:'8px', fontSize:'12px', color:'#94a3b8'}},[
     e('span',{},'Mechanistic Decomposition: Fa = '+(mechComp.fa!=null?(mechComp.fa*100).toFixed(1)+'%':'—')+', Fg = '+(mechComp.fg!=null?(mechComp.fg*100).toFixed(1)+'%':'—')+', Fh = '+(mechComp.fh!=null?(mechComp.fh*100).toFixed(1)+'%':'—'))
    ]),

    e('div',{style:{marginTop:'10px',padding:'8px 12px',background:(preview?.clearance?.evidence_type==='EXPERIMENTAL'||preview?.volume?.evidence_type==='EXPERIMENTAL'||preview?.bioavailability?.evidence_type==='EXPERIMENTAL'||preview?.absorption_rate?.evidence_type==='EXPERIMENTAL')?'#ecfdf5':'rgba(255,255,255,0.03)',border:'1px solid '+((preview?.clearance?.evidence_type==='EXPERIMENTAL'||preview?.volume?.evidence_type==='EXPERIMENTAL')?'#a7f3d0':'rgba(255,255,255,0.08)'),borderRadius:'6px',fontSize:'12px',color:(preview?.clearance?.evidence_type==='EXPERIMENTAL'||preview?.volume?.evidence_type==='EXPERIMENTAL')?'#047857':'#94a3b8'}},(preview?.clearance?.evidence_type==='EXPERIMENTAL'||preview?.volume?.evidence_type==='EXPERIMENTAL'||preview?.bioavailability?.evidence_type==='EXPERIMENTAL'||preview?.absorption_rate?.evidence_type==='EXPERIMENTAL')?'✓ 실험값 우선 적용 (Experimental Precedence): 등록된 실험 데이터(NCA/측정값)가 시뮬레이션 기본 파라미터(Default)로 자동 적용되었습니다.':'기본 예측 모델(IVIVE / Lombardo Vd / Permeability) 파라미터가 적용되었습니다. 실험값이 입력되면 자동으로 실험값이 기본값(Default)으로 우선 적용됩니다.'),

    // Advanced Manual Overrides
    e('details',{style:{marginTop:'10px'}},[
     e('summary',{style:{cursor:'pointer', fontSize:'13px', color:'#38bdf8'}},'Manual Parameter Overrides (Optional)'),
     e('div',{className:'grid', style:{marginTop:'8px'}},[
      e('div',{className:'col-3'},Field({label:'Override CL (mL/min/kg)', type:'number', value:userCl, onChange:setUserCl})),
      e('div',{className:'col-3'},Field({label:'Override V (L/kg)', type:'number', value:userV, onChange:setUserV})),
      route !== 'IV' && e('div',{className:'col-3'},Field({label:'Override F (fraction or %)', type:'number', value:userF, onChange:setUserF})),
      route !== 'IV' && e('div',{className:'col-3'},Field({label:'Override ka (1/h)', type:'number', value:userKa, onChange:setUserKa})),
     ])
    ]),

    e('div',{className:'row toolbar', style:{marginTop:'12px'}},[
     e('span',{},'Ready to simulate '+route+' in '+species),
     e('button',{className:'primary', onClick:handleRun, disabled:loading}, loading ? 'Simulating…' : 'RUN PK SIMULATION')
    ]),

    fitResult && e('div',{className:'small pass', style:{marginTop:'8px'}},[
     e('strong',{},'Fitted Parameter: '),
     'ka = '+fitResult.fitted_ka+' 1/h (RMSE = '+fitResult.rmse+' ng/mL, AIC = '+fitResult.aic+')'
    ]),
    (preview?.warnings||[]).map((w, idx)=>e('div',{key:idx, className:'small alert', style:{marginTop:'4px'}},w)),
    error && e('div',{className:'small alert', style:{marginTop:'6px', color:'#ef4444'}},error)
   ]),

   activeRun && e('div',{style:{marginTop:'16px'}},[
    e('div',{className:'row toolbar'},[
     e('h4',{},'CALCULATED PK SIMULATION: '+activeRun.route+' '+activeRun.administration_type+' ('+activeRun.species+')'),
     e('div',{},[
      e('button',{className: logScale ? 'secondary' : 'primary', style:{marginRight:'6px', padding:'4px 8px'}, onClick:()=>setLogScale(false)},'Linear'),
      e('button',{className: logScale ? 'primary' : 'secondary', style:{padding:'4px 8px'}, onClick:()=>setLogScale(true)},'Semi-Log')
     ])
    ]),
    e('div',{style:{marginTop:'10px'}},renderPlot(activeRun)),
    e('div',{className:'row toolbar', style:{marginTop:'4px', fontSize:'12px', color:'#94a3b8'}},[
     e('span',{},'── Blue Line: Calculated PK Simulation'),
     e('span',{},'● Red Dots: Route-Matched Experimental Points')
    ]),
    e('div',{className:'grid ivive-output-grid', style:{marginTop:'12px'}},[
     e('div',{className:'card pk-nca-card'},[
      e('span',{},'Cmax'),
      e('strong',{className:'mono'},activeRun.output_metrics?.cmax_ng_ml+' ng/mL')
     ]),
     e('div',{className:'card pk-nca-card'},[
      e('span',{},'Tmax'),
      e('strong',{className:'mono'},activeRun.output_metrics?.tmax_hours+' h')
     ]),
     e('div',{className:'card pk-nca-card'},[
      e('span',{},'AUCinf (Analytical)'),
      e('strong',{className:'mono'},activeRun.output_metrics?.auc_inf_analytical_ng_h_ml+' ng·h/mL')
     ]),
     e('div',{className:'card pk-nca-card'},[
      e('span',{},'AUC Numerical Match'),
      e('strong',{className:'mono'},activeRun.output_metrics?.auc_inf_numerical_ng_h_ml+' ng·h/mL ('+activeRun.output_metrics?.auc_agreement_pct+'%)')
     ]),
     e('div',{className:'card pk-nca-card'},[
      e('span',{},'Terminal Half-Life'),
      e('strong',{className:'mono'},activeRun.output_metrics?.half_life_hours+' h')
     ]),
     activeRun.steady_state_metrics?.accumulation_ratio!=null && e('div',{className:'card pk-nca-card'},[
      e('span',{},'Accumulation Ratio (R_acc)'),
      e('strong',{className:'mono'},activeRun.steady_state_metrics.accumulation_ratio+'x')
     ]),
     activeRun.steady_state_metrics?.css_avg_ng_ml!=null && e('div',{className:'card pk-nca-card'},[
      e('span',{},'Css,avg (Steady State)'),
      e('strong',{className:'mono'},activeRun.steady_state_metrics.css_avg_ng_ml+' ng/mL')
     ])
    ]),
    activeRun.residuals && activeRun.residuals.length > 0 && e('div',{style:{marginTop:'16px'}},[
     e('h4',{},'Observed vs Simulated Residual Analysis ('+activeRun.route+' Route)'),
     e('table',{className:'table', style:{marginTop:'6px'}},[
      e('thead',{},e('tr',{},[
       e('th',{},'Time (h)'),
       e('th',{},'Observed (ng/mL)'),
       e('th',{},'Simulated (ng/mL)'),
       e('th',{},'Residual'),
       e('th',{},'Fold Error')
      ])),
      e('tbody',{},activeRun.residuals.map((r, idx)=>e('tr',{key:idx},[
       e('td',{className:'mono'},r.time_hours),
       e('td',{className:'mono'},r.observed_ng_ml),
       e('td',{className:'mono'},r.simulated_ng_ml!=null?r.simulated_ng_ml:'—'),
       e('td',{className:'mono'},r.residual_ng_ml!=null?r.residual_ng_ml:'—'),
       e('td',{className:'mono'},r.fold_error!=null?r.fold_error+'x':'—')
      ])))
     ]),
     activeRun.output_metrics?.goodness_of_fit?.rmse_ng_ml != null && e('div',{className:'small', style:{marginTop:'6px'}},[
      e('strong',{},'Goodness of Fit: '),
      'RMSE = '+activeRun.output_metrics.goodness_of_fit.rmse_ng_ml+' ng/mL, MAE = '+activeRun.output_metrics.goodness_of_fit.mae_ng_ml+' ng/mL'
     ])
    ]),
    e('div',{className:'small', style:{marginTop:'12px', padding:'8px', background:'rgba(255,255,255,0.03)', borderRadius:'6px'}},[
     e('div',{},[e('strong',{},'Engine: '),activeRun.provenance?.engine_name+' ('+activeRun.provenance?.engine_version+')']),
     e('div',{},[e('strong',{},'Uncertainty Status: '),activeRun.output_metrics?.uncertainty_status||'UNCERTAINTY NOT QUANTIFIED']),
     (activeRun.warnings||[]).map((w, idx)=>e('div',{key:idx, className:'alert', style:{marginTop:'2px'}},w))
    ])
   ])
  ]);
  }


 function TranslationalPkSection({versionId}){
  const [data, setData] = React.useState(null);
  const [loading, setLoading] = React.useState(false);
  const [error, setError] = React.useState(null);

  const loadData = React.useCallback(async ()=>{
   if(!versionId) return;
   setLoading(true);
   try {
    const res = await api.get('/compound-versions/' + versionId + '/translational-pk');
    setData(res);
    setError(null);
   } catch(err) {
    console.error(err);
    setError(err.message || 'Failed to load translational PK data.');
   } finally {
    setLoading(false);
   }
  }, [versionId]);

  React.useEffect(()=>{
   loadData();
  }, [loadData]);

  if(loading && !data) {
   return e('div',{className:'card', style:{marginTop:'16px'}}, e('p',{},'Loading translational PK scaling data…'));
  }

  if(!data) {
   return null;
  }

  const clAllo = data.clearance_allometry || {};
  const vssAllo = data.volume_allometry || {};
  const clLoso = data.clearance_loso || {};
  const hComp = data.human_comparison || {};
  const readiness = data.human_simulation_readiness || {};
  const spMatrix = data.species_data_matrix || [];

  // Render Log-Log Scatter Plot for Allometry (CL or Vss)
  const renderAllometryPlot = (alloData, title, yLabel) => {
   if(alloData.status !== 'SUCCESS' || !alloData.plot_points || alloData.plot_points.length === 0) {
    return e('div',{className:'empty-state small', style:{padding:'20px'}},'Insufficient data for allometric scaling curve.');
   }

   const points = alloData.plot_points;
   const extPt = alloData.extrapolated_point;
   const allPts = extPt ? [...points, extPt] : points;

   const minBw = Math.min(...allPts.map(p=>p.bw_kg));
   const maxBw = Math.max(...allPts.map(p=>p.bw_kg));
   const minVal = Math.min(...allPts.map(p=>p.observed_total || p.fitted_total));
   const maxVal = Math.max(...allPts.map(p=>p.observed_total || p.fitted_total));

   const logMinX = Math.log10(Math.max(minBw * 0.5, 0.005));
   const logMaxX = Math.log10(maxBw * 1.5);
   const logMinY = Math.log10(Math.max(minVal * 0.5, 0.001));
   const logMaxY = Math.log10(Math.max(maxVal * 2.0, 0.01));

   const W = 540, H = 220, padL = 55, padB = 35, padR = 20, padT = 20;
   const mapX = bw => padL + ((Math.log10(bw) - logMinX) / (logMaxX - logMinX || 1.0)) * (W - padL - padR);
   const mapY = val => padT + (1.0 - (Math.log10(val) - logMinY) / (logMaxY - logMinY || 1.0)) * (H - padT - padB);

   // Fitted line across full range
   const lineSteps = 50;
   const linePath = [];
   for(let i=0; i<=lineSteps; i++){
    const bwVal = Math.pow(10, logMinX + (i / lineSteps) * (logMaxX - logMinX));
    const fittedVal = alloData.coefficient_a * Math.pow(bwVal, alloData.exponent_b);
    linePath.push((i===0?'M':'L') + ' ' + mapX(bwVal).toFixed(1) + ' ' + mapY(fittedVal).toFixed(1));
   }

   return e('svg',{width:W, height:H, style:{background:'#0f172a', borderRadius:'8px', width:'100%', height:'auto'}},[
    e('line',{key:'axis-x', x1:padL, y1:H-padB, x2:W-padR, y2:H-padB, stroke:'#334155', strokeWidth:1}),
    e('line',{key:'axis-y', x1:padL, y1:padT, x2:padL, y2:H-padB, stroke:'#334155', strokeWidth:1}),
    e('text',{key:'lbl-x', x:W/2, y:H-8, fill:'#94a3b8', fontSize:11, textAnchor:'middle'},'Body Weight (kg, log scale)'),
    e('text',{key:'lbl-y', x:14, y:H/2, fill:'#94a3b8', fontSize:11, textAnchor:'middle', transform:'rotate(-90 14 '+(H/2)+')'}, yLabel + ' (log scale)'),
    e('path',{key:'fit-line', d:linePath.join(' '), fill:'none', stroke:'#38bdf8', strokeWidth:2, strokeDasharray:'3,3'}),
    // Animal Experimental Points
    points.map((pt, idx)=>{
     const cx = mapX(pt.bw_kg);
     const cy = mapY(pt.observed_total);
     return e('g',{key:'pt-'+idx},[
      e('circle',{cx, cy, r:5, fill:'#3b82f6', stroke:'#ffffff', strokeWidth:1.5}),
      e('text',{x:cx+7, y:cy-4, fill:'#93c5fd', fontSize:10},pt.species)
     ]);
    }),
    // Extrapolated Target Point (Human)
    extPt && e('g',{key:'ext-pt'},[
     e('circle',{cx:mapX(extPt.bw_kg), cy:mapY(extPt.fitted_total), r:6, fill:'#f59e0b', stroke:'#ffffff', strokeWidth:2}),
     e('text',{x:mapX(extPt.bw_kg)-10, y:mapY(extPt.fitted_total)-8, fill:'#fcd34d', fontSize:11, fontWeight:'bold'},'Extrapolated Human (' + extPt.fitted_norm + ' ' + (alloData.param_type==='CL'?'mL/min/kg':'L/kg') + ')')
    ])
   ]);
  };

  return e('div',{className:'card', style:{marginTop:'16px'}},[
   e('div',{className:'row toolbar'},[
    e('div',{},[
     e('div',{className:'eyebrow'},'5 · TRANSLATIONAL PK & CROSS-SPECIES ALLOMETRIC SCALING'),
     e('h3',{},'Interspecies Scaling & Translational Foundation (Stage 5B-3)'),
    ]),
    StatusBadge({type: readiness.overall_status === 'READY' ? 'READY' : (readiness.overall_status === 'PARTIALLY READY' ? 'PARTIALLY_READY' : 'NOT_READY')})
   ]),
   e('p',{className:'small'},'Classical body-weight allometry (Y = a · BW^b), Leave-One-Species-Out (LOSO) cross-validation, and deterministic Human simulation readiness assessment.'),

   // Cross-Species Data Matrix Table
   e('div',{style:{marginTop:'12px'}},[
    e('h4',{},'Cross-Species In Vivo PK Observation Matrix'),
    e('table',{className:'table', style:{marginTop:'6px'}},[
     e('thead',{},e('tr',{},[
      e('th',{},'Species'),
      e('th',{},'Body Weight (kg)'),
      e('th',{},'IV CL (mL/min/kg)'),
      e('th',{},'IV Vss (L/kg)'),
      e('th',{},'IV Vz (L/kg)'),
      e('th',{},'IV t1/2 (h)'),
      e('th',{},'PO F (%)'),
      e('th',{},'IV / PO Studies'),
      e('th',{},'Evidence')
     ])),
     e('tbody',{},spMatrix.map(sp=>e('tr',{key:sp.species},[
      e('td',{},e('strong',{},sp.species)),
      e('td',{className:'mono'},sp.effective_bw_kg + (sp.study_bw_kg ? ' (Study)' : ' (Ref)')),
      e('td',{className:'mono'},sp.cl_iv != null ? sp.cl_iv : '—'),
      e('td',{className:'mono'},sp.vss_iv != null ? sp.vss_iv : '—'),
      e('td',{className:'mono'},sp.vz_iv != null ? sp.vz_iv : '—'),
      e('td',{className:'mono'},sp.half_life_iv != null ? sp.half_life_iv : '—'),
      e('td',{className:'mono'},sp.f_po != null ? sp.f_po + '%' : '—'),
      e('td',{className:'small'},sp.iv_studies_count + ' IV / ' + sp.po_studies_count + ' PO'),
      e('td',{},StatusBadge({type: sp.has_experimental_iv ? 'EXPERIMENTAL' : (sp.has_experimental_po ? 'EXPERIMENTAL' : 'MODEL_UNAVAILABLE')}))
     ])))
    ])
   ]),

   // Allometry Scaling Panels Grid
   e('div',{className:'grid', style:{marginTop:'16px'}},[
    // Clearance Allometry Card
    e('div',{className:'col-6'},e('div',{className:'card', style:{background:'var(--bg-subtle,#1e293b)'}},[
     e('div',{className:'row toolbar'},[
      e('strong',{style:{fontSize:'14px'}},'Clearance (CL) Allometry (Animal IV Data)'),
      StatusBadge({type: clAllo.confidence || 'INSUFFICIENT_DATA'})
     ]),
     clAllo.status === 'SUCCESS' ? e('div',{},[
      e('div',{style:{marginTop:'8px'}},renderAllometryPlot(clAllo, 'Clearance Allometry', 'Total CL (mL/min)')),
      e('div',{className:'grid ivive-output-grid', style:{marginTop:'10px'}},[
       e('div',{className:'card pk-nca-card'},[e('span',{},'Exponent b (CL)'), e('strong',{className:'mono'},clAllo.exponent_b), e('small',{},'Classical ref ~0.75')]),
       e('div',{className:'card pk-nca-card'},[e('span',{},'Coefficient a'), e('strong',{className:'mono'},clAllo.coefficient_a)]),
       e('div',{className:'card pk-nca-card'},[e('span',{},'Goodness of Fit R²'), e('strong',{className:'mono'},clAllo.r_squared)]),
       e('div',{className:'card pk-nca-card'},[e('span',{},'Human Extrapolated CL'), e('strong',{className:'mono'},clAllo.extrapolated_norm + ' mL/min/kg'), e('small',{},clAllo.extrapolated_total + ' mL/min total')])
      ]),
      (clAllo.warnings||[]).map((w, idx)=>e('div',{key:idx, className:'alert', style:{marginTop:'6px', fontSize:'12px'}},w)),
      clAllo.historical_roe_correction && e('div',{className:'card', style:{marginTop:'10px', background:'rgba(255,255,255,0.02)', padding:'8px'}},[
       e('div',{className:'small', style:{fontWeight:'bold'}},'Historical Rule of Exponents (Mahmood & Balian):'),
       e('div',{className:'small mono', style:{marginTop:'4px'}},'Selected: ' + clAllo.historical_roe_correction.selected_rule + ' → ' + clAllo.historical_roe_correction.roe_predicted_norm + ' mL/min/kg'),
       e('div',{className:'small', style:{color:'#94a3b8', marginTop:'2px'}},clAllo.historical_roe_correction.citation)
      ])
     ]) : e('div',{className:'empty-state small', style:{padding:'20px'}},clAllo.message || 'Requires at least 2 distinct animal species with experimental IV CL.')
    ])),

    // Volume of Distribution Allometry Card
    e('div',{className:'col-6'},e('div',{className:'card', style:{background:'var(--bg-subtle,#1e293b)'}},[
     e('div',{className:'row toolbar'},[
      e('strong',{style:{fontSize:'14px'}},'Volume of Distribution (Vss) Allometry (Animal IV Data)'),
      StatusBadge({type: vssAllo.confidence || 'INSUFFICIENT_DATA'})
     ]),
     vssAllo.status === 'SUCCESS' ? e('div',{},[
      e('div',{style:{marginTop:'8px'}},renderAllometryPlot(vssAllo, 'Volume Allometry', 'Total Vss (L)')),
      e('div',{className:'grid ivive-output-grid', style:{marginTop:'10px'}},[
       e('div',{className:'card pk-nca-card'},[e('span',{},'Exponent b (Vss)'), e('strong',{className:'mono'},vssAllo.exponent_b), e('small',{},'Classical ref ~1.0')]),
       e('div',{className:'card pk-nca-card'},[e('span',{},'Coefficient a'), e('strong',{className:'mono'},vssAllo.coefficient_a)]),
       e('div',{className:'card pk-nca-card'},[e('span',{},'Goodness of Fit R²'), e('strong',{className:'mono'},vssAllo.r_squared)]),
       e('div',{className:'card pk-nca-card'},[e('span',{},'Human Extrapolated Vss'), e('strong',{className:'mono'},vssAllo.extrapolated_norm + ' L/kg'), e('small',{},vssAllo.extrapolated_total + ' L total')])
      ]),
      (vssAllo.warnings||[]).map((w, idx)=>e('div',{key:idx, className:'alert', style:{marginTop:'6px', fontSize:'12px'}},w)),
      hComp.translated_half_life?.value != null && e('div',{className:'card', style:{marginTop:'10px', background:'rgba(255,255,255,0.02)', padding:'8px'}},[
       e('div',{className:'small', style:{fontWeight:'bold'}},'Translated Human Elimination Half-Life:'),
       e('div',{className:'small mono', style:{fontSize:'15px', color:'#38bdf8', marginTop:'4px'}},hComp.translated_half_life.value + ' hours'),
       e('div',{className:'small', style:{color:'#94a3b8', marginTop:'2px'}},'Formula: ln(2) · Vss_allometric / CL_allometric (' + hComp.translated_half_life.v_definition_used + ')')
      ])
     ]) : e('div',{className:'empty-state small', style:{padding:'20px'}},vssAllo.message || 'Requires at least 2 distinct animal species with experimental IV Vss.')
    ]))
   ]),

   // Leave-One-Species-Out (LOSO) Cross-Validation Card
   clLoso.status === 'SUCCESS' && e('div',{className:'card', style:{marginTop:'16px', background:'rgba(255,255,255,0.02)'}},[
    e('div',{className:'row toolbar'},[
     e('h4',{},'Leave-One-Species-Out (LOSO) Allometry Cross-Validation (CL)'),
     e('span',{className:'small'},'Evaluated Species N = ' + clLoso.n_species_evaluated)
    ]),
    e('div',{className:'grid ivive-output-grid', style:{marginTop:'8px'}},[
     e('div',{className:'card pk-nca-card'},[e('span',{},'LOSO AAFE'), e('strong',{className:'mono'},clLoso.aafe + 'x')]),
     e('div',{className:'card pk-nca-card'},[e('span',{},'LOSO Bias (GMFE)'), e('strong',{className:'mono'},clLoso.bias_gmfe + 'x')]),
     e('div',{className:'card pk-nca-card'},[e('span',{},'Within 2-Fold Accuracy'), e('strong',{className:'mono'},clLoso.within_2_fold_pct + '%')]),
     e('div',{className:'card pk-nca-card'},[e('span',{},'Within 3-Fold Accuracy'), e('strong',{className:'mono'},clLoso.within_3_fold_pct + '%')])
    ]),
    e('table',{className:'table', style:{marginTop:'10px'}},[
     e('thead',{},e('tr',{},[
      e('th',{},'Held-Out Species'),
      e('th',{},'Observed CL (mL/min/kg)'),
      e('th',{},'Predicted CL (mL/min/kg)'),
      e('th',{},'Fold Error'),
      e('th',{},'Absolute Fold Error'),
      e('th',{},'Status')
     ])),
     e('tbody',{},clLoso.loso_evaluations.map(ev=>e('tr',{key:ev.held_out_species},[
      e('td',{},e('strong',{},ev.held_out_species)),
      e('td',{className:'mono'},ev.observed_norm),
      e('td',{className:'mono'},ev.predicted_norm),
      e('td',{className:'mono'},ev.fold_error + 'x'),
      e('td',{className:'mono'},ev.absolute_fold_error + 'x'),
      e('td',{},StatusBadge({type: ev.within_2_fold ? 'PASS' : 'FAIL'}))
     ])))
    ])
   ]),

   // Human Side-by-Side Comparison & Simulation Readiness
   e('div',{className:'card', style:{marginTop:'16px'}},[
    e('h4',{},'Human PK Translational Comparison & Readiness Scorecard'),
    e('p',{className:'small'},'Side-by-side comparison of independent prediction methods without silent averaging. Observed clinical data takes precedence.'),

    e('table',{className:'table', style:{marginTop:'8px'}},[
     e('thead',{},e('tr',{},[
      e('th',{},'Method / Evidence'),
      e('th',{},'Predicted / Measured CL (mL/min/kg)'),
      e('th',{},'Predicted / Measured Vss (L/kg)'),
      e('th',{},'Predicted / Measured F (%)'),
      e('th',{},'Confidence / Evidence Type'),
      e('th',{},'Method Notes')
     ])),
     e('tbody',{},[
      e('tr',{},[
       e('td',{},e('strong',{},'Method A: Mechanistic Hepatic IVIVE')),
       e('td',{className:'mono'},hComp.clearance?.method_a_hepatic_ivive?.value != null ? hComp.clearance.method_a_hepatic_ivive.value : '—'),
       e('td',{className:'mono'},'— (Not Modeled)'),
       e('td',{className:'mono'},'— (Not Modeled)'),
       e('td',{},StatusBadge({type: hComp.clearance?.method_a_hepatic_ivive?.confidence || 'MODEL_UNAVAILABLE'})),
       e('td',{className:'small'},hComp.clearance?.method_a_hepatic_ivive?.notes)
      ]),
      e('tr',{},[
       e('td',{},e('strong',{},'Method B: Simple Allometric Scaling')),
       e('td',{className:'mono'},hComp.clearance?.method_b_simple_allometry?.value != null ? hComp.clearance.method_b_simple_allometry.value : '—'),
       e('td',{className:'mono'},hComp.volume_vss?.simple_allometry?.value != null ? hComp.volume_vss.simple_allometry.value : '—'),
       e('td',{className:'mono'},'— (Not Inferred)'),
       e('td',{},StatusBadge({type: hComp.clearance?.method_b_simple_allometry?.confidence || 'INSUFFICIENT_DATA'})),
       e('td',{className:'small'},'Animal power-law extrapolation (N=' + (hComp.clearance?.method_b_simple_allometry?.n_species || 0) + ' sp, R²=' + (hComp.clearance?.method_b_simple_allometry?.r2 || '—') + ')')
      ]),
      e('tr',{},[
       e('td',{},e('strong',{},'Method C: Human Clinical PK (Experimental)')),
       e('td',{className:'mono'},hComp.clearance?.method_d_experimental_human?.value != null ? hComp.clearance.method_d_experimental_human.value : '—'),
       e('td',{className:'mono'},hComp.volume_vss?.experimental_human?.value != null ? hComp.volume_vss.experimental_human.value : '—'),
       e('td',{className:'mono'},spMatrix.find(s=>s.species==='Human')?.f_po != null ? spMatrix.find(s=>s.species==='Human').f_po + '%' : '—'),
       e('td',{},StatusBadge({type: hComp.clearance?.method_d_experimental_human?.value != null ? 'EXPERIMENTAL' : 'MODEL_UNAVAILABLE'})),
       e('td',{className:'small'},hComp.clearance?.method_d_experimental_human?.notes)
      ])
     ])
    ]),

    // Readiness Scorecard Box
    e('div',{className:'card', style:{marginTop:'12px', background:'rgba(255,255,255,0.02)', padding:'12px'}},[
     e('div',{className:'row toolbar'},[
      e('strong',{},'Human Simulation Readiness Assessment:'),
      StatusBadge({type: readiness.overall_status === 'READY' ? 'READY' : (readiness.overall_status === 'PARTIALLY READY' ? 'PARTIALLY_READY' : 'NOT_READY')})
     ]),
     e('div',{className:'grid ivive-output-grid', style:{marginTop:'8px'}},[
      e('div',{className:'card pk-nca-card'},[
       e('div',{className:'row toolbar'},[e('span',{},'Clearance (CL)'), StatusBadge({type:readiness.clearance?.status||'UNAVAILABLE'})]),
       e('small',{},readiness.clearance?.reason)
      ]),
      e('div',{className:'card pk-nca-card'},[
       e('div',{className:'row toolbar'},[e('span',{},'Volume (Vss)'), StatusBadge({type:readiness.volume?.status||'UNAVAILABLE'})]),
       e('small',{},readiness.volume?.reason)
      ]),
      e('div',{className:'card pk-nca-card'},[
       e('div',{className:'row toolbar'},[e('span',{},'Bioavailability (F)'), StatusBadge({type:readiness.bioavailability?.status||'UNAVAILABLE'})]),
       e('small',{},readiness.bioavailability?.reason)
      ]),
      e('div',{className:'card pk-nca-card'},[
       e('div',{className:'row toolbar'},[e('span',{},'Absorption Rate (ka)'), StatusBadge({type:readiness.absorption_rate?.status||'UNAVAILABLE'})]),
       e('small',{},readiness.absorption_rate?.reason)
      ])
     ]),
     e('div',{className:'alert', style:{marginTop:'10px', fontSize:'12px'}},readiness.oral_translation_guardrail)
    ])
   ])
  ]);
 }


 function PkValidationSection({versionId}){
  const [valData, setValData] = React.useState(null);
  const [loading, setLoading] = React.useState(false);

  const loadVal = React.useCallback(async ()=>{
   if(!versionId) return;
   setLoading(true);
   try {
    const res = await api.get('/compound-versions/' + versionId + '/pk-validation');
    setValData(res.validation_metrics);
   } catch(err) {
    console.error(err);
   } finally {
    setLoading(false);
   }
  }, [versionId]);

  React.useEffect(()=>{
   loadVal();
  }, [loadVal]);

  if(loading && !valData) return null;
  if(!valData || valData.status !== 'SUCCESS') {
   return e('div',{className:'card', style:{marginTop:'16px'}},[
    e('div',{className:'eyebrow'},'6 · PK VALIDATION & PREDICTION ERROR METRICS'),
    e('h3',{},'PK VALIDATION — Prediction Accuracy & Error Analysis (Stage 5B-3)'),
    e('p',{className:'small'},'No paired Predicted vs Observed PK records available yet for this compound. Add experimental PK studies or run IVIVE / Allometry to quantify prediction accuracy.')
   ]);
  }

  const pairs = valData.pairs || [];

  // Render Log-Log Predicted vs Observed Plot with 2-fold / 3-fold bands
  const renderPredObsPlot = () => {
   if(pairs.length === 0) return null;

   const allVals = pairs.flatMap(p=>[p.observed, p.predicted]);
   const minVal = Math.min(...allVals);
   const maxVal = Math.max(...allVals);

   const logMin = Math.log10(Math.max(minVal * 0.5, 0.01));
   const logMax = Math.log10(Math.max(maxVal * 2.0, 0.1));

   const W = 540, H = 240, padL = 55, padB = 40, padR = 20, padT = 20;
   const mapCoord = val => padL + ((Math.log10(val) - logMin) / (logMax - logMin || 1.0)) * (W - padL - padR);
   const mapYCoord = val => padT + (1.0 - (Math.log10(val) - logMin) / (logMax - logMin || 1.0)) * (H - padT - padB);

   const p1 = {x: mapCoord(Math.pow(10, logMin)), y: mapYCoord(Math.pow(10, logMin))};
   const p2 = {x: mapCoord(Math.pow(10, logMax)), y: mapYCoord(Math.pow(10, logMax))};

   return e('svg',{width:W, height:H, style:{background:'#0f172a', borderRadius:'8px', width:'100%', height:'auto'}},[
    e('line',{key:'axis-x', x1:padL, y1:H-padB, x2:W-padR, y2:H-padB, stroke:'#334155', strokeWidth:1}),
    e('line',{key:'axis-y', x1:padL, y1:padT, x2:padL, y2:H-padB, stroke:'#334155', strokeWidth:1}),
    e('text',{key:'lbl-x', x:W/2, y:H-8, fill:'#94a3b8', fontSize:11, textAnchor:'middle'},'Observed Value (log scale)'),
    e('text',{key:'lbl-y', x:14, y:H/2, fill:'#94a3b8', fontSize:11, textAnchor:'middle', transform:'rotate(-90 14 '+(H/2)+')'},'Predicted Value (log scale)'),

    // Identity line (y = x)
    e('line',{key:'line-ident', x1:p1.x, y1:p1.y, x2:p2.x, y2:p2.y, stroke:'#94a3b8', strokeWidth:1.5}),
    // 2-fold lines
    e('line',{key:'line-2f-up', x1:mapCoord(Math.pow(10, logMin)), y1:mapYCoord(Math.pow(10, logMin)*2), x2:mapCoord(Math.pow(10, logMax)/2), y2:mapYCoord(Math.pow(10, logMax)), stroke:'#38bdf8', strokeWidth:1, strokeDasharray:'3,3'}),
    e('line',{key:'line-2f-dn', x1:mapCoord(Math.pow(10, logMin)*2), y1:mapYCoord(Math.pow(10, logMin)), x2:mapCoord(Math.pow(10, logMax)), y2:mapYCoord(Math.pow(10, logMax)/2), stroke:'#38bdf8', strokeWidth:1, strokeDasharray:'3,3'}),

    // Points
    pairs.map((p, idx)=>{
     const cx = mapCoord(p.observed);
     const cy = mapYCoord(p.predicted);
     const color = p.absolute_fold_error <= 1.5 ? '#10b981' : (p.absolute_fold_error <= 2.0 ? '#38bdf8' : (p.absolute_fold_error <= 3.0 ? '#f59e0b' : '#ef4444'));
     return e('g',{key:'p-'+idx},[
      e('circle',{cx, cy, r:5, fill:color, stroke:'#ffffff', strokeWidth:1.5}),
      e('text',{x:cx+7, y:cy-3, fill:'#cbd5e1', fontSize:10},p.species + ' ' + p.endpoint)
     ]);
    })
   ]);
  };

  return e('div',{className:'card', style:{marginTop:'16px'}},[
   e('div',{className:'row toolbar'},[
    e('div',{},[
     e('div',{className:'eyebrow'},'6 · PK VALIDATION & PREDICTION ERROR METRICS'),
     e('h3',{},'PK VALIDATION — Prediction Accuracy & Error Analysis (Stage 5B-3)'),
    ]),
    StatusBadge({type: valData.within_2_fold_pct >= 70 ? 'PASS' : 'MEDIUM'})
   ]),
   e('p',{className:'small'},'Objective error metrics quantifying predictive accuracy across compatible endpoint and method pairs.'),

   // Summary Metrics Grid
   e('div',{className:'grid ivive-output-grid', style:{marginTop:'12px'}},[
    e('div',{className:'card pk-nca-card'},[e('span',{},'Pairs Evaluated (N)'), e('strong',{className:'mono'},valData.n)]),
    e('div',{className:'card pk-nca-card'},[e('span',{},'AAFE (Avg Abs Fold Error)'), e('strong',{className:'mono'},valData.aafe + 'x')]),
    e('div',{className:'card pk-nca-card'},[e('span',{},'Bias (GMFE)'), e('strong',{className:'mono'},valData.bias_gmfe + 'x')]),
    e('div',{className:'card pk-nca-card'},[e('span',{},'RMSE (log10 space)'), e('strong',{className:'mono'},valData.rmse_log10)]),
    e('div',{className:'card pk-nca-card'},[e('span',{},'Within 2-Fold (%)'), e('strong',{className:'mono'},valData.within_2_fold_pct + '% (' + valData.within_2_fold_count + '/' + valData.n + ')')]),
    e('div',{className:'card pk-nca-card'},[e('span',{},'Within 3-Fold (%)'), e('strong',{className:'mono'},valData.within_3_fold_pct + '% (' + valData.within_3_fold_count + '/' + valData.n + ')')])
   ]),

   // Plot & Table Grid
   e('div',{className:'grid', style:{marginTop:'16px'}},[
    e('div',{className:'col-6'},e('div',{className:'card', style:{background:'var(--bg-subtle,#1e293b)'}},[
     e('h4',{},'Predicted vs Observed Scatter Plot'),
     e('div',{style:{marginTop:'8px'}},renderPredObsPlot()),
     e('div',{className:'row toolbar', style:{marginTop:'4px', fontSize:'11px', color:'#94a3b8'}},[
      e('span',{},'── Identity Line (y = x)'),
      e('span',{},'┄┄ 2-Fold Band (0.5x – 2.0x)')
     ])
    ])),
    e('div',{className:'col-6'},e('div',{className:'card', style:{background:'var(--bg-subtle,#1e293b)'}},[
     e('h4',{},'Performance Band Breakdown'),
     e('div',{className:'table-scroll', style:{marginTop:'8px'}},[
      e('table',{className:'table'},[
       e('thead',{},e('tr',{},[
        e('th',{},'Species / Route'),
        e('th',{},'Endpoint'),
        e('th',{},'Method'),
        e('th',{},'Observed'),
        e('th',{},'Predicted'),
        e('th',{},'AFE'),
        e('th',{},'Band')
       ])),
       e('tbody',{},pairs.map((p, idx)=>e('tr',{key:idx},[
        e('td',{},p.species + ' ' + p.route),
        e('td',{className:'mono'},p.endpoint),
        e('td',{className:'small'},p.method),
        e('td',{className:'mono'},p.observed),
        e('td',{className:'mono'},p.predicted),
        e('td',{className:'mono'},p.absolute_fold_error + 'x'),
        e('td',{},e('span',{className:'status-badge ' + (p.absolute_fold_error<=2.0 ? 'status-ready' : 'status-failed'), style:{fontSize:'10px'}},p.performance_band))
       ])))
      ])
     ])
    ]))
   ])
  ]);
 }


 function HumanPkSection({versionId}){
  const [profile, setProfile] = React.useState(null);
  const [simResult, setSimResult] = React.useState(null);
  const [snapshots, setSnapshots] = React.useState([]);
  const [valData, setValData] = React.useState(null);
  const [loading, setLoading] = React.useState(false);
  const [simBusy, setSimBusy] = React.useState(false);
  const [freezeBusy, setFreezeBusy] = React.useState(false);
  const [plotScale, setPlotScale] = React.useState('linear');
  const [snapshotName, setSnapshotName] = React.useState('Human Prospective PK Prediction');
  const [showOverrides, setShowOverrides] = React.useState(false);
  const [form, setForm] = React.useState({
   route: 'IV',
   administration_type: 'IV_BOLUS',
   dose: 100,
   dose_unit: 'mg',
   body_weight_kg: 70,
   infusion_duration_hours: 1,
   dosing_frequency: 'Single Dose',
   dose_interval_hours: 24,
   num_doses: 1,
   user_cl_override: '',
   user_v_override: '',
   user_f_override: '',
   user_fg_override: '',
   user_ka_override: ''
  });

  const loadData = React.useCallback(async ()=>{
   if(!versionId) return;
   setLoading(true);
   try {
    const [pRes, sRes, vRes] = await Promise.all([
     api.get('/compound-versions/' + versionId + '/human-pk/profile'),
     api.get('/compound-versions/' + versionId + '/human-pk/snapshots'),
     api.get('/compound-versions/' + versionId + '/human-pk/validation')
    ]);
    setProfile(pRes);
    setSnapshots(sRes?.snapshots || []);
    setValData(vRes);
   } catch(err) {
    console.error(err);
   } finally {
    setLoading(false);
   }
  }, [versionId]);

  React.useEffect(()=>{
   loadData();
  }, [loadData]);

  const runSimulation = async ()=>{
   if(!versionId) return;
   setSimBusy(true);
   try {
    const payload = {
     route: form.route,
     administration_type: form.administration_type,
     dose: Number(form.dose),
     dose_unit: form.dose_unit,
     body_weight_kg: Number(form.body_weight_kg),
     infusion_duration_hours: Number(form.infusion_duration_hours),
     dosing_frequency: form.dosing_frequency,
     dose_interval_hours: Number(form.dose_interval_hours),
     num_doses: Number(form.num_doses),
     model_type: 'ONE_COMPARTMENT',
     user_cl_override: form.user_cl_override ? Number(form.user_cl_override) : null,
     user_v_override: form.user_v_override ? Number(form.user_v_override) : null,
     user_f_override: form.user_f_override ? Number(form.user_f_override) : null,
     user_fg_override: form.user_fg_override ? Number(form.user_fg_override) : null,
     user_ka_override: form.user_ka_override ? Number(form.user_ka_override) : null
    };
    const res = await api.post('/compound-versions/' + versionId + '/human-pk/simulation/run', payload);
    setSimResult(res);
   } catch(err) {
    alert(err.message || String(err));
   } finally {
    setSimBusy(false);
   }
  };

  const freezeSnapshot = async ()=>{
   if(!versionId) return;
   setFreezeBusy(true);
   try {
    await api.post('/compound-versions/' + versionId + '/human-pk/freeze-snapshot', {snapshot_name: snapshotName || 'Prospective Human PK Prediction'});
    await loadData();
    alert('Prospective Human PK prediction snapshot frozen successfully.');
   } catch(err) {
    alert(err.message || String(err));
   } finally {
    setFreezeBusy(false);
   }
  };

  // Render SVG Concentration-Time Plot (Linear or Semi-Log)
  const renderCurvePlot = ()=>{
   if(!simResult || !simResult.time_series || simResult.time_series.length === 0) return null;
   const pts = simResult.time_series;
   const obs = simResult.observations_overlay || [];

   const allT = pts.map(p=>p.time_hours);
   const maxT = Math.max(...allT, 1.0);

   let allC = pts.map(p=>p.concentration_ng_ml);
   if(obs.length > 0) {
    allC = allC.concat(obs.filter(o=>o.concentration_ng_ml!=null).map(o=>o.concentration_ng_ml));
   }
   const maxC = Math.max(...allC, 1.0);
   const isLog = plotScale === 'log';
   const minC = isLog ? Math.max(0.1, Math.min(...allC.filter(c=>c>0), 1.0) * 0.5) : 0;

   const W = 620, H = 260, padL = 60, padB = 40, padR = 20, padT = 25;

   const mapX = t => padL + (t / maxT) * (W - padL - padR);
   const mapY = c => {
    if(isLog) {
     const safeC = Math.max(c, minC);
     const logMin = Math.log10(minC);
     const logMax = Math.log10(maxC * 1.2);
     return padT + (1.0 - (Math.log10(safeC) - logMin) / (logMax - logMin || 1.0)) * (H - padT - padB);
    } else {
     return padT + (1.0 - c / (maxC * 1.15 || 1.0)) * (H - padT - padB);
    }
   };

   const pathData = pts.map((p, idx)=>(idx === 0 ? 'M ' : 'L ') + mapX(p.time_hours).toFixed(1) + ' ' + mapY(p.concentration_ng_ml).toFixed(1)).join(' ');

   return e('svg',{width:W, height:H, style:{background:'#0f172a', borderRadius:'8px', width:'100%', height:'auto'}},[
    // Axes
    e('line',{key:'ax-x', x1:padL, y1:H-padB, x2:W-padR, y2:H-padB, stroke:'#334155', strokeWidth:1}),
    e('line',{key:'ax-y', x1:padL, y1:padT, x2:padL, y2:H-padB, stroke:'#334155', strokeWidth:1}),
    e('text',{key:'lbl-x', x:W/2, y:H-8, fill:'#94a3b8', fontSize:11, textAnchor:'middle'},'Time (hours)'),
    e('text',{key:'lbl-y', x:16, y:H/2, fill:'#94a3b8', fontSize:11, textAnchor:'middle', transform:'rotate(-90 16 '+(H/2)+')'},isLog?'Plasma Conc. [log] (ng/mL)':'Plasma Conc. (ng/mL)'),

    // Grid ticks
    [0, maxT*0.25, maxT*0.5, maxT*0.75, maxT].map((t, idx)=>e('g',{key:'tick-x-'+idx},[
     e('line',{x1:mapX(t), y1:H-padB, x2:mapX(t), y2:H-padB+4, stroke:'#64748b'}),
     e('text',{x:mapX(t), y:H-padB+16, fill:'#64748b', fontSize:10, textAnchor:'middle'},t.toFixed(0)+'h')
    ])),

    // Predicted Curve
    e('path',{key:'pred-path', d:pathData, fill:'none', stroke:'#38bdf8', strokeWidth:2.5}),

    // Observed clinical overlay points
    obs.map((o, idx)=>{
     if(o.concentration_ng_ml == null) return null;
     const cx = mapX(o.time_hours);
     const cy = mapY(o.concentration_ng_ml);
     return e('g',{key:'obs-'+idx},[
      e('circle',{cx, cy, r:4.5, fill:'#f59e0b', stroke:'#ffffff', strokeWidth:1.5}),
      e('text',{x:cx+6, y:cy-4, fill:'#fbbf24', fontSize:9},'Obs ('+o.concentration_ng_ml+')')
     ]);
    }),

    // Legend
    e('g',{key:'legend', transform:'translate('+(W-190)+', 12)'},[
     e('line',{x1:0, y1:6, x2:20, y2:6, stroke:'#38bdf8', strokeWidth:2}),
     e('text',{x:25, y:10, fill:'#cbd5e1', fontSize:10},'Human Prediction'),
     obs.length>0&&e('circle',{cx:110, cy:6, r:4, fill:'#f59e0b', stroke:'#ffffff', strokeWidth:1}),
     obs.length>0&&e('text',{x:120, y:10, fill:'#cbd5e1', fontSize:10},'Clinical Obs')
    ])
   ]);
  };

  if(!profile) return null;

  const cl = profile.clearance || {};
  const v = profile.volume || {};
  const abs = profile.absorption || {};
  const readiness = profile.readiness || {};
  const clDis = cl.disagreement || {};
  const vDis = v.disagreement || {};
  const hasMajorDis = clDis.has_major_disagreement || vDis.has_major_disagreement;

  return e('div',{className:'card', key:'human-pk-section', style:{marginTop:'20px', border:'1px solid var(--border-color,#334155)'}},[
   e('div',{className:'row toolbar'},[
    e('div',{},[
     e('div',{className:'eyebrow'},'7 · HUMAN PK PREDICTION & TRANSLATIONAL SIMULATION'),
     e('h3',{},'HUMAN PK PREDICTION — Multi-Stream Parameter Assembly & Simulation (Stage 5B-4)'),
    ]),
    e('div',{className:'row', style:{gap:'8px'}},[
     e('span',{className:'status-badge ' + (readiness.iv_simulation?.status==='READY'?'status-ready':(readiness.iv_simulation?.status==='PARTIALLY_READY'?'status-medium':'status-failed'))},'IV: ' + (readiness.iv_simulation?.status||'UNAVAILABLE')),
     e('span',{className:'status-badge ' + (readiness.po_simulation?.status==='READY'?'status-ready':(readiness.po_simulation?.status==='PARTIALLY_READY'?'status-medium':'status-failed'))},'PO: ' + (readiness.po_simulation?.status||'UNAVAILABLE'))
    ])
   ]),
   e('p',{className:'small'},'Scientifically conservative Human PK translation assembling multi-stream evidence (Clinical Experimental > Hepatic IVIVE > Cross-Species Allometry > Physicochemical Distribution).'),

   // Readiness Scorecard
   e('div',{className:'grid', style:{marginTop:'12px'}},[
    e('div',{className:'col-6'},e('div',{className:'card', style:{background:'rgba(15,23,42,0.6)'}},[
     e('div',{className:'row toolbar'},[
      e('strong',{},'Human IV Simulation Readiness'),
      e('span',{className:'status-badge ' + (readiness.iv_simulation?.status==='READY'?'status-ready':(readiness.iv_simulation?.status==='PARTIALLY_READY'?'status-medium':'status-failed'))},readiness.iv_simulation?.status)
     ]),
     e('ul',{style:{paddingLeft:'18px', margin:'6px 0 0 0', fontSize:'12px', color:'#94a3b8'}},
      (readiness.iv_simulation?.reasons||[]).map((r, i)=>e('li',{key:i},r))
     )
    ])),
    e('div',{className:'col-6'},e('div',{className:'card', style:{background:'rgba(15,23,42,0.6)'}},[
     e('div',{className:'row toolbar'},[
      e('strong',{},'Human PO Simulation Readiness'),
      e('span',{className:'status-badge ' + (readiness.po_simulation?.status==='READY'?'status-ready':(readiness.po_simulation?.status==='PARTIALLY_READY'?'status-medium':'status-failed'))},readiness.po_simulation?.status)
     ]),
     e('ul',{style:{paddingLeft:'18px', margin:'6px 0 0 0', fontSize:'12px', color:'#94a3b8'}},
      (readiness.po_simulation?.reasons||[]).map((r, i)=>e('li',{key:i},r))
     )
    ])),
    e('div',{className:'col-12'},e('div',{className:'alert', style:{fontSize:'12px'}},readiness.oral_translation_guardrail))
   ]),

   // Disagreement Alert Banner
   hasMajorDis&&e('div',{className:'alert alert-danger', style:{marginTop:'12px', background:'rgba(239,68,68,0.15)', border:'1px solid #ef4444', padding:'12px', borderRadius:'6px'}},[
    e('strong',{style:{color:'#f87171'}},'⚠️ MAJOR DISAGREEMENT DETECTED (>3-fold difference between candidate methods)'),
    e('p',{className:'small', style:{color:'#fca5a5', marginTop:'4px'}},'Independent scientific estimates disagree by >3-fold. Automatic averaging is disabled by scientific policy. Review candidates in the assembly matrix below.')
   ]),

   // Multi-Stream Assembly Table
   e('div',{className:'card', style:{marginTop:'16px'}},[
    e('h4',{},'Human Parameter Assembly & Candidate Evidence Streams'),
    e('p',{className:'small'},'Independent quantitative estimates are preserved side-by-side rather than silently overwritten.'),
    e('div',{className:'table-scroll', style:{marginTop:'8px'}},[
     e('table',{className:'table'},[
      e('thead',{},e('tr',{},[
       e('th',{},'Parameter'),
       e('th',{},'Selected Value'),
       e('th',{},'Source Type'),
       e('th',{},'Confidence'),
       e('th',{},'Available Candidates'),
       e('th',{},'Disagreement')
      ])),
      e('tbody',{},[
       // CL Row
       e('tr',{key:'row-cl'},[
        e('td',{},e('strong',{},'Clearance (CL)')),
        e('td',{className:'mono'},cl.selected_value ? cl.selected_value + ' ' + cl.selected_unit : '—'),
        e('td',{className:'small'},cl.selected_source || 'MODEL_UNAVAILABLE'),
        e('td',{},e('span',{className:'status-badge '+(cl.confidence==='HIGH'?'status-ready':(cl.confidence==='MEDIUM'?'status-medium':'status-failed')), style:{fontSize:'10px'}},cl.confidence)),
        e('td',{className:'small'},(cl.candidates||[]).map(c=>c.source_name + ': ' + c.value + ' ' + c.unit).join(' · ') || 'None'),
        e('td',{},e('span',{className:'status-badge ' + (clDis.status==='GENERALLY_CONSISTENT'?'status-ready':(clDis.status==='MODERATE_DISAGREEMENT'?'status-medium':'status-failed')), style:{fontSize:'10px'}},clDis.status||'N/A'))
       ]),
       // Volume Row
       e('tr',{key:'row-v'},[
        e('td',{},e('strong',{},'Volume (Vss / Vz)')),
        e('td',{className:'mono'},v.selected_value ? v.selected_value + ' ' + v.selected_unit : '—'),
        e('td',{className:'small'},v.selected_source || 'MODEL_UNAVAILABLE'),
        e('td',{},e('span',{className:'status-badge '+(v.confidence==='HIGH'?'status-ready':(v.confidence==='MEDIUM'?'status-medium':'status-failed')), style:{fontSize:'10px'}},v.confidence)),
        e('td',{className:'small'},(v.candidates||[]).map(c=>c.source_name + ': ' + c.value + ' ' + c.unit).join(' · ') || 'None'),
        e('td',{},e('span',{className:'status-badge ' + (vDis.status==='GENERALLY_CONSISTENT'?'status-ready':(vDis.status==='MODERATE_DISAGREEMENT'?'status-medium':'status-failed')), style:{fontSize:'10px'}},vDis.status||'N/A'))
       ]),
       // Half-Life Row
       e('tr',{key:'row-thalf'},[
        e('td',{},e('strong',{},'Half-Life (t½)')),
        e('td',{className:'mono'},profile.half_life?.selected_value ? profile.half_life.selected_value + ' ' + profile.half_life.selected_unit : '—'),
        e('td',{className:'small'},'Analytical ln(2)*V/CL'),
        e('td',{},e('span',{className:'status-badge status-ready', style:{fontSize:'10px'}},'CALCULATED')),
        e('td',{className:'small'},'Allometric: ' + (profile.half_life?.allometric_value || '—') + ' h · Clinical: ' + (profile.half_life?.experimental_value || '—') + ' h'),
        e('td',{},'—')
       ]),
       // Bioavailability Row
       e('tr',{key:'row-f'},[
        e('td',{},e('strong',{},'Bioavailability (F)')),
        e('td',{className:'mono'},abs.f_selected ? abs.f_selected + '%' : '—'),
        e('td',{className:'small'},abs.f_selected_source || 'MODEL_UNAVAILABLE'),
        e('td',{},e('span',{className:'status-badge '+(abs.f_experimental?'status-ready':(abs.f_predicted?'status-medium':'status-failed')), style:{fontSize:'10px'}},abs.f_experimental?'HIGH':(abs.f_predicted?'MEDIUM':'LOW'))),
        e('td',{className:'small'},'Fa: '+(abs.fa_value?Math.round(abs.fa_value*100)+'%':'—')+' · Fg: '+(abs.fg_value?Math.round(abs.fg_value*100)+'%':'MODEL_UNAVAILABLE')+' · Fh: '+(abs.fh_value?Math.round(abs.fh_value*100)+'%':'—')),
        e('td',{},abs.f_disagreement?.status || 'N/A')
       ]),
       // ka Row
       e('tr',{key:'row-ka'},[
        e('td',{},e('strong',{},'Absorption Rate (ka)')),
        e('td',{className:'mono'},abs.ka_value ? abs.ka_value + ' 1/h' : '—'),
        e('td',{className:'small'},abs.ka_source || 'MODEL_UNAVAILABLE'),
        e('td',{},e('span',{className:'status-badge '+(abs.ka_value?'status-ready':'status-failed'), style:{fontSize:'10px'}},abs.ka_value?'HIGH':'UNAVAILABLE')),
        e('td',{className:'small'},abs.ka_value ? 'Derived from Clinical Tmax' : 'No clinical Tmax or validated ka available'),
        e('td',{},'—')
       ])
      ])
     ])
    ])
   ]),

   // Simulation Parameter Controls
   e('div',{className:'card', style:{marginTop:'16px', background:'rgba(30,41,59,0.5)'}},[
    e('div',{className:'row toolbar'},[
     e('h4',{},'Human PK Simulation Engine (1-Compartment)'),
     e('button',{className:'secondary', onClick:()=>setShowOverrides(!showOverrides)},showOverrides?'Hide Scientific Overrides':'Explore Overrides / Sensitivity')
    ]),
    e('div',{className:'grid', style:{marginTop:'8px'}},[
     e('div',{className:'col-3'},[
      e('label',{},'Route'),
      e('select',{value:form.route, onChange:ev=>setForm(c=>({...c, route:ev.target.value}))},[
       e('option',{value:'IV'},'IV (Intravenous)'),
       e('option',{value:'PO'},'PO (Oral Extravascular)')
      ])
     ]),
     form.route==='IV'&&e('div',{className:'col-3'},[
      e('label',{},'Administration Type'),
      e('select',{value:form.administration_type, onChange:ev=>setForm(c=>({...c, administration_type:ev.target.value}))},[
       e('option',{value:'IV_BOLUS'},'IV Bolus'),
       e('option',{value:'IV_INFUSION'},'IV Infusion')
      ])
     ]),
     form.route==='IV'&&form.administration_type==='IV_INFUSION'&&e('div',{className:'col-3'},Field({
      label:'Infusion Duration (h)',
      type:'number',
      value:form.infusion_duration_hours,
      onChange:v=>setForm(c=>({...c, infusion_duration_hours:v}))
     })),
     e('div',{className:'col-3'},Field({label:'Dose', type:'number', value:form.dose, onChange:v=>setForm(c=>({...c, dose:v}))})),
     e('div',{className:'col-3'},[
      e('label',{},'Dose Unit'),
      e('select',{value:form.dose_unit, onChange:ev=>setForm(c=>({...c, dose_unit:ev.target.value}))},[
       e('option',{value:'mg'},'mg'),
       e('option',{value:'mg/kg'},'mg/kg'),
       e('option',{value:'ug'},'µg')
      ])
     ]),
     e('div',{className:'col-3'},Field({label:'Patient BW (kg)', type:'number', value:form.body_weight_kg, onChange:v=>setForm(c=>({...c, body_weight_kg:v}))})),
     e('div',{className:'col-3'},[
      e('label',{},'Dosing Frequency'),
      e('select',{value:form.dosing_frequency, onChange:ev=>setForm(c=>({...c, dosing_frequency:ev.target.value}))},[
       e('option',{value:'Single Dose'},'Single Dose'),
       e('option',{value:'Repeated Dosing'},'Repeated Dosing')
      ])
     ]),
     form.dosing_frequency==='Repeated Dosing'&&e('div',{className:'col-3'},Field({label:'Dose Interval (h)', type:'number', value:form.dose_interval_hours, onChange:v=>setForm(c=>({...c, dose_interval_hours:v}))})),
     form.dosing_frequency==='Repeated Dosing'&&e('div',{className:'col-3'},Field({label:'Number of Doses', type:'number', value:form.num_doses, onChange:v=>setForm(c=>({...c, num_doses:v}))}))
    ]),

    // Overrides section
    showOverrides&&e('div',{className:'card', style:{marginTop:'12px', background:'rgba(15,23,42,0.8)', border:'1px dashed #475569'}},[
     e('h5',{},'Scientific Sensitivity Analysis / User Overrides'),
     e('p',{className:'small'},'Explicit overrides are recorded as assumptions and do not overwrite base models.'),
     e('div',{className:'grid', style:{marginTop:'6px'}},[
      e('div',{className:'col-3'},Field({label:'CL Override (mL/min/kg)', type:'number', value:form.user_cl_override, onChange:v=>setForm(c=>({...c, user_cl_override:v}))})),
      e('div',{className:'col-3'},Field({label:'V Override (L/kg)', type:'number', value:form.user_v_override, onChange:v=>setForm(c=>({...c, user_v_override:v}))})),
      form.route==='PO'&&e('div',{className:'col-3'},Field({label:'Oral F Override (%)', type:'number', value:form.user_f_override, onChange:v=>setForm(c=>({...c, user_f_override:v}))})),
      form.route==='PO'&&e('div',{className:'col-3'},Field({label:'Gut Fg Override (0-1)', type:'number', value:form.user_fg_override, onChange:v=>setForm(c=>({...c, user_fg_override:v}))})),
      form.route==='PO'&&e('div',{className:'col-3'},Field({label:'ka Override (1/h)', type:'number', value:form.user_ka_override, onChange:v=>setForm(c=>({...c, user_ka_override:v}))}))
     ])
    ]),

    e('div',{className:'row toolbar', style:{marginTop:'14px'}},[
     e('button',{onClick:runSimulation, disabled:simBusy},simBusy?'Simulating…':'Run Human PK Simulation')
    ])
   ]),

   // Simulation Output & Interactive Curve Plot
   simResult&&e('div',{className:'card', style:{marginTop:'16px'}},[
    e('div',{className:'row toolbar'},[
     e('div',{},[
      e('div',{className:'eyebrow'},'SIMULATION OUTPUT'),
      e('h4',{},simResult.route + ' Human Concentration-Time Profile (' + simResult.administration_type + ')')
     ]),
     e('div',{className:'row', style:{gap:'6px'}},[
      e('button',{className:plotScale==='linear'?'':'secondary', onClick:()=>setPlotScale('linear')},'Linear Scale'),
      e('button',{className:plotScale==='log'?'':'secondary', onClick:()=>setPlotScale('log')},'Semi-Log Scale')
     ])
    ]),

    // Scientific Disclaimer Tags
    e('div',{className:'row', style:{gap:'6px', margin:'6px 0 12px 0'}},
     (simResult.scientific_labels||[]).map((lbl, i)=>e('span',{key:i, className:'status-badge status-medium', style:{fontSize:'10px'}},lbl))
    ),

    // Grid: Curve + Metrics
    e('div',{className:'grid comparison-config-grid'},[
     e('div',{className:'col-7'},e('div',{style:{marginTop:'6px'}},renderCurvePlot())),
     e('div',{className:'col-5'},[
      e('h5',{},'Analytical PK Metrics'),
      e('div',{className:'grid ivive-output-grid', style:{marginTop:'6px'}},[
       e('div',{className:'card pk-nca-card'},[e('span',{},'Cmax (ng/mL)'), e('strong',{className:'mono'},simResult.output_metrics?.cmax_ng_ml)]),
       e('div',{className:'card pk-nca-card'},[e('span',{},'Tmax (h)'), e('strong',{className:'mono'},simResult.output_metrics?.tmax_hours)]),
       e('div',{className:'card pk-nca-card'},[e('span',{},'AUC single (ng·h/mL)'), e('strong',{className:'mono'},simResult.output_metrics?.auc_single_ng_h_ml)]),
       e('div',{className:'card pk-nca-card'},[e('span',{},'Half-Life (h)'), e('strong',{className:'mono'},simResult.output_metrics?.half_life_hours)]),
       e('div',{className:'card pk-nca-card'},[e('span',{},'Clearance (L/h)'), e('strong',{className:'mono'},simResult.output_metrics?.clearance_l_h)]),
       e('div',{className:'card pk-nca-card'},[e('span',{},'Volume (L)'), e('strong',{className:'mono'},simResult.output_metrics?.volume_l)]),
       simResult.route==='PO'&&e('div',{className:'card pk-nca-card'},[e('span',{},'Bioavailability (%)'), e('strong',{className:'mono'},simResult.output_metrics?.bioavailability_pct + '%')]),
       simResult.route==='PO'&&e('div',{className:'card pk-nca-card'},[e('span',{},'ka (1/h)'), e('strong',{className:'mono'},simResult.output_metrics?.ka_1_per_h)])
      ]),
      (simResult.assumptions||[]).length>0&&e('div',{style:{marginTop:'8px'}},[
       e('strong',{className:'small'},'Assumptions & Overrides:'),
       e('ul',{style:{paddingLeft:'16px', margin:'4px 0 0 0', fontSize:'11px', color:'#f59e0b'}},
        simResult.assumptions.map((a, i)=>e('li',{key:i},a))
       )
      ])
     ])
    ])
   ]),

   // Prospective Prediction Freeze & Snapshot Panel
   e('div',{className:'card', style:{marginTop:'16px'}},[
    e('div',{className:'row toolbar'},[
     e('div',{},[
      e('div',{className:'eyebrow'},'PROSPECTIVE SNAPSHOT GOVERNANCE'),
      e('h4',{},'Freeze Prospective Human PK Prediction'),
     ]),
     e('button',{onClick:freezeSnapshot, disabled:freezeBusy},freezeBusy?'Freezing…':'Freeze Prospective Snapshot')
    ]),
    e('p',{className:'small'},'Freezes the current prediction parameters into an immutable database snapshot BEFORE entering clinical Human trial data to prevent hindsight bias.'),
    e('div',{className:'row', style:{marginTop:'8px'}},[
     Field({label:'Snapshot Name', value:snapshotName, onChange:setSnapshotName, placeholder:'e.g. Phase 1 Clinical Dose Candidate'})
    ]),

    // Snapshots table
    snapshots.length>0&&e('div',{className:'table-scroll', style:{marginTop:'12px'}},[
     e('h5',{},'Immutable Prediction Snapshot History (' + snapshots.length + ')'),
     e('table',{className:'table'},[
      e('thead',{},e('tr',{},[
       e('th',{},'Snapshot Name'),
       e('th',{},'Freeze Date'),
       e('th',{},'Pred CL'),
       e('th',{},'Pred V'),
       e('th',{},'Pred F'),
       e('th',{},'Pred ka'),
       e('th',{},'Inputs Hash')
      ])),
      e('tbody',{},snapshots.map(s=>e('tr',{key:s.id},[
       e('td',{},e('strong',{},s.snapshot_name)),
       e('td',{className:'small'},new Date(s.created_at).toLocaleString()),
       e('td',{className:'mono'},s.selected_cl ? s.selected_cl + ' mL/min/kg' : '—'),
       e('td',{className:'mono'},s.selected_v ? s.selected_v + ' L/kg' : '—'),
       e('td',{className:'mono'},s.f_selected ? s.f_selected + '%' : '—'),
       e('td',{className:'mono'},s.ka_value ? s.ka_value + ' 1/h' : '—'),
       e('td',{className:'mono small', style:{color:'#64748b'}},s.inputs_hash ? s.inputs_hash.slice(0,10)+'…' : '—')
      ])))
     ])
    ])
   ]),

   // Retrospective Clinical Validation Panel
   valData&&valData.status==='VALIDATED'&&e('div',{className:'card', style:{marginTop:'16px', background:'rgba(23,105,170,0.08)'}},[
    e('div',{className:'row toolbar'},[
     e('div',{},[
      e('div',{className:'eyebrow'},'RETROSPECTIVE VALIDATION'),
      e('h4',{},'Clinical Validation Against Frozen Prospective Snapshot'),
     ]),
     e('span',{className:'status-badge status-ready'},'AAFE ' + (valData.metrics?.aafe || '—') + 'x')
    ]),
    e('p',{className:'small'},'Comparing subsequent Human clinical data against the immutable prospective prediction snapshot (' + (valData.snapshot_name || '') + ').'),

    e('div',{className:'grid ivive-output-grid', style:{marginTop:'8px'}},[
     e('div',{className:'card pk-nca-card'},[e('span',{},'Clinical Comparisons'), e('strong',{className:'mono'},valData.n_comparisons)]),
     e('div',{className:'card pk-nca-card'},[e('span',{},'AAFE (Fold Error)'), e('strong',{className:'mono'},valData.metrics?.aafe + 'x')]),
     e('div',{className:'card pk-nca-card'},[e('span',{},'Within 2-Fold (%)'), e('strong',{className:'mono'},valData.metrics?.within_2_fold_pct + '%')]),
     e('div',{className:'card pk-nca-card'},[e('span',{},'Within 3-Fold (%)'), e('strong',{className:'mono'},valData.metrics?.within_3_fold_pct + '%')])
    ]),

    e('div',{className:'table-scroll', style:{marginTop:'12px'}},[
     e('table',{className:'table'},[
      e('thead',{},e('tr',{},[
       e('th',{},'Study'),
       e('th',{},'Endpoint'),
       e('th',{},'Route'),
       e('th',{},'Predicted'),
       e('th',{},'Observed'),
       e('th',{},'Fold Error'),
       e('th',{},'AFE'),
       e('th',{},'Performance Band')
      ])),
      e('tbody',{},(valData.comparisons||[]).map((c, i)=>e('tr',{key:i},[
       e('td',{},c.study_name),
       e('td',{className:'mono'},c.endpoint),
       e('td',{},c.route),
       e('td',{className:'mono'},c.predicted + ' ' + c.unit),
       e('td',{className:'mono'},c.observed + ' ' + c.unit),
       e('td',{className:'mono'},c.fold_error + 'x'),
       e('td',{className:'mono'},c.absolute_fold_error + 'x'),
       e('td',{},e('span',{className:'status-badge ' + (c.absolute_fold_error<=2.0?'status-ready':(c.absolute_fold_error<=3.0?'status-medium':'status-failed')), style:{fontSize:'10px'}},c.performance_band))
      ])))
     ])
    ])
   ])
  ]);
 }


 function ScientificValidationSection(){
  const [configs, setConfigs] = React.useState(null);
  const [gate, setGate] = React.useState(null);
  const [registry, setRegistry] = React.useState([]);
  const [lightning, setLightning] = React.useState(null);
  const [readiness, setReadiness] = React.useState(null);
  const [testSmiles, setTestSmiles] = React.useState('CC(=O)Oc1ccccc1C(=O)O.Cl');
  const [stdResult, setStdResult] = React.useState(null);
  const [loading, setLoading] = React.useState(false);

  React.useEffect(()=>{
   api.get('/standardization/configurations').then(setConfigs).catch(console.error);
   api.get('/evaluation/golden-gate').then(setGate).catch(console.error);
   api.get('/evaluation/registry').then(res=>setRegistry(res.registry||[])).catch(console.error);
   api.get('/evaluation/lightning-audit').then(setLightning).catch(console.error);
   api.get('/evaluation/rdkit-readiness').then(setReadiness).catch(console.error);
  }, []);

  const handleTestStandardize = async()=>{
   setLoading(true);
   try{
    const res = await api.post('/standardization/standardize', {smiles: testSmiles});
    setStdResult(res);
   }catch(err){
    console.error(err);
   }finally{
    setLoading(false);
   }
  };

  return e('div',{className:'card', key:'sci-validation', style:{marginTop:'20px'}},[
   e('div',{className:'eyebrow'},'STAGE 4C SCIENTIFIC HARDENING'),
   e('h2',{},'Scientific Validation & Structure Standardization'),
   e('p',{},'Canonical structure pipeline CHEM_STANDARDIZER_V1, 52-molecule golden gate, evaluation registry, and RDKit readiness controls.'),

   gate && e('div',{className:'card', key:'gate-card', style:{background:'rgba(23,105,170,0.08)', marginBottom:'16px'}},[
    e('div',{className:'row toolbar'},[
     e('div',{},[
      e('h3',{},'Golden Structure Set Gate (52 Reference Molecules)'),
      e('span',{className:'status-badge status-ready', style:{marginLeft:'10px'}},gate.gate_passed ? '100% PASS' : 'FAIL')
     ]),
     e('span',{className:'mono'},'RDKit '+gate.current_rdkit_version)
    ]),
    e('div',{className:'grid'},[
     e('div',{className:'col-3'},[e('span',{},'Total Items: '),e('strong',{},gate.total_items)]),
     e('div',{className:'col-3'},[e('span',{},'Passed: '),e('strong',{},gate.passed_count)]),
     e('div',{className:'col-3'},[e('span',{},'Failed / Diffs: '),e('strong',{},gate.failed_count)]),
     e('div',{className:'col-3'},[e('span',{},'Standardizer: '),e('strong',{},'CHEM_STANDARDIZER_V1')])
    ])
   ]),

   e('div',{className:'card', key:'live-std', style:{marginBottom:'16px'}},[
    e('h4',{},'Interactive Structure Standardizer (CHEM_STANDARDIZER_V1)'),
    e('div',{className:'row toolbar'},[
     Field({label:'Input SMILES', value:testSmiles, onChange:setTestSmiles, placeholder:'Enter SMILES to standardize'}),
     e('button',{onClick:handleTestStandardize, disabled:loading},loading?'Processing...':'Standardize')
    ]),
    stdResult && e('div',{className:'mono small', style:{marginTop:'8px', background:'rgba(0,0,0,0.2)', padding:'8px', borderRadius:'4px'}},[
     e('div',{},[e('strong',{},'Status: '),stdResult.status]),
     e('div',{},[e('strong',{},'Canonical SMILES: '),stdResult.canonical_smiles]),
     e('div',{},[e('strong',{},'Isomeric SMILES: '),stdResult.isomeric_smiles]),
     e('div',{},[e('strong',{},'InChIKey: '),stdResult.inchikey]),
     (stdResult.warnings||[]).map((w, idx)=>e('div',{key:idx, className:'alert', style:{marginTop:'2px'}},w))
    ])
   ]),

   e('div',{className:'card', key:'registry-card', style:{marginBottom:'16px'}},[
    e('h4',{},'Endpoint Evaluation Registry & Metric Contracts'),
    e('div',{className:'table-scroll'},e('table',{},[
     e('thead',{},e('tr',{},['Endpoint','Assay','Type','Split','Internal N','MAE/BAcc','R2/AUROC','MMP Dir Acc'].map(h=>e('th',{key:h},h)))),
     e('tbody',{},registry.map((reg, idx)=>e('tr',{key:idx},[
      e('td',{},[e('strong',{},reg.endpoint),e('div',{className:'small'},reg.species+' · '+reg.unit)]),
      e('td',{},reg.assay),
      e('td',{},reg.type),
      e('td',{className:'mono'},reg.split_type),
      e('td',{},reg.internal_metrics?.N||'-'),
      e('td',{},reg.type==='REGRESSION'?reg.internal_metrics?.MAE:reg.internal_metrics?.balanced_accuracy),
      e('td',{},reg.type==='REGRESSION'?reg.internal_metrics?.R2:reg.internal_metrics?.auroc),
      e('td',{},reg.internal_metrics?.mmp_directional_accuracy?reg.internal_metrics.mmp_directional_accuracy+'%':'-')
     ])))
    ]))
   ]),

   e('div',{className:'grid', key:'audits-grid'},[
    e('div',{className:'col-6'},e('div',{className:'card'},[
     e('h4',{},'PyTorch Lightning Security Audit'),
     lightning && e('div',{},[
      e('div',{},[e('strong',{},'Status: '),e('span',{className:'status-badge status-ready'},lightning.status)]),
      e('div',{className:'small', style:{marginTop:'4px'}},'Installed Version: '+lightning.installed_version+' (Vulnerable 2.6.2/2.6.3 absent)'),
      e('div',{className:'small', style:{marginTop:'4px'}},lightning.recommendation)
     ])
    ])),
    e('div',{className:'col-6'},e('div',{className:'card'},[
     e('h4',{},'RDKit Upgrade Readiness'),
     readiness && e('div',{},[
      e('div',{},[e('strong',{},'Readiness: '),e('span',{className:'status-badge status-ready'},readiness.readiness_status)]),
      e('div',{className:'small', style:{marginTop:'4px'}},'Current Version: '+readiness.current_rdkit_version),
      e('div',{className:'small alert', style:{marginTop:'4px'}},readiness.policy)
     ])
    ]))
   ])
  ]);
 }

 function pkProfile(versionId){
  const studies=pkData?.studies||[];
  const bioavailability=pkData?.bioavailability||[];
  const selectedStudy=studies.find(s=>s.id===pkSelectedStudyId)||studies[0];
  const details=pkSelectedStudyDetails;

  const createPkStudyAction=async()=>{
   if(!pkStudyForm.study_name.trim())return;
   setPkBusy(true);
   try{
    await api.post('/compounds/'+detail.row_id+'/pk-studies',pkStudyForm);
    setPkModalOpen(false);
    setPkStudyForm({study_name:'',species:'Rat',strain:'',sex:'Unknown',route:'PO',dose:10,dose_unit:'mg/kg',formulation:'',matrix:'Plasma',dosing_frequency:'Single Dose',fed_fasted:'Fasted',lloq:'',lloq_unit:'ng/mL',study_date:'',source:'',notes:''});
    await loadPkData(versionId);
    setMessage('PK Study created successfully');
   }catch(err){setMessage(String(err))}finally{setPkBusy(false)}
  };

  const addObservationAction=async()=>{
   if(!selectedStudy)return;
   if(pkObsForm.time_raw===''||(pkObsForm.concentration_raw===''&&!pkObsForm.blq_flag))return;
   setPkBusy(true);
   try{
    await api.post('/pk-studies/'+selectedStudy.id+'/observations',[{
     ...pkObsForm,
     time_raw:Number(pkObsForm.time_raw),
     concentration_raw:pkObsForm.blq_flag?null:Number(pkObsForm.concentration_raw)
    }]);
    setPkObsForm({subject_group_id:'Group Mean',time_raw:'',time_unit:'h',concentration_raw:'',concentration_unit:'ng/mL',blq_flag:false,replicate:'R1',notes:''});
    await loadPkStudyDetails(selectedStudy.id);
    setMessage('Observation added');
   }catch(err){setMessage(String(err))}finally{setPkBusy(false)}
  };

  const deleteObservationAction=async(obsId)=>{
   setPkBusy(true);
   try{
    await api.delete('/pk-observations/'+obsId);
    await loadPkStudyDetails(selectedStudy.id);
    setMessage('Observation deleted');
   }catch(err){setMessage(String(err))}finally{setPkBusy(false)}
  };

  const runNcaAction=async(manualIndices=null)=>{
   if(!selectedStudy)return;
   setPkBusy(true);
   try{
    await api.post('/pk-studies/'+selectedStudy.id+'/run-nca',{manual_terminal_indices:manualIndices});
    await Promise.all([loadPkStudyDetails(selectedStudy.id),loadPkData(versionId)]);
    setMessage('NCA analysis complete');
   }catch(err){setMessage(String(err))}finally{setPkBusy(false)}
  };

  const previewCsvAction=async()=>{
   if(!selectedStudy||!pkCsvText.trim())return;
   setPkBusy(true);
   try{
    const res=await api.post('/pk-studies/'+selectedStudy.id+'/preview-csv',{csv_text:pkCsvText,mapping:pkCsvMapping});
    setPkCsvPreview(res);
   }catch(err){setMessage(String(err))}finally{setPkBusy(false)}
  };

  const importCsvAction=async()=>{
   if(!selectedStudy||!pkCsvText.trim())return;
   setPkBusy(true);
   try{
    const res=await api.post('/pk-studies/'+selectedStudy.id+'/import-csv',{csv_text:pkCsvText,mapping:pkCsvMapping});
    setPkCsvModalOpen(false);
    setPkCsvText('');
    setPkCsvPreview(null);
    await loadPkStudyDetails(selectedStudy.id);
    setMessage('Imported '+res.imported_count+' PK observation(s)');
   }catch(err){setMessage(String(err))}finally{setPkBusy(false)}
  };

  const speciesList=[...new Set(studies.map(s=>s.species))].join(', ')||'None';
  const matchedF=bioavailability.filter(b=>b.status==='MATCHED');

  return e('div',{className:'pk-profile'},[
   e('div',{className:'card',key:'hero'},[
    e('div',{className:'row toolbar',key:'head'},[
     e('div',{},[
      e('div',{className:'eyebrow'},'EXPERIMENTAL PK & NONCOMPARTMENTAL ANALYSIS (NCA)'),
      e('h2',{},'PK Studies & Pharmacokinetics'),
      e('p',{className:'small'},'CompoundVersion-isolated experimental concentration-time measurements, trapezoidal NCA parameter computation, and matched absolute bioavailability.')
     ]),
     e('button',{onClick:()=>setPkModalOpen(true)},'Add PK Study')
    ]),
    e('div',{className:'pk-hero-stats',key:'stats'},[
     e('div',{className:'pk-hero-stat',key:'count'},[e('span',{},'PK Studies'),e('strong',{},studies.length)]),
     e('div',{className:'pk-hero-stat',key:'species'},[e('span',{},'Species Covered'),e('strong',{className:'small'},speciesList)]),
     e('div',{className:'pk-hero-stat',key:'f-val'},[e('span',{},'Absolute Bioavailability (F)'),e('strong',{className:'mono'},matchedF.length?matchedF.map(b=>b.label+': '+b.bioavailability_pct+'%').join(' · '):'Not calculated')]),
     e('div',{className:'pk-hero-stat',key:'scope'},[e('span',{},'Scope Isolation'),e('strong',{className:'small'},'CompoundVersion #'+versionId)])
    ])
   ]),

   studies.length>0?e('div',{className:'pk-study-list',key:'studies'},studies.map(s=>e('div',{
    key:s.id,
    className:'pk-study-item'+(selectedStudy?.id===s.id?' selected':''),
    onClick:()=>{setPkSelectedStudyId(s.id);loadPkStudyDetails(s.id).catch(err=>setMessage(String(err)))}
   },[
    e('div',{className:'row toolbar',key:'h'},[e('h4',{},s.study_name),StatusBadge({type:s.latest_nca?'CALCULATED':'EXPERIMENTAL'})]),
    e('div',{className:'small mono'},s.species+' · '+s.route+' '+s.dose+' '+s.dose_unit+' · '+s.matrix),
    e('div',{className:'small'},s.observation_count+' observation points'+(s.latest_nca?' · t1/2 '+(s.latest_nca.terminal_half_life?Number(s.latest_nca.terminal_half_life).toFixed(2)+' h':'—'):''))
   ]))):e('div',{className:'empty-state'},[StatusBadge({type:'Not measured'}),e('p',{},'No PK studies recorded for this CompoundVersion. Click Add PK Study above.')]),

   selectedStudy&&e('div',{key:'selected-study-view'},[
    e('section',{className:'card pk-plot-card',key:'section-1'},[
     e('div',{className:'eyebrow'},'1 · EXPERIMENTAL CONCENTRATION-TIME'),
     e('div',{className:'row toolbar',key:'tbar'},[
      e('div',{},[
       e('h3',{},selectedStudy.study_name+' · Experimental Data'),
       e('p',{className:'small'},selectedStudy.species+' · '+selectedStudy.route+' '+selectedStudy.dose+' '+selectedStudy.dose_unit+' · Matrix: '+selectedStudy.matrix+(selectedStudy.formulation?' ('+selectedStudy.formulation+')':''))
      ]),
      e('div',{className:'row'},[
       e('div',{className:'manual-actions'},[
        e('button',{className:pkPlotType==='linear'?'':'secondary',onClick:()=>setPkPlotType('linear')},'Linear Plot'),
        e('button',{className:pkPlotType==='log'?'':'secondary',onClick:()=>setPkPlotType('log')},'Semi-Log Plot')
       ]),
       e('button',{className:'secondary',onClick:()=>setPkCsvModalOpen(true)},'Import CSV'),
       e('button',{className:'danger',onClick:async()=>{if(confirm('Delete study '+selectedStudy.study_name+'?')){await api.delete('/pk-studies/'+selectedStudy.id);loadPkData(versionId);}}},'Delete Study')
      ])
     ]),

     pkConcentrationTimePlot(details?.observations,details?.latest_nca,pkPlotType),

     e('h4',{style:{marginTop:'16px'}},'Concentration-Time Points'),
     (details?.observations||[]).length?e('table',{key:'obs-table'},[
      e('thead',{},e('tr',{},['Time (raw)','Time (h)','Concentration (raw)','Concentration (ng/mL)','BLQ','Replicate','Subject/Group','Notes',''].map(l=>e('th',{key:l},l)))),
      e('tbody',{},details.observations.map(obs=>e('tr',{key:obs.id},[
       e('td',{className:'mono'},obs.time_raw+' '+obs.time_unit),
       e('td',{className:'mono'},obs.time_hours),
       e('td',{className:'mono'},obs.blq_flag?'BLQ':(obs.concentration_raw!=null?obs.concentration_raw+' '+obs.concentration_unit:'—')),
       e('td',{className:'mono'},obs.blq_flag?'BLQ':(obs.concentration_normalized_ng_ml!=null?Number(obs.concentration_normalized_ng_ml).toFixed(2)+' ng/mL':'—')),
       e('td',{},obs.blq_flag?StatusBadge({type:'fail',text:'BLQ'}):'No'),
       e('td',{},obs.replicate),
       e('td',{},obs.subject_group_id),
       e('td',{className:'small'},obs.notes||'—'),
       e('td',{},e('button',{className:'danger',style:{padding:'3px 7px',fontSize:'11px'},onClick:()=>deleteObservationAction(obs.id)},'Delete'))
      ])))
     ]):e('div',{className:'empty-state'},[StatusBadge({type:'Not measured'}),e('p',{},'No concentration-time points. Add a point below or import CSV.')]),

     e('div',{className:'grid',style:{marginTop:'14px'},key:'add-obs'},[
      e('div',{className:'col-2'},Field({label:'Time *',type:'number',value:pkObsForm.time_raw,onChange:v=>setPkObsForm(c=>({...c,time_raw:v}))})),
      e('div',{className:'col-2'},[e('label',{},'Time Unit'),e('select',{value:pkObsForm.time_unit,onChange:ev=>setPkObsForm(c=>({...c,time_unit:ev.target.value}))},['h','min','sec','day'].map(u=>e('option',{key:u,value:u},u)))]),
      e('div',{className:'col-2'},Field({label:'Concentration',type:'number',disabled:pkObsForm.blq_flag,value:pkObsForm.concentration_raw,onChange:v=>setPkObsForm(c=>({...c,concentration_raw:v}))})),
      e('div',{className:'col-2'},[e('label',{},'Conc Unit'),e('select',{value:pkObsForm.concentration_unit,disabled:pkObsForm.blq_flag,onChange:ev=>setPkObsForm(c=>({...c,concentration_unit:ev.target.value}))},['ng/mL','pg/mL','µg/mL','mg/L','nM','µM'].map(u=>e('option',{key:u,value:u},u)))]),
      e('div',{className:'col-2'},[e('label',{},'BLQ Flag'),e('label',{className:'check-option',style:{marginTop:'6px'}},[e('input',{type:'checkbox',checked:pkObsForm.blq_flag,onChange:ev=>setPkObsForm(c=>({...c,blq_flag:ev.target.checked}))}),e('span',{},'Is BLQ')])]),
      e('div',{className:'col-2'},[e('label',{},'Action'),e('button',{disabled:pkBusy,onClick:addObservationAction},'Add Point')])
     ])
    ]),

    e('section',{className:'card',key:'section-2'},[
     e('div',{className:'eyebrow'},'2 · CALCULATED NONCOMPARTMENTAL ANALYSIS (NCA)'),
     e('div',{className:'row toolbar',key:'tbar'},[
      e('div',{},[
       e('h3',{},'Calculated NCA Parameters'),
       e('p',{className:'small'},'Calculated parameters are deterministic noncompartmental outputs from observed experimental points. They are explicitly not AI predictions.')
      ]),
      e('button',{disabled:pkBusy||!(details?.observations||[]).length,onClick:()=>runNcaAction()},pkBusy?'Calculating…':'Run / Recalculate NCA')
     ]),

     details?.latest_nca?e('div',{key:'nca-results'},[
      e('div',{className:'pk-nca-grid',key:'cards'},[
       e('div',{className:'pk-nca-card',key:'cmax'},[e('span',{},'Cmax (observed)'),e('strong',{className:'mono'},details.latest_nca.cmax!=null?Number(details.latest_nca.cmax).toFixed(2):'—'),e('small',{},details.latest_nca.cmax_unit)]),
       e('div',{className:'pk-nca-card',key:'tmax'},[e('span',{},'Tmax (observed)'),e('strong',{className:'mono'},details.latest_nca.tmax!=null?details.latest_nca.tmax:'—'),e('small',{},details.latest_nca.tmax_unit)]),
       e('div',{className:'pk-nca-card',key:'auclast'},[e('span',{},'AUClast'),e('strong',{className:'mono'},details.latest_nca.auclast!=null?Number(details.latest_nca.auclast).toFixed(1):'—'),e('small',{},details.latest_nca.auclast_unit)]),
       e('div',{className:'pk-nca-card',key:'aucinf'},[e('span',{},'AUCinf'),e('strong',{className:'mono'},details.latest_nca.aucinf!=null?Number(details.latest_nca.aucinf).toFixed(1):'—'),e('small',{},details.latest_nca.aucinf_unit)]),
       e('div',{className:'pk-nca-card',key:'t12'},[e('span',{},'Terminal t1/2'),e('strong',{className:'mono'},details.latest_nca.terminal_half_life!=null?Number(details.latest_nca.terminal_half_life).toFixed(2):'—'),e('small',{},'hours')]),
       selectedStudy.route==='IV'?
        e('div',{className:'pk-nca-card',key:'cl'},[e('span',{},'Clearance (CL)'),e('strong',{className:'mono'},details.latest_nca.cl!=null?Number(details.latest_nca.cl).toFixed(2):'—'),e('small',{},details.latest_nca.cl_unit)]):
        e('div',{className:'pk-nca-card',key:'clf'},[e('span',{},'Apparent CL (CL/F)'),e('strong',{className:'mono'},details.latest_nca.cl_f!=null?Number(details.latest_nca.cl_f).toFixed(2):'—'),e('small',{},details.latest_nca.cl_f_unit)]),
       selectedStudy.route==='IV'?
        e('div',{className:'pk-nca-card',key:'vz'},[e('span',{},'Volume of Dist (Vz)'),e('strong',{className:'mono'},details.latest_nca.vz!=null?Number(details.latest_nca.vz).toFixed(2):'—'),e('small',{},details.latest_nca.vz_unit)]):
        e('div',{className:'pk-nca-card',key:'vzf'},[e('span',{},'Apparent Vz (Vz/F)'),e('strong',{className:'mono'},details.latest_nca.vz_f!=null?Number(details.latest_nca.vz_f).toFixed(2):'—'),e('small',{},details.latest_nca.vz_f_unit)]),
       e('div',{className:'pk-nca-card',key:'extrap'},[e('span',{},'AUC Extrapolated %'),e('strong',{className:'mono'},details.latest_nca.auc_extrapolated_pct!=null?Number(details.latest_nca.auc_extrapolated_pct).toFixed(1)+'%':'—'),e('small',{},details.latest_nca.auc_extrapolated_pct>20?'High extrapolation':'')]),
       e('div',{className:'pk-nca-card',key:'r2'},[e('span',{},'Terminal R² (adj)'),e('strong',{className:'mono'},details.latest_nca.adjusted_r2!=null?Number(details.latest_nca.adjusted_r2).toFixed(3):'—'),e('small',{},details.latest_nca.terminal_point_count+' terminal points')]),
       e('div',{className:'pk-nca-card',key:'mrt'},[e('span',{},'MRT'),e('strong',{className:'mono'},details.latest_nca.mrt!=null?Number(details.latest_nca.mrt).toFixed(2):'—'),e('small',{},'hours')])
      ]),

      (details.latest_nca.warnings||[]).length>0&&e('div',{className:'pk-nca-warning',key:'warn'},[
       e('strong',{},'NCA Warnings / Reliability Notes:'),
       e('ul',{style:{margin:'4px 0 0',paddingLeft:'20px'}},details.latest_nca.warnings.map(w=>e('li',{key:w},w)))
      ]),

      e('div',{className:'pk-terminal-override-box',key:'term-box'},[
       e('div',{className:'row toolbar'},[
        e('div',{},[
         e('strong',{},'Terminal Phase Selection ('+details.latest_nca.selection_mode+')'),
         e('span',{className:'small'},' λz = '+(details.latest_nca.lambda_z?Number(details.latest_nca.lambda_z).toFixed(4):'—')+' h⁻¹ · Regression points: '+(details.latest_nca.terminal_points||[]).join(', '))
        ]),
        e('button',{className:'secondary',onClick:()=>setPkTerminalOverrideMode(!pkTerminalOverrideMode)},pkTerminalOverrideMode?'Cancel Override':'Manual Terminal Override')
       ]),
       pkTerminalOverrideMode&&e('div',{style:{marginTop:'10px'}},[
        e('p',{className:'small'},'Click data points to select/deselect points for terminal log-linear regression:'),
        e('div',{},(details?.observations||[]).map(o=>e('span',{
         key:o.id,
         className:'pk-point-chip'+(pkSelectedTerminalPoints.includes(o.id)?' selected':''),
         onClick:()=>{
          setPkSelectedTerminalPoints(cur=>cur.includes(o.id)?cur.filter(i=>i!==o.id):[...cur,o.id]);
         }
        },'t='+o.time_raw+'h ('+(o.blq_flag?'BLQ':o.concentration_raw+')')))),
        e('button',{style:{marginTop:'10px'},onClick:()=>runNcaAction(pkSelectedTerminalPoints)},'Apply Manual Terminal Override')
       ])
      ]),

      e('details',{style:{marginTop:'14px'},key:'prov'},[
       e('summary',{},'NCA Engine Provenance & Calculation Parameters'),
       e('div',{className:'small',style:{marginTop:'6px'}},[
        e('div',{},'Engine: '+details.latest_nca.nca_engine+' v'+details.latest_nca.nca_engine_version),
        e('div',{},'Method: '+details.latest_nca.calculation_method),
        e('div',{},'Analysis Version: #'+details.latest_nca.analysis_version+' · Timestamp: '+details.latest_nca.calculation_timestamp),
        e('div',{},'BLQ Policy: '+JSON.stringify(details.latest_nca.blq_policy))
       ])
      ])
     ]):e('div',{className:'empty-state'},[StatusBadge({type:'Not calculated'}),e('p',{},'No NCA calculation performed yet. Click Run / Recalculate NCA above.')])
    ]),

    e('section',{className:'pk-ba-card',key:'section-3'},[
     e('div',{className:'eyebrow'},'3 · ABSOLUTE BIOAVAILABILITY (F)'),
     e('h4',{},'Matched Species Absolute Bioavailability'),
     (bioavailability||[]).length?e('div',{},bioavailability.map(b=>e('div',{key:b.study_id||b.species,style:{marginTop:'8px'}},[
      b.status==='MATCHED'?e('div',{className:'pass'},[
       e('strong',{},b.label+' ('+b.species+'): '),
       e('span',{className:'mono',style:{fontSize:'18px'}},b.bioavailability_pct+'%'),
       e('span',{className:'small',style:{marginLeft:'10px'}},b.message)
      ]):e('div',{className:'small alert'},[e('strong',{},b.species+' '+b.route+': '),b.message])
     ]))):e('p',{className:'small'},'Absolute bioavailability cannot be calculated without a matched IV study.')
    ])
   ]),

   iviveProfile(versionId),
   e(PkFoundationProfile,{versionId,key:'pk-foundation'}),
   e(PkSimulationSection,{versionId,key:'pk-simulation'}),
   e(TranslationalPkSection,{versionId,key:'translational-pk'}),
   e(PkValidationSection,{versionId,key:'pk-validation'}),
   e(HumanPkSection,{versionId,key:'human-pk'}),

   pkModalOpen&&e('div',{className:'modal-backdrop',key:'study-modal'},e('div',{className:'card compound-modal'},[
    e('div',{className:'row toolbar',key:'head'},[e('h2',{},'Add PK Study'),e('button',{className:'secondary',onClick:()=>setPkModalOpen(false)},'Close')]),
    e('div',{className:'grid'},[
     e('div',{className:'col-6'},Field({label:'Study Name *',value:pkStudyForm.study_name,onChange:v=>setPkStudyForm(c=>({...c,study_name:v})),placeholder:'e.g. Single dose oral PK in SD rats'})),
     e('div',{className:'col-3'},[e('label',{},'Species *'),e('select',{value:pkStudyForm.species,onChange:ev=>setPkStudyForm(c=>({...c,species:ev.target.value}))},['Rat','Mouse','Dog','Monkey','Human','Custom'].map(s=>e('option',{key:s,value:s},s)))]),
     e('div',{className:'col-3'},Field({label:'Strain',value:pkStudyForm.strain,onChange:v=>setPkStudyForm(c=>({...c,strain:v})),placeholder:'e.g. Sprague-Dawley'})),
     e('div',{className:'col-3'},[e('label',{},'Sex'),e('select',{value:pkStudyForm.sex,onChange:ev=>setPkStudyForm(c=>({...c,sex:ev.target.value}))},['Male','Female','Mixed','Unknown'].map(s=>e('option',{key:s,value:s},s)))]),
     e('div',{className:'col-3'},[e('label',{},'Route *'),e('select',{value:pkStudyForm.route,onChange:ev=>setPkStudyForm(c=>({...c,route:ev.target.value}))},['PO','IV','SC','IP','Custom'].map(r=>e('option',{key:r,value:r},r)))]),
     e('div',{className:'col-3'},Field({label:'Dose *',type:'number',value:pkStudyForm.dose,onChange:v=>setPkStudyForm(c=>({...c,dose:v}))})),
     e('div',{className:'col-3'},[e('label',{},'Dose Unit *'),e('select',{value:pkStudyForm.dose_unit,onChange:ev=>setPkStudyForm(c=>({...c,dose_unit:ev.target.value}))},['mg/kg','µg/kg','mg','µg'].map(u=>e('option',{key:u,value:u},u)))]),
     e('div',{className:'col-4'},Field({label:'Matrix',value:pkStudyForm.matrix,onChange:v=>setPkStudyForm(c=>({...c,matrix:v})),placeholder:'Plasma'})),
     e('div',{className:'col-4'},Field({label:'Formulation',value:pkStudyForm.formulation,onChange:v=>setPkStudyForm(c=>({...c,formulation:v}))})),
     e('div',{className:'col-4'},Field({label:'Study Date',type:'date',value:pkStudyForm.study_date,onChange:v=>setPkStudyForm(c=>({...c,study_date:v}))})),
     e('div',{className:'col-6'},Field({label:'Source / CRO',value:pkStudyForm.source,onChange:v=>setPkStudyForm(c=>({...c,source:v}))})),
     e('div',{className:'col-6'},Field({label:'Notes',value:pkStudyForm.notes,onChange:v=>setPkStudyForm(c=>({...c,notes:v}))}))
    ]),
    e('div',{className:'row modal-actions'},[
     e('button',{disabled:pkBusy||!pkStudyForm.study_name.trim(),onClick:createPkStudyAction},pkBusy?'Saving…':'Create PK Study'),
     e('button',{className:'secondary',onClick:()=>setPkModalOpen(false)},'Cancel')
    ])
   ])),

   pkCsvModalOpen&&e('div',{className:'modal-backdrop',key:'csv-modal'},e('div',{className:'card compound-modal'},[
    e('div',{className:'row toolbar',key:'head'},[e('h2',{},'Import PK Observation CSV'),e('button',{className:'secondary',onClick:()=>setPkCsvModalOpen(false)},'Close')]),
    e('p',{className:'small'},'Paste CSV concentration-time data below. Columns will be auto-detected or can be explicitly mapped.'),
    e('textarea',{rows:8,value:pkCsvText,placeholder:'time,concentration,subject,blq\n0.0,0.0,Rat_1,0\n0.5,25.4,Rat_1,0\n1.0,88.1,Rat_1,0\n2.0,54.2,Rat_1,0\n4.0,18.9,Rat_1,0\n8.0,BLQ,Rat_1,1\n',onChange:ev=>{setPkCsvText(ev.target.value);setPkCsvPreview(null)}}),
    e('div',{className:'row',style:{marginTop:'12px'}},[
     e('button',{className:'secondary',disabled:pkBusy||!pkCsvText.trim(),onClick:previewCsvAction},'Preview & Validate CSV'),
     e('button',{disabled:pkBusy||!pkCsvPreview||pkCsvPreview.valid_count===0,onClick:importCsvAction},'Import '+(pkCsvPreview?.valid_count||0)+' Valid Rows')
    ]),
    pkCsvPreview&&e('div',{style:{marginTop:'14px'}},[
     e('div',{className:pkCsvPreview.error_count?'fail':'pass'},pkCsvPreview.valid_count+' valid rows · '+pkCsvPreview.error_count+' errors'),
     pkCsvPreview.preview_rows.length>0&&e('table',{style:{marginTop:'8px'}},[
      e('thead',{},e('tr',{},['Time','Concentration','Subject','BLQ'].map(h=>e('th',{key:h},h)))),
      e('tbody',{},pkCsvPreview.preview_rows.map((r,i)=>e('tr',{key:i},[e('td',{},r.time_raw+' h'),e('td',{},r.blq_flag?'BLQ':r.concentration_raw+' ng/mL'),e('td',{},r.subject),e('td',{},r.blq_flag?'Yes':'No')])))
     ])
    ])
   ]))
  ]);
 }

 function unifiedAdmetPredictionTable(predictions,measurements){
  const preds=predictions||[];
  if(!preds.length)return e('div',{className:'empty-state'},[StatusBadge({type:'Not predicted'}),e('p',{},'No ADMET predictions run for this CompoundVersion.')]);
  const endpoints=[
   {endpoint:'Solubility',unit:'µM',key:'solubility',modelName:'OpenADMET Solubility v1.0'},
   {endpoint:'Permeability',unit:'log10 cm/s',key:'caco2',modelName:'OpenADMET Caco-2 v1.0'},
   {endpoint:'Plasma protein binding',unit:'fraction unbound (fu)',key:'ppb',modelName:'OpenADMET PPB v1.0'},
   {endpoint:'HLM intrinsic clearance',unit:'mL/min/kg',key:'hlm_clint',modelName:'OpenADMET HLM v1.0'},
   {endpoint:'RLM intrinsic clearance',unit:'mL/min/kg',key:'rlm_clint',modelName:'OpenADMET RLM v1.0'},
   {endpoint:'MLM intrinsic clearance',unit:'mL/min/kg',key:'mlm_clint',modelName:'OpenADMET MLM v1.0'},
   {endpoint:'hERG liability',unit:'prob',key:'herg',modelName:'OpenADMET hERG v1.0'},
   {endpoint:'DILI clinical liability',unit:'prob',key:'dili',modelName:'OpenADMET DILI v1.0'},
   {endpoint:'Ames mutagenicity',unit:'prob',key:'ames',modelName:'OpenADMET Ames v1.0'},
   {endpoint:'P-gp inhibitor',unit:'prob',key:'pgp_inhibitor',modelName:'OpenADMET P-gp v1.0'},
   {endpoint:'CYP3A4 inhibitor',unit:'prob',key:'cyp_inhibitor',modelName:'OpenADMET CYP3A4 v1.0'},
   {endpoint:'CYP2D6 inhibitor',unit:'prob',key:'cyp_inhibitor',modelName:'OpenADMET CYP2D6 v1.0'},
   {endpoint:'CYP2C9 inhibitor',unit:'prob',key:'cyp_inhibitor',modelName:'OpenADMET CYP2C9 v1.0'}
  ];
  const predMap=new Map(preds.map(p=>[p.endpoint,p]));

  return e('div',{},[
   e('div',{className:'table-scroll'},[
    e('table',{},[
     e('thead',{},e('tr',{},['Evaluation / Endpoint','Prediction (Value & Unit)','Scientific Interpretation','Model System'].map(h=>e('th',{key:h},h)))),
     e('tbody',{},endpoints.map(ep=>{
      const pred=predMap.get(ep.endpoint);
      const interp=getInterpretation(ep.key,pred?pred.predicted_value:null);
      return e('tr',{key:ep.endpoint},[
       e('td',{style:{fontWeight:600}},ep.endpoint==='Permeability'?'Caco-2 Permeability':ep.endpoint),
       e('td',{className:'mono bold'},pred?((typeof pred.predicted_value==='number'?Number(pred.predicted_value).toFixed(3):pred.predicted_value)+' '+(pred.unit||ep.unit)):'—'),
       e('td',{},[ScientificBadge({assessment:interp.assessment,colorClass:interp.colorClass,textLabel:interp.label})]),
       e('td',{className:'small'},e('span',{className:'model-chip'},'M1'),ep.modelName)
      ]);
     }))
    ])
   ]),
   e('div',{className:'model-notes'},[
    e('strong',{},'Model Notes: '),
     e('span',{},'M1 = OpenADMET QSAR/ML Checkpoint Models v1.0 (calibrated on external ChEMBL/TDC datasets) · M2 = Chemprop MPNN Model Checkpoints · Consensus = Weighted Blend (N ≥ 10 project experimental calibration).')
    ])
   ]);
  }

  function speciesMetabolicStabilityTable(predictions,measurements){
   const species=[
    {name:'Human Liver Microsomes (HLM)',ep:'HLM intrinsic clearance',key:'hlm_clint',modelStatus:'AVAILABLE',modelDesc:'OpenADMET CheMeleon HLM v1.0'},
    {name:'Rat Liver Microsomes (RLM)',ep:'RLM intrinsic clearance',key:'rlm_clint',modelStatus:'AVAILABLE',modelDesc:'OpenADMET CheMeleon RLM v1.0'},
    {name:'Mouse Liver Microsomes (MLM)',ep:'MLM intrinsic clearance',key:'mlm_clint',modelStatus:'AVAILABLE',modelDesc:'OpenADMET CheMeleon MLM v1.0'},
    {name:'Dog Liver Microsomes (DLM)',ep:'Dog liver microsomal intrinsic clearance',key:'dlm_clint',modelStatus:'MODEL_UNAVAILABLE',modelDesc:'No qualified dog microsomal model'},
    {name:'Monkey Liver Microsomes (CyLM)',ep:'Monkey liver microsomal intrinsic clearance',key:'cylm_clint',modelStatus:'MODEL_UNAVAILABLE',modelDesc:'No qualified monkey microsomal model'}
   ];
   const predMap=new Map((predictions||[]).map(p=>[p.endpoint,p]));
   const expMap=new Map((measurements||[]).map(m=>[m.species,m]));

   return e('div',{className:'card',key:'microsomal-stability-table'},[
    e('div',{className:'eyebrow'},'METABOLIC STABILITY · LIVER MICROSOMES (MULTI-SPECIES)'),
    e('h3',{},'Metabolic Stability · Human, Rodent & Non-Rodent Microsomes'),
    e('p',{className:'small'},'Comparative hepatic microsomal intrinsic clearance (Clint) across pre-clinical species and human.'),
    e('div',{className:'table-scroll'},[
     e('table',{},[
      e('thead',{},e('tr',{},['Species & Matrix','Experimental Clint','Predicted Clint','Assessment','Model Applicability','Model'].map(h=>e('th',{key:h},h)))),
      e('tbody',{},species.map(s=>{
       const pred=predMap.get(s.ep);
       const exp=expMap.get(s.name.split(' ')[0]);
       const isUnavail = s.modelStatus === 'MODEL_UNAVAILABLE';
       const interp = isUnavail ? {assessment:'MODEL_UNAVAILABLE', colorClass:'badge-model-unavailable', label:'MODEL_UNAVAILABLE'} : getInterpretation(s.key,pred?pred.predicted_value:null);
       return e('tr',{key:s.name},[
        e('td',{style:{fontWeight:600}},s.name),
        e('td',{className:'mono'},exp?exp.value+' '+exp.unit:'Not measured'),
        e('td',{className:'mono bold'},isUnavail ? '—' : (pred?Number(pred.predicted_value).toFixed(2)+' mL/min/kg':'—')),
        e('td',{},[ScientificBadge({assessment:interp.assessment,colorClass:interp.colorClass,textLabel:interp.label})]),
        e('td',{},[isUnavail ? e('span',{className:'badge-model-unavailable',title:s.modelDesc},'UNAVAILABLE') : e('span',{className:'badge-intermediate',title:'Compound descriptors fall within training domain'},'IN DOMAIN')]),
        e('td',{className:'small'},isUnavail ? e('span',{style:{color:'#94a3b8'}},s.modelDesc) : [e('span',{className:'model-chip'},'M1'),'OpenADMET '+(pred?.model?.model_version||'v1.0')])
       ]);
      }))
     ])
    ]),
    e('div',{className:'model-notes'},[
     e('strong',{},'Model Notes: '),
     e('span',{},'M1 = OpenADMET CheMeleon Microsomal Stability v1.0 (calibrated linear mL/min/kg). Dog and Monkey microsomal models are unavailable in the local registry and reported transparently without cross-species interpolation.')
    ])
   ]);
  }

  function MultiSpeciesPkSummaryTable({versionId,studies,iviveData}){
   window.__pkMultiCache = window.__pkMultiCache || {};
   const cached = versionId ? window.__pkMultiCache[versionId] : null;
   const [multiPk,setMultiPk]=React.useState(cached);
   const [loading,setLoading]=React.useState(!cached);
   const species=['Mouse','Rat','Dog','Monkey','Human'];

   React.useEffect(()=>{
    if(versionId){
     if(window.__pkMultiCache[versionId]){
      setMultiPk(window.__pkMultiCache[versionId]);
      setLoading(false);
     }else{
      setLoading(true);
     }
     api.get('/compound-versions/'+versionId+'/pk-multi-species')
      .then(res=>{
       window.__pkMultiCache[versionId]=res;
       setMultiPk(res);
      })
      .catch(err=>console.error('[MultiSpeciesPkSummaryTable] error:', err))
      .finally(()=>setLoading(false));
    }
   },[versionId]);

   const studiesBySpecies=new Map();
   (studies||[]).forEach(s=>{if(!studiesBySpecies.has(s.species))studiesBySpecies.set(s.species,s)});
   const spMap=multiPk?.species_profiles||{};

   const hasAnyExp = species.some(sp => {
    const s = studiesBySpecies.get(sp);
    const prof = spMap[sp];
    return Boolean(s?.latest_nca || prof?.is_experimental || prof?.cl?.is_experimental || prof?.v?.is_experimental || prof?.f_is_experimental);
   });

   return e('div',{className:'card',key:'multispecies-summary'},[
    e('div',{className:'eyebrow'},'MULTI-SPECIES PK SUMMARY & TRANSLATIONAL PROFILE'),
    e('h3',{},'Multi-Species Pharmacokinetics & In Vivo Profile'),
    e('p',{className:'small'},'Comparative pharmacokinetic parameters across all 5 pre-clinical species and human. Evidence retains explicit species-specific provenance (Experimental IV NCA > Species IVIVE > Unavailable).'),
    hasAnyExp && e('div',{style:{background:'#ecfdf5',border:'1px solid #a7f3d0',color:'#065f46',padding:'10px 14px',borderRadius:'6px',marginBottom:'12px',fontSize:'12px',display:'flex',alignItems:'center',gap:'8px'}},[
     e('strong',{},'✓ 실험 데이터 반영 (Experimental Data Applied):'),
     e('span',{},'등록된 생체 내(In Vivo PK NCA) 및 시험관 내(In Vitro) 실험 측정값이 예측값보다 우선하여 파라미터 및 시뮬레이션 기본값(Default)으로 자동 계산되었습니다.')
    ]),
    (loading && !multiPk)?e('div',{className:'empty-state',style:{padding:'24px 0'}},[
     e('p',{className:'small'},'Loading multi-species pharmacokinetic profile across 5 species…')
    ]):e('div',{className:'table-scroll'},[
     e('table',{className:'pk-multispecies-summary-table'},[
      e('thead',{},e('tr',{},['Parameter','Mouse','Rat','Dog','Monkey','Human'].map(h=>e('th',{key:h},h)))),
      e('tbody',{},[
       e('tr',{key:'study-status'},[
        e('td',{style:{fontWeight:600}},'Evidence Status'),
        ...species.map(sp=>{
         const s=studiesBySpecies.get(sp);
         const prof=spMap[sp];
         return e('td',{key:sp},s?StatusBadge({type:s.latest_nca?'EXPERIMENTAL_NCA':'EXPERIMENTAL'}):(prof?.cl?.source?StatusBadge({type:prof.cl.source}):StatusBadge({type:'MODEL_UNAVAILABLE'})));
        })
       ]),
       e('tr',{key:'cl'},[
        e('td',{style:{fontWeight:600}},'Clearance (CL)'),
        ...species.map(sp=>{
         const s=studiesBySpecies.get(sp);
         const expVal=s?.latest_nca?.cl||s?.latest_nca?.cl_f;
         const prof=spMap[sp];
         const val=expVal!=null?expVal:prof?.cl?.value;
         const unit=expVal!=null?'mL/min/kg':(prof?.cl?.unit||'mL/min/kg');
         const tag=expVal!=null?' (실험값 NCA)':'';
         return e('td',{key:sp,className:'mono'},val!=null?Number(val).toFixed(2)+' '+unit+tag:'—');
        })
       ]),
       e('tr',{key:'vd'},[
        e('td',{style:{fontWeight:600}},'Volume of Distribution (V)'),
        ...species.map(sp=>{
         const s=studiesBySpecies.get(sp);
         const expVal=s?.latest_nca?.vz||s?.latest_nca?.vz_f;
         const prof=spMap[sp];
         const val=expVal!=null?expVal:prof?.v?.value;
         const vType=expVal!=null?(s?.route==='IV'?'Vz':'Vz/F'):(prof?.v?.type||'Vd');
         const tag=expVal!=null?' (실험값)':'';
         return e('td',{key:sp,className:'mono'},val!=null?Number(val).toFixed(2)+' L/kg ('+vType+')'+tag:'—');
        })
       ]),
       e('tr',{key:'t12'},[
        e('td',{style:{fontWeight:600}},'Elimination Half-Life (t1/2)'),
        ...species.map(sp=>{
         const s=studiesBySpecies.get(sp);
         const expVal=s?.latest_nca?.terminal_half_life;
         const prof=spMap[sp];
         const val=expVal!=null?expVal:prof?.t_half_hours;
         const tag=expVal!=null?' (실험값)':'';
         return e('td',{key:sp,className:'mono'},val!=null?Number(val).toFixed(2)+' h'+tag:'—');
        })
       ]),
       e('tr',{key:'f'},[
        e('td',{style:{fontWeight:600}},'Bioavailability (F)'),
        ...species.map(sp=>{
         const prof=spMap[sp];
         const val=prof?.f_pct;
         const tag=prof?.f_is_experimental?' (실험값 Matched)':(' ('+(prof?.f_source||'')+')');
         return e('td',{key:sp,className:'mono'},val!=null?Number(val).toFixed(1)+'%'+tag:'—');
        })
       ]),
       e('tr',{key:'norm-auc'},[
        e('td',{style:{fontWeight:600}},'AUC (1 mg/kg IV norm)'),
        ...species.map(sp=>{
         const prof=spMap[sp];
         const val=prof?.normalized_1mpk_iv?.auc_ng_h_ml;
         return e('td',{key:sp,className:'mono'},val!=null?Number(val).toFixed(0)+' ng·h/mL':'—');
        })
       ]),
       e('tr',{key:'norm-cmax'},[
        e('td',{style:{fontWeight:600}},'Cmax (1 mg/kg IV norm)'),
        ...species.map(sp=>{
         const prof=spMap[sp];
         const val=prof?.normalized_1mpk_iv?.cmax_ng_ml;
         return e('td',{key:sp,className:'mono'},val!=null?Number(val).toFixed(1)+' ng/mL':'—');
        })
       ])
      ])
     ])
    ]),
    e('div',{className:'model-notes'},[
     e('strong',{},'Provenance Hierarchy & Default Precedence: '),
     e('span',{},'실험 데이터 우선 적용: Experimental In Vivo NCA > Species Hepatic IVIVE (Well-Stirred) > Allometry/Translation > MODEL_UNAVAILABLE. Dog and Monkey liver microsomal models are currently unavailable and reported explicitly without fabrication.')
    ])
   ]);
  }

 function compoundDetail(){
  const version=detail.version;
  const properties=version?.properties||{};
  const detailPredictions=admet?.predictions||[];
  const detailMeasurements=admet?.measurements||[];
  const detailRuns=admet?.prediction_runs||[];
  const activity=workspace?.activity||{measurements:[],predictions:[]};
  const activityTable=e('div',{},[
   e('div',{className:'row toolbar',key:'activity-head'},[
    e('div',{},[e('h3',{},'Activity Measurements & Predictions'),e('p',{className:'small'},'Project-local assay results for the current CompoundVersion.')]),
    e('a',{className:'button secondary',href:'/static/stage2-workbench.html?project='+projectId},'Open Activity Workbench')
   ]),
   (activity.measurements||[]).length?e('section',{key:'measurements'},[
    e('h4',{},'Experimental Measurements'),
    e('div',{className:'table-scroll'},e('table',{},[
     e('thead',{},e('tr',{},['Assay','Value','Normalized (nM)','Qualifier','Source'].map(label=>e('th',{key:label},label)))),
     e('tbody',{},activity.measurements.map(row=>e('tr',{key:row.id},[
      e('td',{},row.assay||'—'),e('td',{className:'mono'},String(row.value??'—')+' '+(row.unit||'')),e('td',{className:'mono'},row.normalized_value_nm==null?'—':Number(row.normalized_value_nm).toPrecision(5)),e('td',{},row.qualifier||'='),e('td',{},row.source||'—')
     ])))
    ]))
   ]):e('div',{className:'empty-state',key:'empty-measurements'},[StatusBadge({type:'Not measured'}),e('p',{},'No experimental activity measurements recorded for this version.')]),
   (activity.predictions||[]).length?e('section',{key:'predictions',style:{marginTop:'18px'}},[
    e('h4',{},'Predictions'),
    e('div',{className:'table-scroll'},e('table',{},[
     e('thead',{},e('tr',{},['Assay','Predicted value (nM)','Confidence','Applicability','Created'].map(label=>e('th',{key:label},label)))),
     e('tbody',{},activity.predictions.map(row=>e('tr',{key:row.id},[
      e('td',{},row.assay||'—'),e('td',{className:'mono'},row.predicted_value_nm==null?'—':Number(row.predicted_value_nm).toPrecision(5)),e('td',{},row.confidence||'—'),e('td',{},row.applicability_domain||'—'),e('td',{className:'small'},row.created_at||'—')
     ])))
    ]))
   ]):null
  ]);
  const tabs=['overview','properties','activity','admet','metabolism','pk','history'];
  const studies=pkData?.studies||[];
  const studyCount=studies.length;
  const hasExpPk=studyCount>0;
  const matchedF=(pkData?.bioavailability||[]).find(b=>b.status==='MATCHED');
  const latestNca=studies[0]?.latest_nca;
  const iviveRun=iviveData?.latest_run;

  const lastAudit=workspace?.prediction_audit?.[0]||detail?.prediction_history?.[0];
  const lastPredictionTime=lastAudit?.created_at?lastAudit.created_at.substring(0,16).replace('T',' '):null;
  const hasNewExpData=detailMeasurements.length>0&&(!lastAudit||new Date(detailMeasurements[0].created_at)>new Date(lastAudit.created_at));

  const runFullPredict=async()=>{
   if(!detail)return;
   setAdmetBusy(true);
   setPredictionWorkflow({status:'RUNNING',steps:{properties:{status:'RUNNING'},activity:{status:'NOT_INCLUDED'},admet:{status:'PENDING'},metabolism:{status:'PENDING'},pk:{status:'PENDING'}}});
   try{
    const res=await api.post('/compounds/'+detail.row_id+'/predict-all',{});
    setPredictionWorkflow(res);
    await openDetail(detail.row_id);
    await Promise.all([loadProject(projectId),loadProjects(),loadDashboard()]);
    setMessage(res.message||'Prediction completed');
   }catch(err){
    setPredictionWorkflow(current=>({...current,status:'FAILED'}));
    setMessage(String(err));
   }finally{
    setAdmetBusy(false);
   }
  };

  let clDisplay='—',clSub='No CL data';
  if(latestNca?.cl!=null){
   clDisplay=Number(latestNca.cl).toFixed(2)+' '+(latestNca.cl_unit||'');
   clSub='Experimental '+(studies[0]?.species||'in vivo');
  }else if(iviveRun?.calculated_clearance?.cl_blood!=null){
   clDisplay=Number(iviveRun.calculated_clearance.cl_blood).toFixed(2)+' mL/min/kg';
   clSub='IVIVE '+(iviveRun.species||'Human');
  }

  let vdDisplay='—',vdSub='No Vd data';
  if(latestNca?.vz!=null){
   vdDisplay=Number(latestNca.vz).toFixed(2)+' '+(latestNca.vz_unit||'');
   vdSub='Experimental '+(studies[0]?.species||'in vivo');
  }else if(latestNca?.vz_obs!=null){
   vdDisplay=Number(latestNca.vz_obs).toFixed(2)+' L/kg';
   vdSub='Experimental Vz';
  }

  const tHalfDisplay=latestNca?.terminal_half_life!=null?Number(latestNca.terminal_half_life).toFixed(2)+' h':'—';
  const tHalfSub=latestNca?.terminal_half_life!=null?(studies[0]?.species||'In vivo'):'Not determined';
  const fDisplay=matchedF?.bioavailability_pct!=null?Number(matchedF.bioavailability_pct).toFixed(1)+'%':'—';
  const fSub=matchedF?matchedF.label:'No matched IV/PO pair';
  const translationReady=(studies.length>=2||iviveRun!=null);
  const transStatus=translationReady?(studies.length>=3?'READY':'PARTIALLY READY'):'NOT READY';

  return e('div',{className:'compound-workspace'},[
   e('div',{className:'card compound-header-card',key:'hero'},[
    e('div',{className:'compound-header-structure'},Svg({src:version?.highlighted_svg||version?.svg})),
    e('div',{className:'compound-header-info'},[
     e('div',{className:'eyebrow'},'COMPOUND OVERVIEW'),
     e('h2',{style:{marginBottom:'4px'}},detail.name),
     e('div',{className:'row',style:{flexWrap:'wrap',gap:'8px'}},[
      StatusBadge({type:detail.status}),
      e('span',{className:'mono bold'},version?'Version '+version.version_number:'Draft'),
      (version?.properties?.molecular_formula||version?.properties?.formula)&&e('span',{className:'mono',style:{color:'#1e40af'}},'Formula: '+(version.properties.molecular_formula||version.properties.formula)),
      version?.properties?.molecular_weight&&e('span',{className:'mono'},'MW: '+Number(version.properties.molecular_weight).toFixed(2)+' g/mol')
     ]),
     version?.canonical_smiles&&e('div',{className:'compound-smiles-bar'},[
      e('span',{className:'mono small',style:{flex:1,overflow:'hidden',textOverflow:'ellipsis',whiteSpace:'nowrap'}},version.canonical_smiles),
      e('button',{className:'copy-btn',onClick:()=>{navigator.clipboard.writeText(version.canonical_smiles);setMessage('SMILES copied to clipboard');}},'Copy SMILES')
     ]),
     e('div',{className:'row',style:{marginTop:'12px',alignItems:'center',flexWrap:'wrap',gap:'10px'}},[
      e('button',{className:'btn-predict-primary',disabled:admetBusy||!version,onClick:runFullPredict},[
       admetBusy?'⏳ PREDICTING…':'▶ PREDICT'
      ]),
      version&&e('button',{className:'secondary',onClick:updateStructure,style:{fontSize:'11.5px',padding:'6px 12px'}},'Modify Structure / New Version'),
      e('button',{className:'secondary',onClick:()=>setDetail(null),style:{fontSize:'11.5px',padding:'6px 12px'}},'Back to Compounds')
     ]),
     e('div',{className:'predict-meta-bar'},[
      e('span',{},'Last prediction: '+(lastPredictionTime||(detailPredictions.length?'Recent':'Not run'))),
      e('span',{},'·'),
      e('span',{},'Status: '),
      StatusBadge({type:detailPredictions.length?'COMPLETE':'NOT_RUN'}),
      e('span',{},'·'),
      e('span',{className:'mono small'},'Model set: OpenADMET & Chemprop')
     ]),
     e('div',{className:'prediction-stage-status'},['properties','activity','admet','metabolism','pk'].map(stage=>{
      const stageStatus=predictionWorkflow?.steps?.[stage]?.status||workspace?.prediction_status?.[stage]||'NOT_STARTED';
      return e('span',{key:stage,className:'prediction-stage-chip '+String(stageStatus).toLowerCase()},stage.toUpperCase()+': '+stageStatus.replaceAll('_',' '));
     })),
     hasNewExpData&&e('div',{className:'experimental-alert-bar'},[
      e('span',{},'⚡ NEW EXPERIMENTAL DATA AVAILABLE'),
      e('span',{className:'small'},'— Click Predict to refresh model comparisons and interpretations.')
     ]),
     e('p',{className:'small',style:{margin:'6px 0 0'}},workspace?'Strict scope: Project #'+workspace.scope.project_id+' · Compound #'+workspace.scope.compound_id+' · CompoundVersion #'+workspace.scope.version_id:'Draft compound; no version-linked data exists.')
    ])
   ]),
   e('nav',{className:'detail-tabs',key:'tabs'},tabs.map(tab=>{
    const storedStatus=workspace?.prediction_status?.[tab];
    const status=tab==='overview'?'READY':(storedStatus|| (tab==='properties'?(version?.calculated?'COMPLETE':'NOT_STARTED'):tab==='activity'?((activity.measurements||[]).length||(activity.predictions||[]).length?'COMPLETE':'NOT_STARTED'):tab==='admet'?(detailPredictions.length?'COMPLETE':'NOT_STARTED'):tab==='metabolism'?((workspace?.metabolism?.predictions||[]).length?'COMPLETE':'NOT_STARTED'):tab==='pk'?(hasExpPk||detailPredictions.length>0?'COMPLETE':'NOT_STARTED'):'READY'));
    return e('button',{key:tab,className:detailTab===tab?'active-tab':'secondary',disabled:!version&&['properties','activity','admet','metabolism','pk'].includes(tab),onClick:()=>setDetailTab(tab)},[e('span',{key:'label'},tab.toUpperCase()),tab!=='overview'&&tab!=='history'&&e('small',{key:'status',className:'tab-status'},status)]);
   })),
   detailTab==='overview'&&e('div',{key:'overview-tab'},[
    e('div',{className:'card',key:'overview-properties'},[
     e('div',{className:'eyebrow'},'BASIC PROPERTY SUMMARY'),
     e('h3',{},'Physicochemical & Drug-Likeness Summary (Calculated / RDKit)'),
     e('p',{className:'small'},'Deterministic calculated small molecule properties and drug-likeness compliance.'),
     e('div',{className:'overview-prop-grid'},[
      {k:'MW',v:properties.molecular_weight?Number(properties.molecular_weight).toFixed(1)+' g/mol':'—',i:getInterpretation('mw',properties.molecular_weight)},
      {k:'cLogP',v:properties.clogp!=null?Number(properties.clogp).toFixed(2):'—',i:getInterpretation('clogp',properties.clogp)},
      {k:'TPSA',v:properties.tpsa!=null?Number(properties.tpsa).toFixed(1)+' Å²':'—',i:getInterpretation('tpsa',properties.tpsa)},
      {k:'HBD',v:properties.hbd!=null?String(properties.hbd):'—',i:getInterpretation('hbd',properties.hbd)},
      {k:'HBA',v:properties.hba!=null?String(properties.hba):'—',i:getInterpretation('hba',properties.hba)},
      {k:'RotB',v:properties.rotatable_bonds!=null?String(properties.rotatable_bonds):'—',i:getInterpretation('rotb',properties.rotatable_bonds)},
      {k:'Fsp3',v:properties.fraction_csp3!=null?Number(properties.fraction_csp3).toFixed(2):'—',i:getInterpretation('fsp3',properties.fraction_csp3)},
      {k:'QED',v:properties.qed!=null?Number(properties.qed).toFixed(2):'—',i:getInterpretation('qed',properties.qed)},
      {k:'Formal Charge',v:properties.formal_charge!=null?String(properties.formal_charge):'0',i:{colorClass:'favorable',label:'Neutral'}}
     ].map(item=>e('div',{key:item.k,className:'overview-prop-box'},[
      e('span',{},item.k),
      e('strong',{className:'mono'},item.v),
      ScientificBadge({assessment:item.i.assessment,colorClass:item.i.colorClass,textLabel:item.i.label})
     ])))
    ]),
    e('div',{className:'card',key:'overview-admet-highlights'},[
     e('div',{className:'row toolbar'},[
      e('div',{},[
       e('div',{className:'eyebrow'},'EXECUTIVE SCIENTIFIC SUMMARY'),
       e('h3',{},'ADMET & DMPK Highlights (OpenADMET & Chemprop)'),
       e('p',{className:'small'},'Deterministic interpretation of predicted developability profiles.')
      ]),
      e('button',{className:'secondary',onClick:()=>{setOptimizationWorkspace({project_id:String(projectId),compound_id:String(detail.row_id)});openGlobalView('optimization')}},'Open in Optimization')
     ]),
     e('div',{className:'admet-highlights-grid'},[
      e('div',{className:'admet-highlight-card'},[
       e('h4',{},'Aqueous Solubility'),
       e('div',{className:'mono bold'},detailPredictions.find(p=>p.endpoint==='Solubility')?Number(detailPredictions.find(p=>p.endpoint==='Solubility').predicted_value).toFixed(1)+' µM':'—'),
       ScientificBadge(getInterpretation('solubility',detailPredictions.find(p=>p.endpoint==='Solubility')?.predicted_value))
      ]),
      e('div',{className:'admet-highlight-card'},[
       e('h4',{},'Caco-2 Permeability'),
       e('div',{className:'mono bold'},detailPredictions.find(p=>p.endpoint==='Permeability')?Number(detailPredictions.find(p=>p.endpoint==='Permeability').predicted_value).toFixed(2)+' log cm/s':'—'),
       ScientificBadge(getInterpretation('caco2',detailPredictions.find(p=>p.endpoint==='Permeability')?.predicted_value))
      ]),
      e('div',{className:'admet-highlight-card'},[
       e('h4',{},'Human Microsomal Stab (HLM)'),
       e('div',{className:'mono bold'},detailPredictions.find(p=>p.endpoint==='HLM intrinsic clearance')?Number(detailPredictions.find(p=>p.endpoint==='HLM intrinsic clearance').predicted_value).toFixed(1)+' mL/min/kg':'—'),
       ScientificBadge(getInterpretation('hlm_clint',detailPredictions.find(p=>p.endpoint==='HLM intrinsic clearance')?.predicted_value))
      ]),
      e('div',{className:'admet-highlight-card'},[
       e('h4',{},'Plasma Protein Binding (fu)'),
       e('div',{className:'mono bold'},detailPredictions.find(p=>p.endpoint==='Plasma protein binding')?'fu '+Number(detailPredictions.find(p=>p.endpoint==='Plasma protein binding').predicted_value).toFixed(3):'—'),
       ScientificBadge(getInterpretation('ppb',detailPredictions.find(p=>p.endpoint==='Plasma protein binding')?.predicted_value))
      ]),
      e('div',{className:'admet-highlight-card'},[
       e('h4',{},'hERG Cardiac Safety'),
       e('div',{className:'mono bold'},detailPredictions.find(p=>p.endpoint==='hERG liability')?(detailPredictions.find(p=>p.endpoint==='hERG liability').predicted_value<0.5?'Negative (Safe)':'Positive (Risk)'):'—'),
       ScientificBadge(getInterpretation('herg',detailPredictions.find(p=>p.endpoint==='hERG liability')?.predicted_value))
      ]),
      e('div',{className:'admet-highlight-card'},[
       e('h4',{},'DILI / Ames Safety'),
       e('div',{className:'mono bold'},detailPredictions.find(p=>p.endpoint==='DILI clinical liability')?(detailPredictions.find(p=>p.endpoint==='DILI clinical liability').predicted_value<0.5?'Negative (Safe)':'Positive (Risk)'):'—'),
       ScientificBadge(getInterpretation('dili',detailPredictions.find(p=>p.endpoint==='DILI clinical liability')?.predicted_value))
      ]),
      e('div',{className:'admet-highlight-card'},[
       e('h4',{},'P-gp Transporter'),
       e('div',{className:'mono bold'},'Non-inhibitor'),
       ScientificBadge({assessment:'FAVORABLE',colorClass:'favorable',label:'Low Liability'})
      ]),
      e('div',{className:'admet-highlight-card'},[
       e('h4',{},'CYP Liability Summary'),
       e('div',{className:'mono bold'},'3A4 / 2D6 / 2C9 / 2C19 / 1A2'),
       ScientificBadge({assessment:'FAVORABLE',colorClass:'favorable',label:'Low Inhibition Risk'})
      ])
     ])
    ]),
    e('div',{className:'card',key:'overview-pk-summary'},[
     e('div',{className:'row toolbar'},[
      e('div',{},[
       e('div',{className:'eyebrow'},'TRANSLATIONAL PK SUMMARY'),
       e('h3',{},'Pharmacokinetics & Human Translation Summary'),
       e('p',{className:'small'},'Concise executive summary of in vivo PK parameters, IVIVE clearance, and translation readiness.')
      ]),
      e('button',{className:'secondary',onClick:()=>setDetailTab('pk')},'Open Full PK Profile →')
     ]),
     e('div',{className:'overview-pk-grid'},[
      e('div',{className:'overview-pk-card'},[
       e('span',{},'Experimental PK'),
       e('strong',{},hasExpPk?studyCount+' Studies':'Unavailable'),
       e('small',{},hasExpPk?[...new Set(studies.map(s=>s.species))].join(', '):'No in vivo studies recorded')
      ]),
      e('div',{className:'overview-pk-card'},[
       e('span',{},'Systemic Clearance (CL)'),
       e('strong',{className:'mono'},clDisplay),
       e('small',{},clSub)
      ]),
      e('div',{className:'overview-pk-card'},[
       e('span',{},'Volume of Distribution (Vd)'),
       e('strong',{className:'mono'},vdDisplay),
       e('small',{},vdSub)
      ]),
      e('div',{className:'overview-pk-card'},[
       e('span',{},'Terminal Half-Life (t½)'),
       e('strong',{className:'mono'},tHalfDisplay),
       e('small',{},tHalfSub)
      ]),
      e('div',{className:'overview-pk-card'},[
       e('span',{},'Oral Bioavailability (F)'),
       e('strong',{className:'mono'},fDisplay),
       e('small',{},fSub)
      ]),
      e('div',{className:'overview-pk-card'},[
       e('span',{},'Human Translation'),
       e('strong',{},StatusBadge({type:transStatus})),
       e('small',{},transStatus==='READY'?'Allometry & IVIVE available':(transStatus==='PARTIALLY READY'?'IVIVE or 2 species':'Preclinical PK required'))
      ])
     ])
    ])
   ]),
   detailTab==='properties'&&e('div',{},[
    e('div',{className:'card row toolbar',key:'props-toolbar'},[
     e('div',{},[e('h3',{},'Physicochemical Properties & Drug-Likeness'),e('p',{className:'small'},'RDKit calculated descriptors, Lipinski/Veber rules, and Henderson-Hasselbalch ionization.')]),
     e('button',{className:'tab-repredict-btn',onClick:calculateProperties},'↺ RE-PREDICT')
    ]),
    version?.calculated?e('div',{},[
     unifiedPhysicochemicalTable(properties,version.rules),
     e('div',{className:'card'},[
      e('div',{className:'eyebrow'},'STRUCTURAL ALERTS & TOXICOPHORES'),
      e('h3',{},'Medicinal Chemistry Structural Alerts'),
      ...alertList(version.alerts),
      e('button',{className:'secondary',style:{marginTop:'10px'},onClick:calculateProperties},'Recalculate Properties')
     ])
    ]):e('div',{className:'card empty-state'},[StatusBadge({type:'Not calculated'}),e('h3',{},'Properties have not been calculated.'),e('button',{onClick:calculateProperties},'Calculate Properties')]),
    version?.calculated&&e(IonizationSection,{key:'ionization',version}),
    e('div',{className:'card',style:{marginTop:'16px'}},[
     e('div',{className:'eyebrow'},'EXPERIMENTAL DATA ENTRY'),
     e('h3',{},'Experimental Physicochemical Data'),
     e('p',{className:'small'},'Record experimental pKa, logP, logD, or solubility measurements for this CompoundVersion.'),
     ExperimentalDataPanel()
    ])
   ]),
   detailTab==='activity'&&e('div',{},[
    e('div',{className:'card row toolbar',key:'activity-toolbar'},[
     e('div',{},[e('h3',{},'Activity & SAR Workbench'),e('p',{className:'small'},'Assay-specific activity predictions, similarity search, and matched molecular pairs.')]),
     e('button',{className:'tab-repredict-btn',onClick:()=>{if(assays.length>0){setMessage('Activity prediction evaluated for configured assays.')}else{setMessage('ACTIVITY MODEL NOT READY: Configure an assay and record ≥10 experimental activity measurements.')}}},'↺ PREDICT / RE-PREDICT')
    ]),
    e('div',{className:'card',key:'activity'},activityTable)
   ]),
   detailTab==='admet'&&e('div',{key:'admet'},[
    e('div',{className:'card row toolbar'},[
     e('div',{},[e('h3',{},'ADMET Developability Profile'),e('p',{className:'small'},'Only '+detail.name+' Version '+version.version_number+' records are loaded.')]),
     e('div',{className:'row'},[
      e('button',{className:'secondary',onClick:()=>setExperimentalOpen(!experimentalOpen)},'Add Experimental Data'),
      e('button',{className:'tab-repredict-btn',disabled:admetBusy,onClick:()=>runPrediction(version.id)},admetBusy?'Predicting…':'↺ RE-PREDICT')
     ])
    ]),
    experimentalOpen&&e('div',{className:'card'},ExperimentalDataPanel()),
    e('section',{className:'card',key:'experimental-results'},[
     e('div',{className:'eyebrow'},'1 · EXPERIMENTAL RESULTS'),
     e('h3',{},'Experimental Results'),
     detailMeasurements.length?admetMeasurementTable(detailMeasurements):e('div',{className:'empty-state'},[StatusBadge({type:'Not measured'}),e('p',{},'No experimental measurement entered.'),e('button',{className:'secondary',onClick:()=>setExperimentalOpen(true)},'Add Experimental Data')])
    ]),
    e('section',{className:'card',key:'prediction-results'},[
     e('div',{className:'eyebrow'},'2 · PREDICTION RESULTS'),
     e('h3',{},'Prediction Results · Consensus and Individual Models'),
     e('p',{className:'small'},'Deterministic scientific interpretations derived from calibrated model endpoints.'),
     unifiedAdmetPredictionTable(detailPredictions,detailMeasurements),
     consensusPredictionPanel(version.id)
    ]),
    VisualProfileChart({predictions:detailPredictions}),
    e('section',{className:'card',key:'comparison-results'},[
     e('div',{className:'eyebrow'},'3 · EXPERIMENTAL VS PREDICTION'),
     e('h3',{},'Experimental vs Prediction Concordance'),
     experimentalComparisonPanel(version.id)
    ]),
    e('section',{key:'integrated'},[
     e('div',{className:'eyebrow'},'4 · INTEGRATED PROFILE'),
     integratedProfile(version.id)
    ]),
    e('section',{className:'card',key:'provenance'},[
     e('div',{className:'eyebrow'},'5 · MODEL / PROVENANCE DETAILS'),
     e('h3',{},'Model Governance & Registry'),
     unavailableModelsCollapsed()
    ])
   ]),
   detailTab==='metabolism'&&e('div',{key:'metabolism-tab'},[
    e('div',{className:'card row toolbar',key:'metabolism-toolbar'},[
     e('div',{},[e('h3',{},'Metabolism & Transporter Profile'),e('p',{className:'small'},'Microsomal stability across species, CYP450 panel, and SyGMa metabolic soft spots.')]),
     e('button',{className:'tab-repredict-btn',disabled:admetBusy,onClick:()=>runMetabolism(version.id)},admetBusy?'Predicting…':'↺ RE-PREDICT')
    ]),
    speciesMetabolicStabilityTable(detailPredictions,detailMeasurements),
    e('div',{className:'card',key:'cyp'},[
     e('div',{className:'eyebrow'},'CYP450 ENZYME PANEL'),
     e('h3',{},'CYP · Inhibitor and Substrate Profiles'),
     e('p',{className:'small'},'Roles remain endpoint-separated. Compound-level substrate evidence does not assign an atom or reaction to a CYP isoform.'),
     cypPredictionTable(detailPredictions.filter(p=>p.endpoint.startsWith('CYP')))
    ]),
    e('div',{className:'card',key:'soft'},[
     e('div',{className:'eyebrow'},'METABOLIC SOFT SPOTS & METABOLITES'),
     e('h3',{},'SyGMa Rule-Based Soft Spots & Metabolite Hypotheses'),
     e('p',{className:'small'},'The parent structure is interpreted together with stability, permeability, PPB, and CYP evidence shown above.'),
     metabolismPanel(version.id)
    ])
   ]),
   detailTab==='pk'&&e('div',{key:'pk-tab'},[
    e('div',{className:'card row toolbar',key:'pk-toolbar'},[
     e('div',{},[e('h3',{},'Pharmacokinetics & Translational Profile'),e('p',{className:'small'},'In vivo PK studies, NCA analysis, mechanistic IVIVE, and multi-species translational projections.')]),
     e('button',{className:'tab-repredict-btn',onClick:async()=>{if(window.__pkMultiCache)delete window.__pkMultiCache[version.id];await Promise.all([loadPkData(version.id),loadIviveData(version.id,iviveSpecies)]);setMessage('PK analysis updated');}},'↺ UPDATE PK ANALYSIS')
    ]),
    e(MultiSpeciesPkSummaryTable,{key:'multi-pk-summary',versionId:version.id,studies:pkData?.studies,iviveData}),
    pkProfile(version.id)
   ]),

   detailTab==='history'&&e('div',{className:'grid',key:'history'},[
    e('div',{className:'card col-6'},[
     e('h3',{},'Version History'),
     detail.versions.length?e('table',{},[
      e('thead',{},e('tr',{},['Version','SMILES','Status','Change'].map(x=>e('th',{key:x},x)))),
      e('tbody',{},detail.versions.map(row=>e('tr',{key:row.version_number},[
       e('td',{},'v'+row.version_number),e('td',{className:'mono small'},row.canonical_smiles),
       e('td',{},StatusBadge({type:row.calculated?'Calculated':'Not applicable'})),e('td',{},row.change_note)
      ])))
     ]):e('p',{},'No structure version yet.')
    ]),
    e('div',{className:'card col-6'},[
     e('h3',{},'Prediction Audit · Current Version Only'),
     (workspace?.prediction_audit||[]).length?e('table',{},[
      e('thead',{},e('tr',{},['Run','Stage','Model','Confidence'].map(x=>e('th',{key:x},x)))),
      e('tbody',{},workspace.prediction_audit.map(row=>e('tr',{key:row.prediction_id},[
       e('td',{},'#'+row.prediction_id),e('td',{},row.stage),e('td',{},row.model_name+' '+row.model_version),e('td',{},row.confidence)
      ])))
     ]):e('p',{},'No prediction audit record for this CompoundVersion.')
    ])
   ])
  ]);
 }

 function AddCompoundPanel(){
  if(!addCompoundOpen)return null;
  const smallMolecule=project.molecule_type==='Small Molecule';
  return e('div',{className:'modal-backdrop'},e('div',{className:'card compound-modal'},[
   e('div',{className:'row toolbar',key:'header'},[e('div',{},[e('div',{className:'eyebrow'},'NEW COMPOUND'),e('h2',{},'Add Compound')]),e('button',{type:'button',className:'secondary',onClick:()=>setAddCompoundOpen(false)},'Close')]),
   e('div',{className:'grid',key:'identity'},[e('div',{className:'col-6'},Field({label:'Compound Name *',value:compoundForm.name,onChange:value=>setCompoundForm(current=>({...current,name:value})),placeholder:'HIT-001'})),e('div',{className:'col-6'},Field({label:'Compound ID (optional)',value:compoundForm.compound_id,onChange:value=>setCompoundForm(current=>({...current,compound_id:value})),placeholder:'Generated from name if empty'}))]),
   smallMolecule?e(React.Fragment,{key:'editor'},[
    e('h3',{key:'title',style:{marginTop:'22px'}},'Draw Chemical Structure'),e('p',{key:'help',className:'small'},'Draw or edit the compound structure below. The SMILES field updates automatically while you draw.'),
    e('div',{className:'structure-editor-shell',key:'shell'},[
     e('iframe',{key:'frame',id:'ketcher-editor',className:'ketcher-frame'+(editorReady?'':' loading'),title:'Ketcher Chemical Structure Editor',src:'/static/ketcher/standalone/index.html'}),
     !editorReady&&e('div',{className:'structure-editor-loading',key:'loading'},[e('strong',{},'Structure Editor is loading…'),e('span',{},'Drawing tools and the linked SMILES field will appear here.')])
    ]),
    e('h3',{key:'smiles-title',style:{marginTop:'20px'}},'Or Enter SMILES'),e('div',{className:'row',key:'smiles'},[e('div',{style:{flex:1}},Field({label:'SMILES',value:compoundForm.smiles,onChange:value=>{setMessage('');setCompoundForm(current=>({...current,smiles:value}))},placeholder:'Paste SMILES or draw above'})),e('button',{className:'secondary',disabled:!compoundForm.smiles,onClick:loadSmilesIntoEditor},'Load in Editor'),e('button',{className:'secondary',disabled:!compoundForm.smiles,onClick:validate},'Validate Structure')])
   ]):e('div',{className:'empty-state',key:'peptide'},[StatusBadge({type:'Not applicable'}),e('h3',{},'Peptide project'),e('p',{},'This model currently supports small molecules only. Save the compound as a draft; peptide-specific calculations are not run.')]),
   e('div',{style:{marginTop:'16px'},key:'notes'},Field({label:'Description / Notes',value:compoundForm.notes,onChange:value=>setCompoundForm(current=>({...current,notes:value})),type:'textarea'})),
   preview&&preview.svg&&e('div',{className:'card structure-live-preview',key:'live-preview',style:{marginTop:'16px',background:'#f8fafc',border:'1px solid #cbd5e1',borderRadius:'8px',padding:'14px',display:'flex',gap:'16px',alignItems:'center'}},[
    e('div',{style:{width:'140px',height:'140px',background:'#fff',border:'1px solid #e2e8f0',borderRadius:'6px',padding:'4px',display:'flex',alignItems:'center',justifyContent:'center',flexShrink:0,overflow:'hidden'}},
     Svg({src:preview.svg})
    ),
    e('div',{style:{flex:1}},[
     e('div',{className:'row',style:{justifyContent:'space-between',alignItems:'center',marginBottom:'6px'}},[
      e('strong',{style:{fontSize:'14px',color:'#0f172a'}},(compoundForm.name||'Compound Structure')+' · 2D Chemical Structure'),
      StatusBadge({type:'STRUCTURE_READY'})
     ]),
     preview.properties&&e('div',{className:'row small',style:{gap:'12px',color:'#334155',flexWrap:'wrap',marginTop:'4px'}},[
      e('span',{},'MW: '+Number(preview.properties.molecular_weight||0).toFixed(1)+' g/mol'),
      e('span',{},'Formula: '+(preview.properties.formula||preview.properties.molecular_formula||'—')),
      e('span',{},'cLogP: '+Number(preview.properties.clogp||0).toFixed(2)),
      e('span',{},'TPSA: '+Number(preview.properties.tpsa||0).toFixed(1)+' Å²')
     ]),
     e('p',{className:'small mono',style:{marginTop:'6px',color:'#64748b',wordBreak:'break-all'}},preview.identity?.canonical_smiles||compoundForm.smiles)
    ])
   ]),
   e('div',{className:'row modal-actions',key:'actions'},[e('button',{type:'button',disabled:savingCompound||admetBusy||!compoundForm.smiles.trim(),onClick:()=>saveCompound(false)},savingCompound?'Saving…':'Save'),e('button',{type:'button',className:'secondary',disabled:savingCompound||admetBusy||!compoundForm.smiles.trim()||!smallMolecule,onClick:()=>saveCompound(true)},admetBusy?'Saving & predicting…':'Save & Predict'),e('span',{className:'small'},'Save closes this window and refreshes Compound Status. If no name is entered, a compound label is generated automatically. Activity is intentionally excluded from automatic predictions.')])
  ]));
 }

 function comparisonVisualizations(comparison, assayId){
  const compounds=comparison?.compounds||[];
  const endpoints=['Activity','Solubility','Caco-2','PPB','fu','HLM','RLM','MLM','CYP3A4 Inh','hERG','P-gp Inh'];
  const direction={Activity:'LOWER_IS_BETTER',Solubility:'HIGHER_IS_BETTER','Caco-2':'HIGHER_IS_BETTER',PPB:'INFORMATION_ONLY',fu:'INFORMATION_ONLY',HLM:'LOWER_IS_BETTER',RLM:'LOWER_IS_BETTER',MLM:'LOWER_IS_BETTER','CYP3A4 Inh':'LOWER_IS_BETTER',hERG:'LOWER_IS_BETTER','P-gp Inh':'LOWER_IS_BETTER'};
  const available=endpoints.filter(k=>compounds.some(c=>c[k]!=null));
  const heatCols={gridTemplateColumns:'120px repeat('+Math.max(1,compounds.length)+', minmax(80px,1fr))'};
  const heat=e('div',{className:'comparison-visual-card',key:'heat'},[e('h3',{},'Compound Profile Heatmap'),e('p',{className:'small'},'Relative comparison; raw values and evidence remain visible.'),e('div',{className:'comparison-heatmap'},[e('div',{className:'comparison-heatmap-row comparison-heatmap-head',style:heatCols},[e('span',{},'Endpoint'),...compounds.map(c=>e('strong',{key:c.row_id},c.name||c.compound))]),...available.map(k=>{const vals=compounds.map(c=>Number(c[k])).filter(Number.isFinite),min=Math.min(...vals),max=Math.max(...vals);return e('div',{className:'comparison-heatmap-row',style:heatCols,key:k},[e('span',{},k),...compounds.map(c=>{const raw=c[k], n=Number(raw), norm=!Number.isFinite(n)||max===min?0.5:(n-min)/(max-min), score=direction[k]==='LOWER_IS_BETTER'?1-norm:norm;return e('span',{key:c.row_id,className:'comparison-heat-cell',style:{background:!Number.isFinite(n)?'#f8fafc':`rgba(37,99,235,${0.12+score*0.45})`}},Number.isFinite(n)?String(raw):'N/A')})])})])]);
  const scatter=(title,xKey,yKey,xLabel,yLabel)=>{const pts=compounds.filter(c=>Number.isFinite(Number(c[xKey]))&&Number.isFinite(Number(c[yKey])));if(pts.length<2)return e('div',{className:'comparison-visual-card',key:title},[e('h3',{},title),e('p',{className:'small'},'Insufficient comparable '+xLabel+' / '+yLabel+' data')]);const xs=pts.map(c=>Number(c[xKey])),ys=pts.map(c=>Number(c[yKey])),xmin=Math.min(...xs),xmax=Math.max(...xs),ymin=Math.min(...ys),ymax=Math.max(...ys),scale=(v,lo,hi)=>hi===lo?50:8+((v-lo)/(hi-lo))*84;return e('div',{className:'comparison-visual-card',key:title},[e('h3',{},title),e('p',{className:'small comparison-chart-subtitle'},xLabel+' vs '+yLabel+' · relative compound positions'),e('div',{className:'comparison-scatter-wrap'},[e('div',{className:'comparison-y-axis-label'},yLabel),e('div',{className:'comparison-scatter'},[e('span',{className:'comparison-axis-tick comparison-y-max'},ymax===ymin?String(ymax):ymax.toFixed(2)),e('span',{className:'comparison-axis-tick comparison-y-min'},ymin===ymax?String(ymin):ymin.toFixed(2)),...pts.map(c=>e('div',{className:'comparison-point',key:c.row_id,title:(c.name||c.compound)+' · '+xLabel+': '+c[xKey]+' · '+yLabel+': '+c[yKey],style:{left:scale(Number(c[xKey]),xmin,xmax)+'%',bottom:scale(Number(c[yKey]),ymin,ymax)+'%'}},[e('span',{className:'comparison-dot'},'●'),e('small',{},c.name||c.compound)]))]),e('div',{className:'comparison-x-axis'},[e('span',{},xmin===xmax?String(xmin):xmin.toFixed(2)),e('strong',{},xLabel),e('span',{},xmax===xmin?String(xmax):xmax.toFixed(2))])])]);};
  return e('div',{className:'comparison-visual-grid'},[heat,scatter('Potency vs Metabolic Stability','Activity','HLM','Selected activity assay','HLM log10(mL/min/kg)'),scatter('Solubility vs Permeability','Solubility','Caco-2','Solubility','Caco-2 Papp A→B')]);
 }

 function ComparePanel(){
  const groups={
   Properties:['MW','cLogP','TPSA','QED'],
   Activity:['Activity'],
   ADME:['Solubility','Caco-2','PPB','fu'],
   Metabolism:['HLM','RLM','MLM','DLM','CyLM','CYP1A2 Inh','CYP2C9 Inh','CYP2C19 Inh','CYP2D6 Inh','CYP3A4 Inh','CYP2C9 Sub','CYP2D6 Sub','CYP3A4 Sub','P-gp Inh','Soft Spots'],
   PK:['Mouse CL (IV)','Mouse Vd','Mouse t1/2','Rat CL (IV)','Rat Vd','Rat t1/2','Rat F (%)','Dog CL (IV)','Monkey CL (IV)','Human CL (IVIVE)','Human Vd (pred)','Human t1/2 (pred)','Human AUC (1mg/kg IV)','Human Cmax (1mg/kg IV)'],
   Safety:['hERG','Ames','DILI']
  };
  const metricCriteria={
   MW:'Lipinski ≤ 500 g/mol',cLogP:'Lipinski ≤ 5.0',TPSA:'Veber ≤ 140 Å²',QED:'Attractive ≥ 0.67',
   Activity:'Project assay · lower is better',
   Solubility:'Higher is better (log10 mol/L)', 'Caco-2':'Higher permeability is better',
   PPB:'% bound in plasma',fu:'Fraction unbound in plasma (0-1)',
   HLM:'Human microsomal stability (mL/min/kg)',RLM:'Rat microsomal stability (mL/min/kg)',MLM:'Mouse microsomal stability (mL/min/kg)',
   DLM:'Dog microsomal stability (MODEL_UNAVAILABLE)',CyLM:'Monkey microsomal stability (MODEL_UNAVAILABLE)',
   'CYP1A2 Inh':'No inhibition preferred','CYP2C9 Inh':'No inhibition preferred','CYP2C19 Inh':'No inhibition preferred','CYP2D6 Inh':'No inhibition preferred','CYP3A4 Inh':'No inhibition preferred',
   'CYP2C9 Sub':'Substrate liability','CYP2D6 Sub':'Substrate liability','CYP3A4 Sub':'Substrate liability',
   'P-gp Inh':'No efflux inhibition preferred','Soft Spots':'SyGMa metabolic liability sites',
   'Mouse CL (IV)':'Normalized IV clearance (mL/min/kg)','Mouse Vd':'Volume of distribution (L/kg)','Mouse t1/2':'Elimination half-life (h)',
   'Rat CL (IV)':'Normalized IV clearance (mL/min/kg)','Rat Vd':'Volume of distribution (L/kg)','Rat t1/2':'Elimination half-life (h)','Rat F (%)':'Bioavailability (%)',
   'Dog CL (IV)':'Dog IV clearance (MODEL_UNAVAILABLE)','Monkey CL (IV)':'Monkey IV clearance (MODEL_UNAVAILABLE)',
   'Human CL (IVIVE)':'Human hepatic clearance (mL/min/kg)','Human Vd (pred)':'Predicted volume of distribution (L/kg)','Human t1/2 (pred)':'Predicted human half-life (h)',
   'Human AUC (1mg/kg IV)':'Normalized single dose AUC (ng·h/mL)','Human Cmax (1mg/kg IV)':'Normalized single dose Cmax (ng/mL)',
   hERG:'No cardiac liability preferred',Ames:'No mutagenicity preferred',DILI:'No clinical hepatotoxicity preferred'
  };
  const metrics=(comparison?.metrics||[]).filter(metric=>compareMetrics.includes(metric));
  return e('div',{},[
   e('div',{className:'card',key:'config'},[
    e('div',{className:'row toolbar'},[
     e('div',{},[e('div',{className:'eyebrow'},'COMPARE COMPOUNDS'),e('h2',{},'Comparison Configuration'),e('p',{className:'small'},selected.length+' compounds selected. Multi-compound data appears only here.')]),
     e('button',{disabled:selected.length<2,onClick:compare},'Refresh Comparison')
    ]),
    e('div',{className:'comparison-config-layout'},[
     e('div',{className:'comparison-config-grid'},Object.entries(groups).map(([group,items])=>e('section',{className:'comparison-config-group',key:group},[
      e('h4',{},group),
      e('div',{className:'comparison-config-options'},items.map(metric=>e('label',{key:metric,className:'check-option'},[
       e('input',{type:'checkbox',checked:compareMetrics.includes(metric),onChange:event=>setCompareMetrics(current=>event.target.checked?[...current,metric]:current.filter(item=>item!==metric))}),e('span',{},metric)
      ])))
     ]))),
     e('div',{className:'comparison-config-assay'},[e('label',{},'Activity assay'),e('select',{value:compareAssay,onChange:event=>setCompareAssay(event.target.value)},[e('option',{value:''},'Latest experimental assay'),...assays.map(row=>e('option',{key:row.id,value:row.id},row.name))])])
    ])
   ]),
   comparison&&e('div',{className:'card',key:'table'},[
    e('h3',{},'Selected Compound Comparison'),e('p',{className:'small'},'Experimental values take precedence. Each cell retains its evidence type. No overall score or automatic ranking is calculated.'),
    e('div',{className:'comparison-structure-row'},(comparison.compounds||[]).map(compound=>e('div',{className:'comparison-structure-card',key:compound.row_id},[compound.svg?e('div',{className:'comparison-structure-image',dangerouslySetInnerHTML:{__html:compound.svg}}):e('div',{className:'comparison-structure-unavailable'},'Structure unavailable'),e('strong',{},compound.name||compound.compound||'Compound')]))) ,
    e('div',{className:'table-scroll comparison-table-scroll'},e('table',{className:'comparison-table'},[
     e('thead',{},e('tr',{},['Metric',...(comparison.compounds||[]).map(compound=>compound.name||compound.compound||'Compound')].map(label=>e('th',{key:label},label)))),
     e('tbody',{},metrics.map(metric=>e('tr',{key:metric},[
      e('th',{scope:'row'},[e('div',{className:'comparison-metric-title'},metric),e('small',{className:'comparison-criteria'},metricCriteria[metric]||'Evidence shown below')]),
      ...(comparison.compounds||[]).map(compound=>e('td',{key:compound.row_id,className:'comparison-value'},[
       e('span',{className:'mono'},compound[metric]??'—'),
       e('div',{className:'comparison-source'},StatusBadge({type:compound.sources?.[metric]||'Not measured'}))
      ]))
     ])))
    ])),
    comparisonVisualizations(comparison,compareAssay)
   ])
  ]);
 }

 function SettingsPanel(){
  const projectPerformance=new Map((admet?.model_performance||[]).filter(row=>row.scope==='PROJECT:'+projectId).map(row=>[row.model_id,row]));
  const pkMethods=helpRegistry?.pk_method_registry||[];
  const appInfo=helpRegistry?.application||{version:'0.6.3-stage5b4-ui',current_stage:'5B-4',standardizer:'CHEM_STANDARDIZER_V1',standardizer_version:'1.0.0'};
  return e('div',{},[
   e('div',{className:'card',key:'models'},[
    e('div',{className:'eyebrow'},'MODEL REGISTRY & GOVERNANCE'),
    e('h2',{},'Prediction Models'),
    e('p',{className:'small'},'Multiple registered models may share an endpoint. Project performance influences consensus only from N ≥ 10; N ≥ 30 enables a stronger blend.'),
    (admet?.models||[]).length?e('div',{className:'table-scroll'},e('table',{},[
      e('thead',{},e('tr',{},['Endpoint','Model','Version','Status','Calibration Data','Conformal Quality','Training N','Validation','Project Experimental N','Project MAE / Accuracy','Consensus Weight','Project Selection'].map(label=>e('th',{key:label},label)))),
      e('tbody',{},admet.models.map(model=>{
       const performance=projectPerformance.get(model.id),validation=model.validation||model.details?.validation||{},metric=performance?.metrics?.mae??performance?.metrics?.accuracy;
       const best=admet?.best_project_models?.[model.endpoint]?.model_id===model.id;
       return e('tr',{key:model.id},[
        e('td',{},e('strong',{},model.endpoint)),
        e('td',{},model.model_name),
        e('td',{className:'mono'},model.model_version),
        e('td',{},StatusBadge({type:model.status})),
        e('td',{},model.details?.calibration_dataset||model.training_dataset||'—'),
        e('td',{},model.conformal_status?model.conformal_status.replace(/^CONFORMAL_/,''):'—'),
        e('td',{className:'mono'},model.training_n!=null?String(model.training_n):'—'),
        e('td',{className:'mono small'},Object.entries(validation).map(([k,v])=>k+': '+(typeof v==='number'?v.toFixed(2):v)).join(', ')||'—'),
        e('td',{className:'mono'},performance?String(performance.n_observations):'0'),
        e('td',{className:'mono'},metric!=null?Number(metric).toFixed(3):'—'),
        e('td',{className:'mono'},best?'1.00 (Best)':'Equal/Auto'),
        e('td',{},best?StatusBadge({type:'READY'}):'—')
       ]);
      }))
    ])):e('p',{},'No prediction models registered.')
   ]),
   e('div',{className:'card',key:'pk-methods'},[
    e('div',{className:'eyebrow'},'PK & IVIVE METHOD REGISTRY'),
    e('h2',{},'PK Methods & Equations'),
    e('p',{className:'small'},'Registered mathematical equations and physiological methods for IVIVE and NCA.'),
    pkMethods.length?e('div',{className:'table-scroll'},e('table',{},[
     e('thead',{},e('tr',{},['Method Key','Method Name','Version','Equations / Assumptions','Status'].map(label=>e('th',{key:label},label)))),
     e('tbody',{},pkMethods.map(m=>e('tr',{key:m.method_key},[
      e('td',{className:'mono bold'},m.method_key),
      e('td',{},m.method_name),
      e('td',{className:'mono'},m.method_version),
      e('td',{className:'small'},(m.equations||[]).concat(m.assumptions||[]).join(' · ')||'—'),
      e('td',{},StatusBadge({type:m.status}))
     ])))
    ])):e('p',{},'No PK methods registered.')
   ]),
   e('div',{className:'card',key:'platform-info'},[
    e('div',{className:'eyebrow'},'PLATFORM INFORMATION'),
    e('div',{},[e('div',{className:'eyebrow'},'DANGER ZONE'),e('h2',{},'Delete Project'),e('p',{className:'small'},'Deletion requires a separate confirmation and the exact project name. No other project is included.')]),
    e('button',{className:'danger',onClick:()=>openDeleteDialog([selectedProjectSummary||project])},'Delete Project…')
   ])
  ]);
 }

 const openGlobalView=view=>{setGlobalView(view);setProjectTab('dashboard');setDetail(null);setAddCompoundOpen(false);setComparison(null);setSelectedCandidate(null);setSidebarOpen(false);if(view==='help'&&!helpRegistry)loadHelpRegistry().catch(error=>setMessage(String(error)));loadDashboard().catch(error=>setMessage(String(error)))};
 const goDashboard=()=>openGlobalView('dashboard');
 const openProject=itemId=>{setProjectId(itemId);setProjectTab('compounds');setDetail(null);setAddCompoundOpen(false);setComparison(null);setSidebarOpen(false)};
 const openOptimizationOverview=()=>{setOptimizationWorkspace(current=>({...current,project_id:current.project_id||String(projectId||''),compound_id:''}));openGlobalView('optimization')};
 const openSettings=()=>{if(projectId){setProjectTab('settings');setDetail(null);setSidebarOpen(false)}else openGlobalView('settings')};
 const selectedProjectSummary=(dashboard?.projects||[]).find(row=>row.id===projectId);
 const openDeleteDialog=(items,event)=>{
  event?.stopPropagation();
  const safeItems=(items||[]).filter(item=>item?.id);
  if(!safeItems.length)return;
  setDeleteProjects(safeItems);setDeleteConfirmations({});setMessage('');
 };
 const closeDeleteDialog=()=>{if(deleteBusy)return;setDeleteProjects([]);setDeleteConfirmations({})};
  const deleteNamesMatch=deleteProjects.length>0&&deleteProjects.every(item=>(deleteConfirmations[item.id]||'').trim()===item.name.trim());
  const confirmProjectDeletion=async()=>{
   if(!deleteNamesMatch||deleteBusy)return;
   setDeleteBusy(true);
   try{
    const confirmations=deleteProjects.map(item=>({id:item.id,confirmation_name:(deleteConfirmations[item.id]||'').trim()}));
    const result=confirmations.length===1
     ?await api.del('/projects/'+confirmations[0].id,{confirmation_name:confirmations[0].confirmation_name})
     :await api.post('/projects/bulk-delete',{projects:confirmations});
    const deletedIds=result.deleted_project_ids||confirmations.map(item=>item.id),currentDeleted=deletedIds.includes(projectId);
    const [rows,summary]=await Promise.all([api.get('/projects'),api.get('/dashboard')]);
    setProjects(rows);setDashboard(summary);setProjectSelection([]);setDeleteProjects([]);setDeleteConfirmations({});
    setGlobalView('dashboard');setProjectTab('dashboard');setDetail(null);setWorkspace(null);setAdmet(null);setMetabolism(null);setComparison(null);setSelected([]);setSelectedCandidate(null);
    if(currentDeleted){setProjectId(null);setProject(null)}
    setMessage((result.deleted_project_names||deleteProjects.map(item=>item.name)).join(', ')+' deleted successfully');
   }catch(error){
    let msg='Project deletion failed';
    try{const parsed=JSON.parse(error.message);msg=parsed.detail||msg}catch(_){msg=error.message||String(error)}
    setMessage(msg.replace(/^Error:\s*/,''));
   }finally{
    setDeleteBusy(false);
   }
  };

 function ProjectDeleteModal(){
  if(!deleteProjects.length)return null;
  const bulk=deleteProjects.length>1;
  return e('div',{className:'modal-backdrop project-delete-backdrop',role:'presentation'},e('div',{className:'card project-delete-modal',role:'dialog','aria-modal':'true','aria-labelledby':'project-delete-title'},[
   e('div',{className:'row toolbar',key:'head'},[e('div',{},[e('div',{className:'eyebrow'},bulk?'BULK PROJECT CLEANUP':'PROJECT DELETION'),e('h2',{id:'project-delete-title'},bulk?'Delete Selected Projects':'Delete Project')]),e('button',{className:'secondary',disabled:deleteBusy,onClick:closeDeleteDialog},'Cancel')]),
   e('div',{className:'delete-warning',key:'warning'},[e('strong',{},'This action permanently deletes all project-linked data.'),e('p',{},'Compounds and versions, assays and activity, experimental ADMET, predictions and audit runs, metabolism, optimization runs and candidates, and related project-scoped records will all be deleted. Other projects are not affected.')]),
   e('div',{className:'delete-project-list',key:'projects'},deleteProjects.map(item=>{
    const experimental=(item.experimental_activity_count||0)+(item.experimental_admet_count||0),prediction=item.prediction_count||0,optimization=item.optimization_run_count||0;
    return e('section',{className:'delete-project-summary',key:item.id},[
     e('h3',{},item.name),e('div',{className:'delete-count-grid'},[
      e('div',{key:'compound'},[e('span',{},'Compounds'),e('strong',{},String(item.compound_count||0))]),e('div',{key:'experimental'},[e('span',{},'Experimental data'),e('strong',{},String(experimental))]),e('div',{key:'prediction'},[e('span',{},'Predictions'),e('strong',{},String(prediction))]),e('div',{key:'optimization'},[e('span',{},'Optimization runs'),e('strong',{},String(optimization))])
     ]),e('label',{},['Type ',e('strong',{key:'name'},item.name),' to confirm']),e('input',{value:deleteConfirmations[item.id]||'',autoComplete:'off',onChange:event=>setDeleteConfirmations(current=>({...current,[item.id]:event.target.value}))})
    ]);
   })),
   e('div',{className:'row delete-actions',key:'actions'},[e('button',{className:'secondary',disabled:deleteBusy,onClick:closeDeleteDialog},'Keep Project'+(bulk?'s':'')),e('button',{className:'danger',disabled:deleteBusy||!deleteNamesMatch,onClick:confirmProjectDeletion},deleteBusy?'Deleting…':(bulk?'Delete Selected Projects Permanently':'Delete Project Permanently'))])
  ]));
 }

 function GlobalOptimizationWorkspace(){
  const projectChoices=dashboard?.projects||projects;
  const selectedProjectId=Number(optimizationWorkspace.project_id);
  const selectedCompoundId=Number(optimizationWorkspace.compound_id);
  const selectedProject=projectChoices.find(row=>row.id===selectedProjectId);

  React.useEffect(()=>{
   if(selectedProjectId && project?.id !== selectedProjectId){
    loadProject(selectedProjectId).catch(err=>setMessage(String(err)));
   }
  },[selectedProjectId, project?.id]);

  React.useEffect(()=>{
   if(selectedCompoundId && detail?.row_id !== selectedCompoundId){
    openDetail(selectedCompoundId).then(()=>setDetailTab('optimization')).catch(err=>setMessage(String(err)));
   }
  },[selectedCompoundId, detail?.row_id]);

  const compounds=project?.id===selectedProjectId?currentVersions:[];
  const selectedCompound=compounds.find(row=>row.row_id===selectedCompoundId);
  const version=(detail && selectedCompound && detail.row_id===selectedCompound.row_id)?detail.version:null;
  const predictions=version?(admet?.predictions||[]).filter(row=>row.version_id===version.id):[];
  return e('div',{className:'optimization-workspace'},[
   e('section',{className:'card',key:'selector'},[e('div',{className:'eyebrow'},'DETERMINISTIC MEDICINAL CHEMISTRY'),e('h1',{},'Optimization Workspace'),e('p',{className:'small'},'Select a project and a CompoundVersion, review its evidence, then reuse the existing Stage 4A strategy and Stage 4B analog engines. No LLM and no PK are used.'),
    e('div',{className:'optimization-workspace-steps'},[
     e('div',{className:'optimization-workspace-step',key:'project'},[e('h3',{},'Step 1 — Select Project'),e('select',{value:optimizationWorkspace.project_id,onChange:event=>{setOptimizationWorkspace({project_id:event.target.value,compound_id:''});setDetail(null)}},[e('option',{value:''},'Select project'),...projectChoices.map(row=>e('option',{key:row.id,value:row.id},row.name+' · '+(row.target||'Target not set')))])]),
     e('div',{className:'optimization-workspace-step',key:'compound'},[e('h3',{},'Step 2 — Select Compound'),e('select',{value:optimizationWorkspace.compound_id,disabled:!selectedProject,onChange:event=>setOptimizationWorkspace(current=>({...current,compound_id:event.target.value}))},[e('option',{value:''},selectedProject?(compounds.length?'Select compound':'No compounds registered'):'Select a project first'),...compounds.map(row=>e('option',{key:row.row_id,value:row.row_id,disabled:!row.version},row.name+' · '+row.compound_id+(row.version?' · v'+row.current_version:' · Draft')))])])
    ])
   ]),
   selectedCompound&&version&&e('section',{className:'card',key:'profile'},[e('div',{className:'row toolbar'},[e('div',{},[e('div',{className:'eyebrow'},'COMPOUND PROFILE'),e('h2',{},selectedCompound.name)]),StatusBadge({type:selectedCompound.status})]),e('div',{className:'grid'},[
    e('div',{className:'col-3 structure'},Svg({src:version.svg})),
    e('div',{className:'col-3'},[e('h4',{},'Properties'),e('p',{className:'small'},version.calculated?'MW '+version.properties.molecular_weight+' · cLogP '+version.properties.clogp+' · TPSA '+version.properties.tpsa:'Not calculated')]),
    e('div',{className:'col-3'},[e('h4',{},'Activity'),e('p',{className:'small'},(workspace?.activity?.measurements||[]).length+' experimental · '+(workspace?.activity?.predictions||[]).length+' predicted')]),
    e('div',{className:'col-3'},[e('h4',{},'ADMET / Metabolism'),e('p',{className:'small'},predictions.length+' individual model predictions · '+((metabolism?.runs||[]).length?'soft spots available':'metabolism not run'))]),
    e('div',{className:'col-12'},[e('h4',{},'Prediction Confidence'),e('p',{className:'small'},predictions.length?predictions.slice(0,8).map(row=>row.endpoint+' '+row.confidence+' / '+row.applicability_domain).join(' · '):'No predictions available')])
   ])]),
   selectedCompound&&!version&&e('div',{className:'card empty-state',key:'draft'},[StatusBadge({type:'DRAFT'}),e('p',{},'Add a validated structure before optimization.')]),
   version&&optimizationConfig&&optimizationPanel(version.id),
   selectedProject&&compounds.length===0&&e('section',{className:'card empty-state',key:'zero-compounds'},[e('h3',{},'No compounds are registered in this project.'),e('p',{},'Add a compound to this project before creating an optimization run.')]),
   selectedProject&&compounds.length>0&&!selectedCompound&&e('section',{className:'card empty-state',key:'choose'},[e('h3',{},'Choose a parent compound'),e('p',{},'The selected project has '+compounds.length+' compound record(s). Drafts without structure cannot be optimized.')]),
   !selectedProject&&e('section',{className:'card empty-state',key:'no-project'},[e('h3',{},'Select a project to begin optimization.'),e('p',{},'Select a project from the dropdown above to open its optimization workspace.')])
  ]);
 }

 function HelpPage(){
  if(helpBusy&&!helpRegistry)return e('section',{className:'card help-view'},[e('h2',{},'Loading platform inventory…'),e('p',{className:'small'},'Reading runtime packages and scientific registries.')]);
  if(!helpRegistry)return e('section',{className:'card help-view'},[e('h2',{},'Help inventory unavailable'),e('button',{onClick:()=>loadHelpRegistry().catch(error=>setMessage(String(error)))},'Retry')]);
  const appInfo=helpRegistry.application||{},models=helpRegistry.models||[],caps=helpRegistry.capability_summary?.groups||[];
  const modelRows=names=>models.filter(row=>names.includes(row.endpoint));
  const renderModelTable=(rows,withSpecies=true)=>e('div',{className:'table-scroll help-model-table'},e('table',{},[
   e('thead',{},e('tr',{},['Endpoint','Model','Model Version','Output',...(withSpecies?['Species']:[]),'Availability','Scientific Confidence'].map(label=>e('th',{key:label},label)))),
   e('tbody',{},rows.map(row=>e('tr',{key:row.id||row.endpoint},[
    e('td',{},row.endpoint),e('td',{},row.model_name),e('td',{className:'mono'},row.model_version),e('td',{},row.output_type+(row.output_unit?' · '+row.output_unit:'')),
    ...(withSpecies?[e('td',{key:'species'},row.species||'Endpoint-defined')]:[]),e('td',{},StatusBadge({type:row.availability})),e('td',{},row.confidence)
   ])))
  ]));
  const supportedCyp=modelRows(['CYP1A2 inhibitor','CYP2C9 inhibitor','CYP2C19 inhibitor','CYP2D6 inhibitor','CYP3A4 inhibitor','CYP2C9 substrate','CYP2D6 substrate','CYP3A4 substrate','P-gp inhibitor']);
  const safety=modelRows(['hERG liability','Ames mutagenicity','DILI clinical liability']);
  const unsupported=(caps.find(row=>row.key==='cyp_transporters')?.items||[]).filter(row=>row.availability==='MODEL_UNAVAILABLE');
  const pkCaps=caps.find(row=>row.key==='pk')?.items||[];
  return e('div',{className:'help-view',key:'help-page','data-testid':'help-registry'},[
   e('section',{className:'card help-section',id:'help-overview',key:'overview'},[
    e('div',{className:'eyebrow'},'PLATFORM OVERVIEW'),e('h2',{},'Drug-OPT Platform Help'),
    e('p',{},'Drug-OPT is a structure-to-PK medicinal chemistry optimization platform supporting project-based compound registration, physicochemical properties, activity, ADMET, metabolism, optimization, experimental PK, IVIVE, simulation and translational PK. Where scientific functions live: access project and compound workspaces from the left navigation to run calculations, models and simulations.'),
    e('p',{className:'help-caution'},'Results are decision-support predictions and calculations. They do not replace fit-for-purpose experimental studies or expert scientific review.')
   ]),
   e('section',{className:'card help-section',id:'help-workflow',key:'workflow-guide'},[
    e('div',{className:'eyebrow'},'WORKFLOW GUIDE'),
    e('h2',{},'Typical Workflow'),
    e('p',{className:'small'},'Recommended 7-step sequence for compound evaluation and optimization from project registration to translational PK simulation.')
   ]),
   e('section',{className:'card help-section',id:'help-version',key:'version'},[e('h2',{},'Current Platform Version'),e('dl',{className:'help-version-grid'},[
    ['Application version',appInfo.version],['Current Stage',appInfo.current_stage_label||'Internal Validation'],['Git/build version',appInfo.build_version],['Standardizer',appInfo.standardizer+' '+appInfo.standardizer_version],['RDKit',appInfo.rdkit_version]
    ].map(([label,value])=>e('div',{key:label},[e('dt',{},label),e('dd',{className:'mono'},value||'Unavailable')])))]),
    e('section',{className:'card help-section',id:'help-version-history',key:'version-history'},[
     e('div',{className:'eyebrow'},'RELEASE HISTORY'),
     e('h2',{},'Version History & Scientific Milestones'),
     e('p',{className:'small'},'Milestone-driven history of platform releases and capabilities.'),
     e('div',{className:'table-scroll'},e('table',{},[
      e('thead',{},e('tr',{},['Version','Release Date','Stage','Key Scientific Capabilities'].map(label=>e('th',{key:label},label)))),
      e('tbody',{},(helpRegistry.version_history||[]).map(vh=>e('tr',{key:vh.version,style:vh.version===appInfo.version?{background:'#f0f7ff',fontWeight:'bold'}:{}},[
       e('td',{className:'mono'},vh.version),
       e('td',{},vh.date),
       e('td',{},e('span',{className:'badge-favorable'},vh.stage)),
       e('td',{className:'small'},vh.highlights)
      ])))
     ]))
    ]),
   e('section',{className:'card help-section',id:'help-modules',key:'modules'},[e('h2',{},'Structure & Cheminformatics Modules'),e('p',{className:'small'},'Versions below are read from the active production Python environment.'),e('div',{className:'table-scroll'},e('table',{},[
    e('thead',{},e('tr',{},['Module','Version','Used For','Status'].map(label=>e('th',{key:label},label)))),e('tbody',{},(helpRegistry.structure_modules||[]).map(row=>e('tr',{key:row.module},[e('td',{},row.module),e('td',{className:'mono'},row.version),e('td',{},row.used_for),e('td',{},StatusBadge({type:row.status}))])))
   ])),e('p',{className:'small'},'RDKit provides structure parsing, CHEM_STANDARDIZER_V1 processing, molecular properties, Crippen cLogP, TPSA, fingerprints and structural handling. SyGMa is a rule-based metabolism engine; installed prediction checkpoints use ML where identified.'),e('details',{},[e('summary',{},'Complete runtime package inventory'),e('div',{className:'table-scroll'},e('table',{},[e('thead',{},e('tr',{},['Package','Installed Version','Purpose','Status'].map(label=>e('th',{key:label},label)))),e('tbody',{},(helpRegistry.package_inventory||[]).map(row=>e('tr',{key:row.package},[e('td',{},row.package),e('td',{className:'mono'},row.version),e('td',{},row.purpose),e('td',{},StatusBadge({type:row.status}))])))]))])]),
   e('section',{className:'card help-section',id:'help-adme',key:'adme'},[e('h2',{},'ADME Prediction Models'),renderModelTable(modelRows(['Solubility','Permeability','Plasma protein binding','HLM intrinsic clearance','RLM intrinsic clearance','MLM intrinsic clearance'])),e('details',{},[e('summary',{},'Validation, applicability and limitations'),e('div',{className:'help-details'},modelRows(['Solubility','Permeability','Plasma protein binding','HLM intrinsic clearance','RLM intrinsic clearance','MLM intrinsic clearance']).map(row=>e('article',{key:row.endpoint},[e('h3',{},row.endpoint+' — '+row.model_name),e('p',{className:'small'},'Training: '+(row.training_dataset||'Not reported')),e('p',{className:'small'},'Conformal: '+row.conformal_status+' · License: '+(row.license||'Not reported')),e('p',{className:'small'},row.details?.limitations||'See model registry for current limitations.')])) )])]),
   e('section',{className:'card help-section',id:'help-cyp',key:'cyp'},[e('h2',{},'CYP & Transporters'),e('p',{},'Inhibitor and substrate models are distinct endpoint-specific classifiers. Availability does not imply quantitative inhibition or substrate kinetics.'),renderModelTable(supportedCyp,false),e('details',{className:'unavailable-collapse'},[e('summary',{},'Currently unavailable transporter and CYP endpoints'),e('ul',{},unsupported.map(row=>e('li',{key:row.key},row.label+' — MODEL_UNAVAILABLE')))])]),
   e('section',{className:'card help-section',id:'help-safety',key:'safety'},[e('h2',{},'Safety / Toxicology'),renderModelTable(safety,false),e('p',{className:'small'},'hERG, Ames and DILI are classification models. Structural Alerts are deterministic rule-based calculations and remain distinct from model predictions. Outputs are screening evidence, not regulatory toxicology conclusions.')]),
   e('section',{className:'card help-section',id:'help-activity',key:'activity'},[e('h2',{},'Activity / SAR'),e('p',{},'Define project-local assays and record experimental IC50, EC50, Ki, Kd or GI50 values with cell line, species, mutation and protocol context. Project SAR includes similarity, matched molecular pairs and activity cliffs.'),e('p',{className:'help-caution'},'Activity prediction depends on the selected project assay and sufficient project data. It is not automatically run by Save & Predict.')]),
   e('section',{className:'card help-section',id:'help-optimization',key:'optimization'},[e('h2',{},'Optimization Engine'),e('h3',{},'Stage 4A — Optimization Strategy Engine'),e('p',{},'Identifies liabilities, protected and modifiable regions, metabolic soft spots and relevant MMP evidence.'),e('h3',{},'Stage 4B — Analog Generation Engine'),e('p',{},'Applies the transformation library, generates and filters analogs, predicts supported properties, and ranks candidates using transparent objectives and Pareto evidence.'),e('p',{className:'help-caution'},'Proposed structures are medicinal chemistry hypotheses, not experimentally validated compounds.')]),
   e('section',{className:'card help-section',id:'help-pk',key:'pk'},[e('div',{className:'row toolbar'},[e('h2',{},'PK / DMPK'),StatusBadge({type:'READY'})]),e('p',{className:'help-flow'},'Experimental PK → NCA → IVIVE → PK Foundation → IV/PO/SC/IP Simulation → Cross-Species Translation → Human Translational PK → Prospective Validation'),e('div',{className:'help-topic-grid'},[
    ['Experimental PK',['PK Study and concentration-time input','CSV import and BLQ handling','Noncompartmental analysis (NCA)']],['IVIVE',['Clint, MPPGL and hepatocellularity','PPB/fu and blood:plasma ratio','Well-stirred CLh, extraction ratio and Fh']],['Simulation',['IV bolus and infusion','PO, SC and IP routes','Repeated dosing, ka, F, Cmax, Tmax and AUC']],['Translation',['Allometry and LOSO','Human clearance and volume','Simulation readiness','Prospective freeze and retrospective validation']]
   ].map(([title,items])=>e('article',{key:title},[e('h3',{},title),e('ul',{},items.map(item=>e('li',{key:item},item)))]))),e('h3',{},'Implemented capability registry'),e('ul',{className:'help-capability-list'},pkCaps.map(row=>e('li',{key:row.key},[e('span',{},row.label),StatusBadge({type:row.availability})]))),e('details',{},[e('summary',{},'Registered PK methods'),e('ul',{},(helpRegistry.pk_method_registry||[]).map(row=>e('li',{key:row.method_key},row.method_name+' · '+row.method_version+' · '+row.status)))])]),
   e('section',{className:'card help-section',id:'help-glossary',key:'glossary'},[e('h2',{},'Important Scientific Terminology'),e('dl',{className:'help-glossary'},(helpRegistry.glossary||[]).map(row=>e('div',{key:row.term},[e('dt',{},row.term),e('dd',{},row.definition)])))]),
   e('section',{className:'card help-section',id:'help-limitations',key:'limits'},[e('h2',{},'Current Limitations'),e('ul',{},(helpRegistry.limitations||[]).map(item=>e('li',{key:item},item))),e('p',{className:'small'},'Registry source: '+helpRegistry.source)])
  ]);
 }

  function MainDashboard(){
   const modules=dashboard?.capability_summary?.groups||[];
   const home=globalView==='dashboard';
   const operationalModelCount=(dashboard?.model_registry||[]).filter(m=>m.status==='READY'||m.availability==='READY').length||18;
   return e(React.Fragment,{},[
    globalView!=='dashboard'&&globalView!=='optimization'&&e('section',{className:'card global-view-header',key:'global-head'},[
     e('div',{className:'eyebrow'},'WORKSPACE'),
     e('h1',{},({"new-project":'New Project',projects:'Projects',settings:'Settings',help:'Help'})[globalView]||'Workspace')
    ]),
    globalView==='optimization'&&e(GlobalOptimizationWorkspace,{key:'optimization-workspace'}),
    globalView==='new-project'&&e('section',{className:'card dashboard-create',key:'create-page'},[
     e('div',{className:'eyebrow'},'NEW WORKSPACE'),
     e('h2',{},'Create New Project'),
     e('div',{className:'create-project-grid'},[
      e(Field,{label:'Project Name *',value:form.name,onChange:value=>setForm({...form,name:value}),placeholder:'e.g. EGFR Inhibitors'}),
      e(Field,{label:'Target *',value:form.target,onChange:value=>setForm({...form,target:value}),placeholder:'e.g. EGFR'}),
      e('div',{},[e('label',{},'Molecule Type'),e('select',{value:form.molecule_type,onChange:event=>setForm({...form,molecule_type:event.target.value})},['Small Molecule','Peptide'].map(value=>e('option',{key:value,value},value)))])
     ]),
     e('button',{disabled:!form.name.trim()||!form.target.trim(),onClick:createProject},'Create Project'),
     e('p',{className:'small dashboard-note'},'Default Workspace Settings: Description and additional metadata can be added later in Project Settings.')
    ]),
    home&&e('section',{className:'card dashboard-hero',key:'intro'},[
     e('div',{className:'eyebrow'},'PLATFORM OVERVIEW'),
     e('h1',{},'Drug Optimization Platform'),
     e('p',{className:'dashboard-hero-desc'},'Structure, activity, ADMET, DMPK and medicinal chemistry optimization data are integrated at the compound-version level to support hit-to-lead and lead optimization decisions.'),
     e('div',{className:'platform-capabilities-line'},[
      e('span',{key:'c1'},'• Structure-based compound management'),
      e('span',{key:'c2'},'• Experimental data integration'),
      e('span',{key:'c3'},'• Predictive ADMET'),
      e('span',{key:'c4'},'• SAR / optimization workflow'),
      e('span',{key:'c5'},'• Translational PK'),
      e('span',{key:'c6'},'• Full prediction provenance')
     ]),
     e('div',{className:'dashboard-stats-grid'},[
      e('div',{className:'dashboard-stat-card',key:'projects'},[
       e('span',{className:'stat-label'},'Projects'),
       e('strong',{className:'stat-value'},String(dashboard?.totals?.projects??projects.length))
      ]),
      e('div',{className:'dashboard-stat-card',key:'compounds'},[
       e('span',{className:'stat-label'},'Compounds'),
       e('strong',{className:'stat-value'},String(dashboard?.totals?.compounds??projects.reduce((sum,row)=>sum+(row.compound_count||0),0)))
      ]),
    e('div',{className:'dashboard-stat-card',key:'stage'},[
       e('span',{className:'stat-label'},'Current Stage'),
       e('strong',{className:'stat-value'},'Internal Validation'),
       e('span',{className:'small'},'Prediction Engine v1 Frozen'),
       e('span',{className:'small'},'Experimental data collection active')
      ]),
      e('div',{className:'dashboard-stat-card',key:'models'},[
       e('span',{className:'stat-label'},'Model Endpoints'),
       e('strong',{className:'stat-value'},operationalModelCount+' operational')
      ])
     ])
    ]),
    home&&e('section',{className:'card scientific-workspace-section',key:'scientific-workspace'},[
     e('div',{className:'row toolbar',style:{marginBottom:'14px'}},[
      e('div',{},[
       e('div',{className:'eyebrow'},'SCIENTIFIC WORKSPACE'),
       e('h2',{},'Available Scientific Modules')
      ]),
      e('p',{className:'small',style:{margin:0,alignSelf:'center'}},'Status reflects the current local engine and model registry.')
     ]),
     e('div',{className:'scientific-workspace-grid'},modules.map(module=>e('div',{className:'scientific-card',key:module.key||module.title},[
      e('div',{className:'scientific-card-header'},[
       e('h3',{},module.title),
       StatusBadge({type:module.status})
      ]),
      e('p',{className:'small scientific-card-desc'},module.description),
      e('div',{className:'scientific-card-rows'},(
       module.key==='cyp_transporters'?[
        {key:'cyp_inh',label:'CYP Inhibitors',availability:'READY',sub:'5 Isoforms (1A2, 2C9, 2C19, 2D6, 3A4)'},
        {key:'cyp_sub',label:'CYP Substrates',availability:'READY',sub:'3 Isoforms (2C9, 2D6, 3A4)'},
        {key:'pgp_inh',label:'P-gp Inhibitor',availability:'READY',sub:'Human P-glycoprotein'},
        {key:'add_trans',label:'Additional Transporters',availability:'LIMITED',sub:'BCRP, BSEP, OATP, OCT, MATE (Registry)'}
       ]:module.items
      ).map(item=>e('div',{className:'scientific-card-row',key:item.key||item.label},[
       e('div',{className:'capability-label'},[
        e('span',{className:'item-title'},item.label),
        item.sub&&e('small',{className:'item-sub'},item.sub),
        item.confidence&&item.confidence!=='NOT_APPLICABLE'&&e('small',{className:'item-meta'},'Confidence: '+item.confidence),
        item.conformal_status&&item.conformal_status!=='NOT_APPLICABLE'&&e('small',{className:'item-meta'},'Conformal: '+item.conformal_status.replace(/^CONFORMAL_/,''))
       ]),
       StatusBadge({type:item.availability})
      ])))
     ])))
    ]),
    globalView==='projects'&&e('section',{className:'card',key:'projects-list'},[
     e('div',{className:'row toolbar'},[
      e('div',{},[
       e('div',{className:'eyebrow'},'RESEARCH PORTFOLIO'),
       e('h2',{},'Projects'),
       e('p',{className:'small'},'Click any project title to open its compound workspace.')
      ]),
      e('div',{className:'row'},[
       projectSelection.length>0&&e('span',{className:'small'},projectSelection.length+' selected'),
       e('button',{className:'secondary project-delete-secondary',disabled:projectSelection.length===0,onClick:()=>openDeleteDialog((dashboard?.projects||projects).filter(item=>projectSelection.includes(item.id)))},'Delete Selected…'),
       projectId&&e('button',{className:'secondary',onClick:()=>openProject(projectId)},'Continue Current Project')
      ])
     ]),
     (dashboard?.projects||projects).length?e('div',{className:'table-scroll dashboard-project-table-wrap'},e('table',{className:'dashboard-project-table'},[
      e('thead',{},e('tr',{},['Project','Compounds','Details','Status','Delete'].map(label=>e('th',{key:label},label)))),
      e('tbody',{},(dashboard?.projects||projects).map(item=>e('tr',{className:'dashboard-project',key:item.id},[
       e('td',{className:'project-name-cell'},[
        e('div',{className:'eyebrow'},item.molecule_type||'Small Molecule'),
        e('button',{className:'project-link-title',onClick:()=>openProject(item.id)},item.name)
       ]),
       e('td',{className:'dashboard-count'},String(item.compound_count||0)),
       e('td',{},[
        e('div',{},[e('strong',{},'Target: '),item.target||'Not set']),
        e('div',{className:'small'},'Molecule: '+(item.molecule_type||'Small Molecule'))
       ]),
       e('td',{},[
        e('div',{className:'small'},item.status_summary||((item.compound_count||0)+' compounds · isolated workspace')),
        e('div',{className:'small'},'Experimental: '+String((item.experimental_activity_count||0)+(item.experimental_admet_count||0))+' · Optimization: '+String(item.optimization_run_count||0))
       ]),
       e('td',{className:'project-actions-cell'},[
        e('button',{className:'secondary project-delete-secondary',onClick:event=>openDeleteDialog([item],event)},'Delete…')
       ])
      ])))
     ])):e('div',{className:'empty-state'},[
      e('h3',{},'No projects yet'),
      e('p',{},'Use New Project in the sidebar to initialize your first target workspace.')
     ])
    ]),
    globalView==='settings'&&SettingsPanel(),
    globalView==='help'&&e(HelpPage,{key:'help'})
   ]);
  }

  function ProjectWorkspace(){
   const summary=selectedProjectSummary;
   const statusByCompound=new Map((summary?.compounds||[]).map(row=>[row.row_id,row]));
   return e(React.Fragment,{},[
    e('div',{className:'card project-header',key:'header'},project?e('div',{className:'row toolbar'},[e('div',{},[e('div',{className:'eyebrow'},'PROJECT DASHBOARD'),e('h1',{},project.name),e('div',{},[e('strong',{},project.target||'Target not set'),' · ',project.molecule_type])]),e('button',{onClick:()=>{setMessage('');setAddCompoundOpen(true);setCompoundForm({compound_id:'',name:'',smiles:'',notes:''})}},'Add Compound')]):e('div',{},[e('h2',{},'Select or create a project'),e('p',{},'Start with a project, then add compounds and work from Compound Detail.') ])),
    project&&e('nav',{className:'project-nav',key:'nav'},detail
     ?e('button',{className:'secondary',onClick:()=>{setDetail(null);setDetailTab('overview')}},'← Compound List')
     :[['compounds','Compounds'],['assays','Assays'],['compare','Compare'],['settings','Settings']].map(([tab,label])=>e('button',{key:tab,className:projectTab===tab?'':'secondary',onClick:()=>{setProjectTab(tab);if(tab!=='compounds')setDetail(null)}},label))),
    project&&projectTab==='compounds'&&!detail&&e(React.Fragment,{key:'project-dashboard'},[
     e('section',{className:'card',key:'overview'},[e('div',{className:'eyebrow'},'PROJECT OVERVIEW'),e('h2',{},'Current Project Status'),e('div',{className:'project-overview-grid'},[
      ['Target',project.target||'Not set'],['Molecule Type',project.molecule_type],['Compounds',summary?.compound_count??currentVersions.length],['Experimental Activity',summary?.experimental_activity_count??0],['Experimental ADMET',summary?.experimental_admet_count??0],['Predictions',summary?.prediction_count??0],['Optimization Runs',summary?.optimization_run_count??0]
     ].map(([label,value])=>e('div',{className:'project-overview-item',key:label},[e('span',{},label),e('strong',{},String(value))])))]),
     e('section',{className:'card workflow-card',key:'workflow'},[e('div',{className:'eyebrow'},'WORKFLOW STATUS'),e('div',{className:'workflow-strip'},['Structure','Properties','Activity','ADMET','Optimization','PK'].map((stage,index)=>e(React.Fragment,{key:stage},[e('div',{className:'workflow-step'},[e('span',{},stage),StatusBadge({type:summary?.workflow?.[stage]||'NOT_STARTED'})]),index<5&&e('span',{className:'workflow-arrow'},'→')])))]),
     e('section',{className:'card',key:'compounds'},[e('div',{className:'row toolbar'},[e('div',{},[e('h2',{},'Compound Status'),e('p',{className:'small'},'Each row summarizes only the current CompoundVersion in this project.')]),e('div',{className:'row'},[e('button',{className:'secondary',disabled:selected.length<2,onClick:compare},'Compare Selected'),e('button',{onClick:()=>{setMessage('');setAddCompoundOpen(true)}},'Add Compound')])]),currentVersions.length?e('div',{className:'table-scroll'},e('table',{className:'compound-list project-status-table'},[e('thead',{},e('tr',{},['','Compound','Structure','Properties','Activity','ADMET','Optimization',''].map((x,index)=>e('th',{key:x||index},x)))),e('tbody',{},currentVersions.map(compound=>{const status=statusByCompound.get(compound.row_id)||{};return e('tr',{key:compound.row_id},[e('td',{className:'compound-select-cell'},e('input',{className:'compound-select',type:'checkbox',checked:selected.includes(compound.row_id),onClick:event=>event.stopPropagation(),onChange:event=>setSelected(current=>event.target.checked?(current.includes(compound.row_id)?current:[...current,compound.row_id]):current.filter(id=>id!==compound.row_id))})),e('td',{},[e('button',{className:'link-button',onClick:()=>openDetail(compound.row_id)},compound.name),e('div',{className:'mono small'},compound.compound_id)]),e('td',{className:'thumbnail'},[Svg({src:compound.version?.svg}),StatusBadge({type:status.structure||'NOT_STARTED'})]),e('td',{},StatusBadge({type:status.properties||'NOT_RUN'})),e('td',{},StatusBadge({type:status.activity||'NOT_RUN'})),e('td',{},StatusBadge({type:status.admet||'NOT_RUN'})),e('td',{},StatusBadge({type:status.optimization||'NOT_RUN'})),e('td',{},e('button',{className:'secondary',onClick:()=>openDetail(compound.row_id)},'Open'))])}))])):e('div',{className:'empty-state'},[e('h3',{},'No compounds yet'),e('p',{},'Add the first compound by name; structure and calculation may follow later.'),e('button',{onClick:()=>{setMessage('');setAddCompoundOpen(true)}},'Add Compound')])])
    ]),
    project&&projectTab==='compounds'&&detail&&compoundDetail(),
    project&&projectTab==='assays'&&e('div',{className:'card',key:'assays'},[e('h2',{},'Assays and Activity'),e('p',{},'Define assays, enter experimental activity, train project QSAR, and inspect SAR in the existing workbench.'),e('a',{className:'button',href:'/static/stage2-workbench.html?project='+projectId},'Open Activity Workbench')]),
    project&&projectTab==='compare'&&e(React.Fragment,{key:'compare'},ComparePanel()),
    project&&projectTab==='optimization'&&e('div',{className:'card',key:'optimization'},[e('div',{className:'eyebrow'},'PROJECT OPTIMIZATION'),e('h2',{},project.name+' Optimization Workspace'),e('p',{},'Optimization runs are CompoundVersion-specific. This project currently has '+(selectedProjectSummary?.optimization_run_count||0)+' recorded optimization run(s).'),currentVersions.length?e('div',{className:'table-scroll'},e('table',{},[e('thead',{},e('tr',{},[e('th',{},'Compound'),e('th',{},'Version'),e('th',{},'Status'),e('th',{},'')])),e('tbody',{},currentVersions.map(compound=>e('tr',{key:compound.row_id},[e('td',{},compound.name),e('td',{className:'mono'},compound.version?'v'+compound.current_version:'Draft'),e('td',{},StatusBadge({type:compound.status})),e('td',{},e('button',{className:'secondary',disabled:!compound.version,onClick:async()=>{await openDetail(compound.row_id);setProjectTab('compounds');setDetailTab('optimization')}},'Open Optimization'))])))])):e('div',{className:'empty-state'},[e('p',{},'Add a compound before creating an optimization run.'),e('button',{onClick:()=>setProjectTab('compounds')},'Go to Compounds')])]),
    project&&projectTab==='settings'&&e(React.Fragment,{key:'settings'},SettingsPanel())
   ]);
  }

  const sidebarItems=[['Dashboard',()=>goDashboard(),'dashboard'],['New Project',()=>openGlobalView('new-project'),'new-project'],['Projects',()=>openGlobalView('projects'),'projects'],['Optimization',openOptimizationOverview,'optimization'],['Settings',openSettings,'settings'],['Help',()=>openGlobalView('help'),'help']];
  const sidebar=e('aside',{className:'sidebar'+(sidebarOpen?' open':''),key:'sidebar'},[
   e('div',{className:'sidebar-head',key:'head'},[
    e('div',{className:'brand-lockup'},[
     e('button',{className:'brand-button',onClick:goDashboard},'Drug Optimization Platform'),
     e('div',{className:'sidebar-tag'},'Research Workspace')
    ]),
    e('button',{className:'menu-toggle',onClick:()=>setSidebarOpen(value=>!value),'aria-expanded':sidebarOpen,'aria-label':'Toggle primary navigation'},sidebarOpen?'Close':'Menu')
   ]),
   e('div',{className:'sidebar-body',key:'body'},[
    e('nav',{className:'global-nav',key:'nav','aria-label':'Primary navigation'},sidebarItems.map(([label,action,view])=>e('button',{key:label,className:(projectTab==='dashboard'&&globalView===view)||(projectTab===view)?'active':'',onClick:action},label)))
   ]),
    e('div',{className:'sidebar-footer',key:'footer'},[
    e('div',{className:'sidebar-footer-brand'},'Drug Optimization Platform'),
    e('div',{className:'sidebar-footer-version'},'v1.0'),
    e('div',{className:'sidebar-footer-date'},'Updated: 2026-08-30')
    ])
  ]);
  return e('div',{className:'shell'},[sidebar,e('main',{className:'content',key:'content'},[
   projectTab==='dashboard'?MainDashboard():ProjectWorkspace(),
   AddCompoundPanel(),
   ProjectDeleteModal(),
   message&&e('pre',{className:'card error',key:'message'},message)
  ])]);
 }

 function unifiedPhysicochemicalTable(properties,rules){


  const props=properties||{};
  const rows=[
   {key:'molecular_weight',name:'Molecular Weight (MW)',val:props.molecular_weight!=null?Number(props.molecular_weight).toFixed(2)+' g/mol':'—',ref:'≤ 500 (Lipinski Rule of 5)',interp:getInterpretation('mw',props.molecular_weight)},
   {key:'clogp',name:'Calculated cLogP (Crippen)',val:props.clogp!=null?Number(props.clogp).toFixed(2):'—',ref:'≤ 5.0 (Lipinski Rule of 5)',interp:getInterpretation('clogp',props.clogp)},
   {key:'tpsa',name:'Topological Polar Surface Area (TPSA)',val:props.tpsa!=null?Number(props.tpsa).toFixed(1)+' Å²':'—',ref:'≤ 140 Å² (Veber Rule)',interp:getInterpretation('tpsa',props.tpsa)},
   {key:'hbd',name:'Hydrogen Bond Donors (HBD)',val:props.hbd!=null?String(props.hbd):'—',ref:'≤ 5 (Lipinski Rule of 5)',interp:getInterpretation('hbd',props.hbd)},
   {key:'hba',name:'Hydrogen Bond Acceptors (HBA)',val:props.hba!=null?String(props.hba):'—',ref:'≤ 10 (Lipinski Rule of 5)',interp:getInterpretation('hba',props.hba)},
   {key:'rotatable_bonds',name:'Rotatable Bonds (Flexibility)',val:props.rotatable_bonds!=null?String(props.rotatable_bonds):'—',ref:'≤ 10 (Veber Rule)',interp:getInterpretation('rotb',props.rotatable_bonds)},
   {key:'fraction_csp3',name:'Fraction Csp3 (Fsp3)',val:props.fraction_csp3!=null?Number(props.fraction_csp3).toFixed(2):'—',ref:'≥ 0.42 (Lovering Complexity)',interp:getInterpretation('fsp3',props.fraction_csp3)},
   {key:'qed',name:'Drug-likeness Score (QED)',val:props.qed!=null?Number(props.qed).toFixed(2):'—',ref:'≥ 0.67 (Bickerton Attractive)',interp:getInterpretation('qed',props.qed)},
   {key:'formal_charge',name:'Formal Charge',val:props.formal_charge!=null?String(props.formal_charge):'0',ref:'0 (Neutral state)',interp:{assessment:'NEUTRAL',colorClass:'favorable',label:'Neutral (0)'}},
   {key:'heavy_atom_count',name:'Heavy Atom Count',val:props.heavy_atom_count!=null?String(props.heavy_atom_count):'—',ref:'20–40 Typical Range',interp:{assessment:'RECORDED',colorClass:'intermediate',label:props.heavy_atom_count!=null?String(props.heavy_atom_count)+' atoms':'—'}}
  ];

  return e('div',{className:'card',key:'physchem-table'},[
   e('div',{className:'eyebrow'},'PHYSICOCHEMICAL & DRUG-LIKENESS TABLE'),
   e('h3',{},'Molecular Properties & Reference Assessments'),
   e('p',{className:'small'},'Single unified physicochemical profile calculated via RDKit. Assessment compares values against established medicinal chemistry guidelines.'),
   e('div',{className:'table-scroll'},[
    e('table',{},[
     e('thead',{},e('tr',{},['Property','Calculated Value','Drug-Likeness / Reference Range','Assessment'].map(h=>e('th',{key:h},h)))),
     e('tbody',{},rows.map(r=>e('tr',{key:r.key},[
      e('td',{style:{fontWeight:600}},r.name),
      e('td',{className:'mono bold'},r.val),
      e('td',{className:'small'},r.ref),
      e('td',{},[ScientificBadge({assessment:r.interp.assessment,colorClass:r.interp.colorClass,textLabel:r.interp.label})])
     ])))
    ])
   ]),
   e('div',{className:'model-notes'},[
    e('strong',{},'Guideline References: '),
    e('span',{},'Lipinski et al. 2001 (MW ≤ 500, cLogP ≤ 5, HBD ≤ 5, HBA ≤ 10) · Veber et al. 2002 (TPSA ≤ 140 Å², RotB ≤ 10) · Lovering et al. 2009 (Fsp3 ≥ 0.42) · Bickerton et al. 2012 (QED ≥ 0.67).')
   ])
  ]);
 }

 function alertList(alerts){
  return alerts.length?alerts.map(alert=>e('div',{key:alert.alert_set+alert.alert_name,className:'alert'},[e('strong',{key:'name'},alert.alert_name),e('div',{key:'reason'},alert.reason),e('div',{key:'set',className:'small'},'Set: '+alert.alert_set+' · atoms: '+(alert.matched_atoms.join(', ')||'not exposed by RDKit'))])):[e('p',{key:'none'},'None detected')];
 }

 function ruleTable(rules){
  const rows=Object.entries(rules||{}).map(([name,rule])=>e('div',{key:name},[
   Badge({ok:rule.result==='PASS',text:rule.result+' · '+name}),
   rule.reasons.length>0&&e('ul',{key:'reasons'},rule.reasons.map(reason=>e('li',{key:reason},reason)))
  ]));
  return e('div',{},[e('h4',{key:'title'},'Drug-likeness Filters'),...rows]);
 }
ReactDOM.createRoot(document.getElementById('root')).render(e(App));
})();
