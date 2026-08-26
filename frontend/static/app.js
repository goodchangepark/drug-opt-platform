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
 del:path=>api.req(path,{method:'DELETE'})
};

function Field({label,value,onChange,type='text',placeholder=''}){
 const tag=type==='textarea'?'textarea':'input';
 const props={value:value??'',onChange:event=>onChange(event.target.value),placeholder};
 if(tag==='input')props.type=type;
 return e('div',{},e('label',{},label),e(tag,props));
}
function Svg({src}){if(!src)return e('div',{className:'structure-placeholder'},'No structure saved');return src.startsWith('data:')?e('img',{src,alt:'structure'}):e('span',{className:'structure',dangerouslySetInnerHTML:{__html:src}})}
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
 ['Solubility','ADME'],['Caco-2 Permeability','ADME'],['Plasma Protein Binding (PPB)','ADME'],
 ['Human Microsomal Stability','ADME'],['Rat Microsomal Stability','ADME'],['Hepatocyte Stability','ADME'],['CYP Inhibition','ADME'],['Transporter','ADME'],
 ['Activity','Activity'],['hERG','Safety'],['Ames','Safety'],['DILI','Safety']
];
const EXPERIMENT_PRESETS={
 'Standard Early ADME':['Solubility','Caco-2 Permeability','Plasma Protein Binding (PPB)','Human Microsomal Stability','Rat Microsomal Stability'],
 'DDI Panel':['CYP Inhibition']
};
const EMPTY_EXPERIMENT={value:'',unit:'',species:'Human',measurement:'',assay:'',role:'Inhibition',isoform:'3A4',transporter:'P-gp',matrix:'',pH:'',medium:'',solubility_type:'',source:'User experimental',notes:'',assay_id:''};

function StatusBadge({type}){
 const labels={Experimental:'EXP',Calculated:'CALC',Predicted:'PRED','Not calculated':'NOT CALCULATED','Not measured':'NOT MEASURED','Not predicted':'NOT PREDICTED','Model unavailable':'MODEL UNAVAILABLE','Not applicable':'NOT APPLICABLE',DRAFT:'DRAFT',STRUCTURE_READY:'STRUCTURE READY',CALCULATED:'CALCULATED',READY:'READY',LIMITED:'LIMITED',MODEL_UNAVAILABLE:'MODEL UNAVAILABLE',PLANNED:'PLANNED',PARTIAL:'PARTIAL',NOT_STARTED:'NOT STARTED',NOT_RUN:'NOT RUN',EXPERIMENTAL:'EXPERIMENTAL',PREDICTED:'PREDICTED'};
 return e('span',{className:'status-badge status-'+String(type||'not-applicable').toLowerCase().replace(/[^a-z]+/g,'-')},labels[type]||type||'NOT APPLICABLE');
}

function App(){
 const [projects,setProjects]=useState([]),[projectId,setProjectId]=useState(null),[project,setProject]=useState(null);
 const [dashboard,setDashboard]=useState(null),[sidebarOpen,setSidebarOpen]=useState(false);
 const [form,setForm]=useState({name:'',target:'',molecule_type:'Small Molecule',description:''});
 const [compoundForm,setCompoundForm]=useState({compound_id:'',name:'',smiles:'',notes:''}),[addCompoundOpen,setAddCompoundOpen]=useState(false);
 const [preview,setPreview]=useState(null),[selected,setSelected]=useState([]),[comparison,setComparison]=useState(null),[detail,setDetail]=useState(null),[message,setMessage]=useState('');
 const [projectTab,setProjectTab]=useState('dashboard'),[detailTab,setDetailTab]=useState('overview');
 const [admet,setAdmet]=useState(null),[admetVersionId,setAdmetVersionId]=useState(''),[admetForm,setAdmetForm]=useState({...EMPTY_ADMET_FORM});
 const [admetCsv,setAdmetCsv]=useState(''),[admetCsvPreview,setAdmetCsvPreview]=useState(null),[admetBusy,setAdmetBusy]=useState(false);
 const [metabolism,setMetabolism]=useState(null),[metabolismBusy,setMetabolismBusy]=useState(false),[metabolicTop,setMetabolicTop]=useState(3),[selectedSpotId,setSelectedSpotId]=useState(null);
 const [metaboliteForm,setMetaboliteForm]=useState({...EMPTY_METABOLITE_FORM});
 const [optimizationConfig,setOptimizationConfig]=useState(null),[optimizationRuns,setOptimizationRuns]=useState([]),[optimizationRun,setOptimizationRun]=useState(null),[optimizationBusy,setOptimizationBusy]=useState(false),[assays,setAssays]=useState([]);
 const [proposalRuns,setProposalRuns]=useState([]),[proposalRun,setProposalRun]=useState(null),[proposalView,setProposalView]=useState('top10'),[selectedCandidate,setSelectedCandidate]=useState(null),[proposalBusy,setProposalBusy]=useState(false);
 const [proposalSettings,setProposalSettings]=useState({max_raw_candidates:120,allow_double_transforms:true}),[userAnalog,setUserAnalog]=useState({smiles:'',reason:''});
 const [optimizationForm,setOptimizationForm]=useState({
  assay_id:'',objectives:['Balanced optimization'],custom_objective:'',
  constraints:{potency_max_nm:'',do_not_worsen_fold:'2',clogp_max:'4',tpsa_min:'40',tpsa_max:'100',mw_max:'550',similarity_min:'0.6',logs_min:'-4',caco2_logpapp_min:'-5.5',herg_do_not_increase:true},endpoint_weights:{}
 });
 const [workspace,setWorkspace]=useState(null),[experimentalOpen,setExperimentalOpen]=useState(false),[experimentalSelected,setExperimentalSelected]=useState([]),[experimentalDrafts,setExperimentalDrafts]=useState({});
 const [compareMetrics,setCompareMetrics]=useState(['MW','cLogP','TPSA','QED','Activity','Solubility','Caco-2','PPB','HLM','RLM','hERG','Ames','DILI']),[compareAssay,setCompareAssay]=useState('');
 const [editorReady,setEditorReady]=useState(false);
 const editorSmiles=useRef('');

 const currentVersions=project?.compounds||[];
 const versionLabel=versionId=>{
  const compound=currentVersions.find(item=>item.version?.id===Number(versionId));
  return compound?compound.compound_id+' v'+compound.current_version:'Version not available';
 };
 const endpointName=endpointId=>(admet?.endpoints||[]).find(item=>item.id===endpointId)?.name||endpointId;

 const loadProjects=async()=>{
  const rows=await api.get('/projects');setProjects(rows);
  if(rows.length&&!projectId)setProjectId(rows[0].id);
 };
 const loadDashboard=async()=>{const data=await api.get('/dashboard');setDashboard(data);return data};
 const loadProject=async id=>{
  const data=await api.get('/projects/'+id);setProject(data);
  setAdmetVersionId(current=>data.compounds.some(item=>item.version?.id===Number(current))?current:(data.compounds.find(item=>item.version)?.version.id||''));
  return data;
 };
 const loadAdmet=async(id=projectId)=>{
  if(!id)return null;
  const data=await api.get('/projects/'+id+'/admet');setAdmet(data);return data;
 };
 const loadMetabolism=async(id=projectId)=>{
  if(!id)return null;
  const data=await api.get('/projects/'+id+'/metabolism');setMetabolism(data);return data;
 };
 const loadWorkspace=async versionId=>{
  if(!versionId){setWorkspace(null);setAdmet(null);setMetabolism(null);return null}
  const data=await api.get('/compound-versions/'+versionId+'/workspace');
  if(data.scope.version_id!==Number(versionId))throw new Error('CompoundVersion isolation check failed');
  setWorkspace(data);setAdmet(data.admet);setMetabolism(data.metabolism);return data;
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

 useEffect(()=>{Promise.all([loadProjects(),loadDashboard()]).catch(error=>setMessage(String(error)))},[]);
 useEffect(()=>{
  setProject(null);setDetail(null);setWorkspace(null);setSelected([]);setComparison(null);setAdmet(null);setMetabolism(null);setAdmetCsvPreview(null);setSelectedSpotId(null);setOptimizationConfig(null);setOptimizationRuns([]);setOptimizationRun(null);setAssays([]);setProposalRuns([]);setProposalRun(null);setSelectedCandidate(null);
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

 const loadSmilesIntoEditor=async()=>{
  if(!compoundForm.smiles.trim())return;
  try{
   const editor=document.getElementById('ketcher-editor')?.contentWindow?.ketcher;
   if(!editor)throw new Error('Structure Editor is still loading');
   await editor.setMolecule(compoundForm.smiles.trim());editorSmiles.current=compoundForm.smiles.trim();setMessage('SMILES loaded into Structure Editor');
  }catch(error){setMessage(String(error))}
 };

 const createProject=async()=>{
  try{
   const created=await api.post('/projects',form);setForm({name:'',target:'',molecule_type:'Small Molecule',description:''});
   await Promise.all([loadProjects(),loadDashboard()]);setProjectId(created.id);setProjectTab('compounds');setMessage('Project created');
  }catch(error){setMessage(String(error))}
 };
 const saveProjectSettings=async()=>{
  try{
   const updated=await api.patch('/projects/'+projectId,{name:project.name,target:project.target,molecule_type:project.molecule_type,description:project.description||''});
   setProject(current=>({...current,...updated}));await Promise.all([loadProjects(),loadDashboard()]);setMessage('Project settings saved');
  }catch(error){setMessage(String(error))}
 };
 const validate=async()=>{try{const result=await api.post('/structure/validate',{smiles:compoundForm.smiles});setPreview(result);setMessage('')}catch(error){setPreview(null);setMessage('Invalid structure: '+error.message)}};
 const saveCompound=async calculate=>{
  try{
   const saved=await api.post('/projects/'+projectId+'/compounds',{...compoundForm,calculate});
   setCompoundForm({compound_id:'',name:'',smiles:'',notes:''});setPreview(null);setAddCompoundOpen(false);
   await Promise.all([loadProject(projectId),loadProjects(),loadDashboard()]);setMessage(calculate?'Compound saved and properties calculated':'Compound saved without calculation');
   await openDetail(saved.row_id);
  }catch(error){setMessage(String(error))}
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
  const base={version_id:detail.version.id,value:qualitative?'':row.value,qualitative_value:qualitative?row.value:'',unit:row.unit,species:row.species,matrix:row.matrix,method:row.measurement,source:row.source,notes:row.notes,provenance:{ui_workflow:'Pre-Stage 5 endpoint selector',display_name:name}};
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
  try{
   const result=await api.post('/metabolism/predict/'+versionId,{});const data=(await loadWorkspace(versionId)).metabolism;
   const run=(data?.runs||[]).find(item=>item.version_id===Number(versionId));setSelectedSpotId(run?.spots?.[0]?.id||null);setMessage(result.message);
  }catch(error){setMessage(String(error))}finally{setMetabolismBusy(false)}
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
  return e('table',{},[
   e('thead',{key:'head'},e('tr',{},['CYP','Role','Prediction','Probability','Experimental','Domain','Confidence','Model',''].map(label=>e('th',{key:label},label)))),
   e('tbody',{key:'body'},rows.map(prediction=>{
    const output=prediction.outputs||{},evidence=output.experimental_evidence||[];
    const experimental=evidence.length?evidence.map(item=>item.value+' '+item.unit+' ('+item.comparison+')').join(' · '):'—';
    const liability=output.liability_summary?.flag;
    return e('tr',{key:prediction.id},[
     e('td',{key:'isoform'},output.isoform||prediction.endpoint.split(' ')[0]),
     e('td',{key:'role'},output.role||prediction.endpoint.split(' ')[1]?.toUpperCase()),
     e('td',{key:'class'},[output.classification||'—',liability&&e('div',{key:'flag',className:'fail small'},liability)]),
     e('td',{key:'probability',className:'mono'},Number(output.probability??prediction.predicted_value).toFixed(4)),
     e('td',{key:'experimental'},experimental),
     e('td',{key:'domain'},prediction.applicability_domain),
     e('td',{key:'confidence'},prediction.confidence),
     e('td',{key:'model',className:'small'},prediction.model?.model_name+' '+prediction.model?.model_version),
     e('td',{key:'details'},predictionDetails(prediction))
    ]);
   }))
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
   e('div',{className:'grid',key:'summary'},[group('Strengths',summary.strengths,'strengths'),group('Concerns',summary.concerns,'concerns'),group('Not available',summary.unknown,'')]),
   e('div',{key:'audit',className:profile.provenance_audit?.status==='PASS'?'pass':'fail'},'Provenance audit: '+profile.provenance_audit?.status+' · '+profile.provenance_audit?.checked+' latest endpoint predictions checked')
  ]);
 }

 function admetPredictionTable(rows){
  if(!rows.length)return e('div',{className:'empty-state'},[StatusBadge({type:'Not predicted'}),e('p',{key:'text'},'Prediction not run for this CompoundVersion.'),e('button',{key:'run',className:'secondary',disabled:admetBusy||!detail?.version,onClick:()=>runPrediction(detail.version.id)},'Run Predictions')]);
  return e('table',{},[
   e('thead',{key:'head'},e('tr',{},['Compound','Endpoint','Experimental','Predicted','Confidence','Domain',''].map(label=>e('th',{key:label},label)))),
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
     e('td',{key:'confidence'},prediction.confidence),e('td',{key:'domain'},prediction.applicability_domain),
     e('td',{key:'details'},predictionDetails(prediction))
    ]);
   }))
  ]);
 }

 function ExperimentalDataPanel(){
  if(!detail?.version)return e('div',{className:'empty-state'},[StatusBadge({type:'Not applicable'}),e('p',{},'Add a valid structure before entering structure-linked experimental data.')]);
  const select=(label,value,options,onChange)=>e('div',{},[e('label',{key:'label'},label),e('select',{key:'select',value,onChange:event=>onChange(event.target.value)},options.map(option=>{const item=typeof option==='object'?option:{value:option,label:option};return e('option',{key:item.value,value:item.value},item.label)}))]);
  const input=(name,row,key,label,type='text')=>Field({label,type,value:row[key],onChange:value=>setExperimentValue(name,key,value)});
  const formFor=name=>{
   const row=experimentalDrafts[name]||experimentDefaults(name),fields=[];
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
   ...['Activity','ADME','Safety'].map(category=>e('div',{key:category},[e('h4',{key:'title'},category),e('div',{className:'endpoint-selector',key:'options'},EXPERIMENT_OPTIONS.filter(row=>row[1]===category).map(([name])=>e('label',{key:name,className:'check-option'},[e('input',{type:'checkbox',checked:experimentalSelected.includes(name),onChange:()=>toggleExperiment(name)}),e('span',{},name)]))) ])),
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
     e('div',{className:'col-6 structure',key:'svg'},Svg({src:run.highlighted_svg})),
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
      e('td',{key:'rank'},item.rank),e('td',{key:'smiles',className:'mono small'},item.canonical_smiles),e('td',{key:'transform'},item.transformation),e('td',{key:'atom'},item.source_atom),e('td',{key:'phase'},item.phase),e('td',{key:'confidence'},item.confidence),
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
  const setConstraint=(key,value)=>setOptimizationForm(current=>({...current,constraints:{...current.constraints,[key]:value}}));
  const toggleObjective=name=>setOptimizationForm(current=>({...current,objectives:current.objectives.includes(name)?current.objectives.filter(value=>value!==name):[...current.objectives,name]}));
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
  const constraintField=(key,label,type='number')=>e('div',{className:'col-3',key},Field({label,type,value:optimizationForm.constraints[key],onChange:value=>setConstraint(key,value)}));
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
    e('h3',{key:'title'},'Optimization Run'),
    e('p',{key:'parent',className:'small'},'Parent: '+detail.compound_id+' v'+detail.current_version+' · CompoundVersion #'+versionId),
    e('div',{className:'grid',key:'top'},[
     e('div',{className:'col-4',key:'assay'},[e('label',{},'Selected assay'),e('select',{value:optimizationForm.assay_id,onChange:event=>setOptimizationForm(current=>({...current,assay_id:event.target.value}))},[e('option',{key:'none',value:''},'No assay selected'),...assays.map(assay=>e('option',{key:assay.id,value:assay.id},assay.name+' · '+assay.measurement_type))])]),
     e('div',{className:'col-8',key:'objectives'},[e('label',{},'Optimization objective(s)'),e('div',{className:'objective-grid'},config.objectives.map(name=>e('label',{key:name,className:'check-option'},[e('input',{key:'input',type:'checkbox',checked:optimizationForm.objectives.includes(name),onChange:()=>toggleObjective(name)}),e('span',{key:'label'},name)])))])
    ]),
    optimizationForm.objectives.includes('Custom')&&e('div',{key:'custom',style:{marginTop:'10px'}},Field({label:'Custom objective',value:optimizationForm.custom_objective,onChange:value=>setOptimizationForm(current=>({...current,custom_objective:value}))})),
    e('h4',{key:'constraints-title',style:{marginTop:'18px'}},'Constraints'),
    e('div',{className:'grid',key:'constraints'},[
     constraintField('potency_max_nm','Potency IC50 ≤ (nM)'),constraintField('do_not_worsen_fold','Do not worsen potency > fold'),constraintField('clogp_max','cLogP ≤'),constraintField('mw_max','MW ≤'),
     constraintField('tpsa_min','TPSA minimum Å²'),constraintField('tpsa_max','TPSA maximum Å²'),constraintField('similarity_min','Future analog similarity ≥'),constraintField('logs_min','LogS minimum'),constraintField('caco2_logpapp_min','Caco-2 LogPapp minimum'),
     e('div',{className:'col-3',key:'herg'},e('label',{className:'check-option'},[e('input',{type:'checkbox',checked:!!optimizationForm.constraints.herg_do_not_increase,onChange:event=>setConstraint('herg_do_not_increase',event.target.checked)}),e('span',{},'hERG: do not increase liability')]))
    ]),
    e('p',{key:'precedence',className:'small'},'Experimental evidence takes precedence over prediction. Low-confidence classification alone remains supporting-only. Similarity and do-not-worsen constraints are stored now as hard gates for a future proposal stage; Stage 4A does not create candidates.'),
    e('button',{key:'analyze',disabled:optimizationBusy||!optimizationForm.objectives.length,onClick:()=>analyzeOptimization(versionId)},optimizationBusy?'Analyzing…':'Analyze strategy')
   ]),
   optimizationRuns.length>0&&e('div',{className:'row run-picker',key:'history'},[e('label',{key:'label'},'Saved runs'),e('select',{key:'select',value:run?.id||'',onChange:event=>setOptimizationRun(optimizationRuns.find(item=>item.id===Number(event.target.value)))},optimizationRuns.map(item=>e('option',{key:item.id,value:item.id},'#'+item.id+' · '+item.objectives.join(' + ')+' · '+item.status))) ]),
   run&&e(React.Fragment,{key:'results'},[
    e('div',{className:'card',key:'profile'},[
     e('div',{className:'row toolbar',key:'head'},[e('h3',{},'Current profile'),e('span',{className:'small'},run.engine+' '+run.engine_version)]),
     e('div',{className:'grid',key:'profile-grid'},[
      e('div',{className:'col-4',key:'activity'},[e('h4',{},'Activity'),e('p',{className:'small'},activity.experimental?'Experimental '+activity.experimental.mean_nm+' nM':(activity.predicted?'Predicted '+activity.predicted.value_nm+' nM · '+activity.predicted.confidence:'Not measured / not predicted for selected assay'))]),
      e('div',{className:'col-4',key:'properties'},[e('h4',{},'Properties'),e('p',{className:'small'},['molecular_weight','clogp','tpsa','fraction_csp3'].map(key=>key+' '+(properties[key]?.value??'—')).join(' · ')+' · Calculated / RDKit')]),
      e('div',{className:'col-4',key:'admet'},[e('h4',{},'ADMET / metabolism'),...Object.entries(admetProfile).slice(0,12).map(([name,row])=>e('div',{key:name,className:'small'},name+': '+evidenceValue(row))),Object.keys(admetProfile).length===0&&Empty({children:'No compatible experimental or predicted ADMET evidence.'})])
     ]),
     e('details',{key:'hierarchy'},[e('summary',{key:'summary'},'Evidence hierarchy'),e('ol',{key:'list',className:'small'},(run.evidence.evidence_hierarchy||[]).map(item=>e('li',{key:item.rank},item.type+' · ordinal weight '+item.weight)))])
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

 function compoundDetail(){
  if(!detail)return null;
  const version=detail.version,detailMeasurements=admet?.measurements||[],detailRuns=admet?.prediction_runs||[],detailPredictions=admet?.predictions||[];
  const activity=workspace?.activity||{measurements:[],predictions:[]},properties=version?.properties||{};
  const tabs=['overview','properties','activity','admet','metabolism','optimization','history'];
  const highlights=detailPredictions.filter((row,index,array)=>array.findIndex(item=>item.endpoint===row.endpoint)===index).slice(0,5);
  const activityTable=e('div',{},[
   e('h3',{key:'exp'},'Experimental Activity'),activity.measurements.length?e('table',{key:'exp-table'},[e('thead',{},e('tr',{},['Assay','Measurement','Value','Source'].map(x=>e('th',{key:x},x)))),e('tbody',{},activity.measurements.map(row=>e('tr',{key:row.id},[e('td',{},row.assay),e('td',{},row.measurement_type),e('td',{className:'mono'},row.qualifier+' '+row.value+' '+row.unit),e('td',{},[StatusBadge({type:'Experimental'}),' '+row.source])])))]):e('div',{className:'empty-state'},[StatusBadge({type:'Not measured'}),e('p',{},'No experimental activity measurement entered.')]),
   e('h3',{key:'pred',style:{marginTop:'22px'}},'Activity Prediction'),activity.predictions.length?e('table',{key:'pred-table'},[e('thead',{},e('tr',{},['Assay','Predicted value','Confidence','Domain'].map(x=>e('th',{key:x},x)))),e('tbody',{},activity.predictions.map(row=>e('tr',{key:row.id},[e('td',{},row.assay),e('td',{className:'mono'},row.predicted_value_nm+' nM'),e('td',{},row.confidence),e('td',{},row.applicability_domain)])))]):e('div',{className:'empty-state'},[StatusBadge({type:'Not predicted'}),e('p',{},'No activity prediction run for this CompoundVersion.'),e('a',{className:'button secondary',href:'/static/stage2-workbench.html?project='+projectId},'Open Activity Workbench')])
  ]);
  return e('div',{className:'compound-workspace'},[
   e('div',{className:'card compound-hero',key:'hero'},[e('div',{className:'compound-hero-structure'},Svg({src:version?.highlighted_svg||version?.svg})),e('div',{className:'compound-hero-copy'},[e('div',{className:'eyebrow'},'COMPOUND DETAIL'),e('h2',{},detail.name),e('div',{className:'row'},[StatusBadge({type:detail.status}),e('span',{className:'mono'},detail.compound_id+(version?' · Version '+version.version_number:' · No structure version'))]),e('p',{className:'small'},workspace?'Strict scope: Project #'+workspace.scope.project_id+' · Compound #'+workspace.scope.compound_id+' · CompoundVersion #'+workspace.scope.version_id:'Draft compound; no version-linked data exists.'),e('div',{className:'row'},[version&&e('button',{className:'secondary',onClick:updateStructure},'Modify Structure / New Version'),e('button',{className:'secondary',onClick:()=>setDetail(null)},'Back to Compounds')])])]),
   e('nav',{className:'detail-tabs',key:'tabs'},tabs.map(tab=>e('button',{key:tab,className:detailTab===tab?'':'secondary',disabled:!version&&['properties','activity','admet','metabolism','optimization'].includes(tab),onClick:()=>setDetailTab(tab)},tab.toUpperCase()))),
   detailTab==='overview'&&e('div',{className:'grid',key:'overview'},[
    e('div',{className:'card col-4'},[
     e('h3',{},'Key Properties'),
     ...['clogp','tpsa','qed'].map(key=>e('div',{key,className:'metric-row'},[
      e('span',{},key==='clogp'?'cLogP':key.toUpperCase()),
      properties[key]!=null?e('strong',{className:'mono'},properties[key]):StatusBadge({type:'Not calculated'})
     ])),
     !version?.calculated&&e('button',{onClick:calculateProperties},'Calculate Properties')
    ]),
    e('div',{className:'card col-4'},[e('h3',{},'Activity'),activity.measurements[0]?e('p',{},[StatusBadge({type:'Experimental'}),' ',activity.measurements[0].value+' '+activity.measurements[0].unit+' · '+activity.measurements[0].assay]):activity.predictions[0]?e('p',{},[StatusBadge({type:'Predicted'}),' '+activity.predictions[0].predicted_value_nm+' nM · '+activity.predictions[0].assay]):e('p',{},[StatusBadge({type:'Not measured'}),' No activity data'])]),
    e('div',{className:'card col-12'},[e('h3',{},'ADMET Highlights'),highlights.length?e('div',{className:'highlight-grid'},highlights.map(row=>e('div',{key:row.endpoint,className:'highlight-item'},[e('strong',{},row.endpoint==='Permeability'?'Caco-2 Permeability':row.endpoint),e('div',{className:'mono'},row.predicted_value+' '+row.unit),StatusBadge({type:'Predicted'})]))):e('div',{className:'empty-state'},[StatusBadge({type:'Not predicted'}),e('p',{},'No ADMET predictions run for this CompoundVersion.')])])
   ]),
   detailTab==='properties'&&e('div',{className:'card',key:'properties'},
    version?.calculated?e('div',{className:'grid'},[
     e('div',{className:'col-6'},propertyTable(properties)),
     e('div',{className:'col-6'},[ruleTable(version.rules),e('h4',{},'Structural Alerts'),...alertList(version.alerts),e('button',{className:'secondary',onClick:calculateProperties},'Recalculate Properties')])
    ]):e('div',{className:'empty-state'},[StatusBadge({type:'Not calculated'}),e('h3',{},'Properties have not been calculated.'),e('button',{onClick:calculateProperties},'Calculate Properties')])
   ),
   detailTab==='activity'&&e('div',{className:'card',key:'activity'},activityTable),
   detailTab==='admet'&&e('div',{key:'admet'},[
    e('div',{className:'card row toolbar'},[
     e('div',{},[e('h3',{},'ADMET'),e('p',{className:'small'},'Only '+detail.name+' Version '+version.version_number+' records are loaded.')]),
     e('div',{className:'row'},[e('button',{className:'secondary',onClick:()=>setExperimentalOpen(!experimentalOpen)},'Add Experimental Data'),e('button',{disabled:admetBusy,onClick:()=>runPrediction(version.id)},admetBusy?'Predicting…':'Run Predictions')])
    ]),
    experimentalOpen&&e('div',{className:'card'},ExperimentalDataPanel()),integratedProfile(version.id),
    e('div',{className:'card'},[e('h3',{},'Experimental'),detailMeasurements.length?admetMeasurementTable(detailMeasurements):e('div',{className:'empty-state'},[StatusBadge({type:'Not measured'}),e('p',{},'No experimental measurement entered.'),e('button',{className:'secondary',onClick:()=>setExperimentalOpen(true)},'Add Experimental Data')])]),
    e('div',{className:'card'},[e('h3',{},'Prediction'),e('h4',{},'Absorption · Aqueous Solubility and Caco-2 Permeability'),admetPredictionTable(detailPredictions.filter(row=>['Solubility','Permeability'].includes(row.endpoint))),e('h4',{},'Distribution · Plasma Protein Binding (PPB) and fu'),admetPredictionTable(detailPredictions.filter(row=>row.endpoint==='Plasma protein binding')),e('h4',{},'Microsomal Stability (MS)'),admetPredictionTable(detailPredictions.filter(row=>row.endpoint.endsWith('intrinsic clearance'))),e('h4',{},'CYP Inhibition / Substrate'),cypPredictionTable(detailPredictions.filter(row=>row.endpoint.startsWith('CYP'))),e('h4',{},'Transporters'),transporterPredictionTable(detailPredictions.filter(row=>TRANSPORTER_ENDPOINTS.has(row.endpoint))),unavailableTransporterModels(),e('h4',{},'Safety'),safetyPredictionTable(detailPredictions.filter(row=>SAFETY_ENDPOINTS.has(row.endpoint))),unavailableSafetyModels()])
   ]),
   detailTab==='metabolism'&&e('div',{className:'card',key:'metabolism'},[e('h3',{},'Metabolic Soft Spots and Metabolite Hypotheses'),metabolismPanel(version.id)]),
   detailTab==='optimization'&&(project?.molecule_type==='Small Molecule'?optimizationPanel(version.id):e('div',{className:'card empty-state'},[StatusBadge({type:'Not applicable'}),e('p',{},'This model currently supports small molecules only.')])),
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
   e('div',{className:'row toolbar',key:'header'},[e('div',{},[e('div',{className:'eyebrow'},'NEW COMPOUND'),e('h2',{},'Add Compound')]),e('button',{className:'secondary',onClick:()=>setAddCompoundOpen(false)},'Close')]),
   e('div',{className:'grid',key:'identity'},[e('div',{className:'col-6'},Field({label:'Compound Name *',value:compoundForm.name,onChange:value=>setCompoundForm(current=>({...current,name:value})),placeholder:'HIT-001'})),e('div',{className:'col-6'},Field({label:'Compound ID (optional)',value:compoundForm.compound_id,onChange:value=>setCompoundForm(current=>({...current,compound_id:value})),placeholder:'Generated from name if empty'}))]),
   smallMolecule?e(React.Fragment,{key:'editor'},[
    e('h3',{key:'title',style:{marginTop:'22px'}},'Draw Chemical Structure'),e('p',{key:'help',className:'small'},'Draw or edit the compound structure below. The SMILES field updates automatically while you draw.'),
    e('div',{className:'structure-editor-shell',key:'shell'},[
     e('iframe',{key:'frame',id:'ketcher-editor',className:'ketcher-frame'+(editorReady?'':' loading'),title:'Ketcher Chemical Structure Editor',src:'/static/ketcher/standalone/index.html'}),
     !editorReady&&e('div',{className:'structure-editor-loading',key:'loading'},[e('strong',{},'Structure Editor is loading…'),e('span',{},'Drawing tools and the linked SMILES field will appear here.')])
    ]),
    e('h3',{key:'smiles-title',style:{marginTop:'20px'}},'Or Enter SMILES'),e('div',{className:'row',key:'smiles'},[e('div',{style:{flex:1}},Field({label:'SMILES',value:compoundForm.smiles,onChange:value=>{setCompoundForm(current=>({...current,smiles:value}));setPreview(null)},placeholder:'Paste SMILES or draw above'})),e('button',{className:'secondary',disabled:!compoundForm.smiles,onClick:loadSmilesIntoEditor},'Load in Editor'),e('button',{className:'secondary',disabled:!compoundForm.smiles,onClick:validate},'Validate Structure')])
   ]):e('div',{className:'empty-state',key:'peptide'},[StatusBadge({type:'Not applicable'}),e('h3',{},'Peptide project'),e('p',{},'This model currently supports small molecules only. Save the compound as a draft; peptide-specific calculations are not run.')]),
   e('div',{style:{marginTop:'16px'},key:'notes'},Field({label:'Description / Notes',value:compoundForm.notes,onChange:value=>setCompoundForm(current=>({...current,notes:value})),type:'textarea'})),
   preview&&e('div',{className:'structure-validation',key:'preview'},[StatusBadge({type:'Calculated'}),e('span',{},' Valid structure · '+preview.identity.canonical_smiles)]),
   e('div',{className:'row modal-actions',key:'actions'},[e('button',{className:'secondary',disabled:!compoundForm.name.trim(),onClick:()=>saveCompound(false)},'Save Compound'),e('button',{disabled:!compoundForm.name.trim()||!compoundForm.smiles.trim()||!smallMolecule,onClick:()=>saveCompound(true)},'Save & Calculate'),e('span',{className:'small'},'Save Compound does not require property calculation.')])
  ]));
 }

 function ComparePanel(){
  const groups={Properties:['MW','cLogP','TPSA','QED'],Activity:['Activity'],ADME:['Solubility','Caco-2','PPB','HLM','RLM'],Safety:['hERG','Ames','DILI']};
  const metrics=(comparison?.metrics||[]).filter(metric=>compareMetrics.includes(metric));
  return e('div',{},[
   e('div',{className:'card',key:'config'},[
    e('div',{className:'row toolbar'},[
     e('div',{},[e('div',{className:'eyebrow'},'COMPARE COMPOUNDS'),e('h2',{},'Comparison Configuration'),e('p',{className:'small'},selected.length+' compounds selected. Multi-compound data appears only here.')]),
     e('button',{disabled:selected.length<2,onClick:compare},'Refresh Comparison')
    ]),
    e('div',{className:'grid'},[
     ...Object.entries(groups).map(([group,items])=>e('div',{className:'col-3',key:group},[
      e('h4',{},group),
      ...items.map(metric=>e('label',{key:metric,className:'check-option'},[
       e('input',{type:'checkbox',checked:compareMetrics.includes(metric),onChange:event=>setCompareMetrics(current=>event.target.checked?[...current,metric]:current.filter(item=>item!==metric))}),e('span',{},metric)
      ]))
     ])),
     e('div',{className:'col-4'},[e('label',{},'Activity assay'),e('select',{value:compareAssay,onChange:event=>setCompareAssay(event.target.value)},[e('option',{value:''},'Latest experimental assay'),...assays.map(row=>e('option',{key:row.id,value:row.id},row.name))])])
    ])
   ]),
   comparison&&e('div',{className:'card',key:'table'},[
    e('h3',{},'Selected Compound Comparison'),e('p',{className:'small'},'Experimental values take precedence. Each cell retains its evidence type. No overall score or automatic ranking is calculated.'),
    e('table',{},[
     e('thead',{},e('tr',{},['Compound',...metrics].map(label=>e('th',{key:label},label)))),
     e('tbody',{},comparison.compounds.map(compound=>e('tr',{key:compound.row_id},[
      e('td',{},[e('strong',{},compound.name||compound.compound),e('div',{className:'mono small'},compound.compound)]),
      ...metrics.map(metric=>e('td',{key:metric,className:'mono'},[compound[metric]??'—',e('div',{key:'source'},StatusBadge({type:compound.sources?.[metric]||'Not measured'}))]))
     ])))
    ]),
    e('div',{className:'plot-grid'},pairPlot(comparison,'MW','cLogP'),pairPlot(comparison,'TPSA','cLogP'),qedHistogram(comparison))
   ])
  ]);
 }

 function SettingsPanel(){
  const sample=currentVersions.find(row=>row.version)?.compound_id||'C001',placeholder='compound_id,version_number,endpoint,species,matrix,value,unit,qualifier,replicate,mean,sd,n,method,source,date,notes\n'+sample+',1,Solubility,Human,,12.5,µM,=,R1,,,,shake flask,Study A,2026-08-25,';
  return e('div',{},[
   e('div',{className:'card',key:'identity'},[e('h2',{},'Project Settings'),e('p',{className:'small'},'Indication and mechanism remain preserved in the database but are kept out of the primary creation workflow.'),e('div',{className:'grid'},[e('div',{className:'col-4'},Field({label:'Project Name',value:project.name,onChange:value=>setProject(current=>({...current,name:value}))})),e('div',{className:'col-4'},Field({label:'Target',value:project.target,onChange:value=>setProject(current=>({...current,target:value}))})),e('div',{className:'col-4'},[e('label',{},'Molecule Type'),e('select',{value:project.molecule_type,onChange:event=>setProject(current=>({...current,molecule_type:event.target.value}))},['Small Molecule','Peptide'].map(value=>e('option',{key:value,value},value)))]),e('div',{className:'col-12'},Field({label:'Description',value:project.description||'',onChange:value=>setProject(current=>({...current,description:value})),type:'textarea'}))]),e('button',{disabled:!project.name.trim()||!project.target.trim(),onClick:saveProjectSettings},'Save Project Settings')]),
   e('div',{className:'card',key:'csv'},[e('h2',{},'Experimental ADMET CSV'),e('p',{className:'small'},'Advanced project-wide import/export. Compound Detail remains CompoundVersion-isolated.'),e('a',{className:'button secondary',href:'/api/projects/'+projectId+'/admet/export.csv'},'Export CSV'),e('textarea',{rows:7,value:admetCsv,placeholder,onChange:event=>{setAdmetCsv(event.target.value);setAdmetCsvPreview(null)}}),e('div',{className:'row',style:{marginTop:'10px'}},[e('button',{className:'secondary',disabled:admetBusy||!admetCsv.trim(),onClick:previewAdmet},'Preview CSV'),e('button',{disabled:admetBusy||!admetCsvPreview||admetCsvPreview.errors.length>0||!admetCsvPreview.valid_count,onClick:importAdmet},'Import Valid Rows')]),admetCsvPreview&&e('p',{className:admetCsvPreview.errors.length?'fail':'pass'},admetCsvPreview.valid_count+' valid · '+admetCsvPreview.errors.length+' errors')])
  ]);
 }

 const goDashboard=()=>{setProjectTab('dashboard');setDetail(null);setAddCompoundOpen(false);setComparison(null);setSelectedCandidate(null);setSidebarOpen(false);loadDashboard().catch(error=>setMessage(String(error)))};
 const openProject=itemId=>{setProjectId(itemId);setProjectTab('compounds');setDetail(null);setAddCompoundOpen(false);setComparison(null);setSidebarOpen(false)};
 const navigateScientific=(projectTarget='compounds',detailTarget='')=>{
  if(!projectId){goDashboard();setMessage('Create or select a project to use this module.');return}
  setProjectTab(projectTarget);setSidebarOpen(false);
  if(detail&&detailTarget)setDetailTab(detailTarget);else if(projectTarget!=='compounds')setDetail(null);
 };
 const selectedProjectSummary=(dashboard?.projects||[]).find(row=>row.id===projectId);

 function MainDashboard(){
  const registry=dashboard?.model_registry||[];
  const registryStatus=(endpoint,activeStatus='LIMITED')=>{
   const model=registry.find(row=>row.endpoint===endpoint);
   return model?.active?activeStatus:'MODEL_UNAVAILABLE';
  };
  const modules=[
   {title:'Structure & Chemistry',status:'READY',description:'Version-controlled chemical identity and calculated molecular properties.',items:[['Structure Drawing','READY'],['SMILES Input','READY'],['Compound Versioning','READY'],['Structure Validation','READY'],['Physicochemical Properties','READY'],['Structural Alerts','READY']]},
   {title:'Activity & SAR',status:'READY',description:'Project-local experimental activity, modeling, and structure–activity evidence.',items:[['Experimental Activity','READY'],['IC50 / EC50 / Ki / Kd / GI50','READY'],['Project QSAR','LIMITED'],['Similarity Analysis','READY'],['SAR / MMP / Activity Cliff','READY']]},
   {title:'ADME',status:'READY',description:'Absorption, distribution, and metabolic-liability evidence.',items:[['Solubility',registryStatus('Solubility','READY')],['Caco-2 Permeability',registryStatus('Permeability','READY')],['Plasma Protein Binding',registryStatus('Plasma protein binding','READY')],['fu','READY'],['HLM / RLM / MLM',registryStatus('HLM intrinsic clearance')],['Metabolic Soft Spot','LIMITED'],['Metabolite Hypothesis','LIMITED']]},
   {title:'CYP & Transporters',status:'LIMITED',description:'Endpoint- and role-separated metabolism and transporter classifications.',items:[['CYP1A2 inhibitor',registryStatus('CYP1A2 inhibitor')],['CYP2C9 inhibitor',registryStatus('CYP2C9 inhibitor')],['CYP2C19 inhibitor',registryStatus('CYP2C19 inhibitor')],['CYP2D6 inhibitor',registryStatus('CYP2D6 inhibitor')],['CYP3A4 inhibitor',registryStatus('CYP3A4 inhibitor')],['CYP substrates where available','LIMITED'],['P-gp inhibitor',registryStatus('P-gp inhibitor')]],unavailable:registry.filter(row=>/P-gp substrate|BCRP|BSEP|OATP|OCT|MATE/.test(row.endpoint)&&!row.active).map(row=>row.endpoint)},
   {title:'Safety / Toxicology',status:'LIMITED',description:'Classification evidence and calculated structural safety alerts.',items:[['hERG',registryStatus('hERG liability')],['Ames',registryStatus('Ames mutagenicity')],['DILI',registryStatus('DILI clinical liability')],['Structural Alerts','READY']]},
   {title:'Optimization',status:'READY',description:'Deterministic strategy, analog generation, filtering, and transparent ranking.',items:[['Liability Analysis','READY'],['Protected / Modifiable Regions','READY'],['Medicinal Chemistry Transformations','READY'],['Analog Generation','READY'],['Re-scoring','LIMITED'],['Pareto Optimization','READY'],['Top Candidate Selection','READY']]},
   {title:'PK / DMPK',status:'PLANNED',description:'Stage 5 capabilities are intentionally not active.',items:[['Experimental PK','PLANNED'],['NCA','PLANNED'],['IVIVE','PLANNED'],['PK Simulation','PLANNED']]}
  ];
  const summaries=dashboard?.projects||projects;
  return e(React.Fragment,{},[
   e('section',{className:'card dashboard-hero',key:'intro'},[
    e('div',{className:'eyebrow'},'PLATFORM OVERVIEW'),e('h1',{},'Drug Optimization Platform'),
    e('p',{},'Structure, activity, ADMET and medicinal chemistry optimization data are integrated at the compound-version level to support hit-to-lead and lead optimization decisions.'),
    e('ul',{className:'dashboard-capabilities'},['Structure-based compound management','Experimental data integration','Predictive ADMET','SAR / optimization workflow','Full prediction provenance'].map(item=>e('li',{key:item},item))),
    e('div',{className:'dashboard-stats'},[
     e('div',{className:'dashboard-stat',key:'projects'},[e('span',{},'Projects'),e('strong',{},String(dashboard?.totals?.projects??projects.length))]),
     e('div',{className:'dashboard-stat',key:'compounds'},[e('span',{},'Compounds'),e('strong',{},String(dashboard?.totals?.compounds??projects.reduce((sum,row)=>sum+(row.compound_count||0),0)))]),
     e('div',{className:'dashboard-stat',key:'scope'},[e('span',{},'Default data scope'),e('strong',{className:'dashboard-stat-text'},'CompoundVersion'),e('small',{},'Project-isolated')])
    ])
   ]),
   e('section',{className:'dashboard-section',key:'modules'},[
    e('div',{className:'section-heading'},[e('div',{},[e('div',{className:'eyebrow'},'SCIENTIFIC WORKSPACE'),e('h2',{},'Available Scientific Modules')]),e('p',{className:'small'},'Status reflects the current local engine and model registry.')]),
    e('div',{className:'module-grid'},modules.map(module=>e('article',{className:'module-card',key:module.title},[
     e('div',{className:'module-card-head'},[e('h3',{},module.title),StatusBadge({type:module.status})]),e('p',{className:'small'},module.description),
     e('ul',{className:'module-list'},module.items.map(([label,status])=>e('li',{key:label},[e('span',{},label),StatusBadge({type:status})]))),
     module.unavailable?.length>0&&e('div',{className:'module-unavailable'},['Unavailable: ',module.unavailable.join(' · ')])
    ])))
   ]),
   e('div',{className:'dashboard-split',key:'start'},[
    e('section',{className:'card dashboard-create',key:'create'},[e('div',{className:'eyebrow'},'NEW WORKSPACE'),e('h2',{},'Create New Project'),e('div',{className:'create-project-grid'},[
     e(Field,{label:'Project Name *',value:form.name,onChange:value=>setForm({...form,name:value}),placeholder:'EGFR Exon20ins'}),e(Field,{label:'Target *',value:form.target,onChange:value=>setForm({...form,target:value}),placeholder:'EGFR'}),e('div',{},[e('label',{},'Molecule Type'),e('select',{value:form.molecule_type,onChange:event=>setForm({...form,molecule_type:event.target.value})},['Small Molecule','Peptide'].map(value=>e('option',{key:value,value},value)))])
    ]),e('button',{disabled:!form.name.trim()||!form.target.trim(),onClick:createProject},'Create Project'),e('p',{className:'small dashboard-note'},'Description and additional metadata can be added later in Project Settings.')]),
    e('section',{className:'card quick-start',key:'quick'},[e('div',{className:'eyebrow'},'QUICK START'),e('h2',{},'Typical Workflow'),e('ol',{},['Create Project','Add Compound','Draw Structure','Calculate Properties','Add Experimental Data','Run Predictions','Compare / Optimize'].map(item=>e('li',{key:item},item))),e('p',{className:'small'},'Save and calculation remain separate. Prediction and experimental evidence are never merged.')])
   ]),
   e('section',{className:'card',key:'defaults'},[e('h2',{},'Default Workspace Settings'),e('div',{className:'dashboard-settings'},[
    e('div',{className:'dashboard-setting',key:'type'},[e('span',{},'Default molecule type'),e('strong',{},'Small Molecule')]),e('div',{className:'dashboard-setting',key:'entry'},[e('span',{},'Structure entry'),e('strong',{},'Ketcher or SMILES')]),e('div',{className:'dashboard-setting',key:'calc'},[e('span',{},'Calculation policy'),e('strong',{},'Save first · Calculate on demand')]),e('div',{className:'dashboard-setting',key:'isolation'},[e('span',{},'Data isolation'),e('strong',{},'Project + CompoundVersion')])
   ])]),
   e('section',{className:'card',key:'projects'},[
    e('div',{className:'row toolbar'},[e('div',{},[e('div',{className:'eyebrow'},'RESEARCH PORTFOLIO'),e('h2',{},'Projects'),e('p',{className:'small'},'Project cards summarize recorded evidence without synthetic progress percentages.')]),projectId&&e('button',{className:'secondary',onClick:()=>openProject(projectId)},'Continue Current Project')]),
    summaries.length?e('div',{className:'dashboard-project-grid'},summaries.map(item=>e('article',{className:'dashboard-project',key:item.id,tabIndex:0,onClick:()=>openProject(item.id),onKeyDown:event=>{if(event.key==='Enter'||event.key===' ')openProject(item.id)}},[
     e('div',{className:'dashboard-project-head'},[e('div',{},[e('div',{className:'eyebrow'},item.molecule_type||'Small Molecule'),e('h3',{},item.name)]),e('span',{className:'dashboard-count'},item.compound_count||0)]),
     e('dl',{},[e('div',{key:'target'},[e('dt',{},'Target'),e('dd',{},item.target||'Not set')]),e('div',{key:'experimental'},[e('dt',{},'Experimental records'),e('dd',{},String((item.experimental_activity_count||0)+(item.experimental_admet_count||0)))]),e('div',{key:'optimization'},[e('dt',{},'Optimization runs'),e('dd',{},String(item.optimization_run_count||0))])]),
     e('p',{className:'project-status-summary'},item.status_summary||((item.compound_count||0)+' compounds · experimental and prediction data not started')),e('span',{className:'project-open-link'},'Open Project →')
    ]))):e('div',{className:'empty-state'},[e('h3',{},'No projects yet'),e('p',{},'Use Create New Project above to begin a compound-version-isolated workspace.')])
   ])
  ]);
 }

 function ProjectWorkspace(){
  const summary=selectedProjectSummary;
  const statusByCompound=new Map((summary?.compounds||[]).map(row=>[row.row_id,row]));
  return e(React.Fragment,{},[
   e('div',{className:'card project-header',key:'header'},project?e('div',{className:'row toolbar'},[e('div',{},[e('div',{className:'eyebrow'},'PROJECT DASHBOARD'),e('h1',{},project.name),e('div',{},[e('strong',{},project.target||'Target not set'),' · ',project.molecule_type])]),e('button',{onClick:()=>{setAddCompoundOpen(true);setCompoundForm({compound_id:'',name:'',smiles:'',notes:''})}},'Add Compound')]):e('div',{},[e('h2',{},'Select or create a project'),e('p',{},'Start with a project, then add compounds and work from Compound Detail.') ])),
   project&&e('nav',{className:'project-nav',key:'nav'},[['compounds','Compounds'],['assays','Assays'],['compare','Compare'],['optimization','Optimization Runs'],['settings','Settings']].map(([tab,label])=>e('button',{key:tab,className:projectTab===tab?'':'secondary',onClick:()=>{setProjectTab(tab);if(tab!=='compounds')setDetail(null)}},label))),
   project&&projectTab==='compounds'&&!detail&&e(React.Fragment,{key:'project-dashboard'},[
    e('section',{className:'card',key:'overview'},[e('div',{className:'eyebrow'},'PROJECT OVERVIEW'),e('h2',{},'Current Project Status'),e('div',{className:'project-overview-grid'},[
     ['Target',project.target||'Not set'],['Molecule Type',project.molecule_type],['Compounds',summary?.compound_count??currentVersions.length],['Experimental Activity',summary?.experimental_activity_count??0],['Experimental ADMET',summary?.experimental_admet_count??0],['Predictions',summary?.prediction_count??0],['Optimization Runs',summary?.optimization_run_count??0]
    ].map(([label,value])=>e('div',{className:'project-overview-item',key:label},[e('span',{},label),e('strong',{},String(value))])))]),
    e('section',{className:'card workflow-card',key:'workflow'},[e('div',{className:'eyebrow'},'WORKFLOW STATUS'),e('div',{className:'workflow-strip'},['Structure','Properties','Activity','ADMET','Optimization','PK'].map((stage,index)=>e(React.Fragment,{key:stage},[e('div',{className:'workflow-step'},[e('span',{},stage),StatusBadge({type:summary?.workflow?.[stage]||(stage==='PK'?'PLANNED':'NOT_STARTED')})]),index<5&&e('span',{className:'workflow-arrow'},'→')])))]),
    e('section',{className:'card',key:'compounds'},[e('div',{className:'row toolbar'},[e('div',{},[e('h2',{},'Compound Status'),e('p',{className:'small'},'Each row summarizes only the current CompoundVersion in this project.')]),e('div',{className:'row'},[e('button',{className:'secondary',disabled:selected.length<2,onClick:compare},'Compare Selected'),e('button',{onClick:()=>setAddCompoundOpen(true)},'Add Compound')])]),currentVersions.length?e('div',{className:'table-scroll'},e('table',{className:'compound-list project-status-table'},[e('thead',{},e('tr',{},['','Compound','Structure','Properties','Activity','ADMET','Optimization',''].map((x,index)=>e('th',{key:x||index},x)))),e('tbody',{},currentVersions.map(compound=>{const status=statusByCompound.get(compound.row_id)||{};return e('tr',{key:compound.row_id},[e('td',{},e('input',{type:'checkbox',checked:selected.includes(compound.row_id),onChange:event=>setSelected(event.target.checked?[...selected,compound.row_id]:selected.filter(id=>id!==compound.row_id))})),e('td',{},[e('button',{className:'link-button',onClick:()=>openDetail(compound.row_id)},compound.name),e('div',{className:'mono small'},compound.compound_id)]),e('td',{className:'thumbnail'},[Svg({src:compound.version?.svg}),StatusBadge({type:status.structure||'NOT_STARTED'})]),e('td',{},StatusBadge({type:status.properties||'NOT_RUN'})),e('td',{},StatusBadge({type:status.activity||'NOT_RUN'})),e('td',{},StatusBadge({type:status.admet||'NOT_RUN'})),e('td',{},StatusBadge({type:status.optimization||'NOT_RUN'})),e('td',{},e('button',{className:'secondary',onClick:()=>openDetail(compound.row_id)},'Open'))])}))])):e('div',{className:'empty-state'},[e('h3',{},'No compounds yet'),e('p',{},'Add the first compound by name; structure and calculation may follow later.'),e('button',{onClick:()=>setAddCompoundOpen(true)},'Add Compound')])])
   ]),
   project&&projectTab==='compounds'&&detail&&compoundDetail(),
   project&&projectTab==='assays'&&e('div',{className:'card',key:'assays'},[e('h2',{},'Assays and Activity'),e('p',{},'Define assays, enter experimental activity, train project QSAR, and inspect SAR in the existing workbench.'),e('a',{className:'button',href:'/static/stage2-workbench.html?project='+projectId},'Open Activity Workbench')]),
   project&&projectTab==='compare'&&e(React.Fragment,{key:'compare'},ComparePanel()),
   project&&projectTab==='optimization'&&e('div',{className:'card',key:'optimization'},[e('h2',{},'Optimization Runs'),e('p',{},'Optimization is CompoundVersion-specific. Open a compound and select Optimization to analyze or generate proposals.'),e('button',{className:'secondary',onClick:()=>setProjectTab('compounds')},'Choose a Compound')]),
   project&&projectTab==='settings'&&e(React.Fragment,{key:'settings'},SettingsPanel())
  ]);
 }

 const sidebarGroups=[
  ['WORKSPACE',[['Dashboard',()=>goDashboard()]]],
  ['CHEMISTRY',[['Structure',()=>navigateScientific('compounds','overview')],['Properties',()=>navigateScientific('compounds','properties')]]],
  ['BIOLOGY',[['Activity / SAR',()=>navigateScientific('assays')]]],
  ['ADMET',[['Absorption',()=>navigateScientific('compounds','admet')],['Distribution',()=>navigateScientific('compounds','admet')],['Metabolism',()=>navigateScientific('compounds','metabolism')],['CYP',()=>navigateScientific('compounds','admet')],['Transporters',()=>navigateScientific('compounds','admet')]]],
  ['SAFETY',[['Toxicology',()=>navigateScientific('compounds','admet')]]],
  ['OPTIMIZATION',[['Optimization',()=>navigateScientific(detail?'compounds':'optimization','optimization')]]]
 ];
 const sidebar=e('aside',{className:'sidebar'+(sidebarOpen?' open':''),key:'sidebar'},[
  e('div',{className:'sidebar-head',key:'head'},[e('div',{},[e('button',{className:'brand-button',onClick:goDashboard},'Drug Optimization Platform'),e('div',{className:'tag'},'CompoundVersion research workspace')]),e('button',{className:'menu-toggle',onClick:()=>setSidebarOpen(value=>!value),'aria-expanded':sidebarOpen,'aria-label':'Toggle scientific navigation'},sidebarOpen?'Close':'Menu')]),
  e('div',{className:'sidebar-body',key:'body'},[
   e('nav',{className:'scientific-nav',key:'nav','aria-label':'Scientific navigation'},sidebarGroups.map(([group,items])=>e('section',{className:'nav-group',key:group},[e('h3',{},group),...items.map(([label,action])=>e('button',{key:label,className:label==='Dashboard'&&projectTab==='dashboard'?'active':'',onClick:action},label))]))),
   e('section',{className:'nav-group development-nav',key:'development'},[e('h3',{},'DEVELOPMENT'),e('button',{disabled:true},['PK / DMPK ',StatusBadge({type:'PLANNED'})])]),
   e('section',{className:'current-project',key:'current'},[e('h3',{},'CURRENT PROJECT'),projects.length?e(React.Fragment,{},[e('select',{value:projectId||'',onChange:event=>{setProjectId(Number(event.target.value));setSidebarOpen(false)}},projects.map(item=>e('option',{key:item.id,value:item.id},item.name))),e('div',{className:'current-project-meta'},[(project?.target||projects.find(item=>item.id===projectId)?.target||'Target not set'),' · ',(project?.molecule_type||projects.find(item=>item.id===projectId)?.molecule_type||'Small Molecule')]),e('button',{className:'secondary',onClick:()=>openProject(projectId)},'Open Current Project')]):e('p',{className:'tag'},'Create a project from the dashboard.')])
  ])
 ]);
 return e('div',{className:'shell'},[sidebar,e('main',{className:'content',key:'content'},[
  projectTab==='dashboard'?MainDashboard():ProjectWorkspace(),
  AddCompoundPanel(),
  message&&e('pre',{className:'card error',key:'message'},message)
 ])]);

 function propertyTable(properties){
  const names={molecular_formula:'Molecular formula',molecular_weight:'MW',exact_molecular_weight:'Exact MW',clogp:'cLogP',tpsa:'TPSA',hbd:'HBD',hba:'HBA',rotatable_bonds:'Rotatable bonds',heavy_atom_count:'Heavy atoms',heteroatom_count:'Heteroatoms',ring_count:'Rings',aromatic_ring_count:'Aromatic rings',fraction_csp3:'Fsp3',formal_charge:'Charge',molar_refractivity:'MR',aromatic_proportion:'Aromatic proportion',molecular_flexibility:'Flexibility'};
  return e('div',{},[
   e('h4',{key:'title'},'Physicochemical properties'),
   e('table',{key:'table'},e('tbody',{},Object.entries(names).map(([key,label])=>e('tr',{key},[e('td',{key:'label'},label),e('td',{key:'value',className:'num mono'},String(properties[key]??'-'))]))))
  ]);
 }
 function alertList(alerts){
  return alerts.length?alerts.map(alert=>e('div',{key:alert.alert_set+alert.alert_name,className:'alert'},[e('strong',{key:'name'},alert.alert_name),e('div',{key:'reason'},alert.reason),e('div',{key:'set',className:'small'},'Set: '+alert.alert_set+' · atoms: '+(alert.matched_atoms.join(', ')||'not exposed by RDKit'))])):[e('p',{key:'none'},'None detected')];
 }
 function ruleTable(rules){
  const rows=Object.entries(rules).map(([name,rule])=>e('div',{key:name},[
   Badge({ok:rule.result==='PASS',text:rule.result+' · '+name}),
   rule.reasons.length>0&&e('ul',{key:'reasons'},rule.reasons.map(reason=>e('li',{key:reason},reason)))
  ]));
  return e('div',{},[e('h4',{key:'title'},'Drug-likeness'),...rows]);
 }
 function pairPlot(data,x,y){
  const points=data.compounds.filter(compound=>compound[x]!=null&&compound[y]!=null);if(!points.length)return null;
  const xs=points.map(compound=>compound[x]),ys=points.map(compound=>compound[y]);
  const xMin=Math.min(...xs)*.95,xMax=Math.max(...xs)*1.05,yMin=Math.min(...ys)*.95,yMax=Math.max(...ys)*1.05;
  const xRange=xMax-xMin||1,yRange=yMax-yMin||1;
  return e('div',{className:'card'},[e('h4',{key:'title'},x+' vs '+y),e('div',{key:'plot',style:{height:'160px',position:'relative',borderBottom:'1px solid #ccc',borderLeft:'1px solid #ccc'}},points.map((compound,index)=>{const left=(compound[x]-xMin)/xRange*82+6,top=(yMax-compound[y])/yRange*70+8;return e('div',{key:index,title:compound.compound+' '+x+'='+compound[x]+' '+y+'='+compound[y],onClick:()=>openDetail(compound.row_id),style:{position:'absolute',left:left+'%',top:top+'%',width:'11px',height:'11px',background:'#1769aa',borderRadius:'50%',cursor:'pointer'}})})),e('div',{key:'legend',className:'small'},'Click a point to open the compound · ranges '+xMin.toFixed(1)+'–'+xMax.toFixed(1)+' / '+yMin.toFixed(1)+'–'+yMax.toFixed(1))]);
 }
 function qedHistogram(data){
  const bins=[0,0,0,0,0];data.compounds.forEach(compound=>{if(compound.QED!=null)bins[Math.min(Math.floor(compound.QED/.2),4)]++});const maximum=Math.max(1,...bins);
  return e('div',{className:'card'},[e('h4',{key:'title'},'QED distribution'),...bins.map((count,index)=>e('div',{key:index,style:{display:'flex',alignItems:'center',gap:'8px',margin:'7px 0'}},[e('div',{key:'label',className:'small',style:{width:'52px'}},(index*.2).toFixed(1)+'-'+((index+1)*.2).toFixed(1)),e('div',{key:'bar',style:{flex:1,background:'#e5edf5'}},e('div',{style:{width:(count/maximum*100)+'%',background:'#1769aa',height:'16px'}})),e('span',{key:'count',className:'small'},count)]))]);
 }
}
ReactDOM.createRoot(document.getElementById('root')).render(e(App));
})();
