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
 const [globalView,setGlobalView]=useState('dashboard');
 const [projectSelection,setProjectSelection]=useState([]),[deleteProjects,setDeleteProjects]=useState([]),[deleteConfirmations,setDeleteConfirmations]=useState({}),[deleteBusy,setDeleteBusy]=useState(false);
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
 const [predictionWorkflow,setPredictionWorkflow]=useState(null);
 const [optimizationWorkspace,setOptimizationWorkspace]=useState({project_id:'',compound_id:''});
 const [optimizationForm,setOptimizationForm]=useState({
  assay_id:'',objectives:['Balanced optimization'],custom_objective:'',
  constraints:{potency_max_nm:'',do_not_worsen_fold:'2',clogp_max:'4',tpsa_min:'40',tpsa_max:'100',mw_max:'550',similarity_min:'0.6',logs_min:'-4',caco2_logpapp_min:'-5.5',herg_do_not_increase:true},endpoint_weights:{}
 });
 const [workspace,setWorkspace]=useState(null),[experimentalOpen,setExperimentalOpen]=useState(false),[experimentalSelected,setExperimentalSelected]=useState([]),[experimentalDrafts,setExperimentalDrafts]=useState({});
 const [compareMetrics,setCompareMetrics]=useState(['MW','cLogP','TPSA','QED','Activity','Solubility','Caco-2','PPB','HLM','RLM','hERG','Ames','DILI']),[compareAssay,setCompareAssay]=useState('');
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
 const saveCompound=async predict=>{
  setAdmetBusy(!!predict);setPredictionWorkflow(predict?{status:'RUNNING',steps:{overview:{status:'PENDING'},properties:{status:'PENDING'},admet:{status:'PENDING'},metabolism:{status:'PENDING'}}}:null);
  try{
   const saved=await api.post('/projects/'+projectId+'/compounds',{...compoundForm,calculate:false});
   setCompoundForm({compound_id:'',name:'',smiles:'',notes:''});setPreview(null);setAddCompoundOpen(false);
   let workflow=null;
   if(predict){workflow=await api.post('/compounds/'+saved.row_id+'/predict-workflow',{});setPredictionWorkflow(workflow)}
   await Promise.all([loadProject(projectId),loadProjects(),loadDashboard()]);setMessage(predict?(workflow.message+' Activity was not run.'):'Compound saved without prediction');
   await openDetail(saved.row_id);
  }catch(error){setPredictionWorkflow(current=>current?{...current,status:'FAILED',message:String(error)}:null);setMessage(String(error))}finally{setAdmetBusy(false)}
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
     e('div',{className:'col-4',key:'assay'},[e('label',{},'Selected assay'),e('select',{value:optimizationForm?.assay_id||'',onChange:event=>setOptimizationForm(current=>({...current,assay_id:event.target.value}))},[e('option',{key:'none',value:''},'No assay selected'),...assays.map(assay=>e('option',{key:assay.id,value:assay.id},assay.name+' · '+(assay.measurement_type||'')))])]),
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
  const [adminType, setAdminType] = React.useState('IV_BOLUS');
  const [dose, setDose] = React.useState(5.0);
  const [doseUnit, setDoseUnit] = React.useState('mg/kg');
  const [infusionDur, setInfusionDur] = React.useState(1.0);
  const [frequency, setFrequency] = React.useState('Single Dose');
  const [interval, setInterval] = React.useState(24.0);
  const [numDoses, setNumDoses] = React.useState(3);
  const [modelType, setModelType] = React.useState('ONE_COMPARTMENT');
  const [logScale, setLogScale] = React.useState(false);
  const [preview, setPreview] = React.useState(null);
  const [activeRun, setActiveRun] = React.useState(null);
  const [history, setHistory] = React.useState([]);
  const [loading, setLoading] = React.useState(false);
  const [error, setError] = React.useState(null);

  const loadData = React.useCallback(async()=>{
   if(!versionId) return;
   try{
    const prev = await api.get('/compound-versions/'+versionId+'/pk-simulation/preview?species='+species);
    setPreview(prev);
    const hist = await api.get('/compound-versions/'+versionId+'/pk-simulation/history?species='+species);
    setHistory(hist||[]);
    if(hist && hist.length > 0 && !activeRun){
     setActiveRun(hist[0]);
    }
   }catch(err){
    console.error("Simulation load error:", err);
   }
  },[versionId, species]);

  React.useEffect(()=>{ loadData(); },[loadData]);

  const handleRun = async()=>{
   setLoading(true);
   setError(null);
   try{
    const payload = {
     species,
     administration_type: adminType,
     dose: parseFloat(dose),
     dose_unit: doseUnit,
     infusion_duration_hours: adminType === 'IV_INFUSION' ? parseFloat(infusionDur) : 0.0,
     dosing_frequency: frequency,
     dose_interval_hours: parseFloat(interval),
     num_doses: frequency === 'Repeated Dosing' ? parseInt(numDoses, 10) : 1,
     model_type: modelType
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

  return e('div',{className:'card', style:{marginTop:'16px'}},[
   e('div',{className:'row toolbar'},[
    e('h3',{},'PK SIMULATION — IV Concentration-Time Engine (Stage 5B-1)'),
    StatusBadge({type: activeRun ? activeRun.confidence : (preview?.confidence_ceiling||'MEDIUM')})
   ]),
   e('p',{className:'small'},'Mechanistic mathematical simulation of IV bolus and IV infusion concentration-time profiles using Stage 5A PK parameters.'),

   e('div',{className:'grid', style:{marginTop:'12px'}},[
    e('div',{className:'col-3'},[
     e('label',{},'Species'),
     e('select',{value:species, onChange:ev=>setSpecies(ev.target.value)},['Rat','Mouse','Dog','Monkey','Human'].map(s=>e('option',{key:s,value:s},s)))
    ]),
    e('div',{className:'col-3'},[
     e('label',{},'Administration Type'),
     e('select',{value:adminType, onChange:ev=>setAdminType(ev.target.value)},[
      e('option',{value:'IV_BOLUS'},'IV Bolus'),
      e('option',{value:'IV_INFUSION'},'IV Infusion')
     ])
    ]),
    e('div',{className:'col-3'},Field({label:'Dose', type:'number', value:dose, onChange:setDose})),
    e('div',{className:'col-3'},[
     e('label',{},'Dose Unit'),
     e('select',{value:doseUnit, onChange:ev=>setDoseUnit(ev.target.value)},['mg/kg','µg/kg'].map(u=>e('option',{key:u,value:u},u)))
    ]),
    adminType === 'IV_INFUSION' && e('div',{className:'col-3', key:'inf-dur'},Field({label:'Infusion Duration (h)', type:'number', value:infusionDur, onChange:setInfusionDur})),
    e('div',{className:'col-3'},[
     e('label',{},'Dosing Frequency'),
     e('select',{value:frequency, onChange:ev=>setFrequency(ev.target.value)},['Single Dose','Repeated Dosing'].map(f=>e('option',{key:f,value:f},f)))
    ]),
    frequency === 'Repeated Dosing' && e('div',{className:'col-3', key:'interval'},Field({label:'Dose Interval τ (h)', type:'number', value:interval, onChange:setInterval})),
    frequency === 'Repeated Dosing' && e('div',{className:'col-3', key:'numDoses'},Field({label:'Number of Doses', type:'number', value:numDoses, onChange:setNumDoses})),
    e('div',{className:'col-3'},[
     e('label',{},'Model Type'),
     e('select',{value:modelType, onChange:ev=>setModelType(ev.target.value)},[
      e('option',{value:'ONE_COMPARTMENT'},'1-Compartment Model'),
      e('option',{value:'TWO_COMPARTMENT'},'2-Compartment Model (If fit/parameters available)')
     ])
    ])
   ]),

   e('div',{className:'card', style:{background:'var(--bg-subtle,#1e293b)', marginTop:'12px', padding:'12px'}},[
    e('strong',{style:{fontSize:'14px'}},'Parameter Review Before Run:'),
    e('div',{className:'row toolbar', style:{marginTop:'6px'}},[
     e('span',{},'CL: '+(clPreview?.value!=null ? clPreview.value+' mL/min/kg ('+clPreview.source+')' : 'Unavailable')),
     e('span',{},'Volume: '+(vPreview?.value!=null ? vPreview.value+' L/kg ('+vPreview.type+')' : 'Unavailable')),
     e('button',{className:'primary', onClick:handleRun, disabled:loading}, loading ? 'Simulating...' : 'RUN SIMULATION')
    ]),
    (preview?.warnings||[]).map((w, idx)=>e('div',{key:idx, className:'small alert', style:{marginTop:'4px'}},w)),
    error && e('div',{className:'small alert', style:{marginTop:'6px', color:'#ef4444'}},error)
   ]),

   activeRun && e('div',{style:{marginTop:'16px'}},[
    e('div',{className:'row toolbar'},[
     e('h4',{},'CALCULATED PK SIMULATION: '+(activeRun.administration_type==='IV_INFUSION'?'IV Infusion':'IV Bolus')+' ('+activeRun.species+')'),
     e('div',{},[
      e('button',{className: logScale ? 'secondary' : 'primary', style:{marginRight:'6px', padding:'4px 8px'}, onClick:()=>setLogScale(false)},'Linear'),
      e('button',{className: logScale ? 'primary' : 'secondary', style:{padding:'4px 8px'}, onClick:()=>setLogScale(true)},'Semi-Log')
     ])
    ]),
    e('div',{style:{marginTop:'10px'}},renderPlot(activeRun)),
    e('div',{className:'row toolbar', style:{marginTop:'4px', fontSize:'12px', color:'#94a3b8'}},[
     e('span',{},'── Blue Line: Calculated PK Simulation'),
     e('span',{},'● Red Dots: Experimental Observed Points')
    ]),
    e('div',{className:'grid ivive-output-grid', style:{marginTop:'12px'}},[
     e('div',{className:'card pk-nca-card'},[
      e('span',{},'Cmax / C0'),
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
      e('span',{},'AUC Numerical Cross-Check'),
      e('strong',{className:'mono'},activeRun.output_metrics?.auc_inf_numerical_ng_h_ml+' ng·h/mL ('+activeRun.output_metrics?.auc_agreement_pct+'% match)')
     ]),
     e('div',{className:'card pk-nca-card'},[
      e('span',{},'Terminal Half-Life'),
      e('strong',{className:'mono'},activeRun.output_metrics?.half_life_hours+' h')
     ])
    ]),
    activeRun.residuals && activeRun.residuals.length > 0 && e('div',{style:{marginTop:'16px'}},[
     e('h4',{},'Experimental Observation Overlay & Residual Analysis'),
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

 function compoundDetail(){
  if(!detail)return null;
  const version=detail.version,detailMeasurements=admet?.measurements||[],detailRuns=admet?.prediction_runs||[],detailPredictions=admet?.predictions||[];
  const activity=workspace?.activity||{measurements:[],predictions:[]},properties=version?.properties||{};
  const tabs=['overview','properties','activity','admet','metabolism','pk','history'];
  const highlights=detailPredictions.filter((row,index,array)=>array.findIndex(item=>item.endpoint===row.endpoint)===index).slice(0,5);
  const activityTable=e('div',{},[
   e('h3',{key:'exp'},'Experimental Activity'),activity.measurements.length?e('table',{key:'exp-table'},[e('thead',{},e('tr',{},['Assay','Measurement','Value','Source'].map(x=>e('th',{key:x},x)))),e('tbody',{},activity.measurements.map(row=>e('tr',{key:row.id},[e('td',{},row.assay),e('td',{},row.measurement_type),e('td',{className:'mono'},row.qualifier+' '+row.value+' '+row.unit),e('td',{},[StatusBadge({type:'Experimental'}),' '+row.source])])))]):e('div',{className:'empty-state'},[StatusBadge({type:'Not measured'}),e('p',{},'No experimental activity measurement entered.')]),
   e('h3',{key:'pred',style:{marginTop:'22px'}},'Activity Prediction'),activity.predictions.length?e('table',{key:'pred-table'},[e('thead',{},e('tr',{},['Assay','Predicted value','Confidence','Domain'].map(x=>e('th',{key:x},x)))),e('tbody',{},activity.predictions.map(row=>e('tr',{key:row.id},[e('td',{},row.assay),e('td',{className:'mono'},row.predicted_value_nm+' nM'),e('td',{},row.confidence),e('td',{},row.applicability_domain)])))]):e('div',{className:'empty-state'},[StatusBadge({type:'Not predicted'}),e('h4',{},'Assay configuration required'),e('p',{},'Activity is intentionally excluded from Save & Predict. Configure assay type, conditions, species, cell line, and mutation where applicable.'),e('a',{className:'button secondary',href:'/static/stage2-workbench.html?project='+projectId},'Run Activity Prediction')])
  ]);
  return e('div',{className:'compound-workspace'},[
   e('div',{className:'card compound-hero',key:'hero'},[e('div',{className:'compound-hero-structure'},Svg({src:version?.highlighted_svg||version?.svg})),e('div',{className:'compound-hero-copy'},[e('div',{className:'eyebrow'},'COMPOUND DETAIL'),e('h2',{},detail.name),e('div',{className:'row'},[StatusBadge({type:detail.status}),e('span',{className:'mono'},detail.compound_id+(version?' · Version '+version.version_number:' · No structure version'))]),e('p',{className:'small'},workspace?'Strict scope: Project #'+workspace.scope.project_id+' · Compound #'+workspace.scope.compound_id+' · CompoundVersion #'+workspace.scope.version_id:'Draft compound; no version-linked data exists.'),e('div',{className:'row'},[version&&e('button',{className:'secondary',onClick:updateStructure},'Modify Structure / New Version'),e('button',{className:'secondary',onClick:()=>setDetail(null)},'Back to Compounds')])])]),
   e('nav',{className:'detail-tabs',key:'tabs'},tabs.map(tab=>e('button',{key:tab,className:detailTab===tab?'':'secondary',disabled:!version&&['properties','activity','admet','metabolism','pk'].includes(tab),onClick:()=>setDetailTab(tab)},tab.toUpperCase()))),
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
    e('div',{className:'card col-12'},[e('div',{className:'row toolbar'},[e('h3',{},'ADMET Highlights'),e('button',{className:'secondary',onClick:()=>{setOptimizationWorkspace({project_id:String(projectId),compound_id:String(detail.row_id)});openGlobalView('optimization')}},'Open in Optimization Workspace')]),highlights.length?e('div',{className:'highlight-grid'},highlights.map(row=>e('div',{key:row.endpoint,className:'highlight-item'},[e('strong',{},row.endpoint==='Permeability'?'Caco-2 Permeability':row.endpoint),e('div',{className:'mono'},row.predicted_value+' '+row.unit),e('div',{className:'small'},'Model: '+row.model?.model_name),StatusBadge({type:'Predicted'})]))):e('div',{className:'empty-state'},[StatusBadge({type:'Not predicted'}),e('p',{},'No ADMET predictions run for this CompoundVersion.')])]),
    predictionWorkflow&&predictionWorkflow.compound_id===detail.row_id&&e('div',{className:'card col-12 prediction-workflow-status'},[e('h3',{},'Save & Predict Workflow'),e('div',{className:'workflow-strip'},Object.entries(predictionWorkflow.steps||{}).map(([name,row])=>e('div',{className:'workflow-step',key:name},[e('span',{},name==='admet'?'ADMET / CYP / Transporter / Safety':name),StatusBadge({type:row.status}),row.message&&e('small',{},row.message)])))])
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
    experimentalOpen&&e('div',{className:'card'},ExperimentalDataPanel()),
    e('section',{className:'card',key:'experimental-results'},[e('div',{className:'eyebrow'},'1 · EXPERIMENTAL RESULTS'),e('h3',{},'Experimental Results'),detailMeasurements.length?admetMeasurementTable(detailMeasurements):e('div',{className:'empty-state'},[StatusBadge({type:'Not measured'}),e('p',{},'No experimental measurement entered.'),e('button',{className:'secondary',onClick:()=>setExperimentalOpen(true)},'Add Experimental Data')])]),
    e('section',{className:'card',key:'prediction-results'},[e('div',{className:'eyebrow'},'2 · PREDICTION RESULTS'),e('h3',{},'Prediction Results · Consensus and Individual Models'),e('p',{className:'small'},'Probability and confidence remain distinct. Each model result is preserved; a consensus never overwrites it.'),consensusPredictionPanel(version.id)]),
    e('section',{className:'card',key:'comparison-results'},[e('div',{className:'eyebrow'},'3 · EXPERIMENTAL VS PREDICTION'),e('h3',{},'Experimental vs Prediction'),experimentalComparisonPanel(version.id)]),
    e('section',{key:'integrated'},[e('div',{className:'eyebrow'},'4 · INTEGRATED PROFILE'),integratedProfile(version.id)]),
    e('section',{className:'card',key:'provenance'},[
     e('div',{className:'eyebrow'},'5 · MODEL / PROVENANCE DETAILS'),e('h3',{},'Model Registry and Availability'),unavailableModelsCollapsed(),
     e('details',{style:{marginTop:'12px'}},[
      e('summary',{},'Available model registry entries'),
      e('table',{},[
       e('thead',{},e('tr',{},['Endpoint','Model','Version','Output','Species'].map(label=>e('th',{key:label},label)))),
       e('tbody',{},(admet?.models||[]).filter(model=>model.active).map(model=>e('tr',{key:model.id},[
        e('td',{},model.endpoint),e('td',{},model.model_name),e('td',{},model.model_version),e('td',{},model.output_unit),e('td',{},model.species||'Not specified')
       ])))
      ])
     ])
    ])
   ]),
   detailTab==='metabolism'&&metabolismProfile(version.id),
   detailTab==='pk'&&pkProfile(version.id),
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
   e('div',{className:'row modal-actions',key:'actions'},[e('button',{className:'secondary',disabled:admetBusy||!compoundForm.name.trim(),onClick:()=>saveCompound(false)},'Save'),e('button',{disabled:admetBusy||!compoundForm.name.trim()||!compoundForm.smiles.trim()||!smallMolecule,onClick:()=>saveCompound(true)},admetBusy?'Saving & predicting…':'Save & Predict'),e('span',{className:'small'},'Save stores identity and structure only. Save & Predict runs Properties, ADMET, and Metabolism; Activity requires an assay.')])
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
  const projectPerformance=new Map((admet?.model_performance||[]).filter(row=>row.scope==='PROJECT:'+projectId).map(row=>[row.model_id,row]));
  return e('div',{},[
   e('div',{className:'card',key:'identity'},[e('h2',{},'Project Settings'),e('p',{className:'small'},'Indication and mechanism remain preserved in the database but are kept out of the primary creation workflow.'),e('div',{className:'grid'},[e('div',{className:'col-4'},Field({label:'Project Name',value:project.name,onChange:value=>setProject(current=>({...current,name:value}))})),e('div',{className:'col-4'},Field({label:'Target',value:project.target,onChange:value=>setProject(current=>({...current,target:value}))})),e('div',{className:'col-4'},[e('label',{},'Molecule Type'),e('select',{value:project.molecule_type,onChange:event=>setProject(current=>({...current,molecule_type:event.target.value}))},['Small Molecule','Peptide'].map(value=>e('option',{key:value,value},value)))]),e('div',{className:'col-12'},Field({label:'Description',value:project.description||'',onChange:value=>setProject(current=>({...current,description:value})),type:'textarea'}))]),e('button',{disabled:!project.name.trim()||!project.target.trim(),onClick:saveProjectSettings},'Save Project Settings')]),
   e('div',{className:'card',key:'models'},[
    e('h2',{},'Prediction Models'),e('p',{className:'small'},'Multiple registered models may share an endpoint. Project performance influences consensus only from N ≥ 10; N ≥ 30 enables a stronger blend.'),
    (admet?.models||[]).length?e('div',{className:'table-scroll'},e('table',{},[
     e('thead',{},e('tr',{},['Endpoint','Model','Version','Status','Training N','Validation','Project Experimental N','Project MAE / Accuracy','Consensus Weight','Project Selection'].map(label=>e('th',{key:label},label)))),
     e('tbody',{},admet.models.map(model=>{
      const performance=projectPerformance.get(model.id),validation=model.validation||model.details?.validation||{},metric=performance?.metrics?.mae??performance?.metrics?.accuracy;
      const best=admet?.best_project_models?.[model.endpoint]?.model_id===model.id;
      return e('tr',{key:model.id},[e('td',{},model.endpoint),e('td',{},model.model_name),e('td',{className:'mono'},model.model_version),e('td',{},StatusBadge({type:model.active?'READY':'MODEL_UNAVAILABLE'})),e('td',{className:'mono'},model.details?.training_n??'—'),e('td',{className:'small'},Object.entries(validation).slice(0,2).map(([key,value])=>key+' '+value).join(' · ')||'Not reported'),e('td',{className:'mono'},performance?.n??0),e('td',{className:'mono'},metric==null?'Insufficient experimental data':Number(metric).toFixed(3)),e('td',{className:'mono'},performance?.n>=10?Number(performance.performance_factor).toFixed(3):'Published validation'),e('td',{},best?e('span',{className:'pass'},'Best performing model for this project'):'Insufficient experimental data')]);
     }))
    ])):Empty({children:'Select a project to inspect its model registry.'})
   ]),
   e('div',{className:'card',key:'csv'},[e('h2',{},'Experimental ADMET CSV'),e('p',{className:'small'},'Advanced project-wide import/export. Compound Detail remains CompoundVersion-isolated.'),e('a',{className:'button secondary',href:'/api/projects/'+projectId+'/admet/export.csv'},'Export CSV'),e('textarea',{rows:7,value:admetCsv,placeholder,onChange:event=>{setAdmetCsv(event.target.value);setAdmetCsvPreview(null)}}),e('div',{className:'row',style:{marginTop:'10px'}},[e('button',{className:'secondary',disabled:admetBusy||!admetCsv.trim(),onClick:previewAdmet},'Preview CSV'),e('button',{disabled:admetBusy||!admetCsvPreview||admetCsvPreview.errors.length>0||!admetCsvPreview.valid_count,onClick:importAdmet},'Import Valid Rows')]),admetCsvPreview&&e('p',{className:admetCsvPreview.errors.length?'fail':'pass'},admetCsvPreview.valid_count+' valid · '+admetCsvPreview.errors.length+' errors')]),
   e('div',{className:'card danger-zone',key:'delete'},[e('div',{},[e('div',{className:'eyebrow'},'DANGER ZONE'),e('h2',{},'Delete Project'),e('p',{className:'small'},'Deletion requires a separate confirmation and the exact project name. No other project is included.')]),e('button',{className:'danger',onClick:()=>openDeleteDialog([selectedProjectSummary||project])},'Delete Project…')])
  ]);
 }

 const openGlobalView=view=>{setGlobalView(view);setProjectTab('dashboard');setDetail(null);setAddCompoundOpen(false);setComparison(null);setSelectedCandidate(null);setSidebarOpen(false);loadDashboard().catch(error=>setMessage(String(error)))};
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
 const deleteNamesMatch=deleteProjects.length>0&&deleteProjects.every(item=>deleteConfirmations[item.id]===item.name);
 const confirmProjectDeletion=async()=>{
  if(!deleteNamesMatch)return;
  setDeleteBusy(true);
  try{
   const confirmations=deleteProjects.map(item=>({id:item.id,confirmation_name:deleteConfirmations[item.id]}));
   const result=confirmations.length===1
    ?await api.del('/projects/'+confirmations[0].id,{confirmation_name:confirmations[0].confirmation_name})
    :await api.post('/projects/bulk-delete',{projects:confirmations});
   const deletedIds=result.deleted_project_ids||confirmations.map(item=>item.id),currentDeleted=deletedIds.includes(projectId);
   const [rows,summary]=await Promise.all([api.get('/projects'),api.get('/dashboard')]);
   setProjects(rows);setDashboard(summary);setProjectSelection([]);setDeleteProjects([]);setDeleteConfirmations({});
   setGlobalView('dashboard');setProjectTab('dashboard');setDetail(null);setWorkspace(null);setAdmet(null);setMetabolism(null);setComparison(null);setSelected([]);setSelectedCandidate(null);
   if(currentDeleted){setProjectId(null);setProject(null)}
   setMessage((result.deleted_project_names||deleteProjects.map(item=>item.name)).join(', ')+' deleted');
  }catch(error){setMessage('Project deletion failed: '+error.message)}finally{setDeleteBusy(false)}
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
  const version=detail?.row_id===selectedCompound?.row_id?detail.version:null;
  const predictions=version?(admet?.predictions||[]).filter(row=>row.version_id===version.id):[];
  return e('div',{className:'optimization-workspace'},[
   e('section',{className:'card',key:'selector'},[e('div',{className:'eyebrow'},'DETERMINISTIC MEDICINAL CHEMISTRY'),e('h1',{},'Optimization Workspace'),e('p',{className:'small'},'Select a project and a CompoundVersion, review its evidence, then reuse the existing Stage 4A strategy and Stage 4B analog engines. No LLM and no PK are used.'),
    e('div',{className:'optimization-workspace-steps'},[
     e('div',{className:'optimization-workspace-step',key:'project'},[e('h3',{},'Step 1 — Select Project'),e('select',{value:optimizationWorkspace.project_id,onChange:event=>{setOptimizationWorkspace({project_id:event.target.value,compound_id:''});setDetail(null)}},[e('option',{value:''},'Select project'),...projectChoices.map(row=>e('option',{key:row.id,value:row.id},row.name+' · '+(row.target||'Target not set')))])]),
     e('div',{className:'optimization-workspace-step',key:'compound'},[e('h3',{},'Step 2 — Select Compound'),e('select',{value:optimizationWorkspace.compound_id,disabled:!selectedProject,onChange:event=>setOptimizationWorkspace(current=>({...current,compound_id:event.target.value}))},[e('option',{value:''},selectedProject?'Select compound':'Select a project first'),...compounds.map(row=>e('option',{key:row.row_id,value:row.row_id,disabled:!row.version},row.name+' · '+row.compound_id+(row.version?' · v'+row.current_version:' · Draft')))])])
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
   selectedProject&&!selectedCompound&&e('section',{className:'card empty-state',key:'choose'},[e('h3',{},'Choose a parent compound'),e('p',{},'The selected project has '+compounds.length+' compound record(s). Drafts without structure cannot be optimized.')]),
   !selectedProject&&e('section',{className:'card empty-state',key:'no-project'},[e('h3',{},'No project selected'),e('p',{},'Select a project from the dropdown above to open its optimization workspace.')])
  ]);
 }

 function MainDashboard(){
  const registry=dashboard?.model_registry||[];
  const registryStatus=(endpoint,defaultFallback='MODEL_UNAVAILABLE')=>{
   const model=registry.find(row=>row.endpoint===endpoint);
   if(!model)return defaultFallback;
   return model.active?(model.status||'READY'):'MODEL_UNAVAILABLE';
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
  const home=globalView==='dashboard';
  return e(React.Fragment,{},[
   globalView!=='dashboard'&&globalView!=='optimization'&&e('section',{className:'card global-view-header',key:'global-head'},[e('div',{className:'eyebrow'},'WORKSPACE'),e('h1',{},({"new-project":'New Project',projects:'Projects',settings:'Settings',help:'Help'})[globalView]||'Workspace')]),
   globalView==='optimization'&&e(GlobalOptimizationWorkspace,{key:'optimization-workspace'}),
   home&&e('section',{className:'card dashboard-hero',key:'intro'},[
    e('div',{className:'eyebrow'},'PLATFORM OVERVIEW'),e('h1',{},'Drug Optimization Platform'),
    e('p',{},'Structure, activity, ADMET and medicinal chemistry optimization data are integrated at the compound-version level to support hit-to-lead and lead optimization decisions.'),
    e('ul',{className:'dashboard-capabilities'},['Structure-based compound management','Experimental data integration','Predictive ADMET','SAR / optimization workflow','Full prediction provenance'].map(item=>e('li',{key:item},item))),
    e('div',{className:'dashboard-stats'},[
     e('div',{className:'dashboard-stat',key:'projects'},[e('span',{},'Projects'),e('strong',{},String(dashboard?.totals?.projects??projects.length))]),
     e('div',{className:'dashboard-stat',key:'compounds'},[e('span',{},'Compounds'),e('strong',{},String(dashboard?.totals?.compounds??projects.reduce((sum,row)=>sum+(row.compound_count||0),0)))]),
     e('div',{className:'dashboard-stat',key:'scope'},[e('span',{},'Default data scope'),e('strong',{className:'dashboard-stat-text'},'CompoundVersion'),e('small',{},'Project-isolated')])
    ])
   ]),
   home&&e('section',{className:'dashboard-section',key:'modules'},[
    e('div',{className:'section-heading'},[e('div',{},[e('div',{className:'eyebrow'},'SCIENTIFIC WORKSPACE'),e('h2',{},'Available Scientific Modules')]),e('p',{className:'small'},'Status reflects the current local engine and model registry.')]),
    e('div',{className:'module-grid'},modules.map(module=>e('article',{className:'module-card',key:module.title},[
     e('div',{className:'module-card-head'},[e('h3',{},module.title),StatusBadge({type:module.status})]),e('p',{className:'small'},module.description),
     e('ul',{className:'module-list'},module.items.map(([label,status])=>e('li',{key:label},[e('span',{},label),StatusBadge({type:status})]))),
     module.unavailable?.length>0&&e('div',{className:'module-unavailable'},['Unavailable: ',module.unavailable.join(' · ')])
    ])))
   ]),
   (home||globalView==='new-project')&&e('div',{className:'dashboard-split',key:'start'},[
    e('section',{className:'card dashboard-create',key:'create'},[e('div',{className:'eyebrow'},'NEW WORKSPACE'),e('h2',{},'Create New Project'),e('div',{className:'create-project-grid'},[
     e(Field,{label:'Project Name *',value:form.name,onChange:value=>setForm({...form,name:value}),placeholder:'EGFR Exon20ins'}),e(Field,{label:'Target *',value:form.target,onChange:value=>setForm({...form,target:value}),placeholder:'EGFR'}),e('div',{},[e('label',{},'Molecule Type'),e('select',{value:form.molecule_type,onChange:event=>setForm({...form,molecule_type:event.target.value})},['Small Molecule','Peptide'].map(value=>e('option',{key:value,value},value)))])
    ]),e('button',{disabled:!form.name.trim()||!form.target.trim(),onClick:createProject},'Create Project'),e('p',{className:'small dashboard-note'},'Description and additional metadata can be added later in Project Settings.')]),
    e('section',{className:'card quick-start',key:'quick'},[e('div',{className:'eyebrow'},'QUICK START'),e('h2',{},'Typical Workflow'),e('ol',{},['Create Project','Add Compound','Draw Structure','Calculate Properties','Add Experimental Data','Run Predictions','Compare / Optimize'].map(item=>e('li',{key:item},item))),e('p',{className:'small'},'Save and calculation remain separate. Prediction and experimental evidence are never merged.')])
   ]),
   (home||globalView==='settings')&&e('section',{className:'card',key:'defaults'},[e('h2',{},'Default Workspace Settings'),e('div',{className:'dashboard-settings'},[
    e('div',{className:'dashboard-setting',key:'type'},[e('span',{},'Default molecule type'),e('strong',{},'Small Molecule')]),e('div',{className:'dashboard-setting',key:'entry'},[e('span',{},'Structure entry'),e('strong',{},'Ketcher or SMILES')]),e('div',{className:'dashboard-setting',key:'calc'},[e('span',{},'Calculation policy'),e('strong',{},'Save first · Calculate on demand')]),e('div',{className:'dashboard-setting',key:'isolation'},[e('span',{},'Data isolation'),e('strong',{},'Project + CompoundVersion')])
   ])]),
   (home||globalView==='settings')&&e(ScientificValidationSection,{key:'sci-val-section'}),
   (home||globalView==='projects')&&e('section',{className:'card',key:'projects'},[
    e('div',{className:'row toolbar'},[e('div',{},[e('div',{className:'eyebrow'},'RESEARCH PORTFOLIO'),e('h2',{},'Projects'),e('p',{className:'small'},'Project cards summarize recorded evidence without synthetic progress percentages.')]),e('div',{className:'row'},[projectSelection.length>0&&e('span',{className:'small'},projectSelection.length+' selected'),e('button',{className:'danger',disabled:projectSelection.length===0,onClick:()=>openDeleteDialog(summaries.filter(item=>projectSelection.includes(item.id)))},'Delete Selected'),projectId&&e('button',{className:'secondary',onClick:()=>openProject(projectId)},'Continue Current Project')])]),
    summaries.length?e('div',{className:'dashboard-project-grid'},summaries.map(item=>e('article',{className:'dashboard-project',key:item.id,tabIndex:0,onClick:()=>openProject(item.id),onKeyDown:event=>{if(event.key==='Enter'||event.key===' ')openProject(item.id)}},[
     e('div',{className:'dashboard-project-actions'},[e('label',{className:'project-select',onClick:event=>event.stopPropagation()},[e('input',{type:'checkbox',checked:projectSelection.includes(item.id),onChange:event=>setProjectSelection(current=>event.target.checked?[...current,item.id]:current.filter(id=>id!==item.id))}),e('span',{},'Select')]),e('button',{className:'danger project-delete-button',onClick:event=>openDeleteDialog([item],event)},'Delete…')]),
     e('div',{className:'dashboard-project-head'},[e('div',{},[e('div',{className:'eyebrow'},item.molecule_type||'Small Molecule'),e('h3',{},item.name)]),e('span',{className:'dashboard-count'},item.compound_count||0)]),
     e('dl',{},[e('div',{key:'target'},[e('dt',{},'Target'),e('dd',{},item.target||'Not set')]),e('div',{key:'experimental'},[e('dt',{},'Experimental records'),e('dd',{},String((item.experimental_activity_count||0)+(item.experimental_admet_count||0)))]),e('div',{key:'optimization'},[e('dt',{},'Optimization runs'),e('dd',{},String(item.optimization_run_count||0))])]),
     e('p',{className:'project-status-summary'},item.status_summary||((item.compound_count||0)+' compounds · experimental and prediction data not started')),e('span',{className:'project-open-link'},'Open Project →')
    ]))):e('div',{className:'empty-state'},[e('h3',{},'No projects yet'),e('p',{},'Use Create New Project above to begin a compound-version-isolated workspace.')])
   ]),
   globalView==='help'&&e('section',{className:'card help-view',key:'help'},[
    e('h2',{},'Platform Workflow'),e('p',{},'Create a project, add compounds, calculate properties, add experimental evidence, run available predictions, then compare or optimize compounds.'),
    e('ol',{},['Create or open a Project','Add a Compound by drawing a structure or entering SMILES','Calculate Properties when ready','Add Experimental Data and run available predictions','Compare compounds or open an Optimization run'].map(item=>e('li',{key:item},item))),
    e('h3',{},'Where scientific functions live'),e('p',{},'Structure, Properties, Activity, ADMET, Metabolism, CYP, Transporters, Safety, and Optimization tools are available inside Project and Compound Detail pages. The global sidebar only changes top-level workspace views.'),
    e('h3',{},'Evidence and support'),e('p',{},'Experimental, Calculated, Predicted, Rule-based, Model unavailable, and Planned states remain distinct. PK / DMPK is planned and is not active in the current stage.'),
    e('div',{className:'row'},[e('button',{onClick:()=>openGlobalView('new-project')},'Create a Project'),e('button',{className:'secondary',onClick:()=>openGlobalView('projects')},'View Projects')])
   ])
  ]);
 }

 function ProjectWorkspace(){
  const summary=selectedProjectSummary;
  const statusByCompound=new Map((summary?.compounds||[]).map(row=>[row.row_id,row]));
  return e(React.Fragment,{},[
   e('div',{className:'card project-header',key:'header'},project?e('div',{className:'row toolbar'},[e('div',{},[e('div',{className:'eyebrow'},'PROJECT DASHBOARD'),e('h1',{},project.name),e('div',{},[e('strong',{},project.target||'Target not set'),' · ',project.molecule_type])]),e('button',{onClick:()=>{setAddCompoundOpen(true);setCompoundForm({compound_id:'',name:'',smiles:'',notes:''})}},'Add Compound')]):e('div',{},[e('h2',{},'Select or create a project'),e('p',{},'Start with a project, then add compounds and work from Compound Detail.') ])),
   project&&e('nav',{className:'project-nav',key:'nav'},[['compounds','Compounds'],['assays','Assays'],['compare','Compare'],['settings','Settings']].map(([tab,label])=>e('button',{key:tab,className:projectTab===tab?'':'secondary',onClick:()=>{setProjectTab(tab);if(tab!=='compounds')setDetail(null)}},label))),
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
   project&&projectTab==='optimization'&&e('div',{className:'card',key:'optimization'},[e('div',{className:'eyebrow'},'PROJECT OPTIMIZATION'),e('h2',{},project.name+' Optimization Workspace'),e('p',{},'Optimization runs are CompoundVersion-specific. This project currently has '+(selectedProjectSummary?.optimization_run_count||0)+' recorded optimization run(s).'),currentVersions.length?e('div',{className:'table-scroll'},e('table',{},[e('thead',{},e('tr',{},[e('th',{},'Compound'),e('th',{},'Version'),e('th',{},'Status'),e('th',{},'')])),e('tbody',{},currentVersions.map(compound=>e('tr',{key:compound.row_id},[e('td',{},compound.name),e('td',{className:'mono'},compound.version?'v'+compound.current_version:'Draft'),e('td',{},StatusBadge({type:compound.status})),e('td',{},e('button',{className:'secondary',disabled:!compound.version,onClick:async()=>{await openDetail(compound.row_id);setProjectTab('compounds');setDetailTab('optimization')}},'Open Optimization'))])))])):e('div',{className:'empty-state'},[e('p',{},'Add a compound before creating an optimization run.'),e('button',{onClick:()=>setProjectTab('compounds')},'Go to Compounds')])]),
   project&&projectTab==='settings'&&e(React.Fragment,{key:'settings'},SettingsPanel())
  ]);
 }

 const sidebarItems=[['Dashboard',()=>goDashboard(),'dashboard'],['New Project',()=>openGlobalView('new-project'),'new-project'],['Projects',()=>openGlobalView('projects'),'projects'],['Optimization',openOptimizationOverview,'optimization'],['Settings',openSettings,'settings'],['Help',()=>openGlobalView('help'),'help']];
 const sidebar=e('aside',{className:'sidebar'+(sidebarOpen?' open':''),key:'sidebar'},[
  e('div',{className:'sidebar-head',key:'head'},[e('div',{},[e('button',{className:'brand-button',onClick:goDashboard},'Drug Optimization Platform'),e('div',{className:'tag'},'Research workspace')]),e('button',{className:'menu-toggle',onClick:()=>setSidebarOpen(value=>!value),'aria-expanded':sidebarOpen,'aria-label':'Toggle primary navigation'},sidebarOpen?'Close':'Menu')]),
  e('div',{className:'sidebar-body',key:'body'},[
   e('nav',{className:'global-nav',key:'nav','aria-label':'Primary navigation'},sidebarItems.map(([label,action,view])=>e('button',{key:label,className:(projectTab==='dashboard'&&globalView===view)||(projectTab===view)?'active':'',onClick:action},label)))
  ])
 ]);
 return e('div',{className:'shell'},[sidebar,e('main',{className:'content',key:'content'},[
  projectTab==='dashboard'?MainDashboard():ProjectWorkspace(),
  AddCompoundPanel(),
  ProjectDeleteModal(),
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
