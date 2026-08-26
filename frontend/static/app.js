(function(){
const e=React.createElement, useState=React.useState, useEffect=React.useEffect;
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
function Svg({src}){return src.startsWith('data:')?e('img',{src,alt:'structure'}):e('span',{className:'structure',dangerouslySetInnerHTML:{__html:src}})}
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

function App(){
 const [projects,setProjects]=useState([]),[projectId,setProjectId]=useState(null),[project,setProject]=useState(null);
 const [form,setForm]=useState({name:'',target:'',indication:'',mechanism_modality:'',description:''});
 const [compoundForm,setCompoundForm]=useState({compound_id:'',name:'',smiles:'',notes:''});
 const [preview,setPreview]=useState(null),[selected,setSelected]=useState([]),[comparison,setComparison]=useState(null),[detail,setDetail]=useState(null),[message,setMessage]=useState('');
 const [projectTab,setProjectTab]=useState('compounds'),[detailTab,setDetailTab]=useState('overview');
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

 const currentVersions=project?.compounds||[];
 const versionLabel=versionId=>{
  const compound=currentVersions.find(item=>item.version.id===Number(versionId));
  return compound?compound.compound_id+' v'+compound.current_version:'Unknown version';
 };
 const endpointName=endpointId=>(admet?.endpoints||[]).find(item=>item.id===endpointId)?.name||endpointId;

 const loadProjects=async()=>{
  const rows=await api.get('/projects');setProjects(rows);
  if(rows.length&&!projectId)setProjectId(rows[0].id);
 };
 const loadProject=async id=>{
  const data=await api.get('/projects/'+id);setProject(data);
  setAdmetVersionId(current=>data.compounds.some(item=>item.version.id===Number(current))?current:(data.compounds[0]?.version.id||''));
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

 useEffect(()=>{loadProjects().catch(error=>setMessage(String(error)))},[]);
 useEffect(()=>{
  setProject(null);setDetail(null);setSelected([]);setComparison(null);setAdmet(null);setMetabolism(null);setAdmetCsvPreview(null);setSelectedSpotId(null);setOptimizationConfig(null);setOptimizationRuns([]);setOptimizationRun(null);setAssays([]);setProposalRuns([]);setProposalRun(null);setSelectedCandidate(null);
  if(projectId)loadProject(projectId).catch(error=>setMessage(String(error)));
 },[projectId]);
 useEffect(()=>{
  if(projectId&&(projectTab==='admet'||(detail&&detailTab==='admet')))Promise.all([loadAdmet(),loadMetabolism()]).catch(error=>setMessage(String(error)));
 },[projectId,projectTab,detailTab,detail?.row_id]);
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

 const createProject=async()=>{
  try{
   const created=await api.post('/projects',form);setForm({name:'',target:'',indication:'',mechanism_modality:'',description:''});
   await loadProjects();setProjectId(created.id);setMessage('Project created');
  }catch(error){setMessage(String(error))}
 };
 const validate=async()=>{try{const result=await api.post('/structure/validate',{smiles:compoundForm.smiles});setPreview(result);setMessage('')}catch(error){setPreview(null);setMessage('Invalid structure: '+error.message)}};
 const saveCompound=async()=>{
  try{
   await api.post('/projects/'+projectId+'/compounds',compoundForm);
   setCompoundForm({compound_id:'C'+String((currentVersions.length+2)).padStart(3,'0'),name:'',smiles:'',notes:''});setPreview(null);
   await loadProject(projectId);await loadProjects();setMessage('Compound saved');
  }catch(error){setMessage(String(error))}
 };
 const openDetail=async rowId=>{
  try{setDetail(await api.get('/compounds/'+rowId+'?include_versions=true'));setDetailTab('overview');setMessage('')}catch(error){setMessage(String(error))}
 };
 const updateStructure=async()=>{
  const smiles=prompt('New SMILES (creates a new version)');if(!smiles)return;
  try{await api.patch('/compounds/'+detail.row_id,{smiles,change_note:'Manual structure edit'});await openDetail(detail.row_id);await loadProject(projectId);setAdmet(null);setMessage('Version created')}catch(error){setMessage(String(error))}
 };
 const compare=async()=>{try{setComparison(await api.get('/projects/'+projectId+'/compare?ids='+selected.join(',')));setMessage('')}catch(error){setComparison(null);setMessage(String(error))}};

 const saveAdmet=async versionId=>{
  const targetVersionId=Number(versionId||admetVersionId);
  if(!targetVersionId)return;
  setAdmetBusy(true);
  try{
   await api.post('/projects/'+projectId+'/admet/measurements',{...admetForm,version_id:targetVersionId});
   setAdmetForm(current=>({...current,value:'',mean:'',sd:'',n:'',notes:''}));await loadAdmet();setMessage('Experimental ADMET saved');
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
  try{const result=await api.post('/admet/predict/'+versionId,{});await loadAdmet();setMessage(result.message)}
  catch(error){setMessage(String(error))}finally{setAdmetBusy(false)}
 };
 const runMetabolism=async versionId=>{
  if(!versionId)return;
  setMetabolismBusy(true);
  try{
   const result=await api.post('/metabolism/predict/'+versionId,{});const data=await loadMetabolism();
   const run=(data?.runs||[]).find(item=>item.version_id===Number(versionId));setSelectedSpotId(run?.spots?.[0]?.id||null);setMessage(result.message);
  }catch(error){setMessage(String(error))}finally{setMetabolismBusy(false)}
 };
 const saveExperimentalMetabolite=async versionId=>{
  setMetabolismBusy(true);
  try{
   await api.post('/projects/'+projectId+'/metabolism/experimental',{...metaboliteForm,version_id:Number(versionId)});
   setMetaboliteForm({...EMPTY_METABOLITE_FORM});await loadMetabolism();setMessage('Experimental metabolite saved');
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
    e('td',{key:'result',className:'num mono'},row.value!=null?(row.qualifier||'=')+' '+row.value:(row.mean!=null?'mean '+row.mean:'-')),
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
   e('summary',{key:'summary'},model.endpoint+': MODEL_UNAVAILABLE'),
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
   e('summary',{key:'summary'},model.endpoint+': MODEL_UNAVAILABLE'),
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
   e('div',{className:'grid',key:'summary'},[group('Strengths',summary.strengths,'strengths'),group('Concerns',summary.concerns,'concerns'),group('Unknown',summary.unknown,'')]),
   e('div',{key:'audit',className:profile.provenance_audit?.status==='PASS'?'pass':'fail'},'Provenance audit: '+profile.provenance_audit?.status+' · '+profile.provenance_audit?.checked+' latest endpoint predictions checked')
  ]);
 }

 function admetPredictionTable(rows){
  if(!rows.length)return Empty({children:'No implemented ADMET predictions yet.'});
  return e('table',{},[
   e('thead',{key:'head'},e('tr',{},['Compound','Endpoint','Experimental','Predicted','Confidence','Domain',''].map(label=>e('th',{key:label},label)))),
   e('tbody',{key:'body'},rows.map(prediction=>{
    const comparison=prediction.experimental_comparisons?.[0];
    const experimental=comparison?(comparison.experimental_value+' '+comparison.experimental_unit):'No compatible value';
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
      e('div',{key:'model'},e('strong',{},'Model evidence: '),(selected.model_evidence?.status||'UNKNOWN')+' — '+(selected.model_evidence?.reason||'')),
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
    e('tbody',{key:'body'},experimental.map(item=>e('tr',{key:item.id},[e('td',{key:'type'},item.label),e('td',{key:'smiles',className:'mono small'},item.canonical_smiles||'Unknown'),e('td',{key:'transform'},item.transformation),e('td',{key:'mass'},item.observed_mass==null?'—':item.observed_mass+' '+item.mass_unit),e('td',{key:'source'},item.source||'—'),e('td',{key:'experiment'},item.experiment||'—'),e('td',{key:'notes'},item.notes||'—')])))
   ]):Empty({children:'No experimental metabolites recorded for this CompoundVersion.'})
  ]);
 }

 function proposalCandidatePanel(candidate){
  if(!candidate)return Empty({children:'Select a candidate to inspect its full rescoring snapshot.'});
  const formatCell=cell=>{
   if(!cell)return '—';const value=cell.value==null?'—':(typeof cell.value==='number'?Number(cell.value).toPrecision(5):String(cell.value));
   return value+(cell.unit?' '+cell.unit:'')+' · '+(cell.type||'Unknown')+(cell.confidence?' · '+cell.confidence:'')+(cell.domain?' · '+(typeof cell.domain==='object'?(cell.domain.classification||'UNKNOWN'):cell.domain):'');
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
   e('p',{key:'activity',className:'small'},activity.status==='COMPLETE'?(Number(activity.value_nm).toPrecision(5)+' nM · '+activity.record_type+' · '+activity.confidence+' · '+activity.applicability_domain+' · nearest '+(activity.nearest_neighbors||[]).slice(0,3).map(row=>row.compound_id+' '+row.similarity).join(', ')):('MODEL_UNAVAILABLE — '+(activity.reason||'No selected assay model'))),
   e('h4',{key:'properties-title'},'Stage 1 property changes'),
   e('div',{key:'properties',className:'small'},Object.entries(propertyDelta).map(([key,value])=>e('span',{key,className:value<0?'delta-down':'delta-up'},key+' '+(value>=0?'+':'')+Number(value).toFixed(3)+' '))),
   e('h4',{key:'soft-title'},'Soft spot changes'),
   e('p',{key:'soft',className:'small'},'Parent primary: '+(soft.parent_primary?.transformation||'Unknown')+' · Candidate primary: '+(soft.candidate_primary?.transformation||'None')+' · parent site absent from candidate Top 3: '+(soft.parent_primary_absent_from_candidate_top3?'YES':'NO')+' · new primary liability: '+(soft.new_primary_liability?'YES':'NO')),
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
   const preferred=row?.preferred;if(!preferred)return 'Unknown';
   const value=preferred.classification??preferred.assessment?.category??preferred.value;
   return String(value??'Unknown')+(preferred.unit?' '+preferred.unit:'')+' · '+preferred.type+' · '+(preferred.confidence||'UNKNOWN')+(preferred.applicability_domain?' · '+preferred.applicability_domain:'');
  };
  const constraintField=(key,label,type='number')=>e('div',{className:'col-3',key},Field({label,type,value:optimizationForm.constraints[key],onChange:value=>setConstraint(key,value)}));
  const regionTable=(title,rows,protectedType)=>e('div',{className:'col-6 optimization-region',key:title},[
   e('h3',{key:'title'},title),
   rows?.length?e('table',{key:'table'},[
    e('thead',{key:'head'},e('tr',{},['Atoms / fragment','Reason','Risk','Confidence','Override'].map(label=>e('th',{key:label},label)))),
    e('tbody',{key:'body'},rows.map(row=>e('tr',{key:row.id},[
     e('td',{key:'atoms',className:'mono small'},row.atom_indices?.length?row.atom_indices.join(', '):(row.fragment||'UNKNOWN')),
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
      e('div',{className:'col-4',key:'activity'},[e('h4',{},'Activity'),e('p',{className:'small'},activity.experimental?'Experimental '+activity.experimental.mean_nm+' nM':(activity.predicted?'Predicted '+activity.predicted.value_nm+' nM · '+activity.predicted.confidence:'Unknown for selected assay'))]),
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
     ]):Empty({children:'No deterministic liability threshold was triggered. Unknown evidence remains visible in Current profile.'})
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
  const detailMeasurements=(admet?.measurements||[]).filter(row=>row.version_id===detail.version.id);
  const detailRuns=(admet?.prediction_runs||[]).filter(run=>run.version_id===detail.version.id);
  const detailPredictions=(admet?.predictions||[]).filter(row=>row.version_id===detail.version.id);
  return e('div',{className:'card'},[
   e('div',{className:'row toolbar',key:'header'},[e('h3',{},detail.compound_id+' · v'+detail.current_version),e('div',{className:'row'},[
    ...['overview','admet','optimization'].map(tab=>e('button',{key:tab,className:detailTab===tab?'':'secondary',onClick:()=>setDetailTab(tab)},tab.toUpperCase())),
    e('button',{key:'modify',className:'secondary',onClick:updateStructure},'Modify structure / new version'),e('button',{key:'close',className:'secondary',onClick:()=>setDetail(null)},'Close')
   ])]),
   detailTab==='overview'&&e('div',{className:'grid',key:'overview'},[
    e('div',{className:'col-4 structure',key:'structure'},Svg({src:detail.version.highlighted_svg})),
    e('div',{className:'col-8',key:'identity'},[e('h4',{},'Identity'),e('div',{className:'mono small'},'Canonical: '+detail.version.canonical_smiles),e('div',{className:'mono small'},'Isomeric: '+detail.version.isomeric_smiles),e('div',{className:'mono small'},'InChIKey: '+detail.version.inchikey)]),
    e('div',{className:'col-4',key:'properties'},propertyTable(detail.version.properties)),e('div',{className:'col-4',key:'rules'},ruleTable(detail.version.rules)),
    e('div',{className:'col-4',key:'assessment'},[e('h4',{},'Structural alerts'),...alertList(detail.version.alerts),e('h4',{},'Assessment'),e('ul',{className:'strengths'},detail.version.assessment.strengths.map(text=>e('li',{key:text},text))),e('ul',{className:'concerns'},detail.version.assessment.concerns.map(text=>e('li',{key:text},text))),e('h4',{},'Provenance'),e('div',{className:'small'},detail.version.provenance.type+' · '+detail.version.provenance.engine+' '+detail.version.provenance.engine_version+' · '+detail.version.provenance.methods.join(', '))]),
    e('div',{className:'col-6',key:'versions'},[e('h4',{},'Version history'),e('table',{},[e('thead',{key:'head'},e('tr',{},['v','SMILES','Change'].map(label=>e('th',{key:label},label)))),e('tbody',{key:'body'},detail.versions.map(version=>e('tr',{key:version.version_number},[e('td',{},'v'+version.version_number),e('td',{className:'mono small'},version.canonical_smiles),e('td',{},version.change_note)])))])]),
    e('div',{className:'col-6',key:'audit'},[e('h4',{},'Prediction audit'),e('table',{},[e('thead',{key:'head'},e('tr',{},['Prediction','Model','Confidence'].map(label=>e('th',{key:label},label)))),e('tbody',{key:'body'},detail.prediction_history.map(run=>e('tr',{key:run.prediction_id},[e('td',{},'#'+run.prediction_id),e('td',{},run.model_name+' '+run.model_version),e('td',{},run.confidence)])))])])
   ]),
   detailTab==='admet'&&e('div',{key:'admet'},[
    e('p',{className:'small',key:'scope'},'ADMET records below are attached to '+detail.compound_id+' v'+detail.current_version+' (compound version #'+detail.version.id+').'),
    integratedProfile(detail.version.id),
    e('h4',{key:'add-title'},'Add experimental measurement'),admetFormPanel(detail.version.id),
    e('h4',{key:'experimental-title',style:{marginTop:'22px'}},'Experimental measurements'),e('div',{key:'experimental-table'},admetMeasurementTable(detailMeasurements)),
    e('div',{className:'row toolbar',key:'prediction-title',style:{marginTop:'22px'}},[e('h4',{},'ADMET predictions'),e('button',{disabled:admetBusy,onClick:()=>runPrediction(detail.version.id)},admetBusy?'Predicting…':'Run prediction')]),
    e('h3',{key:'absorption-title'},'Absorption'),
    e('h4',{key:'stage3a-title'},'Aqueous Solubility / LogS & Caco-2 Papp'),
    e('div',{key:'prediction-table'},admetPredictionTable(detailPredictions.filter(row=>['Solubility','Permeability'].includes(row.endpoint)))),
    e('h3',{key:'distribution-section',style:{marginTop:'24px'}},'Distribution'),
    e('h4',{key:'distribution-title'},'Distribution · Human PPB / fu'),
    e('div',{key:'distribution-table'},admetPredictionTable(detailPredictions.filter(row=>row.endpoint==='Plasma protein binding'))),
    e('h3',{key:'metabolism-section',style:{marginTop:'24px'}},'Metabolism'),
    e('h4',{key:'metabolism-title',style:{marginTop:'18px'}},'Metabolism · HLM / RLM / MLM'),
    e('div',{key:'metabolism-table'},admetPredictionTable(detailPredictions.filter(row=>row.endpoint.endsWith('intrinsic clearance')))),
    e('h4',{key:'cyp-title',style:{marginTop:'18px'}},'Metabolism · CYP'),
    e('div',{key:'cyp-table'},cypPredictionTable(detailPredictions.filter(row=>row.endpoint.startsWith('CYP')))),
    e('div',{key:'cyp-unavailable',className:'small'},(admet?.models||[]).filter(model=>model.endpoint.startsWith('CYP')&&!model.active).map(model=>model.endpoint+': MODEL_UNAVAILABLE — '+model.unavailable_reason).join(' · ')),
    e('div',{key:'metabolic-soft-spots'},metabolismPanel(detail.version.id)),
    e('h3',{key:'transporter-section',style:{marginTop:'24px'}},'Transporters'),
    e('h4',{key:'transporter-title'},'P-gp and available endpoints'),
    e('div',{key:'transporter-table'},transporterPredictionTable(detailPredictions.filter(row=>TRANSPORTER_ENDPOINTS.has(row.endpoint)))),
    e('div',{key:'transporter-unavailable'},unavailableTransporterModels()),
    e('h3',{key:'safety-section',style:{marginTop:'24px'}},'Safety'),
    e('div',{key:'safety-table'},safetyPredictionTable(detailPredictions.filter(row=>SAFETY_ENDPOINTS.has(row.endpoint)))),
    e('div',{key:'safety-unavailable'},unavailableSafetyModels()),
    e('h4',{key:'audit-title',style:{marginTop:'22px'}},'Prediction audit'),
    detailRuns.length?e('table',{key:'runs'},[e('thead',{key:'head'},e('tr',{},['Run','Status','Message','Started'].map(label=>e('th',{key:label},label)))),e('tbody',{key:'body'},detailRuns.map(run=>e('tr',{key:run.id},[e('td',{},'#'+run.id),e('td',{},run.status),e('td',{},run.message),e('td',{},new Date(run.started_at).toLocaleString())]))) ]):Empty({children:'No ADMET prediction runs for this compound version.'})
   ]),
   detailTab==='optimization'&&optimizationPanel(detail.version.id)
  ]);
 }

 return e('div',{className:'shell'},[
  e('aside',{className:'sidebar',key:'sidebar'},[e('h1',{},'AI Drug Optimization Platform'),e('div',{className:'tag'},'Stage 4B · Analog Proposal & Ranking'),
   e('h3',{style:{marginTop:'24px'}},'Projects'),e('ul',{className:'projects'},projects.map(item=>e('li',{key:item.id},e('button',{className:'project-link '+(item.id===projectId?'active':''),onClick:()=>setProjectId(item.id)},item.name,e('div',{className:'tag'},(item.target||'No target')+' · '+item.compound_count+' compounds'))))),
   e('div',{style:{marginTop:'28px'}},[
    ...['name','target','indication','mechanism_modality'].map(key=>e('div',{key,style:{marginBottom:'8px'}},e(Field,{label:key.replace(/_/g,' '),value:form[key],onChange:value=>setForm({...form,[key]:value})}))),
    e(Field,{key:'description',label:'description',value:form.description,onChange:value=>setForm({...form,description:value}),type:'textarea'}),
    e('button',{key:'create',style:{marginTop:'8px'},disabled:!form.name,onClick:createProject},'Create project')
   ])
  ]),
  e('main',{className:'content',key:'content'},[
   e('div',{className:'card row toolbar',key:'project-header'},[e('div',{},[e('h2',{},project?project.name:'Select a project'),e('div',{className:'small'},project?[project.target||'Target not set',project.indication||'',project.mechanism_modality||''].filter(Boolean).join(' · '):'Create a project in the sidebar')]),project&&e('div',{className:'small'},'Updated '+new Date(project.updated_at).toLocaleString())]),
   project&&e('div',{className:'card',key:'add-compound'},[e('h3',{},'Add compound'),e('div',{className:'grid'},[
    e('div',{className:'col-3',key:'id'},Field({label:'Compound ID',value:compoundForm.compound_id,onChange:value=>setCompoundForm({...compoundForm,compound_id:value})})),e('div',{className:'col-3',key:'name'},Field({label:'Name',value:compoundForm.name,onChange:value=>setCompoundForm({...compoundForm,name:value})})),
    e('div',{className:'col-6',key:'smiles'},Field({label:'SMILES / editor output',value:compoundForm.smiles,onChange:value=>{setCompoundForm({...compoundForm,smiles:value});setPreview(null)}})),e('div',{className:'col-6',key:'notes'},Field({label:'Notes',value:compoundForm.notes,onChange:value=>setCompoundForm({...compoundForm,notes:value})})),
    e('div',{className:'col-6 row',key:'actions'},[e('button',{className:'secondary',disabled:!compoundForm.smiles,onClick:validate},'Validate & calculate'),e('button',{disabled:!preview||!compoundForm.compound_id,onClick:saveCompound},'Save compound')])
   ])]),
   preview&&e('div',{className:'card',key:'preview'},[e('h3',{},'Live validation'),e('div',{className:'grid col-12'},[e('div',{className:'col-3 structure'},Svg({src:preview.svg})),e('div',{className:'col-9'},[e('div',{className:'mono small'},preview.identity.canonical_smiles),e('table',{},e('tbody',{},[['Formula','molecular_formula'],['MW','molecular_weight'],['cLogP','clogp'],['TPSA','tpsa'],['HBD','hbd'],['HBA','hba'],['QED','qed']].map(([label,key])=>e('tr',{key},[e('td',{},label),e('td',{className:'num mono'},String(preview.properties[key]))]))))]),e('div',{className:'col-12'},[e('strong',{},'Provenance: '),preview.provenance.engine+' '+preview.provenance.engine_version+' · Calculated · '+preview.provenance.methods.join(', ')])])]),
   project&&e('div',{className:'card row toolbar',key:'tabs'},['compounds','admet'].map(tab=>e('button',{key:tab,className:projectTab===tab?'':'secondary',onClick:()=>setProjectTab(tab)},tab==='compounds'?'Compounds / SAR':'ADMET'))),
   project&&projectTab==='compounds'&&e(React.Fragment,{key:'compounds'},[
    e('div',{className:'card',key:'compound-list'},[e('h3',{},'Stage 2 · Assays / Activity / Models'),e('div',{className:'row toolbar'},[e('a',{className:'button secondary',href:'/static/stage2-workbench.html?project='+projectId},'Open workbench')]),e('div',{className:'row toolbar'},[e('h3',{},'Compounds ('+currentVersions.length+')'),e('div',{className:'row'},[e('span',{className:'small'},'Select at least two'),e('button',{className:'secondary',disabled:selected.length<2,onClick:compare},'Compare selected')])]),
     e('table',{},[e('thead',{key:'head'},e('tr',{},['','ID','Name','Version','Canonical SMILES','MW','cLogP','TPSA','QED',''].map((label,index)=>e('th',{key:index,className:['MW','cLogP','TPSA','QED'].includes(label)?'num':''},label)))),e('tbody',{key:'body'},currentVersions.map(compound=>e('tr',{key:compound.row_id},[
      e('td',{key:'select'},e('input',{type:'checkbox',checked:selected.includes(compound.row_id),onChange:event=>setSelected(event.target.checked?[...selected,compound.row_id]:selected.filter(id=>id!==compound.row_id))})),e('td',{key:'id',className:'mono'},compound.compound_id),e('td',{key:'name'},compound.name||'-'),e('td',{key:'version'},'v'+compound.current_version),e('td',{key:'smiles',className:'mono small'},compound.version.canonical_smiles),...['molecular_weight','clogp','tpsa','qed'].map(key=>e('td',{key,className:'num mono'},compound.version.properties[key]??'-')),e('td',{key:'open'},e('button',{className:'secondary',onClick:()=>openDetail(compound.row_id)},'Open'))
     ])))])]),
    comparison&&e('div',{className:'card',key:'comparison'},[e('h3',{},'Comparison'),e('p',{className:'small'},'Activity and ADMET columns are optional. Compatible experimental ADMET values take display priority; predictions remain stored. No overall score or ranking is calculated.'),e('table',{},[e('thead',{key:'head'},e('tr',{},['Compound',...comparison.metrics].map((label,index)=>e('th',{key:index,className:label!=='Compound'?'num':'',title:comparison.metric_units?.[label]||''},label)))),e('tbody',{key:'body'},comparison.compounds.map(compound=>e('tr',{key:compound.row_id},[e('td',{key:'compound',className:'mono'},compound.compound),...comparison.metrics.map(metric=>e('td',{key:metric,className:'num mono',title:comparison.metric_units?.[metric]||''},compound[metric]??'-'))])))]),e('div',{className:'plot-grid',style:{marginTop:'15px'}},pairPlot(comparison,'MW','cLogP'),pairPlot(comparison,'TPSA','cLogP'),qedHistogram(comparison))]),
    compoundDetail()
   ]),
   project&&projectTab==='admet'&&e(React.Fragment,{key:'admet'},projectAdmetTab()),
   message&&e('pre',{className:'card error',key:'message'},message)
  ])
 ]);

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
