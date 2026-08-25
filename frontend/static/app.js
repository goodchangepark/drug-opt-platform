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

 useEffect(()=>{loadProjects().catch(error=>setMessage(String(error)))},[]);
 useEffect(()=>{
  setProject(null);setDetail(null);setSelected([]);setComparison(null);setAdmet(null);setMetabolism(null);setAdmetCsvPreview(null);setSelectedSpotId(null);
  if(projectId)loadProject(projectId).catch(error=>setMessage(String(error)));
 },[projectId]);
 useEffect(()=>{
  if(projectId&&(projectTab==='admet'||(detail&&detailTab==='admet')))Promise.all([loadAdmet(),loadMetabolism()]).catch(error=>setMessage(String(error)));
 },[projectId,projectTab,detailTab,detail?.row_id]);

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
    output.liability_summary&&e('div',{key:'cyp-liability'},[e('strong',{},'CYP liability rule: '),output.liability_summary.flag+' · '+output.liability_summary.rule+' · '+output.liability_summary.basis]),
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
    e('div',{className:'row toolbar',key:'title'},[e('h3',{},'ADMET predictions through Stage 3C'),e('button',{disabled:admetBusy||!admetVersionId,onClick:()=>runPrediction(Number(admetVersionId))},admetBusy?'Predicting…':'Run prediction')]),
    e('p',{key:'scope',className:'small'},'CYP inhibitor and substrate classifiers are isolated endpoints. Probabilities are binary class probabilities and are never converted to IC50. CYP1A2/CYP2C19 substrate models remain unavailable. Stage 3A/3B endpoints retain their definitions.'),
    admetPredictionTable((admet?.predictions||[]).filter(row=>!row.endpoint.startsWith('CYP'))),
    e('h4',{key:'cyp-predictions-title',style:{marginTop:'22px'}},'CYP inhibitor / substrate predictions'),
    cypPredictionTable((admet?.predictions||[]).filter(row=>row.endpoint.startsWith('CYP'))),
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

 function compoundDetail(){
  if(!detail)return null;
  const detailMeasurements=(admet?.measurements||[]).filter(row=>row.version_id===detail.version.id);
  const detailRuns=(admet?.prediction_runs||[]).filter(run=>run.version_id===detail.version.id);
  const detailPredictions=(admet?.predictions||[]).filter(row=>row.version_id===detail.version.id);
  return e('div',{className:'card'},[
   e('div',{className:'row toolbar',key:'header'},[e('h3',{},detail.compound_id+' · v'+detail.current_version),e('div',{className:'row'},[
    ...['overview','admet'].map(tab=>e('button',{key:tab,className:detailTab===tab?'':'secondary',onClick:()=>setDetailTab(tab)},tab.toUpperCase())),
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
    e('h4',{key:'add-title'},'Add experimental measurement'),admetFormPanel(detail.version.id),
    e('h4',{key:'experimental-title',style:{marginTop:'22px'}},'Experimental measurements'),e('div',{key:'experimental-table'},admetMeasurementTable(detailMeasurements)),
    e('div',{className:'row toolbar',key:'prediction-title',style:{marginTop:'22px'}},[e('h4',{},'ADMET predictions'),e('button',{disabled:admetBusy,onClick:()=>runPrediction(detail.version.id)},admetBusy?'Predicting…':'Run prediction')]),
    e('h4',{key:'distribution-title'},'Distribution · Human PPB / fu'),
    e('div',{key:'distribution-table'},admetPredictionTable(detailPredictions.filter(row=>row.endpoint==='Plasma protein binding'))),
    e('h4',{key:'metabolism-title',style:{marginTop:'18px'}},'Metabolism · HLM / RLM / MLM'),
    e('div',{key:'metabolism-table'},admetPredictionTable(detailPredictions.filter(row=>row.endpoint.endsWith('intrinsic clearance')))),
    e('h4',{key:'cyp-title',style:{marginTop:'18px'}},'Metabolism · CYP'),
    e('div',{key:'cyp-table'},cypPredictionTable(detailPredictions.filter(row=>row.endpoint.startsWith('CYP')))),
    e('div',{key:'cyp-unavailable',className:'small'},(admet?.models||[]).filter(model=>model.endpoint.startsWith('CYP')&&!model.active).map(model=>model.endpoint+': MODEL_UNAVAILABLE — '+model.unavailable_reason).join(' · ')),
    e('div',{key:'metabolic-soft-spots'},metabolismPanel(detail.version.id)),
    e('h4',{key:'stage3a-title',style:{marginTop:'18px'}},'Solubility & Caco-2'),
    e('div',{key:'prediction-table'},admetPredictionTable(detailPredictions.filter(row=>['Solubility','Permeability'].includes(row.endpoint)))),
    e('h4',{key:'audit-title',style:{marginTop:'22px'}},'Prediction audit'),
    detailRuns.length?e('table',{key:'runs'},[e('thead',{key:'head'},e('tr',{},['Run','Status','Message','Started'].map(label=>e('th',{key:label},label)))),e('tbody',{key:'body'},detailRuns.map(run=>e('tr',{key:run.id},[e('td',{},'#'+run.id),e('td',{},run.status),e('td',{},run.message),e('td',{},new Date(run.started_at).toLocaleString())]))) ]):Empty({children:'No ADMET prediction runs for this compound version.'})
   ])
  ]);
 }

 return e('div',{className:'shell'},[
  e('aside',{className:'sidebar',key:'sidebar'},[e('h1',{},'AI Drug Optimization Platform'),e('div',{className:'tag'},'Stage 3D · Metabolic Hypotheses'),
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
