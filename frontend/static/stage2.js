(function(){
const e=React.createElement;
const api={async req(path,opt={}){const r=await fetch('/api'+path,{headers:{'Content-Type':'application/json'},...opt});const t=await r.text();let d;try{d=t?JSON.parse(t):null}catch(_){d=t}if(!r.ok)throw new Error(typeof d==='object'?JSON.stringify(d,null,2):d||r.statusText);return d},post:(p,v)=>api.req(p,{method:'POST',body:JSON.stringify(v)}),get:p=>api.req(p)};
window.Stage2={
 async createAssay(projectId,form){return api.post('/projects/'+projectId+'/assays',form)},
 async listAssays(projectId){return api.get('/projects/'+projectId+'/assays')},
 async addMeasurement(assayId,row){return api.post('/assays/'+assayId+'/measurements',row)},
 async trainModel(assayId){return api.post('/assays/'+assayId+'/models/train',{})},
 async predict(assayId,rowId){return api.post('/assays/'+assayId+'/predict/'+rowId,{})},
 async sarTable(projectId,assayId){return api.get('/projects/'+projectId+'/sar?assay_id='+assayId)},
 async cliffs(projectId,assayId){return api.get('/projects/'+projectId+'/cliffs?assay_id='+assayId)},
 exportUrl(projectId,assayId){return '/api/projects/'+projectId+'/sar-export.csv?assay_id='+assayId},
 AssayPanel({project,onChange}){
  const [form,setForm]=React.useState({name:'EGFR Ba/F3 IC50',target:project.target||'',target_type:'Kinase',assay_category:'Cellular proliferation',measurement_type:'IC50',custom_measurement_name:'',unit:'nM',species:'mouse',cell_line:'Ba/F3',mutation_variant:'EGFR D770_N771insSVD',protein_construct:'',substrate:'',atp_concentration:'',incubation_time:'72 h',detection_method:'CellTiter-Glo',experimental_conditions:'',protocol:'',reference_compound:'',reference_structure_smiles:'',reference_activity:'',reference_source:'',notes:''});
  const [assays,setAssays]=React.useState([]),[message,setMessage]=React.useState('');
  const load=()=>api.get('/projects/'+project.id+'/assays').then(setAssays).catch(String);
  React.useEffect(()=>{load()},[project.id]);
  const save=async()=>{try{await this.createAssay(project.id,form);load();setMessage('Assay saved')}catch(x){setMessage(String(x))}};
  return e('div',{className:'card'},[e('h3',{},'Assay settings'),e('div',{className:'grid'},Object.entries(form).map(([k,v])=>e('div',{key:k,className:k.includes('conditions')||['protocol','notes'].includes(k)?'col-6':'col-3'},e('label',{},k.replace(/_/g,' ')),e(k==='protocol'||k==='notes'?'textarea':'input',{value:v,onChange:x=>setForm({...form,[k]:x.target.value})})))),
   e('div',{className:'row'},e('button',{disabled:!form.name,onClick:save},'Save assay'),e('button',{className:'secondary',onClick:load},'Reload'),e('span',{className:'small'},message)),e('table',{},e('thead',{},e('tr',{},['Name','Type','Unit','Species','Cell line','Mutation'])),e('tbody',{},assays.map(a=>e('tr',{key:a.id},[a.name,a.measurement_type,a.unit,a.species,a.cell_line,a.mutation_variant].map((v,i)=>e('td',{key:i},v))))))]);
 },
 ActivityPanel({project,compounds}){
  const [assays,setAssays]=React.useState([]),[assayId,setAssayId]=React.useState(''),[rows,setRows]=React.useState([]),[msg,setMsg]=React.useState('');
  React.useEffect(()=>{api.get('/projects/'+project.id+'/assays').then(a=>{setAssays(a);if(a.length)setAssayId(a[0].id)})},[]);
  const submit=async row=>{try{const r=await api.post('/assays/'+assayId+'/measurements',row);setRows([...rows,r]);setMsg('Experimental measurement saved')}catch(x){setMsg(String(x))}};
  return e('div',{className:'card'},[e('h3',{},'Activity input'),
   e('select',{value:assayId,onChange:x=>setAssayId(x.target.value)},assays.map(a=>e('option',{key:a.id,value:a.id},a.name+' · '+a.measurement_type))),
   ...compounds.map(c=>e('div',{key:c.row_id,className:'row',style:{margin:'7px 0'}},[
    e('strong',{style:{minWidth:'70px'}},c.compound_id),e('input',{placeholder:'value',id:'val-'+c.row_id,style:{maxWidth:'110px'}}),e('input',{defaultValue:c.version.properties.unit||'nM',id:'unit-'+c.row_id,style:{maxWidth:'80px'}}),e('button',{onClick:()=>submit({version_id:c.row_id,value:document.getElementById('val-'+c.row_id).value,unit:document.getElementById('unit-'+c.row_id).value})},'Add replicate')])),e('div',{className:'small'},msg)]);
  }
};
})();
