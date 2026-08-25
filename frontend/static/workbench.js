(function () {
  const e = React.createElement;
  const projectId = new URLSearchParams(location.search).get('project');

  async function req(path, options = {}) {
    const response = await fetch('/api' + path, { headers: { 'Content-Type': 'application/json' }, ...options });
    const text = await response.text();
    let data;
    try { data = text ? JSON.parse(text) : null; } catch (_) { data = text; }
    if (!response.ok) throw new Error(typeof data === 'object' ? JSON.stringify(data, null, 2) : (data || response.statusText));
    return data;
  }

  function Workbench() {
    const [project, setProject] = React.useState(null);
    const [assays, setAssays] = React.useState([]);
    const [assayId, setAssayId] = React.useState('');
    const [form, setForm] = React.useState({ name: 'New IC50 assay', measurement_type: 'IC50', unit: 'nM', species: '', cell_line: '', mutation_variant: '', protocol: '', reference_compound: '' });
    const [model, setModel] = React.useState(null);
    const [cliffs, setCliffs] = React.useState(null);
    const [mmp, setMmp] = React.useState(null);
    const [message, setMessage] = React.useState('');
    const [sar, setSar] = React.useState(null);

    React.useEffect(() => {
      async function load() {
        try {
          const projectData = await req('/projects/' + projectId);
          setProject(projectData);
          const assayData = await req('/projects/' + projectId + '/assays');
          setAssays(assayData);
          if (assayData.length) setAssayId(String(assayData[0].id));
        } catch (error) { setMessage(String(error)); }
      }
      load();
    }, []);

    async function saveAssay() {
      try {
        await req('/projects/' + projectId + '/assays', { method: 'POST', body: JSON.stringify(form) });
        const assayData = await req('/projects/' + projectId + '/assays');
        setAssays(assayData);
        if (assayData.length) setAssayId(String(assayData[assayData.length - 1].id));
        setMessage('Assay saved');
      } catch (error) { setMessage(String(error)); }
    }

    async function addMeasurement(compound) {
      try {
        const result = await req('/assays/' + assayId + '/measurements', {
          method: 'POST',
          body: JSON.stringify({
            version_id: compound.row_id,
            value: document.getElementById('value-' + compound.row_id).value,
            unit: document.getElementById('unit-' + compound.row_id).value
          })
        });
        setMessage(JSON.stringify(result, null, 2));
      } catch (error) { setMessage(String(error)); }
    }

    async function trainModel() {
      try {
        const result = await req('/assays/' + assayId + '/models/train', { method: 'POST' });
        setModel(result); setMessage(JSON.stringify(result, null, 2));
      } catch (error) { setMessage(String(error)); }
    }

    async function predictAll() {
      try {
        for (const compound of project.compounds) {
          try { await req('/assays/' + assayId + '/predict/' + compound.row_id, { method: 'POST' }); } catch (_) {}
        }
        setSar(await req('/projects/' + projectId + '/sar?assay_id=' + assayId));
      } catch (error) { setMessage(String(error)); }
    }

    async function loadCliffs() {
      try { setCliffs(await req('/projects/' + projectId + '/cliffs?assay_id=' + assayId)); } catch (error) { setMessage(String(error)); }
    }

    async function loadMmp() {
      try { setMmp(await req('/projects/' + projectId + '/mmp?assay_id=' + assayId)); } catch (error) { setMessage(String(error)); }
    }

    const controls = e('div', { className: 'card row toolbar' }, [
      e('a', { href: '/' }, 'Back to projects'),
      e('h1', {}, 'Stage 2 Workbench'),
      e('select', { value: assayId, onChange: event => setAssayId(event.target.value) },
        assays.map(assay => e('option', { key: assay.id, value: assay.id }, assay.name + ' #' + assay.id))),
      e('button', { onClick: trainModel }, 'Train model'),
      e('button', { className: 'secondary', onClick: predictAll }, 'Build SAR + predict all'),
      e('button', { className: 'secondary', onClick: loadCliffs }, 'Activity cliffs'),
      e('button', { className: 'secondary', onClick: loadMmp }, 'MMP'),
      e('a', { className: 'button', href: '/api/projects/' + projectId + '/sar-export.csv?assay_id=' + assayId, target: '_blank' }, 'Export CSV')
    ]);

    const children = [
      controls,
      message && e('pre', { className: 'card error' }, message),
      e('div', { className: 'card' }, [
        e('h3', {}, 'Create assay'),
        Object.entries(form).map(([key, value]) => e('div', { key }, [
          e('label', {}, key.replace(/_/g, ' ')),
          e(key === 'protocol' ? 'textarea' : 'input', { value, onChange: event => setForm({ ...form, [key]: event.target.value }) })
        ])),
        e('button', { onClick: saveAssay }, 'Save assay')
      ]),
      e('div', { className: 'card' }, [
        e('h3', {}, 'Experimental input'),
        ...(project?.compounds || []).map(compound => e('div', { key: compound.row_id, className: 'row' }, [
          e('strong', { style: { minWidth: '70px' } }, compound.compound_id),
          e('input', { id: 'value-' + compound.row_id, placeholder: 'value', style: { maxWidth: '120px' } }),
          e('input', { id: 'unit-' + compound.row_id, defaultValue: 'nM', style: { maxWidth: '80px' } }),
          e('button', { onClick: () => addMeasurement(compound) }, 'Add replicate')
        ]))
      ]),
      sar && e('div', { className: 'card' }, [
        e('h3', {}, 'SAR table'),
        e('table', {}, [
          e('thead', {}, e('tr', {}, ['ID', 'N', 'Mean nM', 'pActivity', 'Confidence', 'Domain', 'MW', 'cLogP'].map(label => e('th', { key: label }, label)))),
          e('tbody', {}, sar.compounds.map(compound => e('tr', { key: compound.row_id }, [
            compound.compound_id,
            compound.experimental?.n ?? '-',
            compound.experimental?.mean_nm ?? '-',
            compound.experimental?.pactivity_mean ?? '-',
            compound.predicted?.confidence ?? '-',
            compound.predicted?.applicability_domain ?? '-',
            compound.properties.molecular_weight,
            compound.properties.clogp
          ].map((value, index) => e('td', { key: index }, String(value))))))
        ])
      ]),
      model && e('pre', { className: 'card' }, JSON.stringify(model, null, 2)),
      cliffs && e('pre', { className: 'card' }, JSON.stringify(cliffs, null, 2)),
      mmp && e('pre', { className: 'card' }, JSON.stringify(mmp, null, 2))
    ];

    return e(React.Fragment, null, ...children);
  }

  ReactDOM.createRoot(document.getElementById('root')).render(e(Workbench));
})();
