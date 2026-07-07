import * as FilePond from 'filepond'
import Choices from 'choices.js'
import 'filepond/dist/filepond.min.css'
import 'choices.js/public/assets/styles/choices.min.css'
import dataset from '../datasets/1kg-ont.json'

const el = (id) => document.getElementById(id)

const worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' })

const errorEl = el('wally-error')
const errorMsgEl = el('error-message')
const logEl = el('wally-log')
const resultsEl = el('results-container')

$('#mainTab a, #source-tabs a').on('click', function (e) {
  e.preventDefault()
  $(this).tab('show')
})

function showError(msg) {
  errorMsgEl.textContent = msg
  errorEl.classList.toggle('d-none', !msg)
}
function appendLog(line) { logEl.textContent += line + '\n'; logEl.scrollTop = logEl.scrollHeight }
function showResults(show) { resultsEl.classList.toggle('d-none', !show) }
function activateTab(href) { $(`#mainTab a[href="${href}"]`).tab('show') }
function busy(on) { el('btn-submit').disabled = on }

const fpOpts = { allowMultiple: true, allowProcess: false, instantUpload: false, credits: false }
const pondAln = FilePond.create(el('localAln'), {
  ...fpOpts,
  labelIdle: 'Drag &amp; drop BAM/CRAM + indexes, or <span class="filepond--label-action">browse</span>'
})
const pondRef = FilePond.create(el('localRef'), {
  ...fpOpts,
  labelIdle: 'Drag &amp; drop reference + index, or <span class="filepond--label-action">browse</span>'
})

// Test hook
window.__wally = { pondAln, pondRef, sampleChoices: null }

const rx = {
  alignment: /\.(bam|cram)$/i,
  alnIndex: /\.(bai|crai|csi)$/i,
  reference: /\.(fa|fasta|fna)(\.gz)?$/i
}

// Index files
function hasIndexFor(name, allNames) {
  const base = name.replace(rx.alignment, '')
  return allNames.some((n) => rx.alnIndex.test(n) && (n.startsWith(name) || n.startsWith(base + '.')))
}

function collectLocal() {
  const alnFiles = pondAln.getFiles().map((f) => f.file)
  const refFiles = pondRef.getFiles().map((f) => f.file)
  if (!alnFiles.length) throw new Error('Drop at least one BAM/CRAM (with its index) in the alignments field.')
  if (!refFiles.length) throw new Error('Drop a reference FASTA (with its .fai) in the reference field.')

  const alnNames = alnFiles.map((f) => f.name)
  const alignments = alnFiles.filter((f) => rx.alignment.test(f.name))
  if (!alignments.length) throw new Error('No .bam/.cram found in the alignments field.')
  for (const a of alignments) {
    if (!hasIndexFor(a.name, alnNames)) throw new Error(`Missing index for ${a.name} - add its .bai/.crai file.`)
  }

  const references = refFiles.filter((f) => rx.reference.test(f.name) && !rx.alignment.test(f.name))
  if (references.length === 0) throw new Error('No reference FASTA (.fa/.fa.gz) in the reference field.')
  if (references.length > 1) throw new Error('Multiple reference FASTA files - please provide exactly one.')
  const ref = references[0]
  const refNames = refFiles.map((f) => f.name)
  if (!refNames.some((n) => n === ref.name + '.fai' || n === ref.name.replace(rx.reference, '') + '.fai')) {
    throw new Error(`Missing ${ref.name}.fai - add the reference index (samtools faidx).`)
  }
  if (/\.gz$/i.test(ref.name) && !refNames.some((n) => n.endsWith('.gzi'))) {
    throw new Error(`Missing ${ref.name}.gzi - a bgzipped reference needs its .gzi index.`)
  }

  return {
    mode: 'local',
    files: [...alnFiles, ...refFiles],
    bamNames: alignments.map((f) => f.name),
    refName: ref.name
  }
}

// 1000 Genomes

let sampleChoices = null
try {
  const samples = dataset.builds.hg38.samples 
  sampleChoices = new Choices(el('remoteSample'), {
    removeItemButton: true,
    searchEnabled: true,
    shouldSort: false,
    itemSelectText: '',
    searchResultLimit: 50,
    placeholderValue: 'Search and add samples…',
    choices: samples.map((s) => ({ value: s, label: s }))
  })
  window.__wally.sampleChoices = sampleChoices
} catch (e) {
  showError('Could not build the 1000 Genomes sample list: ' + e.message)
}

function collectRemote() {
  const build = el('remoteBuild').value
  const selected = sampleChoices ? sampleChoices.getValue(true) : []
  const samples = Array.isArray(selected) ? selected : (selected ? [selected] : [])
  if (!samples.length) throw new Error('Please select at least one sample.')
  const b = dataset.builds[build]
  return {
    mode: 'remote',
    refName: b.reference.split('/').pop(),
    refUrl: b.reference,
    samples: samples.map((s) => {
      const cramName = s + b.cramSuffix
      const cramUrl = `${b.cramBase}/${cramName}`
      return { name: s, cramName, cramUrl, craiUrl: cramUrl + '.crai' }
    })
  }
}

// Render

function currentSource() {
  return el('source-remote-tab').classList.contains('active') ? 'remote' : 'local'
}

function parseRegion(s) {
  const m = String(s).trim().match(/^([^:\s]+):([\d,]+)-([\d,]+)$/)
  if (!m) return null
  const start = parseInt(m[2].replace(/,/g, ''), 10)
  const end = parseInt(m[3].replace(/,/g, ''), 10)
  if (!(end > start)) return null
  return { chr: m[1], start, end }
}
const fmtRegion = (v) => `${v.chr}:${v.start}-${v.end}`

function clampView(v) {
  if (v.start < 1) { v.end += 1 - v.start; v.start = 1 }
  if (v.end - v.start < 20) v.end = v.start + 20
  return v
}

let currentView = null

// Render 
function doRender(view) {
  showError('')
  if (!view) { showError('Invalid region. Use e.g. chr1:1000000-1000200'); activateTab('#result-tab'); return }
  let source
  try { source = currentSource() === 'remote' ? collectRemote() : collectLocal() } catch (err) { showError(err.message); activateTab('#result-tab'); return }

  currentView = view
  const region = fmtRegion(view)
  el('region').value = region
  el('regionNav').value = region

  logEl.textContent = ''
  el('results-toolbar').classList.remove('d-none') 
  showResults(true)
  activateTab('#result-tab')
  busy(true)

  worker.postMessage({
    type: 'render',
    ...source,
    region,
    width: 1024,
    height: 1024,
    paired: el('paired').checked,
    clip: el('clip').checked,
    supplementary: el('supp').checked,
    coverage: el('coverage').checked,
    mod: parseInt(el('modType').value, 10) || 0
  })
}

function zoom(factor) {
  if (!currentView) return
  const center = (currentView.start + currentView.end) / 2
  const span = Math.max(20, Math.round((currentView.end - currentView.start) * factor))
  doRender(clampView({ chr: currentView.chr, start: Math.round(center - span / 2), end: Math.round(center + span / 2) }))
}
function shift(frac) {
  if (!currentView) return
  const d = Math.round((currentView.end - currentView.start) * frac)
  doRender(clampView({ chr: currentView.chr, start: currentView.start + d, end: currentView.end + d }))
}

worker.onmessage = (ev) => {
  const msg = ev.data
  if (msg.type === 'log') {
    appendLog(msg.line)
  } else if (msg.type === 'result') {
    const blob = new Blob([msg.png], { type: 'image/png' })
    const url = URL.createObjectURL(blob)
    const img = new Image()
    img.onload = () => {
      const c = el('canvas')
      c.width = img.width
      c.height = img.height
      c.getContext('2d').drawImage(img, 0, 0)
      c.style.transform = ''
      URL.revokeObjectURL(url)
    }
    img.src = url
    showResults(true)
    appendLog(`Rendered in ${msg.elapsed} ms${msg.prefetched ? ' (prefetched window)' : msg.xhr ? ` (${msg.xhr} requests)` : ' (cached)'}.`)
    busy(false)
  } else if (msg.type === 'error') {
    showError(msg.message)
    busy(false)
  }
}

// Initial render
el('btn-submit').addEventListener('click', () => doRender(parseRegion(el('region').value)))

// Results-tab navigation
el('nav-zoomin').addEventListener('click', () => zoom(0.5))
el('nav-zoomout').addEventListener('click', () => zoom(2))
el('nav-left').addEventListener('click', () => shift(-0.2))
el('nav-right').addEventListener('click', () => shift(0.2))
el('nav-go').addEventListener('click', () => doRender(parseRegion(el('regionNav').value)))
el('regionNav').addEventListener('keydown', (e) => { if (e.key === 'Enter') { e.preventDefault(); doRender(parseRegion(el('regionNav').value)) } })
for (const id of ['paired', 'clip', 'supp', 'coverage', 'modType']) {
  el(id).addEventListener('change', () => { if (currentView) doRender(currentView) })
}

// Mouse panning
{
  const canvas = el('canvas')
  let dragging = false
  let startX = 0
  let moved = 0
  canvas.addEventListener('mousedown', (e) => {
    if (!currentView) return
    dragging = true; startX = e.clientX; moved = 0
    canvas.style.cursor = 'grabbing'
    e.preventDefault()
  })
  window.addEventListener('mousemove', (e) => {
    if (!dragging) return
    moved = e.clientX - startX
    canvas.style.transform = `translateX(${moved}px)`
  })
  window.addEventListener('mouseup', () => {
    if (!dragging) return
    dragging = false
    canvas.style.cursor = ''
    if (Math.abs(moved) < 3 || !currentView) { canvas.style.transform = ''; return }
    const frac = moved / (canvas.clientWidth || canvas.width)
    const span = currentView.end - currentView.start
    const d = Math.round(-frac * span)
    doRender(clampView({ chr: currentView.chr, start: currentView.start + d, end: currentView.end + d }))
  })
}

el('btn-example').addEventListener('click', () => {
  showError('')
  $('#source-remote-tab').tab('show')
  el('remoteBuild').value = 'hg38'
  if (sampleChoices) {
    sampleChoices.removeActiveItems()
    sampleChoices.setChoiceByValue(['NA19900', 'NA19720', 'HG02266'])
  }
  el('region').value = 'chr8:33904080-33904200'
})
