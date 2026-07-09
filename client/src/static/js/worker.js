// Web Worker for wally WebAssembly module

import createWally from './wally.js'

let xhrCount = 0
{
  const _get = XMLHttpRequest.prototype.getResponseHeader
  XMLHttpRequest.prototype.getResponseHeader = function (name) {
    const v = _get.call(this, name)
    if (v == null && /^accept-ranges$/i.test(name)) return 'bytes'
    return v
  }
  const _open = XMLHttpRequest.prototype.open
  XMLHttpRequest.prototype.open = function (...a) { xhrCount++; return _open.apply(this, a) }
}

let ready = null
function init() {
  if (!ready) {
    ready = createWally({
      print: (t) => postMessage({ type: 'log', line: t }),
      printErr: (t) => postMessage({ type: 'log', line: t })
    })
  }
  return ready
}

// Emscripten helper
function freshDir(M, path) {
  try {
    for (const name of M.FS.readdir(path)) {
      if (name === '.' || name === '..') continue
	try {
	    M.FS.unlink(path + '/' + name)
	} catch (e) { }
    }
    M.FS.rmdir(path)
  } catch (e) { }
  M.FS.mkdir(path)
}

function mountLocal(M, msg) {
  try {
    M.FS.unmount('/data')
  } catch (e) { }
  try {
    M.FS.rmdir('/data')
  } catch (e) { }
  M.FS.mkdir('/data')
  M.FS.mount(M.FS.filesystems.WORKERFS, { files: msg.files }, '/data')
  return {
    bams: msg.bamNames.map((n) => '/data/' + n).join('\n'),
    genome: '/data/' + msg.refName
  }
}

let remoteKey = null

function mountRemote(M, msg) {
  const key = msg.samples.map((s) => s.cramUrl).join(',') + '|' + msg.refUrl + '|' + (msg.annoUrl || '')
  if (key !== remoteKey) {
    freshDir(M, '/remote')
    freshDir(M, '/ref')
    for (const s of msg.samples) {
      M.FS.createLazyFile('/remote', s.cramName, s.cramUrl, true, false)
      M.FS.createLazyFile('/remote', s.cramName + '.crai', s.craiUrl, true, false)
    }
    M.FS.createLazyFile('/ref', msg.refName, msg.refUrl, true, false)
    M.FS.createLazyFile('/ref', msg.refName + '.fai', msg.refUrl + '.fai', true, false)
    M.FS.createLazyFile('/ref', msg.refName + '.gzi', msg.refUrl + '.gzi', true, false)
    if (msg.annoUrl && msg.annoName) {
      M.FS.createLazyFile('/ref', msg.annoName, msg.annoUrl, true, false)
      M.FS.createLazyFile('/ref', msg.annoName + '.tbi', msg.annoUrl + '.tbi', true, false)
    }
    remoteKey = key
  }
}

function parseRegion(s) {
  const m = String(s).match(/^([^:]+):(\d+)-(\d+)$/)
  return m ? { chr: m[1], start: +m[2], end: +m[3] } : null
}

let slice = { key: null, chr: null, start: 0, end: 0 }

function ensureSlice(M, msg, view) {
  const sliceBams = msg.samples.map((s) => '/slice/' + s.name + '.bam')
  const need = slice.key !== remoteKey || view.chr !== slice.chr ||
    view.start < slice.start || view.end > slice.end
  if (need) {
    const span = Math.max(1, view.end - view.start)
    const margin = Math.min(Math.max(span * 3, 20000), 500000)
    const wStart = Math.max(1, view.start - margin)
    const wEnd = view.end + margin
    freshDir(M, '/slice')
    const remoteBams = msg.samples.map((s) => '/remote/' + s.cramName).join('\n')
    postMessage({ type: 'log', line: `Prefetching ${view.chr}:${wStart}-${wEnd}…` })
    const rc = M.ccall('wally_slice', 'number',
      ['string', 'string', 'string', 'string'],
      [remoteBams, '/ref/' + msg.refName, `${view.chr}:${wStart}-${wEnd}`, sliceBams.join('\n')])
    if (rc !== 0) throw new Error(`prefetch (wally_slice) failed (${rc})`)
    slice = { key: remoteKey, chr: view.chr, start: wStart, end: wEnd }
  }
  return {
    bams: sliceBams.join('\n'),
    genome: '/ref/' + msg.refName,
    prefetched: need
  }
}

self.onmessage = async (ev) => {
  const msg = ev.data
  if (msg.type !== 'render') return
  try {
    const M = await init()

    let bams, genome, bed = '', prefetched = false
    if (msg.mode === 'remote') {
      mountRemote(M, msg)
      const view = parseRegion(msg.region)
      if (!view) throw new Error('Invalid region: ' + msg.region)
      const s = ensureSlice(M, msg, view)
      bams = s.bams; genome = s.genome; prefetched = s.prefetched
      if (msg.anno && msg.annoName) bed = '/ref/' + msg.annoName
    } else {
      const m = mountLocal(M, msg)
      bams = m.bams; genome = m.genome
      if (msg.anno && msg.bedName) bed = '/data/' + msg.bedName
    }

    try { M.FS.unlink('/wallyplot.png') } catch (e) { }

    const xhr0 = xhrCount
    const t0 = performance.now()
    const rc = M.ccall('wally_region', 'number',
      ['string', 'string', 'string', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'string'],
      [bams, genome, msg.region, msg.width, msg.height,
        msg.paired ? 1 : 0, msg.clip ? 1 : 0, msg.supplementary ? 1 : 0, msg.coverage ? 1 : 0,
        msg.delsize || 1000, msg.mod || 0, msg.tlheight || 14, msg.rdheight || 12,
        Number.isFinite(msg.mapq) ? msg.mapq : 1, bed])
    const elapsed = Math.round(performance.now() - t0)

    if (rc !== 0) {
      postMessage({ type: 'error', message: `wally returned ${rc} - check the region, sample and reference genome and that indexes are present.` })
      return
    }
    const png = M.FS.readFile('/wallyplot.png')
    postMessage({ type: 'result', png, elapsed, xhr: xhrCount - xhr0, prefetched }, [png.buffer])
  } catch (e) {
    postMessage({ type: 'error', message: String(e && e.message ? e.message : e) })
  }
}
