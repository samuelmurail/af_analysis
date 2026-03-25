import { api, setStatus } from './api.js';
import { state, renderTable, renderSelectedResidues } from './table.js';
import { renderPlot, renderPaePlot, renderLisPlot, highlightPlotResidues, clearPlotSelection, reapplyPaePlotOverlay, getPlotZoom } from './plot.js';
import { loadStructure, highlightResidues, hoverResidues, unhoverResidues, applyPaeColors, clearPaeColors, clearMolstarSelection, subscribeToMolstarHover } from './molstar.js';
import { initResizableLayout } from './resize.js';

// Track the last rendered plot type so we only restore zoom when staying on the same type.
let _lastPlotType = null;

function collapseCard(cardId) {
  const card = document.getElementById(cardId);
  if (!card) return;
  const titleEl = card.querySelector('.card-title');
  const body = card.querySelector('.card-body');
  if (!titleEl || !body) return;
  body.style.display = 'none';
  titleEl.textContent = '\u25BA ' + titleEl.textContent.replace(/^[\u25BC\u25BA] /, '');
}

function renderCurrentTable(columns, rows) {
  renderTable(columns, rows, async (selectedRow) => {
    state.selectedModel = Number(selectedRow);
    renderCurrentTable(columns, rows);
    await refreshModelPanels();
  });
}

async function refreshModelPanels() {
  const plotType = document.getElementById('plot-type')?.value || 'plddt';

  // Capture zoom before re-rendering, but only restore it when staying on the same plot type.
  const savedZoom = (plotType === _lastPlotType) ? getPlotZoom() : null;
  _lastPlotType = plotType;

  const plddtPayload = await api(`/api/plddt?index=${state.selectedModel}`);
  state.chainIds = plddtPayload.chain_ids || [];
  state.chainLengths = plddtPayload.chain_lengths || [];
  state.chainTypes = plddtPayload.chain_types || [];

  if (plotType === 'pae') {
    try {
      const paePayload = await api(`/api/pae?index=${state.selectedModel}`);
      renderPaePlot(paePayload, makePaeHandlers(), savedZoom);
      if (state.paeSelection) {
        const { xResidues, yResidues } = state.paeSelection;
        reapplyPaePlotOverlay(xResidues, yResidues);
        const el = document.getElementById("selected-residues");
        if (el) el.textContent = `Scored: [${xResidues.join(", ")}] | Aligned: [${yResidues.join(", ")}]`;
      }
    } catch (e) {
      renderPlot(plddtPayload, makePlddtHandlers(), savedZoom);
      document.getElementById('plot-type').value = 'plddt';
    }
  } else if (['lis', 'lia', 'iptm_d0_matrix', 'ipsae_matrix'].includes(plotType)) {
    const LIS_API = { lis: 'lis', lia: 'lia', iptm_d0_matrix: 'iptm_d0', ipsae_matrix: 'ipsae' };
    try {
      const payload = await api(`/api/${LIS_API[plotType]}?index=${state.selectedModel}`);
      renderLisPlot(payload, makeLisHandlers());
    } catch (e) {
      renderPlot(plddtPayload, makePlddtHandlers());
      document.getElementById('plot-type').value = 'plddt';
    }
  } else {
    renderPlot(plddtPayload, makePlddtHandlers(), savedZoom);
  }

  renderSelectedResidues();
  await loadStructure(state.selectedModel);
  if (plotType === 'pae' && state.paeSelection) {
    await applyPaeColors(state.paeSelection.xResidues, state.paeSelection.yResidues);
  }
  subscribeToMolstarHover((globalResidues) => {
    highlightPlotResidues(globalResidues);
  });
}

function makePaeHandlers() {
  return {
    onClick: (residue) => {
      state.paeSelection = { xResidues: [residue], yResidues: [] };
      const el = document.getElementById("selected-residues");
      if (el) el.textContent = `Scored: [${residue}] | Aligned: []`;
      applyPaeColors([residue], []);
    },
    onPaeSelect: ({ xResidues, yResidues }) => {
      state.paeSelection = (xResidues.length > 0 || yResidues.length > 0)
        ? { xResidues, yResidues } : null;
      const el = document.getElementById("selected-residues");
      if (el) el.textContent = `Scored: [${xResidues.join(", ")}] | Aligned: [${yResidues.join(", ")}]`;
      if (xResidues.length > 0 || yResidues.length > 0) {
        applyPaeColors(xResidues, yResidues);
      }
    },
    onDeselect: () => {
      state.paeSelection = null;
      clearPaeColors();
    },
    onHover: ({ xResidue, yResidue }) => {
      const residues = [...new Set([xResidue, yResidue].filter(Number.isInteger))];
      hoverResidues(residues);
    },
    onUnhover: () => {
      unhoverResidues();
    },
  };
}

function _chainIdsToResidues(chainIds) {
  const ids = state.chainIds || [];
  const lengths = state.chainLengths || [];
  const starts = [];
  let offset = 1;
  for (const len of lengths) {
    starts.push(offset);
    offset += len;
  }
  const residues = [];
  for (const chainId of chainIds) {
    const idx = ids.indexOf(chainId);
    if (idx < 0) continue;
    const start = starts[idx];
    const len = lengths[idx];
    for (let r = start; r < start + len; r++) residues.push(r);
  }
  return residues;
}

function makeLisHandlers() {
  return {
    onClick: ({ xChainId, yChainId }) => {
      const residues = _chainIdsToResidues([...new Set([xChainId, yChainId])]);
      state.selectedResidues = residues;
      renderSelectedResidues();
      if (residues.length) highlightResidues(residues);
    },
    onHover: ({ xChainId, yChainId }) => {
      const residues = _chainIdsToResidues([...new Set([xChainId, yChainId])]);
      if (residues.length) hoverResidues(residues);
    },
    onUnhover: () => { unhoverResidues(); },
  };
}

function makePlddtHandlers() {
  return {
    onClick: (residue) => {
      state.selectedResidues = [residue];
      renderSelectedResidues();
      applyPaeColors([residue], []);
    },
    onSelect: (residues) => {
      state.selectedResidues = residues;
      renderSelectedResidues();
      if (residues.length > 0) applyPaeColors(residues, []);
      else clearPaeColors();
    },
    onHover: (residues) => {
      hoverResidues(residues);
    },
    onUnhover: () => {
      unhoverResidues();
    },
  };
}

function _updateLisOption(hasLis, hasLia, hasIptmD0, hasIpsae) {
  const optLis = document.getElementById('opt-lis');
  const optLia = document.getElementById('opt-lia');
  const optIptm = document.getElementById('opt-iptm-d0');
  const optIpsae = document.getElementById('opt-ipsae');
  const plotType = document.getElementById('plot-type');
  if (optLis) optLis.style.display = hasLis ? '' : 'none';
  if (optLia) optLia.style.display = hasLia ? '' : 'none';
  if (optIptm) optIptm.style.display = hasIptmD0 ? '' : 'none';
  if (optIpsae) optIpsae.style.display = hasIpsae ? '' : 'none';
  if (plotType) {
    if (!hasLis && plotType.value === 'lis') plotType.value = 'plddt';
    if (!hasLia && plotType.value === 'lia') plotType.value = 'plddt';
    if (!hasIptmD0 && plotType.value === 'iptm_d0_matrix') plotType.value = 'plddt';
    if (!hasIpsae && plotType.value === 'ipsae_matrix') plotType.value = 'plddt';
  }
}

async function refreshTableAndPanels() {
  const table = await api('/api/table');
  state.tableRows = table.rows;
  state.selectedModel = 0;
  _updateLisOption(!!table.has_lis, !!table.has_lia, !!table.has_iptm_d0_matrix, !!table.has_ipsae_matrix);
  renderCurrentTable(table.columns, table.rows);
  await refreshModelPanels();
}

function initEvents() {
  const loadBtn = document.getElementById('load-btn');
  if (!loadBtn) return;

  loadBtn.addEventListener('click', async () => {
    const directory = document.getElementById('directory').value.trim();
    const format = document.getElementById('format').value;
    if (!directory) {
      setStatus('Please provide a directory path', true);
      return;
    }

    const progressEl = document.getElementById('load-progress');
    const barEl      = document.getElementById('load-progress-bar');
    const labelEl    = document.getElementById('load-progress-label');
    const showProgress = (desc, n, total) => {
      if (!progressEl) return;
      progressEl.style.display = '';
      const pct = total > 0 ? Math.round(100 * n / total) : 0;
      if (barEl)   barEl.style.width = `${pct}%`;
      if (labelEl) labelEl.textContent = total > 0
        ? `${desc}: ${n} / ${total} (${pct}%)`
        : `${desc}…`;
    };
    const hideProgress = () => { if (progressEl) progressEl.style.display = 'none'; };

    setStatus('Loading dataset...');
    hideProgress();

    // Open SSE stream before the POST so no events are missed.
    let sseResolve;
    const sseDone = new Promise(r => { sseResolve = r; });
    const es = new EventSource('/api/progress/stream');
    es.onmessage = (ev) => {
      try {
        const d = JSON.parse(ev.data);
        if (d.desc === '__end__') { es.close(); sseResolve(); }
        else { showProgress(d.desc, d.n, d.total); }
      } catch {}
    };
    es.onerror = () => { es.close(); sseResolve(); };

    try {
      const payload = await api('/api/load', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ directory, format })
      });
      await sseDone;
      hideProgress();
      setStatus(`Dataset loaded (${payload.rows} rows)`);
      collapseCard('load-card');
      await refreshTableAndPanels();
    } catch (err) {
      await sseDone;
      hideProgress();
      setStatus(String(err), true);
    }
  });

  document.getElementById('plot-type')?.addEventListener('change', async () => {
    if (state.selectedModel != null) {
      await refreshModelPanels();
    }
  });

  document.getElementById('color-scheme')?.addEventListener('change', async () => {
    if (state.selectedModel != null) {
      await loadStructure(state.selectedModel);
    }
  });

  document.getElementById('clear-selection-btn')?.addEventListener('click', () => {
    state.selectedResidues = [];
    state.paeSelection = null;
    renderSelectedResidues();
    clearPlotSelection();
    clearPaeColors();
    clearMolstarSelection();
  });

  document.getElementById('compute-btn')?.addEventListener('click', async () => {
    const scores = ['pdockq2', 'LIS', 'LIA', 'iptm_d0'].filter(
      id => document.getElementById(`compute-${id}`)?.checked
    );
    // iptm_d0 implies ipSAE — handled server-side, no separate checkbox needed.
    if (!scores.length) return;
    const statusEl  = document.getElementById('compute-status');
    const progressEl = document.getElementById('compute-progress');
    const barEl     = document.getElementById('progress-bar');
    const labelEl   = document.getElementById('progress-label');

    const showProgress = (desc, n, total) => {
      if (!progressEl) return;
      progressEl.style.display = '';
      const pct = total > 0 ? Math.round(100 * n / total) : 0;
      if (barEl)  barEl.style.width = `${pct}%`;
      if (labelEl) labelEl.textContent = total > 0
        ? `${desc}: ${n} / ${total} (${pct}%)`
        : `${desc}…`;
    };
    const hideProgress = () => { if (progressEl) progressEl.style.display = 'none'; };

    if (statusEl) { statusEl.textContent = 'Computing…'; statusEl.style.color = '#2d3a57'; }
    hideProgress();

    try {
      const added = [];
      for (const score of scores) {
        // Open SSE stream first, then fire the compute request.
        let sseResolve;
        const sseDone = new Promise(r => { sseResolve = r; });
        const es = new EventSource('/api/progress/stream');
        es.onmessage = (ev) => {
          try {
            const d = JSON.parse(ev.data);
            if (d.desc === '__end__') {
              es.close();
              sseResolve();
            } else {
              showProgress(d.desc, d.n, d.total);
            }
          } catch {}
        };
        es.onerror = () => { es.close(); sseResolve(); };

        const lisPaeCutoff   = parseFloat(document.getElementById('lis-pae-cutoff')?.value   ?? 12);
        const liaPaeCutoff   = parseFloat(document.getElementById('lia-pae-cutoff')?.value   ?? 12);
        const liaDistCutoff  = parseFloat(document.getElementById('lia-dist-cutoff')?.value  ?? 8);
        const ipsaePaeCutoff = parseFloat(document.getElementById('ipsae-pae-cutoff')?.value ?? 10);
        const extraParams = score === 'LIS'
          ? { pae_cutoff: lisPaeCutoff }
          : score === 'LIA'
            ? { pae_cutoff: liaPaeCutoff, dist_cutoff: liaDistCutoff }
            : score === 'iptm_d0'
              ? { pae_cutoff: ipsaePaeCutoff }
              : {};
        const res = await api('/api/compute', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ score, ...extraParams }),
        });
        await sseDone;  // wait for last SSE event before moving on
        added.push(...(res.columns_added || []));
      }
      hideProgress();
      if (statusEl) statusEl.textContent = added.length
        ? `Done. New columns: ${added.join(', ')}`
        : 'Done (no new columns added).';
      collapseCard('compute-card');
      await refreshTableAndPanels();
    } catch (err) {
      hideProgress();
      if (statusEl) { statusEl.textContent = String(err); statusEl.style.color = '#a11927'; }
    }
  });
}

// ── Directory browser ────────────────────────────────────────────────────────

let _browseCurrent = null;

async function _browseLoad(path) {
  const data = await api(`/api/browse?path=${encodeURIComponent(path)}`);
  _browseCurrent = data.path;
  document.getElementById('browse-path').textContent = data.path;
  document.getElementById('browse-up').disabled = !data.parent;
  const list = document.getElementById('browse-list');
  list.innerHTML = '';
  for (const dir of data.dirs) {
    const li = document.createElement('li');
    li.textContent = '📁 ' + dir;
    li.style.cssText = 'padding:7px 12px; cursor:pointer; font-size:13px; border-bottom:1px solid #f0f2f7;';
    li.addEventListener('mouseover', () => li.style.background = '#f0f4ff');
    li.addEventListener('mouseout',  () => li.style.background = '');
    li.addEventListener('click', () => _browseLoad(data.path + '/' + dir));
    list.appendChild(li);
  }
  if (!data.dirs.length) {
    const li = document.createElement('li');
    li.textContent = '(no subdirectories)';
    li.style.cssText = 'padding:10px 12px; color:#888; font-size:13px;';
    list.appendChild(li);
  }
}

function initBrowser() {
  const modal = document.getElementById('browse-modal');
  const close = () => { modal.style.display = 'none'; };

  document.getElementById('browse-btn').addEventListener('click', async () => {
    modal.style.display = 'flex';
    const current = document.getElementById('directory').value.trim() || '~';
    await _browseLoad(current);
  });

  document.getElementById('browse-up').addEventListener('click', async () => {
    const pathEl = document.getElementById('browse-path').textContent;
    if (!pathEl) return;
    const parent = pathEl.split('/').slice(0, -1).join('/') || '/';
    await _browseLoad(parent);
  });

  document.getElementById('browse-select').addEventListener('click', () => {
    if (_browseCurrent) document.getElementById('directory').value = _browseCurrent;
    close();
  });

  document.getElementById('browse-cancel').addEventListener('click', close);
  document.getElementById('browse-cancel2').addEventListener('click', close);
  modal.addEventListener('click', (e) => { if (e.target === modal) close(); });
}

// ── Init ─────────────────────────────────────────────────────────────────────

async function main() {
  initResizableLayout();
  initEvents();
  initBrowser();
  renderSelectedResidues();

  // If the server already has data loaded (e.g. --directory CLI arg), populate
  // the directory input and refresh all panels without requiring a manual load.
  try {
    const health = await api('/api/health');
    if (health.loaded) {
      if (health.directory) {
        const dirInput = document.getElementById('directory');
        if (dirInput) dirInput.value = health.directory;
      }
      setStatus(`Dataset loaded — ${health.rows} models.`);
      collapseCard('load-card');
      await refreshTableAndPanels();
    }
  } catch (_) { /* server not ready yet, ignore */ }
}

main().catch((err) => {
  console.error(err);
  setStatus(String(err), true);
});
