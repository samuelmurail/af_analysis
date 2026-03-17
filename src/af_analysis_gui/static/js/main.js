import { api, setStatus } from './api.js';
import { state, renderTable, renderSelectedResidues } from './table.js';
import { renderPlot, renderPaePlot, renderLisPlot, highlightPlotResidues, clearPlotSelection, reapplyPaePlotOverlay } from './plot.js';
import { loadStructure, highlightResidue, highlightResidues, hoverResidues, unhoverResidues, applyPaeColors, clearPaeColors, clearMolstarSelection, subscribeToMolstarHover } from './molstar.js';
import { initResizableLayout } from './resize.js';

function renderCurrentTable(columns, rows) {
  renderTable(columns, rows, async (selectedRow) => {
    state.selectedModel = Number(selectedRow);
    renderCurrentTable(columns, rows);
    await refreshModelPanels();
  });
}

async function refreshModelPanels() {
  const plotType = document.getElementById('plot-type')?.value || 'plddt';

  const plddtPayload = await api(`/api/plddt?index=${state.selectedModel}`);
  state.chainIds = plddtPayload.chain_ids || [];
  state.chainLengths = plddtPayload.chain_lengths || [];

  if (plotType === 'pae') {
    try {
      const paePayload = await api(`/api/pae?index=${state.selectedModel}`);
      renderPaePlot(paePayload, makePaeHandlers());
      if (state.paeSelection) {
        const { xResidues, yResidues } = state.paeSelection;
        reapplyPaePlotOverlay(xResidues, yResidues);
        const el = document.getElementById("selected-residues");
        if (el) el.textContent = `Scored: [${xResidues.join(", ")}] | Aligned: [${yResidues.join(", ")}]`;
      }
    } catch (e) {
      renderPlot(plddtPayload, makePlddtHandlers());
      document.getElementById('plot-type').value = 'plddt';
    }
  } else if (plotType === 'lis') {
    try {
      const lisPayload = await api(`/api/lis?index=${state.selectedModel}`);
      renderLisPlot(lisPayload, makeLisHandlers());
    } catch (e) {
      renderPlot(plddtPayload, makePlddtHandlers());
      document.getElementById('plot-type').value = 'plddt';
    }
  } else if (plotType === 'lia') {
    try {
      const liaPayload = await api(`/api/lia?index=${state.selectedModel}`);
      renderLisPlot(liaPayload, makeLisHandlers());
    } catch (e) {
      renderPlot(plddtPayload, makePlddtHandlers());
      document.getElementById('plot-type').value = 'plddt';
    }
  } else if (plotType === 'iptm_d0_matrix') {
    try {
      const iptmPayload = await api(`/api/iptm_d0?index=${state.selectedModel}`);
      renderLisPlot(iptmPayload, makeLisHandlers());
    } catch (e) {
      renderPlot(plddtPayload, makePlddtHandlers());
      document.getElementById('plot-type').value = 'plddt';
    }
  } else {
    renderPlot(plddtPayload, makePlddtHandlers());
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
      highlightResidue(residue);
    },
    onSelect: (residues) => {
      state.selectedResidues = residues;
      renderSelectedResidues();
      if (residues.length > 0) highlightResidues(residues);
    },
    onHover: (residues) => {
      hoverResidues(residues);
    },
    onUnhover: () => {
      unhoverResidues();
    },
  };
}

function _updateLisOption(hasLis, hasLia, hasIptmD0) {
  const optLis = document.getElementById('opt-lis');
  const optLia = document.getElementById('opt-lia');
  const optIptm = document.getElementById('opt-iptm-d0');
  const plotType = document.getElementById('plot-type');
  if (optLis) optLis.style.display = hasLis ? '' : 'none';
  if (optLia) optLia.style.display = hasLia ? '' : 'none';
  if (optIptm) optIptm.style.display = hasIptmD0 ? '' : 'none';
  if (plotType) {
    if (!hasLis && plotType.value === 'lis') plotType.value = 'plddt';
    if (!hasLia && plotType.value === 'lia') plotType.value = 'plddt';
    if (!hasIptmD0 && plotType.value === 'iptm_d0_matrix') plotType.value = 'plddt';
  }
}

async function refreshTableAndPanels() {
  const table = await api('/api/table');
  state.tableRows = table.rows;
  state.selectedModel = 0;
  _updateLisOption(!!table.has_lis, !!table.has_lia, !!table.has_iptm_d0_matrix);
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

    setStatus('Loading dataset...');
    try {
      const payload = await api('/api/load', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ directory, format })
      });
      setStatus(`Dataset loaded (${payload.rows} rows)`);
      await refreshTableAndPanels();
    } catch (err) {
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
    if (!scores.length) return;
    const statusEl = document.getElementById('compute-status');
    if (statusEl) { statusEl.textContent = 'Computing…'; statusEl.style.color = '#2d3a57'; }
    try {
      const added = [];
      for (const score of scores) {
        const res = await api('/api/compute', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ score }),
        });
        added.push(...(res.columns_added || []));
      }
      if (statusEl) statusEl.textContent = added.length
        ? `Done. New columns: ${added.join(', ')}`
        : 'Done (no new columns added).';
      // Refresh the table so new numeric columns appear.
      await refreshTableAndPanels();
    } catch (err) {
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
}

main().catch((err) => {
  console.error(err);
  setStatus(String(err), true);
});
