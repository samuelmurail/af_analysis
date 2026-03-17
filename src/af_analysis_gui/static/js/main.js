import { api, setStatus } from './api.js';
import { state, renderTable, renderSelectedResidues } from './table.js';
import { renderPlot, renderPaePlot, highlightPlotResidues, clearPlotSelection, reapplyPaePlotOverlay } from './plot.js';
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

async function refreshTableAndPanels() {
  const table = await api('/api/table');
  state.tableRows = table.rows;
  state.selectedModel = 0;
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
