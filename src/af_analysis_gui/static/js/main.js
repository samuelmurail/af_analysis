import { api, setStatus } from './api.js';
import { state, renderTable, renderSelectedResidues } from './table.js';
import { renderPlot, renderPaePlot, renderLisPlot, highlightPlotResidues, clearPlotSelection, reapplyPaePlotOverlay, getPlotZoom } from './plot.js';
import { loadStructure, highlightResidues, hoverResidues, unhoverResidues, applyPaeColors, clearPaeColors, clearMolstarSelection, subscribeToMolstarHover, subscribeToSuperposeHover, loadSuperpose, setClusterColor, getSuperposeState } from './molstar.js';
import { initResizableLayout } from './resize.js';

// Track the last rendered plot type so we only restore zoom when staying on the same type.
let _lastPlotType = null;

// Look up the CLUSTER_PALETTE hex integer for a given row index.
// Returns null if clustering data is not available or the row has no cluster.
function _clusterColorForRow(rowIdx) {
  if (!_clusterData) return null;
  for (const q of _clusterData) {
    const pt = q.points.find(p => p.row === rowIdx);
    if (pt && pt.cluster != null) {
      const clusterNums = [...new Set(q.points.map(p => p.cluster))].sort((a, b) => a - b);
      const ci = clusterNums.indexOf(pt.cluster);
      const hex = CLUSTER_PALETTE[ci % CLUSTER_PALETTE.length];
      return parseInt(hex.replace('#', ''), 16);
    }
  }
  return null;
}

// Build a per-frame colour array (hex ints) for a set of rows based on cluster membership.
function _frameColorsForRows(rows) {
  return rows.map(r => _clusterColorForRow(r) ?? _SUPERPOSE_PALETTE_INT[r % _SUPERPOSE_PALETTE_INT.length]);
}

// Integer version of CLUSTER_PALETTE for fallback use in _frameColorsForRows.
const _SUPERPOSE_PALETTE_INT = [
  0x636efa, 0xef553b, 0x00cc96, 0xab63fa, 0xffa15a,
  0x19d3f3, 0xff6692, 0xb6e880, 0xff97ff, 0xfecb52,
];

// Prepare molstar cluster coloring for the given row before loadStructure.
function _applyClusterColorForRow(rowIdx) {
  if (document.getElementById('color-scheme')?.value !== 'cluster') return;
  const col = _clusterColorForRow(rowIdx);
  setClusterColor(col ?? 0x636efa);
}
// Stored so cluster-plot clicks can re-render the table without refetching.
let _tableColumns = [];
let _tableRows = [];

// Return 'mds' | 'dendrogram' | null based on plot-type dropdown.
function getClusterView() {
  const v = document.getElementById('plot-type')?.value;
  return (v === 'mds' || v === 'dendrogram') ? v : null;
}

function enableToolbarButtons() {
  document.getElementById('toolbar-compute-btn')?.removeAttribute('disabled');
  document.getElementById('toolbar-cluster-btn')?.removeAttribute('disabled');
}

function _hideDialog(id) {
  const el = document.getElementById(id);
  if (el && typeof bootstrap !== 'undefined') bootstrap.Modal.getInstance(el)?.hide();
}

function _resetForNewDataset() {
  _clusterData = null;
  _lastPlotType = null;
  // Hide cluster plot types
  const optMds = document.getElementById('opt-mds');
  const optDendro = document.getElementById('opt-dendro');
  if (optMds) optMds.style.display = 'none';
  if (optDendro) optDendro.style.display = 'none';
  // Reset to pLDDT view
  const plotType = document.getElementById('plot-type');
  if (plotType) plotType.value = 'plddt';
  document.getElementById('plddt-plot').style.display = '';
  document.getElementById('cluster-main-plot').style.display = 'none';
  document.getElementById('cluster-controls').style.display = 'none';
}

// Wire up superpose hover → highlight matching row in table and MDS plot.
function _attachSuperposeHover() {
  subscribeToSuperposeHover((rowIdx) => {
    state.hoveredRow = rowIdx;
    renderCurrentTable(_tableColumns, _tableRows);
    const plotDiv = document.getElementById('cluster-main-plot');
    if (!plotDiv?._fullLayout) return;
    const view = getClusterView() ?? 'mds';
    if (view === 'mds') {
      // Highlight hovered point in MDS scatter plot.
      if (rowIdx != null) {
        const cp = _rowToCurveAndPoint(rowIdx);
        if (cp) Plotly.Fx.hover(plotDiv, [{ curveNumber: cp.curve, pointNumber: cp.point }]);
      } else {
        Plotly.Fx.unhover(plotDiv);
      }
    } else {
      // Highlight hovered leaf in dendrogram with a vertical line shape.
      const leaves   = plotDiv._afDendroLeaves;
      const tickvals = plotDiv._afDendroTickvals;
      if (!leaves || !tickvals) return;
      if (rowIdx != null) {
        const leafIdx = leaves.indexOf(rowIdx);
        if (leafIdx !== -1) {
          const xVal = tickvals[leafIdx];
          const yMax = plotDiv._fullLayout?.yaxis?.range?.[1] ?? 1;
          Plotly.relayout(plotDiv, {
            shapes: [{
              type: 'line', xref: 'x', yref: 'y',
              x0: xVal, x1: xVal, y0: 0, y1: yMax,
              line: { color: '#636efa', width: 2.5, dash: 'dot' },
            }],
          });
        }
      } else {
        Plotly.relayout(plotDiv, { shapes: [] });
      }
    }
  });
}

// Map a dataset row index to { curve, point } indices in the current MDS Plotly chart.
function _rowToCurveAndPoint(rowIdx) {
  if (!_clusterData) return null;
  const queryVal = document.getElementById('cluster-query-select')?.value;
  const q = _clusterData.find(x => x.query === queryVal);
  if (!q) return null;
  const clusterNums = [...new Set(q.points.map(p => p.cluster))].sort((a, b) => a - b);
  for (let ci = 0; ci < clusterNums.length; ci++) {
    const pts = q.points.filter(p => p.cluster === clusterNums[ci]);
    const pi = pts.findIndex(p => p.row === rowIdx);
    if (pi !== -1) return { curve: ci, point: pi };
  }
  return null;
}

// Returns a luminance-contrasted text color ('#000' or '#fff') for a hex bg color.
function _contrastColor(hex) {
  const r = parseInt(hex.slice(1, 3), 16);
  const g = parseInt(hex.slice(3, 5), 16);
  const b = parseInt(hex.slice(5, 7), 16);
  return (r * 299 + g * 587 + b * 114) / 1000 > 128 ? '#000' : '#fff';
}

// Returns a cellStyleFn for renderTable that colors the 'cluster' column with CLUSTER_PALETTE.
function _makeClusterCellStyleFn() {
  if (!_clusterData) return null;
  return (colName, row) => {
    if (colName !== 'cluster') return null;
    const clust = row['cluster'];
    if (clust == null || (typeof clust === 'number' && isNaN(clust))) return null;
    const qData = _clusterData.find(q => q.query === row['query']);
    if (!qData) return null;
    const clusterNums = [...new Set(qData.points.map(p => p.cluster))].sort((a, b) => a - b);
    const ci = clusterNums.indexOf(Number(clust));
    if (ci === -1) return null;
    const color = CLUSTER_PALETTE[ci % CLUSTER_PALETTE.length];
    return `background-color: ${color}; color: ${_contrastColor(color)}; font-weight: bold;`;
  };
}

function renderCurrentTable(columns, rows) {
  _tableColumns = columns;
  _tableRows = rows;
  const cellStyleFn = _makeClusterCellStyleFn();
  renderTable(columns, rows, async (selectedRow) => {
    state.selectedModel = Number(selectedRow);
    state.selectedRows = new Set();  // clear superpose highlight on single-model click
    state.hoveredRow = null;
    // Switch away from cluster coloring when viewing a single model.
    const schemeEl = document.getElementById('color-scheme');
    if (schemeEl?.value === 'cluster') schemeEl.value = 'chain-id';
    renderCurrentTable(columns, rows);
    await refreshModelPanels();
  }, { cellStyleFn });
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
  } else if (plotType !== 'mds' && plotType !== 'dendrogram') {
    renderPlot(plddtPayload, makePlddtHandlers(), savedZoom);
  }

  if (plotType !== 'mds' && plotType !== 'dendrogram') renderSelectedResidues();
  _applyClusterColorForRow(state.selectedModel);
  await loadStructure(state.selectedModel);
  if (plotType === 'pae' && state.paeSelection) {
    await applyPaeColors(state.paeSelection.xResidues, state.paeSelection.yResidues);
  }
  _renderClusterPlot();
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

  // Resize Plotly whenever the plot panel changes size.
  const plotPanelEl = document.getElementById('plot-panel');
  if (plotPanelEl && typeof ResizeObserver !== 'undefined') {
    new ResizeObserver(() => {
      const pEl = document.getElementById('plddt-plot');
      if (pEl?._fullLayout) Plotly.Plots.resize(pEl);
      const cEl = document.getElementById('cluster-main-plot');
      if (cEl?._fullLayout) Plotly.Plots.resize(cEl);
    }).observe(plotPanelEl);
  }


  loadBtn.addEventListener('click', async () => {
    const path = document.getElementById('directory').value.trim();
    const format = document.getElementById('format').value;
    const csvMode = path.toLowerCase().endsWith('.csv');
    const directory = csvMode ? '' : path;
    if (!path) {
      setStatus('Please provide a directory or CSV file path', true);
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
        body: JSON.stringify(csvMode ? { csv: path } : { directory, format }),
      });
      await sseDone;
      hideProgress();
      setStatus(`Dataset loaded (${payload.rows} rows)`);
      _hideDialog('load-dialog');
      enableToolbarButtons();
      _resetForNewDataset();
      await refreshTableAndPanels();
    } catch (err) {
      await sseDone;
      hideProgress();
      setStatus(String(err), true);
    }
  });

  document.getElementById('plot-type')?.addEventListener('change', async () => {
    const pt = document.getElementById('plot-type')?.value;
    const isCluster = pt === 'mds' || pt === 'dendrogram';
    document.getElementById('plddt-plot').style.display = isCluster ? 'none' : '';
    document.getElementById('cluster-main-plot').style.display = isCluster ? '' : 'none';
    document.getElementById('cluster-controls').style.display = isCluster ? '' : 'none';
    if (isCluster) {
      _renderClusterPlot();
    } else if (state.selectedModel != null) {
      await refreshModelPanels();
    }
  });

  // Toggle plot panel visibility.
  document.getElementById('toggle-plot-btn')?.addEventListener('click', () => {
    const plotPanel = document.getElementById('plot-panel');
    const splitterV = document.getElementById('splitter-v');
    const molstarPanel = document.getElementById('molstar-panel');
    const btn = document.getElementById('toggle-plot-btn');
    if (plotPanel.classList.contains('plot-hidden')) {
      plotPanel.classList.remove('plot-hidden');
      if (splitterV) splitterV.style.display = '';
      molstarPanel.style.flex = '';
      if (btn) btn.innerHTML = 'Plot &#9654;';
      // Resize plots after showing
      setTimeout(() => {
        const p = document.getElementById('plddt-plot');
        if (p?._fullLayout) Plotly.Plots.resize(p);
        const c = document.getElementById('cluster-main-plot');
        if (c?._fullLayout) Plotly.Plots.resize(c);
      }, 50);
    } else {
      plotPanel.classList.add('plot-hidden');
      if (splitterV) splitterV.style.display = 'none';
      molstarPanel.style.flex = '1';
      if (btn) btn.innerHTML = '&#9664; Plot';
    }
  });

  document.getElementById('color-scheme')?.addEventListener('change', async () => {
    const sp = getSuperposeState();
    if (sp) {
      // Re-render the active superpose with the new color scheme.
      await loadSuperpose(sp.rows, sp.query, sp.frameColors);
      _attachSuperposeHover();
    } else if (state.selectedModel != null) {
      _applyClusterColorForRow(state.selectedModel);
      await loadStructure(state.selectedModel);
    }
  });

  document.getElementById('export-csv-btn')?.addEventListener('click', () => {
    const a = document.createElement('a');
    a.href = '/api/export_csv';
    a.download = 'af_analysis.csv';
    a.click();
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
    const scores = ['pdockq2', 'LIS_LIA', 'iptm_d0'].filter(
      id => document.getElementById(`compute-${id}`)?.checked
    );
    // LIS_LIA computes both LIS and LIA in one call; iptm_d0 also includes ipSAE server-side.
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
        const liaDistCutoff  = parseFloat(document.getElementById('lia-dist-cutoff')?.value  ?? 8);
        const ipsaePaeCutoff = parseFloat(document.getElementById('ipsae-pae-cutoff')?.value ?? 10);
        const extraParams = score === 'LIS_LIA'
          ? { pae_cutoff: lisPaeCutoff, dist_cutoff: liaDistCutoff }
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
      _hideDialog('compute-dialog');
      await refreshTableAndPanels();
    } catch (err) {
      hideProgress();
      if (statusEl) { statusEl.textContent = String(err); statusEl.style.color = '#a11927'; }
    }
  });

  document.getElementById('cluster-btn')?.addEventListener('click', async () => {
    const threshold   = parseFloat(document.getElementById('cluster-threshold')?.value ?? 2.0);
    const alignSel    = document.getElementById('cluster-align-sel')?.value.trim() || 'backbone';
    const distSel     = document.getElementById('cluster-dist-sel')?.value.trim()  || 'backbone';
    const rmsdScale = !!document.getElementById('cluster-rmsd-scale')?.checked;
    const statusEl = document.getElementById('cluster-status');
    if (statusEl) { statusEl.textContent = 'Clustering…'; statusEl.style.color = '#2d3a57'; }
    _hideDialog('cluster-dialog');
    try {
      const data = await api('/api/cluster', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ threshold, align_selection: alignSel, distance_selection: distSel, rmsd_scale: rmsdScale }),
      });
      if (statusEl) statusEl.textContent = '';
      renderClusterPanel(data.queries);
      await refreshTableAndPanels();
    } catch (err) {
      if (statusEl) {
        statusEl.style.color = '#a11927';
        statusEl.title = err.details || '';
        statusEl.textContent = String(err);
      }
    }
  });
}

// ── Cluster panel ────────────────────────────────────────────────────────────

const CLUSTER_PALETTE = [
  '#636EFA','#EF553B','#00CC96','#AB63FA','#FFA15A',
  '#19D3F3','#FF6692','#B6E880','#FF97FF','#FECB52',
];

let _clusterData = null;
let _clusterListenersAttached = false;

function renderClusterPanel(queries) {
  _clusterData = queries;

  // Show cluster plot type options.
  const optMds = document.getElementById('opt-mds');
  const optDendro = document.getElementById('opt-dendro');
  if (optMds) optMds.style.display = '';
  if (optDendro) optDendro.style.display = '';

  // Populate query dropdown.
  const querySelect = document.getElementById('cluster-query-select');
  if (querySelect) {
    querySelect.innerHTML = '';
    for (const q of queries) {
      const opt = document.createElement('option');
      opt.value = q.query;
      opt.textContent = q.query;
      querySelect.appendChild(opt);
    }
  }

  // Attach query dropdown change listener once.
  if (!_clusterListenersAttached) {
    _clusterListenersAttached = true;
    querySelect?.addEventListener('change', _renderClusterPlot);
  }

  // Switch to MDS view.
  const plotType = document.getElementById('plot-type');
  if (plotType) plotType.value = 'mds';
  document.getElementById('plddt-plot').style.display = 'none';
  document.getElementById('cluster-main-plot').style.display = '';
  document.getElementById('cluster-controls').style.display = '';

  _renderClusterPlot();
}

function _renderClusterPlot() {
  if (!_clusterData) return;
  const view = getClusterView();
  if (!view) return;
  const queryVal = document.getElementById('cluster-query-select')?.value;
  const q = _clusterData.find(x => x.query === queryVal);
  if (!q) return;
  const plotDiv = document.getElementById('cluster-main-plot');
  if (!plotDiv) return;
  if (view === 'mds') {
    _renderMdsPlot(plotDiv, q);
  } else {
    _renderDendrogram(plotDiv, q);
  }
}

function _renderMdsPlot(plotDiv, q) {
  const clusterNums = [...new Set(q.points.map(p => p.cluster))].sort((a, b) => a - b);
  const traces = clusterNums.map((c, ci) => {
    const pts = q.points.filter(p => p.cluster === c);
    return {
      x: pts.map(p => p.x),
      y: pts.map(p => p.y),
      customdata: pts.map(p => p.row),
      text: pts.map(p => `Model ${p.row}<br>Cluster ${p.cluster}`),
      mode: 'markers',
      type: 'scatter',
      name: `Cluster ${c}`,
      marker: { size: 8, color: CLUSTER_PALETTE[ci % CLUSTER_PALETTE.length] },
      hovertemplate: '%{text}<extra></extra>',
    };
  });

  // Add a star marker for the currently selected model (single-view mode).
  if (state.selectedRows.size === 0) {
    const selPt = q.points.find(p => p.row === state.selectedModel);
    if (selPt) {
      traces.push({
        x: [selPt.x], y: [selPt.y],
        mode: 'markers',
        type: 'scatter',
        name: 'Selected',
        showlegend: false,
        hoverinfo: 'skip',
        marker: { size: 14, symbol: 'star', color: '#222', line: { color: '#fff', width: 1 } },
      });
    }
  }

  Plotly.react(plotDiv, traces, {
    autosize: true,
    margin: { l: 50, r: 20, t: 20, b: 50 },
    xaxis: { title: 'MDS 1', zeroline: false },
    yaxis: { title: 'MDS 2', zeroline: false },
    legend: { font: { size: 11 } },
    showlegend: q.n_clusters > 1,
    selections: [],
  }, { responsive: true });

  if (!plotDiv._afClusterClickAttached) {
    plotDiv._afClusterClickAttached = true;
    plotDiv.on('plotly_hover', (ev) => {
      if (getClusterView() !== 'mds') return;
      const pt = ev?.points?.[0];
      if (pt?.customdata == null) return;
      state.hoveredRow = pt.customdata;
      renderCurrentTable(_tableColumns, _tableRows);
    });
    plotDiv.on('plotly_unhover', () => {
      if (getClusterView() !== 'mds') return;
      state.hoveredRow = null;
      renderCurrentTable(_tableColumns, _tableRows);
    });
    plotDiv.on('plotly_click', async (ev) => {
      if (getClusterView() !== 'mds') return;
      if (!ev?.points?.length) return;
      const rowIdx = ev.points[0].customdata;
      state.selectedModel = rowIdx;
      state.selectedRows = new Set();
      state.hoveredRow = null;
      // Switch away from cluster coloring when viewing a single model.
      const schemeEl = document.getElementById('color-scheme');
      if (schemeEl?.value === 'cluster') schemeEl.value = 'chain-id';
      renderCurrentTable(_tableColumns, _tableRows);
      await refreshModelPanels();
    });
    plotDiv.on('plotly_selected', async (ev) => {
      if (getClusterView() !== 'mds') return;
      if (!ev?.points?.length) return;
      const rows = [...new Set(ev.points.map(p => p.customdata))];
      const queryVal = document.getElementById('cluster-query-select')?.value ?? '';
      const schemeEl = document.getElementById('color-scheme');
      if (schemeEl) schemeEl.value = 'cluster';
      state.selectedRows = new Set(rows);
      renderCurrentTable(_tableColumns, _tableRows);
      await loadSuperpose(rows, queryVal, _frameColorsForRows(rows));
      _attachSuperposeHover();
    });
  }
}

// Returns an array (one entry per arm) of Sets of row indices under that arm.
function _computeArmLeaves(d, tickvals) {
  const n = d.icoord.length;
  // Map arm x-center (×100, rounded) → arm index for fast lookup.
  const byCenter = new Map();
  for (let i = 0; i < n; i++) {
    const key = Math.round((d.icoord[i][0] + d.icoord[i][3]) / 2 * 100);
    byCenter.set(key, i);
  }
  const cache = new Array(n).fill(null);
  function get(i) {
    if (cache[i]) return cache[i];
    const ic = d.icoord[i], dc = d.dcoord[i];
    const s = new Set();
    for (const [fx, fy] of [[ic[0], dc[0]], [ic[3], dc[3]]]) {
      if (fy === 0) {
        // Foot touches the baseline → this is a leaf.
        const li = tickvals.findIndex(tv => Math.abs(tv - fx) < 1e-6);
        if (li >= 0) s.add(d.leaves[li]);
      } else {
        // Foot connects to a child arm whose center is at fx.
        const ci = byCenter.get(Math.round(fx * 100));
        if (ci !== undefined) for (const r of get(ci)) s.add(r);
      }
    }
    return (cache[i] = s);
  }
  return Array.from({ length: n }, (_, i) => get(i));
}

function _renderDendrogram(plotDiv, q) {
  if (!q.dendrogram) {
    Plotly.react(plotDiv, [], {
      margin: { l: 50, r: 20, t: 30, b: 60 },
      title: { text: 'No dendrogram data available' },
    }, { responsive: true });
    return;
  }
  const d = q.dendrogram;

  // x-axis tick labels at leaf positions (5, 15, 25, …).
  const tickvals = d.leaves.map((_, i) => (2 * i + 1) * 5);
  const ticktext = d.leaves.map(r => `M${r}`);

  // Build row→cluster map (same sort order as MDS so palette colours match).
  const clusterNums = [...new Set(q.points.map(p => p.cluster))].sort((a, b) => a - b);
  const clusterColor = new Map(
    clusterNums.map((c, ci) => [c, CLUSTER_PALETTE[ci % CLUSTER_PALETTE.length]])
  );
  const rowToCluster = new Map(q.points.map(p => [p.row, p.cluster]));

  // Compute the set of row indices lying under each arm.
  const armLeaves = _computeArmLeaves(d, tickvals);

  // Colour each arm: single-cluster subtrees get the cluster colour, mixed → gray.
  const GRAY = '#adb5bd';
  const armColors = armLeaves.map(leaves => {
    const clusters = new Set([...leaves].map(r => rowToCluster.get(r)));
    clusters.delete(null); clusters.delete(undefined);
    if (clusters.size === 1) return clusterColor.get([...clusters][0]) ?? GRAY;
    return GRAY;
  });

  // One line trace per colour (null-separated polylines within each trace).
  const colorSegs = new Map();
  for (let i = 0; i < d.icoord.length; i++) {
    const c = armColors[i];
    if (!colorSegs.has(c)) colorSegs.set(c, { xs: [], ys: [] });
    const seg = colorSegs.get(c);
    seg.xs.push(...d.icoord[i], null);
    seg.ys.push(...d.dcoord[i], null);
  }
  const armTraces = [...colorSegs.entries()].map(([col, { xs, ys }]) => ({
    x: xs, y: ys,
    mode: 'lines', type: 'scatter',
    line: { color: col, width: 2 },
    hoverinfo: 'skip', showlegend: false,
  }));

  // Threshold line.
  const xMax = d.n * 10;
  const threshTrace = {
    x: [0, xMax], y: [d.threshold, d.threshold],
    mode: 'lines', type: 'scatter',
    line: { color: '#ef4444', width: 1.5, dash: 'dash' },
    name: `Threshold (${d.threshold.toFixed(2)})`,
    showlegend: true,
  };

  // Invisible click-target markers at the peak of every arm.
  // customdata = arm index; hovertemplate shows the leaf count.
  const ctX = [], ctY = [], ctCustom = [], ctText = [];
  for (let i = 0; i < d.icoord.length; i++) {
    ctX.push((d.icoord[i][0] + d.icoord[i][3]) / 2);
    ctY.push(d.dcoord[i][1]);   // top of the U-shape
    ctCustom.push(i);
    const nLeaves = armLeaves[i].size;
    ctText.push(`${nLeaves} model${nLeaves !== 1 ? 's' : ''} — click to superpose`);
  }
  const clickTrace = {
    x: ctX, y: ctY,
    customdata: ctCustom,
    text: ctText,
    mode: 'markers', type: 'scatter',
    marker: { size: 14, opacity: 0, color: '#000' },
    hovertemplate: '%{text}<extra></extra>',
    showlegend: false,
  };

  // Selectable leaf markers at (tickval, 0) — invisible but hittable by lasso/box selection.
  // customdata = row index so plotly_selected can collect the right rows.
  const leafMarkerTrace = {
    x: tickvals, y: tickvals.map(() => 0),
    customdata: d.leaves,
    mode: 'markers', type: 'scatter',
    marker: { size: 10, opacity: 0, color: '#000' },
    hoverinfo: 'skip', showlegend: false,
  };

  // Store leaf ordering so hover can look up x positions.
  plotDiv._afDendroLeaves = d.leaves;
  plotDiv._afDendroTickvals = tickvals;
  // Store leaf trace curve number so plotly_selected can filter out arm click targets.
  plotDiv._afLeafTraceIdx = armTraces.length + 2;

  Plotly.react(plotDiv, [...armTraces, threshTrace, clickTrace, leafMarkerTrace], {
    margin: { l: 50, r: 20, t: 20, b: 70 },
    xaxis: { tickvals, ticktext, tickangle: -45, showgrid: false },
    yaxis: { title: 'Distance', zeroline: false },
    legend: {
      font: { size: 13 },
      x: 1, xanchor: 'right',
      y: 1, yanchor: 'top',
      bgcolor: 'rgba(255,255,255,0.85)',
      bordercolor: '#ccc', borderwidth: 1,
    },
  }, { responsive: true });

  // Store arm leaves so the persistent click handler can look them up.
  plotDiv._afArmLeaves = armLeaves;

  if (!plotDiv._afDendroClickAttached) {
    plotDiv._afDendroClickAttached = true;
    plotDiv.on('plotly_click', async (ev) => {
      if (getClusterView() !== 'dendrogram') return;
      if (!ev?.points?.length) return;
      const pt = ev.points[0];
      if (typeof pt.customdata !== 'number') return;
      const leaves = [...(plotDiv._afArmLeaves?.[pt.customdata] ?? [])]
        .filter(r => typeof r === 'number' && r >= 0);
      if (!leaves.length) return;
      const queryVal = document.getElementById('cluster-query-select')?.value ?? '';
      const schemeEl = document.getElementById('color-scheme');
      if (schemeEl) schemeEl.value = 'cluster';
      state.selectedRows = new Set(leaves);
      renderCurrentTable(_tableColumns, _tableRows);
      await loadSuperpose(leaves, queryVal, _frameColorsForRows(leaves));
      _attachSuperposeHover();
    });
    plotDiv.on('plotly_selected', async (ev) => {
      if (getClusterView() !== 'dendrogram') return;
      if (!ev?.points?.length) return;
      // Collect row indices strictly from the leaf marker trace (not arm click targets).
      const leafIdx = plotDiv._afLeafTraceIdx;
      const rows = [...new Set(
        ev.points
          .filter(p => p.curveNumber === leafIdx && typeof p.customdata === 'number')
          .map(p => p.customdata)
      )];
      if (!rows.length) return;
      const queryVal = document.getElementById('cluster-query-select')?.value ?? '';
      const schemeEl = document.getElementById('color-scheme');
      if (schemeEl) schemeEl.value = 'cluster';
      state.selectedRows = new Set(rows);
      renderCurrentTable(_tableColumns, _tableRows);
      await loadSuperpose(rows, queryVal, _frameColorsForRows(rows));
      _attachSuperposeHover();
    });
  }
}

let _browseCurrent = null;

async function _browseLoad(path) {
  const data = await api(`/api/browse?path=${encodeURIComponent(path)}`);
  _browseCurrent = data.path;
  document.getElementById('browse-path').textContent = data.path;
  document.getElementById('browse-up').disabled = !data.parent;
  const list = document.getElementById('browse-list');
  list.innerHTML = '';

  const addItem = (icon, name, onClick, isSelected) => {
    const li = document.createElement('li');
    li.textContent = icon + ' ' + name;
    li.style.cssText = `padding:7px 12px; cursor:pointer; font-size:13px; border-bottom:1px solid #f0f2f7;${isSelected ? 'background:#e8eeff; font-weight:600;' : ''}`;
    li.addEventListener('mouseover', () => li.style.background = '#f0f4ff');
    li.addEventListener('mouseout',  () => li.style.background = isSelected ? '#e8eeff' : '');
    li.addEventListener('click', onClick);
    list.appendChild(li);
  };

  for (const dir of (data.dirs || [])) {
    addItem('📁', dir, () => _browseLoad(data.path + '/' + dir), false);
  }
  for (const csv of (data.csvs || [])) {
    const fullPath = data.path + '/' + csv;
    addItem('📄', csv, () => {
      _browseCurrent = fullPath;
      // Highlight the selected CSV and update the select button label.
      for (const el of list.querySelectorAll('li')) {
        el.style.background = el.textContent.trim() === '📄 ' + csv ? '#e8eeff' : '';
        el.style.fontWeight  = el.textContent.trim() === '📄 ' + csv ? '600' : '';
      }
    }, false);
  }

  if (!data.dirs?.length && !data.csvs?.length) {
    const li = document.createElement('li');
    li.textContent = '(empty)';
    li.style.cssText = 'padding:10px 12px; color:#888; font-size:13px;';
    list.appendChild(li);
  }
}

function initBrowser() {
  const modal = document.getElementById('browse-modal');
  const close = () => { modal.style.display = 'none'; };

  document.getElementById('browse-btn').addEventListener('click', async () => {
    modal.style.display = 'flex';
    const current = document.getElementById('directory').value.trim();
    if (current) {
      await _browseLoad(current);
    } else {
      // Default to server launch cwd.
      try {
        const h = await api('/api/health');
        await _browseLoad(h.cwd || '');
      } catch (_) {
        await _browseLoad('');
      }
    }
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
      enableToolbarButtons();
      await refreshTableAndPanels();
    }
  } catch (_) { /* server not ready yet, ignore */ }
}

main().catch((err) => {
  console.error(err);
  setStatus(String(err), true);
});
