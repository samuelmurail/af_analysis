import { api, setStatus } from './api.js';
import { state, renderTable, renderSelectedResidues } from './table.js';
import { renderPlot, renderPaePlot } from './plot.js';
import { loadStructure, highlightResidue, highlightResidues, highlightTwoGroups } from './molstar.js';
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
    } catch (e) {
      renderPlot(plddtPayload, makePlddtHandlers());
      document.getElementById('plot-type').value = 'plddt';
    }
  } else {
    renderPlot(plddtPayload, makePlddtHandlers());
  }

  renderSelectedResidues();
  await loadStructure(state.selectedModel);
}

function makePaeHandlers() {
  return {
    onClick: (residue) => {
      const el = document.getElementById("selected-residues");
      if (el) el.textContent = `Scored: [${residue}] | Aligned: []`;
      highlightResidue(residue);
    },
    onPaeSelect: ({ xResidues, yResidues }) => {
      const el = document.getElementById("selected-residues");
      if (el) el.textContent = `Scored: [${xResidues.join(", ")}] | Aligned: [${yResidues.join(", ")}]`;
      if (xResidues.length > 0 || yResidues.length > 0) {
        highlightTwoGroups(xResidues, yResidues);
      }
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
}

async function main() {
  initResizableLayout();
  initEvents();
  renderSelectedResidues();
}

main().catch((err) => {
  console.error(err);
  setStatus(String(err), true);
});
