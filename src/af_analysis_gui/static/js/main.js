import { api, setStatus } from './api.js';
import { state, renderTable, renderSelectedResidues } from './table.js';
import { renderPlot } from './plot.js';
import { loadStructure, highlightResidue, highlightResidues } from './molstar.js';
import { initResizableLayout } from './resize.js';

function renderCurrentTable(columns, rows) {
  renderTable(columns, rows, async (selectedRow) => {
    state.selectedModel = Number(selectedRow);
    renderCurrentTable(columns, rows);
    await refreshModelPanels();
  });
}

async function refreshModelPanels() {
  const plddtPayload = await api(`/api/plddt?index=${state.selectedModel}`);
  state.chainIds = plddtPayload.chain_ids || [];
  state.chainLengths = plddtPayload.chain_lengths || [];
  renderPlot(plddtPayload, {
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
  });
  renderSelectedResidues();
  await loadStructure(state.selectedModel);
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
