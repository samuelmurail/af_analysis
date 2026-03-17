export const state = {
  selectedModel: 0,
  selectedResidues: [],
  paeSelection: null,   // { xResidues, yResidues } | null
  chainIds: [],
  chainLengths: [],
  viewer: null,
  plugin: null,
  tableRows: [],
  sortCol: null,  // column name currently sorted, or null
  sortDir: 1,     // 1 = ascending, -1 = descending
};

export function renderSelectedResidues() {
  const text = state.selectedResidues.length
    ? `Selected residues: ${state.selectedResidues.join(", ")}`
    : "Selected residues: []";
  const el = document.getElementById("selected-residues");
  if (el) el.textContent = text;
}

function _fmtCell(value) {
  if (typeof value === "number" && Number.isFinite(value)) {
    // Use 4 significant figures; strip trailing zeros after the decimal point.
    return parseFloat(value.toPrecision(4)).toString();
  }
  return value ?? "";
}

export function renderTable(columns, rows, onRowClick) {
  const head = document.getElementById("table-head");
  const body = document.getElementById("table-body");
  if (!head || !body) return;

  // Build sortable headers.
  const headerCells = columns.map((c) => {
    const isActive = c === state.sortCol;
    const arrow = isActive ? (state.sortDir === 1 ? ' ▲' : ' ▼') : '';
    return `<th class="sortable${isActive ? ' sort-active' : ''}" data-col="${c}">${c}${arrow}</th>`;
  }).join("");
  head.innerHTML = `<tr>${headerCells}</tr>`;

  // Attach click handlers to each <th>.
  head.querySelectorAll('th.sortable').forEach((th) => {
    th.addEventListener('click', () => {
      const col = th.dataset.col;
      if (state.sortCol === col) {
        state.sortDir *= -1;
      } else {
        state.sortCol = col;
        state.sortDir = 1;
      }
      renderTable(columns, rows, onRowClick);
    });
  });

  // Sort a shallow copy of rows so the original order is preserved.
  const sorted = state.sortCol
    ? [...rows].sort((a, b) => {
        const av = a[state.sortCol];
        const bv = b[state.sortCol];
        const an = Number(av), bn = Number(bv);
        if (!isNaN(an) && !isNaN(bn)) return (an - bn) * state.sortDir;
        return String(av ?? '').localeCompare(String(bv ?? '')) * state.sortDir;
      })
    : rows;

  body.innerHTML = "";
  sorted.forEach((row) => {
    const tr = document.createElement("tr");
    if (Number(row.row) === Number(state.selectedModel)) tr.classList.add("active");
    tr.innerHTML = columns.map((c) => `<td>${_fmtCell(row[c])}</td>`).join("");
    tr.addEventListener("click", () => onRowClick(Number(row.row), columns, rows));
    body.appendChild(tr);
  });
}
