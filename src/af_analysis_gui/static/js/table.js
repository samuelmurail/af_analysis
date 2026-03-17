export const state = {
  selectedModel: 0,
  selectedResidues: [],
  paeSelection: null,   // { xResidues, yResidues } | null
  chainIds: [],
  chainLengths: [],
  viewer: null,
  plugin: null,
  tableRows: []
};

export function renderSelectedResidues() {
  const text = state.selectedResidues.length
    ? `Selected residues: ${state.selectedResidues.join(", ")}`
    : "Selected residues: []";
  const el = document.getElementById("selected-residues");
  if (el) el.textContent = text;
}

export function renderTable(columns, rows, onRowClick) {
  const head = document.getElementById("table-head");
  const body = document.getElementById("table-body");
  if (!head || !body) return;

  head.innerHTML = `<tr>${columns.map((c) => `<th>${c}</th>`).join("")}</tr>`;
  body.innerHTML = "";

  rows.forEach((row) => {
    const tr = document.createElement("tr");
    if (Number(row.row) === Number(state.selectedModel)) tr.classList.add("active");
    tr.innerHTML = columns.map((c) => `<td>${row[c] ?? ""}</td>`).join("");
    tr.addEventListener("click", () => onRowClick(Number(row.row), columns, rows));
    body.appendChild(tr);
  });
}
