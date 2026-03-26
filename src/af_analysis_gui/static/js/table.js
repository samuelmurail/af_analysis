export const state = {
  selectedModel: 0,
  selectedRows: new Set(),  // rows highlighted in superpose mode
  hoveredRow: null,          // row index hovered in Mol* superpose view
  selectedResidues: [],
  paeSelection: null,   // { xResidues, yResidues } | null
  chainIds: [],
  chainLengths: [],
  chainTypes: [],
  viewer: null,
  plugin: null,
  tableRows: [],
  sortCol: null,  // column name currently sorted, or null
  sortDir: 1,     // 1 = ascending, -1 = descending
  collapsedGroups: new Set(['pdockq2', 'iptm_d0', 'ipsae']),  // group ids collapsed by default
};

// Collapsible column groups: columns whose name starts with `prefix`
// (and isn't in `exclude`) are folded into a single summary column.
const COLUMN_GROUPS = [
  {
    id: 'pdockq2',
    prefix: 'pdockq2_',
    exclude: new Set(),
    label: 'pdockq2',
  },
  {
    id: 'iptm_d0',
    prefix: 'ipTM_d0_',
    exclude: new Set(['ipTM_d0_matrix']),
    label: 'ipTM_d0',
  },
  {
    id: 'ipsae',
    prefix: 'ipSAE_',
    exclude: new Set(['ipSAE_matrix']),
    label: 'ipSAE',
  },
];

/** Return the group config for a column name, or null. */
function _groupFor(colName) {
  for (const g of COLUMN_GROUPS) {
    if (colName.startsWith(g.prefix) && !g.exclude.has(colName)) return g;
  }
  return null;
}

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

  // ── collect group member columns present in this dataset ─────────────────
  const groupCols = {};   // groupId → [colName, ...]
  for (const g of COLUMN_GROUPS) {
    groupCols[g.id] = columns.filter(c => c.startsWith(g.prefix) && !g.exclude.has(c));
  }

  // ── build display column list ─────────────────────────────────────────────
  // Each entry is either:
  //   { __group__, label, memberCols }  — collapsed-group placeholder (1 th / 1 td)
  //   string                            — regular column (1 th / 1 td)
  //
  // When a group is expanded, its member columns appear as normal strings.
  // The first member column's <th> carries the collapse button.
  const seenGroupIds = new Set();
  const displayColumns = [];

  for (const col of columns) {
    const g = _groupFor(col);
    if (!g) {
      displayColumns.push(col);
    } else {
      if (!seenGroupIds.has(g.id)) {
        seenGroupIds.add(g.id);
        if (state.collapsedGroups.has(g.id)) {
          // One placeholder for the whole group.
          displayColumns.push({ __group__: g.id, label: g.label, memberCols: groupCols[g.id] });
        } else {
          // All member columns individually.
          for (const mc of groupCols[g.id]) displayColumns.push(mc);
        }
      }
      // Collapsed members are intentionally skipped (seenGroupIds already set).
    }
  }

  // ── header ────────────────────────────────────────────────────────────────
  const firstMemberSeen = new Set();
  const headerCells = displayColumns.map((dc) => {
    if (typeof dc === 'object') {
      // Collapsed group — single summary header with expand button.
      const btn = `<button class="group-toggle" data-group="${dc.__group__}" title="Expand">${dc.label} ▶</button>`;
      return `<th class="group-header">${btn}</th>`;
    }

    const g = _groupFor(dc);
    const isActive = dc === state.sortCol;
    const arrow = isActive ? (state.sortDir === 1 ? ' ▲' : ' ▼') : '';
    const sortAttrs = `class="sortable${isActive ? ' sort-active' : ''}" data-col="${dc}"`;

    if (g && !firstMemberSeen.has(g.id)) {
      // First column of an expanded group — add collapse button.
      firstMemberSeen.add(g.id);
      const btn = `<button class="group-toggle" data-group="${g.id}" title="Collapse">◀</button>`;
      return `<th ${sortAttrs}>${dc}${arrow} ${btn}</th>`;
    }

    return `<th ${sortAttrs}>${dc}${arrow}</th>`;
  });
  head.innerHTML = `<tr>${headerCells.join("")}</tr>`;

  // Group-toggle click handlers.
  head.querySelectorAll('button.group-toggle').forEach((btn) => {
    btn.addEventListener('click', (e) => {
      e.stopPropagation();
      const gid = btn.dataset.group;
      if (state.collapsedGroups.has(gid)) state.collapsedGroups.delete(gid);
      else state.collapsedGroups.add(gid);
      renderTable(columns, rows, onRowClick);
    });
  });

  // Sort click handlers.
  head.querySelectorAll('th.sortable').forEach((th) => {
    th.addEventListener('click', () => {
      const col = th.dataset.col;
      if (state.sortCol === col) state.sortDir *= -1;
      else { state.sortCol = col; state.sortDir = 1; }
      renderTable(columns, rows, onRowClick);
    });
  });

  // ── body ──────────────────────────────────────────────────────────────────
  const sorted = state.sortCol
    ? [...rows].sort((a, b) => {
        const av = a[state.sortCol], bv = b[state.sortCol];
        const an = Number(av), bn = Number(bv);
        if (!isNaN(an) && !isNaN(bn)) return (an - bn) * state.sortDir;
        return String(av ?? '').localeCompare(String(bv ?? '')) * state.sortDir;
      })
    : rows;

  body.innerHTML = "";
  sorted.forEach((row) => {
    const tr = document.createElement("tr");
    const rowNum = Number(row.row);
    if (rowNum === Number(state.selectedModel)) tr.classList.add("active");
    if (state.selectedRows.size > 0 && state.selectedRows.has(rowNum)) tr.classList.add("superpose-row");
    if (state.hoveredRow === rowNum) tr.classList.add("hover-row");
    tr.innerHTML = displayColumns.map((dc) => {
      if (typeof dc === 'object') {
        // Collapsed group: show max of member columns.
        let maxVal = null;
        for (const mc of dc.memberCols) {
          const v = row[mc];
          if (typeof v === 'number' && isFinite(v) && (maxVal === null || v > maxVal)) maxVal = v;
        }
        return `<td class="group-cell">${maxVal !== null ? _fmtCell(maxVal) : ''}</td>`;
      }
      const g = _groupFor(dc);
      return g
        ? `<td class="group-cell">${_fmtCell(row[dc])}</td>`
        : `<td>${_fmtCell(row[dc])}</td>`;
    }).join("");
    tr.addEventListener("click", () => onRowClick(Number(row.row), columns, rows));
    body.appendChild(tr);
  });
}
