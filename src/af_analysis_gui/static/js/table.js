import { api } from './api.js';

// ── Row detail modal ────────────────────────────────────────────────────────
function _openRowDetail(rowIdx) {
  const modal = document.getElementById('row-detail-modal');
  const title = document.getElementById('row-detail-title');
  const tbody = document.getElementById('row-detail-body');
  if (!modal || !tbody) return;
  title.textContent = `Row ${rowIdx} — details`;
  tbody.innerHTML = '<tr><td colspan="2" style="color:#888; padding:12px;">Loading…</td></tr>';
  modal.style.display = 'flex';
  api(`/api/row/${rowIdx}`).then(({ data }) => {
    tbody.innerHTML = Object.entries(data)
      .map(([k, v]) => {
        const display = v === null ? '<em style="color:#aaa">null</em>'
          : typeof v === 'string' && v.startsWith('[array,')
            ? `<em style="color:#888">${_escHtml(v)}</em>`
            : _escHtml(String(v));
        return `<tr><th>${_escHtml(k)}</th><td>${display}</td></tr>`;
      }).join('');
  }).catch(() => {
    tbody.innerHTML = '<tr><td colspan="2" style="color:#a11927;">Failed to load row data.</td></tr>';
  });
}

function _escHtml(s) {
  return s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
}

function _initRowDetailModal() {
  const modal = document.getElementById('row-detail-modal');
  if (!modal || modal._afDetailInited) return;
  modal._afDetailInited = true;
  document.getElementById('row-detail-close')?.addEventListener('click', () => {
    modal.style.display = 'none';
  });
  modal.addEventListener('click', (e) => {
    if (e.target === modal) modal.style.display = 'none';
  });
}

// Module-level storage for the current render options (cellStyleFn, etc.).
// Kept up-to-date on every renderTable call so that sort/collapse re-renders
// always use the latest options rather than stale closures.
let _currentOptions = {};

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
  filters: {},  // { columnName: { min, max } | string }
  hiddenColumns: new Set(),  // user-hidden columns
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

function _escAttr(s) {
  return String(s).replace(/&/g, '&amp;').replace(/"/g, '&quot;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
}

function _matchesFilter(value, filter) {
  if (!filter) return true;
  if (typeof filter === 'string') {
    if (!filter.trim()) return true;
    return String(value ?? '').toLowerCase().includes(filter.trim().toLowerCase());
  }
  // Object: { min, max } for numeric range
  const { min, max } = filter;
  if ((!min || !min.trim()) && (!max || !max.trim())) return true;
  const numVal = parseFloat(value);
  if (isNaN(numVal)) return false;
  if (min && min.trim()) { const n = parseFloat(min); if (!isNaN(n) && numVal < n) return false; }
  if (max && max.trim()) { const n = parseFloat(max); if (!isNaN(n) && numVal > n) return false; }
  return true;
}

// Column type detection cache (reset when rows array reference changes).
let _colTypeCache = {};
let _colTypeCacheRows = null;

function _getColumnType(colName, rows) {
  if (_colTypeCacheRows !== rows) {
    _colTypeCache = {};
    _colTypeCacheRows = rows;
  }
  if (colName in _colTypeCache) return _colTypeCache[colName];
  let numCount = 0, total = 0;
  for (const row of rows) {
    const v = row[colName];
    if (v === null || v === undefined || v === '') continue;
    total++;
    if (typeof v === 'number' && isFinite(v)) numCount++;
  }
  const type = total > 0 && numCount / total > 0.5 ? 'numeric' : 'string';
  _colTypeCache[colName] = type;
  return type;
}

let _filterDebounce = null;

export function renderTable(columns, rows, onRowClick, options = {}) {
  _currentOptions = options;  // always keep the latest options available for re-renders
  const { cellStyleFn = null } = _currentOptions;
  const hiddenColumns = state.hiddenColumns;
  const head = document.getElementById("table-head");
  const body = document.getElementById("table-body");
  if (!head || !body) return;

  // Save focus state for filter restore.
  const focusedFilterCol = document.activeElement?.dataset?.filterCol;
  const focusedFilterBound = document.activeElement?.dataset?.filterBound;
  const focusedCursorPos = document.activeElement?.selectionStart;

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
    if (hiddenColumns.has(col)) continue;
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
          for (const mc of groupCols[g.id]) {
            if (!hiddenColumns.has(mc)) displayColumns.push(mc);
          }
        }
      }
      // Collapsed members are intentionally skipped (seenGroupIds already set).
    }
  }

  _initRowDetailModal();

  // ── header ────────────────────────────────────────────────────────────────
  const firstMemberSeen = new Set();
  const headerCells = ['<th class="row-info-cell"></th>', ...displayColumns.map((dc) => {
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
  })];
  head.innerHTML = `<tr>${headerCells.join("")}</tr>`;

  // Group-toggle click handlers.
  head.querySelectorAll('button.group-toggle').forEach((btn) => {
    btn.addEventListener('click', (e) => {
      e.stopPropagation();
      const gid = btn.dataset.group;
      if (state.collapsedGroups.has(gid)) state.collapsedGroups.delete(gid);
      else state.collapsedGroups.add(gid);
      renderTable(columns, rows, onRowClick, _currentOptions);
    });
  });

  // Sort click handlers.
  head.querySelectorAll('th.sortable').forEach((th) => {
    th.addEventListener('click', () => {
      const col = th.dataset.col;
      if (state.sortCol === col) state.sortDir *= -1;
      else { state.sortCol = col; state.sortDir = 1; }
      renderTable(columns, rows, onRowClick, _currentOptions);
    });
  });

  // ── filter row ────────────────────────────────────────────────────────────
  const filterCells = ['<th></th>', ...displayColumns.map((dc) => {
    const colName = typeof dc === 'object' ? dc.__group__ : dc;
    const colType = _getColumnType(colName, rows);
    if (colType === 'numeric') {
      const filter = state.filters[colName];
      const minVal = (filter && typeof filter === 'object') ? (filter.min || '') : '';
      const maxVal = (filter && typeof filter === 'object') ? (filter.max || '') : '';
      return `<th><div class="filter-range"><input type="text" class="filter-min" data-filter-col="${_escAttr(colName)}" data-filter-bound="min" value="${_escAttr(minVal)}" placeholder="Min"><span class="filter-sep">&ndash;</span><input type="text" class="filter-max" data-filter-col="${_escAttr(colName)}" data-filter-bound="max" value="${_escAttr(maxVal)}" placeholder="Max"></div></th>`;
    } else {
      const val = typeof state.filters[colName] === 'string' ? state.filters[colName] : '';
      return `<th><input type="text" class="filter-input" data-filter-col="${_escAttr(colName)}" value="${_escAttr(val)}" placeholder="Search\u2026"></th>`;
    }
  })];
  const filterRow = document.createElement('tr');
  filterRow.className = 'filter-row';
  filterRow.innerHTML = filterCells.join('');
  head.appendChild(filterRow);

  // Text filter handlers.
  filterRow.querySelectorAll('.filter-input').forEach(input => {
    input.addEventListener('input', () => {
      state.filters[input.dataset.filterCol] = input.value;
      clearTimeout(_filterDebounce);
      _filterDebounce = setTimeout(() => {
        renderTable(columns, rows, onRowClick, _currentOptions);
      }, 200);
    });
  });

  // Numeric min/max filter handlers.
  filterRow.querySelectorAll('.filter-min, .filter-max').forEach(input => {
    input.addEventListener('input', () => {
      const col = input.dataset.filterCol;
      const bound = input.dataset.filterBound;
      let filter = state.filters[col];
      if (!filter || typeof filter !== 'object') filter = { min: '', max: '' };
      filter[bound] = input.value;
      state.filters[col] = filter;
      clearTimeout(_filterDebounce);
      _filterDebounce = setTimeout(() => {
        renderTable(columns, rows, onRowClick, _currentOptions);
      }, 200);
    });
  });

  // ── column picker ─────────────────────────────────────────────────────────
  const pickerMenu = document.getElementById('col-picker-menu');
  const pickerWrap = document.getElementById('col-picker-wrap');
  if (pickerMenu && pickerWrap) {
    pickerWrap.style.display = '';
    pickerMenu.innerHTML = columns
      .filter(c => c !== 'row')
      .map(c => {
        const checked = !state.hiddenColumns.has(c) ? 'checked' : '';
        const id = `pick-${_escAttr(c)}`;
        return `<div class="form-check"><input class="form-check-input" type="checkbox" ${checked} data-pick-col="${_escAttr(c)}" id="${id}"><label class="form-check-label" for="${id}">${_escHtml(c)}</label></div>`;
      }).join('');
    pickerMenu.querySelectorAll('.form-check-input').forEach(cb => {
      cb.addEventListener('change', () => {
        if (cb.checked) state.hiddenColumns.delete(cb.dataset.pickCol);
        else state.hiddenColumns.add(cb.dataset.pickCol);
        renderTable(columns, rows, onRowClick, _currentOptions);
      });
    });
  }

  // ── body ──────────────────────────────────────────────────────────────────
  const sorted = state.sortCol
    ? [...rows].sort((a, b) => {
        const av = a[state.sortCol], bv = b[state.sortCol];
        const an = Number(av), bn = Number(bv);
        if (!isNaN(an) && !isNaN(bn)) return (an - bn) * state.sortDir;
        return String(av ?? '').localeCompare(String(bv ?? '')) * state.sortDir;
      })
    : rows;

  // Apply filters.
  let filtered = sorted;
  for (const [col, expr] of Object.entries(state.filters)) {
    const hasFilter = typeof expr === 'string'
      ? expr.trim()
      : (expr && (expr.min?.trim() || expr.max?.trim()));
    if (hasFilter) {
      filtered = filtered.filter(row => _matchesFilter(row[col], expr));
    }
  }

  // Update table info.
  const infoEl = document.getElementById('table-info');
  if (infoEl) {
    const total = rows.length;
    const shown = filtered.length;
    infoEl.textContent = total === shown ? `${total} models` : `Showing ${shown} of ${total} models`;
  }

  body.innerHTML = "";
  filtered.forEach((row) => {
    const tr = document.createElement("tr");
    const rowNum = Number(row.row);
    if (rowNum === Number(state.selectedModel)) tr.classList.add("active");
    if (state.selectedRows.size > 0 && state.selectedRows.has(rowNum)) tr.classList.add("superpose-row");
    if (state.hoveredRow === rowNum) tr.classList.add("hover-row");
    const infoBtn = document.createElement('td');
    infoBtn.className = 'row-info-cell';
    infoBtn.innerHTML = `<button class="row-info-btn" title="Show all details">i</button>`;
    infoBtn.querySelector('button').addEventListener('click', (e) => {
      e.stopPropagation();
      _openRowDetail(rowNum);
    });
    tr.appendChild(infoBtn);
    const dataCells = document.createElement('template');
    dataCells.innerHTML = displayColumns.map((dc) => {
      if (typeof dc === 'object') {
        // Collapsed group: show max of member columns.
        let maxVal = null;
        for (const mc of dc.memberCols) {
          const v = row[mc];
          if (typeof v === 'number' && isFinite(v) && (maxVal === null || v > maxVal)) maxVal = v;
        }
        return `<td class="group-cell">${maxVal !== null ? _fmtCell(maxVal) : ''}</td>`;
      }
      const style = cellStyleFn ? (cellStyleFn(dc, row) ?? '') : '';
      const styleAttr = style ? ` style="${style}"` : '';
      const g = _groupFor(dc);
      return g
        ? `<td class="group-cell"${styleAttr}>${_fmtCell(row[dc])}</td>`
        : `<td${styleAttr}>${_fmtCell(row[dc])}</td>`;
    }).join("");
    tr.appendChild(dataCells.content);
    tr.addEventListener("click", () => onRowClick(Number(row.row), columns, rows));
    body.appendChild(tr);
  });

  // Restore filter focus.
  if (focusedFilterCol) {
    let input;
    if (focusedFilterBound) {
      input = head.querySelector(`.filter-${focusedFilterBound}[data-filter-col="${CSS.escape(focusedFilterCol)}"]`);
    } else {
      input = head.querySelector(`.filter-input[data-filter-col="${CSS.escape(focusedFilterCol)}"]`);
    }
    if (input) {
      input.focus();
      if (focusedCursorPos != null) input.setSelectionRange(focusedCursorPos, focusedCursorPos);
    }
  }
}
