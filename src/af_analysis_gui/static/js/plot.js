export function renderPlot(plddtPayload, handlers) {
  const trace = {
    x: plddtPayload.residues,
    y: plddtPayload.plddt,
    type: "scatter",
    mode: "lines+markers",
    marker: { size: 5 },
    line: { width: 2 },
    name: "pLDDT"
  };

  const shapes = (plddtPayload.chain_boundaries || []).map((x) => ({
    type: "line",
    x0: x + 0.5,
    x1: x + 0.5,
    y0: 0,
    y1: 100,
    line: { color: "#333", width: 1 }
  }));

  const layout = {
    margin: { l: 40, r: 10, t: 10, b: 40 },
    xaxis: { title: "Residue" },
    yaxis: { title: "pLDDT", range: [0, 100] },
    clickmode: "event+select",
    dragmode: "pan",
    shapes
  };

  const selEl = document.getElementById("selected-residues");
  if (selEl) selEl.style.display = "";

  // Store chain-boundary shapes so highlightPlotResidues can append hover
  // lines to them without needing to read custom properties back from Plotly.
  _plddtChainShapes = shapes;

  Plotly.newPlot("plddt-plot", [trace], layout, { responsive: true });

  const plotDiv = document.getElementById("plddt-plot");
  plotDiv.on("plotly_click", (ev) => {
    if (!ev?.points?.length) return;
    const residue = Number(ev.points[0].x);
    if (!Number.isInteger(residue)) return;
    if (handlers?.onClick) handlers.onClick(residue);
  });

  plotDiv.on("plotly_selected", (ev) => {
    const pts = ev?.points || [];
    const residues = Array.from(new Set(pts.map((p) => Number(p.x)).filter((x) => Number.isInteger(x)))).sort(
      (a, b) => a - b
    );
    if (handlers?.onSelect) handlers.onSelect(residues);
  });
}

// Base shapes for the pLDDT scatter plot (chain boundaries).
// Stored here so highlightPlotResidues never needs to read them back from
// _fullLayout.shapes (where Plotly strips unknown custom properties).
let _plddtChainShapes = [];

// Debounce timer for hover-driven plot updates.
let _hoverHighlightTimer = null;

// Module-level state for the PAE overlay so helpers can access it without
// closing over the stale `paePayload` after a re-render.
let _paeState = null;
let _paeOverlayActive = false; // guard against recursive plotly_deselect

const PAE_GREEN  = "#22c55e";  // scored residues  (x-axis strip, y = 0)
const PAE_ORANGE = "#f97316";  // aligned residues (y-axis strip, x = 0)
const PAE_PINK   = "#ec4899";  // residues in both groups

// Uses SVG shapes (always fully opaque) instead of scatter traces to mark the
// selected residue strips along the top (scored) and left (aligned) edges.
function _applyPaeOverlay(xResidues, yResidues) {
  if (!_paeState) return;
  const { payload, trace, layout } = _paeState;
  const n = payload.residues.length;
  const xSet = new Set(xResidues);
  const ySet = new Set(yResidues);

  // Keep existing chain-boundary lines, then append selection shapes.
  const selectionShapes = [];

  // White background strips so shapes are visible outside the heatmap.
  selectionShapes.push(
    { type: "rect", x0: 0.5, x1: n + 0.5, y0: -0.5, y1: -2,
      fillcolor: "#ffffff", line: { width: 0 }, layer: "above" },
    { type: "rect", x0: -0.5, x1: -2, y0: 0.5, y1: n + 0.5,
      fillcolor: "#ffffff", line: { width: 0 }, layer: "above" },
  );

  // One rectangle per scored residue along the top strip (y ∈ [-0.5, -2]).
  for (const r of xResidues) {
    selectionShapes.push({
      type: "rect", x0: r - 0.5, x1: r + 0.5, y0: -0.5, y1: -2,
      fillcolor: ySet.has(r) ? PAE_PINK : PAE_GREEN,
      line: { width: 0 }, layer: "above",
    });
  }

  // One rectangle per aligned residue along the left strip (x ∈ [-0.5, -2]).
  for (const r of yResidues) {
    selectionShapes.push({
      type: "rect", x0: -0.5, x1: -2, y0: r - 0.5, y1: r + 0.5,
      fillcolor: xSet.has(r) ? PAE_PINK : PAE_ORANGE,
      line: { width: 0 }, layer: "above",
    });
  }

  const overlayLayout = {
    ...layout,
    shapes: [...(layout.shapes || []), ...selectionShapes],
    xaxis: { ...layout.xaxis, range: [-2, n + 0.5], autorange: false },
    yaxis: { ...layout.yaxis, range: [n + 0.5, -2], autorange: false },
  };

  _paeOverlayActive = true;
  Plotly.react("plddt-plot", [trace], overlayLayout);
  _paeOverlayActive = false;
}

function _clearPaeOverlay() {
  if (!_paeState || _paeOverlayActive) return;
  const { trace, layout } = _paeState;
  Plotly.react("plddt-plot", [trace], layout);
}

export function renderPaePlot(paePayload, handlers) {
  const n = paePayload.residues.length;
  const boundaries = paePayload.chain_boundaries || [];

  const trace = {
    z: paePayload.pae,
    x: paePayload.residues,
    y: paePayload.residues,
    type: "heatmap",
    colorscale: [
      [0,   "#0053D6"],
      [0.1, "#0053D6"],
      [0.3, "#65CBF3"],
      [0.5, "#FFDB13"],
      [1.0, "#FF7D45"],
    ],
    reversescale: false,
    zmin: 0,
    zmax: 30,
    colorbar: { title: "PAE (Å)", thickness: 14, len: 0.8 },
    hoverongaps: false,
  };

  // Solid white lines at chain boundaries
  const shapes = [];
  for (const b of boundaries) {
    const pos = b + 0.5;
    shapes.push(
      { type: "line", x0: pos, x1: pos, y0: 0.5, y1: n + 0.5, line: { color: "#fff", width: 2 } },
      { type: "line", x0: 0.5, x1: n + 0.5, y0: pos, y1: pos, line: { color: "#fff", width: 2 } }
    );
  }

  const layout = {
    margin: { l: 50, r: 10, t: 10, b: 50 },
    xaxis: { title: "Scored residue", scaleanchor: "y", scaleratio: 1, constrain: "domain" },
    yaxis: { title: "Aligned residue", autorange: "reversed", constrain: "domain" },
    clickmode: "event+select",
    dragmode: "select",
    selectdirection: "any",
    shapes,
  };

  // Save state so overlay helpers can access payload/trace/layout later.
  _paeState = { payload: paePayload, trace, layout };
  _paeOverlayActive = false;

  const selEl = document.getElementById("selected-residues");
  if (selEl) selEl.style.display = "";

  Plotly.newPlot("plddt-plot", [trace], layout, { responsive: true });

  const plotDiv = document.getElementById("plddt-plot");

  plotDiv.on("plotly_click", (ev) => {
    _clearPaeOverlay();
    if (!ev?.points?.length) return;
    const residue = Number(ev.points[0].x);
    if (!Number.isInteger(residue)) return;
    if (handlers?.onClick) handlers.onClick(residue);
  });

  // Heatmap traces don't populate ev.points on plotly_selected.
  // Use ev.range (coordinate bounds of the drag box) and filter the residue list.
  plotDiv.on("plotly_selected", (ev) => {
    if (!ev?.range) return;
    const [x0, x1] = (ev.range.x || []).map(Number);
    const [y0, y1] = (ev.range.y || []).map(Number);
    const xMin = Math.min(x0, x1), xMax = Math.max(x0, x1);
    const yMin = Math.min(y0, y1), yMax = Math.max(y0, y1);
    const residues = paePayload.residues || [];
    const xResidues = residues.filter(r => r >= xMin && r <= xMax);
    const yResidues = residues.filter(r => r >= yMin && r <= yMax);
    _applyPaeOverlay(xResidues, yResidues);
    if (handlers?.onPaeSelect) handlers.onPaeSelect({ xResidues, yResidues });
  });

  plotDiv.on("plotly_deselect", () => {
    _clearPaeOverlay();
    if (handlers?.onDeselect) handlers.onDeselect();
  });
}

export function resizePlot() {
  const plotDiv = document.getElementById("plddt-plot");
  if (plotDiv && typeof Plotly !== "undefined" && Plotly.Plots?.resize) {
    Plotly.Plots.resize(plotDiv);
  }
}

// Highlight (hover) residues on the active Plotly plot driven by Mol* hover events.
// For pLDDT: draws a vertical dashed orange line at the hovered x position.
// For PAE: triggers Plotly's native crosshair/tooltip at the hovered cell.
// Debounced to 30 ms so rapid Mol* hover events don't flood Plotly.
export function highlightPlotResidues(globalResidues) {
  if (_hoverHighlightTimer) { clearTimeout(_hoverHighlightTimer); }
  _hoverHighlightTimer = setTimeout(() => {
    _doHighlightPlotResidues(globalResidues);
  }, 30);
}

function _doHighlightPlotResidues(globalResidues) {
  const plotDiv = document.getElementById("plddt-plot");
  if (!plotDiv?._fullLayout) return;

  const traceType = plotDiv._fullData?.[0]?.type;

  if (traceType === 'scatter') {
    // _plddtChainShapes is set by renderPlot and never includes hover lines,
    // so we can safely append/remove without relying on Plotly preserving any
    // custom marker property on shapes.
    if (!globalResidues?.length) {
      // Remove hover lines by restoring chain-boundary-only shapes.
      const cur = plotDiv._fullLayout.shapes || [];
      if (cur.length > _plddtChainShapes.length) {
        Plotly.relayout('plddt-plot', { shapes: _plddtChainShapes });
      }
      return;
    }
    const hoverShapes = globalResidues.map(r => ({
      type: 'line',
      x0: r, x1: r, y0: 0, y1: 100,
      line: { color: '#f97316', width: 2, dash: 'dot' },
      layer: 'above',
    }));
    Plotly.relayout('plddt-plot', { shapes: [..._plddtChainShapes, ...hoverShapes] });
  } else if (traceType === 'heatmap') {
    // For the PAE heatmap use Plotly's lightweight native hover to show the
    // crosshair + tooltip without touching the SVG-shape overlay.
    if (!globalResidues?.length) {
      Plotly.Fx.unhover('plddt-plot');
      return;
    }
    const r = globalResidues[0];
    const xs = plotDiv._fullData[0]?.x || [];
    const ys = plotDiv._fullData[0]?.y || [];
    const xi = xs.indexOf(r);
    const yi = ys.indexOf(r);
    if (xi !== -1 && yi !== -1) {
      Plotly.Fx.hover('plddt-plot', [{ curveNumber: 0, pointNumber: xi + yi * xs.length }]);
    }
  }
}
