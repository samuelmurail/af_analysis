export function renderPlot(plddtPayload, handlers, savedZoom = null) {
  const plotlyConfig = { responsive: true };

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
    margin: { l: 40, r: 10, t: 38, b: 40 },
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

  Plotly.newPlot("plddt-plot", [trace], layout, plotlyConfig).then(() => {
    if (savedZoom?.xRange) {
      Plotly.relayout("plddt-plot", { 'xaxis.range': savedZoom.xRange, 'xaxis.autorange': false });
    }
  });
  document.getElementById('plddt-plot')?.classList.remove('pae-active');

  const plotDiv = document.getElementById("plddt-plot");

  plotDiv.on("plotly_click", (ev) => {
    const residue = Number(ev.points[0].x);
    if (!Number.isInteger(residue)) return;
    if (handlers?.onClick) handlers.onClick(residue);
  });

  plotDiv.on("plotly_hover", (ev) => {
    if (!ev?.points?.length) return;
    const residue = Number(ev.points[0].x);
    if (!Number.isInteger(residue)) return;
    if (handlers?.onHover) handlers.onHover([residue]);
  });

  plotDiv.on("plotly_unhover", () => {
    if (handlers?.onUnhover) handlers.onUnhover();
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
let _paeActiveBaseShapes = []; // chain-boundary lines (or chain + selection overlay shapes) currently rendered
let _paeStrip = 2; // strip width (residue units) computed at renderPaePlot time
let _paeHoverLines = [];       // transient vertical+horizontal crosshair added on top of base shapes

const PAE_GREEN  = "#22c55e";  // scored residues  (x-axis strip, y = 0)
const PAE_ORANGE = "#f97316";  // aligned residues (y-axis strip, x = 0)
const PAE_PINK   = "#ec4899";  // residues in both groups

// Uses SVG shapes (always fully opaque) instead of scatter traces to mark the
// selected residue strips along the top (scored) and left (aligned) edges.
function _applyPaeOverlay(xResidues, yResidues) {
  if (!_paeState) return;
  const { payload, trace, layout } = _paeState;
  const plotDiv = document.getElementById("plddt-plot");
  const curXRange = plotDiv?._fullLayout?.xaxis?.range?.slice() || layout.xaxis.range;
  const curYRange = plotDiv?._fullLayout?.yaxis?.range?.slice() || layout.yaxis.range;
  const n = payload.residues.length;
  const xSet = new Set(xResidues);
  const ySet = new Set(yResidues);

  // Re-use the strip width computed at renderPaePlot time so the overlay
  // geometry always matches the pre-reserved gutter in the initial layout.
  const STRIP = _paeStrip;
  const OUTER = -(STRIP + 0.5);  // outer edge of both strips

  // Keep existing chain-boundary lines, then append selection shapes.
  const selectionShapes = [];

  // White background strips so shapes are visible outside the heatmap.
  selectionShapes.push(
    { type: "rect", x0: 0.5, x1: n + 0.5, y0: -0.5, y1: OUTER,
      fillcolor: "#ffffff", line: { width: 0 }, layer: "above" },
    { type: "rect", x0: -0.5, x1: OUTER, y0: 0.5, y1: n + 0.5,
      fillcolor: "#ffffff", line: { width: 0 }, layer: "above" },
  );

  // One rectangle per scored residue along the top strip (y ∈ [-0.5, OUTER]).
  for (const r of xResidues) {
    selectionShapes.push({
      type: "rect", x0: r - 0.5, x1: r + 0.5, y0: -0.5, y1: OUTER,
      fillcolor: ySet.has(r) ? PAE_PINK : PAE_GREEN,
      line: { width: 0 }, layer: "above",
    });
  }

  // One rectangle per aligned residue along the left strip (x ∈ [OUTER, -0.5]).
  for (const r of yResidues) {
    selectionShapes.push({
      type: "rect", x0: -0.5, x1: OUTER, y0: r - 0.5, y1: r + 0.5,
      fillcolor: xSet.has(r) ? PAE_PINK : PAE_ORANGE,
      line: { width: 0 }, layer: "above",
    });
  }

  const overlayLayout = {
    ...layout,
    shapes: [...(layout.shapes || []), ...selectionShapes],
    xaxis: { ...layout.xaxis, range: curXRange, autorange: false },
    yaxis: { ...layout.yaxis, range: curYRange, autorange: false },
  };

  _paeActiveBaseShapes = [...(layout.shapes || []), ...selectionShapes];
  _paeHoverLines = [];
  _paeOverlayActive = true;
  Plotly.react("plddt-plot", [trace], overlayLayout);
  _paeOverlayActive = false;
}

function _clearPaeOverlay() {
  if (!_paeState || _paeOverlayActive) return;
  const { trace, layout } = _paeState;
  const plotDiv = document.getElementById("plddt-plot");
  const curXRange = plotDiv?._fullLayout?.xaxis?.range?.slice() || layout.xaxis.range;
  const curYRange = plotDiv?._fullLayout?.yaxis?.range?.slice() || layout.yaxis.range;
  _paeActiveBaseShapes = layout.shapes || [];
  _paeHoverLines = [];
  const preservedLayout = {
    ...layout,
    xaxis: { ...layout.xaxis, range: curXRange, autorange: false },
    yaxis: { ...layout.yaxis, range: curYRange, autorange: false },
  };
  Plotly.react("plddt-plot", [trace], preservedLayout);
}

export function renderPaePlot(paePayload, handlers, savedZoom = null) {
  // Heatmap traces filter out the built-in "select2d" string from
  // modeBarButtonsToAdd, so we must supply a full custom-button object.
  const plotlyConfig = {
    responsive: true,
    displayModeBar: true,
    modeBarButtonsToAdd: [
      {
        name: "Box Select",
        title: "Box Select",
        icon: Plotly.Icons.selectbox,
        click(gd) { Plotly.relayout(gd, { dragmode: "select" }); },
      },
      {
        name: "Pan",
        title: "Pan",
        icon: Plotly.Icons.pan,
        click(gd) { Plotly.relayout(gd, { dragmode: "pan" }); },
      },
    ],
  };

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

  // Compute strip width once and store it so _applyPaeOverlay uses the same value.
  _paeStrip = Math.max(2, Math.ceil(n * 0.013));
  const PAE_OUTER = -(_paeStrip + 0.5);

  const layout = {
    margin: { l: 50, r: 10, t: 38, b: 50 },
    // Pre-reserve space for the selection strips so the plot never resizes on
    // first selection — range already includes the strip gutter.
    xaxis: { title: "Scored residue", scaleanchor: "y", scaleratio: 1, constrain: "domain",
             range: [PAE_OUTER, n + 0.5], autorange: false },
    yaxis: { title: "Aligned residue", constrain: "domain",
             range: [n + 0.5, PAE_OUTER], autorange: false },
    clickmode: "event+select",
    dragmode: "select",
    selectdirection: "any",
    shapes,
  };

  // Save state so overlay helpers can access payload/trace/layout later.
  _paeState = { payload: paePayload, trace, layout };
  _paeOverlayActive = false;
  _paeActiveBaseShapes = shapes;
  _paeHoverLines = [];

  const selEl = document.getElementById("selected-residues");
  if (selEl) selEl.style.display = "";

  Plotly.newPlot("plddt-plot", [trace], layout, plotlyConfig).then(() => {
    if (savedZoom?.xRange && savedZoom?.yRange) {
      Plotly.relayout("plddt-plot", {
        'xaxis.range': savedZoom.xRange,
        'xaxis.autorange': false,
        'yaxis.range': savedZoom.yRange,
        'yaxis.autorange': false,
      });
    }
  });
  document.getElementById('plddt-plot')?.classList.add('pae-active');

  const plotDiv = document.getElementById("plddt-plot");

  plotDiv.on("plotly_click", (ev) => {
    _clearPaeOverlay();
    if (!ev?.points?.length) return;
    const residue = Number(ev.points[0].x);
    if (!Number.isInteger(residue)) return;
    if (handlers?.onClick) handlers.onClick(residue);
  });

  plotDiv.on("plotly_hover", (ev) => {
    if (!ev?.points?.length) return;
    const xResidue = Number(ev.points[0].x);
    const yResidue = Number(ev.points[0].y);
    if (handlers?.onHover) handlers.onHover({ xResidue, yResidue });
  });

  plotDiv.on("plotly_unhover", () => {
    if (handlers?.onUnhover) handlers.onUnhover();
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

export function renderLisPlot(lisPayload, handlers) {
  const chainIds = lisPayload.chain_ids || [];
  const lis = lisPayload.lis || [];
  const label = lisPayload.label || "LIS";

  const textMatrix = lis.map(row =>
    row.map(v => (typeof v === "number" && Number.isFinite(v)) ? v.toFixed(2) : "N/A")
  );

  const trace = {
    z: lis,
    x: chainIds,
    y: chainIds,
    type: "heatmap",
    colorscale: [
      [0.000, "#7E1A00"],
      [0.125, "#AE4A1A"],
      [0.250, "#D4896A"],
      [0.375, "#EED8C0"],
      [0.500, "#FAFAF8"],
      [0.625, "#C0CEDF"],
      [0.750, "#6898C8"],
      [0.875, "#2860A0"],
      [1.000, "#001880"],
    ],
    zmin: 0,
    zmax: 1,
    colorbar: { title: label, thickness: 14, len: 0.8 },
    text: textMatrix,
    texttemplate: "%{text}",
    textfont: { color: "black", size: 12 },
    hovertemplate: "Chain %{y} → Chain %{x}<br>" + label + ": %{text}<extra></extra>",
    hoverongaps: false,
  };

  const layout = {
    margin: { l: 50, r: 10, t: 38, b: 50 },
    xaxis: { title: "Chain j" },
    yaxis: { title: "Chain i", autorange: "reversed" },
  };

  Plotly.newPlot("plddt-plot", [trace], layout, { responsive: true });
  document.getElementById("plddt-plot")?.classList.remove("pae-active");

  const plotDiv = document.getElementById("plddt-plot");
  plotDiv.on("plotly_click", (ev) => {
    if (!ev?.points?.length) return;
    const xChainId = String(ev.points[0].x);
    const yChainId = String(ev.points[0].y);
    if (handlers?.onClick) handlers.onClick({ xChainId, yChainId });
  });

  plotDiv.on("plotly_hover", (ev) => {
    if (!ev?.points?.length) return;
    const xChainId = String(ev.points[0].x);
    const yChainId = String(ev.points[0].y);
    if (handlers?.onHover) handlers.onHover({ xChainId, yChainId });
  });

  plotDiv.on("plotly_unhover", () => {
    if (handlers?.onUnhover) handlers.onUnhover();
  });
}

export function resizePlot() {
  const plotDiv = document.getElementById("plddt-plot");
  if (plotDiv && typeof Plotly !== "undefined" && Plotly.Plots?.resize) {
    Plotly.Plots.resize(plotDiv);
  }
}

/**
 * Return the current axis zoom ranges of the active plot, or null if no plot.
 * Shape: { xRange: [min, max], yRange: [min, max] }
 */
export function getPlotZoom() {
  const plotDiv = document.getElementById("plddt-plot");
  if (!plotDiv?._fullLayout) return null;
  const xRange = plotDiv._fullLayout.xaxis?.range?.slice();
  const yRange = plotDiv._fullLayout.yaxis?.range?.slice();
  if (!xRange || !yRange) return null;
  return { xRange, yRange };
}

// Re-apply the PAE overlay after a plot re-render (e.g. model switch).
export function reapplyPaePlotOverlay(xResidues, yResidues) {
  _applyPaeOverlay(xResidues, yResidues);
}

// Clear any active Plotly selection (box-select on scatter, or PAE overlay on heatmap).
export function clearPlotSelection() {
  const plotDiv = document.getElementById("plddt-plot");
  if (!plotDiv?._fullLayout) return;
  const traceType = plotDiv._fullData?.[0]?.type;
  if (traceType === 'scatter') {
    Plotly.restyle('plddt-plot', { selectedpoints: [[]] }, [0]);
    Plotly.relayout('plddt-plot', { shapes: _plddtChainShapes });
  } else if (traceType === 'heatmap' && plotDiv.classList.contains('pae-active')) {
    _clearPaeOverlay();
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
    if (!_paeState) return;
    const residues = _paeState.payload.residues || [];
    if (!residues.length) return;
    const r0 = residues[0] - 0.5;
    const r1 = residues[residues.length - 1] + 0.5;

    if (!globalResidues?.length) {
      _paeHoverLines = [];
      if (_paeActiveBaseShapes.length) {
        Plotly.relayout('plddt-plot', { shapes: _paeActiveBaseShapes });
      }
      return;
    }
    _paeHoverLines = globalResidues.flatMap(r => [
      { type: 'line', x0: r, x1: r, y0: r0, y1: r1,
        line: { color: '#f97316', width: 2, dash: 'dot' }, layer: 'above' },
      { type: 'line', x0: r0, x1: r1, y0: r, y1: r,
        line: { color: '#f97316', width: 2, dash: 'dot' }, layer: 'above' },
    ]);
    Plotly.relayout('plddt-plot', { shapes: [..._paeActiveBaseShapes, ..._paeHoverLines] });
  }
}
