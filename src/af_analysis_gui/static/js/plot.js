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

  // Solid white lines at chain boundaries (same style as pLDDT dark lines)
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

  const selEl = document.getElementById("selected-residues");
  if (selEl) selEl.style.display = "";

  Plotly.newPlot("plddt-plot", [trace], layout, { responsive: true });

  const plotDiv = document.getElementById("plddt-plot");

  plotDiv.on("plotly_click", (ev) => {
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
    if (handlers?.onPaeSelect) handlers.onPaeSelect({ xResidues, yResidues });
  });
}

export function resizePlot() {
  const plotDiv = document.getElementById("plddt-plot");
  if (plotDiv && typeof Plotly !== "undefined" && Plotly.Plots?.resize) {
    Plotly.Plots.resize(plotDiv);
  }
}
