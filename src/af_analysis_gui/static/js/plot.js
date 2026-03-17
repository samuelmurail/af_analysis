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
    x0: x,
    x1: x,
    y0: 0,
    y1: 100,
    line: { color: "#333", width: 1 }
  }));

  const layout = {
    margin: { l: 40, r: 10, t: 10, b: 40 },
    xaxis: { title: "Residue" },
    yaxis: { title: "predicted LDDT", range: [0, 100] },
    clickmode: "event+select",
    dragmode: "pan",
    shapes
  };

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

export function resizePlot() {
  const plotDiv = document.getElementById("plddt-plot");
  if (plotDiv && typeof Plotly !== "undefined" && Plotly.Plots?.resize) {
    Plotly.Plots.resize(plotDiv);
  }
}
