import { resizePlot } from "./plot.js";

function _resizeAll() {
  resizePlot();
  const cluster = document.getElementById('cluster-main-plot');
  if (cluster && cluster._fullLayout && typeof Plotly !== 'undefined') {
    Plotly.Plots.resize(cluster);
  }
}

function _initHorizontalSplitter() {
  const splitter = document.getElementById('splitter-h');
  const tableSection = document.getElementById('table-section');
  const workspace = document.getElementById('layout');
  if (!splitter || !tableSection) return;

  let startPos, startSize;

  splitter.addEventListener('mousedown', (e) => {
    e.preventDefault();
    const isLeft = workspace && workspace.classList.contains('workspace--table-left');
    if (isLeft) {
      startPos = e.clientX;
      startSize = tableSection.getBoundingClientRect().width;
      document.body.style.cursor = 'col-resize';
    } else {
      startPos = e.clientY;
      startSize = tableSection.getBoundingClientRect().height;
      document.body.style.cursor = 'row-resize';
    }
    document.body.style.userSelect = 'none';
    document.addEventListener('mousemove', onMove);
    document.addEventListener('mouseup', onUp);
  });

  function onMove(e) {
    const isLeft = workspace && workspace.classList.contains('workspace--table-left');
    if (isLeft) {
      const delta = e.clientX - startPos;
      const newW = Math.max(200, Math.min(window.innerWidth - 400, startSize + delta));
      tableSection.style.width = newW + 'px';
    } else {
      const delta = e.clientY - startPos;
      const newH = Math.max(120, Math.min(window.innerHeight - 280, startSize + delta));
      tableSection.style.height = newH + 'px';
    }
    _resizeAll();
  }

  function onUp() {
    document.body.style.cursor = '';
    document.body.style.userSelect = '';
    document.removeEventListener('mousemove', onMove);
    document.removeEventListener('mouseup', onUp);
    _resizeAll();
  }
}

function _initVerticalSplitter() {
  const splitter = document.getElementById('splitter-v');
  const molstarPanel = document.getElementById('molstar-panel');
  const plotPanel = document.getElementById('plot-panel');
  if (!splitter || !molstarPanel || !plotPanel) return;

  let startX, startMolW, startPlotW;

  splitter.addEventListener('mousedown', (e) => {
    e.preventDefault();
    startX = e.clientX;
    startMolW = molstarPanel.getBoundingClientRect().width;
    startPlotW = plotPanel.getBoundingClientRect().width;
    document.body.style.cursor = 'col-resize';
    document.body.style.userSelect = 'none';
    document.addEventListener('mousemove', onMove);
    document.addEventListener('mouseup', onUp);
  });

  function onMove(e) {
    const totalW = startMolW + startPlotW;
    const delta = e.clientX - startX;
    const newMolW = Math.max(200, Math.min(totalW - 200, startMolW + delta));
    const newPlotW = totalW - newMolW;
    molstarPanel.style.flex = '0 0 ' + newMolW + 'px';
    plotPanel.style.flex = '0 0 ' + newPlotW + 'px';
    _resizeAll();
  }

  function onUp() {
    document.body.style.cursor = '';
    document.body.style.userSelect = '';
    document.removeEventListener('mousemove', onMove);
    document.removeEventListener('mouseup', onUp);
    _resizeAll();
  }
}

export function initResizableLayout() {
  _initHorizontalSplitter();
  _initVerticalSplitter();
}
