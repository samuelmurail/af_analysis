import { resizePlot } from "./plot.js";

export function initResizableLayout() {
  const layout = document.getElementById("layout");
  if (!layout) return;

  const splitters = Array.from(document.querySelectorAll(".splitter"));
  let dragging = null;
  let startX = 0;
  let startWidths = [];

  const getWidths = () => Array.from(layout.children).map((el) => el.getBoundingClientRect().width);

  splitters.forEach((splitter, idx) => {
    splitter.addEventListener("mousedown", (e) => {
      e.preventDefault();
      dragging = { idx };
      startX = e.clientX;
      startWidths = getWidths();
      document.body.style.cursor = "col-resize";
    });
  });

  window.addEventListener("mousemove", (e) => {
    if (!dragging) return;

    const delta = e.clientX - startX;
    const leftMin = 220;
    const midMin = 260;
    const rightMin = 200;

    let left = startWidths[0];
    let mid = startWidths[2];
    let right = startWidths[4];

    if (dragging.idx === 0) {
      left = Math.max(leftMin, left + delta);
      mid = Math.max(midMin, mid - delta);
    } else {
      mid = Math.max(midMin, mid + delta);
      right = Math.max(rightMin, right - delta);
    }

    layout.style.gridTemplateColumns = `${left}px 8px ${mid}px 8px ${right}px`;
    resizePlot();
  });

  window.addEventListener("mouseup", () => {
    if (!dragging) return;
    dragging = null;
    document.body.style.cursor = "";
    resizePlot();
  });
}
