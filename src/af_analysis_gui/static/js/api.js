export async function api(path, options = {}) {
  const res = await fetch(path, options);
  const data = await res.json();
  if (!res.ok) throw new Error(data.error || "Request failed");
  return data;
}

export function setStatus(msg, isError = false) {
  const el = document.getElementById("status");
  if (!el) return;
  el.textContent = msg;
  el.style.color = isError ? "#a11927" : "#2d3a57";
}

export function setMolstarNote(msg) {
  const el = document.getElementById("molstar-note");
  if (!el) return;
  el.textContent = msg || "";
}
