export async function api(path, options = {}) {
  const res = await fetch(path, options);
  const data = await res.json();
  if (!res.ok) {
    if (data.details) console.error(`[api] ${path} failed:\n`, data.details);
    const err = new Error(data.error || "Request failed");
    err.details = data.details || null;
    throw err;
  }
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
