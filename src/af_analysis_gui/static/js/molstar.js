import { api, setMolstarNote } from "./api.js";
import { state } from "./table.js";

// PAE selection color state — read by the 'pae-selection' color theme per-atom.
let _paeSelection = null; // { scored: Set<"chainId:seqId">, aligned: Set<"chainId:seqId"> } | null

// Guard flag: true while hoverResidues/unhoverResidues are in progress.
// Prevents the programmatic lociHighlights call from re-entering the
// subscribeToMolstarHover subscription (feedback loop via Plotly.relayout
// or Plotly.Fx.hover re-firing plotly_hover → hoverResidues → lociHighlights).
let _programmaticHighlight = false;

// The CDN viewer bundle exposes structure exports under window.molstar.lib.structure.
// Script.getStructureSelection is NOT exported by the CDN build, so we build loci
// manually by iterating structure units using StructureElement + StructureProperties.
function getStructLib() {
  return window.molstar?.lib?.structure ?? {};
}

function mapGlobalResidue(globalResidue) {
  const idx = Number(globalResidue) - 1;
  // Prefer the structure-derived map (populated after loadStructure) —
  // it uses the actual label_asym_id / label_seq_id values from the file,
  // which covers polymer chains AND single-residue ion chains.
  if (state.residueMap && idx >= 0 && idx < state.residueMap.length) {
    const e = state.residueMap[idx];
    return { chainId: e.chainId, localResidue: e.seqId };
  }
  // Fallback: use API chain_ids/chain_lengths metadata.
  let r = Number(globalResidue);
  for (let i = 0; i < state.chainLengths.length; i++) {
    const len = Number(state.chainLengths[i] || 0);
    if (r <= len) {
      return { chainId: String(state.chainIds[i] || ""), localResidue: r };
    }
    r -= len;
  }
  return { chainId: null, localResidue: null };
}

// Build a Mol* element-loci for an array of { chainId, seqId } pairs.
// Uses StructureElement + StructureProperties from lib.structure (available in CDN viewer bundle).
// StructureElement.Loci is duck-typed: { kind: 'element-loci', structure, elements: [{unit, indices}] }
function buildMultiResidueLoci(plugin, targets) {
  const sl = getStructLib();
  const StructureElement = sl.StructureElement;
  const StructureProperties = sl.StructureProperties;
  const data = plugin.managers.structure.hierarchy.current.structures[0]?.cell?.obj?.data;

  if (!StructureElement || !StructureProperties || !data) {
    console.warn("[molstar] buildMultiResidueLoci: unavailable —",
      { StructureElement: !!StructureElement, StructureProperties: !!StructureProperties, data: !!data,
        libStructureKeys: Object.keys(sl) });
    return null;
  }

  try {
    const loc = StructureElement.Location.create(data);
    // Group matched UnitIndex values by unit reference.
    const unitIndexMap = new Map();

    for (const unit of data.units) {
      loc.unit = unit;
      for (let i = 0; i < unit.elements.length; i++) {
        loc.element = unit.elements[i];
        const seq   = StructureProperties.residue.label_seq_id(loc);
        const chain = StructureProperties.chain.label_asym_id(loc);
        for (const { chainId, seqId } of targets) {
          if (seq !== seqId) continue;
          if (chainId && chain !== chainId) continue;
          if (!unitIndexMap.has(unit)) unitIndexMap.set(unit, new Set());
          unitIndexMap.get(unit).add(i);
          break; // matched — no need to check remaining targets for this atom
        }
      }
    }

    if (unitIndexMap.size === 0) return null;

    const elements = [];
    for (const [unit, indexSet] of unitIndexMap) {
      elements.push({ unit, indices: new Int32Array([...indexSet].sort((a, b) => a - b)) });
    }
    console.debug("[molstar] buildMultiResidueLoci %d targets → %d units", targets.length, elements.length);
    return { kind: "element-loci", structure: data, elements };
  } catch (e) {
    console.warn("[molstar] buildMultiResidueLoci failed:", e);
    return null;
  }
}

// Log sample residues from the loaded structure to help debug seq-id mismatches.
function debugStructureResidues(plugin) {
  try {
    const sl = getStructLib();
    const StructureElement = sl.StructureElement;
    const StructureProperties = sl.StructureProperties;
    const data = plugin.managers.structure.hierarchy.current.structures[0]?.cell?.obj?.data;
    if (!data || !StructureElement || !StructureProperties) {
      console.warn("[molstar] debug: StructureElement=%o StructureProperties=%o. lib.structure keys: %o",
        !!StructureElement, !!StructureProperties, Object.keys(sl));
      return;
    }
    const loc = StructureElement.Location.create(data);
    const seen = new Set();
    const samples = [];
    for (const unit of data.units) {
      loc.unit = unit;
      for (let i = 0; i < unit.elements.length; i++) {
        loc.element = unit.elements[i];
        const chain = StructureProperties.chain.label_asym_id(loc);
        const seq   = StructureProperties.residue.label_seq_id(loc);
        const key   = `${chain}:${seq}`;
        if (!seen.has(key)) {
          seen.add(key);
          samples.push({ chain, seq });
          if (samples.length >= 15) break;
        }
      }
      if (samples.length >= 15) break;
    }
    console.debug("[molstar] sample residues (label_asym_id:label_seq_id):", samples);
  } catch (e) {
    console.warn("[molstar] debugStructureResidues failed:", e);
  }
}

// Build a global-residue index → {chainId, seqId} map from the actual Mol* structure atoms.
// Iterates all units in file order, collecting the first occurrence of each (chain, seqId) pair.
// This is the ground-truth mapping used by mapGlobalResidue, covering both polymer and ion chains.
function buildResidueMapFromStructure(plugin) {
  try {
    const sl = getStructLib();
    const StructureElement = sl.StructureElement;
    const StructureProperties = sl.StructureProperties;
    const data = plugin.managers.structure.hierarchy.current.structures[0]?.cell?.obj?.data;
    if (!StructureElement || !StructureProperties || !data) return [];
    const loc = StructureElement.Location.create(data);
    const seen = new Set();
    const map = [];
    for (const unit of data.units) {
      loc.unit = unit;
      for (let i = 0; i < unit.elements.length; i++) {
        loc.element = unit.elements[i];
        const chain = StructureProperties.chain.label_asym_id(loc);
        const seq   = StructureProperties.residue.label_seq_id(loc);
        const key   = `${chain}:${seq}`;
        if (!seen.has(key)) {
          seen.add(key);
          map.push({ chainId: chain, seqId: seq });
        }
      }
    }
    console.debug(`[molstar] residueMap: ${map.length} residues from ${new Set(map.map(e => e.chainId)).size} chains`);
    return map;
  } catch (e) {
    console.warn('[molstar] buildResidueMapFromStructure failed:', e);
    return [];
  }
}


// Subscribe to Mol* hover events.
// `callback(globalResidues[])` is called on every pointer-move over the 3D viewer.
// Passes [] when the pointer leaves a residue.
// Only one subscriber is kept alive at a time; calling this again replaces the previous one.
let _hoverSubscription = null;
export function subscribeToMolstarHover(callback) {
  if (!state.plugin) return;
  if (_hoverSubscription) { _hoverSubscription.unsubscribe(); _hoverSubscription = null; }
  const hoverBehavior = state.plugin.behaviors?.interaction?.hover;
  if (!hoverBehavior) {
    console.warn('[molstar] subscribeToMolstarHover: plugin.behaviors.interaction.hover not found — hover feedback disabled');
    return;
  }
  console.debug('[molstar] subscribeToMolstarHover: subscribing to hover events');
  // Follows the pp-editor useMolstarSelection pattern:
  //   1. Dismiss immediately on empty-loci (cursor left the structure).
  //   2. For element-loci: use StructureElement.Loci.getFirstLocation() to extract
  //      chain + seqId efficiently; fall back to manual iteration on CDN builds
  //      that may not export getFirstLocation.
  //   3. For bond-loci (ligand sticks): read aUnit/aIndex from the first bond entry
  //      and construct a minimal location object, as pp-editor does.
  //   4. Look up the pre-built reverse map (built once in loadStructure) instead of
  //      rebuilding it on every hover event.
  _hoverSubscription = hoverBehavior.subscribe((hoverState) => {
    // Ignore events fired as a side-effect of our own programmatic highlights.
    if (_programmaticHighlight) return;
    try {
      const sl = getStructLib();
      const StructureElement = sl.StructureElement;
      const StructureProperties = sl.StructureProperties;

      const loci = hoverState?.current?.loci;

      // pp-editor pattern: clear hover immediately on empty loci.
      if (!loci || loci.kind === 'empty-loci') { callback([]); return; }
      if (!StructureElement || !StructureProperties || !state._reverseResidueMap) { callback([]); return; }

      let chainId, seqId;

      if (loci.kind === 'element-loci') {
        // Try pp-editor's preferred path: getFirstLocation (may not exist in CDN build).
        const loc = StructureElement.Loci?.getFirstLocation?.(loci)
          ?? _manualFirstLocation(StructureElement, loci);
        if (!loc) { callback([]); return; }
        chainId = StructureProperties.chain.label_asym_id(loc);
        seqId   = StructureProperties.residue.label_seq_id(loc);

      } else if (loci.kind === 'bond-loci' && loci.bonds?.length) {
        // pp-editor pattern for bond/stick atoms (covers ligand ball-and-stick rep).
        const bondLoc = loci.bonds[0];
        const aUnit   = bondLoc.aUnit;
        const aIndex  = bondLoc.aIndex;
        const aElem   = aUnit?.elements?.[aIndex];
        const loc = { structure: loci.structure, unit: aUnit, element: aElem ?? aIndex };
        chainId = StructureProperties.chain.label_asym_id(loc);
        seqId   = StructureProperties.residue.label_seq_id(loc);

      } else {
        callback([]);
        return;
      }

      const g = state._reverseResidueMap.get(`${chainId}:${seqId}`);
      callback(g != null ? [g] : []);
    } catch (e) {
      // Hover fires very frequently — only warn on the first failure.
      console.warn('[molstar] subscribeToMolstarHover handler failed:', e);
    }
  });
}

// Fallback for CDN builds where StructureElement.Loci.getFirstLocation is not exported.
// Returns a minimal Location from the first atom of the first element group.
function _manualFirstLocation(StructureElement, loci) {
  if (!loci.elements?.length) return null;
  const { unit, indices } = loci.elements[0];
  if (!indices?.length) return null;
  const data = state.plugin.managers.structure.hierarchy.current.structures[0]?.cell?.obj?.data;
  if (!data) return null;
  const loc = StructureElement.Location.create(data);
  loc.unit = unit;
  loc.element = unit.elements[indices[0]];
  return loc;
}

// Highlight a list of global residue indices in Mol*.
// Falls back from chain+local → no-chain+local → no-chain+global for unmapped residues.
export function highlightResidues(globalResidues) {
  if (!state.plugin || !globalResidues?.length) return false;

  // Build candidate target lists: prefer chain+local, fall back to global seq id.
  const targets = [];
  for (const g of globalResidues) {
    const { chainId, localResidue } = mapGlobalResidue(g);
    if (chainId && localResidue != null) {
      targets.push({ chainId, seqId: localResidue });
    } else {
      targets.push({ chainId: null, seqId: Number(g) });
    }
  }

  const loci = buildMultiResidueLoci(state.plugin, targets);
  if (loci) {
    state.plugin.managers.interactivity.lociSelects.selectOnly({ loci });
    state.plugin.managers.interactivity.lociHighlights.highlightOnly({ loci });
    if (document.getElementById('auto-focus')?.checked) {
      state.plugin.managers.camera.focusLoci(loci, { extraRadius: 2 });
    }
    setMolstarNote(
      globalResidues.length === 1
        ? `Highlighted residue ${globalResidues[0]}`
        : `Highlighted ${globalResidues.length} residues`
    );
    return true;
  }

  // Second pass: without chain constraint (in case chain mapping is off)
  const targetsNoChain = targets.map(t => ({ chainId: null, seqId: t.seqId }));
  const loci2 = buildMultiResidueLoci(state.plugin, targetsNoChain);
  if (loci2) {
    state.plugin.managers.interactivity.lociSelects.selectOnly({ loci: loci2 });
    state.plugin.managers.interactivity.lociHighlights.highlightOnly({ loci: loci2 });
    if (document.getElementById('auto-focus')?.checked) {
      state.plugin.managers.camera.focusLoci(loci2, { extraRadius: 2 });
    }
    setMolstarNote(`Highlighted ${globalResidues.length} residue(s) (no-chain fallback)`);
    return true;
  }

  setMolstarNote(`Residue highlight not found for: ${globalResidues.join(', ')}`);
  return false;
}

export function highlightResidue(globalResidue) {
  return highlightResidues([globalResidue]);
}

// Transiently highlight (glow) residues on pointer-hover from an external source
// (e.g. Plotly). Uses lociHighlights only — no selection state change, no camera move.
export function hoverResidues(globalResidues) {
  if (!state.plugin) return;
  _programmaticHighlight = true;
  try {
    if (!globalResidues?.length) {
      try { state.plugin.managers.interactivity.lociHighlights.clearHighlights(); } catch {}
      return;
    }
    const targets = globalResidues.map(g => {
      const { chainId, localResidue } = mapGlobalResidue(g);
      return (chainId && localResidue != null)
        ? { chainId, seqId: localResidue }
        : { chainId: null, seqId: Number(g) };
    });
    const loci = buildMultiResidueLoci(state.plugin, targets);
    if (loci) {
      state.plugin.managers.interactivity.lociHighlights.highlightOnly({ loci });
    }
  } finally {
    // Clear after a short delay to cover both synchronous and rAF-deferred
    // behaviors.interaction.hover emissions triggered by lociHighlights.
    setTimeout(() => { _programmaticHighlight = false; }, 50);
  }
}

export function unhoverResidues() {
  if (!state.plugin) return;
  _programmaticHighlight = true;
  try {
    state.plugin.managers.interactivity.lociHighlights.clearHighlights();
  } catch {}
  setTimeout(() => { _programmaticHighlight = false; }, 50);
}

// Clear all Mol* loci selections and highlights.
export function clearMolstarSelection() {
  if (!state.plugin) return;
  state.plugin.managers.interactivity.lociSelects.deselectAll();
  state.plugin.managers.interactivity.lociHighlights?.clearHighlights?.();
  setMolstarNote('');
}

// ── PAE-selection colour theme helpers ──────────────────────────────────────

// Rebuild the cartoon (polymer), spacefill (ion) and ball-and-stick (ligand)
// representations using the given color theme name.
async function _recolorRepr(colorName) {
  if (!state.plugin) return;
  try {
    if (state.reprCell) {
      await state.plugin.state.data.build().delete(state.reprCell).commit();
      state.reprCell = null;
    }
    if (state.polymerCell) {
      state.reprCell = await state.plugin.builders.structure.representation.addRepresentation(
        state.polymerCell, { type: 'cartoon', color: colorName }
      );
    }
  } catch (e) {
    console.warn('[molstar] _recolorRepr (polymer) failed:', e);
  }
  try {
    if (state.ionReprCell) {
      await state.plugin.state.data.build().delete(state.ionReprCell).commit();
      state.ionReprCell = null;
    }
    if (state.ionCell) {
      state.ionReprCell = await state.plugin.builders.structure.representation.addRepresentation(
        state.ionCell, { type: 'spacefill', color: colorName }
      );
    }
  } catch (e) {
    console.warn('[molstar] _recolorRepr (ion) failed:', e);
  }
  try {
    if (state.ligandReprCell) {
      await state.plugin.state.data.build().delete(state.ligandReprCell).commit();
      state.ligandReprCell = null;
    }
    if (state.ligandCell) {
      state.ligandReprCell = await state.plugin.builders.structure.representation.addRepresentation(
        state.ligandCell, { type: 'ball-and-stick', color: colorName }
      );
    }
  } catch (e) {
    console.warn('[molstar] _recolorRepr (ligand) failed:', e);
  }
}

// Color the structure with the pae-selection theme:
//   scored residues (x-axis) → green, aligned (y-axis) → orange, both → pink, rest → white.
export async function applyPaeColors(xResidues, yResidues) {
  if (!state.plugin) return;
  const scored  = new Set();
  const aligned = new Set();
  for (const g of xResidues) {
    const { chainId, localResidue } = mapGlobalResidue(g);
    scored.add(chainId ? `${chainId}:${localResidue}` : `*:${Number(g)}`);
  }
  for (const g of yResidues) {
    const { chainId, localResidue } = mapGlobalResidue(g);
    aligned.add(chainId ? `${chainId}:${localResidue}` : `*:${Number(g)}`);
  }
  _paeSelection = { scored, aligned };
  await _recolorRepr('pae-selection');
  // Remove standard selection/highlight visuals so they don't clash.
  try {
    state.plugin.managers.interactivity.lociSelects.deselectAll();
    state.plugin.managers.interactivity.lociHighlights.clearHighlights();
  } catch {}
  if (document.getElementById('auto-focus')?.checked) {
    const allGlobals = [...new Set([...xResidues, ...yResidues])];
    if (allGlobals.length) {
      const targets = allGlobals.map(g => {
        const { chainId, localResidue } = mapGlobalResidue(g);
        return chainId ? { chainId, seqId: localResidue } : { chainId: null, seqId: Number(g) };
      });
      const loci = buildMultiResidueLoci(state.plugin, targets);
      if (loci) state.plugin.managers.camera.focusLoci(loci, { extraRadius: 2 });
    }
  }
  setMolstarNote(`Scored: ${xResidues.length} residue(s) | Aligned: ${yResidues.length} residue(s)`);
}

// Restore the user-selected color scheme and clear PAE selection state.
export async function clearPaeColors() {
  if (!_paeSelection) return;
  _paeSelection = null;
  const colorScheme = document.getElementById('color-scheme')?.value || 'chain-id';
  await _recolorRepr(colorScheme);
}

export async function highlightTwoGroups(xResidues, yResidues) {
  await applyPaeColors(xResidues, yResidues);
}

// AlphaFold 3 pLDDT color scheme — reads B-factor (= pLDDT) from each atom.
// Colors match the AF3 confidence legend: dark-blue / cyan / yellow / orange.
function registerAfPlddt(plugin) {
  const StructureProperties = getStructLib().StructureProperties;
  if (!StructureProperties) {
    console.warn('[molstar] Cannot register af-plddt: StructureProperties unavailable');
    return;
  }

  const DARK_BLUE = 0x0053D6;  // > 90  Very high
  const CYAN     = 0x65CBF3;  // 70-90 Confident
  const YELLOW   = 0xFFDB13;  // 50-70 Low
  const ORANGE   = 0xFF7D45;  // < 50  Very low
  const GREY     = 0x909090;  // fallback

  const provider = {
    name: 'af-plddt',
    label: 'AlphaFold pLDDT',
    category: 'Confidence',
    isApplicable: () => true,
    factory: (_ctx, props) => ({
      granularity: 'groupInstance',
      color: (location) => {
        try {
          const b = StructureProperties.atom.B_iso_or_equiv(location);
          if (b > 90) return DARK_BLUE;
          if (b > 70) return CYAN;
          if (b > 50) return YELLOW;
          return ORANGE;
        } catch {
          return GREY;
        }
      },
      props,
      description: 'AlphaFold pLDDT confidence coloring',
    }),
    defaultValues: {},
    getParams: () => ({}),
  };

  try {
    plugin.representation.structure.themes.colorThemeRegistry.add(provider);
    console.debug('[molstar] Registered af-plddt color theme');
  } catch (e) {
    console.warn('[molstar] Failed to register af-plddt color theme:', e);
  }
}

// Custom color theme: everything white except PAE-selected residues.
function registerPaeSelection(plugin) {
  const StructureProperties = getStructLib().StructureProperties;
  if (!StructureProperties) {
    console.warn('[molstar] Cannot register pae-selection: StructureProperties unavailable');
    return;
  }
  const WHITE  = 0xFFFFFF;
  const GREEN  = 0x22c55e;  // scored  (x-axis)
  const ORANGE = 0xf97316;  // aligned (y-axis)
  const PINK   = 0xec4899;  // both
  const provider = {
    name: 'pae-selection',
    label: 'PAE Selection',
    category: 'Misc',
    isApplicable: () => true,
    factory: (_ctx, props) => ({
      granularity: 'groupInstance',
      color: (location) => {
        if (!_paeSelection) return WHITE;
        try {
          const chain = StructureProperties.chain.label_asym_id(location);
          const seq   = StructureProperties.residue.label_seq_id(location);
          const key   = `${chain}:${seq}`;
          const wkey  = `*:${seq}`;
          const inScored  = _paeSelection.scored.has(key)  || _paeSelection.scored.has(wkey);
          const inAligned = _paeSelection.aligned.has(key) || _paeSelection.aligned.has(wkey);
          if (inScored && inAligned) return PINK;
          if (inScored)  return GREEN;
          if (inAligned) return ORANGE;
        } catch {}
        return WHITE;
      },
      props,
      description: 'PAE residue selection coloring',
    }),
    defaultValues: {},
    getParams: () => ({}),
  };
  try {
    plugin.representation.structure.themes.colorThemeRegistry.add(provider);
    console.debug('[molstar] Registered pae-selection color theme');
  } catch (e) {
    console.warn('[molstar] Failed to register pae-selection theme:', e);
  }
}

export async function ensureViewer() {
  if (state.viewer) return;
  state.viewer = await molstar.Viewer.create("molstar-viewer", {
    layoutIsExpanded: false,
    layoutShowControls: false,
    layoutShowSequence: false,
    layoutShowLog: false,
    layoutShowLeftPanel: false,
    viewportShowExpand: false,
    viewportShowSelectionMode: false,
    viewportShowAnimation: false,
    pdbProvider: "pdbe",
    emdbProvider: "pdbe"
  });
  state.plugin = state.viewer.plugin;
  registerAfPlddt(state.plugin);
  registerPaeSelection(state.plugin);
}

export async function loadStructure(index) {
  const payload = await api(`/api/structure?index=${index}`);
  await ensureViewer();

  _paeSelection = null;
  state.residueMap = null;
  state._reverseResidueMap = null;
  state.polymerCell = null;
  state.reprCell = null;
  state.ionCell = null;
  state.ionReprCell = null;
  state.ligandCell = null;
  state.ligandReprCell = null;
  try {
    await state.plugin.clear();
  } catch (_err) {
    // Ignore clear failures.
  }

  const data = await state.plugin.builders.data.rawData({ data: payload.structure_text, label: "model" });
  let trajectory = null;
  try {
    trajectory = await state.plugin.builders.structure.parseTrajectory(data, payload.structure_format);
  } catch (_err) {
    trajectory = await state.plugin.builders.structure.parseTrajectory(
      data,
      payload.structure_format === "pdb" ? "mmcif" : "pdb"
    );
  }

  const model = await state.plugin.builders.structure.createModel(trajectory);
  const structure = await state.plugin.builders.structure.createStructure(model);
  const polymer = await state.plugin.builders.structure.tryCreateComponentStatic(structure, "polymer");
  if (polymer) {
    const colorScheme = document.getElementById('color-scheme')?.value || 'chain-id';
    state.polymerCell = polymer;
    state.reprCell = await state.plugin.builders.structure.representation.addRepresentation(polymer, {
      type: "cartoon",
      color: colorScheme
    });
    state.plugin.managers.camera.reset();
  }

  const ion = await state.plugin.builders.structure.tryCreateComponentStatic(structure, "ion");
  if (ion) {
    const colorScheme = document.getElementById('color-scheme')?.value || 'chain-id';
    state.ionCell = ion;
    state.ionReprCell = await state.plugin.builders.structure.representation.addRepresentation(ion, {
      type: "spacefill",
      color: colorScheme,
    });
  }

  const ligand = await state.plugin.builders.structure.tryCreateComponentStatic(structure, "ligand");
  if (ligand) {
    const colorScheme = document.getElementById('color-scheme')?.value || 'chain-id';
    state.ligandCell = ligand;
    state.ligandReprCell = await state.plugin.builders.structure.representation.addRepresentation(ligand, {
      type: "ball-and-stick",
      color: colorScheme,
    });
  }

  debugStructureResidues(state.plugin);
  state.residueMap = buildResidueMapFromStructure(state.plugin);

  // Pre-build the reverse map (chainId:seqId → 1-based global index) once here
  // so it is ready at hover time without being rebuilt on every cursor move.
  const rm = new Map();
  state.residueMap.forEach((e, i) => rm.set(`${e.chainId}:${e.seqId}`, i + 1));
  state._reverseResidueMap = rm;

  if (state.selectedResidues.length > 0) {
    highlightResidues(state.selectedResidues);
  }
}
