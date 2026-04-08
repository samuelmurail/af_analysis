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
function getStructLib() {
  return window.molstar?.lib?.structure ?? {};
}

// Mol* Script + StructureSelection — needed for MolScript-based atom queries.
// The CDN build exposes Script at window.molstar.Script (top-level export).
function getMolScriptLib() {
  const mol = window.molstar;
  if (!mol) return {};
  return {
    Script: mol.Script ?? mol.lib?.['mol-script/script']?.Script,
    StructureSelection:
      mol.lib?.structure?.StructureSelection ??
      mol.StructureSelection,
  };
}

// Select individual atoms by their CIF _atom_site.id using MolScript.
// Returns an element-loci covering exactly those atoms, or null on failure.
function selectAtomsByIds(plugin, atomCifIds) {
  const { Script, StructureSelection } = getMolScriptLib();
  if (!Script?.getStructureSelection || !StructureSelection?.toLociWithSourceUnits) {
    return null;
  }
  const data = plugin.managers.structure.hierarchy.current.structures[0]?.cell?.obj?.data;
  if (!data) return null;
  try {
    const ids = atomCifIds.map(Number);
    const sel = Script.getStructureSelection(
      Q => Q.struct.generator.atomGroups({
        'atom-test': Q.core.set.has([
          Q.set(...ids),
          Q.core.type.num([Q.struct.atomProperty.core.id()]),
        ]),
      }),
      data
    );
    return StructureSelection.toLociWithSourceUnits(sel);
  } catch (e) {
    console.warn('[molstar] selectAtomsByIds failed:', e);
    return null;
  }
}

function mapGlobalResidue(globalResidue) {
  const idx = Number(globalResidue) - 1;
  // Prefer the structure-derived map (populated after loadStructure) —
  // it uses the actual label_asym_id / label_seq_id values from the file,
  // which covers polymer chains AND single-residue ion chains.
  if (state.residueMap && idx >= 0 && idx < state.residueMap.length) {
    const e = state.residueMap[idx];
    const result = { chainId: e.chainId, localResidue: e.seqId };
    if (e.ligandAtomIdx !== undefined) result.ligandAtomIdx = e.ligandAtomIdx;
    if (e.atomName !== undefined) result.atomName = e.atomName;
    return result;
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

// Build a Mol* element-loci for an array of targets.
// Each target is either { chainId, seqId } for residue-level selection,
// or { chainId, seqId, atomName } (where atomName = CIF _atom_site.id) for atom-precise ligand selection.
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

  // Atom-level targets (ligand atoms): use Script-based selection by CIF _atom_site.id.
  const atomTargets   = targets.filter(t => t.atomName !== undefined);
  const residueTargets = targets.filter(t => t.atomName === undefined);

  if (atomTargets.length > 0) {
    const cifIds = atomTargets.map(t => Number(t.atomName));
    const scriptLoci = selectAtomsByIds(plugin, cifIds);
    if (scriptLoci && residueTargets.length === 0) return scriptLoci;
    // If we also have residue targets, fall through to build both and merge;
    // for simplicity, the Script loci already covers atom targets accurately.
    if (scriptLoci) {
      // Build residue loci and merge manually.
      const residueLoci = _buildResidueLociOnly(plugin, residueTargets, StructureElement, StructureProperties, data);
      if (!residueLoci) return scriptLoci;
      // Merge: combine elements arrays.
      const merged = [...scriptLoci.elements, ...residueLoci.elements];
      return { kind: 'element-loci', structure: data, elements: merged };
    }
    // Script unavailable — fall through to manual loop for all targets.
  }

  return _buildResidueLociOnly(plugin, targets, StructureElement, StructureProperties, data);
}

// Internal helper: build element-loci by iterating units.
// Handles both residue-level targets ({ chainId, seqId }) and
// atom-level targets ({ chainId, seqId, atomName }) where atomName is
// the CIF _atom_site.id stored during loadStructure expansion.
function _buildResidueLociOnly(plugin, targets, StructureElement, StructureProperties, data) {
  try {
    const loc = StructureElement.Location.create(data);
    const unitIndexMap = new Map();

    // Separate residue-level from atom-level targets.
    const residueTargets = targets.filter(t => t.atomName === undefined);
    const atomTargets    = targets.filter(t => t.atomName !== undefined);

    // Build per-chain Set of CIF atom ids for atom-level targets.
    // atomName holds String(_SP.atom.id?.(loc) ?? loc.element) from the expansion.
    const atomSetByChain = new Map();
    for (const t of atomTargets) {
      if (!atomSetByChain.has(t.chainId)) atomSetByChain.set(t.chainId, new Set());
      atomSetByChain.get(t.chainId).add(t.atomName);
    }

    for (const unit of data.units) {
      loc.unit = unit;
      for (let i = 0; i < unit.elements.length; i++) {
        loc.element = unit.elements[i];
        const seq   = StructureProperties.residue.label_seq_id(loc);
        const chain = StructureProperties.chain.label_asym_id(loc);
        let matched = false;

        // Residue-level matching.
        for (const { chainId, seqId } of residueTargets) {
          if (seq !== seqId) continue;
          if (chainId && chain !== chainId) continue;
          matched = true;
          break;
        }

        // Atom-level matching: use the same CIF atom id as stored in atomName.
        if (!matched && atomSetByChain.has(chain)) {
          const cifId = String(StructureProperties.atom.id?.(loc) ?? loc.element);
          if (atomSetByChain.get(chain).has(cifId)) matched = true;
        }

        if (matched) {
          if (!unitIndexMap.has(unit)) unitIndexMap.set(unit, new Set());
          unitIndexMap.get(unit).add(i);
        }
      }
    }

    if (unitIndexMap.size === 0) return null;
    const elements = [];
    for (const [unit, indexSet] of unitIndexMap) {
      elements.push({ unit, indices: new Int32Array([...indexSet].sort((a, b) => a - b)) });
    }
    return { kind: 'element-loci', structure: data, elements };
  } catch (e) {
    console.warn("[molstar] _buildResidueLociOnly failed:", e);
    return null;
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


// Subscribe to Mol* hover events in superpose mode.
// `callback(rowIndex)` is called with the dataset row index of the hovered model,
// or `null` when the pointer leaves the structure.
// Only used during superpose; calling this again replaces the previous subscriber.
let _superposeHoverSub = null;
export function subscribeToSuperposeHover(callback) {
  if (_superposeHoverSub) { _superposeHoverSub.unsubscribe(); _superposeHoverSub = null; }
  if (!state.plugin) return;
  const hoverBehavior = state.plugin.behaviors?.interaction?.hover;
  if (!hoverBehavior) return;
  _superposeHoverSub = hoverBehavior.subscribe((hoverState) => {
    if (_programmaticHighlight) return;
    try {
      const loci = hoverState?.current?.loci;
      if (!loci || loci.kind === 'empty-loci') { callback(null); return; }
      if (!_superposeState) { callback(null); return; }
      // Extract the model number from the unit (Mol* uses 1-based MODEL numbers from PDB).
      let unit = null;
      if (loci.kind === 'element-loci' && loci.elements?.length) {
        unit = loci.elements[0].unit;
      } else if (loci.kind === 'bond-loci' && loci.bonds?.length) {
        unit = loci.bonds[0].aUnit;
      }
      if (!unit) { callback(null); return; }
      // modelNum is 1-based (from PDB MODEL record); map to 0-based frame index.
      const modelNum = unit.model?.modelNum ?? unit.model?.header?.modelNum;
      if (modelNum == null) { callback(null); return; }
      const frameIdx = modelNum - 1;
      const rowIdx = _superposeState.rows[frameIdx] ?? null;
      callback(rowIdx);
    } catch (_) { callback(null); }
  });
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
      callback(g ?? []);
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
  // For ligand atoms, include ligandAtomIdx for atom-precise highlighting.
  const targets = [];
  for (const g of globalResidues) {
    const { chainId, localResidue, ligandAtomIdx, atomName } = mapGlobalResidue(g);
    if (chainId && localResidue != null) {
      const t = { chainId, seqId: localResidue };
      if (ligandAtomIdx !== undefined) t.ligandAtomIdx = ligandAtomIdx;
      if (atomName !== undefined) t.atomName = atomName;
      targets.push(t);
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
      const { chainId, localResidue, ligandAtomIdx, atomName } = mapGlobalResidue(g);
      if (chainId && localResidue != null) {
        const t = { chainId, seqId: localResidue };
        if (ligandAtomIdx !== undefined) t.ligandAtomIdx = ligandAtomIdx;
        if (atomName !== undefined) t.atomName = atomName;
        return t;
      }
      return { chainId: null, seqId: Number(g) };
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
  // 'cluster' uses the built-in uniform theme with the current cluster colour.
  const colorParams = colorName === 'cluster'
    ? { color: 'uniform', colorParams: { value: _clusterColor } }
    : { color: colorName };
  if (!state.plugin) return;
  try {
    if (state.reprCell) {
      await state.plugin.state.data.build().delete(state.reprCell).commit();
      state.reprCell = null;
    }
    if (state.polymerCell) {
      state.reprCell = await state.plugin.builders.structure.representation.addRepresentation(
        state.polymerCell, { type: 'cartoon', ...colorParams }
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
        state.ionCell, { type: 'spacefill', ...colorParams }
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
        state.ligandCell, { type: 'ball-and-stick', ...colorParams }
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
    const { chainId, localResidue, atomName } = mapGlobalResidue(g);
    if (chainId && atomName !== undefined) {
      scored.add(`${chainId}:${localResidue}:${atomName}`);
    } else {
      scored.add(chainId ? `${chainId}:${localResidue}` : `*:${Number(g)}`);
    }
  }
  for (const g of yResidues) {
    const { chainId, localResidue, atomName } = mapGlobalResidue(g);
    if (chainId && atomName !== undefined) {
      aligned.add(`${chainId}:${localResidue}:${atomName}`);
    } else {
      aligned.add(chainId ? `${chainId}:${localResidue}` : `*:${Number(g)}`);
    }
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
  // 'cluster' scheme requires the current color to already be set via setClusterColor.
  await _recolorRepr(colorScheme);
}

// ── Cluster-uniform colour theme ─────────────────────────────────────────────
// A solid-colour theme whose colour is set externally by setClusterColor().
// Each row may belong to a different cluster, so main.js calls setClusterColor
// with the correct palette entry before calling loadStructure.

let _clusterColor = 0x636efa;  // default = first CLUSTER_PALETTE entry

export function setClusterColor(hexInt) {
  _clusterColor = (typeof hexInt === 'number') ? hexInt : 0x636efa;
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
// For ligands, atom-precise keys (chain:seq:cifId) are stored in _paeSelection,
// so 'groupInstance' granularity (called once per atom in ball-and-stick) gives
// per-atom coloring within a single-residue ligand (e.g. ATP).
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
      // 'groupInstance' calls color() once per atom for ball-and-stick/spacefill
      // representations (each atom is its own group in those geometries).
      granularity: 'groupInstance',
      color: (location) => {
        if (!_paeSelection) return WHITE;
        try {
          const chain  = StructureProperties.chain.label_asym_id(location);
          const seq    = StructureProperties.residue.label_seq_id(location);
          // CIF _atom_site.id via StructureProperties.atom.id (if available),
          // otherwise fall back to the raw ElementIndex.
          const cifId  = String(StructureProperties.atom.id?.(location) ?? location.element);
          const resKey  = `${chain}:${seq}`;
          const atomKey = `${chain}:${seq}:${cifId}`;
          const wkey    = `*:${seq}`;
          const inScored  = _paeSelection.scored.has(resKey)  || _paeSelection.scored.has(atomKey)  || _paeSelection.scored.has(wkey);
          const inAligned = _paeSelection.aligned.has(resKey) || _paeSelection.aligned.has(atomKey) || _paeSelection.aligned.has(wkey);
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
  _superposeState = null;  // leaving superpose mode
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
  // Resolve color params: 'cluster' uses built-in uniform theme with current cluster colour.
  const _resolveColorParams = () => {
    const s = document.getElementById('color-scheme')?.value || 'chain-id';
    return s === 'cluster'
      ? { color: 'uniform', colorParams: { value: _clusterColor } }
      : { color: s };
  };
  const polymer = await state.plugin.builders.structure.tryCreateComponentStatic(structure, "polymer");
  if (polymer) {
    state.polymerCell = polymer;
    state.reprCell = await state.plugin.builders.structure.representation.addRepresentation(polymer, {
      type: "cartoon",
      ..._resolveColorParams()
    });
    state.plugin.managers.camera.reset();
  }

  const ion = await state.plugin.builders.structure.tryCreateComponentStatic(structure, "ion");
  if (ion) {
    state.ionCell = ion;
    state.ionReprCell = await state.plugin.builders.structure.representation.addRepresentation(ion, {
      type: "spacefill",
      ..._resolveColorParams(),
    });
  }

  const ligand = await state.plugin.builders.structure.tryCreateComponentStatic(structure, "ligand");
  if (ligand) {
    state.ligandCell = ligand;
    state.ligandReprCell = await state.plugin.builders.structure.representation.addRepresentation(ligand, {
      type: "ball-and-stick",
      ..._resolveColorParams(),
    });
  }

  state.residueMap = buildResidueMapFromStructure(state.plugin);

  // For ligand chains the PAE/pLDDT arrays have one entry per heavy atom,
  // but buildResidueMapFromStructure only stores one entry per unique residue.
  // Expand ligand entries so the residueMap index aligns with pLDDT/PAE positions.
  const _cIds   = state.chainIds   || [];
  const _cLens  = state.chainLengths || [];
  const _cTypes = state.chainTypes  || [];
  if (_cTypes.includes('ligand')) {
    const perChain = {};
    for (const e of state.residueMap) {
      if (!perChain[e.chainId]) perChain[e.chainId] = [];
      perChain[e.chainId].push(e);
    }
    const expanded = [];
    // For atom-level ligand highlighting we need per-atom entries with ligandAtomIdx.
    // Build these directly from the Mol* structure so the ordering matches the CIF file
    // (= PAE/pLDDT atom order).
    const _sl = getStructLib();
    const _SE = _sl.StructureElement;
    const _SP = _sl.StructureProperties;
    const _structData = state.plugin.managers.structure.hierarchy.current.structures[0]?.cell?.obj?.data;
    for (let i = 0; i < _cIds.length; i++) {
      const entries = perChain[_cIds[i]] || [];
      if (_cTypes[i] === 'ligand') {
        // Build one entry per heavy atom using Mol* atom order (= PAE order).
        if (_SE && _SP && _structData) {
          const loc = _SE.Location.create(_structData);
          const ligAtoms = [];
          for (const unit of _structData.units) {
            loc.unit = unit;
            for (let k = 0; k < unit.elements.length; k++) {
              loc.element = unit.elements[k];
              if (_SP.chain.label_asym_id(loc) !== _cIds[i]) continue;
              ligAtoms.push({
                chainId: _cIds[i],
                seqId: _SP.residue.label_seq_id(loc),
                ligandAtomIdx: ligAtoms.length,
                // Store CIF _atom_site.id (via StructureProperties.atom.id if available,
                // otherwise fall back to ElementIndex). Used for Script-based selection
                // and for the 'element'-granularity pae-selection color theme.
                atomName: String(_SP.atom.id?.(loc) ?? loc.element),
              });
            }
          }
          if (ligAtoms.length > 0) {
            expanded.push(...ligAtoms);
            continue;
          }
        }
        // Fallback: repeat the single-residue entry without atom-precision
        const atomCount = Number(_cLens[i] || entries.length || 1);
        const template = entries[0] || { chainId: _cIds[i], seqId: 1 };
        for (let j = 0; j < atomCount; j++) expanded.push({ ...template });
      } else {
        expanded.push(...entries);
      }
    }
    state.residueMap = expanded;
  }

  // Pre-build the reverse map (chainId:seqId → array of 1-based global indices).
  // Storing arrays means hovering a ligand in Mol* highlights all its PAE/pLDDT rows.
  const rm = new Map();
  state.residueMap.forEach((e, i) => {
    const key = `${e.chainId}:${e.seqId}`;
    const existing = rm.get(key);
    if (existing === undefined) rm.set(key, [i + 1]);
    else existing.push(i + 1);
  });
  state._reverseResidueMap = rm;

  if (state.selectedResidues.length > 0) {
    highlightResidues(state.selectedResidues);
  }
}

// ── Superpose multiple aligned models ────────────────────────────────────────

const _SUPERPOSE_PALETTE = [
  0x636efa, 0xef553b, 0x00cc96, 0xab63fa, 0xffa15a,
  0x19d3f3, 0xff6692, 0xb6e880, 0xff97ff, 0xfecb52,
];

// Remembers the last superpose call so re-coloring can replay it.
let _superposeState = null;  // { rows, query, frameColors } | null

export function getSuperposeState() { return _superposeState; }

// frameColors: optional array of hex ints (one per frame) used when colorScheme === 'cluster'.
export async function loadSuperpose(rows, query, frameColors = null) {
  if (!rows || !rows.length) return;
  _superposeState = { rows, query, frameColors };

  const payload = await api(
    `/api/superpose?query=${encodeURIComponent(query)}&rows=${rows.join(',')}`
  );
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
  try { await state.plugin.clear(); } catch (_) {}

  const data = await state.plugin.builders.data.rawData({
    data: payload.structure_text,
    label: 'superpose',
  });
  let trajectory;
  try {
    trajectory = await state.plugin.builders.structure.parseTrajectory(data, payload.structure_format);
  } catch (_) {
    trajectory = await state.plugin.builders.structure.parseTrajectory(data, 'mmcif');
  }

  const colorScheme = document.getElementById('color-scheme')?.value || 'chain-id';
  // 'chain-id' and 'cluster' use a per-frame uniform colour so models are distinguishable
  // even when all frames share identical chain IDs (as in MDAnalysis-written PDB files).
  // 'cluster' uses a per-frame solid colour (Mol* built-in 'uniform' theme) so each model
  // gets its cluster palette entry.  All other schemes (chain-id, af-plddt, …) are passed
  // directly to Mol* — each model is loaded as an independent structure so chain-id, etc.
  // work normally within each frame.
  const useUniform = colorScheme === 'cluster';

  const nFrames = payload.n_frames ?? rows.length;
  for (let i = 0; i < nFrames; i++) {
    try {
      const model = await state.plugin.builders.structure.createModel(trajectory, { modelIndex: i });
      const structure = await state.plugin.builders.structure.createStructure(model);
      // Cluster colour for this frame (only used when useUniform is true).
      const frameHex = frameColors?.[i] ?? _SUPERPOSE_PALETTE[i % _SUPERPOSE_PALETTE.length];
      const polymer = await state.plugin.builders.structure.tryCreateComponentStatic(structure, 'polymer');
      if (polymer) {
        const reprParams = useUniform
          ? { type: 'cartoon', color: 'uniform', colorParams: { value: frameHex } }
          : { type: 'cartoon', color: colorScheme };
        await state.plugin.builders.structure.representation.addRepresentation(polymer, reprParams);
      }
      const ligand = await state.plugin.builders.structure.tryCreateComponentStatic(structure, 'ligand');
      if (ligand) {
        const ligParams = useUniform
          ? { type: 'ball-and-stick', color: 'uniform', colorParams: { value: frameHex } }
          : { type: 'ball-and-stick', color: colorScheme };
        await state.plugin.builders.structure.representation.addRepresentation(ligand, ligParams);
      }
    } catch (_) { /* skip frame if model index is out of range */ }
  }
  state.plugin.managers.camera.reset();
  setMolstarNote(`Superposing ${nFrames} model${nFrames > 1 ? 's' : ''}`);
}
