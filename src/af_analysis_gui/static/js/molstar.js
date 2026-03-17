import { api, setMolstarNote } from "./api.js";
import { state } from "./table.js";

// The CDN viewer bundle exposes structure exports under window.molstar.lib.structure.
// Script.getStructureSelection is NOT exported by the CDN build, so we build loci
// manually by iterating structure units using StructureElement + StructureProperties.
function getStructLib() {
  return window.molstar?.lib?.structure ?? {};
}

function mapGlobalResidue(globalResidue) {
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

// Highlight two groups with distinct visual styles:
//   xResidues → persistent selection (Mol* selection colour, typically orange)
//   yResidues → hover-style highlight (Mol* highlight colour, typically green)
export function highlightTwoGroups(xResidues, yResidues) {
  if (!state.plugin) return;

  function buildTargets(globals) {
    const targets = [];
    for (const g of globals) {
      const { chainId, localResidue } = mapGlobalResidue(g);
      targets.push(chainId && localResidue != null
        ? { chainId, seqId: localResidue }
        : { chainId: null, seqId: Number(g) });
    }
    return targets;
  }

  const xTargets = buildTargets(xResidues);
  const yTargets = buildTargets(yResidues);

  const xLoci = xTargets.length ? buildMultiResidueLoci(state.plugin, xTargets) : null;
  const yLoci = yTargets.length ? buildMultiResidueLoci(state.plugin, yTargets) : null;

  if (xLoci) {
    state.plugin.managers.interactivity.lociSelects.selectOnly({ loci: xLoci });
    if (document.getElementById('auto-focus')?.checked) {
      state.plugin.managers.camera.focusLoci(xLoci, { extraRadius: 2 });
    }
  } else {
    state.plugin.managers.interactivity.lociSelects.deselectAll();
  }

  if (yLoci) {
    state.plugin.managers.interactivity.lociHighlights.highlightOnly({ loci: yLoci });
  } else {
    state.plugin.managers.interactivity.lociHighlights.clearHighlights();
  }

  setMolstarNote(
    `Scored: ${xResidues.length} residue(s) | Aligned: ${yResidues.length} residue(s)`
  );
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
}

export async function loadStructure(index) {
  const payload = await api(`/api/structure?index=${index}`);
  await ensureViewer();

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
    await state.plugin.builders.structure.representation.addRepresentation(polymer, {
      type: "cartoon",
      color: colorScheme
    });
    state.plugin.managers.camera.reset();
  }

  debugStructureResidues(state.plugin);

  if (state.selectedResidues.length > 0) {
    highlightResidues(state.selectedResidues);
  }
}
