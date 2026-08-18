"use strict";

const stages = [
  {
    kicker: "Stage 1 · Input specification",
    title: "The target structure stays fixed",
    summary: "RFdiffusion receives the coordinates of insulin-receptor residues A1–A150. Residues A59, A83 and A91 indicate the desired surface region.",
    has: "Target residue identities, backbone coordinates and hotspot labels.",
    lacks: "A binder backbone, binder sequence, affinity or experimental measurement.",
    takeaway: "A hotspot is spatial guidance—not a force field and not a known energetic hot spot.",
    progress: 0
  },
  {
    kicker: "Stage 2 · Initialize the generated chain",
    title: "Binder residue frames begin as noise",
    summary: "The requested 70–100 binder residues are initialized with noisy positions and orientations. They are not physical atoms diffusing through solvent.",
    has: "A fixed target plus randomly initialized frames for the generated chain.",
    lacks: "Protein-like connectivity, a stable fold or amino-acid identities.",
    takeaway: "Generative noise is a mathematical starting distribution, not a high-temperature molecular configuration.",
    progress: 0
  },
  {
    kicker: "Stage 3 · Reverse diffusion",
    title: "Geometry becomes protein-like step by step",
    summary: "At each noise level, the network uses the current noisy binder, the fixed target and hotspot conditioning to estimate a cleaner set of residue frames.",
    has: "A changing geometric state, target context and a learned structural prior.",
    lacks: "A Boltzmann energy trajectory or a guarantee of experimental binding.",
    takeaway: "Move the slider: target coordinates do not change, while binder translations and rotations become coherent.",
    progress: 72
  },
  {
    kicker: "Stage 4 · Inverse folding",
    title: "ProteinMPNN assigns amino-acid identities",
    summary: "The generated backbone is converted into a spatial graph. ProteinMPNN samples residues compatible with each local backbone environment and the designed interface.",
    has: "A fixed complex backbone and the identities of fixed target residues.",
    lacks: "Experimental folding, solubility, specificity or a measured dissociation constant.",
    takeaway: "Several sequences can be sampled for one backbone; placeholder glycines are not the designed sequence.",
    progress: 100
  },
  {
    kicker: "Stage 5 · Structure prediction",
    title: "Ask whether sequence recovers the design",
    summary: "A predictor receives the designed sequence and target, then produces a new complex model. Agreement with the intended geometry is a self-consistency test.",
    has: "Designed sequences, a predicted complex and confidence estimates.",
    lacks: "Direct evidence that purified molecules bind in solution.",
    takeaway: "The prediction should be compared with the design, not treated as an experimental structure.",
    progress: 100
  },
  {
    kicker: "Stage 6 · Filtering and decision",
    title: "Rank designs without confusing scores with truth",
    summary: "We combine prespecified structural checks: interface uncertainty, design recovery, contacts, clashes, buried area and polar satisfaction.",
    has: "A reproducible shortlist and explicit reasons for rejecting designs.",
    lacks: "Binding ground truth until an experiment such as BLI or SPR is performed.",
    takeaway: "Filtering spends experimental effort more intelligently; it does not replace the experiment.",
    progress: 100
  }
];

const codeHelp = {
  runner: { kind: "PROGRAM", title: "Start RFdiffusion inference", body: "This Python entry point resolves the configuration, selects a compatible checkpoint and runs reverse diffusion. It generates structures; it does not train the network.", example: "One command launches the backbone-generation stage." },
  prefix: { kind: "OUTPUT", title: "Choose the output prefix", body: "Generated PDB files and trajectory metadata are named beneath this prefix. A production experiment should place raw outputs in a run-specific directory.", example: "Expect names based on example_outputs/design_ppi." },
  pdb: { kind: "TARGET INPUT", title: "Load the insulin-receptor coordinates", body: "The PDB supplies the atomic structure used as fixed conditioning. Chain IDs and residue numbers must match every contig and hotspot reference exactly.", example: "A59 means chain A, PDB residue 59—not sequence index 59 after arbitrary preprocessing." },
  contig: { kind: "STRUCTURE SPECIFICATION", title: "Keep one chain and generate another", body: "A1–A150 is retained from the target. /0 creates a chain break. RFdiffusion samples a new binder length between 70 and 100 residues.", example: "The output is a non-covalent target–binder complex, not a fusion protein." },
  hotspots: { kind: "SPATIAL CONDITIONING", title: "Direct generation to a target patch", body: "These three target residues tell the model where the generated chain should form its interface. The model is expected to make additional contacts around them.", example: "Hotspot labels guide interface location; they are not restraints with energies in kcal/mol." },
  number: { kind: "SAMPLING", title: "Generate one candidate", body: "This first executable test intentionally stops after one backbone. It verifies the installation, records one complete trajectory and keeps the scientific claim narrow.", example: "One design is a smoke test, not a search and not a benchmark." },
  seed: { kind: "REPRODUCIBILITY", title: "Name this design with index 42", body: "RFdiffusion appends this design index to the output name. With deterministic mode enabled, the implementation also seeds Python, NumPy and PyTorch from this value.", example: "The output is design_ppi_42.pdb and can be regenerated from the frozen setup." },
  deterministic: { kind: "REPRODUCIBILITY", title: "Use deterministic seeded execution", body: "The runner asks PyTorch for deterministic algorithms and sets all relevant random seeds before the design. Exact cross-hardware reproducibility can still depend on CUDA and library behavior.", example: "Deterministic is a computational control, not evidence that the generated design is unique." },
  trajectory: { kind: "OUTPUT", title: "Preserve both trajectory views", body: "RFdiffusion writes the sampled X(t−1) state and its per-step clean pX0 estimate. Both contain 50 coordinate blocks for this run.", example: "Neither file is a molecular-dynamics trajectory or a free-energy path." },
  ca: { kind: "DENOISER SETTING", title: "Remove additional Cα translation noise", body: "This official example sets the scale of extra translational noise in the reverse updates to zero to favor design quality over additional stochastic diversity.", example: "It does not mean that the binder is initialized without noise." },
  frame: { kind: "DENOISER SETTING", title: "Remove additional frame-rotation noise", body: "The analogous setting controls rotational-frame noise during inference updates. Initial conditions and sampling can still differ between designs.", example: "Translation and rotation are separate parts of each residue frame." }
};

const tabs = Array.from(document.querySelectorAll(".stage-tabs button"));
const slider = document.getElementById("diffusionSlider");
const sliderOutput = document.getElementById("sliderOutput");
const playButton = document.getElementById("playButton");
const resetButton = document.getElementById("resetButton");
const pointsGroup = document.getElementById("binderPoints");
const linesGroup = document.getElementById("binderLines");
const overlay = document.getElementById("predictedOverlay");
const contacts = document.getElementById("interfaceContacts");
const timeLabel = document.getElementById("timeLabel");

let currentStage = 0;
let playing = false;
let playTimer = null;

const aminoAcids = "MDELRKAFEEVAKYGLDPQITN".split("");
const noisyPoints = [
  [605,92],[493,76],[689,135],[547,154],[737,197],[577,225],[655,251],[521,285],[708,314],[610,351],[747,382],[553,404],[674,430],[455,393],[474,317],[437,249],[512,194],[449,123],[643,173],[591,292],[699,221],[518,354]
];
const finalPoints = [
  [403,159],[420,147],[442,145],[461,158],[466,178],[453,194],[431,195],
  [414,212],[414,232],[426,249],[448,253],[469,243],[480,261],[476,284],
  [458,300],[438,300],[425,315],[428,339],[444,356],[467,359],[489,345],[504,324]
];

function svgElement(name, attrs = {}) {
  const node = document.createElementNS("http://www.w3.org/2000/svg", name);
  Object.entries(attrs).forEach(([key, value]) => node.setAttribute(key, value));
  return node;
}

const pointNodes = noisyPoints.map((point, index) => {
  const group = svgElement("g");
  const circle = svgElement("circle", { cx: point[0], cy: point[1], r: 8 });
  const text = svgElement("text", { x: point[0], y: point[1] });
  text.textContent = "G";
  group.append(circle, text);
  pointsGroup.appendChild(group);
  return { group, circle, text };
});

const lineNodes = finalPoints.slice(1).map(() => {
  const line = svgElement("line");
  linesGroup.appendChild(line);
  return line;
});

function smoothStep(x) { return x * x * (3 - 2 * x); }

function updateBinder(progress, stage = currentStage) {
  const fraction = smoothStep(progress / 100);
  const visible = stage > 0;
  pointNodes.forEach((node, index) => {
    const start = noisyPoints[index];
    const end = finalPoints[index];
    const wiggle = Math.sin(index * 2.1 + fraction * 13) * (1 - fraction) * 10;
    const x = start[0] + (end[0] - start[0]) * fraction + wiggle;
    const y = start[1] + (end[1] - start[1]) * fraction - wiggle * .45;
    node.circle.setAttribute("cx", x);
    node.circle.setAttribute("cy", y);
    node.text.setAttribute("x", x);
    node.text.setAttribute("y", y);
    node.group.style.opacity = visible ? "1" : "0";
    node.text.textContent = stage >= 3 ? aminoAcids[index] : "G";
    node.text.style.opacity = stage >= 3 ? "1" : "0";
    node.circle.setAttribute("r", stage >= 3 ? "10" : "8");
  });

  lineNodes.forEach((line, index) => {
    const a = pointNodes[index].circle;
    const b = pointNodes[index + 1].circle;
    line.setAttribute("x1", a.getAttribute("cx"));
    line.setAttribute("y1", a.getAttribute("cy"));
    line.setAttribute("x2", b.getAttribute("cx"));
    line.setAttribute("y2", b.getAttribute("cy"));
    line.style.opacity = visible ? String(Math.max(0, (fraction - .2) / .8)) : "0";
  });

  slider.value = String(progress);
  sliderOutput.value = `${Math.round(progress)}%`;
  sliderOutput.textContent = `${Math.round(progress)}%`;
  const remaining = Math.max(0, Math.round(50 * (1 - progress / 100)));
  timeLabel.textContent = progress < 100 ? `t ≈ ${remaining} · reverse diffusion` : "t = 0 · generated backbone";
  overlay.style.opacity = stage >= 4 ? ".9" : "0";
  contacts.style.opacity = stage >= 2 ? "1" : "0";
}

function renderStage(index, useDefaultProgress = true) {
  currentStage = index;
  const stage = stages[index];
  tabs.forEach((tab, tabIndex) => tab.setAttribute("aria-selected", String(tabIndex === index)));
  document.getElementById("stageKicker").textContent = stage.kicker;
  document.getElementById("stageTitle").textContent = stage.title;
  document.getElementById("stageSummary").textContent = stage.summary;
  document.getElementById("modelHas").textContent = stage.has;
  document.getElementById("modelLacks").textContent = stage.lacks;
  document.getElementById("stageTakeaway").textContent = stage.takeaway;
  slider.disabled = index === 0 || index > 3;
  if (useDefaultProgress) updateBinder(stage.progress, index);
  else updateBinder(Number(slider.value), index);
}

function stopPlayback() {
  playing = false;
  if (playTimer) window.clearInterval(playTimer);
  playTimer = null;
  playButton.textContent = "▶ Play all";
}

function playAll() {
  if (playing) { stopPlayback(); return; }
  playing = true;
  playButton.textContent = "❚❚ Pause";
  renderStage(0);
  let phase = 0;
  let progress = 0;
  let hold = 0;
  playTimer = window.setInterval(() => {
    if (!playing) return;
    hold += 1;
    if (phase === 0 && hold >= 5) { phase = 1; hold = 0; renderStage(1); }
    else if (phase === 1 && hold >= 5) { phase = 2; hold = 0; renderStage(2); progress = 0; }
    else if (phase === 2) {
      progress = Math.min(100, progress + 2);
      updateBinder(progress, 2);
      if (progress === 100) { phase = 3; hold = 0; }
    } else if (phase >= 3 && hold >= 8) {
      renderStage(phase);
      hold = 0;
      phase += 1;
      if (phase > 5) {
        window.setTimeout(() => { renderStage(5); stopPlayback(); }, 1100);
      }
    }
  }, 90);
}

tabs.forEach(tab => tab.addEventListener("click", () => { stopPlayback(); renderStage(Number(tab.dataset.stage)); }));
slider.addEventListener("input", () => {
  stopPlayback();
  if (currentStage < 1 || currentStage > 3) renderStage(2);
  updateBinder(Number(slider.value), currentStage);
});
playButton.addEventListener("click", playAll);
resetButton.addEventListener("click", () => { stopPlayback(); renderStage(0); });

document.querySelectorAll(".code-line").forEach(line => {
  line.addEventListener("click", () => {
    document.querySelectorAll(".code-line").forEach(item => item.classList.remove("active"));
    line.classList.add("active");
    const help = codeHelp[line.dataset.code];
    document.getElementById("codeKind").textContent = help.kind;
    document.getElementById("codeTitle").textContent = help.title;
    document.getElementById("codeBody").textContent = help.body;
    document.querySelector("#codeExample span").textContent = help.example;
  });
});

renderStage(0);

// ---------------------------------------------------------------------------
// Actual RFdiffusion trajectory viewer
// ---------------------------------------------------------------------------

const realRun = window.RFDIFFUSION_REAL_RUN;

if (realRun) {
  const realCanvas = document.getElementById("realTrajectoryCanvas");
  const realContext = realCanvas.getContext("2d");
  const realSlider = document.getElementById("realTrajectorySlider");
  const realPlayButton = document.getElementById("realPlayButton");
  const cameraResetButton = document.getElementById("cameraResetButton");
  const realTimestep = document.getElementById("realTimestep");
  const realStepOutput = document.getElementById("realStepOutput");
  const realConfidence = document.getElementById("realConfidence");
  const trajectoryExplanation = document.getElementById("trajectoryExplanation");
  const confidenceChart = document.getElementById("confidenceChart");

  document.getElementById("statBinderLength").textContent = String(realRun.binder.length);
  document.getElementById("statRuntime").textContent = `${realRun.metadata.modelReportedSeconds.toFixed(2)} s`;

  const targetCoords = realRun.target.ca;
  const finalBinderCoords = realRun.binder.finalCa;
  const allFinalCoords = targetCoords.concat(finalBinderCoords);
  const center = [0, 1, 2].map(axis => allFinalCoords.reduce((sum, xyz) => sum + xyz[axis], 0) / allFinalCoords.length);
  const centeredRadius = Math.max(...allFinalCoords.map(xyz => Math.hypot(xyz[0] - center[0], xyz[1] - center[1], xyz[2] - center[2])));
  const camera = { yaw: -0.62, pitch: 0.34, zoom: 1.0 };
  let realFrame = 0;
  let trajectoryKind = "xtMinus1";
  let realPlaying = false;
  let realAnimation = null;
  let dragStart = null;

  function resetCamera() {
    camera.yaw = -0.62;
    camera.pitch = 0.34;
    camera.zoom = 1.0;
    drawRealFrame();
  }

  function resizeRealCanvas() {
    const ratio = Math.min(window.devicePixelRatio || 1, 2);
    const bounds = realCanvas.getBoundingClientRect();
    const width = Math.max(320, Math.round(bounds.width * ratio));
    const height = Math.max(280, Math.round(bounds.height * ratio));
    if (realCanvas.width !== width || realCanvas.height !== height) {
      realCanvas.width = width;
      realCanvas.height = height;
    }
  }

  function projectCoordinate(xyz) {
    const x = xyz[0] - center[0];
    const y = xyz[1] - center[1];
    const z = xyz[2] - center[2];
    const cy = Math.cos(camera.yaw), sy = Math.sin(camera.yaw);
    const cp = Math.cos(camera.pitch), sp = Math.sin(camera.pitch);
    const xYaw = cy * x + sy * z;
    const zYaw = -sy * x + cy * z;
    const yPitch = cp * y - sp * zYaw;
    const zPitch = sp * y + cp * zYaw;
    const scale = Math.min(realCanvas.width, realCanvas.height) * .41 * camera.zoom / centeredRadius;
    return {
      x: realCanvas.width * .5 + xYaw * scale,
      y: realCanvas.height * .5 - yPitch * scale,
      z: zPitch,
      scale
    };
  }

  function drawTrace(coords, stroke, width, alpha) {
    realContext.save();
    realContext.beginPath();
    coords.forEach((xyz, index) => {
      const point = projectCoordinate(xyz);
      if (index === 0) realContext.moveTo(point.x, point.y);
      else realContext.lineTo(point.x, point.y);
    });
    realContext.strokeStyle = stroke;
    realContext.globalAlpha = alpha;
    realContext.lineWidth = width * (window.devicePixelRatio || 1);
    realContext.lineJoin = "round";
    realContext.lineCap = "round";
    realContext.stroke();
    realContext.restore();
  }

  function drawDots(coords, color, radius, stride = 1, alpha = 1) {
    realContext.save();
    realContext.fillStyle = color;
    realContext.globalAlpha = alpha;
    coords.forEach((xyz, index) => {
      if (index % stride) return;
      const point = projectCoordinate(xyz);
      realContext.beginPath();
      realContext.arc(point.x, point.y, radius * (window.devicePixelRatio || 1), 0, Math.PI * 2);
      realContext.fill();
    });
    realContext.restore();
  }

  function drawActualContacts(binderCoords, completion) {
    if (completion < .75) return;
    realContext.save();
    realContext.setLineDash([7, 8].map(value => value * (window.devicePixelRatio || 1)));
    realContext.strokeStyle = "#ffc978";
    realContext.globalAlpha = Math.min(1, (completion - .75) * 4) * .72;
    realContext.lineWidth = 1.25 * (window.devicePixelRatio || 1);
    realRun.contacts.forEach(contact => {
      const a = projectCoordinate(targetCoords[contact.targetIndex]);
      const b = projectCoordinate(binderCoords[contact.binderIndex]);
      realContext.beginPath();
      realContext.moveTo(a.x, a.y);
      realContext.lineTo(b.x, b.y);
      realContext.stroke();
    });
    realContext.restore();
  }

  function drawHotspots() {
    realRun.target.hotspots.forEach(hotspot => {
      const point = projectCoordinate(hotspot.coord);
      const ratio = window.devicePixelRatio || 1;
      realContext.save();
      realContext.shadowColor = "rgba(230,106,61,.8)";
      realContext.shadowBlur = 11 * ratio;
      realContext.fillStyle = "#f0784c";
      realContext.beginPath();
      realContext.arc(point.x, point.y, 5.2 * ratio, 0, Math.PI * 2);
      realContext.fill();
      realContext.shadowBlur = 0;
      realContext.fillStyle = "#ffd5bf";
      realContext.font = `${11 * ratio}px ui-monospace, monospace`;
      realContext.fillText(hotspot.label, point.x + 8 * ratio, point.y - 8 * ratio);
      realContext.restore();
    });
  }

  function drawRealFrame() {
    resizeRealCanvas();
    const ratio = window.devicePixelRatio || 1;
    const coords = realRun.binder.trajectories[trajectoryKind][realFrame];
    const completion = realFrame / 49;
    const gradient = realContext.createRadialGradient(
      realCanvas.width * .5, realCanvas.height * .5, 10,
      realCanvas.width * .5, realCanvas.height * .5, realCanvas.width * .62
    );
    gradient.addColorStop(0, "#233b35");
    gradient.addColorStop(1, "#101a17");
    realContext.fillStyle = gradient;
    realContext.fillRect(0, 0, realCanvas.width, realCanvas.height);

    drawTrace(targetCoords, "#61b8d3", 3.2, .74);
    drawDots(targetCoords, "#92d7e7", 1.45, 3, .8);
    drawActualContacts(coords, completion);
    drawTrace(coords, "#a783ec", 4.1, .9);
    drawDots(coords, "#dfc1ff", 2.0, 2, .95);
    drawHotspots();

    realContext.fillStyle = "rgba(228,236,232,.65)";
    realContext.font = `${12 * ratio}px ui-monospace, monospace`;
    realContext.fillText("drag: rotate  ·  wheel: zoom", 18 * ratio, (realCanvas.height / ratio - 18) * ratio);
  }

  function updateConfidenceCursor() {
    const oldCursor = confidenceChart.querySelector(".chart-cursor");
    if (oldCursor) oldCursor.remove();
    const x = 28 + realFrame * (372 / 49);
    const cursor = svgElement("line", { x1: x, x2: x, y1: 12, y2: 103, class: "chart-cursor" });
    confidenceChart.appendChild(cursor);
  }

  function renderConfidenceChart() {
    confidenceChart.replaceChildren();
    const values = realRun.binder.confidence.map(item => item.binderMean);
    const points = values.map((value, index) => `${28 + index * (372 / 49)},${103 - value * 86}`).join(" ");
    confidenceChart.appendChild(svgElement("line", { x1: 28, x2: 400, y1: 103, y2: 103, class: "chart-axis" }));
    confidenceChart.appendChild(svgElement("line", { x1: 28, x2: 28, y1: 17, y2: 103, class: "chart-axis" }));
    confidenceChart.appendChild(svgElement("polyline", { points, class: "chart-line" }));
    const startText = svgElement("text", { x: 28, y: 122, class: "chart-label" });
    startText.textContent = "t=50";
    const endText = svgElement("text", { x: 400, y: 122, class: "chart-label", "text-anchor": "end" });
    endText.textContent = "t=1";
    const highText = svgElement("text", { x: 23, y: 21, class: "chart-label", "text-anchor": "end" });
    highText.textContent = "1";
    confidenceChart.append(startText, endText, highText);
    updateConfidenceCursor();
  }

  function setRealFrame(index) {
    realFrame = Math.max(0, Math.min(49, index));
    realSlider.value = String(realFrame);
    const diagnostic = realRun.binder.confidence[realFrame];
    realTimestep.textContent = `t = ${diagnostic.timestep}`;
    realStepOutput.value = `${realFrame + 1} / 50`;
    realStepOutput.textContent = `${realFrame + 1} / 50`;
    realConfidence.textContent = diagnostic.binderMean.toFixed(3);
    updateConfidenceCursor();
    drawRealFrame();
  }

  function stopRealPlayback() {
    realPlaying = false;
    if (realAnimation) window.clearInterval(realAnimation);
    realAnimation = null;
    realPlayButton.textContent = "▶ Play real trajectory";
  }

  function playRealTrajectory() {
    if (realPlaying) { stopRealPlayback(); return; }
    if (realFrame === 49) setRealFrame(0);
    realPlaying = true;
    realPlayButton.textContent = "❚❚ Pause";
    realAnimation = window.setInterval(() => {
      if (realFrame >= 49) { stopRealPlayback(); return; }
      setRealFrame(realFrame + 1);
    }, 115);
  }

  realSlider.addEventListener("input", () => { stopRealPlayback(); setRealFrame(Number(realSlider.value)); });
  realPlayButton.addEventListener("click", playRealTrajectory);
  cameraResetButton.addEventListener("click", resetCamera);
  document.querySelectorAll('input[name="trajectory"]').forEach(input => input.addEventListener("change", event => {
    trajectoryKind = event.target.value;
    trajectoryExplanation.innerHTML = trajectoryKind === "xtMinus1"
      ? "<strong>X<sub>t−1</sub></strong> is the sampled state passed into the next reverse step. Early frames are deliberately noisy."
      : "<strong>pX<sub>0</sub></strong> is the network's clean-structure estimate at each noise level—not the state actually propagated to the next step.";
    drawRealFrame();
  }));

  realCanvas.addEventListener("pointerdown", event => {
    dragStart = { x: event.clientX, y: event.clientY, yaw: camera.yaw, pitch: camera.pitch };
    realCanvas.setPointerCapture(event.pointerId);
  });
  realCanvas.addEventListener("pointermove", event => {
    if (!dragStart) return;
    camera.yaw = dragStart.yaw + (event.clientX - dragStart.x) * .009;
    camera.pitch = Math.max(-1.45, Math.min(1.45, dragStart.pitch + (event.clientY - dragStart.y) * .009));
    drawRealFrame();
  });
  realCanvas.addEventListener("pointerup", () => { dragStart = null; });
  realCanvas.addEventListener("pointercancel", () => { dragStart = null; });
  realCanvas.addEventListener("wheel", event => {
    event.preventDefault();
    camera.zoom = Math.max(.45, Math.min(2.8, camera.zoom * Math.exp(-event.deltaY * .001)));
    drawRealFrame();
  }, { passive: false });
  window.addEventListener("resize", drawRealFrame);

  const selectedCandidate = realRun.proteinMpnn.candidates.find(candidate => candidate.sample === realRun.proteinMpnn.selectedSample);
  document.getElementById("selectedRealSequence").textContent = selectedCandidate.sequence;
  const mpnnRows = document.getElementById("mpnnCandidateRows");
  realRun.proteinMpnn.candidates.forEach(candidate => {
    const row = document.createElement("tr");
    if (candidate.sample === realRun.proteinMpnn.selectedSample) row.className = "selected-candidate";
    [
      `#${candidate.sample}`,
      candidate.score.toFixed(4),
      candidate.globalScore.toFixed(4),
      candidate.sample === realRun.proteinMpnn.selectedSample ? "yes" : "—"
    ].forEach(value => {
      const cell = document.createElement("td");
      cell.textContent = value;
      row.appendChild(cell);
    });
    mpnnRows.appendChild(row);
  });

  const atomCanvas = document.getElementById("atomCanvas");
  const atomContext = atomCanvas.getContext("2d");
  const atomCamera = { yaw: -.62, pitch: .34, zoom: 1.28 };
  const residueSlider = document.getElementById("residueSlider");
  const showResidueLabels = document.getElementById("showResidueLabels");
  const residueSequenceMap = document.getElementById("residueSequenceMap");
  const elementColors = { N: "#315bd4", C: "#a9b0b3", O: "#df4141" };
  const aminoAcidNames = {
    A: "alanine", R: "arginine", N: "asparagine", D: "aspartate", C: "cysteine",
    Q: "glutamine", E: "glutamate", G: "glycine", H: "histidine", I: "isoleucine",
    L: "leucine", K: "lysine", M: "methionine", F: "phenylalanine", P: "proline",
    S: "serine", T: "threonine", W: "tryptophan", Y: "tyrosine", V: "valine"
  };
  const residueClasses = {
    A: "hydrophobic", V: "hydrophobic", I: "hydrophobic", L: "hydrophobic", M: "hydrophobic", F: "hydrophobic", W: "hydrophobic", Y: "hydrophobic",
    S: "polar", T: "polar", N: "polar", Q: "polar",
    D: "acidic", E: "acidic",
    K: "basic", R: "basic", H: "basic",
    G: "special", P: "special", C: "special"
  };
  const residueColors = { hydrophobic: "#d0aa70", polar: "#72bfa0", acidic: "#df7468", basic: "#72a7dd", special: "#bd91dc" };
  let selectedResidue = 1;
  let atomDragStart = null;
  let atomHitPoints = [];

  function atomBonds(atoms) {
    const byResidue = new Map();
    atoms.forEach(atom => {
      if (!byResidue.has(atom.residue)) byResidue.set(atom.residue, new Map());
      byResidue.get(atom.residue).set(atom.atom, atom);
    });
    const bonds = [];
    const residues = Array.from(byResidue.keys()).sort((a, b) => a - b);
    residues.forEach((residue, index) => {
      const current = byResidue.get(residue);
      [["N", "CA"], ["CA", "C"], ["C", "O"]].forEach(([first, second]) => {
        if (current.has(first) && current.has(second)) bonds.push([current.get(first), current.get(second)]);
      });
      if (index < residues.length - 1) {
        const next = byResidue.get(residues[index + 1]);
        if (current.has("C") && next.has("N")) bonds.push([current.get("C"), next.get("N")]);
      }
    });
    return bonds;
  }

  const targetAtomBonds = atomBonds(realRun.target.backboneAtoms);
  const binderAtomBonds = atomBonds(realRun.binder.backboneAtoms);

  function resizeAtomCanvas() {
    const ratio = Math.min(window.devicePixelRatio || 1, 2);
    const bounds = atomCanvas.getBoundingClientRect();
    const width = Math.max(320, Math.round(bounds.width * ratio));
    const height = Math.max(280, Math.round(bounds.height * ratio));
    if (atomCanvas.width !== width || atomCanvas.height !== height) {
      atomCanvas.width = width;
      atomCanvas.height = height;
    }
  }

  function projectAtomCoordinate(xyz) {
    const x = xyz[0] - center[0], y = xyz[1] - center[1], z = xyz[2] - center[2];
    const cy = Math.cos(atomCamera.yaw), sy = Math.sin(atomCamera.yaw);
    const cp = Math.cos(atomCamera.pitch), sp = Math.sin(atomCamera.pitch);
    const xYaw = cy * x + sy * z;
    const zYaw = -sy * x + cy * z;
    const yPitch = cp * y - sp * zYaw;
    const zPitch = sp * y + cp * zYaw;
    const scale = Math.min(atomCanvas.width, atomCanvas.height) * .43 * atomCamera.zoom / centeredRadius;
    return { x: atomCanvas.width * .5 + xYaw * scale, y: atomCanvas.height * .5 - yPitch * scale, z: zPitch };
  }

  function drawAtomBonds(bonds, color, alpha, width) {
    const ratio = window.devicePixelRatio || 1;
    atomContext.save();
    atomContext.strokeStyle = color;
    atomContext.globalAlpha = alpha;
    atomContext.lineWidth = width * ratio;
    atomContext.lineCap = "round";
    bonds.forEach(([first, second]) => {
      const a = projectAtomCoordinate(first.coord), b = projectAtomCoordinate(second.coord);
      atomContext.beginPath();
      atomContext.moveTo(a.x, a.y);
      atomContext.lineTo(b.x, b.y);
      atomContext.stroke();
    });
    atomContext.restore();
  }

  function drawAtomSphere(atom, chain, alpha) {
    const ratio = window.devicePixelRatio || 1;
    const point = projectAtomCoordinate(atom.coord);
    const radius = (atom.atom === "CA" ? 3.3 : 2.8) * ratio;
    const selected = chain === "B" && atom.residue === selectedResidue;
    atomContext.save();
    atomContext.globalAlpha = alpha;
    if (selected) {
      atomContext.shadowColor = residueColors[residueClasses[atom.oneLetter]];
      atomContext.shadowBlur = 10 * ratio;
    }
    atomContext.fillStyle = elementColors[atom.element] || "#d8d8d8";
    atomContext.strokeStyle = selected ? "#fff4d8" : "rgba(255,255,255,.35)";
    atomContext.lineWidth = (selected ? 1.5 : .55) * ratio;
    atomContext.beginPath();
    atomContext.arc(point.x, point.y, selected ? radius * 1.32 : radius, 0, Math.PI * 2);
    atomContext.fill();
    atomContext.stroke();
    atomContext.restore();
    if (chain === "B") atomHitPoints.push({ x: point.x, y: point.y, residue: atom.residue, atom: atom.atom });
  }

  function drawResidueLabels() {
    const ratio = window.devicePixelRatio || 1;
    realRun.binder.backboneAtoms.filter(atom => atom.atom === "CA").forEach(atom => {
      if (!showResidueLabels.checked && atom.residue !== selectedResidue) return;
      const point = projectAtomCoordinate(atom.coord);
      const className = residueClasses[atom.oneLetter];
      const isSelected = atom.residue === selectedResidue;
      atomContext.save();
      atomContext.fillStyle = residueColors[className];
      atomContext.globalAlpha = isSelected ? 1 : .82;
      atomContext.beginPath();
      atomContext.roundRect(point.x + 5 * ratio, point.y - 13 * ratio, (isSelected ? 24 : 15) * ratio, (isSelected ? 20 : 14) * ratio, 3 * ratio);
      atomContext.fill();
      atomContext.fillStyle = "#14201c";
      atomContext.font = `800 ${(isSelected ? 12 : 9) * ratio}px ui-monospace, monospace`;
      atomContext.textAlign = "center";
      atomContext.textBaseline = "middle";
      atomContext.fillText(atom.oneLetter, point.x + (isSelected ? 17 : 12.5) * ratio, point.y + (isSelected ? -3 : -6) * ratio);
      atomContext.restore();
    });
  }

  function drawSelectedAtomNames() {
    const ratio = window.devicePixelRatio || 1;
    realRun.binder.backboneAtoms.filter(atom => atom.residue === selectedResidue).forEach(atom => {
      const point = projectAtomCoordinate(atom.coord);
      atomContext.save();
      atomContext.fillStyle = "rgba(244,248,246,.92)";
      atomContext.font = `700 ${9 * ratio}px ui-monospace, monospace`;
      atomContext.fillText(atom.atom === "CA" ? "Cα" : atom.atom, point.x + 5 * ratio, point.y + 11 * ratio);
      atomContext.restore();
    });
  }

  function drawAtomScene() {
    resizeAtomCanvas();
    atomHitPoints = [];
    const background = atomContext.createRadialGradient(atomCanvas.width * .5, atomCanvas.height * .5, 10, atomCanvas.width * .5, atomCanvas.height * .5, atomCanvas.width * .65);
    background.addColorStop(0, "#253a35");
    background.addColorStop(1, "#0e1815");
    atomContext.fillStyle = background;
    atomContext.fillRect(0, 0, atomCanvas.width, atomCanvas.height);
    drawAtomBonds(targetAtomBonds, "#81918c", .18, 1.1);
    drawAtomBonds(binderAtomBonds, "#d5deda", .64, 1.35);
    const atomsByDepth = [
      ...realRun.target.backboneAtoms.map(atom => ({ atom, chain: "A", alpha: .28, z: projectAtomCoordinate(atom.coord).z })),
      ...realRun.binder.backboneAtoms.map(atom => ({ atom, chain: "B", alpha: .98, z: projectAtomCoordinate(atom.coord).z }))
    ].sort((first, second) => first.z - second.z);
    atomsByDepth.forEach(item => drawAtomSphere(item.atom, item.chain, item.alpha));
    drawResidueLabels();
    drawSelectedAtomNames();
  }

  function selectResidue(residue) {
    selectedResidue = Math.max(1, Math.min(realRun.binder.length, residue));
    const identity = realRun.binder.assignedResidues[selectedResidue - 1];
    const className = residueClasses[identity.oneLetter];
    document.getElementById("residueLetter").textContent = identity.oneLetter;
    document.getElementById("residueLetter").style.background = residueColors[className];
    document.getElementById("residueName").textContent = `${identity.residueName} · ${aminoAcidNames[identity.oneLetter]}`;
    document.getElementById("residuePosition").textContent = `chain B · residue ${selectedResidue} of ${realRun.binder.length}`;
    document.getElementById("residueClass").textContent = className === "special" ? "glycine / proline / cysteine" : className;
    residueSlider.value = String(selectedResidue);
    document.getElementById("residueSliderOutput").value = String(selectedResidue);
    document.getElementById("residueSliderOutput").textContent = String(selectedResidue);
    residueSequenceMap.querySelectorAll("button").forEach((button, index) => button.setAttribute("aria-pressed", String(index + 1 === selectedResidue)));
    drawAtomScene();
  }

  realRun.binder.assignedResidues.forEach(identity => {
    const button = document.createElement("button");
    const className = residueClasses[identity.oneLetter];
    button.type = "button";
    button.textContent = identity.oneLetter;
    button.dataset.residue = String(identity.residue);
    button.style.background = residueColors[className];
    button.setAttribute("aria-label", `Residue ${identity.residue}: ${aminoAcidNames[identity.oneLetter]} (${identity.oneLetter})`);
    button.setAttribute("aria-pressed", String(identity.residue === 1));
    button.addEventListener("click", () => selectResidue(identity.residue));
    residueSequenceMap.appendChild(button);
  });

  residueSlider.addEventListener("input", () => selectResidue(Number(residueSlider.value)));
  showResidueLabels.addEventListener("change", drawAtomScene);
  document.getElementById("atomCameraReset").addEventListener("click", () => {
    atomCamera.yaw = -.62; atomCamera.pitch = .34; atomCamera.zoom = 1.28; drawAtomScene();
  });
  atomCanvas.addEventListener("pointerdown", event => {
    atomDragStart = { x: event.clientX, y: event.clientY, yaw: atomCamera.yaw, pitch: atomCamera.pitch, moved: false };
    atomCanvas.setPointerCapture(event.pointerId);
  });
  atomCanvas.addEventListener("pointermove", event => {
    if (!atomDragStart) return;
    const dx = event.clientX - atomDragStart.x, dy = event.clientY - atomDragStart.y;
    if (Math.hypot(dx, dy) > 3) atomDragStart.moved = true;
    if (!atomDragStart.moved) return;
    atomCamera.yaw = atomDragStart.yaw + dx * .009;
    atomCamera.pitch = Math.max(-1.45, Math.min(1.45, atomDragStart.pitch + dy * .009));
    drawAtomScene();
  });
  atomCanvas.addEventListener("pointerup", event => {
    if (atomDragStart && !atomDragStart.moved) {
      const bounds = atomCanvas.getBoundingClientRect();
      const x = (event.clientX - bounds.left) * atomCanvas.width / bounds.width;
      const y = (event.clientY - bounds.top) * atomCanvas.height / bounds.height;
      const nearest = atomHitPoints.map(point => ({ ...point, distance: Math.hypot(point.x - x, point.y - y) })).sort((a, b) => a.distance - b.distance)[0];
      if (nearest && nearest.distance < 15 * atomCanvas.width / bounds.width) selectResidue(nearest.residue);
    }
    atomDragStart = null;
  });
  atomCanvas.addEventListener("pointercancel", () => { atomDragStart = null; });
  atomCanvas.addEventListener("wheel", event => {
    event.preventDefault();
    atomCamera.zoom = Math.max(.55, Math.min(4.5, atomCamera.zoom * Math.exp(-event.deltaY * .001)));
    drawAtomScene();
  }, { passive: false });
  window.addEventListener("resize", drawAtomScene);

  renderConfidenceChart();
  setRealFrame(0);
  selectResidue(1);
}

// ---------------------------------------------------------------------------
// Ten-seed structural-diversity viewer
// ---------------------------------------------------------------------------

const diversityRun = window.RFDIFFUSION_DIVERSITY;

if (diversityRun) {
  const diversityCanvas = document.getElementById("diversityCanvas");
  const diversityContext = diversityCanvas.getContext("2d");
  const showAllDesigns = document.getElementById("showAllDesigns");
  const seedButtonContainer = document.getElementById("diversitySeedButtons");
  const metricSelect = document.getElementById("diversityMetric");
  const matrixSvg = document.getElementById("diversityMatrix");
  const pairHighlight = document.getElementById("pairHighlight");
  const matrixDefinition = document.getElementById("matrixDefinition");
  const diversityRows = document.getElementById("diversityRows");
  const palette = ["#b88cea", "#62c6a7", "#e7a154", "#e47171", "#6baadd", "#d47bc1", "#90bc57", "#d3bd55", "#6ec6ce", "#a6a0e9"];
  const targetCoordsDiversity = diversityRun.target.ca;
  const allDiversityCoords = targetCoordsDiversity.concat(...diversityRun.designs.map(design => design.alignedBinderCa));
  const diversityCenter = [0, 1, 2].map(axis => targetCoordsDiversity.reduce((sum, xyz) => sum + xyz[axis], 0) / targetCoordsDiversity.length);
  const diversityRadius = Math.max(...allDiversityCoords.map(xyz => Math.hypot(xyz[0] - diversityCenter[0], xyz[1] - diversityCenter[1], xyz[2] - diversityCenter[2])));
  const diversityCamera = { yaw: -0.62, pitch: 0.34, zoom: 1.0 };
  let selectedDesignIndex = 0;
  let diversityDragStart = null;

  function resizeDiversityCanvas() {
    const ratio = Math.min(window.devicePixelRatio || 1, 2);
    const bounds = diversityCanvas.getBoundingClientRect();
    const width = Math.max(320, Math.round(bounds.width * ratio));
    const height = Math.max(280, Math.round(bounds.height * ratio));
    if (diversityCanvas.width !== width || diversityCanvas.height !== height) {
      diversityCanvas.width = width;
      diversityCanvas.height = height;
    }
  }

  function projectDiversity(xyz) {
    const x = xyz[0] - diversityCenter[0];
    const y = xyz[1] - diversityCenter[1];
    const z = xyz[2] - diversityCenter[2];
    const cy = Math.cos(diversityCamera.yaw), sy = Math.sin(diversityCamera.yaw);
    const cp = Math.cos(diversityCamera.pitch), sp = Math.sin(diversityCamera.pitch);
    const xYaw = cy * x + sy * z;
    const zYaw = -sy * x + cy * z;
    const yPitch = cp * y - sp * zYaw;
    const zPitch = sp * y + cp * zYaw;
    const scale = Math.min(diversityCanvas.width, diversityCanvas.height) * .42 * diversityCamera.zoom / diversityRadius;
    return { x: diversityCanvas.width * .5 + xYaw * scale, y: diversityCanvas.height * .5 - yPitch * scale, z: zPitch };
  }

  function drawDiversityTrace(coordinates, color, width, alpha) {
    const ratio = window.devicePixelRatio || 1;
    diversityContext.save();
    diversityContext.beginPath();
    coordinates.forEach((xyz, index) => {
      const point = projectDiversity(xyz);
      if (index === 0) diversityContext.moveTo(point.x, point.y);
      else diversityContext.lineTo(point.x, point.y);
    });
    diversityContext.strokeStyle = color;
    diversityContext.globalAlpha = alpha;
    diversityContext.lineWidth = width * ratio;
    diversityContext.lineJoin = "round";
    diversityContext.lineCap = "round";
    diversityContext.stroke();
    diversityContext.restore();
  }

  function drawDiversityHotspots() {
    const ratio = window.devicePixelRatio || 1;
    diversityRun.target.hotspots.forEach(residue => {
      const point = projectDiversity(targetCoordsDiversity[residue - 1]);
      diversityContext.save();
      diversityContext.fillStyle = "#f0784c";
      diversityContext.shadowColor = "#f0784c";
      diversityContext.shadowBlur = 8 * ratio;
      diversityContext.beginPath();
      diversityContext.arc(point.x, point.y, 4.2 * ratio, 0, Math.PI * 2);
      diversityContext.fill();
      diversityContext.restore();
    });
  }

  function drawDiversityLabel(design, color) {
    const ratio = window.devicePixelRatio || 1;
    const center = [0, 1, 2].map(axis => design.alignedBinderCa.reduce((sum, xyz) => sum + xyz[axis], 0) / design.alignedBinderCa.length);
    const point = projectDiversity(center);
    diversityContext.save();
    diversityContext.fillStyle = color;
    diversityContext.font = `700 ${11 * ratio}px ui-monospace, monospace`;
    diversityContext.fillText(`seed ${design.seed}`, point.x + 7 * ratio, point.y - 7 * ratio);
    diversityContext.restore();
  }

  function drawDiversityScene() {
    resizeDiversityCanvas();
    const background = diversityContext.createRadialGradient(
      diversityCanvas.width * .5, diversityCanvas.height * .5, 20,
      diversityCanvas.width * .5, diversityCanvas.height * .5, diversityCanvas.width * .65
    );
    background.addColorStop(0, "#233a35");
    background.addColorStop(1, "#0f1916");
    diversityContext.fillStyle = background;
    diversityContext.fillRect(0, 0, diversityCanvas.width, diversityCanvas.height);
    drawDiversityTrace(targetCoordsDiversity, "#63bcd6", 3.2, .78);

    if (showAllDesigns.checked) {
      diversityRun.designs.forEach((design, index) => {
        if (index !== selectedDesignIndex) drawDiversityTrace(design.alignedBinderCa, palette[index], 1.65, .26);
      });
    }
    const selected = diversityRun.designs[selectedDesignIndex];
    drawDiversityTrace(selected.alignedBinderCa, "#d2adff", 4.6, .98);
    drawDiversityHotspots();
    drawDiversityLabel(selected, "#ecdfff");
  }

  function selectDiversityDesign(index) {
    selectedDesignIndex = index;
    const design = diversityRun.designs[index];
    seedButtonContainer.querySelectorAll("button").forEach((button, buttonIndex) => button.setAttribute("aria-pressed", String(index === buttonIndex)));
    diversityRows.querySelectorAll("tr").forEach((row, rowIndex) => row.dataset.selected = String(index === rowIndex));
    document.getElementById("selectedDiversitySeed").textContent = String(design.seed);
    document.getElementById("selectedLength").textContent = `${design.length} residues`;
    document.getElementById("selectedRg").textContent = `${design.radiusOfGyrationAngstrom.toFixed(2)} Å`;
    document.getElementById("selectedContacts").textContent = `${design.contactResiduesWithin10A} target residues`;
    document.getElementById("selectedDiversityConfidence").textContent = design.binderConfidence.toFixed(4);
    document.getElementById("selectedHotspots").textContent = ["A59", "A83", "A91"].map(label => design.hotspotMinCaDistanceAngstrom[label].toFixed(2)).join(" / ") + " Å";
    document.getElementById("selectedClashes").textContent = String(design.backbonePairsBelow2A);
    drawDiversityScene();
  }

  diversityRun.designs.forEach((design, index) => {
    const button = document.createElement("button");
    button.type = "button";
    button.textContent = String(design.seed);
    button.setAttribute("aria-label", `Show seed ${design.seed}`);
    button.setAttribute("aria-pressed", String(index === 0));
    button.addEventListener("click", () => selectDiversityDesign(index));
    seedButtonContainer.appendChild(button);

    const row = document.createElement("tr");
    row.dataset.selected = String(index === 0);
    const rowValues = [
      design.seed,
      design.length,
      design.radiusOfGyrationAngstrom.toFixed(2),
      design.contactResiduesWithin10A,
      design.hotspotMinCaDistanceAngstrom.A59.toFixed(2),
      design.hotspotMinCaDistanceAngstrom.A83.toFixed(2),
      design.hotspotMinCaDistanceAngstrom.A91.toFixed(2),
      design.binderConfidence.toFixed(4),
      design.backbonePairsBelow2A
    ];
    rowValues.forEach(value => {
      const cell = document.createElement("td");
      cell.textContent = String(value);
      row.appendChild(cell);
    });
    row.addEventListener("click", () => selectDiversityDesign(index));
    diversityRows.appendChild(row);
  });

  const matrixMetadata = {
    intrinsicTraceRmsdAngstrom: { label: "intrinsic trace RMSD", unit: "Å", definition: "64-point resampled trace RMSD after independent structural superposition.", reverse: false },
    poseChamferDistanceAngstrom: { label: "target-aligned pose Chamfer", unit: "Å", definition: "Symmetric nearest-neighbor Cα distance after target alignment.", reverse: false },
    centerOfMassDistanceAngstrom: { label: "binder center displacement", unit: "Å", definition: "Distance between binder Cα centers after target alignment.", reverse: false },
    contactJaccardSimilarity: { label: "target-contact Jaccard", unit: "", definition: "Overlap of target-residue sets having a binder Cα within 10 Å; higher means more similar.", reverse: true }
  };

  function blendColor(fraction) {
    const stops = [[214,238,229], [242,203,117], [217,104,84]];
    const scaled = Math.max(0, Math.min(1, fraction)) * 2;
    const lower = Math.min(1, Math.floor(scaled));
    const amount = scaled - lower;
    return `rgb(${stops[lower].map((value, index) => Math.round(value + (stops[lower + 1][index] - value) * amount)).join(",")})`;
  }

  function renderDiversityMatrix() {
    matrixSvg.replaceChildren();
    const key = metricSelect.value;
    const matrix = diversityRun.matrices[key];
    const metadata = matrixMetadata[key];
    matrixDefinition.textContent = metadata.definition;
    const offDiagonal = matrix.flatMap((row, i) => row.filter((_, j) => i !== j));
    const minimum = Math.min(...offDiagonal), maximum = Math.max(...offDiagonal);
    const cellSize = 32, left = 70, top = 48;
    diversityRun.seeds.forEach((seed, index) => {
      const xLabel = svgElement("text", { x: left + index * cellSize + cellSize / 2, y: top - 10, class: "matrix-label", "text-anchor": "middle" });
      xLabel.textContent = seed;
      const yLabel = svgElement("text", { x: left - 12, y: top + index * cellSize + 20, class: "matrix-label", "text-anchor": "end" });
      yLabel.textContent = seed;
      matrixSvg.append(xLabel, yLabel);
    });
    matrix.forEach((row, first) => row.forEach((value, second) => {
      const rawFraction = maximum === minimum ? 0 : (value - minimum) / (maximum - minimum);
      const fraction = metadata.reverse ? 1 - rawFraction : rawFraction;
      const cell = svgElement("rect", {
        x: left + second * cellSize,
        y: top + first * cellSize,
        width: cellSize,
        height: cellSize,
        rx: 3,
        fill: first === second ? "#d7d1dc" : blendColor(fraction),
        class: "matrix-cell",
        tabindex: 0
      });
      const exact = `${metadata.label}: ${value.toFixed(4)}${metadata.unit ? ` ${metadata.unit}` : ""}`;
      const title = svgElement("title");
      title.textContent = `Seeds ${diversityRun.seeds[first]} and ${diversityRun.seeds[second]} — ${exact}`;
      cell.appendChild(title);
      const inspect = () => {
        pairHighlight.textContent = first === second
          ? `Seed ${diversityRun.seeds[first]} compared with itself: ${value.toFixed(4)}${metadata.unit ? ` ${metadata.unit}` : ""}.`
          : `Seeds ${diversityRun.seeds[first]} ↔ ${diversityRun.seeds[second]}: ${exact}.`;
      };
      cell.addEventListener("click", inspect);
      cell.addEventListener("keydown", event => { if (event.key === "Enter" || event.key === " ") inspect(); });
      matrixSvg.appendChild(cell);
    }));
  }

  metricSelect.addEventListener("change", renderDiversityMatrix);
  showAllDesigns.addEventListener("change", drawDiversityScene);
  document.getElementById("diversityCameraReset").addEventListener("click", () => {
    diversityCamera.yaw = -.62;
    diversityCamera.pitch = .34;
    diversityCamera.zoom = 1;
    drawDiversityScene();
  });
  diversityCanvas.addEventListener("pointerdown", event => {
    diversityDragStart = { x: event.clientX, y: event.clientY, yaw: diversityCamera.yaw, pitch: diversityCamera.pitch };
    diversityCanvas.setPointerCapture(event.pointerId);
  });
  diversityCanvas.addEventListener("pointermove", event => {
    if (!diversityDragStart) return;
    diversityCamera.yaw = diversityDragStart.yaw + (event.clientX - diversityDragStart.x) * .009;
    diversityCamera.pitch = Math.max(-1.45, Math.min(1.45, diversityDragStart.pitch + (event.clientY - diversityDragStart.y) * .009));
    drawDiversityScene();
  });
  diversityCanvas.addEventListener("pointerup", () => { diversityDragStart = null; });
  diversityCanvas.addEventListener("pointercancel", () => { diversityDragStart = null; });
  diversityCanvas.addEventListener("wheel", event => {
    event.preventDefault();
    diversityCamera.zoom = Math.max(.45, Math.min(2.8, diversityCamera.zoom * Math.exp(-event.deltaY * .001)));
    drawDiversityScene();
  }, { passive: false });
  window.addEventListener("resize", drawDiversityScene);

  renderDiversityMatrix();
  selectDiversityDesign(0);
}
