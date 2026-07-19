"use strict";
const $ = (s, r = document) => r.querySelector(s);
const el = (t, c, h) => { const e = document.createElement(t); if (c) e.className = c;
  if (h !== undefined) e.innerHTML = h; return e; };
const pct = x => (x == null ? "—" : Math.round(x * 100) + "%");
const CALL = { "likely to fail": "fail", "likely to work": "work", "no-call": "nocall" };
const EVLABEL = { known_resistance_determinant: ["known determinant", "known"],
  statistical_association: ["statistical assoc.", "stat"],
  no_known_resistance_signal: ["no known signal", ""] };

let STATE = { genomes: [], selected: null, meta: {} };

/* ---------- boot ---------- */
async function boot() {
  initTheme(); initTabs();
  STATE.meta = await fetch("/api/meta").then(r => r.json()).catch(() => ({}));
  $("#version").textContent = "v" + (STATE.meta.version || "") + " · " + (STATE.meta.species || "");
  if (STATE.meta.features_synthetic) $("#synthBadge").style.display = "";
  $("#aboutSafety").innerHTML = "<b>⚕️ " + (STATE.meta.safety_notice || "") + "</b>";
  const idx = await fetch("/api/genomes").then(r => r.json()).catch(() => ({ genomes: [] }));
  STATE.genomes = idx.genomes || [];
  renderList();
  if (STATE.genomes.length) selectGenome(filteredList()[0]?.genome_id);
  else $("#reportPane").innerHTML =
    `<div class="empty"><div>No precomputed genomes yet.<br>
     Run <code>scripts/precompute_reports.py</code>, or query a genome id via the API.</div></div>`;
  $("#genomeSearch").addEventListener("input", renderList);
  $("#partitionFilter").addEventListener("change", renderList);
}

/* ---------- sidebar ---------- */
function filteredList() {
  const q = $("#genomeSearch").value.trim().toLowerCase();
  const part = $("#partitionFilter").value;
  return STATE.genomes.filter(g =>
    (!part || g.partition === part) &&
    (!q || g.genome_id.toLowerCase().includes(q) || (g.mlst_group || "").toLowerCase().includes(q)));
}
function renderList() {
  const list = $("#genomeList"); list.innerHTML = "";
  const items = filteredList();
  if (!items.length) { list.appendChild(el("li", "muted", "<div style='padding:14px'>No matches</div>")); return; }
  for (const g of items) {
    const s = g.summary || {};
    const li = el("li", "gitem" + (g.genome_id === STATE.selected ? " active" : ""));
    li.innerHTML =
      `<div class="gid">${g.genome_id}</div>
       <div class="meta"><span>${g.mlst_group || "ST?"}</span><span>·</span>
         <span>${g.partition}</span><span>·</span><span>${g.n_proteins} ORFs</span></div>
       <div class="chips">
         <span class="chip f">${s["likely to fail"] || 0} fail</span>
         <span class="chip w">${s["likely to work"] || 0} work</span>
         <span class="chip n">${s["no-call"] || 0} no-call</span></div>`;
    li.onclick = () => selectGenome(g.genome_id);
    list.appendChild(li);
  }
}

/* ---------- report ---------- */
async function selectGenome(gid) {
  if (!gid) return;
  STATE.selected = gid; renderList();
  $("#reportPane").innerHTML = `<div class="empty"><div class="spinner"></div></div>`;
  const rep = await fetch("/api/report/" + encodeURIComponent(gid)).then(r => r.json());
  renderReport(rep);
}

function renderReport(rep) {
  const pane = $("#reportPane"); pane.innerHTML = "";
  const s = rep.summary || {};
  const head = el("div");
  head.innerHTML =
    `<div class="report-head"><h1>${rep.genome_id}</h1>
      <div class="tags">
        <span class="tag">${rep.species}</span>
        <span class="tag">${rep.n_proteins_predicted} ORFs</span>
        <span class="tag">${rep.scope_ok ? "in scope" : "OUT OF SCOPE"}</span>
        <span class="tag">QC ${rep.qc && rep.qc.passed ? "pass" : "fail"}</span>
      </div></div>`;
  pane.appendChild(head);
  if (rep.features_synthetic)
    pane.appendChild(el("div", "safety", `<span>⚠️</span><span><b>Synthetic M1 features.</b>
      These calls demonstrate the pipeline only — swap in real AMRFinderPlus output for real predictions.</span>`));

  const tiles = el("div", "tiles");
  tiles.innerHTML =
    `<div class="tile f"><div class="n">${s["likely to fail"] || 0}</div><div class="l">Likely to fail</div></div>
     <div class="tile w"><div class="n">${s["likely to work"] || 0}</div><div class="l">Likely to work</div></div>
     <div class="tile n"><div class="n">${s["no-call"] || 0}</div><div class="l">No-call</div></div>
     <div class="tile p"><div class="n">${rep.n_proteins_predicted}</div><div class="l">Proteins scanned</div></div>`;
  pane.appendChild(tiles);

  const D = rep.decisions || [];
  const fail = D.filter(d => d.call === "likely to fail");
  const work = D.filter(d => d.call === "likely to work");
  const ncModel = D.filter(d => d.call === "no-call" && (d.tier === "A" || d.tier === "B"));
  const ncOther = D.filter(d => d.call === "no-call" && d.tier === "C");

  addSection(pane, "Likely to fail", fail);
  addSection(pane, "Likely to work", work);
  addSection(pane, "No-call (modelled panel)", ncModel);
  if (ncOther.length) {
    pane.appendChild(el("div", "section-title", "Not covered — insufficient laboratory evidence to model"));
    pane.appendChild(el("div", "muted",
      ncOther.map(d => d.drug_display).join(" · ")));
  }
}

function addSection(pane, title, decisions) {
  if (!decisions.length) return;
  pane.appendChild(el("div", "section-title", `${title} · ${decisions.length}`));
  const grid = el("div", "cards");
  decisions.forEach(d => grid.appendChild(drugCard(d)));
  pane.appendChild(grid);
}

function drugCard(d) {
  const kind = CALL[d.call];
  const c = el("div", "card " + kind);
  const [evtxt, evcls] = EVLABEL[d.evidence_category] || [d.evidence_category, ""];
  const tstatus = d.target_status;
  const tcls = tstatus === "present" ? "target-present"
             : tstatus === "absent" ? "target-absent" : "target-na";
  const support = (d.supporting_determinants || []).length
    ? `<div class="evidence">Determinant: ${d.supporting_determinants.map(g => `<code>${g}</code>`).join(", ")}</div>` : "";
  let targetEv = "";
  if (d.target_evidence && d.target_evidence.length) {
    const cited = d.target_evidence.map(e => e.gene).join(", ");
    const hits = d.target_evidence.map(e =>
      `<div class="hitline">${e.gene}: ${Math.round((e.identity||0)*100)}% id · ${(e.contig||"").split("|").pop()}</div>`).join("");
    targetEv = `<details class="det"><summary>target evidence (${cited})</summary>${hits}</details>`;
  }
  c.innerHTML =
    `<div class="top">
       <div><div class="drug">${d.drug_display}</div><div class="class">${d.drug_class} · tier ${d.tier}</div></div>
       <span class="call ${kind}">${d.call}</span>
     </div>
     <div class="meter ${kind}"><span style="width:${Math.round((d.confidence||0)*100)}%"></span></div>
     <div class="conf-row"><span>confidence ${pct(d.confidence)}</span>
       <span>${d.p_resistant != null ? "p(R)=" + d.p_resistant.toFixed(2) : ""}</span></div>
     <div class="evrow">
       <span class="pill ${evcls}">${evtxt}</span>
       <span class="pill ${tcls}">target: ${tstatus}</span>
     </div>
     ${support}${targetEv}
     <div class="why">${d.rationale || ""}</div>`;
  return c;
}

/* ---------- performance ---------- */
let perfLoaded = false;
async function loadPerformance() {
  if (perfLoaded) return; perfLoaded = true;
  const pane = $("#perfPane");
  const data = await fetch("/api/metrics").then(r => r.json());
  const M = data.metrics || [];
  const synth = data.synthetic_features
    ? `<span class="badge-synth">SYNTHETIC — illustrative only</span>` : "";
  pane.innerHTML = `<h2>Performance on the grouped hidden test ${synth}</h2>
    <p class="muted">Whole MLST lineages held out from training. Metrics per antibiotic below;
      confidence quality shown by Brier score and the reliability curves.</p>
    <div class="grid2">
      <div class="imgcard"><img src="/reports/reliability.png" alt="reliability"/></div>
      <div class="imgcard"><img src="/reports/performance.png" alt="performance"/></div>
    </div>`;
  const cols = [["drug","Drug"],["tier","Tier"],["n_test","N"],["auroc","AUROC"],
    ["pr_auc","PR-AUC"],["balanced_accuracy","Bal.acc"],["recall_resistant","Recall R"],
    ["recall_susceptible","Recall S"],["brier","Brier"],["no_call_rate","No-call"]];
  const t = el("table", "metrics");
  t.innerHTML = "<thead><tr>" + cols.map(c => `<th>${c[1]}</th>`).join("") + "</tr></thead>";
  const tb = el("tbody");
  M.sort((a,b) => (b.auroc||0)-(a.auroc||0)).forEach(m => {
    const tr = el("tr");
    tr.innerHTML = cols.map(([k]) => {
      let v = m[k];
      if (k === "drug") return `<td>${v}</td>`;
      if (k === "tier") return `<td class="tier${v}">${v}</td>`;
      if (v == null) return `<td class="muted">—</td>`;
      if (k === "n_test") return `<td>${v}</td>`;
      let cls = ""; if (["auroc","pr_auc","balanced_accuracy","recall_resistant","recall_susceptible"].includes(k))
        cls = v >= .85 ? "good" : v >= .7 ? "mid" : "bad";
      if (k === "brier") cls = v <= .1 ? "good" : v <= .2 ? "mid" : "bad";
      return `<td class="${cls}">${(+v).toFixed(k==="no_call_rate"?3:2)}</td>`;
    }).join("");
    tb.appendChild(tr);
  });
  t.appendChild(tb); pane.appendChild(t);
}

/* ---------- tabs & theme ---------- */
function initTabs() {
  document.querySelectorAll(".tab").forEach(btn => btn.onclick = () => {
    document.querySelectorAll(".tab").forEach(b => b.classList.remove("active"));
    document.querySelectorAll(".view").forEach(v => v.classList.remove("active"));
    btn.classList.add("active");
    $("#view-" + btn.dataset.view).classList.add("active");
    if (btn.dataset.view === "performance") loadPerformance();
  });
}
function initTheme() {
  const saved = localStorage.getItem("gfw-theme") || "dark";
  document.documentElement.setAttribute("data-theme", saved);
  $("#themeToggle").onclick = () => {
    const cur = document.documentElement.getAttribute("data-theme");
    const next = cur === "dark" ? "light" : "dark";
    document.documentElement.setAttribute("data-theme", next);
    localStorage.setItem("gfw-theme", next);
  };
}
boot();
