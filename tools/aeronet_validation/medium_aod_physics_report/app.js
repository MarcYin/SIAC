const app = document.getElementById("app");

const state = {
  data: null,
  view: "overview",
  selectedId: null,
  detail: "evidence",
  search: "",
  transition: "all",
  direction: "all",
  outcome: "all",
  sort: "error",
};

const transitionLabels = {
  gain: "Gain vs fresh",
  loss: "Loss vs fresh",
  stable_hit: "Stable hit",
  stable_miss: "Stable miss",
};

function comparisonConfig() {
  return state.data?.comparison || {
    id: "conflict",
    label: "Prior conflict",
    short_label: "conflict",
    reference_id: "fresh",
    reference_label: "Fresh SIAC",
    reference_short_label: "fresh",
  };
}

function validViews() {
  const views = ["overview", "cases", "experiments", "reproduction", "audit"];
  if (state.data?.surface) views.splice(1, 0, "surface");
  return views;
}

const escapeHtml = (value) => String(value ?? "")
  .replaceAll("&", "&amp;")
  .replaceAll("<", "&lt;")
  .replaceAll(">", "&gt;")
  .replaceAll('"', "&quot;")
  .replaceAll("'", "&#039;");

const finite = (value) => value !== null && value !== undefined && Number.isFinite(Number(value));

function fmt(value, digits = 3) {
  return finite(value) ? Number(value).toFixed(digits) : "NA";
}

function signed(value, digits = 3) {
  if (!finite(value)) return "NA";
  const number = Number(value);
  return `${number >= 0 ? "+" : ""}${number.toFixed(digits)}`;
}

function pct(value, digits = 1) {
  return finite(value) ? `${(100 * Number(value)).toFixed(digits)}%` : "NA";
}

function cloudPct(value, digits = 1) {
  if (!finite(value)) return "NA";
  const number = Number(value);
  return `${(number <= 1 ? number * 100 : number).toFixed(digits)}%`;
}

function metric(label, value, sub = "", tone = "") {
  return `<div class="metric"><div class="metric-label">${escapeHtml(label)}</div><div class="metric-value ${escapeHtml(tone)}">${escapeHtml(value)}</div>${sub ? `<div class="metric-sub">${escapeHtml(sub)}</div>` : ""}</div>`;
}

function statusPill(value, label = null) {
  return `<span class="pill ${escapeHtml(value)}"><span class="dot"></span>${escapeHtml(label || transitionLabels[value] || value)}</span>`;
}

function topbar() {
  const links = [
    ["overview", "Overview"],
    ...(state.data.surface ? [["surface", "Surface harmonization"]] : []),
    ["cases", "Case explorer"],
    ["experiments", "Experiments"],
    ["reproduction", "Algorithm"],
    ["audit", "Audit & sources"],
  ];
  return `<header class="topbar">
    <div class="title"><h1>${escapeHtml(state.data.title)}</h1><p>${escapeHtml(state.data.subtitle)}; locked development evidence</p></div>
    <nav class="nav" aria-label="Report views">${links.map(([id, label]) => `<button data-view="${id}" class="${state.view === id ? "active" : ""}">${label}</button>`).join("")}</nav>
  </header>`;
}

function reportShell(body) {
  return `${topbar()}${body}<dialog class="image-dialog" id="image-dialog"><div class="dialog-bar"><span id="dialog-label"></span><button class="icon-button" id="dialog-close" aria-label="Close image" title="Close">×</button></div><img id="dialog-image" alt="Expanded diagnostic"></dialog>`;
}

function variantRows(source = state.data.variants) {
  return [...source].sort((a, b) => {
    if (a.complete !== b.complete) return a.complete ? -1 : 1;
    return (b.metrics.within_ee_rate ?? -1) - (a.metrics.within_ee_rate ?? -1);
  });
}

function experimentTable(rows = variantRows()) {
  return `<div class="table-wrap"><table class="table rank-table"><thead><tr>
    <th>Experiment</th><th>Family</th><th class="right">OK</th><th class="right">Within EE</th><th class="right">Rate</th><th class="right">Bias</th><th class="right">MAE</th><th class="right">RMSE</th><th class="right">Gains</th><th class="right">Losses</th><th>Status</th>
  </tr></thead><tbody>${rows.map((row) => {
    const metricRow = row.metrics;
    const transitions = row.transitions_vs_fresh || {};
    const status = Object.entries(row.status_counts || {}).map(([key, count]) => `${key} ${count}`).join(", ") || "none";
    return `<tr><td class="variant-name"><strong>${escapeHtml(row.label)}</strong><span>${escapeHtml(row.description)}</span><span class="mono wrap">${escapeHtml(row.path)}</span></td><td>${escapeHtml(row.family)}</td><td class="right mono">${metricRow.n}/36</td><td class="right mono">${metricRow.within_ee}/${metricRow.n}</td><td class="right mono">${pct(metricRow.within_ee_rate)}</td><td class="right mono">${signed(metricRow.bias, 4)}</td><td class="right mono">${fmt(metricRow.mae, 4)}</td><td class="right mono">${fmt(metricRow.rmse, 4)}</td><td class="right mono gain">${transitions.gain || 0}</td><td class="right mono loss">${transitions.loss || 0}</td><td><span class="pill ${row.complete ? "met" : "locked"}"><span class="dot"></span>${row.complete ? "Complete" : escapeHtml(status)}</span></td></tr>`;
  }).join("")}</tbody></table></div>`;
}

function overviewView() {
  const data = state.data;
  const target = data.target;
  const comparison = comparisonConfig();
  const fresh = data.variants.find((row) => row.id === comparison.reference_id);
  const fullMedium = data.cohort.archived_metrics.medium;
  const fullAll = data.cohort.archived_metrics.all;
  const targetTone = target.met ? "met" : "not-met";
  return `<main class="page">
    <div class="page-head"><div><h2>Development result</h2><p>Medium AERONET AOD 0.25-0.85; retrieval cloud fraction below 20%; site-group holdout not opened</p></div><div class="page-head-actions"><a class="button-link" href="downloads/cases.csv">Cases CSV</a><a class="button-link" href="downloads/variants.csv">Experiments CSV</a><a class="button-link" href="data/report.json">Report JSON</a></div></div>
    <section class="summary-band">
      <div class="page-head"><div><h2>${target.met ? "Development gate met" : "Development gate not met"}</h2><p>The best complete physical candidate is ${escapeHtml(target.best_variant)}. It changes the fresh replay by the paired transitions shown below; the locked holdout remains unscored.</p></div><div class="status-line ${targetTone}"><span class="dot"></span>${target.met ? "ACCEPTED" : "BELOW 87% GATE"}</div></div>
      <div class="metric-strip">
        ${metric("Best development", `${target.best_hits}/36`, pct(target.best_rate), targetTone)}
        ${metric("Required", `${target.development_required_hits}/36`, "at least 87.0%")}
        ${metric("Gap", `${target.gap_cases} cases`, "to development gate", target.gap_cases > 0 ? "not-met" : "met")}
        ${metric("Fresh replay", `${fresh.metrics.within_ee}/36`, pct(fresh.metrics.within_ee_rate))}
        ${metric("Archived medium", `${fullMedium.within_ee}/${fullMedium.n}`, pct(fullMedium.within_ee_rate))}
        ${metric("Holdout", "LOCKED", `${data.cohort.holdout_count} total cases`, "locked")}
      </div>
      <div class="target-track" aria-label="Best development performance against target"><span style="width:${Math.min(100, 100 * target.best_rate)}%"></span><i style="left:${100 * target.within_ee_rate}%"></i></div>
    </section>

    <section class="section"><div class="section-head"><div><h3>Retrieval and physical-evidence overview</h3><p>All 36 locked medium-development cases. Expected-error regions and signed error scales are identical within each comparison.</p></div></div><div class="figure-grid">
      <figure class="figure"><img class="zoomable" src="assets/retrieval-scatter.png" alt="AOD retrieval comparison"><figcaption>${escapeHtml(data.overview_captions?.retrieval_scatter || "AERONET versus retrieved AOD for the principal complete candidates.")}</figcaption></figure>
      <figure class="figure"><img class="zoomable" src="assets/variant-ranking.png" alt="Experiment performance ranking"><figcaption>Complete experiment performance on the development subset; the vertical line is 87%.</figcaption></figure>
      <figure class="figure"><img class="zoomable" src="assets/evidence-matrix.png" alt="Signed error matrix by case and experiment"><figcaption>${escapeHtml(data.overview_captions?.evidence_matrix || "Signed error in expected-error units across retrieval and physical estimators.")}</figcaption></figure>
      <figure class="figure"><img class="zoomable" src="assets/spatial-error.png" alt="Spatial distribution of AOD changes"><figcaption>${escapeHtml(data.overview_captions?.spatial_error || "Geographic distribution of current-code error and candidate changes.")}</figcaption></figure>
    </div></section>

    <section class="section"><div class="section-head"><div><h3>Experiment ledger</h3><p>Sorted by complete development performance. Gains and losses are paired against the 36-case fresh current-code replay.</p></div><span class="mono muted">full archived low-cloud: ${fullAll.within_ee}/${fullAll.n} (${pct(fullAll.within_ee_rate)})</span></div>${experimentTable()}</section>
  </main>`;
}

function surfaceView() {
  const surface = state.data.surface;
  const terrain = surface.terrain_ablation;
  const terrainRetrieval = surface.terrain_retrieval_ablation;
  const figures = (surface.figures || []).map((item) => `<figure class="figure ${item.wide ? "wide-figure" : ""}"><img class="zoomable" src="${escapeHtml(item.path)}" alt="${escapeHtml(item.title)}"><figcaption><strong>${escapeHtml(item.title)}</strong><br>${escapeHtml(item.caption)}</figcaption></figure>`).join("");
  const metricRows = (surface.band_metrics || []).map((row) => `<tr><td>${escapeHtml(row.training_set)}</td><td>${escapeHtml(row.method)}</td><td>${escapeHtml(row.band)}</td><td class="right mono">${fmt(row.scene_bias, 6)}</td><td class="right mono">${fmt(row.scene_mae, 6)}</td><td class="right mono">${fmt(row.scene_rmse, 6)}</td><td class="right mono">${fmt(row.pixel_mae, 6)}</td><td class="right mono">${fmt(row.pixel_rmse, 6)}</td></tr>`).join("");
  const siteRows = (surface.site_rows || []).map((row) => `<tr><td>${escapeHtml(row.site)}</td><td class="right mono">${row.scene_count}</td><td class="right mono">${row.sample_count}</td><td class="right mono">${fmt(row.maiac_aot_median, 3)}</td><td class="right mono">${fmt(row.delta_aot_median, 3)}</td><td class="right mono">${fmt(row.identity_visible_mae, 5)}</td><td class="right mono">${fmt(row.harmonized_visible_mae, 5)}</td></tr>`).join("");
  const terrainRows = (terrain?.band_rows || []).map((row) => `<tr><td>${escapeHtml(row.band)}</td><td class="right mono">${fmt(row.baseline_scene_bias, 6)}</td><td class="right mono">${fmt(row.terrain_scene_bias, 6)}</td><td class="right mono">${signed(row.scene_bias_delta, 6)}</td><td class="right mono">${fmt(row.baseline_scene_mae, 6)}</td><td class="right mono">${fmt(row.terrain_scene_mae, 6)}</td><td class="right mono ${Number(row.scene_mae_delta) <= 0 ? "gain" : "loss"}">${signed(row.scene_mae_delta, 6)}</td></tr>`).join("");
  const terrainSection = terrain ? `<section class="section"><div class="section-head"><div><h3>Terrain-feature ablation</h3><p>Both models use the same sampled exact pairs and held-site folds. The only added inputs are GLO-30 local elevation, slope, and solar-incidence geometry. Negative MAE delta favours the terrain model.</p></div><div class="link-row"><a href="downloads/terrain-surface-band-metrics.csv">Band CSV</a><a href="downloads/terrain-surface-scenes.csv">Scene CSV</a><a href="downloads/terrain-surface-sites.csv">Site CSV</a></div></div><div class="table-wrap"><table class="table rank-table"><thead><tr><th>Band</th><th class="right">Baseline bias</th><th class="right">Terrain bias</th><th class="right">Bias delta</th><th class="right">Baseline scene MAE</th><th class="right">Terrain scene MAE</th><th class="right">MAE delta</th></tr></thead><tbody>${terrainRows}</tbody></table></div></section>` : "";
  const terrainRaw = terrain ? rawPanel("Terrain-conditioned ablation", { pair_audit: terrain.pair_audit, model: terrain.model, scene_count: terrain.scene_count, site_count: terrain.site_count, retrieval: terrainRetrieval ? { case_count: terrainRetrieval.case_count, transitions: terrainRetrieval.transitions } : null }) : "";
  const terrainRetrievalLink = terrainRetrieval ? '<a class="button-link" href="downloads/terrain-retrieval-cases.csv">Terrain AOD CSV</a>' : "";
  return `<main class="page"><div class="page-head"><div><h2>L2A to current-RT surface harmonization</h2><p>Exact same-day operational L2A and MAIAC-conditioned L1C surface pairs; held-site predictions and no AERONET inputs.</p></div><div class="page-head-actions"><a class="button-link" href="assets/pair-examples/pair-gallery-metadata.json">Pair image manifest</a><a class="button-link" href="downloads/terrain-wvp-diagnostic.csv">Terrain/WVP CSV</a>${terrainRetrievalLink}<a class="button-link" href="downloads/exact-pair-scatter-metrics.csv">Pair scatter CSV</a><a class="button-link" href="downloads/surface-band-metrics.csv">Band metrics CSV</a><a class="button-link" href="downloads/surface-scenes.csv">Scene metrics CSV</a><a class="button-link" href="downloads/surface-sites.csv">Site metrics CSV</a></div></div>
    <div class="metric-strip">
      ${metric("Exact paired pixels", surface.pair_audit.sample_count.toLocaleString(), `${surface.pair_audit.successful_scenes} successful scenes`)}
      ${metric("Clean training pixels", surface.clean_training.sample_count.toLocaleString(), `${surface.clean_training.scene_count} clean of ${surface.application_audit?.scene_count ?? surface.clean_training.scene_count} applied scenes`)}
      ${metric("Clean training sites", surface.clean_training.site_count, `${surface.clean_training.missing_site_count} sites without clean pairs`)}
      ${metric("AOD gate", `<= ${fmt(surface.clean_training.maiac_aot_max, 2)}`, `${surface.domain_guard_audit?.mapping_applied_scenes ?? "-"}/${surface.domain_guard_audit?.scene_count ?? "-"} scenes mapped in bounded control`)}
      ${metric("AERONET features", "NONE", "uses_aeronet=false", "hit")}
      ${metric("Cross-validation", "HELD SITE", `${surface.clean_training.fold_count} fixed folds`)}
    </div>
    <section class="section"><div class="section-head"><div><h3>Surface evidence</h3><p>Scene-level out-of-fold errors are the relevant unit for a prior composite. Pixel metrics remain visible for completeness.</p></div></div><div class="figure-grid">${figures}</div></section>
    <section class="section"><div class="section-head"><div><h3>Per-band validation</h3><p>Identity is unmodified L2A; legacy is the current visible AOD offset; harmonized is capped ridge residual prediction.</p></div></div><div class="table-wrap"><table class="table rank-table"><thead><tr><th>Training domain</th><th>Method</th><th>Band</th><th class="right">Scene bias</th><th class="right">Scene MAE</th><th class="right">Scene RMSE</th><th class="right">Pixel MAE</th><th class="right">Pixel RMSE</th></tr></thead><tbody>${metricRows}</tbody></table></div></section>
    ${terrainSection}
    <section class="section"><div class="section-head"><div><h3>Clean-domain site coverage</h3><p>Every retained site is shown. These values describe surface-reflectance closure only and are not AOD calibration scores.</p></div></div><div class="table-wrap"><table class="table"><thead><tr><th>Site</th><th class="right">Scenes</th><th class="right">Pixels</th><th class="right">MAIAC AOD</th><th class="right">MAIAC - Sen2Cor AOD</th><th class="right">Raw visible MAE</th><th class="right">Harmonized visible MAE</th></tr></thead><tbody>${siteRows}</tbody></table></div></section>
    <section class="section"><div class="section-head"><div><h3>Model and target provenance</h3><p>The exported JSON is directly executable without scikit-learn.</p></div></div><div class="raw-grid">${rawPanel("Selected spatial-pair manifest", surface.pair_examples)}${rawPanel("Terrain/WVP diagnostic", surface.terrain_wvp_diagnostic)}${terrainRaw}${rawPanel("Target RT", surface.target_rt)}${rawPanel("Model", surface.model)}${rawPanel("Pair audit", surface.pair_audit)}${rawPanel("Primary all-scene application", surface.application_audit)}${rawPanel("Bounded domain-guard application", surface.domain_guard_audit)}${rawPanel("Clean-domain exclusions", surface.clean_training.missing_sites)}</div></section>
  </main>`;
}

function filteredCases() {
  const comparison = comparisonConfig();
  const query = state.search.trim().toLowerCase();
  const rows = state.data.cases.filter((item) => {
    if (query && !`${item.site} ${item.matchup_id}`.toLowerCase().includes(query)) return false;
    if (state.transition !== "all" && item.transition !== state.transition) return false;
    if (state.direction !== "all" && item.direction !== state.direction) return false;
    const hit = item.candidate_by_id[comparison.reference_id].within_ee;
    if (state.outcome === "hit" && !hit) return false;
    if (state.outcome === "miss" && hit) return false;
    return true;
  });
  return rows.sort((a, b) => {
    if (state.sort === "truth") return b.truth - a.truth;
    if (state.sort === "cloud") return (b.cloud_fraction ?? -1) - (a.cloud_fraction ?? -1);
    if (state.sort === "support") return (a.valid_fraction ?? 2) - (b.valid_fraction ?? 2);
    if (state.sort === "site") return a.site.localeCompare(b.site);
    if (state.sort === "conflict") return Math.abs(b.candidate_by_id[comparison.id].error_over_ee ?? 0) - Math.abs(a.candidate_by_id[comparison.id].error_over_ee ?? 0);
    return Math.abs(b.fresh_error_over_ee) - Math.abs(a.fresh_error_over_ee);
  });
}

function selectedCase(rows = state.data.cases) {
  return state.data.cases.find((item) => item.matchup_id === state.selectedId) || rows[0] || state.data.cases[0];
}

function caseSidebar(rows) {
  const comparison = comparisonConfig();
  const list = rows.map((item) => `<button class="case-row ${item.matchup_id === state.selectedId ? "active" : ""}" data-case="${escapeHtml(item.matchup_id)}">
    <span class="case-mark ${escapeHtml(item.transition)}"></span>
    <span class="case-name"><strong>${escapeHtml(item.site)}</strong><span>AERONET ${fmt(item.truth)}; ${escapeHtml(comparison.reference_short_label)} ${fmt(item.fresh)}; ${escapeHtml(comparison.short_label)} ${fmt(item.candidate_by_id[comparison.id].value)}</span></span>
    <span class="case-score"><strong>${signed(item.fresh_error_over_ee, 2)}x EE</strong><span>${escapeHtml(transitionLabels[item.transition])}</span></span>
  </button>`).join("");
  return `<aside class="case-sidebar"><div class="filters">
    <label for="search">Matchup or site<input id="search" type="search" value="${escapeHtml(state.search)}" placeholder="Search"></label>
    <div class="filter-grid">
      <label for="transition">${escapeHtml(comparison.short_label)} transition<select id="transition"><option value="all">All transitions</option>${Object.entries(transitionLabels).map(([id, label]) => `<option value="${id}" ${state.transition === id ? "selected" : ""}>${label}</option>`).join("")}</select></label>
      <label for="outcome">${escapeHtml(comparison.reference_short_label)} outcome<select id="outcome"><option value="all">Hit and miss</option><option value="miss" ${state.outcome === "miss" ? "selected" : ""}>Outside EE</option><option value="hit" ${state.outcome === "hit" ? "selected" : ""}>Within EE</option></select></label>
      <label for="direction">${escapeHtml(comparison.reference_short_label)} direction<select id="direction"><option value="all">Under and over</option><option value="under" ${state.direction === "under" ? "selected" : ""}>Underestimate</option><option value="over" ${state.direction === "over" ? "selected" : ""}>Overestimate</option></select></label>
      <label for="sort">Sort<select id="sort"><option value="error" ${state.sort === "error" ? "selected" : ""}>${escapeHtml(comparison.reference_short_label)} |error / EE|</option><option value="conflict" ${state.sort === "conflict" ? "selected" : ""}>${escapeHtml(comparison.short_label)} |error / EE|</option><option value="truth" ${state.sort === "truth" ? "selected" : ""}>AERONET AOD</option><option value="cloud" ${state.sort === "cloud" ? "selected" : ""}>Cloud fraction</option><option value="support" ${state.sort === "support" ? "selected" : ""}>Lowest support</option><option value="site" ${state.sort === "site" ? "selected" : ""}>Site</option></select></label>
    </div><div class="filter-count">${rows.length} of ${state.data.cases.length} development cases</div>
  </div><div class="case-list">${list || '<div class="empty">No cases match the filters.</div>'}</div></aside>`;
}

function candidateComparisonSvg(item) {
  const rows = item.candidate_list.filter((row) => finite(row.value));
  const truth = item.truth;
  const ee = item.ee;
  const maximum = Math.max(1.0, truth + ee, ...rows.map((row) => Number(row.value))) * 1.06;
  const width = 1180;
  const labelWidth = 300;
  const right = 45;
  const rowHeight = 29;
  const top = 40;
  const height = top + rows.length * rowHeight + 38;
  const plotWidth = width - labelWidth - right;
  const x = (value) => labelWidth + Math.max(0, Math.min(maximum, value)) / maximum * plotWidth;
  const grid = [0, .2, .4, .6, .8, 1].map((fraction) => {
    const value = maximum * fraction;
    return `<line x1="${x(value)}" x2="${x(value)}" y1="22" y2="${height - 28}" stroke="#d7dddf"/><text x="${x(value)}" y="${height - 9}" text-anchor="middle" fill="#647075" font-size="10">${value.toFixed(2)}</text>`;
  }).join("");
  const marks = rows.map((row, index) => {
    const y = top + index * rowHeight;
    const colour = row.within_ee ? "#176d68" : "#a34436";
    return `<text x="10" y="${y + 4}" fill="#182125" font-size="10">${escapeHtml(row.label)}</text><line x1="${x(truth)}" x2="${x(row.value)}" y1="${y}" y2="${y}" stroke="#aeb8bc"/><circle cx="${x(row.value)}" cy="${y}" r="4.5" fill="${colour}"/><text x="${Math.min(width - 38, x(row.value) + 8)}" y="${y + 4}" fill="#5d696e" font-size="9">${Number(row.value).toFixed(3)}</text>`;
  }).join("");
  return `<svg class="comparison-svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="AOD estimates for ${escapeHtml(item.site)}"><rect x="${x(Math.max(0, truth - ee))}" y="22" width="${x(truth + ee) - x(Math.max(0, truth - ee))}" height="${height - 50}" fill="#edf1f1"/><line x1="${x(truth)}" x2="${x(truth)}" y1="18" y2="${height - 28}" stroke="#182125" stroke-width="1.5"/><text x="${x(truth) + 5}" y="15" fill="#182125" font-size="10">AERONET ${truth.toFixed(3)}</text>${grid}${marks}</svg>`;
}

function candidateTable(item) {
  return `<div class="table-wrap"><table class="table"><thead><tr><th>Output</th><th>Source ID</th><th class="right">AOD550</th><th class="right">Signed error</th><th class="right">Error / EE</th><th>Within EE</th></tr></thead><tbody>${item.candidate_list.map((row) => `<tr><td>${escapeHtml(row.label)}</td><td class="mono">${escapeHtml(row.source)}</td><td class="right mono">${fmt(row.value, 6)}</td><td class="right mono">${signed(row.error, 6)}</td><td class="right mono">${signed(row.error_over_ee, 3)}</td><td>${row.within_ee === null ? "NA" : row.within_ee ? '<span class="hit">Yes</span>' : '<span class="miss">No</span>'}</td></tr>`).join("")}</tbody></table></div>`;
}

function evidenceFigure(item, kind, title, caption) {
  const source = item[kind];
  if (!source) return "";
  return `<figure class="figure wide-figure"><div class="section-head"><div><h3>${escapeHtml(title)}</h3><p>${escapeHtml(caption)}</p></div><a class="button-link" href="${escapeHtml(source)}" target="_blank" rel="noopener">Open PNG</a></div><img class="zoomable" src="${escapeHtml(source)}" alt="${escapeHtml(title)} for ${escapeHtml(item.site)}" loading="eager"><figcaption>${escapeHtml(caption)}</figcaption></figure>`;
}

function evidenceTab(item) {
  const terrainLabel = item.terrain_evidence_label || "Terrain-conditioned";
  return `<section class="case-section evidence-stack">
    ${evidenceFigure(item, "harmonized_spatial_image", "Harmonized spatial evidence", "Harmonized surface prior, solved AOD, per-band minima, residuals, and station context on the 60 m grid.")}
    ${evidenceFigure(item, "harmonized_diagnostic_image", "Harmonized spectral and cost evidence", "Harmonized total and per-band likelihoods, signed residuals, reflectance spectrum, and AOD maps.")}
    ${evidenceFigure(item, "terrain_spatial_image", `${terrainLabel} spatial evidence`, `${terrainLabel} L2A history, solved AOD, per-band minima, residuals, and station context on the same 60 m grid.`)}
    ${evidenceFigure(item, "terrain_diagnostic_image", `${terrainLabel} spectral and cost evidence`, `${terrainLabel} total and per-band likelihoods, signed residuals, reflectance spectrum, and AOD maps. This is an independent replay, not a visual transformation of the clean-day arm.`)}
    ${evidenceFigure(item, "spatial_image", "Spatial evidence", "TOA, predicted surface, MAIAC prior, solved AOD, no-backstop replay, per-band minima, residuals, and station context on the 60 m grid.")}
    ${evidenceFigure(item, "diagnostic_image", "Spectral and cost evidence", "Scene and station-window likelihood profiles, per-band residuals, reflectance spectra, truth-node closure, and scalar estimates.")}
    ${evidenceFigure(item, "surface_image", "Surface-model evidence", "Current and pooled B02 surface fields, model difference, historical ET20 OOF bias/spread, uncertainty fields, and solver support.")}
    <div class="link-row"><a href="${escapeHtml(item.result_url)}" target="_blank" rel="noopener">Reference result JSON</a><a href="${escapeHtml(item.cost_url)}" target="_blank" rel="noopener">Reference cost cube NPZ</a>${item.harmonized_result_url ? `<a href="${escapeHtml(item.harmonized_result_url)}" target="_blank" rel="noopener">Harmonized result JSON</a>` : ""}${item.harmonized_cost_url ? `<a href="${escapeHtml(item.harmonized_cost_url)}" target="_blank" rel="noopener">Harmonized cost cube NPZ</a>` : ""}${item.terrain_result_url ? `<a href="${escapeHtml(item.terrain_result_url)}" target="_blank" rel="noopener">${escapeHtml(terrainLabel)} result JSON</a>` : ""}${item.terrain_cost_url ? `<a href="${escapeHtml(item.terrain_cost_url)}" target="_blank" rel="noopener">${escapeHtml(terrainLabel)} cost cube NPZ</a>` : ""}</div>
  </section>`;
}

function comparisonTab(item) {
  return `<section class="case-section"><div class="section-head"><div><h3>All scalar AOD estimates</h3><p>One row per tested physical output or diagnostic minimum. Shading is the AERONET expected-error interval.</p></div></div>${candidateComparisonSvg(item)}<div style="height:16px"></div>${candidateTable(item)}</section>`;
}

function flattenObject(value, prefix = "", depth = 0) {
  if (!value || typeof value !== "object" || Array.isArray(value) || depth >= 4) {
    return [[prefix || "value", Array.isArray(value) || typeof value === "object" ? JSON.stringify(value) : value]];
  }
  return Object.entries(value).flatMap(([key, item]) => flattenObject(item, prefix ? `${prefix}.${key}` : key, depth + 1));
}

function rawPanel(title, value) {
  const rows = flattenObject(value).sort((a, b) => a[0].localeCompare(b[0]));
  return `<section class="raw-panel"><h4>${escapeHtml(title)}</h4><table class="kv-table"><tbody>${rows.length ? rows.map(([key, item]) => `<tr><th>${escapeHtml(key)}</th><td>${escapeHtml(item === null ? "null" : item)}</td></tr>`).join("") : '<tr><td class="muted">No saved fields</td></tr>'}</tbody></table></section>`;
}

function rawTab(item) {
  const terrainLabel = item.terrain_evidence_label || "Terrain";
  const scalar = {
    matchup_id: item.matchup_id,
    truth: item.truth,
    ee_low: item.ee_low,
    ee_high: item.ee_high,
    fresh: item.fresh,
    fresh_error: item.fresh_error,
    fresh_error_over_ee: item.fresh_error_over_ee,
    transition: item.transition,
    lon: item.lon,
    lat: item.lat,
    cloud_fraction: item.cloud_fraction,
    scene_cloud_cover: item.scene_cloud_cover,
    aeronet_count: item.aeronet_count,
    aeronet_std: item.aeronet_std,
    angstrom: item.angstrom,
    elevation_m: item.elevation_m,
    surface_min: item.surface_min,
    band_min_spread: item.band_min_spread,
    cost_per_band: item.cost_per_band,
    curvature: item.curvature,
    valid_fraction: item.valid_fraction,
  };
  const optionalPanel = (title, value) => value && Object.keys(value).length ? rawPanel(title, value) : "";
  return `<section class="case-section"><div class="section-head"><div><h3>Saved numerical diagnostics</h3><p>Unabridged scalar fields from the reference and comparison retrievals, atmospheric prior, surface prior, history application, and solver.</p></div></div><div class="raw-grid">${rawPanel("Case and validation", scalar)}${rawPanel("Reference solver", item.solver)}${rawPanel("Comparison solver", item.harmonized_solver)}${optionalPanel(`${terrainLabel} solver`, item.terrain_solver)}${rawPanel("Atmospheric prior", item.atmo_prior)}${rawPanel("Reference surface BOA prior", item.prior_boa)}${rawPanel("Comparison surface BOA prior", item.harmonized_prior_boa)}${optionalPanel(`${terrainLabel} surface BOA prior`, item.terrain_prior_boa)}${rawPanel("Reference surface uncertainty", item.prior_unc)}${rawPanel("Comparison surface uncertainty", item.harmonized_prior_unc)}${optionalPanel(`${terrainLabel} surface uncertainty`, item.terrain_prior_unc)}${rawPanel("Harmonization history", item.harmonization_history)}${optionalPanel(`${terrainLabel} mapping history`, item.terrain_history)}${rawPanel("Anchor iteration", item.anchor_iterate)}${rawPanel("Figure extraction", item.diagnostic)}</div><div class="link-row" style="margin-top:18px"><a href="${escapeHtml(item.result_url)}" target="_blank" rel="noopener">Reference result JSON</a><a href="${escapeHtml(item.cost_url)}" target="_blank" rel="noopener">Reference cost cube NPZ</a>${item.harmonized_result_url ? `<a href="${escapeHtml(item.harmonized_result_url)}" target="_blank" rel="noopener">Comparison result JSON</a>` : ""}${item.harmonized_cost_url ? `<a href="${escapeHtml(item.harmonized_cost_url)}" target="_blank" rel="noopener">Comparison cost cube NPZ</a>` : ""}${item.terrain_result_url ? `<a href="${escapeHtml(item.terrain_result_url)}" target="_blank" rel="noopener">${escapeHtml(terrainLabel)} result JSON</a>` : ""}${item.terrain_cost_url ? `<a href="${escapeHtml(item.terrain_cost_url)}" target="_blank" rel="noopener">${escapeHtml(terrainLabel)} cost cube NPZ</a>` : ""}</div></section>`;
}

function caseMain(item, rows) {
  const comparison = comparisonConfig();
  const index = Math.max(0, rows.findIndex((row) => row.matchup_id === item.matchup_id));
  const previous = rows.length ? rows[(index - 1 + rows.length) % rows.length] : item;
  const next = rows.length ? rows[(index + 1) % rows.length] : item;
  const reference = item.candidate_by_id[comparison.reference_id];
  const candidate = item.candidate_by_id[comparison.id];
  const body = state.detail === "comparison" ? comparisonTab(item) : state.detail === "raw" ? rawTab(item) : evidenceTab(item);
  return `<main class="case-main"><header class="case-header">
    <div class="case-head-top"><div><h2>${escapeHtml(item.site)}</h2><div class="matchup">${escapeHtml(item.matchup_id)}</div></div><div class="case-actions"><button class="icon-button" data-case="${escapeHtml(previous.matchup_id)}" title="Previous case" aria-label="Previous case">←</button><button class="icon-button" data-case="${escapeHtml(next.matchup_id)}" title="Next case" aria-label="Next case">→</button></div></div>
    <div class="metric-strip case-metrics">
      ${metric("AERONET", fmt(item.truth, 4), `EE ${fmt(item.ee_low)} to ${fmt(item.ee_high)}`)}
      ${metric(comparison.reference_label, fmt(reference.value, 4), `${signed(reference.error_over_ee, 2)} x EE`, reference.within_ee ? "hit" : "miss")}
      ${metric(comparison.label, fmt(candidate.value, 4), `${signed(candidate.error_over_ee, 2)} x EE`, candidate.within_ee ? "hit" : "miss")}
      ${metric("Staged MAIAC", fmt(item.maiac, 4), `sigma ${fmt(item.maiac_unc, 4)}`)}
      ${metric("Surface minimum", fmt(item.surface_min, 4), `band spread ${fmt(item.band_min_spread, 4)}`)}
      ${metric("Cloud / support", cloudPct(item.cloud_fraction), `${pct(item.valid_fraction)} valid`)}
    </div>
    <div class="case-tags">${statusPill(item.transition)}<span><strong>Direction</strong> ${escapeHtml(item.direction)}</span><span><strong>Cost / band</strong> ${fmt(item.cost_per_band, 5)}</span><span><strong>Curvature</strong> ${fmt(item.curvature, 5)}</span><span><strong>AERONET n</strong> ${fmt(item.aeronet_count, 0)}</span><span><strong>Angstrom</strong> ${fmt(item.angstrom, 3)}</span></div>
  </header><nav class="tabs" aria-label="Case details"><button data-detail="evidence" class="${state.detail === "evidence" ? "active" : ""}">Spatial, spectral & surface</button><button data-detail="comparison" class="${state.detail === "comparison" ? "active" : ""}">All outputs</button><button data-detail="raw" class="${state.detail === "raw" ? "active" : ""}">Raw diagnostics</button></nav><div class="case-body">${body}</div></main>`;
}

function casesView() {
  const rows = filteredCases();
  const item = selectedCase(rows);
  if (!rows.some((row) => row.matchup_id === item.matchup_id) && rows.length) state.selectedId = rows[0].matchup_id;
  const selected = selectedCase(rows);
  return `<div class="case-layout">${caseSidebar(rows)}${caseMain(selected, rows)}</div>`;
}

function experimentsView() {
  const inventory = state.data.experiment_inventory || state.data.variants;
  const complete = inventory.filter((row) => row.complete);
  const partial = inventory.filter((row) => !row.complete);
  const best = complete.reduce((left, right) => (right.metrics.within_ee_rate > left.metrics.within_ee_rate ? right : left));
  return `<main class="page"><div class="page-head"><div><h2>Physical experiment ledger</h2><p>Exhaustive saved-directory inventory. It includes end-to-end runs, cost-cube replays, partial targeted diagnostics, and mixed-version comparisons; a high score alone does not make those categories equivalent.</p></div><div class="page-head-actions"><a class="button-link" href="downloads/experiment-inventory.csv">Full inventory CSV</a><a class="button-link" href="downloads/variants.csv">Curated variants CSV</a></div></div>
    <div class="metric-strip">
      ${metric("Saved arms", inventory.length, `${complete.length} with 36 outputs`)}
      ${metric("Highest saved", `${best.metrics.within_ee}/36`, best.label)}
      ${metric("Highest saved rate", pct(best.metrics.within_ee_rate), "category-neutral inventory", best.metrics.within_ee_rate >= .87 ? "met" : "not-met")}
      ${metric("Fresh current code", `${state.data.variants.find((row) => row.id === "fresh").metrics.within_ee}/36`, "cost-cube replay")}
      ${metric("Incomplete arms", partial.length, "explicitly retained")}
      ${metric("Holdout", "NOT SCORED", "folds 2 and 3", "locked")}
    </div>
    <section class="section"><div class="section-head"><div><h3>All discovered variants</h3><p>Status counts distinguish missing or invalid algorithmic outputs. Family and inventory kind remain visible so mixed diagnostics are not presented as production-equivalent runs.</p></div></div>${experimentTable(variantRows(inventory))}</section>
    <section class="section"><div class="section-head"><div><h3>Performance ranking image</h3><p>Complete variants only; zero-based bar axis and fixed 87% reference.</p></div></div><figure class="figure"><img class="zoomable" src="assets/variant-ranking.png" alt="Complete experiment ranking"><figcaption>Within-EE rate on 36 locked medium-development cases.</figcaption></figure></section>
  </main>`;
}

function reproductionView() {
  const reproduction = state.data.reproduction;
  const implementation = reproduction.implementation_files.map((item) => `<a href="${escapeHtml(item.url)}" target="_blank" rel="noopener"><strong>${escapeHtml(item.path)}</strong><span>sha256 ${escapeHtml(item.sha256)}; ${item.bytes} bytes</span></a>`).join("");
  const formula = reproduction.formula || `EE(A) = 0.05 + 0.15 A
within_EE = |AOD_SIAC - AOD_AERONET| <= EE(AOD_AERONET)

J_surface(tau, p) = sum_b [(rho_BOA,b(tau,p) - rho_surface,b(p)) / sigma_b(p)]^2
J_pool(tau, p) = median of J_surface(tau,q) for q in the centred 20 x 20 window at p
J_total(tau, p) = J_pool(tau,p) + [(tau - tau_MAIAC(p)) / sigma_MAIAC(p)]^2
tau_hat(p) = argmin_tau J_total(tau,p)
scene_output = mean of finite tau_hat over the retrieval AOI

Conflict diagnostic: ${reproduction.conflict_rule}`;
  const calibration = reproduction.calibration_formula || `B02 correction = -0.0003 + 0.0243 * anchor_AOD
B03 correction = -0.0006 + 0.0235 * anchor_AOD
B04 correction = -0.0011 + 0.0223 * anchor_AOD`;
  return `<main class="page"><div class="page-head"><div><h2>Production reproduction contract</h2><p>Input, calibration dependencies, surface construction, radiative transfer, objective, output extraction, and acceptance split.</p></div></div>
    <section class="section" style="margin-top:0"><div class="section-head"><div><h3>Operational constraints</h3></div></div><div class="constraints">${reproduction.constraints.map((item) => `<span>${escapeHtml(item)}</span>`).join("")}</div></section>
    <section class="section"><div class="section-head"><div><h3>End-to-end pipeline</h3><p>The current method is reproduced exactly; experimental deviations are recorded in the experiment ledger.</p></div></div><div class="pipeline">${reproduction.pipeline.map((item, index) => `<article><span class="step">${index + 1}</span><h3>${escapeHtml(item.stage)}</h3><p>${escapeHtml(item.detail)}</p></article>`).join("")}</div></section>
    <section class="section"><div class="section-head"><div><h3>Objective and acceptance</h3><p>Symbols below describe the saved implementation at each 60 m solver pixel p and visible band b.</p></div></div>
      <pre class="formula">${escapeHtml(formula)}</pre>
      <p class="muted" style="font-size:10px;line-height:1.55">${escapeHtml(reproduction.acceptance)}</p>
    </section>
    <section class="section"><div class="section-head"><div><h3>${escapeHtml(reproduction.calibration_title || "Current calibration dependency")}</h3><p>${escapeHtml(reproduction.calibration_description || "The coefficient-free OOF and pooled-surface arms remove the three visible-band terms below; their measured results remain in the experiment ledger.")}</p></div></div><pre class="formula">${escapeHtml(calibration)}</pre></section>
    <section class="section"><div class="section-head"><div><h3>Resolved experiment parameters</h3><p>Machine-readable values used by the fresh 36-case replay.</p></div></div>${rawPanel("Parameters", reproduction.parameters)}</section>
    <section class="section"><div class="section-head"><div><h3>Exact implementation snapshot</h3><p>Public copies of the radiometry, surface provider, predictor, solver, LUT backend, Rust pooling kernel, preset, and submission environment. Hashes are computed at report build time.</p></div><a class="button-link" href="downloads/implementation/manifest.json">Manifest</a></div><div class="source-grid">${implementation}</div></section>
  </main>`;
}

function auditView() {
  const data = state.data;
  const sourceRows = Object.entries(data.sources).map(([key, href]) => `<a href="${escapeHtml(href)}" target="_blank" rel="noopener"><strong>${escapeHtml(key.replaceAll("_", " "))}</strong><span>${escapeHtml(href)}</span></a>`).join("");
  const unresolvedJobs = data.jobs.filter((job) => job.final_state === "UNRESOLVED");
  const repairedJobs = data.jobs.filter((job) => job.final_state === "REPAIRED");
  const jobRows = data.jobs.map((job) => `<tr><td class="mono">${escapeHtml(job.job_id)}</td><td class="mono">${escapeHtml(Object.entries(job.states).map(([key, count]) => `${key} ${count}`).join(", ") || "not returned")}</td><td class="mono">${escapeHtml(Object.entries(job.exit_codes).map(([key, count]) => `${key} ${count}`).join(", ") || "not returned")}</td><td class="mono">${escapeHtml((job.repair_job_ids || []).join(", ") || "-")}</td><td class="${job.final_state === "UNRESOLVED" ? "miss" : "hit"}">${escapeHtml(job.final_state)}</td></tr>`).join("");
  const variantStatus = data.variants.map((row) => `<tr><td>${escapeHtml(row.label)}</td><td class="mono">${escapeHtml(Object.entries(row.status_counts).map(([key, count]) => `${key} ${count}`).join(", "))}</td><td class="mono">${row.metrics.n}/36</td><td>${row.complete ? '<span class="hit">Complete</span>' : '<span class="locked">Partial diagnostic</span>'}</td></tr>`).join("");
  return `<main class="page"><div class="page-head"><div><h2>Run and evidence audit</h2><p>Artifact completeness, result status, batch accounting, split lock, exclusions, and direct source files.</p></div><span class="mono muted">generated ${escapeHtml(data.generated_at)}</span></div>
    <div class="metric-strip">
      ${metric("Cases", data.cases.length, "36 expected")}
      ${metric("Case figures", data.audit_counts?.case_figures ?? 3 * data.cases.length, "spatial + spectral + surface")}
      ${metric("Cost cubes", data.audit_counts?.cost_cubes ?? data.cases.length, `${data.cases.length} expected`)}
      ${metric("Unresolved jobs", unresolvedJobs.length, `${repairedJobs.length} historical arrays repaired`, unresolvedJobs.length ? "miss" : "hit")}
      ${metric("Excluded extremes", data.cohort.excluded_extreme_ids.length, "AOD above 0.9")}
      ${metric("Holdout", "LOCKED", data.cohort.holdout_status, "locked")}
    </div>
    <section class="section"><div class="section-head"><div><h3>Result completeness by experiment</h3><p>Partial diagnostic arms are visible but excluded from complete-candidate ranking.</p></div></div><div class="table-wrap"><table class="table"><thead><tr><th>Experiment</th><th>Saved statuses</th><th>Comparable OK</th><th>Classification</th></tr></thead><tbody>${variantStatus}</tbody></table></div></section>
    <section class="section"><div class="section-head"><div><h3>Scheduler accounting</h3><p>Historical task states are retained. REPAIRED means every failed task was rerun by the listed submission and the required result directory is complete; unresolved failures remain zero.</p></div></div><div class="table-wrap"><table class="table"><thead><tr><th>Job</th><th>Historical states</th><th>Exit codes</th><th>Repair job</th><th>Final state</th></tr></thead><tbody>${jobRows}</tbody></table></div></section>
    ${data.limitations?.length ? `<section class="section"><div class="section-head"><div><h3>Limits and unresolved scientific questions</h3><p>These are evidence boundaries, not hidden decisions.</p></div></div><div class="constraints">${data.limitations.map((item) => `<span>${escapeHtml(item)}</span>`).join("")}</div></section>` : ""}
    <section class="section"><div class="section-head"><div><h3>Split and deferred extreme cases</h3><p>Site-group split seed ${data.cohort.split_seed}; holdout folds ${data.cohort.holdout_folds.join(", ")}. The three deferred records are not part of the medium-development analysis.</p></div></div><div class="constraints">${data.cohort.excluded_extreme_ids.map((item) => `<span class="mono">${escapeHtml(item)}</span>`).join("")}</div></section>
    <section class="section"><div class="section-head"><div><h3>Source index</h3><p>Direct paths to the bounded report extracts, saved retrievals, cost cubes, and cohort list.</p></div></div><div class="source-grid">${sourceRows}</div></section>
  </main>`;
}

function render() {
  let body;
  if (state.view === "surface") body = surfaceView();
  else if (state.view === "cases") body = casesView();
  else if (state.view === "experiments") body = experimentsView();
  else if (state.view === "reproduction") body = reproductionView();
  else if (state.view === "audit") body = auditView();
  else body = overviewView();
  app.innerHTML = reportShell(body);
  bindEvents();
}

function bindEvents() {
  document.querySelectorAll("[data-view]").forEach((button) => button.addEventListener("click", () => {
    state.view = button.dataset.view;
    window.location.hash = state.view;
    render();
    window.scrollTo({ top: 0, behavior: "auto" });
  }));
  document.querySelectorAll("[data-case]").forEach((button) => button.addEventListener("click", () => {
    state.selectedId = button.dataset.case;
    render();
    if (window.innerWidth < 760) window.scrollTo({ top: 0, behavior: "auto" });
  }));
  document.querySelectorAll("[data-detail]").forEach((button) => button.addEventListener("click", () => {
    state.detail = button.dataset.detail;
    render();
  }));
  const controls = ["transition", "direction", "outcome", "sort"];
  controls.forEach((id) => document.getElementById(id)?.addEventListener("change", (event) => {
    state[id] = event.target.value;
    const rows = filteredCases();
    if (rows.length && !rows.some((item) => item.matchup_id === state.selectedId)) state.selectedId = rows[0].matchup_id;
    render();
  }));
  const search = document.getElementById("search");
  search?.addEventListener("input", (event) => {
    const position = event.target.selectionStart;
    state.search = event.target.value;
    const rows = filteredCases();
    if (rows.length && !rows.some((item) => item.matchup_id === state.selectedId)) state.selectedId = rows[0].matchup_id;
    render();
    const replacement = document.getElementById("search");
    replacement?.focus();
    replacement?.setSelectionRange(position, position);
  });
  const dialog = document.getElementById("image-dialog");
  document.querySelectorAll(".zoomable").forEach((image) => image.addEventListener("click", () => {
    document.getElementById("dialog-image").src = image.src;
    document.getElementById("dialog-label").textContent = image.alt;
    dialog.showModal();
  }));
  document.getElementById("dialog-close")?.addEventListener("click", () => dialog.close());
  dialog?.addEventListener("click", (event) => {
    if (event.target === dialog) dialog.close();
  });
}

window.addEventListener("hashchange", () => {
  const requested = window.location.hash.replace("#", "");
  if (validViews().includes(requested) && requested !== state.view) {
    state.view = requested;
    render();
    window.scrollTo({ top: 0, behavior: "auto" });
  }
});

async function start() {
  try {
    const response = await fetch("data/report.json", { cache: "no-store" });
    if (!response.ok) throw new Error(`HTTP ${response.status}`);
    state.data = await response.json();
    const comparison = comparisonConfig();
    transitionLabels.gain = `Gain vs ${comparison.reference_short_label}`;
    transitionLabels.loss = `Loss vs ${comparison.reference_short_label}`;
    const requested = window.location.hash.replace("#", "");
    if (validViews().includes(requested)) state.view = requested;
    state.selectedId = [...state.data.cases].sort((a, b) => Math.abs(b.fresh_error_over_ee) - Math.abs(a.fresh_error_over_ee))[0]?.matchup_id;
    render();
  } catch (error) {
    app.innerHTML = `<div class="empty"><strong>Report data could not be loaded.</strong><div class="mono" style="margin-top:8px">${escapeHtml(error.message)}</div></div>`;
  }
}

start();
