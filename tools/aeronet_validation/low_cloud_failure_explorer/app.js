const app = document.getElementById("app");

const state = {
  data: null,
  view: "cases",
  detailTab: "spatial",
  selectedId: null,
  search: "",
  group: "all",
  direction: "all",
  severity: "all",
  transition: "all",
  unresolvedOnly: false,
  sort: "severity",
  galleryMode: "spatial",
};

const groupColors = {
  A: "#1769aa",
  B: "#6d4c41",
  C: "#7b1fa2",
  D: "#00838f",
  E: "#ad6800",
  F: "#546e7a",
  G: "#9e3d64",
};

const escapeHtml = (value) => String(value ?? "")
  .replaceAll("&", "&amp;")
  .replaceAll("<", "&lt;")
  .replaceAll(">", "&gt;")
  .replaceAll('"', "&quot;")
  .replaceAll("'", "&#039;");

const fmt = (value, digits = 3) => {
  if (value === null || value === undefined || !Number.isFinite(Number(value))) return "NA";
  return Number(value).toFixed(digits);
};

const pct = (value, digits = 1) => {
  if (value === null || value === undefined || !Number.isFinite(Number(value))) return "NA";
  return `${(100 * Number(value)).toFixed(digits)}%`;
};

const signed = (value, digits = 3) => {
  if (value === null || value === undefined || !Number.isFinite(Number(value))) return "NA";
  const number = Number(value);
  return `${number >= 0 ? "+" : ""}${number.toFixed(digits)}`;
};

function filteredCases() {
  const query = state.search.trim().toLowerCase();
  const rows = state.data.cases.filter((item) => {
    if (query && !`${item.matchup_id} ${item.site}`.toLowerCase().includes(query)) return false;
    if (state.group !== "all" && item.mechanism_code !== state.group) return false;
    if (state.direction !== "all" && item.direction !== state.direction) return false;
    if (state.severity !== "all" && item.severity !== state.severity) return false;
    if (state.transition !== "all" && item.baseline_transition !== state.transition) return false;
    if (state.unresolvedOnly && item.tested_one_prior_recoverable) return false;
    return true;
  });
  return rows.sort((a, b) => {
    if (state.sort === "truth") return b.truth_aod - a.truth_aod;
    if (state.sort === "spread") return b.band_spread - a.band_spread;
    if (state.sort === "cloud") return b.cloud_fraction - a.cloud_fraction;
    if (state.sort === "cost") return b.cost_per_band - a.cost_per_band;
    if (state.sort === "backstop") return Math.abs(b.backstop_shift) - Math.abs(a.backstop_shift);
    if (state.sort === "site") return a.site.localeCompare(b.site);
    return a.severity_rank - b.severity_rank;
  });
}

function selectedCase() {
  const rows = state.data.cases;
  return rows.find((item) => item.matchup_id === state.selectedId) || rows[0];
}

function groupOptions() {
  return Object.entries(state.data.group_definitions)
    .map(([code, item]) => `<option value="${code}" ${state.group === code ? "selected" : ""}>${code}. ${escapeHtml(item.label)}</option>`)
    .join("");
}

function optionValues(field) {
  return [...new Set(state.data.cases.map((item) => item[field]))];
}

function sidebarHtml(rows) {
  const severityOptions = optionValues("severity").map((value) => `<option value="${escapeHtml(value)}" ${state.severity === value ? "selected" : ""}>${escapeHtml(value)}</option>`).join("");
  const transitionOptions = optionValues("baseline_transition").map((value) => `<option value="${escapeHtml(value)}" ${state.transition === value ? "selected" : ""}>${escapeHtml(value)}</option>`).join("");
  const list = rows.map((item) => `
    <button class="case-row ${item.matchup_id === state.selectedId ? "active" : ""}" data-case="${escapeHtml(item.matchup_id)}">
      <span class="case-rank">${String(item.severity_rank).padStart(2, "0")}</span>
      <span class="case-name">
        <strong>${escapeHtml(item.site)}</strong>
        <span>AERONET ${fmt(item.truth_aod)} · current ${fmt(item.retrieved_aod)}</span>
      </span>
      <span class="case-score">
        <span class="group-dot" style="background:${groupColors[item.mechanism_code] || "#667"}">${item.mechanism_code}</span>
        <strong>${fmt(item.ee_ratio, 2)}× EE</strong>
      </span>
    </button>`).join("");
  return `
    <aside class="case-sidebar">
      <div class="filter-panel">
        <label class="wide" for="search">Matchup or site</label>
        <input class="wide" id="search" type="search" value="${escapeHtml(state.search)}" placeholder="Search">
        <div class="filter-grid">
          <div><label for="group">Rule-based group</label><select id="group"><option value="all">All groups</option>${groupOptions()}</select></div>
          <div><label for="direction">Direction</label><select id="direction"><option value="all">Both</option><option ${state.direction === "Underread" ? "selected" : ""}>Underread</option><option ${state.direction === "Overread" ? "selected" : ""}>Overread</option></select></div>
          <div><label for="severity">Error / EE</label><select id="severity"><option value="all">All ranges</option>${severityOptions}</select></div>
          <div><label for="sort">Sort</label><select id="sort"><option value="severity" ${state.sort === "severity" ? "selected" : ""}>Error / EE</option><option value="backstop" ${state.sort === "backstop" ? "selected" : ""}>Backstop shift</option><option value="truth" ${state.sort === "truth" ? "selected" : ""}>AERONET AOD</option><option value="spread" ${state.sort === "spread" ? "selected" : ""}>Band spread</option><option value="cost" ${state.sort === "cost" ? "selected" : ""}>Cost / band</option><option value="cloud" ${state.sort === "cloud" ? "selected" : ""}>Cloud fraction</option><option value="site" ${state.sort === "site" ? "selected" : ""}>Site</option></select></div>
          <div class="wide"><label for="transition">Historical transition</label><select id="transition"><option value="all">All transitions</option>${transitionOptions}</select></div>
          <label class="check-row wide"><input id="unresolved" type="checkbox" ${state.unresolvedOnly ? "checked" : ""}> Unresolved by tested one-prior outputs</label>
        </div>
        <div class="filter-summary">${rows.length} of ${state.data.failure_count} cases</div>
      </div>
      <div class="case-list">${list || '<div class="empty-state">No cases match these filters.</div>'}</div>
    </aside>`;
}

function metric(label, value, className = "") {
  return `<div class="metric"><div class="metric-label">${escapeHtml(label)}</div><div class="metric-value ${className}">${escapeHtml(value)}</div></div>`;
}

function outputSvg(item) {
  const rows = item.candidates.filter((candidate) => candidate.value !== null && Number.isFinite(candidate.value));
  const truth = item.truth_aod;
  const ee = item.ee_threshold;
  const maxValue = Math.max(1, truth + ee, ...rows.map((row) => row.value)) * 1.08;
  const labelWidth = 220;
  const right = 30;
  const width = 920;
  const rowHeight = 28;
  const top = 34;
  const height = top + rows.length * rowHeight + 34;
  const plotWidth = width - labelWidth - right;
  const x = (value) => labelWidth + Math.max(0, Math.min(maxValue, value)) / maxValue * plotWidth;
  const grid = [0, .25, .5, .75, 1].map((fraction) => {
    const value = maxValue * fraction;
    return `<line x1="${x(value)}" x2="${x(value)}" y1="20" y2="${height - 26}" stroke="var(--line)"/><text x="${x(value)}" y="${height - 8}" text-anchor="middle" fill="var(--muted)" font-size="10">${value.toFixed(2)}</text>`;
  }).join("");
  const marks = rows.map((row, index) => {
    const y = top + index * rowHeight;
    const colour = row.within_ee ? "var(--positive)" : "var(--muted)";
    return `<text x="8" y="${y + 4}" fill="var(--ink)" font-size="10">${escapeHtml(row.label)}</text><line x1="${x(truth)}" x2="${x(row.value)}" y1="${y}" y2="${y}" stroke="var(--line-strong)"/><circle cx="${x(row.value)}" cy="${y}" r="4" fill="${colour}"/><text x="${Math.min(width - 34, x(row.value) + 7)}" y="${y + 4}" fill="var(--muted)" font-size="9">${row.value.toFixed(3)}</text>`;
  }).join("");
  return `<svg class="output-svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="Saved scalar AOD outputs for ${escapeHtml(item.site)}"><rect x="${x(truth - ee)}" y="20" width="${x(truth + ee) - x(truth - ee)}" height="${height - 46}" fill="var(--surface-2)"/><line x1="${x(truth)}" x2="${x(truth)}" y1="18" y2="${height - 26}" stroke="var(--ink)" stroke-width="1.4"/><text x="${x(truth) + 4}" y="14" fill="var(--ink)" font-size="10">AERONET</text>${grid}${marks}</svg>`;
}

function outputsTable(item) {
  return `<table class="output-table"><thead><tr><th>Type</th><th>Saved output</th><th>AOD</th><th>Within EE</th></tr></thead><tbody>${item.candidates.map((row) => `<tr><td>${escapeHtml(row.kind)}</td><td>${escapeHtml(row.label)}</td><td class="mono">${fmt(row.value)}</td><td class="${row.within_ee ? "status-hit" : "status-miss"}">${row.within_ee === null ? "NA" : row.within_ee ? "Yes" : "No"}</td></tr>`).join("")}</tbody></table>`;
}

function kvRows(entries) {
  return entries.map(([label, value]) => `<tr><th>${escapeHtml(label)}</th><td class="mono wrap">${escapeHtml(value)}</td></tr>`).join("");
}

function spatialTab(item) {
  return `<section class="evidence-section"><div class="section-heading"><h3>Spatial evidence grid</h3><span>60 m solver grid · dashed circle is the 1.5 km station window</span></div><figure class="evidence-figure"><img class="zoomable" src="${item.spatial_image}" alt="Spatial evidence mosaic for ${escapeHtml(item.site)}" loading="eager"><figcaption><span>TOA, surface prior, solver CAMS backstop, operational and no-backstop AOD, backstop-effect map, per-band minima, and residuals.</span><a href="${item.spatial_image}" target="_blank" rel="noopener">Open PNG</a></figcaption></figure></section>`;
}

function spectralTab(item) {
  return `<section class="evidence-section"><div class="section-heading"><h3>Spectral, cost, and scalar-output evidence</h3><span>${item.cube.aod_nodes} AOD nodes · ${item.curves.curve_support_count} pixels · ${escapeHtml(item.curves.curve_support_scope)}</span></div><figure class="evidence-figure"><img class="zoomable" src="${item.diagnostic_image}" alt="Spectral and cost diagnostics for ${escapeHtml(item.site)}" loading="eager"><figcaption><span>Scene and selected-support total cost, per-band costs and residuals, reflectance spectrum, node closure, and scalar outputs.</span><a href="${item.diagnostic_image}" target="_blank" rel="noopener">Open PNG</a></figcaption></figure></section>`;
}

function outputsTab(item) {
  return `<section class="evidence-section"><div class="section-heading"><h3>Saved outputs and diagnostic anchors</h3><span>EE band is shaded; green means the saved scalar lies inside EE</span></div>${outputSvg(item)}<div style="height:18px"></div>${outputsTable(item)}</section>`;
}

function rawTab(item) {
  const extraction = item.retrieval_extraction || {};
  const solver = item.canonical_solver || {};
  const cube = item.cube || {};
  return `<section class="evidence-section"><div class="section-heading"><h3>Raw diagnostic values and source files</h3><span>Canonical scene-mean result is kept separate from the dump-only rerun</span></div><div class="data-grid">
    <div class="data-panel"><h3>Retrieval and extraction</h3><table class="kv-table">${kvRows([
      ["AERONET AOD", fmt(item.truth_aod, 6)], ["Current scene mean", fmt(item.retrieved_aod, 6)], ["Station pixel", fmt(extraction.station, 6)], ["1.5 km median", fmt(extraction.winmed, 6)], ["Dump-only rerun", fmt(item.diagnostic_retrieved, 6)], ["Signed error", signed(item.error, 6)], ["EE threshold", fmt(item.ee_threshold, 6)], ["Error / EE", fmt(item.ee_ratio, 4)], ["OOF audit AOD", fmt(item.oof_aod, 6)], ["Historical R2", fmt(item.baseline_aod, 6)]
    ])}</table></div>
    <div class="data-panel"><h3>Cost and spatial support</h3><table class="kv-table">${kvRows([
      ["Cost mode", solver.surface_cost_mode], ["Final cost / band", fmt(item.cost_per_band, 6)], ["Global curve minimum", fmt(item.curve_min_aod, 6)], ["Relative second-node delta", fmt(item.curve_relative_second_delta, 6)], ["Band argmin spread", fmt(item.band_spread, 6)], ["Valid support", pct(item.valid_support_fraction, 2)], ["Station-window valid pixels", item.curves.station_window_valid_count], ["Curve support", item.curves.curve_support_scope], ["Curve support pixels", item.curves.curve_support_count], ["Pooled station-window median", fmt(item.curves.pooled_local_median, 6)], ["Pooled station-window mean", fmt(item.curves.pooled_local_mean, 6)], ["Pooled station-window uncertainty", fmt(item.curves.pooled_local_unc_median, 6)], ["Cube shape", cube.shape?.join(" × ")], ["Pool window / minimum", `${cube.pool_window} / ${cube.pool_min_count}`]
    ])}</table></div>
    <div class="data-panel"><h3>Spectral scalars</h3><table class="kv-table">${kvRows([
      ["Surface-prior anchor CAMS", fmt(item.surface_anchor_cams_aod, 6)], ["Solver CAMS backstop", fmt(item.solver_cams_aod, 6)], ["Solver CAMS sigma", fmt(item.solver_cams_unc, 6)], ["No-backstop scene mean", fmt(item.no_backstop_aod, 6)], ["Backstop shift", signed(item.backstop_shift, 6)], ["B02 cost minimum", fmt(item.band_b02_min_aod, 6)], ["B03 cost minimum", fmt(item.band_b03_min_aod, 6)], ["B04 cost minimum", fmt(item.band_b04_min_aod, 6)], ["Prior B02", fmt(item.prior_boa.B02, 6)], ["Prior B03", fmt(item.prior_boa.B03, 6)], ["Prior B04", fmt(item.prior_boa.B04, 6)], ["Prior uncertainty B02", fmt(item.prior_unc.B02, 6)], ["Prior uncertainty B03", fmt(item.prior_unc.B03, 6)], ["Prior uncertainty B04", fmt(item.prior_unc.B04, 6)]
    ])}</table></div>
    <div class="data-panel"><h3>Rule-based filter fields</h3><table class="kv-table">${kvRows([
      ["Group", `${item.mechanism_code}. ${item.mechanism}`], ["Direction", item.direction], ["Severity range", item.severity], ["Historical transition", item.baseline_transition], ["Any band inside EE", item.any_band_within_ee], ["All bands below EE", item.all_bands_low], ["All bands above EE", item.all_bands_high], ["Surface-anchor CAMS below / above EE", `${item.cams_low} / ${item.cams_high}`], ["Global curve inside EE", item.curve_min_within_ee], ["Tested-output recoverable", item.tested_one_prior_recoverable]
    ])}</table></div>
  </div><div class="source-links"><a href="${item.canonical_result_url}" target="_blank" rel="noopener">Canonical result JSON</a><a href="${item.diagnostic_result_url}" target="_blank" rel="noopener">Dump-only result JSON</a><a href="${item.cost_cube_url}" target="_blank" rel="noopener">Cost cube NPZ</a></div></section>`;
}

function detailBody(item) {
  if (state.detailTab === "spectral") return spectralTab(item);
  if (state.detailTab === "outputs") return outputsTab(item);
  if (state.detailTab === "raw") return rawTab(item);
  return spatialTab(item);
}

function caseMainHtml(item, rows) {
  const index = Math.max(0, rows.findIndex((row) => row.matchup_id === item.matchup_id));
  const previous = rows[(index - 1 + rows.length) % rows.length];
  const next = rows[(index + 1) % rows.length];
  return `<main class="case-main">
    <header class="case-header">
      <div class="case-header-top"><div><h2>${escapeHtml(item.site)}</h2><div class="matchup">${escapeHtml(item.matchup_id)}</div></div><div class="case-actions"><button class="icon-button" data-case="${escapeHtml(previous?.matchup_id)}" title="Previous case" aria-label="Previous case">←</button><button class="icon-button" data-case="${escapeHtml(next?.matchup_id)}" title="Next case" aria-label="Next case">→</button></div></div>
      <div class="metric-row">
        ${metric("AERONET", fmt(item.truth_aod))}${metric("Current", fmt(item.retrieved_aod))}${metric("Signed error", signed(item.error), item.error < 0 ? "negative" : "positive")}${metric("Absolute error / EE", `${fmt(item.ee_ratio, 2)}×`, "negative")}${metric("Cloud", pct(item.cloud_fraction))}${metric("Valid support", pct(item.valid_support_fraction))}
      </div>
      <div class="signature-row"><span class="tag"><strong>Group ${item.mechanism_code}</strong>${escapeHtml(item.mechanism)}</span><span class="tag"><strong>Direction</strong>${item.direction}</span><span class="tag"><strong>Solver CAMS</strong>${fmt(item.solver_cams_aod)} ± ${fmt(item.solver_cams_unc)}</span><span class="tag"><strong>Backstop shift</strong>${signed(item.backstop_shift)}</span><span class="tag"><strong>Historical</strong>${escapeHtml(item.baseline_transition)}</span><span class="tag"><strong>One-prior union</strong>${item.tested_one_prior_recoverable ? "has a truth-selected hit" : "no tested hit"}</span></div>
    </header>
    <nav class="detail-tabs" aria-label="Case evidence views">
      <button data-tab="spatial" class="${state.detailTab === "spatial" ? "active" : ""}">Spatial</button>
      <button data-tab="spectral" class="${state.detailTab === "spectral" ? "active" : ""}">Spectral & cost</button>
      <button data-tab="outputs" class="${state.detailTab === "outputs" ? "active" : ""}">Saved outputs</button>
      <button data-tab="raw" class="${state.detailTab === "raw" ? "active" : ""}">Raw diagnostics</button>
    </nav>
    <div class="detail-body">${detailBody(item)}</div>
  </main>`;
}

function casesView() {
  const rows = filteredCases();
  let item = selectedCase();
  if (rows.length && !rows.some((row) => row.matchup_id === item.matchup_id)) {
    item = rows[0];
    state.selectedId = item.matchup_id;
  }
  return `<div class="explorer-layout">${sidebarHtml(rows)}${rows.length ? caseMainHtml(item, rows) : '<main class="case-main"><div class="empty-state">No cases match these filters.</div></main>'}</div>`;
}

function overviewView() {
  const stats = state.data.stats;
  const groups = Object.entries(state.data.group_definitions).map(([code, item]) => `<tr><td><span class="group-dot" style="background:${groupColors[code]}">${code}</span></td><td>${escapeHtml(item.label)}</td><td class="mono">${item.count}</td><td class="mono">${item.underreads} / ${item.overreads}</td><td>${escapeHtml(item.evidence)}</td></tr>`).join("");
  const sources = state.data.sources;
  return `<main class="page-view"><header class="page-view-header"><h2>Cross-case evidence</h2><p>All ${state.data.cohort_count} retrievals provide context; detailed matrices contain every outside-EE case.</p></header>
    <div class="metric-row overview-metrics">${metric("Status OK", stats.status_ok)}${metric("Within EE", stats.within_ee)}${metric("Within-EE rate", pct(stats.within_ee_rate))}${metric("Outside EE", stats.outside_ee, "negative")}${metric("Underreads", stats.underreads)}${metric("Overreads", stats.overreads)}</div>
    <div class="overview-grid"><figure class="evidence-figure"><img class="zoomable" src="assets/cross_case_overview.png" alt="Cross-case AOD diagnostics"><figcaption>Retrieval, normalized error, band spread, and atmospheric-anchor relationships.</figcaption></figure><figure class="evidence-figure"><img class="zoomable" src="assets/evidence_matrix.png" alt="Failure evidence matrix"><figcaption>Every failure ordered by normalized error, with signed AOD evidence and context-field percentiles.</figcaption></figure>
    <section class="data-panel"><h3>Rule-based evidence groups</h3><div class="neutral-note">These groups are reproducible filters over saved diagnostics. They are not physical diagnoses or recommended decisions.</div><table class="group-table"><thead><tr><th>Code</th><th>Filter label</th><th>Cases</th><th>Under / over</th><th>Observed rule evidence</th></tr></thead><tbody>${groups}</tbody></table></section>
    <section class="data-panel"><h3>Scope, definitions, and source snapshots</h3><table class="kv-table">${kvRows([["Snapshot generated", state.data.generated_at], ["Cohort", `${state.data.cohort_count} Sentinel-2/AERONET matchups with image cloud fraction below 20%`], ["Current method", state.data.method], ["Current output", "Scene-mean pooled AOD"], ["Expected error", "absolute retrieval error <= 0.05 + 0.15 x AERONET AOD"], ["Surface-prior anchor CAMS", "Point/scene value used while constructing the S2 SWIR/NIR-anchored surface prior"], ["Solver CAMS backstop", "Exact atmospheric-prior center used in the final Gaussian AOD penalty"], ["No-backstop replay", "Same cost cube, pooling, and scene mean with only the Gaussian backstop disabled"], ["Diagnostic rerun", "Same operational extraction with exact cost-cube dumping enabled"], ["Evaluation reference", "AERONET AOD; shown for investigation and not available to an operational selector"]])}</table><div class="source-links"><a href="${sources.failure_csv}" target="_blank" rel="noopener">Failure CSV</a><a href="${sources.failure_analysis}" target="_blank" rel="noopener">Failure analysis JSON</a><a href="${sources.current_results}" target="_blank" rel="noopener">Current results</a><a href="${sources.diagnostic_results}" target="_blank" rel="noopener">Diagnostic reruns</a><a href="${sources.cost_cubes}" target="_blank" rel="noopener">Cost cubes</a></div></section></div></main>`;
}

function galleryFilters(rows) {
  const severityOptions = optionValues("severity").map((value) => `<option value="${escapeHtml(value)}" ${state.severity === value ? "selected" : ""}>${escapeHtml(value)}</option>`).join("");
  return `<div class="gallery-toolbar">
    <div class="segmented" role="group" aria-label="Gallery evidence type"><button data-gallery-mode="spatial" class="${state.galleryMode === "spatial" ? "active" : ""}">Spatial mosaics</button><button data-gallery-mode="diagnostic" class="${state.galleryMode === "diagnostic" ? "active" : ""}">Spectral & cost</button></div>
    <label><span>Matchup or site</span><input id="search" type="search" value="${escapeHtml(state.search)}" placeholder="Search"></label>
    <label><span>Rule-based group</span><select id="group"><option value="all">All groups</option>${groupOptions()}</select></label>
    <label><span>Direction</span><select id="direction"><option value="all">Both</option><option ${state.direction === "Underread" ? "selected" : ""}>Underread</option><option ${state.direction === "Overread" ? "selected" : ""}>Overread</option></select></label>
    <label><span>Error / EE</span><select id="severity"><option value="all">All ranges</option>${severityOptions}</select></label>
    <label><span>Sort</span><select id="sort"><option value="severity" ${state.sort === "severity" ? "selected" : ""}>Error / EE</option><option value="backstop" ${state.sort === "backstop" ? "selected" : ""}>Backstop shift</option><option value="truth" ${state.sort === "truth" ? "selected" : ""}>AERONET AOD</option><option value="spread" ${state.sort === "spread" ? "selected" : ""}>Band spread</option><option value="cost" ${state.sort === "cost" ? "selected" : ""}>Cost / band</option><option value="cloud" ${state.sort === "cloud" ? "selected" : ""}>Cloud fraction</option><option value="site" ${state.sort === "site" ? "selected" : ""}>Site</option></select></label>
    <label class="check-row"><input id="unresolved" type="checkbox" ${state.unresolvedOnly ? "checked" : ""}> Unresolved in tested one-prior outputs</label>
    <span class="gallery-count">${rows.length} / ${state.data.failure_count}</span>
  </div>`;
}

function galleryView() {
  const rows = filteredCases();
  const imageKey = state.galleryMode === "diagnostic" ? "diagnostic_image" : "spatial_image";
  const evidenceLabel = state.galleryMode === "diagnostic" ? "Spectral and cost" : "Spatial";
  const cards = rows.map((item) => `<article class="gallery-item"><button data-open-case="${escapeHtml(item.matchup_id)}"><img src="${item[imageKey]}" alt="${evidenceLabel} evidence thumbnail for ${escapeHtml(item.site)}" loading="lazy"><div class="gallery-meta"><strong>${escapeHtml(item.site)}</strong><span>Group ${item.mechanism_code} · AERONET ${fmt(item.truth_aod)} · current ${fmt(item.retrieved_aod)} · ${fmt(item.ee_ratio, 2)}× EE</span></div></button></article>`).join("");
  return `<main class="page-view"><header class="page-view-header"><h2>All-case evidence gallery</h2><p>Side-by-side evidence for all outside-EE cases matching the active filters.</p></header>${galleryFilters(rows)}<div class="gallery-grid">${cards || '<div class="empty-state">No cases match these filters.</div>'}</div></main>`;
}

function topbar() {
  return `<header class="topbar"><div class="title-lockup"><h1>${escapeHtml(state.data.title)}</h1><p>${state.data.failure_count} outside-EE cases · ${escapeHtml(state.data.method)}</p></div><nav class="view-nav" aria-label="Workspace views"><button data-view="cases" class="${state.view === "cases" ? "active" : ""}">Case explorer</button><button data-view="overview" class="${state.view === "overview" ? "active" : ""}">Cross-case</button><button data-view="gallery" class="${state.view === "gallery" ? "active" : ""}">Evidence gallery</button></nav></header>`;
}

function render() {
  const body = state.view === "overview" ? overviewView() : state.view === "gallery" ? galleryView() : casesView();
  app.innerHTML = `${topbar()}${body}<dialog class="image-dialog"><div class="dialog-bar"><strong></strong><button class="icon-button" data-close-dialog title="Close" aria-label="Close">×</button></div><img alt="Expanded diagnostic figure"></dialog>`;
  bindEvents();
}

function updateFilter(key, value) {
  state[key] = value;
  render();
  if (key === "search") {
    const search = document.getElementById("search");
    search?.focus();
    search?.setSelectionRange(value.length, value.length);
  }
}

function openImage(src, alt) {
  const dialog = document.querySelector(".image-dialog");
  dialog.querySelector("strong").textContent = alt;
  const image = dialog.querySelector("img");
  image.src = src;
  image.alt = alt;
  dialog.showModal();
}

function selectCase(matchupId) {
  if (!matchupId) return;
  state.selectedId = matchupId;
  state.view = "cases";
  history.replaceState(null, "", `#case=${encodeURIComponent(matchupId)}`);
  render();
}

function bindEvents() {
  document.querySelectorAll("[data-view]").forEach((button) => button.addEventListener("click", () => { state.view = button.dataset.view; render(); }));
  document.querySelectorAll("[data-case]").forEach((button) => button.addEventListener("click", () => selectCase(button.dataset.case)));
  document.querySelectorAll("[data-open-case]").forEach((button) => button.addEventListener("click", () => selectCase(button.dataset.openCase)));
  document.querySelectorAll("[data-gallery-mode]").forEach((button) => button.addEventListener("click", () => { state.galleryMode = button.dataset.galleryMode; render(); }));
  document.querySelectorAll("[data-tab]").forEach((button) => button.addEventListener("click", () => { state.detailTab = button.dataset.tab; render(); }));
  document.querySelectorAll(".zoomable").forEach((image) => image.addEventListener("click", () => openImage(image.src, image.alt)));
  document.querySelector("[data-close-dialog]")?.addEventListener("click", () => document.querySelector(".image-dialog")?.close());
  document.getElementById("search")?.addEventListener("input", (event) => updateFilter("search", event.target.value));
  document.getElementById("group")?.addEventListener("change", (event) => updateFilter("group", event.target.value));
  document.getElementById("direction")?.addEventListener("change", (event) => updateFilter("direction", event.target.value));
  document.getElementById("severity")?.addEventListener("change", (event) => updateFilter("severity", event.target.value));
  document.getElementById("transition")?.addEventListener("change", (event) => updateFilter("transition", event.target.value));
  document.getElementById("sort")?.addEventListener("change", (event) => updateFilter("sort", event.target.value));
  document.getElementById("unresolved")?.addEventListener("change", (event) => updateFilter("unresolvedOnly", event.target.checked));
}

document.addEventListener("keydown", (event) => {
  if (state.view !== "cases" || ["INPUT", "SELECT", "TEXTAREA"].includes(document.activeElement?.tagName)) return;
  if (event.key !== "ArrowLeft" && event.key !== "ArrowRight") return;
  const rows = filteredCases();
  const index = rows.findIndex((row) => row.matchup_id === state.selectedId);
  if (index < 0 || !rows.length) return;
  const delta = event.key === "ArrowLeft" ? -1 : 1;
  selectCase(rows[(index + delta + rows.length) % rows.length].matchup_id);
});

fetch("data/cases.json")
  .then((response) => {
    if (!response.ok) throw new Error(`HTTP ${response.status}`);
    return response.json();
  })
  .then((data) => {
    state.data = data;
    const match = location.hash.match(/^#case=(.+)$/);
    const requested = match ? decodeURIComponent(match[1]) : null;
    state.selectedId = data.cases.some((item) => item.matchup_id === requested) ? requested : data.cases[0]?.matchup_id;
    render();
  })
  .catch((error) => {
    app.innerHTML = `<div class="empty-state">Unable to load the diagnostic snapshot: ${escapeHtml(error.message)}</div>`;
  });
