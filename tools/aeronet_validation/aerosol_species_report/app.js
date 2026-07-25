const app = document.getElementById("app");

const state = {
  data: null,
  view: "overview",
  detailTab: "spatial",
  selectedId: null,
  search: "",
  transition: "all",
  group: "all",
  baseline: "all",
  sort: "species_effect",
  galleryMode: "spatial",
};

const transitionLabels = {
  gain: "Gain",
  loss: "Loss",
  both_hit: "Both within EE",
  both_miss: "Both outside EE",
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

const signed = (value, digits = 3) => {
  if (value === null || value === undefined || !Number.isFinite(Number(value))) return "NA";
  const number = Number(value);
  return `${number >= 0 ? "+" : ""}${number.toFixed(digits)}`;
};

const pct = (value, digits = 1) => {
  if (value === null || value === undefined || !Number.isFinite(Number(value))) return "NA";
  return `${(100 * Number(value)).toFixed(digits)}%`;
};

const seconds = (value) => {
  if (value === null || value === undefined || !Number.isFinite(Number(value))) return "NA";
  const number = Number(value);
  return number >= 120 ? `${(number / 60).toFixed(1)} min` : `${number.toFixed(1)} s`;
};

function metric(label, value, className = "") {
  return `<div class="metric"><div class="metric-label">${escapeHtml(label)}</div><div class="metric-value ${className}">${escapeHtml(value)}</div></div>`;
}

function kvRows(entries) {
  return entries.map(([label, value]) => `<tr><th>${escapeHtml(label)}</th><td class="mono wrap">${escapeHtml(value)}</td></tr>`).join("");
}

function transitionClass(value) {
  return `transition-${String(value).replaceAll("_", "-")}`;
}

function filteredCases() {
  const query = state.search.trim().toLowerCase();
  return state.data.cases.filter((item) => {
    if (query && !`${item.matchup_id} ${item.site}`.toLowerCase().includes(query)) return false;
    if (state.transition !== "all" && item.transition_cci_vs_fixed !== state.transition) return false;
    if (state.group !== "all" && item.baseline_failure_group !== state.group) return false;
    const baselineHit = item.methods.lut_continental.within_ee;
    if (state.baseline === "hit" && !baselineHit) return false;
    if (state.baseline === "miss" && baselineHit) return false;
    return true;
  }).sort((a, b) => {
    if (state.sort === "truth") return b.truth_aod - a.truth_aod;
    if (state.sort === "cloud") return b.cloud_fraction - a.cloud_fraction;
    if (state.sort === "site") return a.site.localeCompare(b.site);
    if (state.sort === "confidence") return a.species.selection_confidence_median - b.species.selection_confidence_median;
    if (state.sort === "fixed_error") return b.methods.sixs_continental.abs_error - a.methods.sixs_continental.abs_error;
    return a.cci_minus_fixed_abs_error - b.cci_minus_fixed_abs_error;
  });
}

function selectedCase() {
  return state.data.cases.find((item) => item.matchup_id === state.selectedId) || state.data.cases[0];
}

function filterControls() {
  const groups = [...new Set(state.data.cases.map((item) => item.baseline_failure_group))].sort();
  return `<div class="filter-grid">
    <div class="wide"><label for="search">Matchup or site</label><input id="search" type="search" value="${escapeHtml(state.search)}" placeholder="Search"></div>
    <div><label for="transition">CCI-3 vs fixed 6S</label><select id="transition"><option value="all">All transitions</option>${Object.entries(transitionLabels).map(([value, label]) => `<option value="${value}" ${state.transition === value ? "selected" : ""}>${label}</option>`).join("")}</select></div>
    <div><label for="group">Prior baseline group</label><select id="group"><option value="all">All groups</option>${groups.map((value) => `<option value="${value}" ${state.group === value ? "selected" : ""}>${value}</option>`).join("")}</select></div>
    <div><label for="baseline">LUT baseline status</label><select id="baseline"><option value="all">Hit and miss</option><option value="hit" ${state.baseline === "hit" ? "selected" : ""}>Within EE</option><option value="miss" ${state.baseline === "miss" ? "selected" : ""}>Outside EE</option></select></div>
    <div><label for="sort">Sort</label><select id="sort"><option value="species_effect" ${state.sort === "species_effect" ? "selected" : ""}>CCI absolute-error effect</option><option value="fixed_error" ${state.sort === "fixed_error" ? "selected" : ""}>Fixed-6S absolute error</option><option value="confidence" ${state.sort === "confidence" ? "selected" : ""}>Species separation</option><option value="truth" ${state.sort === "truth" ? "selected" : ""}>AERONET AOD</option><option value="cloud" ${state.sort === "cloud" ? "selected" : ""}>Cloud fraction</option><option value="site" ${state.sort === "site" ? "selected" : ""}>Site</option></select></div>
  </div>`;
}

function methodMetrics() {
  const methods = ["lut_continental", "sixs_continental", "sixs_cci3"];
  return `<div class="method-metrics">${methods.map((method) => {
    const row = state.data.metrics[method];
    const ci = row.strict_rate_ci95;
    return `<section class="method-block" style="--method:${state.data.method_colors[method]}">
      <h3>${escapeHtml(state.data.method_short_labels[method])}</h3>
      <strong>${row.hits}/${row.expected}</strong><span>within EE, ${pct(row.strict_rate)}</span>
      <dl><dt>95% Wilson CI</dt><dd>${ci ? `${pct(ci[0])} to ${pct(ci[1])}` : "NA"}</dd><dt>Gap to above 87%</dt><dd>${row.additional_hits_to_above_87} hits</dd><dt>RMSE</dt><dd>${fmt(row.rmse)}</dd><dt>MAE</dt><dd>${fmt(row.mae)}</dd><dt>Bias</dt><dd>${signed(row.bias)}</dd><dt>R2</dt><dd>${fmt(row.r2)}</dd><dt>Median runtime</dt><dd>${seconds(row.median_runtime_s)}</dd></dl>
    </section>`;
  }).join("")}</div>`;
}

function comparisonTable() {
  const definitions = [
    ["sixs_continental_vs_lut", "6S fixed vs LUT", "RT backend/control change"],
    ["sixs_cci3_vs_sixs_continental", "CCI-3 vs fixed 6S", "Added variable species within native 6S"],
    ["sixs_cci3_vs_lut", "CCI-3 vs LUT", "End-to-end species arm vs current control"],
  ];
  return `<table class="data-table"><thead><tr><th>Pair</th><th>Interpretation</th><th>Gains</th><th>Losses</th><th>Net hits</th><th>McNemar exact p</th><th>Mean delta |error|</th><th>Bootstrap 95% CI</th></tr></thead><tbody>${definitions.map(([key, label, interpretation]) => {
    const row = state.data.comparisons[key];
    const boot = row.bootstrap_abs_error_delta;
    return `<tr><td><strong>${escapeHtml(label)}</strong></td><td>${escapeHtml(interpretation)}</td><td class="status-hit mono">${row.gains}</td><td class="status-miss mono">${row.losses}</td><td class="mono">${signed(row.net_hits, 0)}</td><td class="mono">${fmt(row.mcnemar_exact_p, 4)}</td><td class="mono">${signed(boot.mean, 4)}</td><td class="mono">${boot.ci95 ? `${signed(boot.ci95[0], 4)} to ${signed(boot.ci95[1], 4)}` : "NA"}</td></tr>`;
  }).join("")}</tbody></table>`;
}

function policyReplayTable() {
  const rows = Object.entries(state.data.final_pass_policy_replays);
  return `<table class="data-table"><thead><tr><th>Final-pass policy</th><th>AERONET-free selector</th><th>Within EE</th><th>Rate</th><th>RMSE</th><th>MAE</th><th>Bias</th></tr></thead><tbody>${rows.map(([, row]) => `<tr><td><strong>${escapeHtml(row.label)}</strong></td><td>${row.operational_without_aeronet ? "Yes" : "No - diagnostic oracle"}</td><td class="mono">${row.hits}/${row.count}</td><td class="mono">${pct(row.rate)}</td><td class="mono">${fmt(row.rmse)}</td><td class="mono">${fmt(row.mae)}</td><td class="mono">${signed(row.bias)}</td></tr>`).join("")}</tbody></table>`;
}

function directSmokeTable() {
  const ids = new Set(state.data.direct_libradtran_smoke.matchup_ids);
  const rows = state.data.cases.filter((item) => ids.has(item.matchup_id));
  const metrics = state.data.metrics.libradtran_continental_smoke;
  const exactMetrics = state.data.metrics.sixs_cci_exact_smoke;
  const comparison = state.data.direct_libradtran_smoke.comparison_vs_lut;
  const exactComparison = state.data.exact_cci_smoke.comparisons?.vs_sixs_continental;
  return `<div class="metric-row smoke-metrics">${metric("Direct within EE", `${metrics.hits}/${metrics.expected}`)}${metric("Direct RMSE", fmt(metrics.rmse))}${metric("Direct vs LUT gains", comparison?.gains ?? "NA", "positive")}${metric("Direct vs LUT losses", comparison?.losses ?? "NA", "negative")}${metric("Exact CCI within EE", `${exactMetrics.hits}/${exactMetrics.expected}`)}${metric("Exact CCI RMSE", fmt(exactMetrics.rmse))}${metric("Exact vs fixed gains", exactComparison?.gains ?? "NA", "positive")}${metric("Exact vs fixed losses", exactComparison?.losses ?? "NA", "negative")}</div><table class="data-table"><thead><tr><th>Site</th><th>AERONET</th><th>LUT</th><th>6S fixed</th><th>CCI-3</th><th>Exact CCI</th><th>Direct libRadTran</th><th>Exact / direct within EE</th></tr></thead><tbody>${rows.map((item) => {
    const direct = item.methods.libradtran_continental_smoke;
    const exact = item.methods.sixs_cci_exact_smoke;
    return `<tr><td><strong>${escapeHtml(item.site)}</strong><span class="subtext">${escapeHtml(item.matchup_id)}</span></td><td class="mono">${fmt(item.truth_aod)}</td><td class="mono">${fmt(item.methods.lut_continental.aod)}</td><td class="mono">${fmt(item.methods.sixs_continental.aod)}</td><td class="mono">${fmt(item.methods.sixs_cci3.aod)}</td><td class="mono">${fmt(exact.aod)}</td><td class="mono">${fmt(direct.aod)}</td><td><span class="${exact.within_ee ? "status-hit" : "status-miss"}">${exact.within_ee ? "Exact: yes" : "Exact: no"}</span><br><span class="${direct.within_ee ? "status-hit" : "status-miss"}">${direct.within_ee ? "Direct: yes" : "Direct: no"}</span></td></tr>`;
  }).join("")}</tbody></table>`;
}

function stratumTable(key) {
  const rows = state.data.stratification[key];
  return `<table class="data-table compact"><thead><tr><th>Stratum</th><th>Cases</th><th>LUT</th><th>6S fixed</th><th>CCI-3</th></tr></thead><tbody>${rows.map((row) => `<tr><td><strong>${escapeHtml(row.label)}</strong>${row.description ? `<span class="subtext">${escapeHtml(row.description)}</span>` : ""}</td><td class="mono">${row.count}</td>${["lut_continental", "sixs_continental", "sixs_cci3"].map((method) => `<td class="mono">${row.methods[method].hits}/${row.count} (${pct(row.methods[method].rate)})</td>`).join("")}</tr>`).join("")}</tbody></table>`;
}

function overviewView() {
  const speciesPair = state.data.comparisons.sixs_cci3_vs_sixs_continental;
  const smoke = state.data.direct_libradtran_smoke;
  return `<main class="page-view">
    <header class="page-heading"><div><h2>Full-cohort paired results</h2><p>${escapeHtml(state.data.cohort_definition)}. Every primary method uses the same 152-case denominator.</p></div><div class="snapshot">Generated ${escapeHtml(state.data.generated_at)}</div></header>
    ${methodMetrics()}
    <div class="metric-row overview-metrics">${metric("CCI-3 gains", speciesPair.gains, "positive")}${metric("CCI-3 losses", speciesPair.losses, "negative")}${metric("CCI-3 net hits", signed(speciesPair.net_hits, 0), speciesPair.net_hits >= 0 ? "positive" : "negative")}${metric("Paired common cases", speciesPair.common_valid)}${metric("Direct libRadTran smoke", `${smoke.case_count} cases`)}${metric("Target", ">87% within EE")}</div>
    <section class="report-section"><div class="section-heading"><h3>Retrieval and paired-effect overview</h3><span>All 152 paired outputs</span></div><figure class="evidence-figure"><img class="zoomable" src="assets/global/performance_overview.png" alt="Full-cohort retrieval and paired effects"><figcaption>Outside-EE points are outlined in red. Negative paired absolute-error deltas indicate improvement.</figcaption></figure></section>
    <section class="report-section"><div class="section-heading"><h3>Paired comparisons</h3><span>Gains and losses are expected-error transitions on identical matchup IDs</span></div>${comparisonTable()}</section>
    <section class="report-section"><div class="section-heading"><h3>Species-selection policy diagnostics</h3><span>Same saved CCI final-pass cubes; oracle rows are diagnostic ceilings and are not operational algorithms</span></div>${policyReplayTable()}<div class="neutral-note">The nearest-mixture and scene-wide rows do not replay pass 1. Any promising policy requires a complete end-to-end rerun before it is compared as an experiment arm.</div></section>
    <section class="report-section"><div class="section-heading"><h3>Exact-CCI and direct-libRadTran smoke controls</h3><span>Six predefined cases; neither smoke score replaces a full-cohort result</span></div><figure class="evidence-figure"><img class="zoomable" src="assets/global/smoke_controls.png" alt="Six-case RT and aerosol-species smoke controls"><figcaption>Retrieval outputs, normalized error, anchor-iteration movement, and end-to-end runtime for the predefined cases.</figcaption></figure><div class="smoke-table-block">${directSmokeTable()}</div></section>
    <section class="report-section"><div class="section-heading"><h3>Stratified and spatial evidence</h3><span>Group H is the prior LUT-baseline within-EE set</span></div><figure class="evidence-figure"><img class="zoomable" src="assets/global/stratified_spatial.png" alt="Stratified and spatial species evidence"><figcaption>Rates are shown by AERONET AOD, cloud fraction, prior failure-evidence group, and matchup location.</figcaption></figure><div class="table-grid"><div><h4>AERONET AOD</h4>${stratumTable("truth")}</div><div><h4>Image cloud fraction</h4>${stratumTable("cloud")}</div></div><div class="wide-table"><h4>Prior baseline failure-evidence groups</h4>${stratumTable("baseline_groups")}</div></section>
    <section class="report-section"><div class="section-heading"><h3>CCI candidate behavior</h3><span>Candidate rank denotes distance from the monthly Aerosol_cci mixture, not one fixed species</span></div><figure class="evidence-figure"><img class="zoomable" src="assets/global/species_selection.png" alt="CCI candidate composition and selection behavior"><figcaption>Selection shares, component mixtures, cost separation, spatial dominance, and paired retrieval effect are kept distinct.</figcaption></figure></section>
    <section class="report-section"><div class="section-heading"><h3>Scope and audit notes</h3><span>No operational method decision is made in this report</span></div><ul class="audit-list">${state.data.audit_notes.map((note) => `<li>${escapeHtml(note)}</li>`).join("")}</ul></section>
  </main>`;
}

function sidebar(rows) {
  return `<aside class="case-sidebar"><div class="filter-panel">${filterControls()}<div class="filter-summary">${rows.length} of ${state.data.cohort_count} cases</div></div><div class="case-list">${rows.map((item, index) => {
    const fixed = item.methods.sixs_continental;
    const cci = item.methods.sixs_cci3;
    return `<button class="case-row ${item.matchup_id === state.selectedId ? "active" : ""}" data-case="${escapeHtml(item.matchup_id)}"><span class="case-rank">${String(index + 1).padStart(3, "0")}</span><span class="case-name"><strong>${escapeHtml(item.site)}</strong><span>AERONET ${fmt(item.truth_aod)} | fixed ${fmt(fixed.aod)} | CCI ${fmt(cci.aod)}</span></span><span class="case-score"><span class="transition-chip ${transitionClass(item.transition_cci_vs_fixed)}">${escapeHtml(transitionLabels[item.transition_cci_vs_fixed])}</span><strong class="${item.cci_minus_fixed_abs_error <= 0 ? "status-hit" : "status-miss"}">${signed(item.cci_minus_fixed_abs_error)} |err|</strong></span></button>`;
  }).join("") || '<div class="empty-state">No cases match these filters.</div>'}</div></aside>`;
}

function methodOutputTable(item) {
  const rows = Object.entries(item.methods).filter(([, row]) => row.status !== "MISSING");
  return `<table class="data-table"><thead><tr><th>Method</th><th>Status</th><th>Pass 1</th><th>Final AOD</th><th>Error</th><th>Error / EE</th><th>Within EE</th><th>Runtime</th></tr></thead><tbody>${rows.map(([method, row]) => `<tr><td><strong>${escapeHtml(state.data.method_short_labels[method])}</strong></td><td>${escapeHtml(row.status)}</td><td class="mono">${fmt(row.pass1_aod)}</td><td class="mono">${fmt(row.aod)}</td><td class="mono">${signed(row.error)}</td><td class="mono">${signed(row.error_ee, 2)}</td><td class="${row.within_ee ? "status-hit" : "status-miss"}">${row.within_ee ? "Yes" : "No"}</td><td class="mono">${seconds(row.runtime_s)}</td></tr>`).join("")}</tbody></table>`;
}

function spatialTab(item) {
  return `<section class="evidence-section"><div class="section-heading"><h3>Exact final-pass spatial evidence</h3><span>60 m solver grid; cross and dashed circle mark the station and 1.5 km window</span></div><figure class="evidence-figure"><img class="zoomable" src="${item.spatial_image}" alt="Spatial aerosol-species evidence for ${escapeHtml(item.site)}" loading="eager"><figcaption><span>TOA, surface prior, support, atmospheric backstop, fixed and candidate AOD maps, selected species, cost separation, and species fractions.</span><a href="${item.spatial_image}" target="_blank" rel="noopener">Open PNG</a></figcaption></figure></section>`;
}

function spectralTab(item) {
  return `<section class="evidence-section"><div class="section-heading"><h3>Spectral, cost, iteration, and mixture evidence</h3><span>${item.cost_grid.aod_nodes} AOD nodes; ${item.cost_grid.curve_support_count} ${escapeHtml(item.cost_grid.curve_scope)} pixels</span></div><figure class="evidence-figure"><img class="zoomable" src="${item.diagnostic_image}" alt="Spectral and cost aerosol-species evidence for ${escapeHtml(item.site)}" loading="eager"><figcaption><span>Cost profiles, mixture composition, selected-pixel share, spectral closure, pass-1/final outputs, and final surface-prior spectra.</span><a href="${item.diagnostic_image}" target="_blank" rel="noopener">Open PNG</a></figcaption></figure></section>`;
}

function mixtureBars(candidate) {
  const colors = { dust: "#bd7b2c", sea_salt: "#3b83bd", fine_strong: "#9b3f52", fine_weak: "#4b8c61" };
  return `<div class="mixture-bar" title="${escapeHtml(candidate.mixture_text)}">${Object.entries(candidate.mixture).map(([name, value]) => `<span style="width:${100 * value}%;background:${colors[name]}" aria-label="${escapeHtml(name)} ${pct(value)}"></span>`).join("")}</div>`;
}

function speciesTab(item) {
  const species = item.species;
  return `<section class="evidence-section"><div class="section-heading"><h3>Per-scene CCI candidates and exact selection replay</h3><span>Reported minus replayed scene mean: ${signed(species.replay_delta, 9)}</span></div><div class="metric-row">${metric("Selected pixels", species.selected_pixel_count)}${metric("Reported CCI AOD", fmt(species.reported_scene_mean_aod, 6))}${metric("Replayed CCI AOD", fmt(species.replayed_scene_mean_aod, 6))}${metric("Median separation", fmt(species.selection_confidence_median, 4))}${metric("P10 separation", fmt(species.selection_confidence_p10, 4))}${metric("CCI - fixed map mean", signed(species.cci_minus_fixed_map_mean, 4))}</div><div class="target-mixture"><strong>Monthly 1-degree CCI climatology target</strong><span>${escapeHtml(species.climatology_target_text)}</span>${mixtureBars({ mixture: species.climatology_target, mixture_text: species.climatology_target_text })}</div><table class="data-table species-table"><thead><tr><th>Candidate</th><th>Physical mixture</th><th>Target L1 distance</th><th>Component fractions</th><th>Selected pixels</th><th>Selected share</th><th>Candidate mean AOD</th><th>Median surface cost</th></tr></thead><tbody>${species.candidates.map((candidate) => `<tr><td><strong>CCI-${candidate.index + 1}</strong></td><td>${escapeHtml(candidate.mixture_text)}</td><td class="mono">${fmt(candidate.climatology_l1_distance, 3)}</td><td>${mixtureBars(candidate)}</td><td class="mono">${candidate.selected_count}</td><td class="mono">${pct(candidate.selected_fraction, 2)}</td><td class="mono">${fmt(candidate.scene_mean_aod, 6)}</td><td class="mono">${fmt(candidate.median_surface_cost, 5)}</td></tr>`).join("")}</tbody></table><div class="neutral-note">Candidate selection uses the minimum pooled observation-only surface cost at each valid pixel. A candidate's AOD solution still includes the common atmospheric backstop during node selection.</div></section>`;
}

function rawTab(item) {
  const fixed = item.methods.sixs_continental;
  const cci = item.methods.sixs_cci3;
  const sourceLinks = Object.entries(item.source_urls).map(([method, url]) => `<a href="${url}" target="_blank" rel="noopener">${escapeHtml(state.data.method_short_labels[method])} JSON</a>`).join("");
  const costLinks = Object.entries(item.cost_urls).map(([label, url]) => `<a href="${url}" target="_blank" rel="noopener">${escapeHtml(label.replaceAll("_", " "))} NPZ</a>`).join("");
  return `<section class="evidence-section"><div class="section-heading"><h3>Raw outputs, solver diagnostics, and source artifacts</h3><span>Values below are from the saved operational result JSON and exact cost cubes</span></div>${methodOutputTable(item)}<div class="data-grid">
    <div class="data-panel"><h3>Matchup and evaluation</h3><table class="kv-table">${kvRows([["Matchup ID", item.matchup_id], ["Site", item.site], ["Longitude / latitude", `${fmt(item.lon, 6)} / ${fmt(item.lat, 6)}`], ["AERONET AOD", fmt(item.truth_aod, 8)], ["Expected-error threshold", fmt(item.ee_threshold, 8)], ["Cloud fraction", pct(item.cloud_fraction, 2)], ["Invalid fraction", pct(item.invalid_fraction, 2)], ["Mask fallback", item.mask_fallback ? "Recorded zero-support fallback" : "Normal cloud/water mask"], ["Prior baseline group", `${item.baseline_failure_group}. ${item.baseline_failure_label}`], ["CCI vs fixed transition", transitionLabels[item.transition_cci_vs_fixed]], ["CCI - fixed absolute error", signed(item.cci_minus_fixed_abs_error, 8)]])}</table></div>
    <div class="data-panel"><h3>Cost grid and atmospheric backstop</h3><table class="kv-table">${kvRows([["Grid shape", item.cost_grid.shape.join(" x ")], ["AOD axis", `${item.cost_grid.aod_nodes} nodes; ${fmt(item.cost_grid.aod_min, 3)} to ${fmt(item.cost_grid.aod_max, 2)}`], ["Pool window / minimum", `${item.cost_grid.pool_window} / ${item.cost_grid.pool_min_count}`], ["Curve scope", item.cost_grid.curve_scope], ["Curve support pixels", item.cost_grid.curve_support_count], ["Station-window valid pixels", item.cost_grid.station_window_valid_count], ["Backstop median AOD", fmt(item.cost_grid.solver_backstop_median, 8)], ["Backstop median sigma", fmt(item.cost_grid.solver_backstop_unc_median, 8)]])}</table></div>
    <div class="data-panel"><h3>Fixed continental 6S solver</h3><table class="kv-table">${kvRows([["Final cost / band", fmt(fixed.solver.surface_static_cost_per_band, 8)], ["Cost-curve minimum AOD", fmt(fixed.solver.surface_cost_curve_min_aot, 8)], ["B02 minimum AOD", fmt(fixed.solver.surface_band_B02_argmin_aot, 8)], ["B03 minimum AOD", fmt(fixed.solver.surface_band_B03_argmin_aot, 8)], ["B04 minimum AOD", fmt(fixed.solver.surface_band_B04_argmin_aot, 8)], ["Band-minimum spread", fmt(fixed.solver.surface_band_argmin_spread, 8)], ["Solved-pixel fraction", pct(fixed.solver.surface_solved_pixel_fraction, 3)], ["Station / window AOD", `${fmt(fixed.station_aod, 6)} / ${fmt(fixed.window_median_aod, 6)}`]])}</table></div>
    <div class="data-panel"><h3>CCI-selected 6S solver</h3><table class="kv-table">${kvRows([["RT branch", cci.solver.surface_rt_branch], ["Diagnostic candidate", cci.solver.surface_rt_branch_diagnostic_candidate], ["Final cost / band", fmt(cci.solver.surface_static_cost_per_band, 8)], ["Cost-curve minimum AOD", fmt(cci.solver.surface_cost_curve_min_aot, 8)], ["B02 minimum AOD", fmt(cci.solver.surface_band_B02_argmin_aot, 8)], ["B03 minimum AOD", fmt(cci.solver.surface_band_B03_argmin_aot, 8)], ["B04 minimum AOD", fmt(cci.solver.surface_band_B04_argmin_aot, 8)], ["Band-minimum spread", fmt(cci.solver.surface_band_argmin_spread, 8)], ["Solved-pixel fraction", pct(cci.solver.surface_solved_pixel_fraction, 3)], ["Station / window AOD", `${fmt(cci.station_aod, 6)} / ${fmt(cci.window_median_aod, 6)}`]])}</table></div>
  </div><div class="source-links">${sourceLinks}${costLinks}</div></section>`;
}

function detailBody(item) {
  if (state.detailTab === "spectral") return spectralTab(item);
  if (state.detailTab === "species") return speciesTab(item);
  if (state.detailTab === "raw") return rawTab(item);
  return spatialTab(item);
}

function caseMain(item, rows) {
  const index = Math.max(0, rows.findIndex((row) => row.matchup_id === item.matchup_id));
  const previous = rows[(index - 1 + rows.length) % rows.length];
  const next = rows[(index + 1) % rows.length];
  const baseline = item.methods.lut_continental;
  const fixed = item.methods.sixs_continental;
  const cci = item.methods.sixs_cci3;
  return `<main class="case-main"><header class="case-header"><div class="case-header-top"><div><h2>${escapeHtml(item.site)}</h2><div class="matchup">${escapeHtml(item.matchup_id)}</div></div><div class="case-actions"><button class="icon-button" data-case="${escapeHtml(previous?.matchup_id)}" title="Previous case" aria-label="Previous case">&larr;</button><button class="icon-button" data-case="${escapeHtml(next?.matchup_id)}" title="Next case" aria-label="Next case">&rarr;</button></div></div><div class="metric-row">${metric("AERONET", fmt(item.truth_aod))}${metric("LUT", fmt(baseline.aod), baseline.within_ee ? "positive" : "negative")}${metric("6S fixed", fmt(fixed.aod), fixed.within_ee ? "positive" : "negative")}${metric("CCI-3", fmt(cci.aod), cci.within_ee ? "positive" : "negative")}${metric("CCI - fixed |error|", signed(item.cci_minus_fixed_abs_error), item.cci_minus_fixed_abs_error <= 0 ? "positive" : "negative")}${metric("Cloud", pct(item.cloud_fraction))}</div><div class="signature-row"><span class="tag ${transitionClass(item.transition_cci_vs_fixed)}"><strong>CCI vs fixed</strong>${escapeHtml(transitionLabels[item.transition_cci_vs_fixed])}</span><span class="tag"><strong>Baseline group</strong>${escapeHtml(item.baseline_failure_group)}</span>${item.mask_fallback ? '<span class="tag"><strong>Mask</strong>Zero-support fallback</span>' : ""}<span class="tag"><strong>Backstop</strong>${fmt(item.cost_grid.solver_backstop_median)} +/- ${fmt(item.cost_grid.solver_backstop_unc_median)}</span><span class="tag"><strong>Candidate separation</strong>${fmt(item.species.selection_confidence_median, 4)}</span><span class="tag"><strong>Valid pixels</strong>${item.species.selected_pixel_count}</span></div></header><nav class="detail-tabs" aria-label="Case evidence views"><button data-tab="spatial" class="${state.detailTab === "spatial" ? "active" : ""}">Spatial</button><button data-tab="spectral" class="${state.detailTab === "spectral" ? "active" : ""}">Spectral and cost</button><button data-tab="species" class="${state.detailTab === "species" ? "active" : ""}">Species candidates</button><button data-tab="raw" class="${state.detailTab === "raw" ? "active" : ""}">Raw diagnostics</button></nav><div class="detail-body">${detailBody(item)}</div></main>`;
}

function casesView() {
  const rows = filteredCases();
  let item = selectedCase();
  if (rows.length && !rows.some((row) => row.matchup_id === item.matchup_id)) {
    item = rows[0];
    state.selectedId = item.matchup_id;
  }
  return `<div class="explorer-layout">${sidebar(rows)}${rows.length ? caseMain(item, rows) : '<main class="case-main"><div class="empty-state">No cases match these filters.</div></main>'}</div>`;
}

function galleryView() {
  const rows = filteredCases();
  const key = state.galleryMode === "diagnostic" ? "diagnostic_image" : "spatial_image";
  return `<main class="page-view"><header class="page-heading"><div><h2>All-case evidence gallery</h2><p>Spatial and spectral/cost evidence across every matching case.</p></div></header><div class="gallery-toolbar"><div class="segmented"><button data-gallery-mode="spatial" class="${state.galleryMode === "spatial" ? "active" : ""}">Spatial</button><button data-gallery-mode="diagnostic" class="${state.galleryMode === "diagnostic" ? "active" : ""}">Spectral and cost</button></div>${filterControls()}<span class="gallery-count">${rows.length} / ${state.data.cohort_count}</span></div><div class="gallery-grid">${rows.map((item) => `<article class="gallery-item"><button data-open-case="${escapeHtml(item.matchup_id)}"><img src="${item[key]}" alt="Evidence thumbnail for ${escapeHtml(item.site)}" loading="lazy"><div class="gallery-meta"><strong>${escapeHtml(item.site)}</strong><span>AERONET ${fmt(item.truth_aod)} | fixed ${fmt(item.methods.sixs_continental.aod)} | CCI ${fmt(item.methods.sixs_cci3.aod)}</span><span class="transition-chip ${transitionClass(item.transition_cci_vs_fixed)}">${escapeHtml(transitionLabels[item.transition_cci_vs_fixed])}</span></div></button></article>`).join("") || '<div class="empty-state">No cases match these filters.</div>'}</div></main>`;
}

function methodsView() {
  const recipe = state.data.fixed_recipe;
  const source = state.data.source_urls;
  return `<main class="page-view"><header class="page-heading"><div><h2>Method contract and source snapshot</h2><p>The controlled variables and changing variables are stated explicitly.</p></div></header><section class="report-section"><div class="section-heading"><h3>Experiment arms</h3><span>Only the primary three arms use the complete 152-case cohort</span></div><table class="data-table"><thead><tr><th>Method</th><th>Radiative transfer</th><th>Aerosol species treatment</th><th>Scope</th></tr></thead><tbody>${state.data.method_contract.map((row) => `<tr><td><strong>${escapeHtml(state.data.method_labels[row.method])}</strong></td><td>${escapeHtml(row.rt)}</td><td>${escapeHtml(row.species)}</td><td>${escapeHtml(row.scope)}</td></tr>`).join("")}</tbody></table></section><section class="report-section"><div class="section-heading"><h3>Fixed retrieval recipe</h3><span>These settings are held constant across the primary arms</span></div><table class="kv-table recipe-table">${kvRows(Object.entries(recipe).map(([key, value]) => [key.replaceAll("_", " "), value]))}</table></section><section class="report-section"><div class="section-heading"><h3>Raw source directories</h3><span>Exact JSON and NPZ artifacts are directly linked</span></div><div class="source-links">${Object.entries(source).map(([label, url]) => `<a href="${url}" target="_blank" rel="noopener">${escapeHtml(label.replaceAll("_", " "))}</a>`).join("")}</div></section><section class="report-section"><div class="section-heading"><h3>Direct libRadTran smoke scope</h3><span>This is an implementation control, not a full-cohort score</span></div><table class="kv-table recipe-table">${kvRows([["Completed cases", state.data.direct_libradtran_smoke.case_count], ["Matchup IDs", state.data.direct_libradtran_smoke.matchup_ids.join(", ")], ["Purpose", "Check direct uvspec behavior against the precomputed continental Zarr control on selected cases"]])}</table></section></main>`;
}

function topbar() {
  return `<header class="topbar"><div class="title-lockup"><h1>${escapeHtml(state.data.title)}</h1><p>${state.data.cohort_count} low-cloud matchups | one surface-prior type | fixed 6S versus CCI aerosol mixtures</p></div><nav class="view-nav"><button data-view="overview" class="${state.view === "overview" ? "active" : ""}">Overview</button><button data-view="cases" class="${state.view === "cases" ? "active" : ""}">Case explorer</button><button data-view="gallery" class="${state.view === "gallery" ? "active" : ""}">Evidence gallery</button><button data-view="methods" class="${state.view === "methods" ? "active" : ""}">Methods and sources</button></nav></header>`;
}

function render() {
  let body = overviewView();
  if (state.view === "cases") body = casesView();
  if (state.view === "gallery") body = galleryView();
  if (state.view === "methods") body = methodsView();
  app.innerHTML = `${topbar()}${body}<dialog class="image-dialog"><div class="dialog-bar"><strong></strong><button class="icon-button" data-close-dialog title="Close" aria-label="Close">&times;</button></div><img alt="Expanded diagnostic figure"></dialog>`;
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
  document.querySelectorAll("[data-tab]").forEach((button) => button.addEventListener("click", () => { state.detailTab = button.dataset.tab; render(); }));
  document.querySelectorAll("[data-gallery-mode]").forEach((button) => button.addEventListener("click", () => { state.galleryMode = button.dataset.galleryMode; render(); }));
  document.querySelectorAll(".zoomable").forEach((image) => image.addEventListener("click", () => {
    const dialog = document.querySelector(".image-dialog");
    dialog.querySelector("strong").textContent = image.alt;
    dialog.querySelector("img").src = image.src;
    dialog.querySelector("img").alt = image.alt;
    dialog.showModal();
  }));
  document.querySelector("[data-close-dialog]")?.addEventListener("click", () => document.querySelector(".image-dialog")?.close());
  document.getElementById("search")?.addEventListener("input", (event) => updateFilter("search", event.target.value));
  document.getElementById("transition")?.addEventListener("change", (event) => updateFilter("transition", event.target.value));
  document.getElementById("group")?.addEventListener("change", (event) => updateFilter("group", event.target.value));
  document.getElementById("baseline")?.addEventListener("change", (event) => updateFilter("baseline", event.target.value));
  document.getElementById("sort")?.addEventListener("change", (event) => updateFilter("sort", event.target.value));
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

fetch("data/report.json")
  .then((response) => {
    if (!response.ok) throw new Error(`HTTP ${response.status}`);
    return response.json();
  })
  .then((data) => {
    state.data = data;
    const match = location.hash.match(/^#case=(.+)$/);
    const requested = match ? decodeURIComponent(match[1]) : null;
    if (data.cases.some((item) => item.matchup_id === requested)) {
      state.selectedId = requested;
      state.view = "cases";
    } else {
      state.selectedId = data.cases[0]?.matchup_id;
    }
    render();
  })
  .catch((error) => {
    app.innerHTML = `<div class="empty-state">Unable to load the experiment snapshot: ${escapeHtml(error.message)}</div>`;
  });
