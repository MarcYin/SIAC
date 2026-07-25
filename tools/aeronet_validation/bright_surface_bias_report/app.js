(() => {
  "use strict";

  const state = {
    report: null,
    view: "overview",
    query: "",
    brightnessGroup: "all",
    outcome: "all",
    selectedId: null,
  };

  const $ = (selector, root = document) => root.querySelector(selector);
  const escapeHtml = (value) => String(value ?? "").replace(/[&<>'"]/g, (char) => ({
    "&": "&amp;", "<": "&lt;", ">": "&gt;", "'": "&#39;", '"': "&quot;",
  }[char]));
  const finite = (value) => value !== null && value !== "" && Number.isFinite(Number(value));
  const fmt = (value, digits = 3) => finite(value) ? Number(value).toFixed(digits) : "NA";
  const signed = (value, digits = 3) => finite(value) ? `${Number(value) >= 0 ? "+" : ""}${Number(value).toFixed(digits)}` : "NA";
  const pct = (value, digits = 1) => finite(value) ? `${Number(value).toFixed(digits)}%` : "NA";
  const fraction = (value, digits = 1) => finite(value) ? pct(100 * Number(value), digits) : "NA";
  const outcomeClass = (value) => value === "hit" ? "hit" : value === "under" ? "under" : "over";
  const metric = (label, value, className = "") => `<div class="metric"><div class="metric-label">${escapeHtml(label)}</div><div class="metric-value ${className}">${escapeHtml(value)}</div></div>`;

  function filteredCases() {
    const query = state.query.trim().toLowerCase();
    return state.report.cases.filter((item) => {
      const matchesQuery = !query || `${item.site} ${item.matchup_id}`.toLowerCase().includes(query);
      return matchesQuery
        && (state.brightnessGroup === "all" || item.brightness_group === state.brightnessGroup)
        && (state.outcome === "all" || item.outcome === state.outcome);
    });
  }

  function appHeader() {
    return `<header class="topbar"><div class="title-block"><h1>${escapeHtml(state.report.title)}</h1><p>${escapeHtml(state.report.subtitle)}</p></div><nav class="nav" aria-label="Report views"><button data-view="overview" class="${state.view === "overview" ? "active" : ""}">Evidence overview</button><button data-view="cases" class="${state.view === "cases" ? "active" : ""}">All 36 cases</button></nav></header>`;
  }

  function overview() {
    const summary = state.report.summary;
    const current = summary.current;
    const groups = summary.brightness_groups;
    const cutpoints = summary.brightness_cutpoints;
    const corr = summary.correlations_with_current_error;
    const extraction = summary.extraction_modes;
    return `<main class="page"><header class="page-heading"><h2>Evidence overview</h2><p>The bright-surface signal is real, but it does not identify a single correction that can be promoted without losing performance elsewhere.</p></header>
      <section class="metric-grid">${metric("Current within EE", `${current.within_ee}/${current.n} (${pct(current.within_ee_percent)})`, "hit")}${metric("Low brightness", `${groups.low.within_ee}/${groups.low.n} (${pct(groups.low.within_ee_percent)})`, "hit")}${metric("Middle brightness", `${groups.mid.within_ee}/${groups.mid.n} (${pct(groups.mid.within_ee_percent)})`)}${metric("High brightness", `${groups.high.within_ee}/${groups.high.n} (${pct(groups.high.within_ee_percent)})`, "under")}${metric("High-group bias", signed(groups.high.bias), "under")}${metric("High-group MAE", fmt(groups.high.mae))}</section>
      <section class="overview-grid"><section class="panel"><h3>Brightness, error, and information</h3><figure class="figure"><img class="zoomable" src="${state.report.figures.brightness_error}" alt="Current AOD error and retrieval-information summaries by operating-prior brightness"><figcaption>Each point is one Sentinel-2/AERONET matchup; labels mark only outside-EE cases. The left y-axis is current retrieved AOD minus AERONET, not MAIAC AOD. Bright cases have weaker local AOD information and more disagreement between band-level AOD minima.</figcaption></figure></section>
      <section class="panel"><h3>Exact pair evidence</h3><figure class="figure"><img class="zoomable" src="${state.report.figures.pair_bias}" alt="Exact same-day L2A to current RT pair bias by brightness"><figcaption>Each point is one matchup's historical, same-day L2A/current-RT pair archive. The primary slope uses the mean of L2A and current-RT brightness, avoiding the direct shared-axis coupling of difference versus L2A.</figcaption></figure></section>
      <section class="panel"><h3>Negative correction screen</h3><figure class="figure"><img class="zoomable" src="${state.report.figures.screen}" alt="Performance of the raw L2A pair post-prior correction screen"><figcaption>The raw-pair affine mapping is held-site and AERONET-free, but it moves most bright under-retrievals in the wrong direction. It is diagnostic only, not promoted.</figcaption></figure></section>
      <section class="panel"><h3>Observed findings</h3><ul>${summary.findings.map((item) => `<li>${escapeHtml(item)}</li>`).join("")}</ul><h3 style="margin-top:20px">Interpretation limits</h3><ul>${summary.limitations.map((item) => `<li>${escapeHtml(item)}</li>`).join("")}</ul></section>
      <section class="panel"><h3>Extraction check</h3><p>Station values are missing where the actual station pixel is invalid. Neither available local extraction improves the saved scene-mean score.</p><div class="table-wrap"><table><thead><tr><th>Saved extraction</th><th>Valid cases</th><th>Within EE</th><th>Bias</th><th>MAE</th></tr></thead><tbody><tr><td>Canonical scene mean</td><td>${extraction.scene_mean.n}</td><td>${extraction.scene_mean.within_ee}/${extraction.scene_mean.n}</td><td class="mono">${signed(extraction.scene_mean.bias)}</td><td class="mono">${fmt(extraction.scene_mean.mae)}</td></tr><tr><td>Station pixel</td><td>${extraction.station.n}</td><td>${extraction.station.within_ee}/${extraction.station.n}</td><td class="mono">${signed(extraction.station.bias)}</td><td class="mono">${fmt(extraction.station.mae)}</td></tr><tr><td>1.5 km median</td><td>${extraction.station_window_median.n}</td><td>${extraction.station_window_median.within_ee}/${extraction.station_window_median.n}</td><td class="mono">${signed(extraction.station_window_median.bias)}</td><td class="mono">${fmt(extraction.station_window_median.mae)}</td></tr></tbody></table></div></section>
      <section class="panel full-width"><h3>Brightness strata and diagnostic associations</h3><div class="table-wrap"><table><thead><tr><th>Stratum</th><th>Operating prior luminance</th><th>Within EE</th><th>Bias</th><th>MAE</th><th>Brightness/error correlation</th><th>Band-minimum spread/error correlation</th><th>AOD-information/error correlation</th></tr></thead><tbody><tr><td>Low</td><td>&lt;= ${fmt(cutpoints.low_to_mid, 4)}</td><td>${groups.low.within_ee}/${groups.low.n}</td><td class="mono">${signed(groups.low.bias)}</td><td class="mono">${fmt(groups.low.mae)}</td><td class="mono">${fmt(corr.brightness)}</td><td class="mono">${fmt(corr.band_minimum_spread)}</td><td class="mono">${fmt(corr.aod_information)}</td></tr><tr><td>Middle</td><td>${fmt(cutpoints.low_to_mid, 4)} to ${fmt(cutpoints.mid_to_high, 4)}</td><td>${groups.mid.within_ee}/${groups.mid.n}</td><td class="mono">${signed(groups.mid.bias)}</td><td class="mono">${fmt(groups.mid.mae)}</td><td colspan="3" rowspan="2">Correlations are cohort-wide, shown here once to avoid suggesting they are causal within a brightness group.</td></tr><tr><td>High</td><td>&gt; ${fmt(cutpoints.mid_to_high, 4)}</td><td>${groups.high.within_ee}/${groups.high.n}</td><td class="mono">${signed(groups.high.bias)}</td><td class="mono">${fmt(groups.high.mae)}</td></tr></tbody></table></div></section>
      <section class="panel full-width"><h3>Scope and reproducibility</h3><div class="note">${escapeHtml(state.report.scope.target)} The pair model is historical through ${escapeHtml(state.report.scope.pair_training_cutoff)} and excludes each held site during fitting. The detailed case page retains all intermediate values and raw source links.</div><div class="source-links"><a href="${state.report.sources.cases_csv}" target="_blank" rel="noopener">Download case metrics CSV</a><a href="${state.report.sources.report_json}" target="_blank" rel="noopener">Report JSON</a><a href="${state.report.sources.current_results}" target="_blank" rel="noopener">Current results</a><a href="${state.report.sources.exact_pairs}" target="_blank" rel="noopener">Exact pair archive</a><a href="${state.report.sources.previous_harmonization_report}" target="_blank" rel="noopener">Previous harmonization report</a></div></section></section></main>`;
  }

  function caseRow(item) {
    const selected = item.matchup_id === state.selectedId;
    return `<tr data-case="${escapeHtml(item.matchup_id)}" class="${selected ? "selected" : ""}"><td class="site">${escapeHtml(item.site)}</td><td>${escapeHtml(item.brightness_group)}</td><td class="mono">${fmt(item.brightness, 4)}</td><td class="mono">${fmt(item.truth)}</td><td class="mono">${fmt(item.retrieved)}</td><td class="mono ${outcomeClass(item.outcome)}">${signed(item.error)}</td><td class="outcome ${outcomeClass(item.outcome)}">${escapeHtml(item.outcome)}</td><td class="mono">${fmt(item.aod_information)}</td><td class="mono">${fmt(item.band_minimum_spread)}</td><td class="mono">${fmt(item.pair.l2a_minus_current_rt_vs_mean_slope)}</td></tr>`;
  }

  function caseDetail(item) {
    if (!item) return "";
    const pair = item.pair || {};
    const model = item.pair_model || {};
    const coeffs = model.coefficients_target_minus_l2a || {};
    const context = item.atmospheric_context || {};
    const savedMins = item.saved_solver_band_minima || {};
    const savedResiduals = item.saved_solver_band_final_residual || {};
    const savedPrior = item.saved_prior_boa || {};
    const savedUncertainty = item.saved_prior_uncertainty || {};
    const decisions = item.prior_conflict_decision || "NA";
    return `<section class="case-detail"><div class="case-title"><div><h3>${escapeHtml(item.site)} <span class="outcome ${outcomeClass(item.outcome)}">${escapeHtml(item.outcome)}</span></h3><p>${escapeHtml(item.matchup_id)}</p><div class="legend"><span class="tag"><strong>Brightness:</strong> ${escapeHtml(item.brightness_group)}</span><span class="tag"><strong>Current conflict rule:</strong> ${escapeHtml(decisions)}</span><span class="tag"><strong>Pair samples:</strong> ${escapeHtml(String(pair.samples ?? "NA"))}</span></div></div><button data-close-detail title="Close selected case" aria-label="Close selected case">×</button></div>
      <div class="metric-grid case-metrics">${metric("AERONET", fmt(item.truth))}${metric("Current AOD", fmt(item.retrieved), item.outcome === "hit" ? "hit" : outcomeClass(item.outcome))}${metric("Error / EE", `${fmt(Math.abs(item.error) / item.ee, 2)}x`, outcomeClass(item.outcome))}${metric("Surface-only minimum", fmt(item.local_surface_min))}${metric("MAIAC", fmt(item.maiac_aod))}${metric("Valid support", fraction(item.valid_support_fraction))}</div>
      <div class="detail-grid"><figure class="figure"><img class="zoomable" src="${escapeHtml(item.figure)}" alt="Detailed spatial, spectral, cost, and pair diagnostics for ${escapeHtml(item.site)}"><figcaption><span>TOA and prior RGB, station-window truth-node residual map, joint and per-band cost curves, signed residual curves, and midpoint difference diagnostics for the exact pairs and target BOA.</span> <a href="${escapeHtml(item.figure)}" target="_blank" rel="noopener">Open PNG</a></figcaption></figure>
        <div class="detail-side"><section class="panel"><h3>Surface and AOD evidence</h3><table class="key-value"><tbody><tr><th>Prior brightness median / p90</th><td>${fmt(item.brightness, 5)} / ${fmt(item.brightness_p90, 5)}</td></tr><tr><th>BOA minus prior at truth</th><td>${signed(item.truth_common_residual, 5)}</td></tr><tr><th>Common residual z</th><td>${signed(item.truth_common_residual_z, 3)}</td></tr><tr><th>Spectral residual</th><td>${fmt(item.truth_spectral_residual, 5)}</td></tr><tr><th>Spatial residual SD</th><td>${fmt(item.truth_spatial_residual_std, 5)}</td></tr><tr><th>AOD information / resolution</th><td>${fmt(item.aod_information)} / ${fmt(item.aod_resolution_proxy)}</td></tr><tr><th>Saved B02/B03/B04 minima</th><td>${fmt(savedMins.B02)} / ${fmt(savedMins.B03)} / ${fmt(savedMins.B04)}</td></tr><tr><th>Saved B02/B03/B04 final residual</th><td>${signed(savedResiduals.B02, 5)} / ${signed(savedResiduals.B03, 5)} / ${signed(savedResiduals.B04, 5)}</td></tr><tr><th>Saved B02/B03/B04 prior</th><td>${fmt(savedPrior.B02, 5)} / ${fmt(savedPrior.B03, 5)} / ${fmt(savedPrior.B04, 5)}</td></tr><tr><th>Saved B02/B03/B04 uncertainty</th><td>${fmt(savedUncertainty.B02, 5)} / ${fmt(savedUncertainty.B03, 5)} / ${fmt(savedUncertainty.B04, 5)}</td></tr><tr><th>Truth cost penalty</th><td>${fmt(item.surface_truth_penalty)}</td></tr></tbody></table></section>
        <section class="panel"><h3>Exact same-day L2A/current-RT pair</h3><table class="key-value"><tbody><tr><th>Difference versus midpoint slope</th><td>${signed(pair.l2a_minus_current_rt_vs_mean_slope, 5)}</td></tr><tr><th>Midpoint Q4 minus Q1</th><td>${signed(pair.l2a_minus_current_rt_q4_minus_q1, 5)}</td></tr><tr><th>Midpoint correlation</th><td>${signed(pair.l2a_minus_current_rt_vs_mean_correlation, 3)}</td></tr><tr><th>Pair median</th><td>${signed(pair.l2a_minus_current_rt_median, 5)}</td></tr><tr><th>Difference versus L2A slope</th><td>${signed(pair.l2a_minus_current_rt_slope, 5)}</td></tr><tr><th>Target residual versus midpoint slope</th><td>${signed(item.target_residual_vs_midpoint_slope, 5)}</td></tr><tr><th>Target residual versus midpoint correlation</th><td>${signed(item.target_residual_vs_midpoint_correlation, 3)}</td></tr><tr><th>Target residual versus prior slope</th><td>${signed(item.target_residual_vs_brightness_slope, 5)}</td></tr></tbody></table><div class="note">The L2A-axis and prior-axis slopes share a term with their differences. Use the midpoint fields and plots for interpretation.</div></section>
        <section class="panel"><h3>Support and atmospheric context</h3><table class="key-value"><tbody><tr><th>Cloud / invalid fraction</th><td>${fraction(item.cloud_fraction)} / ${fraction(item.invalid_fraction)}</td></tr><tr><th>Valid AOT / valid support</th><td>${fraction(item.valid_aot_fraction)} / ${fraction(item.valid_support_fraction)}</td></tr><tr><th>Scene / station / 1.5 km AOD</th><td>${fmt(item.retrieved)} / ${fmt(item.retrieved_station)} / ${fmt(item.retrieved_winmed)}</td></tr><tr><th>MAIAC AOD / sigma</th><td>${fmt(item.maiac_aod)} / ${fmt(item.maiac_uncertainty)}</td></tr><tr><th>Prior conflict z</th><td>${fmt(item.prior_conflict_z)}</td></tr><tr><th>TCWV / elevation</th><td>${fmt(context.tcwv_median)} / ${fmt(context.elevation_median)}</td></tr><tr><th>Total ozone</th><td>${fmt(context.tco3_median)}</td></tr><tr><th>AOD quality score</th><td>${fmt(item.aod_quality_score)}</td></tr></tbody></table></section>
        <section class="panel"><h3>Held-site affine diagnostic screen</h3><table class="key-value"><tbody><tr><th>Historical cutoff</th><td>${escapeHtml(model.training_cutoff || "NA")}</td></tr><tr><th>Raw-pair screen AOD</th><td>${fmt(item.screen_mapped_pair_calibrated)}</td></tr><tr><th>Raw-pair product screen AOD</th><td>${fmt(item.screen_mapped_pair_product)}</td></tr><tr><th>Independent-pixel baseline</th><td>${fmt(item.screen_baseline_calibrated)}</td></tr><tr><th>B02 target-L2A (a,b)</th><td>${fmt(coeffs.B02?.intercept, 5)}, ${fmt(coeffs.B02?.slope, 5)}</td></tr><tr><th>B03 target-L2A (a,b)</th><td>${fmt(coeffs.B03?.intercept, 5)}, ${fmt(coeffs.B03?.slope, 5)}</td></tr><tr><th>B04 target-L2A (a,b)</th><td>${fmt(coeffs.B04?.intercept, 5)}, ${fmt(coeffs.B04?.slope, 5)}</td></tr></tbody></table><div class="note">The screen intentionally does not reproduce spatial pooling. Its role is to establish correction direction, not to select an operational result.</div></section>
        <div class="source-links"><a href="${escapeHtml(item.urls.current_result)}" target="_blank" rel="noopener">Current result JSON</a><a href="${escapeHtml(item.urls.fresh_result)}" target="_blank" rel="noopener">Fresh baseline JSON</a><a href="${escapeHtml(item.urls.cost_cube)}" target="_blank" rel="noopener">Cost cube NPZ</a><a href="${escapeHtml(item.urls.exact_pair)}" target="_blank" rel="noopener">Exact pair NPZ</a>${item.urls.direct_harmonized_result ? `<a href="${escapeHtml(item.urls.direct_harmonized_result)}" target="_blank" rel="noopener">Raw harmonized JSON</a>` : ""}</div></div></div></section>`;
  }

  function cases() {
    const rows = filteredCases();
    const selected = state.report.cases.find((item) => item.matchup_id === state.selectedId);
    return `<main class="page"><header class="explorer-header"><div><h2>All-case evidence explorer</h2><p>Select any case to inspect the full spatial, spectral, cost, and upstream-pair record. No automatic causal label is assigned.</p></div></header><section class="toolbar"><label><span>Search site or matchup</span><input id="case-search" type="search" value="${escapeHtml(state.query)}" placeholder="Search"></label><label><span>Brightness stratum</span><select id="brightness-filter"><option value="all">All</option><option value="low" ${state.brightnessGroup === "low" ? "selected" : ""}>Low</option><option value="mid" ${state.brightnessGroup === "mid" ? "selected" : ""}>Middle</option><option value="high" ${state.brightnessGroup === "high" ? "selected" : ""}>High</option></select></label><label><span>Current result</span><select id="outcome-filter"><option value="all">All</option><option value="hit" ${state.outcome === "hit" ? "selected" : ""}>Within EE</option><option value="under" ${state.outcome === "under" ? "selected" : ""}>Under-retrieval</option><option value="over" ${state.outcome === "over" ? "selected" : ""}>Over-retrieval</option></select></label><div class="case-count">${rows.length} / ${state.report.cases.length} cases</div></section><div class="table-wrap"><table class="case-table"><thead><tr><th>Site</th><th>Stratum</th><th>Prior brightness</th><th>AERONET</th><th>Current</th><th>Error</th><th>Outcome</th><th>AOD information</th><th>Band min spread</th><th>Pair midpoint slope</th></tr></thead><tbody>${rows.map(caseRow).join("") || `<tr><td colspan="10">No cases match these filters.</td></tr>`}</tbody></table></div>${caseDetail(selected)}</main>`;
  }

  function render() {
    document.getElementById("app").innerHTML = `<div class="app-shell">${appHeader()}${state.view === "overview" ? overview() : cases()}</div>`;
    bind();
  }

  function showImage(src, alt) {
    const dialog = document.createElement("dialog");
    dialog.className = "image-dialog";
    dialog.innerHTML = `<div class="dialog-bar"><strong>${escapeHtml(alt)}</strong><button aria-label="Close enlarged image">×</button></div><img src="${escapeHtml(src)}" alt="${escapeHtml(alt)}">`;
    $("button", dialog).addEventListener("click", () => dialog.close());
    dialog.addEventListener("close", () => dialog.remove());
    document.body.append(dialog);
    dialog.showModal();
  }

  function bind() {
    document.querySelectorAll("[data-view]").forEach((button) => button.addEventListener("click", () => { state.view = button.dataset.view; render(); }));
    document.querySelectorAll(".zoomable").forEach((image) => image.addEventListener("click", () => showImage(image.src, image.alt)));
    const search = $("#case-search");
    if (search) search.addEventListener("input", (event) => { state.query = event.target.value; render(); });
    const brightness = $("#brightness-filter");
    if (brightness) brightness.addEventListener("change", (event) => { state.brightnessGroup = event.target.value; render(); });
    const outcome = $("#outcome-filter");
    if (outcome) outcome.addEventListener("change", (event) => { state.outcome = event.target.value; render(); });
    document.querySelectorAll("[data-case]").forEach((row) => row.addEventListener("click", () => { state.selectedId = row.dataset.case; render(); document.querySelector(".case-detail")?.scrollIntoView({ behavior: "smooth", block: "start" }); }));
    $("[data-close-detail]")?.addEventListener("click", () => { state.selectedId = null; render(); });
  }

  fetch("data/report.json")
    .then((response) => {
      if (!response.ok) throw new Error(`Report data failed to load (${response.status})`);
      return response.json();
    })
    .then((report) => { state.report = report; render(); })
    .catch((error) => { document.getElementById("app").innerHTML = `<div class="error-state">${escapeHtml(error.message)}</div>`; });
})();
