const app = document.getElementById("app");

const state = {
  data: null,
  view: "overview",
  selectedId: null,
  detailTab: "predictions",
  search: "",
  transition: "all",
  hitStatus: "all",
  fold: "all",
  seen: "all",
  diagnostics: "all",
  sort: "final_error",
  galleryMode: "spatial",
  galleryTransition: "all",
  experimentCategory: "all",
};

const colours = {
  baseline: "#68757c",
  raw: "#b46b00",
  final: "#1769aa",
  target: "#222c31",
  gain: "#1769aa",
  loss: "#b24a3b",
  remaining_miss: "#8a5a13",
  retained_hit: "#68757c",
};

const transitionLabels = {
  gain: "Gain",
  loss: "Loss",
  remaining_miss: "Remaining miss",
  retained_hit: "Retained hit",
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

const percent = (value, digits = 1) => {
  if (value === null || value === undefined || !Number.isFinite(Number(value))) return "NA";
  return `${Number(value).toFixed(digits)}%`;
};

const pctFraction = (value, digits = 1) => {
  if (value === null || value === undefined || !Number.isFinite(Number(value))) return "NA";
  return `${(100 * Number(value)).toFixed(digits)}%`;
};

const metricRate = (metric) => metric?.within_ee_percent ?? 0;

function metric(label, value, tone = "") {
  return `<div class="metric"><div class="metric-label">${escapeHtml(label)}</div><div class="metric-value ${tone}">${escapeHtml(value)}</div></div>`;
}

function statusPill(value, label = null) {
  return `<span class="status-pill ${escapeHtml(value)}"><span class="status-dot"></span>${escapeHtml(label || transitionLabels[value] || value)}</span>`;
}

function topbar() {
  const views = [
    ["overview", "Performance"],
    ["cases", "All cases"],
    ["evidence", "Evidence gallery"],
    ["robustness", "Robustness"],
    ["experiments", "Experiment ledger"],
    ["method", "Method & sources"],
  ];
  return `<header class="topbar">
    <div class="title-lockup"><h1>${escapeHtml(state.data.title)}</h1><p>${escapeHtml(state.data.subtitle)}</p></div>
    <nav class="view-nav" aria-label="Dashboard views">${views.map(([value, label]) => `<button data-view="${value}" class="${state.view === value ? "active" : ""}">${label}</button>`).join("")}</nav>
  </header>`;
}

function groupedRateChart(rows, options = {}) {
  const width = 920;
  const rowHeight = 42;
  const top = 34;
  const bottom = 36;
  const labelWidth = options.labelWidth || 210;
  const right = 48;
  const height = top + rows.length * rowHeight + bottom;
  const plotWidth = width - labelWidth - right;
  const x = (value) => labelWidth + Math.max(0, Math.min(100, value)) / 100 * plotWidth;
  const ticks = [0, 25, 50, 75, 87, 100].map((value) => `<line x1="${x(value)}" x2="${x(value)}" y1="20" y2="${height - bottom + 4}" class="chart-grid ${value === 87 ? "target-grid" : ""}"/><text x="${x(value)}" y="${height - 10}" text-anchor="middle" class="chart-tick">${value}%</text>`).join("");
  const marks = rows.map((row, index) => {
    const y = top + index * rowHeight;
    const count = row.count !== undefined ? `n=${row.count}` : "";
    return `<text x="${labelWidth - 12}" y="${y + 4}" text-anchor="end" class="chart-label">${escapeHtml(row.label)}</text><text x="${labelWidth - 12}" y="${y + 17}" text-anchor="end" class="chart-sub">${count}</text>
      <line x1="${x(0)}" x2="${x(row.baseline)}" y1="${y - 7}" y2="${y - 7}" stroke="${colours.baseline}" stroke-width="7"><title>Physical: ${row.baseline.toFixed(1)}%</title></line>
      <circle cx="${x(row.baseline)}" cy="${y - 7}" r="4" fill="${colours.baseline}"/><text x="${Math.min(width - 30, x(row.baseline) + 7)}" y="${y - 3}" class="chart-value">${row.baseline.toFixed(1)}</text>
      <line x1="${x(0)}" x2="${x(row.final)}" y1="${y + 9}" y2="${y + 9}" stroke="${colours.final}" stroke-width="7"><title>Final: ${row.final.toFixed(1)}%</title></line>
      <circle cx="${x(row.final)}" cy="${y + 9}" r="4" fill="${colours.final}"/><text x="${Math.min(width - 30, x(row.final) + 7)}" y="${y + 13}" class="chart-value final-value">${row.final.toFixed(1)}</text>`;
  }).join("");
  return `<svg class="chart-svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="${escapeHtml(options.aria || "Within expected error comparison")}">${ticks}${marks}<g transform="translate(${labelWidth},10)"><line x1="0" x2="18" y1="0" y2="0" stroke="${colours.baseline}" stroke-width="6"/><text x="24" y="4" class="chart-legend">Physical</text><line x1="105" x2="123" y1="0" y2="0" stroke="${colours.final}" stroke-width="6"/><text x="129" y="4" class="chart-legend">Final</text><line x1="205" x2="223" y1="0" y2="0" class="target-grid"/><text x="229" y="4" class="chart-legend">87% target</text></g></svg>`;
}

function singleRateChart(rows, options = {}) {
  const width = 920;
  const rowHeight = 32;
  const top = 34;
  const bottom = 36;
  const labelWidth = options.labelWidth || 210;
  const right = 48;
  const height = top + rows.length * rowHeight + bottom;
  const plotWidth = width - labelWidth - right;
  const x = (value) => labelWidth + Math.max(0, Math.min(100, value)) / 100 * plotWidth;
  const ticks = [0, 25, 50, 75, 87, 100].map((value) => `<line x1="${x(value)}" x2="${x(value)}" y1="20" y2="${height - bottom + 4}" class="chart-grid ${value === 87 ? "target-grid" : ""}"/><text x="${x(value)}" y="${height - 10}" text-anchor="middle" class="chart-tick">${value}%</text>`).join("");
  const marks = rows.map((row, index) => {
    const y = top + index * rowHeight;
    const count = row.count !== undefined ? `n=${row.count}` : "";
    return `<text x="${labelWidth - 12}" y="${y + 4}" text-anchor="end" class="chart-label">${escapeHtml(row.label)}</text><text x="${labelWidth - 12}" y="${y + 17}" text-anchor="end" class="chart-sub">${count}</text>
      <line x1="${x(0)}" x2="${x(row.rate)}" y1="${y}" y2="${y}" stroke="${colours.final}" stroke-width="7"><title>Within EE: ${row.rate.toFixed(1)}%</title></line>
      <circle cx="${x(row.rate)}" cy="${y}" r="4" fill="${colours.final}"/><text x="${Math.min(width - 30, x(row.rate) + 7)}" y="${y + 4}" class="chart-value final-value">${row.rate.toFixed(1)}</text>`;
  }).join("");
  return `<svg class="chart-svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="${escapeHtml(options.aria || "Within expected error rates")}">${ticks}${marks}<g transform="translate(${labelWidth},10)"><line x1="0" x2="18" y1="0" y2="0" stroke="${colours.final}" stroke-width="6"/><text x="24" y="4" class="chart-legend">Measured within EE</text><line x1="145" x2="163" y1="0" y2="0" class="target-grid"/><text x="169" y="4" class="chart-legend">87% target</text></g></svg>`;
}

function offsetChart(rows) {
  const width = 920;
  const height = 330;
  const margin = { left: 58, right: 26, top: 28, bottom: 46 };
  const xMin = Math.min(...rows.map((row) => row.offset));
  const xMax = Math.max(...rows.map((row) => row.offset));
  const yMin = 65;
  const yMax = 96;
  const x = (value) => margin.left + (value - xMin) / (xMax - xMin) * (width - margin.left - margin.right);
  const y = (value) => height - margin.bottom - (value - yMin) / (yMax - yMin) * (height - margin.top - margin.bottom);
  const series = [
    ["Development folds 0-3", "development", "#1769aa"],
    ["Untouched fold 4", "confirmation", "#b46b00"],
    ["All 152 descriptive", "all", "#59666d"],
  ];
  const grid = [70, 75, 80, 85, 87, 90, 95].map((value) => `<line x1="${margin.left}" x2="${width - margin.right}" y1="${y(value)}" y2="${y(value)}" class="chart-grid ${value === 87 ? "target-grid" : ""}"/><text x="${margin.left - 10}" y="${y(value) + 4}" text-anchor="end" class="chart-tick">${value}%</text>`).join("");
  const paths = series.map(([label, key, colour]) => {
    const points = rows.map((row) => `${x(row.offset)},${y(row[key].within_ee_percent)}`).join(" ");
    const dots = rows.map((row) => `<circle cx="${x(row.offset)}" cy="${y(row[key].within_ee_percent)}" r="2.3" fill="${colour}"><title>${label}; offset ${row.offset.toFixed(4)}; ${row[key].hits}/${row[key].count} (${row[key].within_ee_percent.toFixed(2)}%)</title></circle>`).join("");
    return `<polyline points="${points}" fill="none" stroke="${colour}" stroke-width="2"/>${dots}`;
  }).join("");
  const selected = rows.find((row) => Math.abs(row.offset - state.data.method.global_log_offset) < 1e-9);
  const selectedLine = selected ? `<line x1="${x(selected.offset)}" x2="${x(selected.offset)}" y1="${margin.top}" y2="${height - margin.bottom}" stroke="#222c31" stroke-dasharray="4 4"/><text x="${x(selected.offset) + 5}" y="${margin.top + 11}" class="chart-value">selected ${selected.offset.toFixed(4)}</text>` : "";
  const xTicks = rows.filter((_, index) => index % 4 === 0).map((row) => `<text x="${x(row.offset)}" y="${height - 16}" text-anchor="middle" class="chart-tick">${row.offset.toFixed(3)}</text>`).join("");
  const legend = series.map(([label, , colour], index) => `<line x1="${margin.left + index * 210}" x2="${margin.left + 18 + index * 210}" y1="12" y2="12" stroke="${colour}" stroke-width="3"/><text x="${margin.left + 25 + index * 210}" y="16" class="chart-legend">${label}</text>`).join("");
  return `<svg class="chart-svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="Within expected error across global log offsets">${grid}${paths}${selectedLine}${xTicks}${legend}<text x="${width / 2}" y="${height - 2}" text-anchor="middle" class="chart-label">Global shifted-log offset</text></svg>`;
}

function featureImportanceChart(rows) {
  const selected = rows.slice(0, 18);
  const width = 920;
  const rowHeight = 26;
  const labelWidth = 390;
  const right = 50;
  const height = 24 + selected.length * rowHeight + 22;
  const maxValue = Math.max(...selected.map((row) => row.importance));
  const x = (value) => labelWidth + value / maxValue * (width - labelWidth - right);
  const marks = selected.map((row, index) => {
    const y = 25 + index * rowHeight;
    return `<text x="${labelWidth - 10}" y="${y + 4}" text-anchor="end" class="chart-label mono-label">${escapeHtml(row.feature)}</text><rect x="${labelWidth}" y="${y - 7}" width="${x(row.importance) - labelWidth}" height="13" fill="#1769aa"><title>${row.feature}: ${row.importance.toFixed(5)}</title></rect><text x="${x(row.importance) + 6}" y="${y + 4}" class="chart-value">${row.importance.toFixed(3)}</text>`;
  }).join("");
  return `<svg class="chart-svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="ExtraTrees feature importance">${marks}</svg>`;
}

function aodComparisonSvg(item) {
  const entries = [
    ["AERONET", item.truth, colours.target],
    ["Physical SIAC", item.baseline, colours.baseline],
    ["ExtraTrees raw", item.unadjusted, colours.raw],
    ["Final", item.adjusted, colours.final],
    ["CAMS 550", item.cams_aod, "#8b5e83"],
    ["MAIAC backstop", item.maiac_aod, "#4f7d72"],
  ].filter(([, value]) => value !== null && value !== undefined && Number.isFinite(Number(value)));
  const width = 900;
  const labelWidth = 175;
  const right = 55;
  const top = 35;
  const rowHeight = 34;
  const height = top + entries.length * rowHeight + 36;
  const maxValue = Math.max(1.0, ...entries.map(([, value]) => Number(value)), item.truth + item.expected_error) * 1.08;
  const x = (value) => labelWidth + Math.max(0, Number(value)) / maxValue * (width - labelWidth - right);
  const grid = [0, .25, .5, .75, 1].map((fraction) => `<line x1="${x(maxValue * fraction)}" x2="${x(maxValue * fraction)}" y1="22" y2="${height - 29}" class="chart-grid"/><text x="${x(maxValue * fraction)}" y="${height - 10}" text-anchor="middle" class="chart-tick">${(maxValue * fraction).toFixed(2)}</text>`).join("");
  const marks = entries.map(([label, value, colour], index) => {
    const y = top + index * rowHeight;
    return `<text x="${labelWidth - 12}" y="${y + 4}" text-anchor="end" class="chart-label">${escapeHtml(label)}</text><line x1="${x(item.truth)}" x2="${x(value)}" y1="${y}" y2="${y}" stroke="#c2c9cd"/><circle cx="${x(value)}" cy="${y}" r="6" fill="${colour}"><title>${label}: ${Number(value).toFixed(5)}</title></circle><text x="${Math.min(width - 32, x(value) + 9)}" y="${y + 4}" class="chart-value">${Number(value).toFixed(3)}</text>`;
  }).join("");
  return `<div class="chart-scroll"><svg class="chart-svg comparison-svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="AOD values for ${escapeHtml(item.site)}"><rect x="${x(Math.max(0, item.truth - item.expected_error))}" y="22" width="${x(item.truth + item.expected_error) - x(Math.max(0, item.truth - item.expected_error))}" height="${height - 51}" class="ee-band"/><line x1="${x(item.truth)}" x2="${x(item.truth)}" y1="18" y2="${height - 29}" stroke="#222c31" stroke-width="1.4"/>${grid}${marks}</svg></div>`;
}

function overviewView() {
  const data = state.data;
  const all = data.metrics.all;
  const transition = data.transition_counts;
  const cohortRows = [
    ["All target cases", data.metrics.all],
    ["Development folds 0-3", data.metrics.development],
    ["Untouched fold 4", data.metrics.confirmation],
    ["External-site seen", data.metrics.seen_sites],
    ["External-site unseen", data.metrics.unseen_sites],
  ].map(([label, row]) => ({ label, count: row.count, baseline: metricRate(row.baseline), final: metricRate(row.final) }));
  const truthRows = data.segments.truth_aod.map((row) => ({ label: row.label, count: row.count, baseline: metricRate(row.baseline), final: metricRate(row.final) }));
  const cloudRows = data.segments.cloud_cover.map((row) => ({ label: row.label, count: row.count, baseline: metricRate(row.baseline), final: metricRate(row.final) }));
  return `<main class="page-view">
    <section class="summary-band">
      <div class="summary-title"><div><h2>Fixed 152-case low-cloud benchmark</h2><p>${escapeHtml(data.benchmark.expected_error)} · ${escapeHtml(data.benchmark.cloud_rule)}</p></div><span class="run-status"><span class="status-dot"></span>152/152 retrievals terminal OK</span></div>
      <div class="metric-row hero-metrics">
        ${metric("Physical within EE", `${all.baseline.hits}/${all.count}`, "muted-value")}
        ${metric("Final within EE", `${all.final.hits}/${all.count}`, "positive")}
        ${metric("Final rate", percent(all.final.within_ee_percent, 2), all.final.hits >= data.benchmark.target_hits ? "positive" : "negative")}
        ${metric("Target", `${data.benchmark.target_hits}/${data.benchmark.count}`, "")}
        ${metric("RMSE", fmt(all.final.rmse, 4), "")}
        ${metric("Bias", signed(all.final.bias, 4), "")}
      </div>
      <div class="progress-track" aria-label="Final within expected error rate"><span style="width:${all.final.within_ee_percent}%"></span><i style="left:${data.benchmark.target_rate}%" title="87% target"></i></div>
    </section>

    <section class="transition-strip">
      <div>${statusPill("gain")}<strong>${transition.gain}</strong><span>physical misses moved inside EE</span></div>
      <div>${statusPill("loss")}<strong>${transition.loss}</strong><span>physical hits moved outside EE</span></div>
      <div>${statusPill("remaining_miss")}<strong>${transition.remaining_miss}</strong><span>outside EE after calibration</span></div>
      <div>${statusPill("retained_hit")}<strong>${transition.retained_hit}</strong><span>inside EE before and after</span></div>
    </section>

    <section class="report-section"><div class="section-header"><div><h2>Retrieval relationship</h2><p>Same 152 cases and axes; grey region is expected error.</p></div><a class="source-button" href="assets/performance_scatter.png" target="_blank" rel="noopener">Open PNG</a></div><img class="wide-figure zoomable" src="assets/performance_scatter.png" alt="Physical and final AOD against AERONET"></section>
    <section class="report-section"><div class="section-header"><div><h2>Case-level error movement</h2><p>Signed error normalized by each case's expected-error width; cases retain a common ordering.</p></div><a class="source-button" href="assets/error_transition.png" target="_blank" rel="noopener">Open PNG</a></div><img class="wide-figure zoomable" src="assets/error_transition.png" alt="Case-level normalized error before and after calibration"></section>

    <section class="report-section"><div class="section-header"><div><h2>Cohort performance</h2><p>Physical and final within-EE rates; the vertical reference is 87%.</p></div></div>${groupedRateChart(cohortRows, { aria: "Within expected error by validation cohort" })}</section>
    <div class="two-column-sections">
      <section class="report-section"><div class="section-header"><div><h2>AERONET AOD strata</h2><p>Fixed truth-AOD bins; denominators shown per row.</p></div></div>${groupedRateChart(truthRows, { labelWidth: 165, aria: "Within expected error by AERONET AOD" })}</section>
      <section class="report-section"><div class="section-header"><div><h2>Scene cloud strata</h2><p>Catalogue cloud cover within the under-20% cohort.</p></div></div>${groupedRateChart(cloudRows, { labelWidth: 165, aria: "Within expected error by scene cloud cover" })}</section>
    </div>
    <section class="disclosure-section"><div class="section-header"><div><h2>Validation disclosures</h2><p>Evidence that constrains interpretation of the aggregate score.</p></div></div><ol>${data.disclosures.map((text) => `<li>${escapeHtml(text)}</li>`).join("")}</ol></section>
  </main>`;
}

function filteredCases() {
  const query = state.search.trim().toLowerCase();
  const rows = state.data.cases.filter((item) => {
    if (query && !`${item.site} ${item.matchup_id}`.toLowerCase().includes(query)) return false;
    if (state.transition !== "all" && item.transition !== state.transition) return false;
    if (state.hitStatus === "hit" && !item.adjusted_within_ee) return false;
    if (state.hitStatus === "miss" && item.adjusted_within_ee) return false;
    if (state.fold !== "all" && String(item.site_fold) !== state.fold) return false;
    if (state.seen === "seen" && !item.external_site_seen) return false;
    if (state.seen === "unseen" && item.external_site_seen) return false;
    if (state.diagnostics === "yes" && !item.diagnostic) return false;
    if (state.diagnostics === "no" && item.diagnostic) return false;
    return true;
  });
  return rows.sort((a, b) => {
    if (state.sort === "improvement") return (Math.abs(b.baseline_error_over_ee) - Math.abs(b.adjusted_error_over_ee)) - (Math.abs(a.baseline_error_over_ee) - Math.abs(a.adjusted_error_over_ee));
    if (state.sort === "seed") return b.seed300_adjusted_std - a.seed300_adjusted_std;
    if (state.sort === "truth") return b.truth - a.truth;
    if (state.sort === "cloud") return (b.cloud_cover ?? -1) - (a.cloud_cover ?? -1);
    if (state.sort === "site") return a.site.localeCompare(b.site);
    return Math.abs(b.adjusted_error_over_ee) - Math.abs(a.adjusted_error_over_ee);
  });
}

function selectedCase() {
  return state.data.cases.find((item) => item.matchup_id === state.selectedId) || state.data.cases[0];
}

function caseSidebar(rows) {
  const list = rows.map((item) => `<button class="case-row ${item.matchup_id === state.selectedId ? "active" : ""}" data-case="${escapeHtml(item.matchup_id)}">
    <span class="case-state ${item.transition}"></span>
    <span class="case-name"><strong>${escapeHtml(item.site)}</strong><span>truth ${fmt(item.truth)} · final ${fmt(item.adjusted)}</span></span>
    <span class="case-score"><strong>${fmt(Math.abs(item.adjusted_error_over_ee), 2)}×</strong><span>fold ${item.site_fold}${item.diagnostic ? " · evidence" : ""}</span></span>
  </button>`).join("");
  return `<aside class="case-sidebar"><div class="filter-panel">
    <label for="case-search">Site or matchup</label><input id="case-search" type="search" value="${escapeHtml(state.search)}" placeholder="Search">
    <div class="filter-grid">
      <label><span>Transition</span><select id="case-transition"><option value="all">All</option>${Object.entries(transitionLabels).map(([value, label]) => `<option value="${value}" ${state.transition === value ? "selected" : ""}>${label}</option>`).join("")}</select></label>
      <label><span>Final status</span><select id="case-hit"><option value="all">All</option><option value="hit" ${state.hitStatus === "hit" ? "selected" : ""}>Within EE</option><option value="miss" ${state.hitStatus === "miss" ? "selected" : ""}>Outside EE</option></select></label>
      <label><span>Site fold</span><select id="case-fold"><option value="all">All folds</option>${[0,1,2,3,4].map((value) => `<option value="${value}" ${state.fold === String(value) ? "selected" : ""}>Fold ${value}</option>`).join("")}</select></label>
      <label><span>External sites</span><select id="case-seen"><option value="all">Seen + unseen</option><option value="seen" ${state.seen === "seen" ? "selected" : ""}>Seen</option><option value="unseen" ${state.seen === "unseen" ? "selected" : ""}>Unseen</option></select></label>
      <label><span>Evidence</span><select id="case-diagnostics"><option value="all">With + without</option><option value="yes" ${state.diagnostics === "yes" ? "selected" : ""}>Has spatial/cost</option><option value="no" ${state.diagnostics === "no" ? "selected" : ""}>No cost cube</option></select></label>
      <label><span>Sort</span><select id="case-sort"><option value="final_error" ${state.sort === "final_error" ? "selected" : ""}>Final |error| / EE</option><option value="improvement" ${state.sort === "improvement" ? "selected" : ""}>Error reduction</option><option value="seed" ${state.sort === "seed" ? "selected" : ""}>Seed spread</option><option value="truth" ${state.sort === "truth" ? "selected" : ""}>AERONET AOD</option><option value="cloud" ${state.sort === "cloud" ? "selected" : ""}>Cloud cover</option><option value="site" ${state.sort === "site" ? "selected" : ""}>Site</option></select></label>
    </div><div class="filter-count">${rows.length} of ${state.data.cases.length} cases</div>
  </div><div class="case-list">${list || '<div class="empty-state">No matching cases.</div>'}</div></aside>`;
}

function predictionsTab(item) {
  const rows = [
    ["Physical SIAC", item.baseline, item.baseline_error, item.baseline_within_ee],
    ["ExtraTrees raw", item.unadjusted, item.unadjusted_error, item.unadjusted_within_ee],
    ["Final offset", item.adjusted, item.adjusted_error, item.adjusted_within_ee],
  ];
  return `<section class="case-section"><div class="section-header"><div><h3>Scalar AOD comparison</h3><p>Expected-error range ${fmt(item.truth - item.expected_error)} to ${fmt(item.truth + item.expected_error)}.</p></div></div>${aodComparisonSvg(item)}
    <div class="table-wrap"><table class="data-table"><thead><tr><th>Output</th><th>AOD</th><th>Signed error</th><th>Error / EE</th><th>Within EE</th></tr></thead><tbody>${rows.map(([label, value, error, hit]) => `<tr><td><strong>${label}</strong></td><td class="mono">${fmt(value, 6)}</td><td class="mono">${signed(error, 6)}</td><td class="mono">${signed(error / item.expected_error, 3)}</td><td>${hit ? '<span class="hit-text">Yes</span>' : '<span class="miss-text">No</span>'}</td></tr>`).join("")}</tbody></table></div>
    <div class="context-grid">
      <div><span>CAMS AOD 550</span><strong>${fmt(item.cams_aod, 4)}</strong></div><div><span>MAIAC backstop</span><strong>${fmt(item.maiac_aod, 4)}</strong></div><div><span>Band argmin spread</span><strong>${fmt(item.band_argmin_spread, 4)}</strong></div><div><span>Solver cost / band</span><strong>${fmt(item.solver_cost, 4)}</strong></div><div><span>Scene cloud</span><strong>${item.cloud_cover === null ? "NA" : percent(item.cloud_cover)}</strong></div><div><span>Seed-300 range</span><strong>${fmt(item.seed300_adjusted_min)}–${fmt(item.seed300_adjusted_max)}</strong></div>
    </div></section>`;
}

function imageTab(item, kind) {
  const diagnostic = item.diagnostic;
  if (!diagnostic) return `<section class="case-section"><div class="empty-state">No diagnostic cost cube was generated for this retained case. Scalar predictions and operational features remain available.</div></section>`;
  const spatial = kind === "spatial";
  const source = spatial ? diagnostic.spatial_image : diagnostic.diagnostic_image;
  const title = spatial ? "Spatial evidence" : "Spectral and cost evidence";
  const caption = spatial ? "TOA, surface prior, atmospheric backstop, AOD maps, per-band minima, and truth-node residuals." : "Scene and station-window cost curves, per-band residuals, reflectance spectrum, closure, and scalar outputs.";
  return `<section class="case-section"><div class="section-header"><div><h3>${title}</h3><p>${escapeHtml(diagnostic.source)}</p></div><a class="source-button" href="${source}" target="_blank" rel="noopener">Open PNG</a></div><figure class="evidence-figure"><img class="zoomable" src="${source}" alt="${title} for ${escapeHtml(item.site)}"><figcaption>${caption}</figcaption></figure><div class="source-links"><a href="${diagnostic.canonical_result_url}" target="_blank" rel="noopener">Physical result JSON</a><a href="${diagnostic.diagnostic_result_url}" target="_blank" rel="noopener">Diagnostic result JSON</a><a href="${diagnostic.cost_cube_url}" target="_blank" rel="noopener">Cost cube NPZ</a></div></section>`;
}

function featureRows(item, selectedOnly = true) {
  const importance = new Map(state.data.robustness.model_feature_importance.map((row) => [row.feature, row.importance]));
  const names = selectedOnly ? state.data.robustness.model_feature_importance.map((row) => row.feature) : Object.keys(item.operational_features).sort();
  return names.map((name) => ({ name, value: item.operational_features[name], importance: importance.get(name) })).filter((row) => row.value !== undefined);
}

function featuresTab(item) {
  const rows = featureRows(item, true);
  return `<section class="case-section"><div class="section-header"><div><h3>Selected operational features</h3><p>${rows.length} available values from the 35-feature model subset; importance is global model dependence, not a per-case contribution.</p></div></div><div class="table-wrap"><table class="data-table"><thead><tr><th>Feature</th><th>Case value</th><th>Global tree importance</th></tr></thead><tbody>${rows.map((row) => `<tr><td class="mono wrap">${escapeHtml(row.name)}</td><td class="mono">${fmt(row.value, 6)}</td><td class="mono">${fmt(row.importance, 5)}</td></tr>`).join("")}</tbody></table></div></section>`;
}

function kvTable(title, object) {
  const rows = Object.entries(object || {}).filter(([, value]) => ["string", "number", "boolean"].includes(typeof value));
  return `<div class="kv-panel"><h4>${escapeHtml(title)}</h4><table class="kv-table"><tbody>${rows.map(([key, value]) => `<tr><th>${escapeHtml(key)}</th><td class="mono wrap">${typeof value === "number" ? fmt(value, 7) : escapeHtml(value)}</td></tr>`).join("")}</tbody></table></div>`;
}

function rawTab(item) {
  return `<section class="case-section"><div class="section-header"><div><h3>Raw operational receipts</h3><p>Values shown without interpretation or case routing.</p></div></div><div class="raw-grid">${kvTable("Solver", item.solver)}${kvTable("Retrieval extraction", item.retrieval_extraction)}${kvTable("Surface prior BOA", item.prior_boa)}${kvTable("Surface uncertainty", item.prior_unc)}${kvTable("Anchor iteration", item.anchor_iterate)}</div><details class="all-features"><summary>All ${Object.keys(item.operational_features).length} extracted operational features</summary><div class="table-wrap"><table class="data-table"><tbody>${featureRows(item, false).map((row) => `<tr><td class="mono wrap">${escapeHtml(row.name)}</td><td class="mono">${fmt(row.value, 8)}</td></tr>`).join("")}</tbody></table></div></details></section>`;
}

function caseDetail(item, rows) {
  const index = Math.max(0, rows.findIndex((row) => row.matchup_id === item.matchup_id));
  const previous = rows[(index - 1 + rows.length) % rows.length];
  const next = rows[(index + 1) % rows.length];
  const body = state.detailTab === "spatial" ? imageTab(item, "spatial") : state.detailTab === "spectral" ? imageTab(item, "spectral") : state.detailTab === "features" ? featuresTab(item) : state.detailTab === "raw" ? rawTab(item) : predictionsTab(item);
  return `<main class="case-main"><header class="case-header"><div class="case-header-top"><div><div class="case-title-line"><h2>${escapeHtml(item.site)}</h2>${statusPill(item.transition)}</div><div class="matchup">${escapeHtml(item.matchup_id)}</div></div><div class="case-actions"><button class="icon-button" data-case="${escapeHtml(previous?.matchup_id)}" title="Previous case" aria-label="Previous case">←</button><button class="icon-button" data-case="${escapeHtml(next?.matchup_id)}" title="Next case" aria-label="Next case">→</button></div></div>
    <div class="metric-row case-metrics">${metric("AERONET", fmt(item.truth))}${metric("Physical", fmt(item.baseline), item.baseline_within_ee ? "" : "negative")}${metric("ExtraTrees raw", fmt(item.unadjusted), item.unadjusted_within_ee ? "" : "negative")}${metric("Final", fmt(item.adjusted), item.adjusted_within_ee ? "positive" : "negative")}${metric("Final error / EE", `${fmt(Math.abs(item.adjusted_error_over_ee), 2)}×`, item.adjusted_within_ee ? "positive" : "negative")}${metric("Cloud", item.cloud_cover === null ? "NA" : percent(item.cloud_cover))}</div>
    <div class="case-tags"><span>Fold <strong>${item.site_fold}</strong></span><span>External training site <strong>${item.external_site_seen ? "seen" : "unseen"}</strong></span><span>Seed std <strong>${fmt(item.seed300_adjusted_std, 5)}</strong></span><span>Diagnostic <strong>${item.diagnostic ? item.diagnostic.source : "not generated"}</strong></span></div></header>
    <nav class="detail-tabs" aria-label="Case evidence views"><button data-tab="predictions" class="${state.detailTab === "predictions" ? "active" : ""}">Predictions</button><button data-tab="spatial" class="${state.detailTab === "spatial" ? "active" : ""}">Spatial</button><button data-tab="spectral" class="${state.detailTab === "spectral" ? "active" : ""}">Spectral & cost</button><button data-tab="features" class="${state.detailTab === "features" ? "active" : ""}">Model features</button><button data-tab="raw" class="${state.detailTab === "raw" ? "active" : ""}">Raw receipts</button></nav><div class="case-body">${body}</div></main>`;
}

function casesView() {
  const rows = filteredCases();
  let item = selectedCase();
  if (rows.length && !rows.some((row) => row.matchup_id === item.matchup_id)) {
    item = rows[0];
    state.selectedId = item.matchup_id;
  }
  return `<div class="explorer-layout">${caseSidebar(rows)}${rows.length ? caseDetail(item, rows) : '<main class="case-main"><div class="empty-state">No matching cases.</div></main>'}</div>`;
}

function evidenceView() {
  const query = state.search.trim().toLowerCase();
  const rows = state.data.cases.filter((item) => item.diagnostic && (!query || `${item.site} ${item.matchup_id}`.toLowerCase().includes(query)) && (state.galleryTransition === "all" || item.transition === state.galleryTransition));
  const cards = rows.map((item) => {
    const image = state.galleryMode === "spatial" ? item.diagnostic.spatial_image : item.diagnostic.diagnostic_image;
    return `<article class="gallery-item"><button data-case="${escapeHtml(item.matchup_id)}"><img src="${image}" alt="${state.galleryMode} evidence for ${escapeHtml(item.site)}" loading="lazy"><div class="gallery-meta"><div><strong>${escapeHtml(item.site)}</strong><span>${escapeHtml(item.matchup_id)}</span></div>${statusPill(item.transition)}</div><div class="gallery-values"><span>truth ${fmt(item.truth)}</span><span>physical ${fmt(item.baseline)}</span><span>final ${fmt(item.adjusted)}</span><span>${fmt(Math.abs(item.adjusted_error_over_ee), 2)}× EE</span></div></button></article>`;
  }).join("");
  return `<main class="page-view gallery-view"><div class="page-heading"><div><h2>Spatial, spectral, and cost evidence</h2><p>${state.data.diagnostics.case_count} diagnostic examples, including all ${state.data.diagnostics.transition_case_count} gains, losses, and remaining misses.</p></div></div><div class="gallery-toolbar"><label><span>Site or matchup</span><input id="gallery-search" type="search" value="${escapeHtml(state.search)}" placeholder="Search"></label><label><span>Transition</span><select id="gallery-transition"><option value="all">All transitions</option>${Object.entries(transitionLabels).map(([value, label]) => `<option value="${value}" ${state.galleryTransition === value ? "selected" : ""}>${label}</option>`).join("")}</select></label><div class="segmented" role="group" aria-label="Evidence image"><button data-gallery-mode="spatial" class="${state.galleryMode === "spatial" ? "active" : ""}">Spatial</button><button data-gallery-mode="diagnostic" class="${state.galleryMode === "diagnostic" ? "active" : ""}">Spectral & cost</button></div><div class="gallery-count">${rows.length} cases</div></div><div class="gallery-grid">${cards}</div></main>`;
}

function metricTableRows(rows) {
  return rows.map((row) => `<tr><td><strong>${escapeHtml(row.label)}</strong><span class="subtext">n=${row.count}</span></td><td class="mono">${row.baseline.hits}/${row.count}</td><td class="mono">${percent(row.baseline.within_ee_percent, 2)}</td><td class="mono">${row.final.hits}/${row.count}</td><td class="mono">${percent(row.final.within_ee_percent, 2)}</td><td class="mono">${fmt(row.final.rmse, 4)}</td><td class="mono">${fmt(row.final.mae, 4)}</td><td class="mono">${signed(row.final.bias, 4)}</td></tr>`).join("");
}

function robustnessView() {
  const data = state.data;
  const cohortRows = [
    ["All target", data.metrics.all], ["Development folds 0-3", data.metrics.development], ["Untouched fold 4", data.metrics.confirmation], ["External-site seen", data.metrics.seen_sites], ["External-site unseen", data.metrics.unseen_sites],
  ].map(([label, row]) => ({ label, count: row.count, baseline: row.baseline, final: row.final }));
  const seedRows = data.robustness.seed_tree.map((row) => `<tr><td class="mono wrap">${escapeHtml(row.variant)}</td><td class="mono">${row.target.hits}/${row.target.count}</td><td class="mono">${percent(row.target.within_ee_percent, 2)}</td><td class="mono">${row.confirmation.hits}/${row.confirmation.count}</td><td class="mono">${row.unseen_sites.hits}/${row.unseen_sites.count}</td><td class="mono">${row.external_holdout.hits}/${row.external_holdout.count}</td><td class="mono">${fmt(row.target.rmse, 4)}</td></tr>`).join("");
  const nestedRows = data.robustness.nested_folds.map((row) => `<tr><td>Fold ${row.held_out_fold}</td><td class="mono">${signed(row.selected_offset, 4)}</td><td class="mono">${row.development.candidate.hits}/${row.development.count}</td><td class="mono">${row.confirmation.candidate.hits}/${row.confirmation.count}</td><td class="mono">${percent(row.confirmation.candidate.within_ee_percent, 2)}</td><td class="mono">${fmt(row.confirmation.candidate.rmse, 4)}</td></tr>`).join("");
  const ablations = Object.entries(data.robustness.feature_ablation.variants).map(([name, row]) => ({ name, ...row }));
  const ablationRows = ablations.map((row) => `<tr><td class="mono wrap">${escapeHtml(row.name)}</td><td class="mono">${row.feature_count}</td><td class="mono">${row.target_all.candidate.hits}/152</td><td class="mono">${row.target_confirmation.candidate.hits}/26</td><td class="mono">${row.target_unseen_sites.candidate.hits}/53</td><td class="mono">${row.external_holdout.candidate.hits}/123</td><td class="mono">${fmt(row.target_all.candidate.rmse, 4)}</td></tr>`).join("");
  return `<main class="page-view"><div class="page-heading"><div><h2>Robustness and dependence audit</h2><p>Held-out cohorts, offset selection, random seeds, tree counts, feature removals, and serialized-model reproduction.</p></div></div>
    <section class="report-section"><div class="section-header"><div><h2>Validation cohorts</h2><p>All denominators and rates are retained.</p></div></div><div class="table-wrap"><table class="data-table"><thead><tr><th>Cohort</th><th>Physical hits</th><th>Physical rate</th><th>Final hits</th><th>Final rate</th><th>Final RMSE</th><th>Final MAE</th><th>Final bias</th></tr></thead><tbody>${metricTableRows(cohortRows)}</tbody></table></div></section>
    <section class="report-section"><div class="section-header"><div><h2>Global offset selection</h2><p>Each point applies one offset to every case; folds 0-3 choose -0.0125 and fold 4 is untouched.</p></div></div>${offsetChart(data.robustness.offset_curve)}</section>
    <section class="report-section"><div class="section-header"><div><h2>Seed and tree-count replay</h2><p>Every full-schema replay meets the 133/152 target; external holdout is reported independently.</p></div><span class="run-status"><span class="status-dot"></span>All variants ≥133</span></div><div class="table-wrap"><table class="data-table"><thead><tr><th>Variant</th><th>Target hits</th><th>Target rate</th><th>Fold-4 hits</th><th>Unseen-site hits</th><th>External holdout</th><th>Target RMSE</th></tr></thead><tbody>${seedRows}</tbody></table></div></section>
    <section class="report-section"><div class="section-header"><div><h2>Nested offset replay</h2><p>The ExtraTrees recipe is fixed; each row selects only the scalar offset without its held-out fold.</p></div></div><div class="metric-row compact-metrics">${metric("Aggregate", `${data.metrics.nested_offset.candidate.hits}/${data.metrics.nested_offset.count}`, "positive")}${metric("Rate", percent(data.metrics.nested_offset.candidate.within_ee_percent, 2), "positive")}${metric("RMSE", fmt(data.metrics.nested_offset.candidate.rmse, 4))}${metric("MAE", fmt(data.metrics.nested_offset.candidate.mae, 4))}</div><div class="table-wrap"><table class="data-table"><thead><tr><th>Held-out fold</th><th>Selected offset</th><th>Development hits</th><th>Held-out hits</th><th>Held-out rate</th><th>Held-out RMSE</th></tr></thead><tbody>${nestedRows}</tbody></table></div></section>
    <section class="report-section"><div class="section-header"><div><h2>Operational feature dependence</h2><p>Absolute-error ExtraTrees importance from the exported 1,500-tree model; this is not causal attribution.</p></div></div>${featureImportanceChart(data.robustness.model_feature_importance)}</section>
    <section class="report-section"><div class="section-header"><div><h2>Feature-removal sensitivity</h2><p>Fixed 300-tree diagnostic ablations; none were selected using target labels. The no-MAIAC point estimate is not seed-stable at the target threshold.</p></div></div><div class="table-wrap"><table class="data-table"><thead><tr><th>Feature schema</th><th>Features</th><th>Target</th><th>Fold 4</th><th>Unseen sites</th><th>External holdout</th><th>Target RMSE</th></tr></thead><tbody>${ablationRows}</tbody></table></div></section>
    <section class="report-section"><div class="section-header"><div><h2>Serialized artifact reproduction</h2><p>Cross-host ExtraTrees floating-point variation is bounded and changes no within-EE classification.</p></div><a class="source-button" href="${data.sources.model_manifest}" target="_blank" rel="noopener">Manifest JSON</a></div><div class="metric-row compact-metrics">${metric("Model SHA-256", data.method.model_sha256.slice(0, 12) + "…")}${metric("Max AOD delta", fmt(data.robustness.artifact_reproduction.max_abs_prediction_delta, 6))}${metric("Tolerance", fmt(data.robustness.artifact_reproduction.prediction_tolerance, 4))}${metric("Hit disagreements", data.robustness.artifact_reproduction.within_ee_classification_disagreements, "positive")}</div></section>
  </main>`;
}

function experimentsView() {
  const categories = [...new Set(state.data.experiments.map((row) => row.category))];
  const rows = state.data.experiments.filter((row) => state.experimentCategory === "all" || row.category === state.experimentCategory);
  const chartRows = rows.filter((row) => row.count === 152).map((row) => ({ label: row.method, count: row.count, rate: 100 * row.hits / row.count }));
  return `<main class="page-view"><div class="page-heading"><div><h2>Experiment ledger</h2><p>Measured candidates, diagnostic ceilings, and selection outcomes; source receipts are linked per row.</p></div><label class="inline-filter"><span>Category</span><select id="experiment-category"><option value="all">All categories</option>${categories.map((value) => `<option value="${escapeHtml(value)}" ${state.experimentCategory === value ? "selected" : ""}>${escapeHtml(value)}</option>`).join("")}</select></label></div>
    <section class="report-section"><div class="section-header"><div><h2>Target-cohort comparison</h2><p>Methods with the common 152-case denominator; each bar is one measured within-EE rate, not a recommendation ranking.</p></div></div>${singleRateChart(chartRows, { labelWidth: 350, aria: "Experiment within expected error rates" })}</section>
    <section class="report-section"><div class="table-wrap"><table class="data-table ledger-table"><thead><tr><th>Category</th><th>Method</th><th>Scope</th><th>Target result</th><th>External holdout</th><th>Operational</th><th>Status</th><th>Evidence note</th><th>Source</th></tr></thead><tbody>${rows.map((row) => `<tr><td>${escapeHtml(row.category)}</td><td><strong>${escapeHtml(row.method)}</strong></td><td>${escapeHtml(row.scope)}</td><td class="mono">${row.hits}/${row.count} (${percent(100 * row.hits / row.count, 1)})</td><td class="mono">${row.holdout_hits === undefined ? "NA" : `${row.holdout_hits}/${row.holdout_count} (${percent(100 * row.holdout_hits / row.holdout_count, 1)})`}</td><td>${row.operational ? "Yes" : "No"}</td><td><span class="ledger-status">${escapeHtml(row.status)}</span></td><td>${escapeHtml(row.note)}</td><td><a href="${row.source}" target="_blank" rel="noopener">Receipt</a></td></tr>`).join("")}</tbody></table></div></section>
  </main>`;
}

function methodView() {
  const data = state.data;
  const constraints = data.method.inference_constraints;
  return `<main class="page-view method-view"><div class="page-heading"><div><h2>Method, provenance, and files</h2><p>Frozen operational recipe and validation boundaries for this published snapshot.</p></div><a class="source-button" href="${data.method.model_url}" target="_blank" rel="noopener">Model artifact</a></div>
    <section class="method-flow">
      <div><span>1</span><h3>Physical retrieval</h3><p>${escapeHtml(data.method.physical_retrieval)}</p></div>
      <div><span>2</span><h3>Single surface prior</h3><p>${escapeHtml(data.method.surface_prior)}</p></div>
      <div><span>3</span><h3>Uniform calibrator</h3><p>${escapeHtml(data.method.calibrator.name)} · ${data.method.selected_feature_count} selected operational features</p></div>
      <div><span>4</span><h3>Global offset</h3><p>${signed(data.method.global_log_offset, 4)} in shifted-log AOD space</p></div>
    </section>
    <section class="report-section"><div class="section-header"><div><h2>Exact inference expression</h2><p>Clipping occurs between the model prediction and global offset.</p></div></div><pre class="formula"><code>${escapeHtml(data.method.formula)}</code></pre><div class="constraint-row"><span>One surface prior <strong>${constraints.one_surface_prior ? "yes" : "no"}</strong></span><span>Case routing <strong>${constraints.case_routing ? "yes" : "no"}</strong></span><span>Geography features <strong>${constraints.geography_features ? "yes" : "no"}</strong></span><span>AERONET at inference <strong>${constraints.aeronet_at_inference ? "yes" : "no"}</strong></span><span>Target label at inference <strong>${constraints.target_label_at_inference ? "yes" : "no"}</strong></span></div></section>
    <section class="report-section"><div class="section-header"><div><h2>Training and selection boundaries</h2><p>Counts refer to this frozen snapshot.</p></div></div><div class="data-grid"><div class="kv-panel"><h4>External model fit</h4><table class="kv-table"><tbody><tr><th>Records</th><td class="mono">${data.method.training.count}</td></tr><tr><th>Sites</th><td class="mono">${data.method.training.site_count}</td></tr><tr><th>Cloud domain</th><td class="mono">&lt; ${data.method.training.max_scene_cloud_cover_percent}%</td></tr><tr><th>Target labels in fit</th><td class="mono">${data.method.training.target_labels_used_for_fit ? "yes" : "no"}</td></tr><tr><th>Feature schema</th><td class="mono">${data.method.feature_schema_count}</td></tr><tr><th>Selected features</th><td class="mono">${data.method.selected_feature_count}</td></tr></tbody></table></div><div class="kv-panel"><h4>Target-domain validation</h4><table class="kv-table"><tbody><tr><th>Development</th><td class="mono">folds 0-3 · 126 cases</td></tr><tr><th>Untouched confirmation</th><td class="mono">fold 4 · 26 cases</td></tr><tr><th>Offset nested replay</th><td class="mono">5 held-out site folds</td></tr><tr><th>External holdout</th><td class="mono">123 records</td></tr><tr><th>External-site unseen</th><td class="mono">53 target cases</td></tr><tr><th>Random seed replay</th><td class="mono">5 seeds</td></tr></tbody></table></div></div></section>
    <section class="report-section"><div class="section-header"><div><h2>Source artifacts</h2><p>Published JSON, CSV, model, and diagnostic receipts.</p></div></div><div class="source-index">${Object.entries(data.sources).map(([label, url]) => `<a href="${url}" target="_blank" rel="noopener"><strong>${escapeHtml(label.replaceAll("_", " "))}</strong><span>${escapeHtml(url)}</span></a>`).join("")}<a href="${data.method.model_url}" target="_blank" rel="noopener"><strong>serialized model</strong><span>${escapeHtml(data.method.model_url)}</span></a></div></section>
    <section class="disclosure-section"><div class="section-header"><div><h2>Interpretation boundaries</h2><p>No inference is made beyond these recorded checks.</p></div></div><ol>${data.disclosures.map((text) => `<li>${escapeHtml(text)}</li>`).join("")}</ol></section>
  </main>`;
}

function render() {
  let body;
  if (state.view === "cases") body = casesView();
  else if (state.view === "evidence") body = evidenceView();
  else if (state.view === "robustness") body = robustnessView();
  else if (state.view === "experiments") body = experimentsView();
  else if (state.view === "method") body = methodView();
  else body = overviewView();
  app.innerHTML = `${topbar()}${body}<dialog class="image-dialog"><div class="dialog-bar"><strong></strong><button class="icon-button" data-close-dialog title="Close" aria-label="Close">×</button></div><img alt="Expanded evidence figure"></dialog>`;
  bindEvents();
}

function setView(view) {
  state.view = view;
  history.replaceState(null, "", view === "overview" ? "#overview" : `#${view}`);
  render();
}

function selectCase(matchupId) {
  if (!matchupId) return;
  state.selectedId = matchupId;
  state.view = "cases";
  history.replaceState(null, "", `#case=${encodeURIComponent(matchupId)}`);
  render();
}

function updateState(key, value, focusId = null) {
  state[key] = value;
  render();
  if (focusId) {
    const input = document.getElementById(focusId);
    input?.focus();
    input?.setSelectionRange?.(value.length, value.length);
  }
}

function openImage(source, alt) {
  const dialog = document.querySelector(".image-dialog");
  dialog.querySelector("strong").textContent = alt;
  const image = dialog.querySelector("img");
  image.src = source;
  image.alt = alt;
  dialog.showModal();
}

function bindEvents() {
  document.querySelectorAll("[data-view]").forEach((button) => button.addEventListener("click", () => setView(button.dataset.view)));
  document.querySelectorAll("[data-case]").forEach((button) => button.addEventListener("click", () => selectCase(button.dataset.case)));
  document.querySelectorAll("[data-tab]").forEach((button) => button.addEventListener("click", () => { state.detailTab = button.dataset.tab; render(); }));
  document.querySelectorAll("[data-gallery-mode]").forEach((button) => button.addEventListener("click", () => { state.galleryMode = button.dataset.galleryMode; render(); }));
  document.querySelectorAll(".zoomable").forEach((image) => image.addEventListener("click", () => openImage(image.src, image.alt)));
  document.querySelector("[data-close-dialog]")?.addEventListener("click", () => document.querySelector(".image-dialog")?.close());
  document.getElementById("case-search")?.addEventListener("input", (event) => updateState("search", event.target.value, "case-search"));
  document.getElementById("case-transition")?.addEventListener("change", (event) => updateState("transition", event.target.value));
  document.getElementById("case-hit")?.addEventListener("change", (event) => updateState("hitStatus", event.target.value));
  document.getElementById("case-fold")?.addEventListener("change", (event) => updateState("fold", event.target.value));
  document.getElementById("case-seen")?.addEventListener("change", (event) => updateState("seen", event.target.value));
  document.getElementById("case-diagnostics")?.addEventListener("change", (event) => updateState("diagnostics", event.target.value));
  document.getElementById("case-sort")?.addEventListener("change", (event) => updateState("sort", event.target.value));
  document.getElementById("gallery-search")?.addEventListener("input", (event) => updateState("search", event.target.value, "gallery-search"));
  document.getElementById("gallery-transition")?.addEventListener("change", (event) => updateState("galleryTransition", event.target.value));
  document.getElementById("experiment-category")?.addEventListener("change", (event) => updateState("experimentCategory", event.target.value));
}

document.addEventListener("keydown", (event) => {
  if (state.view !== "cases" || ["INPUT", "SELECT", "TEXTAREA"].includes(document.activeElement?.tagName)) return;
  if (!["ArrowLeft", "ArrowRight"].includes(event.key)) return;
  const rows = filteredCases();
  const index = rows.findIndex((row) => row.matchup_id === state.selectedId);
  if (index < 0 || !rows.length) return;
  const delta = event.key === "ArrowLeft" ? -1 : 1;
  selectCase(rows[(index + delta + rows.length) % rows.length].matchup_id);
});

fetch("data/dashboard.json")
  .then((response) => {
    if (!response.ok) throw new Error(`HTTP ${response.status}`);
    return response.json();
  })
  .then((data) => {
    state.data = data;
    const caseMatch = location.hash.match(/^#case=(.+)$/);
    if (caseMatch) {
      const requested = decodeURIComponent(caseMatch[1]);
      state.selectedId = data.cases.some((item) => item.matchup_id === requested) ? requested : data.cases[0].matchup_id;
      state.view = "cases";
    } else {
      const requestedView = location.hash.replace(/^#/, "");
      state.view = ["overview", "cases", "evidence", "robustness", "experiments", "method"].includes(requestedView) ? requestedView : "overview";
      state.selectedId = data.cases.find((item) => item.transition === "remaining_miss")?.matchup_id || data.cases[0].matchup_id;
    }
    render();
  })
  .catch((error) => {
    app.innerHTML = `<div class="empty-state">Unable to load the published snapshot: ${escapeHtml(error.message)}</div>`;
  });
