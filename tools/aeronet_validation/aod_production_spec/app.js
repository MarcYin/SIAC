(() => {
  "use strict";

  const byId = (id) => document.getElementById(id);

  function initializeToc() {
    const links = [...document.querySelectorAll(".toc a[href^='#']")];
    const sections = links
      .map((link) => document.querySelector(link.getAttribute("href")))
      .filter(Boolean);
    if (!links.length || !sections.length || !("IntersectionObserver" in window)) return;

    const activate = (id) => {
      links.forEach((link) => {
        const active = link.getAttribute("href") === `#${id}`;
        link.classList.toggle("active", active);
        if (active) link.setAttribute("aria-current", "location");
        else link.removeAttribute("aria-current");
      });
    };
    const visible = new Map();
    const observer = new IntersectionObserver(
      (entries) => {
        entries.forEach((entry) => visible.set(entry.target.id, entry));
        const candidate = [...visible.values()]
          .filter((entry) => entry.isIntersecting)
          .sort((left, right) => left.boundingClientRect.top - right.boundingClientRect.top)[0];
        if (candidate) activate(candidate.target.id);
      },
      { rootMargin: "-72px 0px -68% 0px", threshold: [0, 0.05] },
    );
    sections.forEach((section) => observer.observe(section));
    activate(sections[0].id);
  }

  function initializeFeatureFilter() {
    const search = byId("feature-search");
    const group = byId("feature-group");
    const count = byId("feature-count");
    const rows = [...document.querySelectorAll("[data-feature-row]")];
    if (!search || !group || !count || !rows.length) return;

    const update = () => {
      const query = search.value.trim().toLowerCase();
      const selectedGroup = group.value;
      let visible = 0;
      rows.forEach((row) => {
        const matchesQuery = !query || row.dataset.search.includes(query);
        const matchesGroup = selectedGroup === "all" || row.dataset.group === selectedGroup;
        row.hidden = !(matchesQuery && matchesGroup);
        if (!row.hidden) visible += 1;
      });
      count.textContent = `${visible} of ${rows.length} features`;
    };
    search.addEventListener("input", update);
    group.addEventListener("change", update);
    update();
  }

  function initializeSourceFilter() {
    const search = byId("source-search");
    const count = byId("source-count");
    const rows = [...document.querySelectorAll("[data-source-row]")];
    if (!search || !count || !rows.length) return;

    const update = () => {
      const query = search.value.trim().toLowerCase();
      let visible = 0;
      rows.forEach((row) => {
        row.hidden = Boolean(query) && !row.dataset.search.includes(query);
        if (!row.hidden) visible += 1;
      });
      count.textContent = `${visible} of ${rows.length} files`;
    };
    search.addEventListener("input", update);
    update();
  }

  async function copyText(text) {
    if (navigator.clipboard && window.isSecureContext) {
      await navigator.clipboard.writeText(text);
      return;
    }
    const area = document.createElement("textarea");
    area.value = text;
    area.style.position = "fixed";
    area.style.opacity = "0";
    document.body.appendChild(area);
    area.select();
    document.execCommand("copy");
    area.remove();
  }

  function initializeCopyButtons() {
    document.querySelectorAll("[data-copy]").forEach((button) => {
      button.addEventListener("click", async () => {
        const code = button.parentElement?.querySelector("code");
        if (!code) return;
        const prior = button.textContent;
        try {
          await copyText(code.textContent);
          button.textContent = "Copied";
        } catch (_error) {
          button.textContent = "Failed";
        }
        window.setTimeout(() => { button.textContent = prior; }, 1200);
      });
    });
  }

  function initializeChecklist() {
    const checklist = document.querySelector("[data-checklist]");
    if (!checklist) return;
    const inputs = [...checklist.querySelectorAll("input[type='checkbox']")];
    const key = "siac-aod-production-spec-checklist-v1";
    try {
      const saved = JSON.parse(window.localStorage.getItem(key) || "[]");
      inputs.forEach((input, index) => { input.checked = Boolean(saved[index]); });
    } catch (_error) {
      // Storage may be disabled; the checklist remains fully usable in memory.
    }
    inputs.forEach((input) => {
      input.addEventListener("change", () => {
        try {
          window.localStorage.setItem(key, JSON.stringify(inputs.map((item) => item.checked)));
        } catch (_error) {
          // No persistent storage available.
        }
      });
    });
  }

  initializeToc();
  initializeFeatureFilter();
  initializeSourceFilter();
  initializeCopyButtons();
  initializeChecklist();
})();
