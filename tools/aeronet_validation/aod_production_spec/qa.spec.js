const { expect, test } = require("@playwright/test");

const reportUrl = `https://gws-access.jasmin.ac.uk/public/nceo_isp/siac_refactor/reports/aod-production-reproduction-spec-20260713/?qa=${Date.now()}`;

function captureRuntimeFailures(page) {
  const failures = [];
  page.on("pageerror", (error) => failures.push(`pageerror: ${error.message}`));
  page.on("console", (message) => {
    if (message.type() === "error") failures.push(`console: ${message.text()}`);
  });
  page.on("requestfailed", (request) => {
    const reason = request.failure()?.errorText || "unknown";
    if (reason !== "net::ERR_ABORTED") failures.push(`request: ${request.url()} (${reason})`);
  });
  return failures;
}

async function expectNoDocumentOverflow(page) {
  const dimensions = await page.evaluate(() => ({
    viewport: document.documentElement.clientWidth,
    scroll: document.documentElement.scrollWidth,
  }));
  expect(dimensions.scroll).toBeLessThanOrEqual(dimensions.viewport + 1);
}

test.describe("published AOD production specification", () => {
  test.use({ colorScheme: "light" });

  test("desktop content, filters, downloads, and runtime", async ({ page, request }) => {
    await page.setViewportSize({ width: 1440, height: 1000 });
    const failures = captureRuntimeFailures(page);
    await page.goto(reportUrl, { waitUntil: "domcontentloaded" });

    await expect(page.locator("h1")).toHaveText("Low-cloud Sentinel-2 AOD retrieval");
    await expect(page.locator(".metric-strip")).toContainText("88.16%");
    await expect(page.locator(".metric-strip")).toContainText("134 / 152");
    await expect(page.locator("#surface")).toContainText("Single S2 L2A surface prior");
    await expect(page.locator("#rt")).toContainText("Radiative-transfer model");
    await expect(page.locator("#solver")).toContainText("68 float32 nodes");
    await expect(page.locator("#calibration")).toContainText("1500 trees");
    await expectNoDocumentOverflow(page);
    await page.screenshot({ path: "/tmp/aod-production-spec-desktop.png", fullPage: false });

    await page.locator("#feature-search").fill("tau");
    await expect(page.locator("#feature-count")).toHaveText("5 of 35 features");
    await page.locator("#feature-search").fill("");
    await page.locator("#feature-group").selectOption({ label: "CAMS" });
    await expect(page.locator("#feature-count")).toHaveText("8 of 35 features");
    await page.locator("#feature-group").selectOption("all");

    await page.locator("#source-search").fill("surface_driven");
    await expect(page.locator("[data-source-row]:visible")).toHaveCount(2);
    await page.locator("#source-search").fill("");

    const downloadLinks = await page.locator("a[href^='downloads/']").evaluateAll((links) =>
      [...new Set(links.map((link) => link.href))],
    );
    expect(downloadLinks.length).toBeGreaterThan(20);
    for (const url of downloadLinks) {
      const response = await request.fetch(url, { method: "HEAD" });
      expect(response.status(), url).toBeLessThan(400);
    }

    await page.locator("#validation").scrollIntoViewIfNeeded();
    await expect(page.locator("#validation")).toContainText("83/123");
    await page.screenshot({ path: "/tmp/aod-production-spec-validation-desktop.png", fullPage: false });
    expect(failures).toEqual([]);
  });

  test("mobile layout remains readable without document overflow", async ({ page }) => {
    await page.setViewportSize({ width: 390, height: 844 });
    const failures = captureRuntimeFailures(page);
    await page.goto(reportUrl, { waitUntil: "domcontentloaded" });
    await expect(page.locator("h1")).toBeVisible();
    await expect(page.locator(".metric-strip")).toContainText("134 / 152");
    await expectNoDocumentOverflow(page);
    await page.screenshot({ path: "/tmp/aod-production-spec-mobile.png", fullPage: false });

    await page.locator("#features").scrollIntoViewIfNeeded();
    await page.locator("#feature-search").fill("maiac");
    await expect(page.locator("#feature-count")).toHaveText("6 of 35 features");
    await expectNoDocumentOverflow(page);
    await page.screenshot({ path: "/tmp/aod-production-spec-features-mobile.png", fullPage: false });
    expect(failures).toEqual([]);
  });
});
