const { expect, test } = require("@playwright/test");

const dashboardUrl = `https://gws-access.jasmin.ac.uk/public/nceo_isp/siac_refactor/reports/aod-final-performance-dashboard-20260713/?qa=${Date.now()}`;

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
  page.on("response", (response) => {
    if (response.status() >= 400 && response.url().startsWith(dashboardUrl.split("?")[0])) {
      failures.push(`response: ${response.status()} ${response.url()}`);
    }
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

async function expectVisibleImagesLoaded(page, count) {
  await expect.poll(async () => page.locator(".gallery-item img").evaluateAll(
    (images, expected) => images.slice(0, expected).length === expected
      && images.slice(0, expected).every((image) => image.complete && image.naturalWidth > 0),
    count,
  )).toBe(true);
}

test.describe("published AOD dashboard", () => {
  test.use({ colorScheme: "light" });

  test("desktop views and interactions", async ({ page }) => {
    await page.setViewportSize({ width: 1440, height: 1000 });
    const failures = captureRuntimeFailures(page);
    await page.goto(dashboardUrl, { waitUntil: "domcontentloaded" });
    await expect(page.locator(".summary-band")).toContainText("134/152");
    await expect(page.locator(".summary-band")).toContainText("88.16%");
    await expectNoDocumentOverflow(page);
    await page.screenshot({ path: "/tmp/aod-dashboard-overview-desktop.png", fullPage: true });

    await page.locator('[data-view="cases"]').click();
    await expect(page.locator(".filter-count")).toHaveText("152 of 152 cases");
    await page.locator("#case-transition").selectOption("loss");
    await expect(page.locator(".filter-count")).toHaveText("5 of 152 cases");
    await expect(page.locator(".case-title-line .status-pill.loss")).toBeVisible();
    await page.locator('[data-tab="spatial"]').click();
    const spatialImage = page.locator(".evidence-figure img");
    await expect(spatialImage).toBeVisible();
    await spatialImage.evaluate((image) => image.complete || new Promise((resolve) => image.addEventListener("load", resolve, { once: true })));
    await page.screenshot({ path: "/tmp/aod-dashboard-case-desktop.png", fullPage: false });

    await page.locator('[data-view="evidence"]').click();
    await expect(page.locator(".gallery-count")).toHaveText("52 cases");
    await expect(page.locator(".gallery-item")).toHaveCount(52);
    await page.locator('[data-gallery-mode="diagnostic"]').click();
    await expect(page.locator('[data-gallery-mode="diagnostic"]')).toHaveClass(/active/);
    await expectVisibleImagesLoaded(page, 6);
    await page.screenshot({ path: "/tmp/aod-dashboard-gallery-desktop.png", fullPage: false });

    await page.locator('[data-view="robustness"]').click();
    await expect(page.locator("main")).toContainText("All variants ≥133");
    await expect(page.locator("main")).toContainText("134/152");
    await page.screenshot({ path: "/tmp/aod-dashboard-robustness-desktop.png", fullPage: true });

    await page.locator('[data-view="experiments"]').click();
    await expect(page.locator("svg[aria-label='Experiment within expected error rates']")).toContainText("Measured within EE");
    await expect(page.locator("svg[aria-label='Experiment within expected error rates']")).not.toContainText("Physical");
    await page.locator("#experiment-category").selectOption({ index: 1 });
    await expect(page.locator(".ledger-table tbody tr").first()).toBeVisible();

    await page.locator('[data-view="method"]').click();
    await expect(page.locator("main")).toContainText("One surface prior yes");
    await expect(page.locator("main")).toContainText("Case routing no");
    await expectNoDocumentOverflow(page);
    await page.screenshot({ path: "/tmp/aod-dashboard-method-desktop.png", fullPage: true });
    expect(failures).toEqual([]);
  });

  test("mobile views do not overflow the document", async ({ page }) => {
    await page.setViewportSize({ width: 375, height: 812 });
    const failures = captureRuntimeFailures(page);
    await page.goto(dashboardUrl, { waitUntil: "domcontentloaded" });
    await expect(page.locator(".summary-band")).toContainText("134/152");
    await expectNoDocumentOverflow(page);
    await page.screenshot({ path: "/tmp/aod-dashboard-overview-mobile.png", fullPage: false });

    await page.locator('[data-view="cases"]').click();
    await expect(page.locator(".filter-count")).toHaveText("152 of 152 cases");
    await expectNoDocumentOverflow(page);
    await page.screenshot({ path: "/tmp/aod-dashboard-case-mobile.png", fullPage: false });

    await page.locator('[data-view="evidence"]').click();
    await expect(page.locator(".gallery-count")).toHaveText("52 cases");
    await expectVisibleImagesLoaded(page, 1);
    await expectNoDocumentOverflow(page);
    await page.screenshot({ path: "/tmp/aod-dashboard-gallery-mobile.png", fullPage: false });
    expect(failures).toEqual([]);
  });
});
