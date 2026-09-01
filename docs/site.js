document.documentElement.classList.add("js");

const themeSelect = document.querySelector("[data-theme-select]");
const colorScheme = window.matchMedia("(prefers-color-scheme: dark)");

const applyTheme = mode => {
  const resolved = mode === "system" ? (colorScheme.matches ? "dark" : "light") : mode;
  document.documentElement.dataset.theme = resolved;
  document.documentElement.dataset.themeMode = mode;
  document.querySelector('meta[name="theme-color"]')?.setAttribute(
    "content",
    resolved === "dark" ? "#071a1c" : "#f2f1e9"
  );
};

const savedTheme = localStorage.getItem("gpu-eqtl-theme") || "system";
if (themeSelect) themeSelect.value = savedTheme;
applyTheme(savedTheme);

themeSelect?.addEventListener("change", () => {
  const mode = themeSelect.value;
  localStorage.setItem("gpu-eqtl-theme", mode);
  applyTheme(mode);
});

colorScheme.addEventListener("change", () => {
  if ((localStorage.getItem("gpu-eqtl-theme") || "system") === "system") applyTheme("system");
});

const analysisContent = {
  eqtl: {
    kicker: "Single-variant association",
    title: "Test every variant against every expression trait.",
    description: "Use FP64 by default, let the backend and tile sizes resolve automatically, and report only associations crossing the requested threshold.",
    unit: "variant × phenotype",
    input: "none",
    result: "effect, t, log10 p",
    command: `java -Xmx12g --enable-native-access=ALL-UNNAMED \`
  -jar target\\gpu-eqtl-2.0.0-SNAPSHOT-all.jar \`
  --genotype cohort.vcf.gz --genotype-format vcf \`
  --expression expression.csv --covariates covariates.csv \`
  --genotype-id-column participant_id \`
  --expression-id-column sample_id \`
  --fixed-covariates Age,Sex,PC1,PC2 \`
  --analysis eqtl --precision fp64 \`
  --threshold pval 1e-4 --output results.csv`
  },
  burden: {
    kicker: "Collapsed rare-variant signal",
    title: "Test weighted burdens in native sliding windows.",
    description: "Create chromosome-local windows directly from canonical variant coordinates. Window size and stride replace a regions file for regular grids.",
    unit: "window × phenotype",
    input: "coordinates in variant IDs",
    result: "burden effect, t, log10 p",
    command: `java -Xmx12g --enable-native-access=ALL-UNNAMED \`
  -jar target\\gpu-eqtl-2.0.0-SNAPSHOT-all.jar \`
  --genotype cohort.vcf.gz --genotype-format vcf \`
  --expression expression.csv --covariates covariates.csv \`
  --fixed-covariates Age,Sex,PC1,PC2 \`
  --analysis burden --window-size 50000 \`
  --window-stride 10000 --min-mac 0 \`
  --precision fp64 --output burden-results.csv`
  },
  skat: {
    kicker: "Variance-component set test",
    title: "Detect mixed-direction effects with SKAT.",
    description: "Residualize declared-effect dosages, calculate weighted score covariance in FP64, and retain the named p-value method in the audit output.",
    unit: "set × phenotype",
    input: "windows or variant-set TSV",
    result: "Q statistic, p-value method",
    command: `java -Xmx12g --enable-native-access=ALL-UNNAMED \`
  -jar target\\gpu-eqtl-2.0.0-SNAPSHOT-all.jar \`
  --genotype cohort.vcf.gz --genotype-format vcf \`
  --expression expression.csv --covariates covariates.csv \`
  --fixed-covariates Age,Sex,PC1,PC2 \`
  --analysis skat --window-size 50000 \`
  --window-stride 10000 --min-mac 0 \`
  --precision fp64 --output skat-results.csv`
  },
  "skat-o": {
    kicker: "Adaptive omnibus set test",
    title: "Balance burden and SKAT signals with SKAT-O.",
    description: "Evaluate the configured rho grid and report a seeded, correlated-null adjusted omnibus p-value—not an unadjusted minimum component p-value.",
    unit: "set × phenotype",
    input: "windows or variant-set TSV",
    result: "adjusted omnibus p-value",
    command: `java -Xmx12g --enable-native-access=ALL-UNNAMED \`
  -jar target\\gpu-eqtl-2.0.0-SNAPSHOT-all.jar \`
  --genotype cohort.vcf.gz --genotype-format vcf \`
  --expression expression.csv --covariates covariates.csv \`
  --fixed-covariates Age,Sex,PC1,PC2 \`
  --analysis skat-o --window-size 50000 \`
  --window-stride 10000 --skat-o-simulations 10000 \`
  --precision fp64 --output skat-o-results.csv`
  }
};

const header = document.querySelector("[data-header]");
const menuButton = document.querySelector("[data-menu-button]");
const menu = document.querySelector("[data-menu]");

const updateHeader = () => header?.classList.toggle("scrolled", window.scrollY > 18);
updateHeader();
window.addEventListener("scroll", updateHeader, { passive: true });

menuButton?.addEventListener("click", () => {
  const open = menuButton.getAttribute("aria-expanded") === "true";
  menuButton.setAttribute("aria-expanded", String(!open));
  menu?.classList.toggle("open", !open);
});

menu?.querySelectorAll("a").forEach(link => link.addEventListener("click", () => {
  menuButton?.setAttribute("aria-expanded", "false");
  menu?.classList.remove("open");
}));

const tabs = document.querySelectorAll("[data-analysis]");
const command = document.querySelector("[data-command]");
const fields = {
  kicker: document.querySelector("[data-analysis-kicker]"),
  title: document.querySelector("[data-analysis-title]"),
  description: document.querySelector("[data-analysis-description]"),
  unit: document.querySelector("[data-analysis-unit]"),
  input: document.querySelector("[data-analysis-input]"),
  result: document.querySelector("[data-analysis-result]")
};

tabs.forEach(tab => tab.addEventListener("click", () => {
  const key = tab.dataset.analysis;
  const content = analysisContent[key];
  if (!content) return;

  tabs.forEach(candidate => candidate.setAttribute("aria-selected", String(candidate === tab)));
  Object.entries(fields).forEach(([name, node]) => {
    if (node) node.textContent = content[name];
  });
  if (command) command.textContent = content.command;
}));

const copyButton = document.querySelector("[data-copy-button]");
copyButton?.addEventListener("click", async () => {
  if (!command) return;
  const label = copyButton.querySelector("span");
  try {
    await navigator.clipboard.writeText(command.textContent);
    if (label) label.textContent = "Copied";
  } catch {
    if (label) label.textContent = "Select text";
  }
  window.setTimeout(() => { if (label) label.textContent = "Copy"; }, 1600);
});

const revealItems = document.querySelectorAll(".reveal");
if ("IntersectionObserver" in window && !window.matchMedia("(prefers-reduced-motion: reduce)").matches) {
  const observer = new IntersectionObserver(entries => {
    entries.forEach(entry => {
      if (!entry.isIntersecting) return;
      entry.target.classList.add("visible");
      observer.unobserve(entry.target);
    });
  }, { threshold: .12 });
  revealItems.forEach(item => observer.observe(item));
} else {
  revealItems.forEach(item => item.classList.add("visible"));
}
