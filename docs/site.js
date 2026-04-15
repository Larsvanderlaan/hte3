const hteMethodData = {
  cate: {
    title: "Conditional Average Treatment Effects",
    subtitle: "Estimate how treatment effects vary across effect modifiers.",
    methods: {
      dr: {
        label: "DR-learner",
        badge: "Binary/categorical treatment",
        title: "Doubly robust pseudo-outcome regression for CATE",
        summary:
          "Builds a doubly robust pseudo-outcome from the propensity and outcome models, then regresses that signal on the effect modifiers.",
        good: [
          "You want a strong default CATE workflow.",
          "Treatment is binary or categorical and you can estimate both nuisance models.",
          "You want a clear practitioner-facing path from task construction to prediction."
        ],
        tradeoffs: [
          "Requires both propensity and outcome nuisance estimates.",
          "Not currently the production path for continuous treatment.",
          "More moving pieces than a simple T-learner baseline."
        ],
        meta: {
          Target: "CATE",
          Supports: "Binary and categorical treatment",
          Output: "Difference in conditional means",
          Call: 'fit_cate(task, method = "dr", ...)'
        },
        note:
          "Uses both propensity and outcome nuisance estimates to form a doubly robust pseudo-outcome."
      },
      r: {
        label: "R-learner",
        badge: "Continuous treatment supported",
        title: "Residual-on-residual learning for effect heterogeneity",
        summary:
          "Residualizes the outcome and treatment first, then learns the heterogeneity surface from those residuals. In the production API, this is the CATE method that extends to continuous treatment when treatment is numerically coded, using the partially linear effect-model view with an A times tau(X) term.",
        good: [
          "Treatment is continuous or numerically coded.",
          "You want a CATE surface without choosing treated and control labels.",
          "You prefer a residualization-based view of effect learning."
        ],
        tradeoffs: [
          "Requires numeric treatment coding.",
          "For continuous treatment, relies on the partially linear A times tau(X) effect-model assumption.",
          "Can become unstable when treatment residuals are near zero.",
          "Does not use the same contrast-style interface as DR, EP, or T."
        ],
        meta: {
          Target: "CATE",
          Supports: "Binary or continuous treatment",
          Output: "Conditional effect surface",
          Call: 'fit_cate(task, method = "r", ...)'
        },
        note:
          "Residualizes treatment and outcome before learning the effect surface."
      },
      ep: {
        label: "EP-learner",
        badge: "Flexible sieve debiasing",
        title: "Efficient plug-in learner for flexible CATE surfaces",
        summary:
          "Starts from arm-specific outcome regressions and debiases them with a sieve basis over the effect modifiers. In practice, this combines the intuition of T-learning with a debiasing step that can recover richer surfaces than a plain baseline model.",
        good: [
          "You want a flexible alternative to DR-learner.",
          "You are willing to tune or cross-validate sieve complexity.",
          "You expect smooth but nonlinear structure in the modifier space."
        ],
        tradeoffs: [
          "More compute and tuning than DR or T.",
          "Sieve settings matter more in moderate or high dimension.",
          "Best used when you are comfortable moving one level deeper into the modeling stack."
        ],
        meta: {
          Target: "CATE",
          Supports: "Binary and categorical treatment",
          Output: "Debiased conditional mean contrast",
          Call: 'fit_cate(task, method = "ep", ...)'
        },
        note:
          "Combines arm-specific outcome regression with sieve-based debiasing."
      },
      t: {
        label: "T-learner",
        badge: "Simple baseline",
        title: "Separate treatment-arm models, then subtract",
        summary:
          "Fits outcome models by treatment arm and subtracts the fitted means to get the effect surface. This is a common benchmark for comparison against more robust methods.",
        good: [
          "You want a simple baseline for comparison.",
          "You expect separate arm-specific outcome models to capture the signal.",
          "You want the most direct mapping from outcome modeling to CATE predictions."
        ],
        tradeoffs: [
          "Not doubly robust.",
          "Can carry arm-specific modeling bias directly into the effect surface.",
          "Often better used as a benchmark than as the only final analysis."
        ],
        meta: {
          Target: "CATE",
          Supports: "Binary and categorical treatment",
          Output: "Difference of arm-specific fitted means",
          Call: 'fit_cate(task, method = "t", ...)'
        },
        note:
          "Fits outcome models separately by treatment arm and subtracts them."
      }
    },
    defaultMethod: "r"
  },
  crr: {
    title: "Conditional Risk Ratios",
    subtitle: "Estimate relative treatment effects when the outcome is nonnegative.",
    methods: {
      ep: {
        label: "EP-learner",
        badge: "Sieve debiasing",
        title: "Flexible debiasing for conditional risk-ratio estimation",
        summary:
          "Starts from arm-specific outcome regressions, debiases them with a sieve basis, and returns a conditional risk-ratio surface.",
        good: [
          "Outcome is binary or otherwise nonnegative.",
          "You want the main production-path CRR method.",
          "You need a flexible ratio surface rather than a simple weighted baseline."
        ],
        tradeoffs: [
          "Needs nonnegative outcomes.",
          "More tuning and compute than IPW or T.",
          "Sieve complexity can matter in larger modifier spaces."
        ],
        meta: {
          Target: "CRR",
          Supports: "Binary and categorical treatment",
          Output: "Conditional risk ratio",
          Call: 'fit_crr(task, method = "ep", ...)'
        },
        note:
          "Uses sieve-based debiasing on top of arm-specific outcome regression."
      },
      ipw: {
        label: "IPW",
        badge: "Weighting baseline",
        title: "Inverse-probability weighted CRR learning",
        summary:
          "Uses inverse-probability weights to learn how treated-versus-control risk changes across modifier values. This is a benchmark that places more weight on the propensity model than on outcome regression.",
        good: [
          "You want a weighting-based CRR benchmark.",
          "Outcome is nonnegative and you trust the propensity model more than the outcome model.",
          "You want something simpler than EP while staying in the CRR framework."
        ],
        tradeoffs: [
          "Sensitive to unstable propensities and limited overlap.",
          "Usually noisier than EP when outcome modeling is informative.",
          "Only applies to binary or categorical treatment settings."
        ],
        meta: {
          Target: "CRR",
          Supports: "Binary and categorical treatment",
          Output: "Weighted conditional risk ratio",
          Call: 'fit_crr(task, method = "ipw", ...)'
        },
        note:
          "Places more emphasis on weighting than on outcome regression."
      },
      t: {
        label: "T-learner",
        badge: "Simple CRR baseline",
        title: "Arm-specific outcome models on the risk-ratio scale",
        summary:
          "Fits separate outcome models by treatment arm and then takes the log or ratio transformation needed for conditional risk-ratio estimation. This is best treated as a benchmark rather than the primary analysis.",
        good: [
          "You want an easy-to-explain benchmark.",
          "Most of your signal is expected to come from the outcome model.",
          "You are comparing simple and flexible CRR pipelines side by side."
        ],
        tradeoffs: [
          "Not doubly robust.",
          "Needs nonnegative outcomes.",
          "Can inherit bias from weak arm-specific outcome models."
        ],
        meta: {
          Target: "CRR",
          Supports: "Binary and categorical treatment",
          Output: "Ratio of arm-specific fitted risks",
          Call: 'fit_crr(task, method = "t", ...)'
        },
        note:
          "Computes the ratio from separate arm-specific outcome models."
      }
    },
    defaultMethod: "ep"
  }
};

function renderMethodExplorer() {
  const explorer = document.querySelector("[data-module='hte-method-explorer']");
  if (!explorer) {
    return;
  }

  const estimandButtons = explorer.querySelectorAll("[data-estimand]");
  const methodButtonRow = explorer.querySelector("[data-method-buttons]");
  const titleNode = explorer.querySelector("[data-method-title]");
  const badgeNode = explorer.querySelector("[data-method-badge]");
  const summaryNode = explorer.querySelector("[data-method-summary]");
  const goodNode = explorer.querySelector("[data-method-good]");
  const tradeoffNode = explorer.querySelector("[data-method-tradeoffs]");
  const noteNode = explorer.querySelector("[data-method-note]");
  const metaNode = explorer.querySelector("[data-method-meta]");
  const codeNode = explorer.querySelector("[data-method-code]");
  const estimandTitleNode = explorer.querySelector("[data-estimand-title]");
  const estimandSubtitleNode = explorer.querySelector("[data-estimand-subtitle]");

  let activeEstimand = "cate";
  let activeMethod = hteMethodData[activeEstimand].defaultMethod;

  function fillList(node, items) {
    node.innerHTML = "";
    items.forEach((item) => {
      const li = document.createElement("li");
      li.textContent = item;
      node.appendChild(li);
    });
  }

  function fillMeta(node, meta) {
    node.innerHTML = "";
    Object.entries(meta).forEach(([label, value]) => {
      const card = document.createElement("div");
      card.className = "metric-card";
      const labelNode = document.createElement("span");
      labelNode.textContent = label;
      const valueNode = document.createElement("strong");
      valueNode.textContent = value;
      card.appendChild(labelNode);
      card.appendChild(valueNode);
      node.appendChild(card);
    });
  }

  function renderMethodButtons() {
    const estimand = hteMethodData[activeEstimand];
    methodButtonRow.innerHTML = "";

    Object.entries(estimand.methods).forEach(([key, method]) => {
      const button = document.createElement("button");
      button.type = "button";
      button.className = `tab-button${key === activeMethod ? " is-active" : ""}`;
      button.textContent = method.label;
      button.setAttribute("data-method", key);
      button.addEventListener("click", () => {
        activeMethod = key;
        renderExplorer();
      });
      methodButtonRow.appendChild(button);
    });
  }

  function renderExplorer() {
    const estimand = hteMethodData[activeEstimand];
    const method = estimand.methods[activeMethod];

    estimandButtons.forEach((button) => {
      button.classList.toggle("is-active", button.dataset.estimand === activeEstimand);
    });

    renderMethodButtons();

    estimandTitleNode.textContent = estimand.title;
    estimandSubtitleNode.textContent = estimand.subtitle;
    titleNode.textContent = method.title;
    badgeNode.textContent = method.badge;
    summaryNode.textContent = method.summary;
    noteNode.textContent = method.note;
    codeNode.textContent = method.meta.Call;

    fillList(goodNode, method.good);
    fillList(tradeoffNode, method.tradeoffs);
    fillMeta(metaNode, method.meta);
  }

  estimandButtons.forEach((button) => {
    button.addEventListener("click", () => {
      activeEstimand = button.dataset.estimand;
      activeMethod = hteMethodData[activeEstimand].defaultMethod;
      renderExplorer();
    });
  });

  renderExplorer();
}

function revealOnScroll() {
  const items = document.querySelectorAll("[data-reveal]");
  if (!items.length || typeof IntersectionObserver === "undefined") {
    items.forEach((item) => item.classList.add("is-revealed"));
    return;
  }

  const observer = new IntersectionObserver(
    (entries) => {
      entries.forEach((entry) => {
        if (entry.isIntersecting) {
          entry.target.classList.add("is-revealed");
          observer.unobserve(entry.target);
        }
      });
    },
    { threshold: 0.12 }
  );

  items.forEach((item) => observer.observe(item));
}

document.addEventListener("DOMContentLoaded", () => {
  renderMethodExplorer();
  revealOnScroll();
});
