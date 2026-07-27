# GPT Context

## Repository identity

This is the public `ANC-Research/ACDMCANC` repository.

It contains a paper-specific MATLAB implementation of asynchronous-communication distributed multichannel active noise control. Do not treat it as a general ANC library, the coprocessor-assisted intermittent-communication project, IMC feedback ANC, or output-constrained ANC.

Because this repository is public, do not add private repository names, unpublished methods, private data descriptions, internal results, credentials, or confidential future plans.

## Required reading before making changes

Before modifying code, read the latest version of:

1. `PROJECT_SCOPE.md`;
2. the root `README.md`;
3. `AC_DMCANC.m`;
4. `AC_DCMANC_tst.m`;
5. the directly affected comparison implementation;
6. the paper definition or existing implementation evidence relevant to any algorithm change.

Always use the latest requested branch. Do not rely only on previous conversations, cached code, or similarly named distributed ANC implementations from another repository.

## Implementation authority

- `AC_DMCANC.m` is the authoritative implementation for the asynchronous and synchronous communication methods.
- `AC_DCMANC_tst.m` is the authoritative reproduction script and parameter set for the included experiment.
- `DMANC_CompensateSP.m` is the MGDFxLMS / gradient-transmission comparison used by this project.
- `McANC_FxLMS_SIMO.m` is the centralized 1 × 6 × 6 comparison used by this project.
- The included `.mat` path files define the current public simulation plant.
- Keep unrelated files and comparison methods read-only unless the requested change affects them.

## Task isolation

At the start of each task, identify:

- whether the target is asynchronous communication, synchronous communication, local WCFxLMS-style adaptation, trigger logic, compensation-filter identification, weight-difference fusion, a comparison method, or plotting;
- the exact method and entry point;
- the controlling node count, sample rate, filter lengths, step sizes, `alpha`, trigger interval, paths, and input signals;
- whether the task changes equations, communication behaviour, parameters, array generality, data dependencies, validation, or only presentation;
- which residual-error, communication-event, and comparison outputs are expected to change.

Do not assume that asynchronous and synchronous methods should receive identical code changes.

## Current structural assumptions

The current reproduction is explicitly implemented for:

```text
1 reference × 6 secondary sources × 6 error sensors
```

Important shapes include:

```text
SecP:      (6, 6, Ls)
Dis:       (6, N)
yc:        (6, N)
C:         (6, 6, Lc)
iscomm:    (6, N)
```

Many controller and error calculations are expanded node by node. Do not claim support for arbitrary `K` merely because a constructor accepts `node_num`. Generalization requires replacing and validating every hard-coded six-node operation.

## Core algorithm invariants

Preserve these unless the user explicitly requests an algorithm change:

```text
e_k(n) = Dis_k(n) - y_k(n)
```

```text
W_k(n+1) = W_k(n)
           + mu * xf_k(n) * e_k(n)
           + alpha_k * mu * (Wcenter_k(n) - W_k(n))
```

Also preserve:

- sample-by-sample recursive adaptation;
- full cross-secondary-path coupling in physical error generation;
- diagonal secondary-path filtering for the local update;
- separate local controller and center-point states;
- weight-difference calculation `Nabla = Wc - Wcsubopt`;
- compensation filtering of other-node differences before fusion;
- independent node communication events in the asynchronous method;
- distinct synchronous and asynchronous result arrays;
- equal experiment conditions when comparing methods.

Do not replace the center-point attraction term, compensation-filtered fusion, or event trigger with a conventional centralized FxLMS update without explicitly defining a new method.

## Communication-trigger changes

When changing the trigger:

1. state the current and proposed residual metric;
2. identify the evaluation interval and its sample-rate dependency;
3. explain the comparison direction and initialization;
4. state whether decisions are independent per node or synchronized;
5. verify behaviour when residual metrics are equal, non-finite, or noisy;
6. report communication counts per node and total communication events;
7. compare residual-noise performance against the unchanged trigger;
8. update documentation and plots that interpret `iscomm`.

Do not describe fewer trigger events as an improvement without reporting the associated convergence and final noise-reduction trade-off.

## Local constraint and fusion changes

When changing `alpha`, `Wcsubopt`, `Nabla`, compensation filters, or fusion:

- map each mathematical quantity to the exact MATLAB property and update lines;
- distinguish local adaptation from communication-time aggregation;
- verify lengths `Lw`, `Lc`, and `Lw + Lc - 1`;
- check the selected tail segment after compensation filtering;
- verify which node state is updated at an asynchronous event;
- preserve or explicitly redefine the center-point reset step;
- compare controller norms, finite values, communication counts, and residual error;
- update both asynchronous and synchronous paths only when mathematically intended.

A weight constraint or center-point penalty should not be described as a general proof that every possible networked ANC configuration is stable.

## Experiment and data rules

- Preserve the included simulation-path assets and relative paths.
- Do not commit private acoustic recordings, private measured paths, or machine-specific absolute paths.
- Treat generated plots and saved workspace files as experiment artifacts, not source-code authority.
- Record random seeds when reproducibility is added or changed.
- When adding a new network or plant case, document node count, path dimensions, sample rate, duration, filter lengths, communication assumptions, and all algorithm parameters.
- Do not overwrite existing public `.mat` assets with unrelated data.

## MATLAB change discipline

- Make the smallest coherent change that satisfies the request.
- Avoid broad stylistic rewrites of expanded six-node reference code during an algorithm fix.
- Preserve class and method interfaces unless an interface change is required.
- Check row/column orientation and `reshape` order explicitly.
- Keep algorithm changes separate from plotting changes where practical.
- Do not rename the `AC_DCMANC_tst.m` file merely to correct spelling without checking references.
- Report every modified, added, renamed, or deleted file.
- Add a concise entry to `CHANGELOG.md` for significant user-visible changes.
- Update the root README when setup, algorithms, parameters, dependencies, or reproduction instructions change.

## Cross-repository restrictions

Do not import code, parameters, or conclusions from other ANC repositories unless the user explicitly requests a documented dependency or comparison.

In particular, do not merge this implementation with `IC-DMCANC-CPA`. The two repositories represent different communication architectures and must remain independently traceable.

Similarly named WCFxLMS, MWD, FxLMS, compensation-filter, or communication routines in different repositories are independent implementations until their equations and conventions have been compared.

## Verification

Use the strongest verification available for the change.

For static inspection, check:

- MATLAB syntax, class names, method names, and call signatures;
- path-file names and dimensions;
- controller, gradient, compensation-filter, and history-array sizes;
- initialization of `NR_pre`, `NR_now`, `Wc`, `Wcsubopt`, and `Nabla`;
- sample-rate dependence of the trigger interval;
- finite values in control filters, errors, compensation filters, and residual metrics;
- asynchronous updates affecting only intended nodes;
- comparison methods receiving the same disturbance, reference, paths, and step sizes.

For MATLAB execution, run:

```matlab
AC_DCMANC_tst
```

As applicable, inspect:

- residual-noise curves for each node and their mean;
- centralized, MGDFxLMS, synchronous, and asynchronous comparisons;
- communication-event traces and counts;
- local and center-point controller norms;
- compensation-filter identification errors;
- convergence and final noise-reduction trade-offs;
- regenerated figures and saved variables.

MATLAB toolbox availability and the 120-second sample-by-sample simulation may prevent execution in the current environment. Clearly separate code inspection from MATLAB results that must be produced locally.

## Response requirements

For every completed modification, provide:

- files changed;
- the purpose of each change;
- the affected method and communication phase;
- original and updated equation or trigger mapping, when applicable;
- parameter, path, and array-dimension impact;
- verification performed and results;
- MATLAB or local-runtime checks still required.
