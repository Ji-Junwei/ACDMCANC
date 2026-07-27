# Project Scope

## Purpose

`ACDMCANC` is the public MATLAB reference repository for **Asynchronous Communication Distributed Multichannel Active Noise Control (AC-DMCANC)**.

The project studies how a distributed multichannel ANC system can reduce communication while retaining useful local adaptation under heterogeneous or bandwidth-limited network conditions.

This is a paper-specific Distributed ANC repository. It is not a general ANC algorithm collection.

## Authoritative implementation

The proposed asynchronous and synchronous communication methods are implemented in:

```text
AC_DMCANC.m
```

The primary experiment entry point is:

```text
AC_DCMANC_tst.m
```

The current filename uses `DCMANC` in the test-script name and `DMCANC` in the class name. Do not rename files solely for stylistic consistency without checking references and reproduction instructions.

## Current system model

The included experiment implements a specialized:

```text
1 reference × 6 secondary sources × 6 error sensors
```

with:

```text
Wc:   6 × (Lw + Lc - 1) in AC_DMCANC internal state
SecP: 6 × 6 × Ls
Dis:  6 × N
yc:   6 × N
C:    6 × 6 × Lc
```

The experiment loads:

```text
simulation path/PrimaryPath_1x6.mat
simulation path/SecondaryPath_6x6.mat
```

The current test configuration uses:

```text
Fs       = 16000 Hz
T        = 120 s
wLen     = 512
sLen     = 256
Numnode  = 6
cLen     = 33
muw      = 1e-6
muc      = 1e-5
alpha    = [800 800 800 800 800 800]
```

These are experiment parameters, not universal defaults for all DMCANC systems.

## Algorithm components

### Compensation-filter identification

`AC_DMCANC.CompensateSecP` identifies compensation filters for cross-secondary-path contributions.

For each off-diagonal path pair, white noise is passed through the cross path and a Filtered-x LMS system uses the corresponding diagonal secondary path. The final compensation filter is stored in:

```text
C(m, k, :)
```

The compensation-filter stage is part of the current weight-difference fusion implementation and must not be removed or replaced silently.

### Local WCFxLMS-style update

During local adaptation, node `k` uses its local error and diagonal secondary path. The implemented update has the form:

```text
W_k(n+1) = W_k(n)
           + mu * xf_k(n) * e_k(n)
           + alpha_k * mu * (Wcenter_k(n) - W_k(n))
```

where:

- `Wc` is the current local control-filter state;
- `Wcsubopt` is the center-point or constrained reference state;
- `alpha` controls attraction to that center point;
- `xf` is formed using the diagonal secondary path;
- `e` uses the full physical 6 × 6 secondary-path plant.

Changing the sign, scaling, center-point definition, filtered-reference construction, or role of `alpha` is an algorithm change.

### Communication trigger

The current implementation evaluates a residual-noise ratio every:

```text
t = 0.3 s
```

using the hard-coded `16000 * t` interval in `AC_DMCANC.m`.

For each node, communication is requested when the current interval metric no longer improves relative to the stored previous value. Communication events are recorded in:

```text
iscomm(node, sample)
```

Any change to the metric, comparison direction, averaging interval, initialization, hysteresis, or trigger synchronization must be documented explicitly.

### Weight-difference fusion

At a communication event, the implementation forms:

```text
Nabla = Wc - Wcsubopt
```

The requesting node combines its own weight difference with compensation-filtered weight differences from the other nodes, updates its center point, and resets its local controller to that updated center point.

This repository describes the exchange and fusion as compact weight-difference communication / Mixed Weight Difference. Changes to transmitted quantities, compensation filtering, node participation, aggregation order, or center-point replacement are communication-algorithm changes.

### Synchronous comparison

`AC_DMCANC.m` also contains a synchronous-communication comparison method:

```text
SC_DMCANC_166
```

The asynchronous and synchronous methods must remain distinguishable. A change to one method must not be copied to the other without checking the intended communication schedule and experiment comparison.

## Comparison implementations

| File | Role |
|---|---|
| `AC_DMCANC.m` | Proposed asynchronous method and synchronous communication comparison. |
| `DMANC_CompensateSP.m` | MGDFxLMS / gradient-transmission distributed comparison with compensation filters. |
| `McANC_FxLMS_SIMO.m` | Centralized 1 × 6 × 6 multichannel FxLMS comparison. |
| `AC_DCMANC_tst.m` | Reproduction script comparing centralized, MGDFxLMS, synchronous DMCANC, and asynchronous DMCANC. |

The comparison implementations are authoritative only for this repository's reproduction. They must not be treated as the general versions maintained in another repository.

## In scope

- Asynchronous communication distributed multichannel ANC.
- Event-triggered communication based on residual-noise performance.
- WCFxLMS-style local constrained adaptation.
- Compact weight-difference exchange and compensation-filtered fusion.
- Synchronous-communication comparison.
- Centralized and MGDFxLMS comparison baselines used by the experiment.
- The included six-node simulation paths and reproduction script.
- Communication-event, residual-error, and noise-reduction evaluation.
- MATLAB documentation and clearly scoped bug fixes.

## Out of scope

- General-purpose single-channel or centralized ANC libraries.
- Intermittent-communication coprocessor assistance implemented in `IC-DMCANC-CPA`.
- IMC feedback ANC.
- Output-constrained MOV-FxLMS.
- Private acoustic data, credentials, machine-specific absolute paths, or confidential results.
- Claims about arbitrary node counts without generalizing and validating the current explicitly expanded six-node code.
- Communication latency, packet loss, quantization, topology changes, or hardware-network protocols unless introduced as a documented extension.
- Python, C/C++, embedded, fixed-point, or firmware ports unless explicitly added as a separate implementation track.

## Core invariants

Unless the user explicitly requests an algorithm change, preserve:

- residual-error convention `e = Dis - y`;
- sample-by-sample adaptive updating;
- one reference, six controllers, and six error sensors in the current reproduction;
- full 6 × 6 physical secondary-path coupling when computing error signals;
- diagonal secondary-path filtering for each node's local update;
- the center-point attraction term controlled by `alpha`;
- compensation-filter generation before distributed fusion;
- independent node communication decisions in the asynchronous method;
- separate asynchronous and synchronous outputs;
- equal input signals, paths, step sizes, and comparison windows when comparing methods.

Do not describe the center-point penalty or event trigger alone as a formal proof of closed-loop stability.

## MATLAB dependencies

The current scripts use MATLAB functionality that may require relevant toolboxes, including:

- `dsp.FilteredXLMSFilter`;
- `awgn`;
- `fir1` and `filter`;
- `smooth` and plotting functions.

Toolbox-dependent execution results must be separated from static code inspection when the required MATLAB environment is unavailable.

## Validation expectations

Changes should use the strongest applicable checks:

1. inspect `AC_DMCANC.m` and `AC_DCMANC_tst.m` together;
2. identify whether the change affects asynchronous, synchronous, compensation-filter, centralized, or MGDFxLMS behaviour;
3. verify array sizes for `Wc`, `Wcsubopt`, `Nabla`, `SecP`, `C`, `Dis`, `yc`, `e`, and `iscomm`;
4. check the trigger interval and sample-rate dependency;
5. verify finite controller weights, errors, residual metrics, and compensation filters;
6. verify that non-communication operation continues local adaptation;
7. verify that communication updates only the intended requesting node in the asynchronous method;
8. compare communication counts and residual-noise curves against synchronous and centralized baselines;
9. rerun `AC_DCMANC_tst.m` and regenerate plots when MATLAB is available;
10. document changes to equations, trigger logic, parameters, simulation paths, or figure generation;
11. report every modified, added, renamed, or deleted file.

When MATLAB execution is unavailable, clearly separate code-level validation from reproduction results that still require local MATLAB and toolbox execution.
