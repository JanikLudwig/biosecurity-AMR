---
name: user-deterministic-no-llm-agents
description: "User wants the AMR decision path deterministic (no LLM agents) and target-presence to gate \"likely to work\""
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 9318c317-d76e-47bb-b922-2fc433f264f8
  modified: 2026-07-19T00:48:19.136Z
---

For Genome Firewall, the user wants the scientific decision path **deterministic and auditable —
no autonomous LLM agents** in it. The "drug-target agent" they originally described is a
**deterministic module** (pyrodigal ORFs → pyhmmer target search), not an LLM.

The core requirement (they added it to `# design notes.md`): a **"likely to work"** call must
**demonstrate the drug's molecular target is present in the genome** — proving the model isn't
just predicting success from an *absence* of resistance data. If the target can't be proven
present → **no-call**, never a false "works".

**Why:** biosecurity decision-support must be explainable/trustworthy; an LLM guessing whether a
target gene exists is exactly the failure mode to avoid.

**How to apply:** keep M2 (`gfw/targets/`) and M4 (`gfw/decide.py`) deterministic; any LLM is
optional narration over already-computed facts, off by default. See [[project-genome-firewall]].
