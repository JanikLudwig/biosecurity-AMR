import type { AnalysisReport } from "@/lib/genome-firewall-api";

export type DemoDecisionTone = "work" | "fail" | "no-call";

export interface DemoGenome {
  id: string;
  label: string;
  title: string;
  summary: string;
  decisions: Array<{
    drug: string;
    result: string;
    tone: DemoDecisionTone;
  }>;
}

export const DEMO_GENOMES: DemoGenome[] = [
  {
    id: "1280.51926",
    label: "Mixed response",
    title: "Directional calls with a target safeguard",
    summary:
      "Cefoxitin and ciprofloxacin pass every decision gate. Erythromycin is withheld because its molecular target could not be verified.",
    decisions: [
      { drug: "Cefoxitin", result: "Likely to work", tone: "work" },
      { drug: "Ciprofloxacin", result: "Likely to work", tone: "work" },
      { drug: "Erythromycin", result: "No-call", tone: "no-call" },
    ],
  },
  {
    id: "1280.51872",
    label: "Resistance evidence",
    title: "Known ciprofloxacin resistance signal",
    summary:
      "A high calibrated resistance probability and known resistance evidence produce a likely-to-fail ciprofloxacin result.",
    decisions: [
      { drug: "Cefoxitin", result: "No-call", tone: "no-call" },
      { drug: "Ciprofloxacin", result: "Likely to fail", tone: "fail" },
      { drug: "Erythromycin", result: "No-call", tone: "no-call" },
    ],
  },
  {
    id: "1280.51741",
    label: "Safety gates",
    title: "Out-of-distribution feature profile",
    summary:
      "All three directional model signals are withheld because this genome falls outside the training feature distribution.",
    decisions: [
      { drug: "Cefoxitin", result: "No-call", tone: "no-call" },
      { drug: "Ciprofloxacin", result: "No-call", tone: "no-call" },
      { drug: "Erythromycin", result: "No-call", tone: "no-call" },
    ],
  },
];

export async function loadDemoAnalysis(id: string, signal?: AbortSignal): Promise<AnalysisReport> {
  const response = await fetch(`/demo-results/${encodeURIComponent(id)}.report.json`, { signal });
  if (!response.ok) throw new Error(`Saved demo analysis is unavailable (${response.status})`);
  const report = (await response.json()) as AnalysisReport;
  return { ...report, analysis_id: `demo-${id}` };
}
