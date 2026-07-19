import { createFileRoute, Link } from "@tanstack/react-router";
import { useEffect, useMemo, useState } from "react";
import { z } from "zod";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { AlertTriangle, ArrowLeft, Dna, Loader2 } from "lucide-react";
import {
  GenomeFirewallApiError,
  getGenomeReport,
  getPerformanceMetrics,
  type AntibioticDecision,
  type GenomeReport,
  type PerformanceMetric,
  type PredictionCall,
  type TargetEvidence,
} from "@/lib/genome-firewall-api";

export const Route = createFileRoute("/dashboard")({
  validateSearch: z.object({ genomeId: z.string().optional() }),
  component: DashboardScreen,
  head: () => ({ meta: [{ title: "Analysis - GyraseX" }] }),
});

type Status = "work" | "fail" | "nocall";

function toStatus(call: PredictionCall): Status {
  if (call === "likely to work") return "work";
  if (call === "likely to fail") return "fail";
  return "nocall";
}

function statusStyles(status: Status): { bg: string; fg: string; label: string } {
  if (status === "work") return { bg: "var(--success)", fg: "#fff", label: "Likely to work" };
  if (status === "fail") return { bg: "var(--danger)", fg: "#fff", label: "Likely to fail" };
  return { bg: "var(--warning)", fg: "#1a1a1a", label: "No-call" };
}

function percent(value: number | undefined | null, digits = 1): string {
  if (value == null || !Number.isFinite(value)) return "Not available";
  return `${(value * 100).toFixed(digits)}%`;
}

function decimal(value: number | undefined | null, digits = 3): string {
  if (value == null || !Number.isFinite(value)) return "Not available";
  return value.toFixed(digits);
}

function evidenceLabel(category: string): string {
  const labels: Record<string, string> = {
    known_resistance_determinant: "Known resistance determinant",
    statistical_association: "Statistical association",
    no_known_resistance_signal: "No known resistance signal",
  };
  return labels[category] || category.replaceAll("_", " ");
}

function DashboardScreen() {
  const { genomeId } = Route.useSearch();
  const [report, setReport] = useState<GenomeReport | null>(null);
  const [metrics, setMetrics] = useState<PerformanceMetric[]>([]);
  const [selectedId, setSelectedId] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState(Boolean(genomeId));

  useEffect(() => {
    if (!genomeId) {
      setLoading(false);
      setReport(null);
      setMetrics([]);
      return;
    }

    const controller = new AbortController();
    setLoading(true);
    setError(null);
    Promise.all([
      getGenomeReport(genomeId, controller.signal),
      getPerformanceMetrics(controller.signal),
    ])
      .then(([nextReport, nextMetrics]) => {
        setReport(nextReport);
        setMetrics(nextMetrics);
        setSelectedId(nextReport.decisions[0]?.drug ?? null);
      })
      .catch((reason: unknown) => {
        if ((reason as DOMException)?.name !== "AbortError") {
          setError(reason instanceof GenomeFirewallApiError ? reason.message : "The analysis report could not be loaded.");
        }
      })
      .finally(() => {
        if (!controller.signal.aborted) setLoading(false);
      });

    return () => controller.abort();
  }, [genomeId]);

  const metricsByDrug = useMemo(() => new Map(metrics.map((metric) => [metric.drug, metric])), [metrics]);
  const selected = report?.decisions.find((decision) => decision.drug === selectedId) || report?.decisions[0];
  const anyPositive = Boolean(report?.decisions.some((decision) => decision.call === "likely to work"));

  if (!genomeId) {
    return <DashboardMessage title="No sample selected" message="Choose a registered genome before opening its analysis." />;
  }

  if (loading) {
    return (
      <div className="flex min-h-screen items-center justify-center bg-background">
        <div className="flex items-center gap-3 font-mono-tabular text-xs tracking-wider uppercase text-muted-foreground">
          <Loader2 className="h-4 w-4 animate-spin" /> Loading report
        </div>
      </div>
    );
  }

  if (error || !report || !selected) {
    return <DashboardMessage title="Analysis unavailable" message={error || "No antibiotic decisions were returned for this sample."} />;
  }

  return (
    <div className="flex min-h-screen flex-col bg-background text-foreground">
      <header className="flex items-center justify-between border-b border-border bg-card px-6 py-3">
        <div className="flex items-center gap-6">
          <Link to="/" className="flex items-center gap-2 text-sm text-muted-foreground hover:text-foreground">
            <ArrowLeft className="h-4 w-4" />
            <span className="font-mono-tabular text-[11px] tracking-wider uppercase">Back</span>
          </Link>
          <div className="flex items-center gap-2">
            <Dna className="h-4 w-4" />
            <span className="font-mono-tabular text-[11px] font-semibold tracking-[0.2em] uppercase">
              GyraseX / Analysis
            </span>
          </div>
        </div>
        <div className="flex items-center gap-6 font-mono-tabular text-[11px] tracking-wider uppercase text-muted-foreground">
          <div>Sample <span className="text-foreground">{report.genome_id}</span></div>
          <div>Species <span className="italic text-foreground">{report.species}</span></div>
          <div>Drugs <span className="text-foreground">{report.decisions.length}</span></div>
        </div>
      </header>

      {report.features_synthetic && (
        <Alert className="rounded-none border-x-0 border-t-0 border-border bg-muted/50 px-6">
          <AlertTriangle className="h-4 w-4" style={{ color: "var(--warning)" }} />
          <AlertTitle className="font-mono-tabular text-[11px] tracking-wider uppercase">Synthetic M1 features</AlertTitle>
          <AlertDescription>Predictions and validation metrics are illustrative until real AMRFinderPlus features replace the placeholder matrix.</AlertDescription>
        </Alert>
      )}

      <div className="grid flex-1 grid-cols-[340px_minmax(0,1fr)]">
        <aside className="flex h-[calc(100vh-52px)] flex-col border-r border-border bg-card">
          <div className="border-b border-border px-5 py-4">
            <div className="font-mono-tabular text-[10px] tracking-[0.2em] uppercase text-muted-foreground">Directory</div>
            <h2 className="mt-1 text-base font-semibold tracking-tight">List of Antibiotics</h2>
            <div className="mt-2 font-mono-tabular text-[10px] text-muted-foreground">
              {report.summary["likely to work"] || 0} work / {report.summary["likely to fail"] || 0} fail / {report.summary["no-call"] || 0} no-call
            </div>
          </div>
          <div className="flex-1 overflow-y-auto">
            {report.decisions.map((decision) => {
              const status = toStatus(decision.call);
              const style = statusStyles(status);
              const active = decision.drug === selected.drug;
              return (
                <button
                  key={decision.drug}
                  onClick={() => setSelectedId(decision.drug)}
                  className={`flex w-full items-center gap-3 border-b border-border px-5 py-4 text-left transition-colors hover:bg-accent ${active ? "bg-accent" : ""}`}
                >
                  <div className="h-8 w-[3px] shrink-0 rounded-sm" style={{ background: active ? style.bg : "transparent" }} />
                  <div className="min-w-0 flex-1">
                    <div className="truncate text-sm font-medium">{decision.drug_display}</div>
                    <div className="mt-0.5 font-mono-tabular text-[11px] text-muted-foreground">
                      {percent(decision.confidence)} confidence
                    </div>
                  </div>
                  <span
                    className="shrink-0 rounded-sm px-2 py-1 font-mono-tabular text-[10px] font-medium tracking-wider uppercase"
                    style={{ background: style.bg, color: style.fg }}
                  >
                    {style.label}
                  </span>
                </button>
              );
            })}
          </div>
        </aside>

        <section className="flex h-[calc(100vh-52px)] flex-col overflow-hidden">
          <DetailHeader selected={selected} />
          <ReportWarnings report={report} />
          <div className="flex-1 overflow-hidden">
            <Tabs defaultValue="bio" className="flex h-full flex-col">
              <div className="border-b border-border bg-card px-8">
                <TabsList className="h-auto gap-1 rounded-none border-none bg-transparent p-0">
                  <TabsTrigger value="bio" className="rounded-none border-b-2 border-transparent px-4 py-3 font-mono-tabular text-[11px] tracking-[0.2em] uppercase data-[state=active]:border-foreground data-[state=active]:bg-transparent data-[state=active]:shadow-none">
                    Bio
                  </TabsTrigger>
                  <TabsTrigger value="model" className="rounded-none border-b-2 border-transparent px-4 py-3 font-mono-tabular text-[11px] tracking-[0.2em] uppercase data-[state=active]:border-foreground data-[state=active]:bg-transparent data-[state=active]:shadow-none">
                    Model
                  </TabsTrigger>
                </TabsList>
              </div>
              <TabsContent value="bio" className="mt-0 flex-1 overflow-y-auto p-8">
                <BioTab selected={selected} anyPositive={anyPositive} />
              </TabsContent>
              <TabsContent value="model" className="mt-0 flex-1 overflow-y-auto p-8">
                <ModelTab selected={selected} metric={metricsByDrug.get(selected.drug)} />
              </TabsContent>
            </Tabs>
          </div>
        </section>
      </div>

      <footer className="border-t border-border bg-card px-6 py-2 text-center text-[11px] text-muted-foreground">
        {report.safety_notice}
      </footer>
    </div>
  );
}



function ReportWarnings({ report }: { report: GenomeReport }) {
  const messages = [
    !report.scope_ok ? "Genome is outside the supported species scope." : null,
    report.qc?.passed === false ? "Assembly quality control did not pass." : null,
    ...report.warnings,
  ].filter((message): message is string => Boolean(message));

  if (!messages.length) return null;

  return (
    <Alert className="rounded-none border-x-0 border-t-0 border-border bg-muted/50 px-8">
      <AlertTriangle className="h-4 w-4" style={{ color: "var(--warning)" }} />
      <AlertTitle className="font-mono-tabular text-[11px] tracking-wider uppercase">Input and pipeline warnings</AlertTitle>
      <AlertDescription>
        <ul className="list-disc space-y-1 pl-4">
          {messages.map((message) => <li key={message}>{message}</li>)}
        </ul>
      </AlertDescription>
    </Alert>
  );
}
function DashboardMessage({ title, message }: { title: string; message: string }) {
  return (
    <div className="flex min-h-screen items-center justify-center bg-background px-4">
      <div className="max-w-md rounded-md border border-border bg-card p-6 text-center">
        <AlertTriangle className="mx-auto h-5 w-5" style={{ color: "var(--warning)" }} />
        <h1 className="mt-3 text-lg font-semibold">{title}</h1>
        <p className="mt-2 text-sm text-muted-foreground">{message}</p>
        <Link to="/" className="mt-5 inline-flex rounded-md bg-foreground px-4 py-2 font-mono-tabular text-xs tracking-wider text-background">
          Return to sample selection
        </Link>
      </div>
    </div>
  );
}

function DetailHeader({ selected }: { selected: AntibioticDecision }) {
  const style = statusStyles(toStatus(selected.call));
  return (
    <div className="border-b border-border bg-card px-8 py-6">
      <div className="grid grid-cols-[minmax(0,1fr)_auto] items-end gap-8">
        <div className="min-w-0">
          <div className="font-mono-tabular text-[10px] tracking-[0.2em] uppercase text-muted-foreground">
            Antibiotic / Tier {selected.tier}
          </div>
          <h1 className="mt-1 truncate text-4xl font-semibold tracking-tight">{selected.drug_display}</h1>
          <div className="mt-2 text-sm text-muted-foreground">{selected.drug_class}</div>
          <div className="mt-3 inline-flex items-center rounded-sm px-3 py-1.5 font-mono-tabular text-xs font-semibold tracking-[0.15em] uppercase" style={{ background: style.bg, color: style.fg }}>
            {style.label}
          </div>
        </div>
        <div className="text-right">
          <div className="font-mono-tabular text-[10px] tracking-[0.2em] uppercase text-muted-foreground">Confidence Score</div>
          <div className="font-mono-tabular text-5xl font-semibold tracking-tight" style={{ color: style.bg }}>
            {percent(selected.confidence, 2)}
          </div>
          {selected.p_resistant != null && <div className="mt-2 font-mono-tabular text-[11px] text-muted-foreground">P(resistant) {percent(selected.p_resistant, 1)}</div>}
        </div>
      </div>
    </div>
  );
}

function BioTab({ selected, anyPositive }: { selected: AntibioticDecision; anyPositive: boolean }) {
  const showAmr = selected.supporting_determinants.length > 0;
  const showTargets = selected.target_evidence.length > 0;
  const nothing = !showAmr && !showTargets;
  const status = toStatus(selected.call);

  return (
    <div className="space-y-6">
      {status === "nocall" && (
        <Alert className="rounded-md border-border bg-muted/50">
          <AlertTriangle className="h-4 w-4" style={{ color: "var(--warning)" }} />
          <AlertTitle className="font-mono-tabular text-[11px] tracking-wider uppercase">No-call</AlertTitle>
          <AlertDescription>{selected.rationale}</AlertDescription>
        </Alert>
      )}

      {!anyPositive && (
        <Alert className="rounded-md border-border" style={{ background: "color-mix(in oklab, var(--warning) 12%, transparent)" }}>
          <AlertTriangle className="h-4 w-4" />
          <AlertTitle className="font-mono-tabular text-[11px] tracking-wider uppercase">No positive classifications</AlertTitle>
          <AlertDescription>No standard antibiotic passed the success threshold. Consider bacteriophage matching for the identified strain profile.</AlertDescription>
        </Alert>
      )}

      {showAmr && (
        <BioTable
          title="AMR Markers"
          subtitle="M1 determinants supporting the prediction. The current M5 report does not include mutation coordinates."
          accent="var(--danger)"
          columns={["AMR feature", "Prediction context", "Coordinate", "Evidence"]}
          rows={selected.supporting_determinants.map((marker) => [
            marker,
            selected.call,
            "Not emitted by API",
            evidenceLabel(selected.evidence_category),
          ])}
        />
      )}

      {!showAmr && selected.evidence_category === "statistical_association" && (
        <Alert className="rounded-md border-border bg-muted/50">
          <AlertTriangle className="h-4 w-4" style={{ color: "var(--warning)" }} />
          <AlertTitle className="font-mono-tabular text-[11px] tracking-wider uppercase">Statistical resistance signal</AlertTitle>
          <AlertDescription>The model found an association but the report contains no catalogued AMR determinant for this antibiotic.</AlertDescription>
        </Alert>
      )}

      {showTargets && (
        <BioTable
          title="Drug Target Mapping"
          subtitle="M2 target evidence physically detected in the genome"
          accent="var(--success)"
          columns={["Target protein", "ORF / contig", "Coverage", "Identity"]}
          rows={selected.target_evidence.map(targetRow)}
        />
      )}

      {!showTargets && selected.target_status !== "present" && (
        <div className="rounded-md border border-dashed border-border bg-card p-6 text-sm text-muted-foreground">
          Target status: <span className="font-mono-tabular text-foreground">{selected.target_status}</span>.
          {selected.target_status === "not_applicable"
            ? " This drug has no single detectable protein target in the current M2 reference panel."
            : " No qualifying target hit was included in the report."}
        </div>
      )}

      {nothing && status === "nocall" && (
        <div className="rounded-md border border-dashed border-border bg-card p-8 text-center text-sm text-muted-foreground">
          No AMR markers or drug targets were returned for this No-call decision.
        </div>
      )}
    </div>
  );
}

function targetRow(target: TargetEvidence): string[] {
  return [
    target.gene,
    target.orf_id || target.contig || "Not emitted by API",
    percent(target.ref_coverage),
    percent(target.identity),
  ];
}

function BioTable({
  title,
  subtitle,
  accent,
  columns,
  rows,
}: {
  title: string;
  subtitle: string;
  accent: string;
  columns: string[];
  rows: string[][];
}) {
  return (
    <div className="overflow-hidden rounded-md border border-border bg-card">
      <div className="flex items-center gap-3 border-b border-border px-5 py-4">
        <div className="h-4 w-1 rounded-sm" style={{ background: accent }} />
        <div className="min-w-0 flex-1">
          <div className="text-sm font-semibold tracking-tight">{title}</div>
          <div className="mt-0.5 text-xs text-muted-foreground">{subtitle}</div>
        </div>
        <div className="font-mono-tabular text-[10px] tracking-wider uppercase text-muted-foreground">{rows.length} rows</div>
      </div>
      <div className="overflow-x-auto">
        <table className="w-full font-mono-tabular text-xs">
          <thead>
            <tr className="border-b border-border bg-muted/40">
              {columns.map((column) => (
                <th key={column} className="px-5 py-2.5 text-left text-[10px] font-medium tracking-[0.15em] uppercase text-muted-foreground">
                  {column}
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {rows.map((row, rowIndex) => (
              <tr key={rowIndex} className="border-b border-border last:border-0 hover:bg-accent/50">
                {row.map((cell, cellIndex) => <td key={cellIndex} className="px-5 py-3 text-foreground">{cell}</td>)}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

function ModelTab({ selected, metric }: { selected: AntibioticDecision; metric?: PerformanceMetric }) {
  const values = [
    { label: "Brier Score", value: decimal(metric?.brier), hint: "lower is better" },
    { label: "AUROC", value: decimal(metric?.auroc), hint: "0-1" },
    { label: "PR-AUC", value: decimal(metric?.pr_auc), hint: "0-1" },
    { label: "Balanced Accuracy", value: decimal(metric?.balanced_accuracy), hint: "0-1" },
    { label: "Recall", value: decimal(metric?.recall_resistant), hint: "resistant class" },
    { label: "F1 Score", value: decimal(metric?.f1), hint: "harmonic mean" },
  ];

  return (
    <div className="space-y-6">
      {!metric && (
        <Alert className="rounded-md border-border bg-muted/50">
          <AlertTriangle className="h-4 w-4" style={{ color: "var(--warning)" }} />
          <AlertTitle className="font-mono-tabular text-[11px] tracking-wider uppercase">No benchmark metrics</AlertTitle>
          <AlertDescription>No held-out metric row was emitted for {selected.drug_display}. This is expected for some structural No-call drugs.</AlertDescription>
        </Alert>
      )}

      <div className="grid grid-cols-2 gap-3 md:grid-cols-3">
        {values.map((value) => (
          <div key={value.label} className="rounded-md border border-border bg-card p-4">
            <div className="font-mono-tabular text-[10px] tracking-[0.2em] uppercase text-muted-foreground">{value.label}</div>
            <div className="mt-1 font-mono-tabular text-2xl font-semibold tracking-tight">{value.value}</div>
            <div className="mt-0.5 text-[11px] text-muted-foreground">{value.hint}</div>
          </div>
        ))}
      </div>

      <div className="grid grid-cols-1 gap-4 lg:grid-cols-2">
        <ReliabilityPlot data={metric?.reliability} />
        <FeatureImportancePlaceholder />
      </div>
    </div>
  );
}

function ReliabilityPlot({ data }: { data?: PerformanceMetric["reliability"] }) {
  const reliability = data;
  const points = reliability
    ? reliability.mean_predicted
        .map((prediction, index) => [prediction, reliability.observed_frequency[index]] as const)
        .filter(([prediction, observed]) => Number.isFinite(prediction) && Number.isFinite(observed))
    : [];

  return (
    <div className="rounded-md border border-border bg-card p-5">
      <div className="mb-4">
        <div className="text-sm font-semibold tracking-tight">Reliability Plot</div>
        <div className="mt-0.5 text-xs text-muted-foreground">Real held-out calibration values from the per-drug metric report</div>
      </div>
      {points.length ? (
        <svg viewBox="0 0 200 200" className="h-64 w-full">
          {Array.from({ length: 11 }).map((_, index) => (
            <g key={index}>
              <line x1={20 + index * 16} y1={10} x2={20 + index * 16} y2={180} stroke="currentColor" strokeOpacity={0.06} />
              <line x1={20} y1={10 + index * 17} x2={180} y2={10 + index * 17} stroke="currentColor" strokeOpacity={0.06} />
            </g>
          ))}
          <line x1={20} y1={180} x2={180} y2={10} stroke="currentColor" strokeOpacity={0.35} strokeDasharray="3 3" />
          <polyline fill="none" stroke="var(--success)" strokeWidth={1.8} points={points.map(([x, y]) => `${20 + x * 160},${180 - y * 170}`).join(" ")} />
          {points.map(([x, y], index) => <circle key={index} cx={20 + x * 160} cy={180 - y * 170} r={2.5} fill="var(--success)" />)}
          <text x={100} y={196} textAnchor="middle" fontSize="8" fill="currentColor" opacity={0.6}>predicted probability</text>
        </svg>
      ) : (
        <div className="flex h-64 items-center justify-center rounded-md border border-dashed border-border text-center text-sm text-muted-foreground">
          Reliability data was not emitted for this drug.
        </div>
      )}
    </div>
  );
}

function FeatureImportancePlaceholder() {
  return (
    <div className="rounded-md border border-border bg-card p-5">
      <div className="mb-4">
        <div className="text-sm font-semibold tracking-tight">Feature Importance</div>
        <div className="mt-0.5 text-xs text-muted-foreground">Model interpretability / k-mer association values</div>
      </div>
      <div className="flex h-64 items-center justify-center rounded-md border border-dashed border-border bg-muted/30 px-8 text-center text-sm text-muted-foreground">
        The current M3/M5 API does not emit feature-attribution or SHAP vectors. This panel intentionally shows no mock values.
      </div>
    </div>
  );
}
