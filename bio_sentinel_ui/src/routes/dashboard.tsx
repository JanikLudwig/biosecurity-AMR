import { createFileRoute, Link } from "@tanstack/react-router";
import { useEffect, useMemo, useState } from "react";
import { z } from "zod";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { Button } from "@/components/ui/button";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import {
  type AnalysisDecision,
  type AnalysisReport,
  type DrugEvaluation,
  type ModelsResponse,
  displayDrugName,
  fetchAnalysis,
  fetchModels,
  loadStoredAnalysis,
} from "@/lib/genome-firewall-api";
import {
  AlertTriangle,
  Activity,
  ArrowLeft,
  CheckCircle2,
  CircleHelp,
  Dna,
  FileCheck2,
  GitMerge,
  Loader2,
  ScanSearch,
  Target,
} from "lucide-react";

export const Route = createFileRoute("/dashboard")({
  validateSearch: z.object({ analysisId: z.string().optional() }),
  component: DashboardScreen,
  head: () => ({ meta: [{ title: "Analysis — GyraseX" }] }),
});

type StatusStyle = { bg: string; fg: string; label: string };

function statusStyles(call: AnalysisDecision["call"]): StatusStyle {
  if (call === "likely_to_work")
    return { bg: "var(--success)", fg: "#fff", label: "Likely to work" };
  if (call === "likely_to_fail")
    return { bg: "var(--danger)", fg: "#fff", label: "Likely to fail" };
  return { bg: "var(--warning)", fg: "#1a1a1a", label: "No-call" };
}

function readable(value: string): string {
  return value.replaceAll("_", " ").replace(/\b\w/g, (letter) => letter.toUpperCase());
}

function pct(value: number | null | undefined, digits = 1): string {
  return value == null ? "—" : `${(value * 100).toFixed(digits)}%`;
}

function DashboardScreen() {
  const { analysisId } = Route.useSearch();
  const [stored, setStored] = useState<{ report: AnalysisReport; sampleName: string } | null>(null);
  const [models, setModels] = useState<ModelsResponse | null>(null);
  const [selectedId, setSelectedId] = useState<string | null>(null);
  const [ready, setReady] = useState(false);
  const [modelsError, setModelsError] = useState<string | null>(null);

  useEffect(() => {
    const controller = new AbortController();
    const cached = loadStoredAnalysis();
    if (cached && (!analysisId || cached.report.analysis_id === analysisId)) {
      setStored(cached);
      setSelectedId(cached.report.decisions[0]?.antibiotic ?? null);
      setReady(true);
    } else if (analysisId) {
      fetchAnalysis(analysisId, controller.signal)
        .then((report) => {
          setStored({ report, sampleName: report.genome_id });
          setSelectedId(report.decisions[0]?.antibiotic ?? null);
        })
        .catch((error: unknown) => {
          if ((error as DOMException)?.name !== "AbortError") {
            setModelsError(error instanceof Error ? error.message : "Analysis unavailable");
          }
        })
        .finally(() => {
          if (!controller.signal.aborted) setReady(true);
        });
    } else {
      setReady(true);
    }
    return () => controller.abort();
  }, [analysisId]);

  useEffect(() => {
    if (!stored) return;
    const controller = new AbortController();
    fetchModels(controller.signal)
      .then(setModels)
      .catch((error: unknown) => {
        if ((error as DOMException)?.name !== "AbortError") {
          setModelsError(error instanceof Error ? error.message : "Model metadata unavailable");
        }
      });
    return () => controller.abort();
  }, [stored]);

  if (!ready) {
    return (
      <FullPageMessage
        icon={<Loader2 className="h-5 w-5 animate-spin" />}
        title="Loading analysis"
      />
    );
  }
  if (!stored) {
    return (
      <FullPageMessage
        icon={<Dna className="h-5 w-5" />}
        title="No analysis in this browser session"
        body="Upload an assembled S. aureus FASTA to generate a report."
        action={
          <Button asChild>
            <Link to="/">Start an analysis</Link>
          </Button>
        }
      />
    );
  }

  const { report, sampleName } = stored;
  const selected =
    report.decisions.find((decision) => decision.antibiotic === selectedId) ?? report.decisions[0];
  if (!selected) {
    return (
      <FullPageMessage
        icon={<AlertTriangle className="h-5 w-5" />}
        title="The analysis contains no antibiotic decisions"
      />
    );
  }

  return (
    <div className="flex min-h-screen flex-col bg-background text-foreground">
      <header className="flex flex-wrap items-center justify-between gap-3 border-b border-border bg-card px-4 py-3 sm:px-6">
        <div className="flex items-center gap-5">
          <Link
            to="/"
            className="flex items-center gap-2 text-sm text-muted-foreground hover:text-foreground"
          >
            <ArrowLeft className="h-4 w-4" />
            <span className="font-mono-tabular text-[11px] tracking-wider uppercase">
              New analysis
            </span>
          </Link>
          <div className="flex items-center gap-2">
            <Dna className="h-4 w-4" />
            <span className="font-mono-tabular text-[11px] font-semibold tracking-[0.12em]">
              GyraseX
            </span>
          </div>
        </div>
        <div className="flex flex-wrap gap-x-5 gap-y-1 font-mono-tabular text-[10px] tracking-wider uppercase text-muted-foreground">
          <div>
            Genome ID <span className="text-foreground">{sampleName}</span>
          </div>
          <div>
            Species <span className="italic text-foreground">{report.species_scope}</span>
          </div>
          <div>
            Drugs <span className="text-foreground">{report.decisions.length}</span>
          </div>
        </div>
      </header>

      <div className="grid flex-1 lg:grid-cols-[320px_minmax(0,1fr)]">
        <aside className="border-b border-r-0 border-border bg-card lg:border-b-0 lg:border-r">
          <div className="border-b border-border px-5 py-4">
            <div className="font-mono-tabular text-[10px] tracking-[0.2em] uppercase text-muted-foreground">
              Antibiotic directory
            </div>
            <h2 className="mt-1 text-base font-semibold tracking-tight">Calibrated decisions</h2>
          </div>
          <div className="grid sm:grid-cols-3 lg:block">
            {report.decisions.map((decision) => {
              const style = statusStyles(decision.call);
              const active = decision.antibiotic === selected.antibiotic;
              return (
                <button
                  key={decision.antibiotic}
                  onClick={() => setSelectedId(decision.antibiotic)}
                  className={`flex w-full items-center gap-3 border-b border-border px-5 py-4 text-left transition-colors hover:bg-accent ${active ? "bg-accent" : ""}`}
                >
                  <div
                    className="h-9 w-[3px] shrink-0 rounded-sm"
                    style={{ background: active ? style.bg : "transparent" }}
                  />
                  <div className="min-w-0 flex-1">
                    <div className="truncate text-sm font-medium">
                      {displayDrugName(decision.antibiotic)}
                    </div>
                    <div className="mt-0.5 font-mono-tabular text-[11px] text-muted-foreground">
                      P(resistant) {pct(decision.probability_resistant)}
                    </div>
                  </div>
                  <span
                    className="shrink-0 rounded-sm px-2 py-1 font-mono-tabular text-[9px] font-medium tracking-wider uppercase"
                    style={{ background: style.bg, color: style.fg }}
                  >
                    {style.label}
                  </span>
                </button>
              );
            })}
          </div>
        </aside>

        <section className="min-w-0">
          <DetailHeader selected={selected} models={models} />
          <Tabs defaultValue="bio" className="w-full">
            <div className="border-b border-border bg-card px-5 sm:px-8">
              <TabsList className="h-auto max-w-full gap-1 overflow-x-auto rounded-none border-none bg-transparent p-0">
                <TabsTrigger
                  value="bio"
                  className="rounded-none border-b-2 border-transparent px-4 py-3 font-mono-tabular text-[11px] tracking-[0.2em] uppercase data-[state=active]:border-foreground data-[state=active]:bg-transparent data-[state=active]:shadow-none"
                >
                  Bio evidence
                </TabsTrigger>
                <TabsTrigger
                  value="model"
                  className="rounded-none border-b-2 border-transparent px-4 py-3 font-mono-tabular text-[11px] tracking-[0.2em] uppercase data-[state=active]:border-foreground data-[state=active]:bg-transparent data-[state=active]:shadow-none"
                >
                  Model &amp; QC
                </TabsTrigger>
                <TabsTrigger
                  value="pipeline"
                  className="rounded-none border-b-2 border-transparent px-4 py-3 font-mono-tabular text-[11px] tracking-[0.2em] whitespace-nowrap uppercase data-[state=active]:border-foreground data-[state=active]:bg-transparent data-[state=active]:shadow-none"
                >
                  How this was decided
                </TabsTrigger>
              </TabsList>
            </div>
            <TabsContent value="bio" className="m-0 p-5 sm:p-8">
              <BioTab selected={selected} report={report} models={models} />
            </TabsContent>
            <TabsContent value="model" className="m-0 p-5 sm:p-8">
              <ModelTab
                selected={selected}
                report={report}
                models={models}
                modelsError={modelsError}
              />
            </TabsContent>
            <TabsContent value="pipeline" className="m-0 p-5 sm:p-8">
              <PipelineTab selected={selected} report={report} models={models} />
            </TabsContent>
          </Tabs>
        </section>
      </div>
    </div>
  );
}

function FullPageMessage({
  icon,
  title,
  body,
  action,
}: {
  icon: React.ReactNode;
  title: string;
  body?: string;
  action?: React.ReactNode;
}) {
  return (
    <div className="flex min-h-screen items-center justify-center bg-background p-6">
      <div className="max-w-md rounded-md border border-border bg-card p-8 text-center">
        <div className="mx-auto mb-3 flex h-10 w-10 items-center justify-center rounded-md bg-muted">
          {icon}
        </div>
        <h1 className="text-xl font-semibold">{title}</h1>
        {body && <p className="mt-2 text-sm text-muted-foreground">{body}</p>}
        {action && <div className="mt-5">{action}</div>}
      </div>
    </div>
  );
}

function DetailHeader({
  selected,
  models,
}: {
  selected: AnalysisDecision;
  models: ModelsResponse | null;
}) {
  const style = statusStyles(selected.call);
  const targetPresent = selected.target_status === "present";
  const thresholds = models?.bundle.models[selected.antibiotic]?.decision_thresholds;
  return (
    <div className="border-b border-border bg-card px-5 py-6 sm:px-8">
      <div className="grid items-end gap-5 sm:grid-cols-[minmax(0,1fr)_auto]">
        <div className="min-w-0">
          <div className="font-mono-tabular text-[10px] tracking-[0.2em] uppercase text-muted-foreground">
            Antibiotic
          </div>
          <h1 className="mt-1 truncate text-3xl font-semibold tracking-tight sm:text-4xl">
            {displayDrugName(selected.antibiotic)}
          </h1>
          <div className="mt-3 flex flex-wrap items-center gap-2">
            <span
              className="rounded-sm px-3 py-1.5 font-mono-tabular text-xs font-semibold tracking-[0.15em] uppercase"
              style={{ background: style.bg, color: style.fg }}
            >
              {style.label}
            </span>
            <DecisionInfo selected={selected} thresholds={thresholds} />
            <span
              className="rounded-sm px-3 py-1.5 font-mono-tabular text-[10px] font-semibold tracking-wider uppercase"
              style={{
                background: targetPresent ? "var(--success)" : "var(--muted)",
                color: targetPresent ? "#fff" : "var(--muted-foreground)",
              }}
            >
              Target {readable(selected.target_status)}
            </span>
          </div>
        </div>
        <div className="sm:text-right">
          <div className="font-mono-tabular text-[10px] tracking-[0.2em] uppercase text-muted-foreground">
            Resistance probability
          </div>
          <div
            className="font-mono-tabular text-4xl font-semibold tracking-tight sm:text-5xl"
            style={{ color: style.bg }}
          >
            {pct(selected.probability_resistant, 2)}
          </div>
          <div className="mt-1 font-mono-tabular text-[10px] text-muted-foreground">
            {thresholds ? (
              <>
                Work ≤ {pct(thresholds.lower, 2)} ·{" "}
                {thresholds.upper == null
                  ? "Fail threshold unavailable"
                  : `Fail ≥ ${pct(thresholds.upper, 2)}`}
              </>
            ) : (
              "Loading decision thresholds…"
            )}
          </div>
        </div>
      </div>
    </div>
  );
}

type DecisionThresholds = ModelsResponse["bundle"]["models"][string]["decision_thresholds"];

function DecisionInfo({
  selected,
  thresholds,
}: {
  selected: AnalysisDecision;
  thresholds?: DecisionThresholds;
}) {
  const style = statusStyles(selected.call);
  const reason =
    selected.call === "no_call"
      ? selected.reasons.length
        ? selected.reasons.map(readable).join(" · ")
        : "One or more safety gates did not pass."
      : selected.call === "likely_to_work"
        ? "The susceptible model threshold and all safety gates passed, including target presence and absence of relevant known resistance evidence."
        : "The resistant model threshold passed; resistance evidence is reported separately as known or statistical evidence.";

  return (
    <TooltipProvider delayDuration={150}>
      <Tooltip>
        <TooltipTrigger asChild>
          <button
            type="button"
            aria-label={`Explain ${style.label} decision`}
            className="inline-flex h-7 w-7 items-center justify-center rounded-full border border-border bg-background text-muted-foreground transition-colors hover:text-foreground focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring"
          >
            <CircleHelp className="h-4 w-4" />
          </button>
        </TooltipTrigger>
        <TooltipContent
          side="bottom"
          align="start"
          className="w-[340px] max-w-[calc(100vw-2rem)] border border-border bg-card p-4 text-card-foreground shadow-xl"
        >
          <div className="font-mono-tabular text-[10px] font-semibold tracking-[0.16em] uppercase text-muted-foreground">
            Calibrated decision regions
          </div>
          {thresholds ? (
            <div className="mt-3 space-y-2 text-xs">
              <ThresholdRow
                active={selected.model_signal === "susceptible_signal"}
                color="var(--success)"
                label="Likely to work"
                value={`P(resistant) ≤ ${pct(thresholds.lower, 2)}`}
              />
              <ThresholdRow
                active={selected.model_signal === "no_call"}
                color="var(--warning)"
                label="No-call"
                value={
                  thresholds.upper == null
                    ? `P(resistant) > ${pct(thresholds.lower, 2)}`
                    : `${pct(thresholds.lower, 2)} < P(resistant) < ${pct(thresholds.upper, 2)}`
                }
              />
              <ThresholdRow
                active={selected.model_signal === "resistant_signal"}
                color="var(--danger)"
                label="Likely to fail"
                value={
                  thresholds.upper == null
                    ? "No safe threshold calibrated"
                    : `P(resistant) ≥ ${pct(thresholds.upper, 2)}`
                }
              />
              <div className="mt-3 border-t border-border pt-3 leading-5">
                <span className="font-semibold">This sample:</span>{" "}
                {pct(selected.probability_resistant, 2)} P(resistant) →{" "}
                {readable(selected.model_signal)} → final call {style.label}.
              </div>
              <p className="leading-5 text-muted-foreground">{reason}</p>
              <div className="rounded-sm border border-border bg-muted/40 p-3 leading-5">
                <div className="font-semibold">How these thresholds were learned</div>
                <p className="mt-1 text-muted-foreground">
                  They are selected from a separate calibration set, not chosen by hand. Each
                  candidate boundary had to make at least 5 calls with no more than 10% wrong-class
                  cases. Among valid, non-overlapping pairs, the pipeline chose the pair covering
                  the most calibration cases.
                </p>
                <p className="mt-1 text-muted-foreground">
                  The 5-call and 10% limits are an explicit prototype policy, not a clinical
                  standard; they must be validated on larger external data before real use.
                </p>
                <ul className="mt-2 space-y-1 text-muted-foreground">
                  <li>
                    Work boundary: {thresholds.susceptible_calls} calibration calls;{" "}
                    {pct(thresholds.susceptible_call_error, 1)} were resistant.
                  </li>
                  <li>
                    Fail boundary: {thresholds.resistant_calls} calibration calls;{" "}
                    {pct(thresholds.resistant_call_error, 1)} were susceptible.
                  </li>
                </ul>
              </div>
              <p className="leading-5 text-muted-foreground">
                Thresholds create the model signal. Assembly quality, unfamiliar genomes, known
                resistance markers, and target presence may still change a directional signal to a
                no-call.
              </p>
            </div>
          ) : (
            <div className="mt-3 flex items-center gap-2 text-xs text-muted-foreground">
              <Loader2 className="h-3.5 w-3.5 animate-spin" /> Loading calibrated thresholds…
            </div>
          )}
        </TooltipContent>
      </Tooltip>
    </TooltipProvider>
  );
}

function ThresholdRow({
  active,
  color,
  label,
  value,
}: {
  active: boolean;
  color: string;
  label: string;
  value: string;
}) {
  return (
    <div
      className={`grid grid-cols-[8px_minmax(0,1fr)_auto] items-center gap-2 rounded-sm px-2 py-1.5 ${active ? "bg-muted" : ""}`}
    >
      <span className="h-2 w-2 rounded-full" style={{ background: color }} />
      <span className={active ? "font-semibold" : ""}>{label}</span>
      <span className="font-mono-tabular text-[10px] text-muted-foreground">{value}</span>
    </div>
  );
}

function BioTab({
  selected,
  report,
  models,
}: {
  selected: AnalysisDecision;
  report: AnalysisReport;
  models: ModelsResponse | null;
}) {
  const relevant = useMemo(
    () => new Set([...selected.all_relevant_elements, ...selected.supporting_elements]),
    [selected],
  );
  const amr = report.workflows.M1.evidence.filter((item) => relevant.has(item.element_symbol));
  const target = report.workflows.M2.drugs[selected.antibiotic];
  const thresholds = models?.bundle.models[selected.antibiotic]?.decision_thresholds;
  const verifiedTargetEvidence = target?.evidence.filter((item) => item.present) ?? [];
  const targetMethod = target?.probe_kind.includes("rrna") ? "BLASTN" : "PyHMMER";
  const targetNotVerified = selected.reasons.includes("molecular_target_not_verified");
  const probabilityUncertain = selected.reasons.includes(
    "calibrated_probability_inside_no_call_region",
  );

  return (
    <div className="space-y-6">
      {selected.call === "no_call" && (
        <Alert className="rounded-md border-warning/40 bg-warning/10">
          <AlertTriangle className="h-4 w-4 text-warning" />
          <AlertTitle className="font-mono-tabular text-[11px] tracking-wider uppercase">
            No-call is intentional
          </AlertTitle>
          <AlertDescription>
            {targetNotVerified ? (
              <>
                P(resistant) is {pct(selected.probability_resistant, 2)}, which passes the{" "}
                likely-to-work probability threshold of {pct(thresholds?.lower, 2)}. The final call
                is still withheld because {targetMethod} did not verify the required molecular
                target in this assembly.
              </>
            ) : probabilityUncertain ? (
              <>
                P(resistant) is {pct(selected.probability_resistant, 2)}. A likely-to-work call
                requires it at or below {pct(thresholds?.lower, 2)}.{" "}
                {thresholds?.upper == null
                  ? "No safe likely-to-fail threshold was calibrated for this drug."
                  : `Likely-to-fail requires it at or above ${pct(thresholds.upper, 2)}.`}{" "}
                This value is inside the calibrated no-call band.
              </>
            ) : (
              <>One or more safety checks prevented a directional call: {selected.reasons.map(readable).join(" · ")}.</>
            )}
          </AlertDescription>
        </Alert>
      )}

      <div className="grid gap-3 sm:grid-cols-2">
        <SummaryCell
          label="Known resistance markers"
          value={`${amr.length} relevant marker${amr.length === 1 ? "" : "s"} detected`}
        />
        <SummaryCell
          label="Drug target check"
          value={
            selected.target_status === "present"
              ? `Passed · ${selected.targets_detected.length} target${selected.targets_detected.length === 1 ? "" : "s"} detected`
              : `Verification failed · no qualifying ${targetMethod} hit`
          }
        />
      </div>

      {selected.evidence_category === "statistical_association_only" && (
        <Alert className="rounded-md border-warning/40 bg-warning/10">
          <AlertTriangle className="h-4 w-4 text-warning" />
          <AlertTitle className="font-mono-tabular text-[11px] tracking-wider uppercase">
            Model-only resistance signal
          </AlertTitle>
          <AlertDescription>
            The model found a combination of features associated with resistance, but no known
            resistance gene or mutation was identified for this antibiotic. This is statistical
            evidence, not a proven biological cause.
          </AlertDescription>
        </Alert>
      )}

      {amr.length > 0 ? (
        <BioTable
          title="Known resistance markers"
          subtitle="Detected by AMRFinderPlus (workflow M1). These genes or mutations are relevant to this antibiotic."
          accent="var(--danger)"
          columns={["Element", "Feature", "AMR subclass", "Method", "Coverage / identity"]}
          rows={amr.map((item) => [
            item.element_symbol,
            item.feature_key,
            item.amr_subclass || "—",
            item.method,
            `${item.coverage.toFixed(1)}% / ${item.identity.toFixed(1)}%`,
          ])}
        />
      ) : (
        <EmptyEvidence
          title="No relevant known resistance marker detected"
          body="This is not proof that the antibiotic will work. The statistical model, target gate, lineage check, and calibrated thresholds still apply."
        />
      )}

      {selected.target_status === "present" && verifiedTargetEvidence.length > 0 ? (
        <BioTable
          title="Drug target detected"
          subtitle={`${selected.target_label}. A qualifying ${targetMethod} match passed workflow M2.`}
          accent="var(--success)"
          columns={
            target?.probe_kind.includes("rrna")
              ? ["Target", "Method", "Result"]
              : ["Target", "Reference", "Contig / ORF", "Coverage", "Identity"]
          }
          rows={verifiedTargetEvidence.map((item) =>
            target?.probe_kind.includes("rrna")
              ? [item.symbol, item.method ?? "BLASTN", "Qualifying nucleotide match"]
              : [
                  item.symbol,
                  item.reference ?? "—",
                  item.orf_id ?? item.contig ?? "—",
                  pct(item.reference_coverage),
                  pct(item.identity),
                ],
          )}
        />
      ) : (
        <EmptyEvidence
          title="Target verification failed"
          body={`${targetMethod} completed, but no match passed the target-presence criteria. This may reflect a fragmented or missing target sequence in the assembly; it is not evidence that the target is biologically absent.`}
        />
      )}
    </div>
  );
}

function SummaryCell({ label, value }: { label: string; value: string }) {
  return (
    <div className="rounded-md border border-border bg-card p-4">
      <div className="font-mono-tabular text-[10px] tracking-[0.16em] uppercase text-muted-foreground">
        {label}
      </div>
      <div className="mt-1 text-sm font-medium">{value}</div>
    </div>
  );
}

function EmptyEvidence({ title, body }: { title: string; body: string }) {
  return (
    <div className="rounded-md border border-dashed border-border bg-card p-6">
      <div className="flex items-center gap-2 text-sm font-semibold">
        <CheckCircle2 className="h-4 w-4 text-muted-foreground" />
        {title}
      </div>
      <p className="mt-2 text-xs leading-5 text-muted-foreground">{body}</p>
    </div>
  );
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
          <div className="mt-0.5 text-xs leading-5 text-muted-foreground">{subtitle}</div>
        </div>
        <div className="font-mono-tabular text-[10px] tracking-wider uppercase text-muted-foreground">
          {rows.length} rows
        </div>
      </div>
      <div className="overflow-x-auto">
        <table className="w-full min-w-[720px] font-mono-tabular text-xs">
          <thead>
            <tr className="border-b border-border bg-muted/40">
              {columns.map((column) => (
                <th
                  key={column}
                  className="px-5 py-2.5 text-left text-[10px] font-medium tracking-[0.15em] uppercase text-muted-foreground"
                >
                  {column}
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {rows.map((row, rowIndex) => (
              <tr
                key={rowIndex}
                className="border-b border-border last:border-0 hover:bg-accent/50"
              >
                {row.map((cell, cellIndex) => (
                  <td key={cellIndex} className="px-5 py-3 text-foreground">
                    {cell}
                  </td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

function PipelineTab({
  selected,
  report,
  models,
}: {
  selected: AnalysisDecision;
  report: AnalysisReport;
  models: ModelsResponse | null;
}) {
  const relevantMarkers = report.workflows.M1.evidence.filter((item) =>
    selected.all_relevant_elements.includes(item.element_symbol),
  );
  const target = report.workflows.M2.drugs[selected.antibiotic];
  const thresholds = models?.bundle.models[selected.antibiotic]?.decision_thresholds;
  const decisionStyle = statusStyles(selected.call);
  const decisionReason =
    selected.call === "likely_to_work"
      ? "Low resistance probability, no conflicting known marker, and the drug target was verified."
      : selected.call === "likely_to_fail"
        ? selected.supporting_elements.length
          ? `High resistance probability with supporting marker${selected.supporting_elements.length === 1 ? "" : "s"}: ${selected.supporting_elements.join(", ")}.`
          : "High resistance probability from a statistical pattern; no known causal marker was detected."
        : selected.reasons.length
          ? selected.reasons.map(readable).join(" · ")
          : "One or more safety checks did not support a directional call.";

  return (
    <div className="space-y-6">
      <div>
        <div className="font-mono-tabular text-[10px] tracking-[0.2em] uppercase text-muted-foreground">
          Decision path for {displayDrugName(selected.antibiotic)}
        </div>
        <h2 className="mt-1 text-xl font-semibold tracking-tight">How this result was produced</h2>
        <p className="mt-2 max-w-3xl text-sm leading-6 text-muted-foreground">
          The resistance model and target detector run as separate branches. Target evidence never
          changes the model probability; both branches meet only at the final decision step.
        </p>
      </div>

      <div className="grid items-center gap-4 lg:grid-cols-[minmax(0,1fr)_64px_minmax(0,1fr)]">
        <div className="space-y-3 rounded-md border border-border bg-muted/20 p-3">
          <div className="font-mono-tabular text-[10px] tracking-[0.16em] uppercase text-muted-foreground">
            Resistance branch
          </div>
          <PipelineStage
            icon={<ScanSearch className="h-4 w-4" />}
            module="M1"
            title="Read resistance markers"
            status="Complete"
            color="var(--success)"
            summary="AMRFinderPlus converted detected genes and mutations into model features."
            details={[
              `${report.workflows.M1.recognized_features.length} recognized feature${report.workflows.M1.recognized_features.length === 1 ? "" : "s"} in this genome`,
              `${relevantMarkers.length} relevant to ${displayDrugName(selected.antibiotic)}`,
            ]}
          />
          <FlowArrow label="Features flow into the predictor" />
          <PipelineStage
            icon={<Activity className="h-4 w-4" />}
            module="M3"
            title="Estimate resistance probability"
            status={selected.model_signal === "no_call" ? "Uncertain region" : "Directional signal"}
            color={selected.model_signal === "no_call" ? "var(--warning)" : "var(--success)"}
            summary={`The calibrated logistic-regression model estimated ${pct(selected.probability_resistant, 2)} probability of resistance.`}
            details={[
              thresholds
                ? `Work ≤ ${pct(thresholds.lower, 2)} · ${thresholds.upper == null ? "fail threshold unavailable" : `fail ≥ ${pct(thresholds.upper, 2)}`}`
                : "Loading calibrated thresholds",
            ]}
          />
        </div>

        <div className="flex flex-col items-center justify-center gap-2 text-center text-muted-foreground">
          <GitMerge className="h-6 w-6" />
          <span className="font-mono-tabular text-[9px] tracking-wider uppercase">Join at M4</span>
        </div>

        <div className="space-y-3 rounded-md border border-border bg-muted/20 p-3">
          <div className="font-mono-tabular text-[10px] tracking-[0.16em] uppercase text-muted-foreground">
            Independent target branch
          </div>
          <PipelineStage
            icon={<Target className="h-4 w-4" />}
            module="M2"
            title="Verify the drug target"
            status={selected.target_status === "present" ? "Target present" : "Not verified"}
            color={selected.target_status === "present" ? "var(--success)" : "var(--warning)"}
            summary={
              selected.target_status === "present"
                ? `${selected.target_label} was detected, so the target gate passed.`
                : "The required molecular target could not be verified."
            }
            details={[
              target?.detected.length
                ? `Detected: ${target.detected.join(", ")}`
                : "No qualifying target hit",
              report.workflows.M2.method,
            ]}
          />
        </div>
      </div>

      <FlowArrow label="The two independent branches meet" />

      <div className="mx-auto max-w-3xl space-y-3">
        <PipelineStage
          icon={<GitMerge className="h-4 w-4" />}
          module="M4"
          title="Apply thresholds and safety checks"
          status={decisionStyle.label}
          color={decisionStyle.bg}
          summary={decisionReason}
          details={[
            `Assembly QC: ${report.qc.passed ? "passed" : "failed"}`,
            `Training-population check: ${readable(report.lineage.status)}`,
          ]}
        />
        <FlowArrow label="Decision becomes a report" />
        <PipelineStage
          icon={<FileCheck2 className="h-4 w-4" />}
          module="M5"
          title="Create the decision report"
          status="Report ready"
          color="var(--success)"
          summary="The final call, probability, supporting evidence, QC and provenance are packaged for review."
          details={[
            `Analysis ${report.analysis_id}`,
            "Standard laboratory susceptibility testing is still required",
          ]}
        />
      </div>
    </div>
  );
}

function PipelineStage({
  icon,
  module,
  title,
  status,
  color,
  summary,
  details,
}: {
  icon: React.ReactNode;
  module: string;
  title: string;
  status: string;
  color: string;
  summary: string;
  details: string[];
}) {
  return (
    <section
      className="rounded-md border border-border bg-card p-4"
      aria-label={`${module} ${title}`}
    >
      <div className="flex flex-wrap items-start justify-between gap-3">
        <div className="flex items-center gap-3">
          <div className="flex h-8 w-8 items-center justify-center rounded-md bg-muted">{icon}</div>
          <div>
            <div className="font-mono-tabular text-[10px] tracking-[0.16em] uppercase text-muted-foreground">
              {module}
            </div>
            <h3 className="text-sm font-semibold">{title}</h3>
          </div>
        </div>
        <span
          className="rounded-sm px-2 py-1 font-mono-tabular text-[9px] font-semibold tracking-wider uppercase"
          style={{ background: color, color: color === "var(--warning)" ? "#1a1a1a" : "#fff" }}
        >
          {status}
        </span>
      </div>
      <p className="mt-3 text-sm leading-5">{summary}</p>
      <ul className="mt-3 space-y-1 text-xs leading-5 text-muted-foreground">
        {details.map((detail) => (
          <li key={detail} className="flex gap-2">
            <span aria-hidden>·</span>
            <span>{detail}</span>
          </li>
        ))}
      </ul>
    </section>
  );
}

function FlowArrow({ label }: { label: string }) {
  return (
    <div
      className="flex items-center justify-center gap-2 py-0.5 text-muted-foreground"
      aria-label={label}
    >
      <span className="font-mono-tabular text-lg leading-none">↓</span>
      <span className="font-mono-tabular text-[9px] tracking-wider uppercase">{label}</span>
    </div>
  );
}

function ModelTab({
  selected,
  report,
  models,
  modelsError,
}: {
  selected: AnalysisDecision;
  report: AnalysisReport;
  models: ModelsResponse | null;
  modelsError: string | null;
}) {
  const evaluation: DrugEvaluation | undefined =
    models?.bundle.models[selected.antibiotic]?.evaluation;
  const thresholds = models?.bundle.models[selected.antibiotic]?.decision_thresholds;
  const reliability = models?.bundle.reliability?.[selected.antibiotic] ?? [];
  const metrics = evaluation?.test_metrics;
  return (
    <div className="space-y-6">
      {modelsError && (
        <Alert>
          <AlertTriangle className="h-4 w-4" />
          <AlertTitle>Evaluation metadata unavailable</AlertTitle>
          <AlertDescription>{modelsError}</AlertDescription>
        </Alert>
      )}
      {!models && !modelsError && (
        <div className="flex items-center gap-2 text-sm text-muted-foreground">
          <Loader2 className="h-4 w-4 animate-spin" /> Loading held-out evaluation metadata…
        </div>
      )}

      {metrics && (
        <>
          <div>
            <div className="font-mono-tabular text-[10px] tracking-[0.2em] uppercase text-muted-foreground">
              Held-out evaluation
            </div>
            <div className="mt-1 text-sm font-semibold">
              {readable(models!.bundle.evaluation_status)}
            </div>
            <p className="mt-1 text-xs text-muted-foreground">
              Grouped development split across {models!.bundle.split_genomes} genomes ·{" "}
              {models!.feature_count} AMR features · calibration{" "}
              {readable(evaluation!.calibration_status)}
            </p>
          </div>
          <div className="grid grid-cols-2 gap-3 md:grid-cols-4">
            <MetricCell
              label="Brier score"
              value={metrics.brier.toFixed(3)}
              hint="lower is better"
            />
            <MetricCell label="AUROC" value={metrics.auroc.toFixed(3)} hint="0–1" />
            <MetricCell label="PR-AUC" value={metrics.pr_auc.toFixed(3)} hint="class imbalance" />
            <MetricCell
              label="Balanced accuracy"
              value={metrics.balanced_accuracy.toFixed(3)}
              hint="held-out"
            />
            <MetricCell
              label="Resistant recall"
              value={metrics.resistant_recall.toFixed(3)}
              hint="likely-to-fail cases"
            />
            <MetricCell
              label="Susceptible recall"
              value={metrics.susceptible_recall.toFixed(3)}
              hint="likely-to-work cases"
            />
            <MetricCell label="F1 score" value={metrics.f1.toFixed(3)} hint="resistant class" />
            <MetricCell
              label="No-call rate"
              value={pct(evaluation!.no_call_rate)}
              hint="held-out test"
            />
          </div>
        </>
      )}

      {reliability.length > 0 && <ReliabilityPlot points={reliability} />}

      <div className="grid gap-4 xl:grid-cols-2">
        <AuditCard
          title="Decision audit"
          rows={[
            ["P(resistant)", pct(selected.probability_resistant, 2)],
            ["Likely-to-work threshold", `< ${pct(thresholds?.lower, 2)} P(resistant)`],
            [
              "Likely-to-fail threshold",
              thresholds?.upper == null
                ? "Not safely calibrated"
                : `> ${pct(thresholds.upper, 2)} P(resistant)`,
            ],
            [
              "Feature similarity",
              `${pct(selected.feature_similarity)} (floor ${pct(selected.feature_similarity_floor)})`,
            ],
            [
              "Why this result",
              selected.reasons.length
                ? selected.reasons.map(readable).join(" · ")
                : "Probability and safety checks passed",
            ],
          ]}
        />
        <AuditCard
          title="Genome & lineage QC"
          rows={[
            ["Assembly QC", report.qc.passed ? "Passed" : "Failed"],
            ["Genome length", `${report.qc.genome_length.toLocaleString()} bp`],
            ["Contigs / N50", `${report.qc.contigs} / ${report.qc.contig_n50.toLocaleString()} bp`],
            ["Lineage status", readable(report.lineage.status)],
            ["Nearest training genome", report.lineage.nearest_training_genome ?? "—"],
            ["Maximum training ANI", pct(report.lineage.maximum_training_ani, 3)],
          ]}
        />
      </div>

      <div className="rounded-md border border-border bg-card p-5 text-xs leading-5 text-muted-foreground">
        <div className="font-mono-tabular text-[10px] tracking-[0.2em] uppercase">Provenance</div>
        <p className="mt-2">
          AMRFinderPlus {report.provenance.amrfinder_version} · database{" "}
          {report.provenance.amrfinder_database} · analysis {report.analysis_id}
        </p>
        <p className="mt-1">
          Metrics describe the provisional grouped development split, not clinical validity. No SHAP
          or biological-causality claim is generated by this interface.
        </p>
      </div>
    </div>
  );
}

function ReliabilityPlot({ points }: { points: ModelsResponse["bundle"]["reliability"][string] }) {
  const coordinates = points
    .filter((point) => point.samples > 0)
    .map((point) => [point.mean_probability_resistant, point.observed_resistant_fraction] as const);
  return (
    <div className="rounded-md border border-border bg-card p-5">
      <div className="text-sm font-semibold">Held-out reliability</div>
      <div className="mt-0.5 text-xs text-muted-foreground">
        Mean predicted resistance probability versus observed resistant fraction; point size
        reflects bin count.
      </div>
      <svg
        viewBox="0 0 420 220"
        className="mt-4 h-64 w-full"
        role="img"
        aria-label="Reliability plot"
      >
        {[0, 0.25, 0.5, 0.75, 1].map((tick) => (
          <g key={tick}>
            <line
              x1={40 + tick * 350}
              y1="15"
              x2={40 + tick * 350}
              y2="185"
              stroke="currentColor"
              strokeOpacity="0.08"
            />
            <line
              x1="40"
              y1={185 - tick * 170}
              x2="390"
              y2={185 - tick * 170}
              stroke="currentColor"
              strokeOpacity="0.08"
            />
            <text
              x={40 + tick * 350}
              y="202"
              textAnchor="middle"
              fontSize="9"
              fill="currentColor"
              opacity="0.55"
            >
              {tick}
            </text>
          </g>
        ))}
        <line
          x1="40"
          y1="185"
          x2="390"
          y2="15"
          stroke="currentColor"
          strokeOpacity="0.35"
          strokeDasharray="4 4"
        />
        <polyline
          fill="none"
          stroke="var(--success)"
          strokeWidth="2"
          points={coordinates.map(([x, y]) => `${40 + x * 350},${185 - y * 170}`).join(" ")}
        />
        {points.map((point) => (
          <circle
            key={point.probability_bin}
            cx={40 + point.mean_probability_resistant * 350}
            cy={185 - point.observed_resistant_fraction * 170}
            r={Math.min(8, 3 + Math.sqrt(point.samples))}
            fill="var(--success)"
            fillOpacity="0.8"
          >
            <title>{`${point.samples} samples · predicted ${pct(point.mean_probability_resistant)} · observed ${pct(point.observed_resistant_fraction)}`}</title>
          </circle>
        ))}
        <text x="215" y="217" textAnchor="middle" fontSize="9" fill="currentColor" opacity="0.6">
          Predicted resistance probability
        </text>
        <text
          x="12"
          y="105"
          textAnchor="middle"
          fontSize="9"
          fill="currentColor"
          opacity="0.6"
          transform="rotate(-90 12 105)"
        >
          Observed resistant fraction
        </text>
      </svg>
    </div>
  );
}

function MetricCell({ label, value, hint }: { label: string; value: string; hint: string }) {
  return (
    <div className="rounded-md border border-border bg-card p-4">
      <div className="font-mono-tabular text-[10px] tracking-[0.16em] uppercase text-muted-foreground">
        {label}
      </div>
      <div className="mt-1 font-mono-tabular text-2xl font-semibold">{value}</div>
      <div className="mt-0.5 text-[11px] text-muted-foreground">{hint}</div>
    </div>
  );
}

function AuditCard({ title, rows }: { title: string; rows: string[][] }) {
  return (
    <div className="overflow-hidden rounded-md border border-border bg-card">
      <div className="border-b border-border px-5 py-4 text-sm font-semibold">{title}</div>
      <dl>
        {rows.map(([label, value]) => (
          <div
            key={label}
            className="grid grid-cols-[minmax(0,1fr)_minmax(0,1.3fr)] gap-4 border-b border-border px-5 py-3 text-xs last:border-0"
          >
            <dt className="text-muted-foreground">{label}</dt>
            <dd className="break-words text-right font-mono-tabular">{value}</dd>
          </div>
        ))}
      </dl>
    </div>
  );
}
