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
  fetchRawM1,
  loadStoredAnalysis,
  rawWorkflowUrl,
} from "@/lib/genome-firewall-api";
import {
  AlertTriangle,
  ArrowLeft,
  CheckCircle2,
  CircleHelp,
  Dna,
  Download,
  Loader2,
  ShieldCheck,
} from "lucide-react";

export const Route = createFileRoute("/dashboard")({
  validateSearch: z.object({ analysisId: z.string().optional() }),
  component: DashboardScreen,
  head: () => ({ meta: [{ title: "Analysis — Genome Firewall" }] }),
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
  const [rawM1, setRawM1] = useState<string | null>(null);
  const [rawError, setRawError] = useState<string | null>(null);
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
    fetchRawM1(stored.report.analysis_id, controller.signal)
      .then(setRawM1)
      .catch((error: unknown) => {
        if ((error as DOMException)?.name !== "AbortError") {
          setRawError(error instanceof Error ? error.message : "Raw M1 output unavailable");
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
            <ShieldCheck className="h-4 w-4" />
            <span className="font-mono-tabular text-[11px] font-semibold tracking-[0.2em] uppercase">
              Genome Firewall
            </span>
          </div>
        </div>
        <div className="flex flex-wrap gap-x-5 gap-y-1 font-mono-tabular text-[10px] tracking-wider uppercase text-muted-foreground">
          <div>
            Sample <span className="text-foreground">{sampleName}</span>
          </div>
          <div>
            Scope <span className="italic text-foreground">{report.species_scope}</span>
          </div>
          <div>
            Drugs <span className="text-foreground">{report.decisions.length}</span>
          </div>
        </div>
      </header>

      <Alert className="m-3 rounded-md border-warning/40 bg-warning/10 sm:mx-6">
        <AlertTriangle className="h-4 w-4 text-warning" />
        <AlertTitle className="font-mono-tabular text-[11px] tracking-wider uppercase">
          Research decision support only
        </AlertTitle>
        <AlertDescription>{report.warning}</AlertDescription>
      </Alert>

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
                      {pct(decision.confidence)} confidence
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
              <TabsList className="h-auto gap-1 rounded-none border-none bg-transparent p-0">
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
                  value="raw"
                  className="rounded-none border-b-2 border-transparent px-4 py-3 font-mono-tabular text-[11px] tracking-[0.2em] uppercase data-[state=active]:border-foreground data-[state=active]:bg-transparent data-[state=active]:shadow-none"
                >
                  Raw outputs
                </TabsTrigger>
              </TabsList>
            </div>
            <TabsContent value="bio" className="m-0 p-5 sm:p-8">
              <BioTab selected={selected} report={report} models={models} />
            </TabsContent>
            <TabsContent value="raw" className="m-0 p-5 sm:p-8">
              <RawOutputs report={report} rawM1={rawM1} rawError={rawError} />
            </TabsContent>
            <TabsContent value="model" className="m-0 p-5 sm:p-8">
              <ModelTab
                selected={selected}
                report={report}
                models={models}
                modelsError={modelsError}
              />
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
            Calibrated confidence
          </div>
          <div
            className="font-mono-tabular text-4xl font-semibold tracking-tight sm:text-5xl"
            style={{ color: style.bg }}
          >
            {pct(selected.confidence, 2)}
          </div>
          <div className="mt-1 font-mono-tabular text-[10px] text-muted-foreground">
            P(resistant) {pct(selected.probability_resistant, 2)}
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
              <p className="leading-5 text-muted-foreground">
                Thresholds create the model signal. QC, lineage, feature novelty, known resistance
                evidence, and target presence can still turn a directional signal into a no-call.
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

  return (
    <div className="space-y-6">
      {selected.call === "no_call" && (
        <Alert className="rounded-md border-warning/40 bg-warning/10">
          <AlertTriangle className="h-4 w-4 text-warning" />
          <AlertTitle className="font-mono-tabular text-[11px] tracking-wider uppercase">
            No-call is intentional
          </AlertTitle>
          <AlertDescription>
            P(resistant) is {pct(selected.probability_resistant, 2)}. A likely-to-work call requires
            it to be below {pct(thresholds?.lower, 2)}.{" "}
            {thresholds?.upper == null
              ? "No safe likely-to-fail threshold was calibrated for this drug."
              : `Likely-to-fail requires it above ${pct(thresholds.upper, 2)}.`}{" "}
            This value is inside the calibrated no-call band. No known marker is not proof of
            susceptibility.
          </AlertDescription>
        </Alert>
      )}

      <div className="grid gap-3 sm:grid-cols-3">
        <SummaryCell label="Evidence category" value={readable(selected.evidence_category)} />
        <SummaryCell
          label="M1 resistance evidence"
          value={`${amr.length} relevant feature${amr.length === 1 ? "" : "s"}`}
        />
        <SummaryCell
          label="M2 molecular target"
          value={`${readable(selected.target_status)} · ${selected.targets_detected.length} detected`}
        />
      </div>

      {amr.length > 0 ? (
        <BioTable
          title="M1 · AMRFinderPlus evidence"
          subtitle="Known resistance genes or mutations relevant to this antibiotic; association does not by itself prove causality."
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

      {target?.evidence.length ? (
        <BioTable
          title="M2 · Drug target verification"
          subtitle={`${selected.target_label}. ${report.workflows.M2.method}`}
          accent="var(--success)"
          columns={["Target", "Reference / method", "Contig / ORF", "Coverage", "Identity"]}
          rows={target.evidence.map((item) => [
            item.symbol,
            item.reference ?? item.method ?? target.probe_kind,
            item.orf_id ?? item.contig ?? "—",
            pct(item.reference_coverage),
            pct(item.identity),
          ])}
        />
      ) : (
        <EmptyEvidence
          title="No target evidence available"
          body="The system will not make a likely-to-work call solely from the absence of resistance markers."
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

function RawOutputs({
  report,
  rawM1,
  rawError,
}: {
  report: AnalysisReport;
  rawM1: string | null;
  rawError: string | null;
}) {
  return (
    <div className="space-y-6">
      <Alert className="rounded-md border-border bg-muted/30">
        <Dna className="h-4 w-4" />
        <AlertTitle className="font-mono-tabular text-[11px] tracking-wider uppercase">
          Uninterpreted workflow output
        </AlertTitle>
        <AlertDescription>
          These records are provided for expert inspection and reproducibility. They do not replace
          the calibrated decision report or laboratory testing.
        </AlertDescription>
      </Alert>

      <RawPanel
        title="M1 · Raw AMRFinderPlus TSV"
        subtitle="Unmodified tab-separated output produced by AMRFinderPlus."
        href={rawWorkflowUrl(report.analysis_id, "m1")}
        content={rawM1}
        error={rawError}
      />
      <RawPanel
        title="M2 · Raw target-detection JSON"
        subtitle="Complete PyHMMER protein-homology and BLASTN RNA-target output."
        href={rawWorkflowUrl(report.analysis_id, "m2")}
        content={JSON.stringify(report.workflows.M2, null, 2)}
      />
    </div>
  );
}

function RawPanel({
  title,
  subtitle,
  href,
  content,
  error,
}: {
  title: string;
  subtitle: string;
  href: string;
  content: string | null;
  error?: string | null;
}) {
  return (
    <div className="overflow-hidden rounded-md border border-border bg-card">
      <div className="flex flex-wrap items-center justify-between gap-3 border-b border-border px-5 py-4">
        <div>
          <div className="text-sm font-semibold">{title}</div>
          <div className="mt-0.5 text-xs text-muted-foreground">{subtitle}</div>
        </div>
        <Button asChild variant="secondary" size="sm">
          <a href={href} target="_blank" rel="noreferrer">
            <Download className="mr-1 h-3.5 w-3.5" /> Open raw output
          </a>
        </Button>
      </div>
      {error ? (
        <div className="p-5 text-sm text-danger">{error}</div>
      ) : content == null ? (
        <div className="flex items-center gap-2 p-5 text-sm text-muted-foreground">
          <Loader2 className="h-4 w-4 animate-spin" /> Loading raw output…
        </div>
      ) : (
        <pre className="max-h-[420px] overflow-auto whitespace-pre p-5 font-mono-tabular text-[11px] leading-5 text-foreground">
          {content}
        </pre>
      )}
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
            ["Model signal", readable(selected.model_signal)],
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
            ["Unknown AMR features", String(selected.unknown_features.length)],
            ["Evidence category", readable(selected.evidence_category)],
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
