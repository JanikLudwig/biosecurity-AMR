import { createFileRoute, useNavigate } from "@tanstack/react-router";
import { useRef, useState } from "react";
import { DnaBackground } from "@/components/dna-background";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { analyzeGenome, sampleIdFromFileName, storeAnalysis } from "@/lib/genome-firewall-api";
import {
  AlertTriangle,
  CheckCircle2,
  Dna,
  FileText,
  Loader2,
  Upload,
} from "lucide-react";

export const Route = createFileRoute("/")({
  component: UploadScreen,
  head: () => ({
    meta: [{ title: "GyraseX — Antibiotic response prediction" }],
  }),
});

type Phase = "idle" | "analyzing" | "error";
const ACCEPTED_SUFFIXES = [".fna", ".fa", ".fasta"];

function validateFile(file: File): string | null {
  const lower = file.name.toLowerCase();
  if (!ACCEPTED_SUFFIXES.some((suffix) => lower.endsWith(suffix))) {
    return "Choose an assembled .fna, .fa, or .fasta genome. Read files such as FASTQ are outside this prototype's scope.";
  }
  if (file.size > 25 * 1024 * 1024) return "The FASTA exceeds the 25 MiB API limit.";
  if (file.size === 0) return "The selected FASTA is empty.";
  return null;
}

function UploadScreen() {
  const navigate = useNavigate();
  const fileInput = useRef<HTMLInputElement>(null);
  const [sampleName, setSampleName] = useState("");
  const [file, setFile] = useState<File | null>(null);
  const [phase, setPhase] = useState<Phase>("idle");
  const [dragOver, setDragOver] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const selectFile = (candidate?: File) => {
    if (!candidate) return;
    const validationError = validateFile(candidate);
    setError(validationError);
    setFile(validationError ? null : candidate);
    if (!validationError) {
      setSampleName((current) => current || sampleIdFromFileName(candidate.name));
    }
    setPhase(validationError ? "error" : "idle");
  };

  const startAnalysis = async () => {
    if (!file || !sampleName.trim()) return;
    setError(null);
    setPhase("analyzing");
    try {
      const report = await analyzeGenome(file);
      storeAnalysis(report, sampleName.trim());
      await navigate({
        to: "/dashboard",
        search: { analysisId: report.analysis_id },
      });
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : "Analysis failed unexpectedly.");
      setPhase("error");
    }
  };

  return (
    <div className="relative min-h-screen overflow-hidden bg-background">
      <DnaBackground />

      <header className="relative z-10 flex items-center justify-between border-b border-border/60 px-5 py-5 sm:px-8">
        <div className="flex items-center gap-2">
          <Dna className="h-5 w-5 text-foreground" />
          <span className="font-mono-tabular text-sm font-semibold tracking-[0.12em]">
            GyraseX
          </span>
        </div>
        <div className="hidden font-mono-tabular text-[11px] tracking-wider uppercase text-muted-foreground sm:block">
          Defensive AMR decision support
        </div>
      </header>

      <main className="relative z-10 mx-auto flex min-h-[calc(100vh-72px)] max-w-3xl flex-col items-center justify-center gap-5 px-5 py-12">
        <div className="w-full animate-fade-in-up rounded-md border border-border bg-card/95 p-6 shadow-[0_1px_0_rgba(0,0,0,0.02),0_20px_50px_-20px_rgba(15,23,42,0.15)] backdrop-blur-sm sm:p-8">
          <div className="mb-6 flex items-start justify-between gap-6">
            <div>
              <div className="font-mono-tabular text-[10px] tracking-[0.2em] uppercase text-muted-foreground">
                Genome analysis
              </div>
              <h1 className="mt-1 text-2xl font-semibold tracking-tight">
                Predict antibiotic response
              </h1>
              <p className="mt-2 max-w-xl text-sm leading-6 text-muted-foreground">
                Upload one quality-checked, assembled{" "}
                <span className="italic">Staphylococcus aureus</span> genome. The pipeline combines
                AMRFinderPlus features with a target-presence check.
              </p>
            </div>
            <Dna className="mt-1 hidden h-8 w-8 text-muted-foreground sm:block" />
          </div>

          <label className="block">
            <span className="font-mono-tabular text-[11px] tracking-wider uppercase text-muted-foreground">
              Sample label
            </span>
            <Input
              value={sampleName}
              onChange={(event) => setSampleName(event.target.value)}
              placeholder="e.g. SAMPLE-000342"
              disabled={phase === "analyzing"}
              className="mt-1.5 h-11 rounded-md font-mono-tabular text-sm"
            />
          </label>

          <div
            onDragOver={(event) => {
              event.preventDefault();
              setDragOver(true);
            }}
            onDragLeave={() => setDragOver(false)}
            onDrop={(event) => {
              event.preventDefault();
              setDragOver(false);
              selectFile(event.dataTransfer.files?.[0]);
            }}
            className={`mt-4 rounded-md border border-dashed p-6 text-center transition-colors ${
              dragOver ? "border-foreground bg-accent" : "border-border bg-muted/40"
            }`}
          >
            <div className="mx-auto mb-3 flex h-10 w-10 items-center justify-center rounded-md bg-background">
              {file ? <FileText className="h-4 w-4" /> : <Upload className="h-4 w-4" />}
            </div>
            <div className="text-sm text-muted-foreground">
              {file ? (
                <>
                  <span className="font-mono-tabular text-foreground">{file.name}</span>
                  <span className="ml-2 font-mono-tabular text-[11px]">
                    ({(file.size / 1024 / 1024).toFixed(2)} MiB)
                  </span>
                </>
              ) : (
                <>
                  Drag an assembled genome here ·{" "}
                  <span className="font-mono-tabular">.fna / .fa / .fasta</span>
                </>
              )}
            </div>
            <input
              ref={fileInput}
              type="file"
              accept=".fna,.fa,.fasta"
              className="hidden"
              onChange={(event) => selectFile(event.target.files?.[0])}
              disabled={phase === "analyzing"}
            />
            <Button
              type="button"
              variant="secondary"
              onClick={() => fileInput.current?.click()}
              disabled={phase === "analyzing"}
              className="mt-3 h-9 rounded-md font-mono-tabular text-[11px] tracking-wider uppercase"
            >
              Choose FASTA
            </Button>
          </div>

          <Button
            onClick={startAnalysis}
            disabled={!file || !sampleName.trim() || phase === "analyzing"}
            className="mt-5 h-11 w-full rounded-md font-mono-tabular text-xs tracking-[0.2em] uppercase"
          >
            {phase === "analyzing" ? (
              <>
                <Loader2 className="mr-2 h-4 w-4 animate-spin" /> Running M1 + M2 analysis
              </>
            ) : (
              "Analyze genome"
            )}
          </Button>
        </div>

        {phase === "analyzing" && (
          <div className="w-full animate-fade-in-up rounded-md border border-border bg-card/95 p-5 backdrop-blur-sm">
            <div className="flex items-center gap-3">
              <Loader2 className="h-4 w-4 animate-spin" />
              <div className="flex-1">
                <div className="font-mono-tabular text-[11px] tracking-wider uppercase text-muted-foreground">
                  Pipeline active
                </div>
                <div className="text-sm font-medium">
                  Quality checking, annotating resistance, and verifying drug targets…
                </div>
              </div>
            </div>
            <div className="mt-3 h-1 overflow-hidden rounded-full bg-muted">
              <div className="h-full w-2/3 animate-pulse bg-foreground" />
            </div>
            <p className="mt-3 text-xs text-muted-foreground">
              A bacterial genome can take a little while to process locally. Keep this page open.
            </p>
          </div>
        )}

        {error && (
          <Alert className="w-full animate-fade-in-up rounded-md border-danger/40 bg-card/95">
            <AlertTriangle className="h-4 w-4 text-danger" />
            <AlertTitle>Analysis could not start</AlertTitle>
            <AlertDescription>{error}</AlertDescription>
          </Alert>
        )}

        <div className="grid w-full gap-3 text-xs text-muted-foreground sm:grid-cols-2">
          <div className="flex gap-2 rounded-md border border-border bg-card/70 p-3">
            <CheckCircle2 className="mt-0.5 h-4 w-4 shrink-0 text-success" />
            <span>
              Supported scope: assembled, isolated <span className="italic">S. aureus</span> genomes
              only.
            </span>
          </div>
          <div className="flex gap-2 rounded-md border border-border bg-card/70 p-3">
            <AlertTriangle className="mt-0.5 h-4 w-4 shrink-0 text-warning" />
            <span>
              Every prediction must be confirmed by standard laboratory susceptibility testing.
            </span>
          </div>
        </div>
      </main>
    </div>
  );
}
