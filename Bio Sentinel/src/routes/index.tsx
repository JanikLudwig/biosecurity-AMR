import { createFileRoute, useNavigate } from "@tanstack/react-router";
import { useEffect, useMemo, useState, type ChangeEvent, type DragEvent } from "react";
import { DnaBackground } from "@/components/dna-background";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { Upload, Loader2, CheckCircle2, Dna, AlertTriangle } from "lucide-react";
import {
  GenomeFirewallApiError,
  getGenomeReport,
  listGenomes,
  sampleIdFromFileName,
  type GenomeIndexEntry,
  type GenomeReport,
} from "@/lib/genome-firewall-api";

export const Route = createFileRoute("/")({
  component: UploadScreen,
});

type Phase = "idle" | "loading" | "revealed";

function UploadScreen() {
  const navigate = useNavigate();
  const [sampleName, setSampleName] = useState("");
  const [fileName, setFileName] = useState<string | null>(null);
  const [phase, setPhase] = useState<Phase>("idle");
  const [dragOver, setDragOver] = useState(false);
  const [genomes, setGenomes] = useState<GenomeIndexEntry[]>([]);
  const [catalogError, setCatalogError] = useState<string | null>(null);
  const [report, setReport] = useState<GenomeReport | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    const controller = new AbortController();
    listGenomes(controller.signal)
      .then(setGenomes)
      .catch((reason: unknown) => {
        if ((reason as DOMException)?.name !== "AbortError") {
          setCatalogError("The prepared-sample catalog is unavailable. You can still enter a known genome ID.");
        }
      });
    return () => controller.abort();
  }, []);

  const knownIds = useMemo(() => new Set(genomes.map((genome) => genome.genome_id)), [genomes]);

  const selectFile = (file: File | undefined) => {
    if (!file) return;
    setFileName(file.name);
    setSampleName((current) => current || sampleIdFromFileName(file.name));
    setError(null);
  };

  const onFileChange = (event: ChangeEvent<HTMLInputElement>) => {
    selectFile(event.target.files?.[0]);
  };

  const onDrop = (event: DragEvent<HTMLDivElement>) => {
    event.preventDefault();
    setDragOver(false);
    selectFile(event.dataTransfer.files?.[0]);
  };

  const startAnalysis = async () => {
    const genomeId = sampleName.trim() || (fileName ? sampleIdFromFileName(fileName) : "");
    if (!genomeId) return;

    setPhase("loading");
    setError(null);
    setReport(null);
    try {
      const nextReport = await getGenomeReport(genomeId);
      setReport(nextReport);
      setPhase("revealed");
    } catch (reason) {
      setPhase("idle");
      if (reason instanceof GenomeFirewallApiError) {
        setError(reason.status === 404
          ? `No registered genome report was found for "${genomeId}".`
          : reason.message);
      } else {
        setError("The genome prediction could not be loaded.");
      }
    }
  };

  const candidateId = sampleName.trim() || (fileName ? sampleIdFromFileName(fileName) : "");
  const isKnownSample = candidateId ? knownIds.has(candidateId) : false;

  return (
    <div className="relative min-h-screen overflow-hidden bg-background">
      <DnaBackground />

      <header className="relative z-10 flex items-center justify-between px-8 py-6">
        <div className="flex items-center gap-2">
          <Dna className="h-5 w-5 text-foreground" />
          <span className="font-mono-tabular text-sm font-semibold tracking-[0.2em] uppercase">
            GyraseX
          </span>
        </div>
        <div className="font-mono-tabular text-[11px] tracking-wider uppercase text-muted-foreground">
          Antibiotic susceptibility engine
        </div>
      </header>

      <main className="relative z-10 mx-auto flex min-h-[calc(100vh-88px)] max-w-2xl flex-col items-center justify-center gap-6 px-6 pb-24">
        <div className="w-full animate-fade-in-up rounded-md border border-border bg-card/95 p-6 shadow-[0_1px_0_rgba(0,0,0,0.02),0_20px_50px_-20px_rgba(15,23,42,0.15)] backdrop-blur-sm">
          <div className="mb-4">
            <div className="font-mono-tabular text-[10px] tracking-[0.2em] uppercase text-muted-foreground">
              Step 01 - Sample
            </div>
            <h1 className="mt-1 text-xl font-semibold tracking-tight">Initialize sequence</h1>
          </div>

          <label className="block">
            <span className="font-mono-tabular text-[11px] tracking-wider uppercase text-muted-foreground">
              Genome ID / Sample Name
            </span>
            <Input
              value={sampleName}
              onChange={(event) => setSampleName(event.target.value)}
              placeholder="e.g. 1280.10033"
              list="prepared-genomes"
              disabled={phase !== "idle"}
              className="mt-1.5 h-11 rounded-md font-mono-tabular text-sm"
            />
            <datalist id="prepared-genomes">
              {genomes.map((genome) => (
                <option key={genome.genome_id} value={genome.genome_id}>
                  {genome.mlst_group ? `${genome.mlst_group} - ${genome.partition || "prepared"}` : "prepared sample"}
                </option>
              ))}
            </datalist>
          </label>

          <div
            onDragOver={(event) => {
              event.preventDefault();
              setDragOver(true);
            }}
            onDragLeave={() => setDragOver(false)}
            onDrop={onDrop}
            className={`mt-4 rounded-md border border-dashed p-6 text-center transition-colors ${
              dragOver ? "border-foreground bg-accent" : "border-border bg-muted/40"
            }`}
          >
            <div className="mx-auto mb-3 flex h-10 w-10 items-center justify-center rounded-md bg-background">
              <Upload className="h-4 w-4" />
            </div>
            <div className="text-sm text-muted-foreground">
              {fileName ? (
                <span className="font-mono-tabular text-foreground">{fileName}</span>
              ) : (
                <>Select a registered FASTA filename or enter its genome ID</>
              )}
            </div>
            <label className="mt-3 inline-block">
              <input
                type="file"
                accept=".fasta,.fa,.fna,.gz"
                className="hidden"
                onChange={onFileChange}
                disabled={phase !== "idle"}
              />
              <span className="inline-flex cursor-pointer items-center gap-2 rounded-md bg-foreground px-4 py-2 text-xs font-medium tracking-wide text-background hover:opacity-90">
                <Upload className="h-3.5 w-3.5" />
                Genome seq Upload
              </span>
            </label>
          </div>

          <p className="mt-3 text-xs leading-relaxed text-muted-foreground">
            This deployment reads registered genomes through the GyraseX report API.
            Selecting a file fills its genome ID; raw file upload is not sent to the pipeline yet.
          </p>

          {candidateId && genomes.length > 0 && (
            <div className="mt-2 font-mono-tabular text-[10px] tracking-wide text-muted-foreground">
              {isKnownSample ? "Prepared report available" : "A live lookup will be attempted for this genome ID"}
            </div>
          )}
          {catalogError && <div className="mt-2 text-xs text-muted-foreground">{catalogError}</div>}

          <Button
            onClick={startAnalysis}
            disabled={!candidateId || phase !== "idle"}
            className="mt-5 h-11 w-full rounded-md font-mono-tabular text-xs tracking-[0.2em] uppercase"
          >
            Analyze Sequence
          </Button>
        </div>

        {error && (
          <Alert className="w-full rounded-md border-border bg-card/95 backdrop-blur-sm">
            <AlertTriangle className="h-4 w-4" style={{ color: "var(--danger)" }} />
            <AlertTitle className="font-mono-tabular text-[11px] tracking-wider uppercase">Analysis unavailable</AlertTitle>
            <AlertDescription>{error}</AlertDescription>
          </Alert>
        )}

        {phase !== "idle" && (
          <div className="w-full animate-fade-in-up rounded-md border border-border bg-card/95 p-5 backdrop-blur-sm">
            <div className="flex items-center gap-3">
              {phase === "loading" ? (
                <Loader2 className="h-4 w-4 animate-spin" />
              ) : (
                <CheckCircle2 className="h-4 w-4" style={{ color: "var(--success)" }} />
              )}
              <div className="flex-1">
                <div className="font-mono-tabular text-[11px] tracking-wider uppercase text-muted-foreground">
                  {phase === "loading" ? "Pipeline" : "Complete"}
                </div>
                <div className="text-sm font-medium">
                  {phase === "loading" ? "Loading prediction..." : "Prediction ready"}
                </div>
              </div>
              <div className="font-mono-tabular text-[11px] text-muted-foreground">
                {phase === "loading" ? "M1 - M5" : `${report?.decisions.length || 0} drugs`}
              </div>
            </div>
            <div className="mt-3 h-1 overflow-hidden rounded-full bg-muted">
              <div
                className="h-full transition-all duration-500"
                style={{ width: phase === "loading" ? "70%" : "100%", background: "var(--foreground)" }}
              />
            </div>
          </div>
        )}

        {phase === "revealed" && report && (
          <div className="w-full animate-fade-in-up rounded-md border border-border bg-card/95 p-5 backdrop-blur-sm">
            <div className="flex items-end justify-between gap-4">
              <div>
                <div className="font-mono-tabular text-[10px] tracking-[0.2em] uppercase text-muted-foreground">
                  Species name
                </div>
                <div className="mt-1 text-2xl font-semibold italic tracking-tight">{report.species}</div>
              </div>
              <Button
                onClick={() => navigate({ to: "/dashboard", search: { genomeId: report.genome_id } })}
                className="h-10 rounded-md font-mono-tabular text-[11px] tracking-[0.2em] uppercase"
              >
                Open Analysis -&gt;
              </Button>
            </div>
          </div>
        )}
      </main>
    </div>
  );
}
