/**
 * Frontend-owned adapter for the public GyraseX JSON contract.
 * It intentionally knows nothing about Python internals, keeping this branch
 * straightforward to merge with independent pipeline work.
 */

export type PredictionCall = "likely to work" | "likely to fail" | "no-call";

export interface GenomeIndexEntry {
  genome_id: string;
  mlst_group?: string;
  partition?: string;
  n_proteins?: number;
  summary?: Partial<Record<PredictionCall, number>>;
}

export interface TargetEvidence {
  gene: string;
  present: boolean;
  identity?: number;
  ref_coverage?: number;
  evalue?: number;
  bitscore?: number;
  orf_id?: string | null;
  contig?: string | null;
}

export interface AntibioticDecision {
  drug: string;
  drug_display: string;
  drug_class: string;
  tier: string;
  call: PredictionCall;
  confidence: number;
  evidence_category: string;
  target_status: "present" | "absent" | "not_applicable" | string;
  p_resistant: number | null;
  supporting_determinants: string[];
  target_evidence: TargetEvidence[];
  rationale: string;
  no_call_reason: string | null;
}

export interface GenomeReport {
  version: string;
  genome_id: string;
  species: string;
  scope_ok: boolean;
  qc: { passed?: boolean; flags?: string[] };
  n_proteins_predicted: number;
  features_synthetic: boolean;
  m1?: {
    feature_row_available?: boolean;
    nonzero_feature_count?: number;
    source?: string;
  };
  summary: Record<PredictionCall, number>;
  decisions: AntibioticDecision[];
  warnings: string[];
  safety_notice: string;
}

export interface ReliabilityData {
  mean_predicted: number[];
  observed_frequency: number[];
  count: number[];
}

export interface PerformanceMetric {
  drug: string;
  tier: string;
  n_test: number;
  auroc?: number;
  pr_auc?: number;
  balanced_accuracy?: number;
  recall_resistant?: number;
  f1?: number;
  brier?: number;
  reliability?: ReliabilityData;
}

interface GenomeIndexResponse {
  genomes: GenomeIndexEntry[];
  count: number;
}

interface MetricsResponse {
  metrics: PerformanceMetric[];
  synthetic_features: boolean;
}

const apiBase = (import.meta.env.VITE_GENOME_FIREWALL_API_BASE_URL || "/api").replace(/\/$/, "");

export class GenomeFirewallApiError extends Error {
  status?: number;

  constructor(message: string, status?: number) {
    super(message);
    this.name = "GenomeFirewallApiError";
    this.status = status;
  }
}

async function request<T>(path: string, signal?: AbortSignal): Promise<T> {
  let response: Response;
  try {
    response = await fetch(`${apiBase}${path}`, { signal });
  } catch (error) {
    // React development mode can intentionally cancel the first of two effect
    // runs. Preserve that signal so callers can ignore it instead of showing a
    // false API-connection error.
    if (typeof error === "object" && error !== null && "name" in error && error.name === "AbortError") {
      throw error;
    }
    throw new GenomeFirewallApiError(
      "Unable to reach the GyraseX API. Check the API base URL or local proxy configuration.",
    );
  }

  if (!response.ok) {
    let message = `GyraseX API request failed (${response.status}).`;
    try {
      const body = (await response.json()) as { detail?: string };
      if (body.detail) message = body.detail;
    } catch {
      // The HTTP status remains useful if the backend did not send JSON.
    }
    throw new GenomeFirewallApiError(message, response.status);
  }

  return response.json() as Promise<T>;
}

export async function listGenomes(signal?: AbortSignal): Promise<GenomeIndexEntry[]> {
  const response = await request<GenomeIndexResponse>("/genomes", signal);
  return response.genomes || [];
}

export function getGenomeReport(genomeId: string, signal?: AbortSignal): Promise<GenomeReport> {
  return request<GenomeReport>(`/report/${encodeURIComponent(genomeId)}`, signal);
}

export async function getPerformanceMetrics(signal?: AbortSignal): Promise<PerformanceMetric[]> {
  const response = await request<MetricsResponse>("/metrics", signal);
  return response.metrics || [];
}

export function sampleIdFromFileName(fileName: string): string {
  return fileName.replace(/\.(fasta|fna|fa)(\.gz)?$/i, "");
}
