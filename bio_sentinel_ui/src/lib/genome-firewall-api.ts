export type DecisionCall = "likely_to_work" | "likely_to_fail" | "no_call";

export interface AnalysisDecision {
  antibiotic: string;
  call: DecisionCall;
  confidence: number;
  probability_resistant: number;
  model_signal: string;
  evidence_category: string;
  supporting_elements: string[];
  all_relevant_elements: string[];
  target_status: string;
  target_label: string;
  targets_detected: string[];
  feature_similarity: number;
  feature_similarity_floor: number;
  unknown_features: string[];
  lineage_status: string;
  maximum_training_ani: number | null;
  minimum_training_ani: number | null;
  nearest_training_genome: string | null;
  reasons: string[];
}

export interface M1Evidence {
  feature_key: string;
  element_symbol: string;
  element_name: string;
  evidence_category: string;
  amr_class: string;
  amr_subclass: string;
  method: string;
  coverage: number;
  identity: number;
}

export interface M2Evidence {
  symbol: string;
  present: boolean;
  reference?: string;
  contig?: string;
  orf_id?: string;
  reference_coverage?: number;
  identity?: number;
  method?: string;
}

export interface TargetResult {
  status: string;
  detected: string[];
  required_count: number;
  probe_kind: string;
  evidence: M2Evidence[];
}

export interface AnalysisReport {
  analysis_id: string;
  schema_version: string;
  generated_at: string;
  genome_id: string;
  species_scope: string;
  warning: string;
  defensive_use_only: boolean;
  qc: {
    passed: boolean;
    genome_length: number;
    contigs: number;
    contig_n50: number;
    ambiguous_fraction: number;
    reasons: string[];
  };
  provenance: {
    amrfinder_version: string;
    amrfinder_database: string;
    drug_registry_version: string;
    model_directory: string;
  };
  workflows: {
    M1: {
      name: string;
      recognized_features: string[];
      unknown_features: string[];
      evidence: M1Evidence[];
    };
    M2: {
      workflow: string;
      method: string;
      predicted_proteins: number;
      drugs: Record<string, TargetResult>;
    };
  };
  lineage: {
    status: string;
    maximum_training_ani: number | null;
    minimum_training_ani: number | null;
    nearest_training_genome: string | null;
  };
  decisions: AnalysisDecision[];
}

export interface DrugEvaluation {
  calibration_status: string;
  no_call_rate: number;
  test_metrics: {
    auroc: number;
    balanced_accuracy: number;
    brier: number;
    f1: number;
    pr_auc: number;
    resistant_recall: number;
    susceptible_recall: number;
  };
  called_test_metrics: {
    accuracy: number;
    called_count: number;
    coverage: number;
  };
  class_counts: Record<string, { resistant: number; susceptible: number }>;
}

export interface ModelsResponse {
  species: string;
  antibiotics: string[];
  feature_count: number;
  feature_schema_sha256: string;
  warning: string;
  bundle: {
    evaluation_status: string;
    split_genomes: number;
    models: Record<
      string,
      {
        evaluation: DrugEvaluation;
        decision_thresholds: {
          lower: number | null;
          upper: number | null;
          susceptible_calls: number;
          resistant_calls: number;
          susceptible_call_error: number | null;
          resistant_call_error: number | null;
        };
      }
    >;
    reliability: Record<
      string,
      Array<{
        probability_bin: number;
        samples: number;
        mean_probability_resistant: number;
        observed_resistant_fraction: number;
      }>
    >;
  };
}

const API_BASE = (import.meta.env.VITE_API_BASE_URL ?? "http://127.0.0.1:8000").replace(/\/$/, "");
const ANALYSIS_KEY = "genome-firewall:analysis-v1";
const SAMPLE_KEY = "genome-firewall:sample-name";

async function apiError(response: Response): Promise<Error> {
  try {
    const body = (await response.json()) as { detail?: string | unknown[] };
    const detail = typeof body.detail === "string" ? body.detail : JSON.stringify(body.detail);
    return new Error(detail || `Request failed (${response.status})`);
  } catch {
    return new Error(`Request failed (${response.status})`);
  }
}

export async function analyzeGenome(file: File): Promise<AnalysisReport> {
  const form = new FormData();
  form.append("fasta", file);
  const response = await fetch(`${API_BASE}/api/v1/analyses`, {
    method: "POST",
    body: form,
  });
  if (!response.ok) throw await apiError(response);
  return (await response.json()) as AnalysisReport;
}

export async function fetchAnalysis(
  analysisId: string,
  signal?: AbortSignal,
): Promise<AnalysisReport> {
  const response = await fetch(`${API_BASE}/api/v1/analyses/${encodeURIComponent(analysisId)}`, {
    signal,
  });
  if (!response.ok) throw await apiError(response);
  return (await response.json()) as AnalysisReport;
}

export async function fetchModels(signal?: AbortSignal): Promise<ModelsResponse> {
  const response = await fetch(`${API_BASE}/api/v1/models`, { signal });
  if (!response.ok) throw await apiError(response);
  return (await response.json()) as ModelsResponse;
}

export function storeAnalysis(report: AnalysisReport, sampleName: string): void {
  sessionStorage.setItem(ANALYSIS_KEY, JSON.stringify(report));
  sessionStorage.setItem(SAMPLE_KEY, sampleName);
}

export function loadStoredAnalysis(): { report: AnalysisReport; sampleName: string } | null {
  const raw = sessionStorage.getItem(ANALYSIS_KEY);
  if (!raw) return null;
  try {
    return {
      report: JSON.parse(raw) as AnalysisReport,
      sampleName: sessionStorage.getItem(SAMPLE_KEY) || "Unnamed sample",
    };
  } catch {
    sessionStorage.removeItem(ANALYSIS_KEY);
    return null;
  }
}

export function displayDrugName(value: string): string {
  return value.charAt(0).toUpperCase() + value.slice(1);
}

export function sampleIdFromFileName(fileName: string): string {
  return fileName.replace(/\.(fasta|fna|fa)$/i, "");
}
