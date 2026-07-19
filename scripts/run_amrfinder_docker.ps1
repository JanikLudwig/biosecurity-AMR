<#
.SYNOPSIS
Wrapper script to run AMRFinderPlus using Docker backend via Python runner.

.DESCRIPTION
This script acts as a thin wrapper around the Python CLI `python -m genome_firewall.run_amrfinder`.
It enforces the Docker backend and passes along optional parameters.
#>

[CmdletBinding()]
Param(
    [string]$InputDir,
    [string]$OutputDir,
    [string]$LogDir,
    [int]$Threads,
    [int]$Workers,
    [string]$Image,
    [string]$Organism,
    [string]$Manifest,
    [switch]$Force,
    [int]$Limit,
    [switch]$FailFast,
    [switch]$Plus
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

if (-not (Get-Command "python" -ErrorAction SilentlyContinue)) {
    Write-Error "Python is not available in PATH. Please activate the conda environment."
}

$pythonArgs = @("-m", "genome_firewall.run_amrfinder", "--backend", "docker")

if ($InputDir) { $pythonArgs += "--input-dir", $InputDir }
if ($OutputDir) { $pythonArgs += "--output-dir", $OutputDir }
if ($LogDir) { $pythonArgs += "--log-dir", $LogDir }
if ($Threads) { $pythonArgs += "--threads", $Threads.ToString() }
if ($Workers) { $pythonArgs += "--workers", $Workers.ToString() }
if ($Image) { $pythonArgs += "--image", $Image }
if ($Organism) { $pythonArgs += "--organism", $Organism }
if ($Manifest) { $pythonArgs += "--manifest", $Manifest }
if ($Force) { $pythonArgs += "--force" }
if ($Limit) { $pythonArgs += "--limit", $Limit.ToString() }
if ($FailFast) { $pythonArgs += "--fail-fast" }
if ($Plus) { $pythonArgs += "--plus" }

Write-Host "Executing: python $($pythonArgs -join ' ')"
$process = Start-Process -FilePath "python" -ArgumentList $pythonArgs -NoNewWindow -Wait -PassThru

if ($process.ExitCode -ne 0) {
    Write-Error "Python runner failed with exit code $($process.ExitCode)"
}
exit $process.ExitCode
