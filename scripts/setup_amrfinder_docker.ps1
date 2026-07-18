<#
.SYNOPSIS
Validates Docker and AMRFinderPlus setup for Genome Firewall.

.DESCRIPTION
Checks if Docker is running, pulls the pinned AMRFinderPlus image,
and verifies if the required organism (Staphylococcus_aureus) is supported.
Reports are written to the reports/ directory.
#>

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$PSScriptRoot = Split-Path -Parent -Path $MyInvocation.MyCommand.Definition
$RepoRoot = Resolve-Path (Join-Path $PSScriptRoot "..")
$ReportsDir = Join-Path $RepoRoot "reports"

if (-not (Test-Path $ReportsDir)) {
    New-Item -ItemType Directory -Path $ReportsDir | Out-Null
}

$TargetImage = "ncbi/amr:4.2.7-2026-05-15.1"
$ExpectedOrganism = "Staphylococcus_aureus"

function Invoke-NativeDocker {
    param(
        [Parameter(Mandatory=$true, ValueFromRemainingArguments=$true)]
        [string[]]$CommandArgs
    )
    $oldPreference = $global:ErrorActionPreference
    $global:ErrorActionPreference = "Continue"
    try {
        $output = & docker $CommandArgs 2>&1
        $exitCode = $LASTEXITCODE
        $outputStr = ($output | Out-String).Trim()
        if ($exitCode -ne 0) {
            throw "Native docker command failed with exit code $exitCode.`nCommand: docker $($CommandArgs -join ' ')`nOutput: $outputStr"
        }
        return $outputStr
    } finally {
        $global:ErrorActionPreference = $oldPreference
    }
}

Write-Host "Checking if docker is available..."
if (-not (Get-Command "docker" -ErrorAction SilentlyContinue)) {
    Write-Error "Docker is not available in PATH. Please install Docker Desktop and ensure it is in your PATH."
}

Write-Host "Checking Docker daemon status..."
try {
    $dockerInfo = docker info 2>&1
    if ($LASTEXITCODE -ne 0) {
        Write-Error "Docker daemon is not running. Please start Docker Desktop.`nDetails: $dockerInfo"
    }
} catch {
    Write-Error "Docker daemon is not running or accessible. Please start Docker Desktop."
}

Write-Host "Ensuring image $TargetImage is available locally..."
$imageCheck = docker images -q $TargetImage
if (-not $imageCheck) {
    Write-Host "Image not found locally. Pulling $TargetImage ..."
    docker pull $TargetImage
    if ($LASTEXITCODE -ne 0) {
        Write-Error "Failed to pull image $TargetImage."
    }
} else {
    Write-Host "Image $TargetImage is already present locally."
}

Write-Host "Fetching AMRFinderPlus version from container..."
$versionOutput = Invoke-NativeDocker -CommandArgs @("run", "--rm", "--entrypoint", "amrfinder", $TargetImage, "-V")
$versionFile = Join-Path $ReportsDir "amrfinder_version.txt"
$versionOutput | Out-File -FilePath $versionFile -Encoding utf8
Write-Host "Version output saved to $versionFile"

Write-Host "Fetching supported organisms from container..."
$organismsOutput = Invoke-NativeDocker -CommandArgs @("run", "--rm", "--entrypoint", "amrfinder", $TargetImage, "-l")
$organismsFile = Join-Path $ReportsDir "amrfinder_organisms.txt"
$organismsOutput | Out-File -FilePath $organismsFile -Encoding utf8
Write-Host "Organisms output saved to $organismsFile"

$organismsText = $organismsOutput -join "`n"
if ($organismsText -notmatch "\b$ExpectedOrganism\b") {
    Write-Error "The required organism '$ExpectedOrganism' was not found in the supported organisms list."
} else {
    Write-Host "Validation successful: Organism '$ExpectedOrganism' is supported."
}

Write-Host "Collecting Docker container metadata..."
$timestamp = (Get-Date).ToString("yyyy-MM-ddTHH:mm:ssK")
$imageId = docker image inspect $TargetImage --format '{{.Id}}'
$repoDigests = docker image inspect $TargetImage --format '{{json .RepoDigests}}'
$osType = docker info --format '{{.OSType}}'
$arch = docker info --format '{{.Architecture}}'

$containerReport = @"
Setup Timestamp: $timestamp
Image Name: ncbi/amr
Image Tag: 4.2.7-2026-05-15.1
Image ID: $imageId
Repo Digests: $repoDigests
Docker Server OSType: $osType
Docker Architecture: $arch
"@

$containerFile = Join-Path $ReportsDir "amrfinder_container.txt"
$containerReport | Out-File -FilePath $containerFile -Encoding utf8
Write-Host "Container metadata saved to $containerFile"

Write-Host "Docker setup and validation completed successfully."
