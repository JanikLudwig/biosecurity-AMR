#!/usr/bin/env bash
set -e

# Hinweis: Dieses Skript ist für Linux, macOS, WSL oder Git Bash gedacht.
# Es wird für Windows-Nutzer nicht zwingend vorausgesetzt.


# Verzeichnisse anlegen falls nicht vorhanden
mkdir -p data/raw/
mkdir -p reports/

echo "Starte Download der BV-BRC Metadaten (Release Files)..."

# URLs
AMR_URL="ftps://ftp.bvbrc.org/RELEASE_NOTES/PATRIC_genome_AMR.txt"
METADATA_URL="ftps://ftp.bvbrc.org/RELEASE_NOTES/genome_metadata"

# Zieldateien
AMR_TARGET="data/raw/PATRIC_genome_AMR.txt"
METADATA_TARGET="data/raw/genome_metadata"

function download_ftps {
    local url=$1
    local target=$2
    local part="${target}.part"

    if [ -f "$target" ]; then
        echo "Datei $target existiert bereits. Ãœberspringe Download."
        return 0
    fi

    echo "Lade $url nach $part herunter..."
    # --ftp-ssl requires TLS, --retry for retries, -L to follow redirects (though ftp doesn't), -o output
    curl --ftp-ssl -s -S -f --retry 3 --retry-connrefused -o "$part" "$url"
    
    if [ $? -eq 0 ]; then
        mv "$part" "$target"
        echo "Erfolgreich nach $target verschoben."
    else
        echo "Fehler beim Download von $url."
        rm -f "$part"
        exit 1
    fi
}

download_ftps "$AMR_URL" "$AMR_TARGET"
download_ftps "$METADATA_URL" "$METADATA_TARGET"

echo "Downloads abgeschlossen."
