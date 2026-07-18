import hashlib
import logging
import urllib.request
import urllib.error
import urllib.parse
from pathlib import Path
import time
import ssl

logger = logging.getLogger(__name__)

def compute_sha256(file_path: Path) -> str:
    """Computes SHA-256 for a given file."""
    sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        for block in iter(lambda: f.read(65536), b""):
            sha256.update(block)
    return sha256.hexdigest()

def robust_ftps_download(url: str, output_path: Path, max_retries: int = 3, timeout: int = 30) -> bool:
    """
    Lade eine Datei robust per FTPS (oder FTP/HTTP) herunter.
    Unterstützt idempotenten Download, Timeout, Retries und atomare Umbenennung (.part).
    """
    if output_path.exists() and output_path.stat().st_size > 0:
        logger.info(f"Datei existiert bereits: {output_path}")
        return True

    part_file = output_path.with_suffix(".part")
    
    # Erlaube FTPS über unbestätigte SSL Zertifikate, da manche ftps-server Probleme machen können, 
    # aber BV-BRC sollte SSL können. Wir nutzen den Default Context.
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE

    for attempt in range(1, max_retries + 1):
        try:
            logger.info(f"Download (Versuch {attempt}): {url}")
            # BV-BRC blockt evtl. User-Agents nicht, aber zur Sicherheit
            req = urllib.request.Request(url, headers={'User-Agent': 'GenomeFirewall/0.1'})
            
            with urllib.request.urlopen(req, context=ctx, timeout=timeout) as response:
                with open(part_file, "wb") as f_out:
                    # Chunks lesen
                    while True:
                        chunk = response.read(8192)
                        if not chunk:
                            break
                        f_out.write(chunk)
                        
            # Wenn fertig, umbenennen (atomar)
            part_file.replace(output_path)
            logger.info(f"Erfolgreich heruntergeladen: {output_path}")
            return True
            
        except (urllib.error.URLError, OSError) as e:
            logger.warning(f"Download-Fehler: {e}")
            if part_file.exists():
                part_file.unlink()
            if attempt < max_retries:
                time.sleep(2 ** attempt)  # Exponential Backoff
            else:
                logger.error(f"Alle {max_retries} Versuche fehlgeschlagen für {url}.")
                return False

    return False
