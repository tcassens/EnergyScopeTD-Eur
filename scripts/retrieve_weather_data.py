"""Command line helper for downloading ERA5 cutout files.

The module provides a small CLI that accepts a single year and saves the
corresponding PyPSA-Eur cutout from the Sigma2 archive into a local directory.
Downloading is handled with a simple progress indicator,
and existing downloads are skipped.
"""

import argparse
import os
import sys
from urllib.parse import urlparse
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError

FILE_URL_TEMPLATE = (
    "https://data.archive.sigma2.no/dataset/"
    "5c1be6cd-71a1-45a7-8b78-5f18bbbe5728/download/"
    "europe-era5_{year}.nc"
)
DEFAULT_OUTPUT_DIR = "cutouts"

# Function to convert bytes to a human-readable format
def human_readable_size(num_bytes: int) -> str:
    step = 1024.0
    units = ["B", "KB", "MB", "GB", "TB", "PB"]
    num = float(num_bytes)
    for unit in units:
        if num < step or unit == units[-1]:
            return f"{num:.1f} {unit}"
        num /= step

# Function to create a safe filename from a URL
def safe_filename_from_url(url: str) -> str:
    path = urlparse(url).path
    name = os.path.basename(path)
    return name or "downloaded_file"

# Function to download a file with progress indication
def download_file(url: str, out_dir: str):
    os.makedirs(out_dir, exist_ok=True)
    filename = safe_filename_from_url(url)
    out_path = os.path.join(out_dir, filename)

    if os.path.exists(out_path):
        print(f"[+] Skipping (already exists): {out_path}")
        return

    print(f"[*] Downloading:\n    {url}\n    -> {out_path}")

    req = Request(url, headers={"User-Agent": "pypsa-eur-retrieve-weather-data"})
    try:
        response = urlopen(req)
    except HTTPError as err:
        raise RuntimeError(f"HTTP error {err.code}: {err.reason}") from err
    except URLError as err:
        raise RuntimeError(f"Failed to reach server: {err.reason}") from err

    with response:
        total_size = int(response.headers.get("Content-Length") or 0)
        if total_size:
            print(f"    Size: {human_readable_size(total_size)}")

        downloaded = 0
        chunk_size = 1024 * 1024
        progress_width = 40

        def render_progress(current: int):
            if total_size:
                filled = int(progress_width * current / total_size)
                bar = "#" * filled + "-" * (progress_width - filled)
                percent = (current / total_size) * 100
                status = (
                    f"[{bar}] {percent:6.2f}% "
                    f"({human_readable_size(current)}/"
                    f"{human_readable_size(total_size)})"
                )
            else:
                status = f"Downloaded {human_readable_size(current)}"
            print(f"\r    {status}", end="", flush=True)

        with open(out_path, "wb") as f:
            while True:
                chunk = response.read(chunk_size)
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)
                render_progress(downloaded)

        if downloaded:
            print()
    print(f"[+] Done: {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Download ERA5 PyPSA-Eur cutout file for a given year."
    )
    parser.add_argument("year", type=int, help="Year to download (e.g. 2017)")
    parser.add_argument(
        "-o",
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    args = parser.parse_args()

    url = FILE_URL_TEMPLATE.format(year=args.year)

    try:
        download_file(url, args.output_dir)
    except Exception as e:
        print(f"[!] Error downloading {url}: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
