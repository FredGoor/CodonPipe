"""Download the coding-sequence FASTA files CodonPipe ships with.

The genomes are not stored in the repository: 34 CDS FASTA files come to about
230 MB, which is not what version control is for. They are listed in
``genomes.tsv`` by organism and INSDC accession, and this script retrieves them
from the NCBI Nucleotide service into this directory under the file names the
interface expects.

Usage, from the repository root::

    python "Preloaded genomes/fetch_genomes.py"                 # everything missing
    python "Preloaded genomes/fetch_genomes.py" --list          # show what would be fetched
    python "Preloaded genomes/fetch_genomes.py" Salmonella      # only matching organisms
    python "Preloaded genomes/fetch_genomes.py" --email you@example.org

NCBI asks automated clients for a contact address; pass one with --email if you
run large batches. One genome, Chrysiobacterium nematophagum, is an in-house
assembly with no public accession and is shipped with the repository.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import time
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
MANIFEST = os.path.join(HERE, "genomes.tsv")
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
LOCAL_ONLY = "in-house assembly"


def read_manifest():
    with open(MANIFEST, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fetch_one(accessions, email=None, tries=3):
    """Return the CDS FASTA of one or more replicons as a single string."""
    query = {"db": "nuccore", "id": ",".join(accessions),
             "rettype": "fasta_cds_na", "retmode": "text", "tool": "codonpipe"}
    if email:
        query["email"] = email
    url = EFETCH + "?" + urllib.parse.urlencode(query)
    for attempt in range(1, tries + 1):
        try:
            with urllib.request.urlopen(url, timeout=600) as response:
                text = response.read().decode("utf-8", "replace")
            if text.lstrip().startswith(">"):
                return text
            raise RuntimeError(text.strip()[:200] or "empty response")
        except Exception as exc:
            if attempt == tries:
                raise
            print(f"    attempt {attempt} failed ({exc}); retrying", flush=True)
            time.sleep(5 * attempt)
    raise RuntimeError("unreachable")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("pattern", nargs="*", help="only organisms containing these words")
    parser.add_argument("--email", default=os.environ.get("NCBI_EMAIL", ""),
                        help="contact address sent with each request")
    parser.add_argument("--list", action="store_true", help="list without downloading")
    parser.add_argument("--force", action="store_true", help="re-download existing files")
    args = parser.parse_args(argv)

    rows = read_manifest()
    if args.pattern:
        needle = [p.lower() for p in args.pattern]
        rows = [r for r in rows if any(p in r["organism"].lower() for p in needle)]

    todo, skipped = [], 0
    for row in rows:
        target = os.path.join(HERE, row["file"])
        if row["accessions"].startswith(LOCAL_ONLY):
            continue
        if os.path.exists(target) and not args.force:
            skipped += 1
            continue
        todo.append(row)

    print(f"{len(rows)} genomes listed, {skipped} already present, {len(todo)} to fetch")
    if args.list or not todo:
        for row in todo:
            print(f"  {row['organism']:<40s} {row['accessions']}")
        return 0

    failures = []
    for i, row in enumerate(todo, start=1):
        print(f"[{i}/{len(todo)}] {row['organism']}", flush=True)
        try:
            text = fetch_one(row["accessions"].split(";"), email=args.email or None)
        except Exception as exc:
            print(f"    FAILED: {exc}")
            failures.append(row["organism"])
            continue
        with open(os.path.join(HERE, row["file"]), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
        n = text.count(">")
        expected = int(row["CDS"]) if row["CDS"].isdigit() else 0
        flag = "" if not expected or abs(n - expected) <= max(5, 0.02 * expected) else \
               f"  WARNING: {n} CDS, manifest says {expected}"
        print(f"    {n} CDS written{flag}")
        time.sleep(1.0)

    if failures:
        print("\nnot retrieved: " + ", ".join(failures))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
