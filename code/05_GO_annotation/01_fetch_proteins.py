import requests
import time
import sys
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

NCBI_BATCH_SIZE = 200
UNIPROT_BATCH_SIZE = 200

# NCBI: <= 3 req/s without API key, <= 10 req/s with API key
NCBI_SLEEP = 0.34
UNIPROT_SLEEP = 0.34

NCBI_API_KEY = None  # put the key here

INPUT_FILE = "protein_ids.txt"
OUTPUT_FASTA = "proteins.fa"
FAILED_IDS_FILE = "failed_ids.txt"


def build_session():
    retry = Retry(
        total=5,
        connect=5,
        read=5,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry)
    s = requests.Session()
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    s.headers.update({"User-Agent": "protein-fetcher/1.0"})
    return s


session = build_session()


def count_fasta_entries(text: str) -> int:
    return sum(1 for line in text.splitlines() if line.startswith(">"))


def show_status(fetched, total, failed):
    msg = f"\rProteins fetched ({fetched}/{total}) | Failed ({failed})"
    sys.stdout.write(msg)
    sys.stdout.flush()


def fetch_ncbi(ids):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "protein",
        "id": ",".join(ids),
        "rettype": "fasta",
        "retmode": "text",
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    r = session.get(url, params=params, timeout=120)
    r.raise_for_status()
    return r.text


def fetch_uniprot(ids):
    url = "https://rest.uniprot.org/uniprotkb/stream"
    query = "(" + " OR ".join(f"accession:{i}" for i in ids) + ")"
    params = {
        "format": "fasta",
        "query": query,
        "download": "true",
    }

    r = session.get(url, params=params, timeout=120)
    r.raise_for_status()
    return r.text


def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


with open(INPUT_FILE) as f:
    proteins = [line.strip() for line in f if line.strip()]

ncbi_ids = []
uniprot_ids = []

for p in proteins:
    if p.startswith("sp|") or p.startswith("tr|"):
        parts = p.split("|")
        if len(parts) > 1:
            uniprot_ids.append(parts[1])
        else:
            uniprot_ids.append(p)
    else:
        ncbi_ids.append(p)

total = len(ncbi_ids) + len(uniprot_ids)
fetched = 0
failed = 0
failed_ids = []
fasta_blocks = []

show_status(fetched, total, failed)

# NCBI
for chunk in chunks(ncbi_ids, NCBI_BATCH_SIZE):
    try:
        text = fetch_ncbi(chunk)
        n = count_fasta_entries(text)

        if text.strip():
            fasta_blocks.append(text.strip())

        fetched += n

        # IDs requested but not returned
        missing = len(chunk) - n
        if missing > 0:
            failed += missing
            # exact missing IDs are not trivial to map from efetch response alone
            # so we keep the whole chunk as "possibly incomplete" only if needed
            # comment this in if you want chunk-level tracking:
            # failed_ids.extend(chunk)

    except requests.RequestException:
        failed += len(chunk)
        failed_ids.extend(chunk)

    show_status(fetched, total, failed)
    time.sleep(NCBI_SLEEP)

# UniProt
for chunk in chunks(uniprot_ids, UNIPROT_BATCH_SIZE):
    try:
        text = fetch_uniprot(chunk)
        n = count_fasta_entries(text)

        if text.strip():
            fasta_blocks.append(text.strip())

        fetched += n

        missing = len(chunk) - n
        if missing > 0:
            failed += missing
            # optional exact tracking would need header parsing against accession IDs

    except requests.RequestException:
        failed += len(chunk)
        failed_ids.extend(chunk)

    show_status(fetched, total, failed)
    time.sleep(UNIPROT_SLEEP)

with open(OUTPUT_FASTA, "w") as f:
    if fasta_blocks:
        f.write("\n".join(fasta_blocks) + "\n")

if failed_ids:
    with open(FAILED_IDS_FILE, "w") as f:
        for pid in failed_ids:
            f.write(pid + "\n")

sys.stdout.write("\nDone.\n")