# Fetching from UniProt — project-scoped cheat sheet

What this project actually needs from UniProt: two endpoints, two HTTP calls per organism. Extracted from `R/uniprot.R`.

## What we use

| Job                          | Endpoint                                            |
|------------------------------|-----------------------------------------------------|
| Resolve `taxid → UPID`       | `https://rest.uniprot.org/proteomes/search`         |
| Download a proteome's FASTA  | `https://rest.uniprot.org/uniprotkb/stream`         |

Always use `httr2`, not `httr`. The R `quarto` and `uniprot` packages are *not* needed.

---

## 1. Proteome lookup (taxid → UPID)

```r
req <- request("https://rest.uniprot.org/proteomes/search") |>
  req_url_query(
    query  = paste0("taxonomy_id:", taxid),    # see warning below
    format = "json",
    fields = "upid,organism,organism_id,protein_count",
    size   = 25L
  ) |>
  req_user_agent("project-name/0.1 (R; httr2)") |>
  req_retry(max_tries = 4, backoff = function(i) 1.5^i)

resp <- req_perform(req)
body <- resp_body_json(resp)
results <- body$results        # list of zero or more proteome hits
```

### Use `taxonomy_id:` not `organism_id:`

- `organism_id:<taxid>` matches **exactly** that taxid. Misses Reference proteomes that UniProt indexes under a strain-level child of TEMPURA's species taxid.
- `taxonomy_id:<taxid>` does **NCBI tree descent** — matches the taxid or any descendant. Use this.

In the reference run, switching from `organism_id:` to `taxonomy_id:` halved the no-match rate and recovered Reference-tier proteomes.

### Valid `fields` values

Only these are valid as `fields` for the proteomes endpoint:

```
upid, organism, organism_id, components, mnemonic, lineage,
busco, cpd, genome_assembly, genome_representation, protein_count
```

`proteome_type` is **not** in this list — but it always comes back in the JSON automatically as `proteomeType`. Don't ask for it explicitly; the request will be rejected.

### JSON shape per result

```r
results[[1]]$id                       # the UPID, e.g. "UP000619545"
results[[1]]$proteomeType             # "Reference and representative proteome" | ...
results[[1]]$taxonomy$taxonId         # actual organism taxid (may be a strain child)
results[[1]]$taxonomy$scientificName  # organism name
results[[1]]$proteinCount             # n proteins in the proteome
```

### Proteome types — what to keep

| `proteomeType`                          | Use? | Rank |
|-----------------------------------------|------|------|
| `Reference and representative proteome` | yes  | 1    |
| `Reference proteome`                    | yes  | 1    |
| `Representative proteome`               | yes  | 2    |
| `Other proteome`                        | yes  | 3    |
| `Redundant proteome`                    | no   | —    |
| `Excluded proteome` / `Excluded`        | no   | —    |

When a single taxid query returns multiple proteomes (typical when a strain has both a Reference and a species-level Other), pick by **type rank**, ties broken by `proteinCount` descending.

In current (2026) UniProt data, plain `Reference proteome` is rare — most curated genomes are now `Reference and representative proteome` or `Representative proteome`. Accept the broader set.

---

## 2. FASTA download

```r
req <- request("https://rest.uniprot.org/uniprotkb/stream") |>
  req_url_query(
    query      = paste0("proteome:", upid),
    format     = "fasta",
    compressed = "true"
  ) |>
  req_user_agent("project-name/0.1") |>
  req_retry(max_tries = 3, backoff = function(i) 2 * i) |>
  req_timeout(300)

resp <- req_perform(req, path = "data/proteomes/<taxid>.fasta.gz")
```

- **`path = ...` streams the body straight to disk** without holding it in memory. Important for large proteomes.
- **`compressed=true`** makes the response *body* a gzip file (HTTP `Content-Encoding` is identity). Save the bytes verbatim — do **not** decompress. The file is gzip-magic on the first two bytes (`1f 8b`).
- File extension `.fasta.gz` matches the contents.

### Restartability

Before the call, check `file.exists(dest) && file.size(dest) > 0`. If yes, skip — the file is already on disk. No HEAD request, no caching headers, just trust the local file. This makes re-runs cost nothing.

---

## 3. Politeness and reliability

Three things to set on every request:

| Setting                          | Why                                                  |
|----------------------------------|------------------------------------------------------|
| `req_user_agent("project/0.1")`  | UniProt asks for one; helps them help you in errors  |
| `req_retry(max_tries = 3-4)`     | Handles transient 429 / 5xx                          |
| `Sys.sleep(0.3)` between calls   | Don't hammer the API                                 |

`req_retry()` honors `Retry-After` headers automatically. The `backoff` function above is exponential.

---

## 4. Cache the lookup, never the FASTA

```
data/processed/uniprot_proteome_lookup.tsv   # append-only API cache
data/proteomes/<taxid>.fasta.gz              # downloaded FASTA, one per organism
```

- The lookup TSV stores **one row per (query_taxid, returned_proteome)**, plus a status row when the API returned no match or errored. Re-runs read this first and only query the API for new taxids.
- FASTA files are themselves the cache — there is no separate index.

Cache schema for lookups:

```
query_taxid     character    the taxid we asked about
upid            character    NA on no_match/error
proteome_type   character    NA on no_match/error
organism_taxid  character    actual taxid UniProt returned (often a strain)
organism_name   character
protein_count   integer
status          character    "ok" | "no_match" | "error"
error_msg       character    NA on success
fetched_at      character    ISO 8601 UTC
```

---

## 5. Pitfalls already burned

- **`organism_id:` instead of `taxonomy_id:`** — silent loss of Reference proteomes indexed under strain children. *Symptom*: zero hits for `Reference proteome` type across the whole run.
- **`fields=proteome_type`** — request is rejected with `Invalid fields parameter value 'proteome_type'`. The field exists in responses but cannot be requested by name.
- **Decompressing the FASTA download** — `compressed=true` makes the *body* a gzip file. If your client decompresses on the way in (some HTTP libraries do), you get plaintext FASTA written to a `.fasta.gz` file. `httr2` does the right thing only when you pass `path = ...` (streams bytes verbatim).
- **Per-second rate limit** — without `req_retry()`, occasional 429s halt the run. With it, they're invisible.
- **Stale cache after query change** — if you change query semantics (`organism_id:` → `taxonomy_id:`), delete the cache before re-running. Cached entries from the wrong query are not re-checked.

---

## Reference — minimal end-to-end snippet

```r
library(httr2)

look_up <- function(taxid) {
  request("https://rest.uniprot.org/proteomes/search") |>
    req_url_query(query  = paste0("taxonomy_id:", taxid),
                  format = "json",
                  fields = "upid,organism,organism_id,protein_count",
                  size   = 25L) |>
    req_user_agent("project/0.1") |>
    req_retry(max_tries = 4) |>
    req_perform() |>
    resp_body_json()
}

download <- function(upid, dest_path) {
  request("https://rest.uniprot.org/uniprotkb/stream") |>
    req_url_query(query = paste0("proteome:", upid),
                  format = "fasta",
                  compressed = "true") |>
    req_user_agent("project/0.1") |>
    req_retry(max_tries = 3) |>
    req_timeout(300) |>
    req_perform(path = dest_path)
}
```

Two functions. Everything else in `R/uniprot.R` is bookkeeping (caching, picking the best proteome, recording attrition).
