#!/usr/bin/env python3
"""Build a two-generation exhaustive and third-generation targeted backward-citation graph for CalorType set 2 using the public NIH iCite API."""
from __future__ import annotations

import csv
import gzip
import json
import math
import random
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator

BASE_URL = "https://icite.od.nih.gov/api/pubs"
BATCH_SIZE = 200
OUTPUT_DIR = Path("output")
SEED_FILE = Path("tmp_calortype/seed_set2.json")
USER_AGENT = "CalorType-crossref/0.1 (research use; contact: haikg96@gmail.com)"
MINIMAL_FIELDS = ["pmid", "year", "title", "authors", "journal", "doi", "is_research_article", "relative_citation_ratio", "nih_percentile", "human", "animal", "molecular_cellular", "is_clinical", "citation_count", "citations_per_year", "last_modified"]
STRONG_PHRASES = {
    "temperature-sensitive": 10, "temperature sensitive": 10, "thermosensitive": 10,
    "thermo-sensitive": 10, "heat-sensitive": 9, "heat sensitive": 9,
    "cold-sensitive": 9, "cold sensitive": 9, "conditional lethal": 9,
    "nonpermissive temperature": 9, "non-permissive temperature": 9,
    "permissive temperature": 8, "temperature-dependent mutant": 9,
    "temperature dependent mutant": 9, "heat-intolerant": 8,
    "heat intolerant": 8, "cold-intolerant": 8, "cold intolerant": 8,
    "febrile seizure": 8, "fever-sensitive": 8, "fever sensitive": 8,
}
TEMP_TERMS = ["temperature", "thermal", "thermo", "heat", "cold", "fever", "febrile", "hypertherm", "hypotherm", "nonpermissive", "permissive"]
VARIANT_TERMS = ["mutation", "mutant", "variant", "allele", "missense", "substitution", "amino acid", "genotype", "phenotype", "channelopathy", "protein stability", "destabil", "folding", "conditional", "loss of function", "gain of function"]
MECHANISM_TERMS = ["thermolabile", "thermostability", "heat shock", "protein folding", "proteostasis", "stability", "misfold", "denatur", "aggregation", "trafficking", "degradation", "conformation"]

def chunks(items: list[str], n: int) -> Iterator[list[str]]:
    for i in range(0, len(items), n):
        yield items[i:i+n]

def normalize_pmid(value: Any) -> str | None:
    if value is None:
        return None
    s = str(value).strip()
    if s.endswith(".0"):
        s = s[:-2]
    return s if s.isdigit() else None

def extract_records(payload: Any) -> list[dict[str, Any]]:
    if isinstance(payload, list):
        return [x for x in payload if isinstance(x, dict)]
    if isinstance(payload, dict):
        for key in ("data", "results", "publications"):
            value = payload.get(key)
            if isinstance(value, list):
                return [x for x in value if isinstance(x, dict)]
        if "pmid" in payload:
            return [payload]
    return []

def request_json(pmids: list[str], fields: list[str] | None, include_refs: bool) -> list[dict[str, Any]]:
    params = {"pmids": ",".join(pmids)}
    if fields:
        params["fl"] = ",".join(fields)
    if include_refs:
        params["refs"] = "true"
    url = BASE_URL + "?" + urllib.parse.urlencode(params)
    last_error: Exception | None = None
    for attempt in range(8):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as resp:
                payload = json.loads(resp.read().decode("utf-8"))
            return extract_records(payload)
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last_error = exc
            status = getattr(exc, "code", None)
            if status not in (None, 429, 500, 502, 503, 504):
                raise
            sleep_s = min(90.0, (2 ** attempt) + random.random())
            print(f"Retry {attempt+1}/8 after {type(exc).__name__}: {exc}; sleeping {sleep_s:.1f}s", flush=True)
            time.sleep(sleep_s)
    raise RuntimeError(f"iCite request failed after retries: {last_error}")

def fetch_records(pmids: Iterable[str], *, fields: list[str] | None, include_refs: bool, label: str) -> dict[str, dict[str, Any]]:
    ids = sorted({p for x in pmids if (p := normalize_pmid(x))})
    result: dict[str, dict[str, Any]] = {}
    batches = list(chunks(ids, BATCH_SIZE))
    for index, batch in enumerate(batches, start=1):
        records = request_json(batch, fields=fields, include_refs=include_refs)
        for record in records:
            pmid = normalize_pmid(record.get("pmid"))
            if pmid:
                record["pmid"] = pmid
                refs = record.get("references") or record.get("citedPmids") or []
                cited_by = record.get("cited_by") or record.get("citedByPmids") or []
                record["references"] = [p for x in refs if (p := normalize_pmid(x))]
                record["cited_by"] = [p for x in cited_by if (p := normalize_pmid(x))]
                result[pmid] = record
        if index == 1 or index % 10 == 0 or index == len(batches):
            print(f"{label}: batch {index}/{len(batches)}; records={len(result)}/{len(ids)}", flush=True)
        time.sleep(0.06)
    return result

def author_text(authors: Any) -> str:
    if not authors:
        return ""
    names: list[str] = []
    if isinstance(authors, list):
        for a in authors:
            if isinstance(a, str):
                names.append(a)
            elif isinstance(a, dict):
                full = a.get("fullName") or a.get("full_name")
                if not full:
                    full = " ".join(x for x in [str(a.get("firstName") or "").strip(), str(a.get("lastName") or "").strip()] if x)
                names.append(str(full).strip())
    elif isinstance(authors, str):
        return authors
    return "; ".join(x for x in names if x)

def record_title(record: dict[str, Any]) -> str:
    return str(record.get("title") or "").strip()

def relevance(record: dict[str, Any], genes: list[str]) -> tuple[int, list[str]]:
    title = record_title(record)
    text = re.sub(r"\s+", " ", title.lower())
    score = 0
    reasons: list[str] = []
    for phrase, points in STRONG_PHRASES.items():
        if phrase in text:
            score += points
            reasons.append(phrase)
            break
    has_temp = any(term in text for term in TEMP_TERMS)
    has_variant = any(term in text for term in VARIANT_TERMS)
    has_mechanism = any(term in text for term in MECHANISM_TERMS)
    if has_temp and has_variant:
        score += 6
        reasons.append("temperature+variant terminology")
    elif has_temp:
        score += 3
        reasons.append("temperature terminology")
    if has_mechanism and has_variant:
        score += 3
        reasons.append("stability/folding+variant terminology")
    elif has_mechanism:
        score += 1
        reasons.append("stability/folding terminology")
    upper_title = title.upper()
    matched_genes: list[str] = []
    for gene in genes:
        g = gene.strip()
        if len(g) < 3 or not re.fullmatch(r"[A-Za-z0-9_-]+", g):
            continue
        if re.search(rf"(?<![A-Z0-9]){re.escape(g.upper())}(?![A-Z0-9])", upper_title):
            matched_genes.append(g)
            if len(matched_genes) >= 5:
                break
    if matched_genes:
        score += 5
        reasons.append("gene:" + ",".join(matched_genes))
    if record.get("is_research_article") is True:
        score += 1
        reasons.append("primary research")
    try:
        c = float(record.get("citation_count") or 0)
        if c >= 100:
            score += 2
            reasons.append("highly cited")
        elif c >= 20:
            score += 1
            reasons.append("cited")
    except (TypeError, ValueError):
        pass
    return score, reasons

def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding="utf-8")

def write_jsonl_gz(path: Path, records: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        for record in records:
            fh.write(json.dumps(record, ensure_ascii=False, separators=(",", ":")) + "\n")

def write_csv_gz(path: Path, rows: Iterable[dict[str, Any]], fieldnames: list[str]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            count += 1
    return count

def metadata_row(record: dict[str, Any]) -> dict[str, Any]:
    pmid = normalize_pmid(record.get("pmid")) or ""
    doi = str(record.get("doi") or "").strip()
    return {
        "pmid": pmid,
        "title": record_title(record),
        "year": record.get("year") or record.get("pubYear") or "",
        "journal": record.get("journal") or record.get("journalNameIso") or "",
        "authors": author_text(record.get("authors")),
        "doi": doi,
        "is_research_article": record.get("is_research_article", record.get("iCiteArticle", "")),
        "citation_count": record.get("citation_count", record.get("citedByPmidCount", "")),
        "relative_citation_ratio": record.get("relative_citation_ratio", record.get("rcr", "")),
        "nih_percentile": record.get("nih_percentile", record.get("nihRcrPercentile", "")),
        "human": record.get("human", ""), "animal": record.get("animal", ""),
        "molecular_cellular": record.get("molecular_cellular", record.get("molCell", "")),
        "is_clinical": record.get("is_clinical", record.get("isClinicalArticle", "")),
        "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "",
        "doi_url": f"https://doi.org/{doi}" if doi else "",
    }

def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    seed_payload = json.loads(SEED_FILE.read_text(encoding="utf-8"))
    seed_entries = seed_payload.get("seeds") or []
    raw_seed_ids = seed_payload.get("pmids") or [x.get("pmid") for x in seed_entries]
    seed_ids = [normalize_pmid(x) for x in raw_seed_ids]
    seed_ids = [x for x in seed_ids if x]
    genes = [str(x) for x in seed_payload.get("genes", []) if x]
    seed_input = {str(x["pmid"]): x for x in seed_entries if x.get("pmid")}
    print(f"Starting CalorType citation expansion from {len(seed_ids)} seed PMIDs", flush=True)
    seeds = fetch_records(seed_ids, fields=None, include_refs=True, label="seeds")
    missing_seeds = sorted(set(seed_ids) - set(seeds))
    print(f"Seed coverage: {len(seeds)}/{len(seed_ids)}; missing={len(missing_seeds)}", flush=True)

    edges1: list[dict[str, Any]] = []
    origin1: dict[str, set[str]] = defaultdict(set)
    for seed_id, record in seeds.items():
        for ref in record.get("references", []):
            if ref != seed_id:
                edges1.append({"source_pmid": seed_id, "target_pmid": ref, "generation": 1})
                origin1[ref].add(seed_id)
    gen1_ids = sorted(set(origin1) - set(seed_ids))
    print(f"Generation 1: edges={len(edges1):,}; unique new PMIDs={len(gen1_ids):,}", flush=True)
    gen1 = fetch_records(gen1_ids, fields=None, include_refs=True, label="generation 1 full")
    unresolved1 = sorted(set(gen1_ids) - set(gen1))

    edges2: list[dict[str, Any]] = []
    parents2: dict[str, set[str]] = defaultdict(set)
    origin2: dict[str, set[str]] = defaultdict(set)
    for parent_id, record in gen1.items():
        inherited = origin1.get(parent_id, set())
        for ref in record.get("references", []):
            if ref != parent_id:
                edges2.append({"source_pmid": parent_id, "target_pmid": ref, "generation": 2})
                parents2[ref].add(parent_id)
                if inherited and len(origin2[ref]) < 50:
                    origin2[ref].update(list(inherited)[:50 - len(origin2[ref])])
    seen01 = set(seed_ids) | set(gen1_ids)
    gen2_ids = sorted(set(parents2) - seen01)
    print(f"Generation 2: edges={len(edges2):,}; unique new PMIDs={len(gen2_ids):,}", flush=True)
    gen2 = fetch_records(gen2_ids, fields=MINIMAL_FIELDS, include_refs=False, label="generation 2 metadata")
    unresolved2 = sorted(set(gen2_ids) - set(gen2))

    scores1 = {pmid: relevance(record, genes) for pmid, record in gen1.items()}
    scores2 = {pmid: relevance(record, genes) for pmid, record in gen2.items()}
    frontier: set[str] = {pmid for pmid, (score, _) in scores2.items() if score >= 4}
    for gen1_parent, parent_record in gen1.items():
        if scores1.get(gen1_parent, (0, []))[0] < 6:
            continue
        refs = [r for r in parent_record.get("references", []) if r in gen2]
        refs.sort(key=lambda r: (scores2.get(r, (0, []))[0], float(gen2.get(r, {}).get("citation_count") or 0)), reverse=True)
        frontier.update(refs[:5])
    def frontier_key(pmid: str) -> tuple[float, float, float]:
        rec = gen2.get(pmid, {})
        return (float(scores2.get(pmid, (0, []))[0]), math.log1p(float(rec.get("citation_count") or 0)), float(rec.get("year") or 0))
    frontier_ids = sorted(frontier, key=frontier_key, reverse=True)[:5000]
    print(f"Generation 3 frontier: {len(frontier_ids):,} generation-2 parents selected", flush=True)
    frontier_full = fetch_records(frontier_ids, fields=None, include_refs=True, label="generation 2 frontier full")

    edges3: list[dict[str, Any]] = []
    parents3: dict[str, set[str]] = defaultdict(set)
    origin3: dict[str, set[str]] = defaultdict(set)
    for parent_id, record in frontier_full.items():
        inherited = origin2.get(parent_id, set())
        for ref in record.get("references", []):
            if ref != parent_id:
                edges3.append({"source_pmid": parent_id, "target_pmid": ref, "generation": 3})
                parents3[ref].add(parent_id)
                if inherited and len(origin3[ref]) < 50:
                    origin3[ref].update(list(inherited)[:50 - len(origin3[ref])])
    seen012 = seen01 | set(gen2_ids)
    gen3_ids = sorted(set(parents3) - seen012)
    print(f"Generation 3 targeted: edges={len(edges3):,}; unique new PMIDs={len(gen3_ids):,}", flush=True)
    gen3 = fetch_records(gen3_ids, fields=MINIMAL_FIELDS, include_refs=False, label="generation 3 metadata")
    unresolved3 = sorted(set(gen3_ids) - set(gen3))
    scores3 = {pmid: relevance(record, genes) for pmid, record in gen3.items()}

    all_meta: dict[str, dict[str, Any]] = {}
    all_meta.update(seeds); all_meta.update(gen1); all_meta.update(gen2); all_meta.update(gen3)
    def edge_rows(edges: list[dict[str, Any]]) -> Iterator[dict[str, Any]]:
        for edge in edges:
            s, t = edge["source_pmid"], edge["target_pmid"]
            yield {**edge, "source_title": record_title(all_meta.get(s, {})), "target_title": record_title(all_meta.get(t, {})), "source_pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{s}/", "target_pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{t}/"}
    edge_fields = ["generation", "source_pmid", "source_title", "target_pmid", "target_title", "source_pubmed_url", "target_pubmed_url"]
    n_edge1 = write_csv_gz(OUTPUT_DIR / "edges_generation_1.csv.gz", edge_rows(edges1), edge_fields)
    n_edge2 = write_csv_gz(OUTPUT_DIR / "edges_generation_2.csv.gz", edge_rows(edges2), edge_fields)
    n_edge3 = write_csv_gz(OUTPUT_DIR / "edges_generation_3_targeted.csv.gz", edge_rows(edges3), edge_fields)
    write_jsonl_gz(OUTPUT_DIR / "seed_records.jsonl.gz", seeds.values())
    write_jsonl_gz(OUTPUT_DIR / "works_generation_1.jsonl.gz", gen1.values())
    write_jsonl_gz(OUTPUT_DIR / "works_generation_2.jsonl.gz", gen2.values())
    write_jsonl_gz(OUTPUT_DIR / "works_generation_3_targeted.jsonl.gz", gen3.values())

    queue: list[dict[str, Any]] = []
    generation_data = [(1, gen1, scores1, origin1, defaultdict(set)), (2, gen2, scores2, origin2, parents2), (3, gen3, scores3, origin3, parents3)]
    for generation, recs, scores, origins, parents in generation_data:
        for pmid, rec in recs.items():
            score, reasons = scores.get(pmid, (0, []))
            origin_ids = sorted(origins.get(pmid, set()))
            parent_ids = sorted(parents.get(pmid, set()))
            row = metadata_row(rec)
            row.update({"generation": generation, "relevance_score": score, "relevance_reasons": "; ".join(reasons), "seed_origin_count": len(origin_ids), "seed_pmids_sample": ";".join(origin_ids[:10]), "seed_titles_sample": " | ".join(str((seed_input.get(x) or seeds.get(x) or {}).get("title") or "") for x in origin_ids[:5]), "direct_parent_count": len(parent_ids), "direct_parent_pmids_sample": ";".join(parent_ids[:10])})
            queue.append(row)
    queue.sort(key=lambda r: (int(r["relevance_score"] or 0), float(r["citation_count"] or 0), int(r["year"] or 0)), reverse=True)
    queue_fields = ["generation", "relevance_score", "relevance_reasons", "pmid", "title", "year", "journal", "authors", "doi", "is_research_article", "citation_count", "relative_citation_ratio", "nih_percentile", "human", "animal", "molecular_cellular", "is_clinical", "seed_origin_count", "seed_pmids_sample", "seed_titles_sample", "direct_parent_count", "direct_parent_pmids_sample", "pubmed_url", "doi_url"]
    queue_count = write_csv_gz(OUTPUT_DIR / "review_queue.csv.gz", queue, queue_fields)

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(), "project": seed_payload.get("project"), "set": seed_payload.get("set"),
        "method": {"direction": "backward citations / cited references", "source": "NIH iCite Open Citation Collection API", "generation_1": "exhaustive for PMID-resolved references from seed papers", "generation_2": "exhaustive for PMID-resolved references from generation-1 papers", "generation_3": "targeted: expanded generation-2 papers with relevance score >=4, plus up to five influential references per highly relevant generation-1 parent; frontier capped at 5,000 parents", "deduplication": "PMID; each work assigned to its earliest newly discovered generation", "limitations": ["References without PubMed identifiers are outside iCite and are not represented.", "Generation 3 is targeted rather than exhaustive to control citation-network explosion.", "Relevance scoring is title/metadata based and is intended for triage, not final inclusion."]},
        "counts": {"seed_input": len(seed_ids), "seed_resolved": len(seeds), "seed_missing": len(missing_seeds), "generation_1_edges": n_edge1, "generation_1_unique_new": len(gen1_ids), "generation_1_resolved": len(gen1), "generation_1_unresolved": len(unresolved1), "generation_2_edges": n_edge2, "generation_2_unique_new": len(gen2_ids), "generation_2_resolved": len(gen2), "generation_2_unresolved": len(unresolved2), "generation_3_frontier_parents": len(frontier_ids), "generation_3_edges": n_edge3, "generation_3_unique_new": len(gen3_ids), "generation_3_resolved": len(gen3), "generation_3_unresolved": len(unresolved3), "review_queue_rows": queue_count},
        "missing_pmids": {"seeds": missing_seeds, "generation_1": unresolved1, "generation_2": unresolved2[:10000], "generation_3": unresolved3[:10000], "generation_2_truncated": len(unresolved2) > 10000, "generation_3_truncated": len(unresolved3) > 10000},
        "relevance_score_distribution": {"generation_1": dict(sorted({str(k): v for k, v in Counter(score for score, _ in scores1.values()).items()}.items(), key=lambda kv: int(kv[0]), reverse=True)), "generation_2": dict(sorted({str(k): v for k, v in Counter(score for score, _ in scores2.values()).items()}.items(), key=lambda kv: int(kv[0]), reverse=True)), "generation_3": dict(sorted({str(k): v for k, v in Counter(score for score, _ in scores3.values()).items()}.items(), key=lambda kv: int(kv[0]), reverse=True))}
    }
    write_json(OUTPUT_DIR / "summary.json", summary)
    write_json(OUTPUT_DIR / "frontier_generation_2_pmids.json", frontier_ids)
    print(json.dumps(summary["counts"], indent=2), flush=True)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
