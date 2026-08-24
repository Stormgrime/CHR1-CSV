#!/usr/bin/env python3
"""Build a multi-generation backward-citation graph for CalorType set 2.

The script uses NIH iCite/Open Citation Collection for PubMed-to-PubMed
reference links and NCBI E-utilities for abstracts. It expands all seeds two
complete generations, then expands a relevance-ranked subset of generation 2
for a selective third generation.
"""

from __future__ import annotations

import csv
import gzip
import html
import io
import json
import math
import os
import re
import shutil
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, Sequence, Set, Tuple

ROOT = Path(__file__).resolve().parents[1]
INPUT_TSV = ROOT / "data" / "calortype_set2_seed_map.tsv"
OUT_DIR = ROOT / "output" / "calortype_set2_crossrefs"
BUNDLE = ROOT / "output" / "calortype_set2_crossrefs_bundle.zip"

ICITE_ENDPOINT = "https://icite.od.nih.gov/api/pubs"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
USER_AGENT = "CalorTypeCitationSnowball/1.0"

ICITE_BATCH = 100
EFETCH_BATCH = 180
MAX_ABSTRACT_FETCH = 60_000
MAX_G3_FRONTIER = 3_000
MAX_G3_ABSTRACT_FETCH = 50_000
REQUEST_TIMEOUT = 120

TEMP_STRONG_PATTERNS = [
    r"temperature[- ]sensitive",
    r"temperature sensitivity",
    r"temperature[- ]dependent",
    r"thermo[- ]?sensitive",
    r"thermolabile",
    r"heat[- ]sensitive",
    r"heat[- ]labile",
    r"cold[- ]sensitive",
    r"thermal instability",
    r"non[- ]?permissive",
    r"permissive temperature",
    r"restrictive temperature",
    r"febrile",
    r"hyperthermi\w*",
    r"malignant hyperthermia",
    r"heat intolerance",
    r"cold[- ]induced",
    r"heat[- ]induced",
    r"thermosensitive",
]
TEMP_GENERAL_PATTERNS = [
    r"\btemperature\w*\b",
    r"\bthermal\w*\b",
    r"\bheat\b",
    r"\bheated\b",
    r"\bcold\b",
    r"\bfever\w*\b",
    r"\bfebrile\b",
    r"\bhypertherm\w*\b",
    r"\bhypotherm\w*\b",
    r"\bthermo\w*\b",
    r"heat shock",
    r"heat stress",
]
VARIANT_PATTERNS = [
    r"\bmutant\w*\b",
    r"\bmutation\w*\b",
    r"\bvariant\w*\b",
    r"\ballele\w*\b",
    r"\bpolymorphism\w*\b",
    r"\bsubstitution\w*\b",
    r"\bmissense\b",
    r"\bnonsense\b",
    r"\bdeletion\w*\b",
    r"\bconditional\b",
    r"\bgenotyp\w*\b",
    r"\bdefect\w*\b",
    r"loss[- ]of[- ]function",
    r"gain[- ]of[- ]function",
]
MECHANISM_PATTERNS = [
    r"\bstabilit\w*\b",
    r"\bunstable\b",
    r"\bfold\w*\b",
    r"\bmisfold\w*\b",
    r"\bdegrad\w*\b",
    r"\btraffick\w*\b",
    r"\bconformation\w*\b",
    r"\bactivity\b",
    r"\bfunction\w*\b",
    r"\bchannel\w*\b",
    r"\benzyme\w*\b",
    r"\bprotein\w*\b",
    r"\bexpression\b",
]
REVIEW_PATTERNS = [r"\breview\b", r"meta-analysis", r"systematic review"]

STRONG_RX = [re.compile(p, re.I) for p in TEMP_STRONG_PATTERNS]
TEMP_RX = [re.compile(p, re.I) for p in TEMP_GENERAL_PATTERNS]
VARIANT_RX = [re.compile(p, re.I) for p in VARIANT_PATTERNS]
MECH_RX = [re.compile(p, re.I) for p in MECHANISM_PATTERNS]
REVIEW_RX = [re.compile(p, re.I) for p in REVIEW_PATTERNS]

GENE_STOPWORDS = {
    "alpha", "beta", "gamma", "protein", "gene", "l", "ha", "3c",
    "cct", "fas", "ns5", "aca", "tpi", "src", "mag", "dcase",
}


def log(message: str) -> None:
    print(f"[{datetime.now(timezone.utc).isoformat(timespec='seconds')}] {message}", flush=True)


def chunks(items: Sequence[str], size: int) -> Iterator[List[str]]:
    for start in range(0, len(items), size):
        yield list(items[start:start + size])


def pmid_string(value: object) -> str:
    try:
        return str(int(value))
    except (TypeError, ValueError):
        return ""


def clean_text(value: object) -> str:
    if value is None:
        return ""
    text = html.unescape(str(value))
    text = re.sub(r"<[^>]+>", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def request_bytes(request: urllib.request.Request, retries: int = 7) -> bytes:
    delay = 1.0
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=REQUEST_TIMEOUT) as response:
                return response.read()
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            last_error = exc
            retryable = not isinstance(exc, urllib.error.HTTPError) or exc.code in {408, 429, 500, 502, 503, 504}
            if not retryable or attempt == retries:
                break
            log(f"Request failed ({exc}); retry {attempt}/{retries} after {delay:.1f}s")
            time.sleep(delay)
            delay = min(delay * 1.8, 30.0)
    raise RuntimeError(f"Request failed after {retries} attempts: {last_error}")


def fetch_icite(pmids: Iterable[str], include_refs: bool, phase: str) -> Dict[str, dict]:
    ids = sorted({pmid_string(p) for p in pmids if pmid_string(p)}, key=int)
    result: Dict[str, dict] = {}
    batches = list(chunks(ids, ICITE_BATCH))
    log(f"iCite {phase}: {len(ids):,} PMIDs in {len(batches):,} batches; refs={include_refs}")

    def fetch_one(batch: List[str]) -> None:
        params = {"pmids": ",".join(batch)}
        if include_refs:
            params["refs"] = "true"
        url = ICITE_ENDPOINT + "?" + urllib.parse.urlencode(params, safe=",")
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
        try:
            payload = json.loads(request_bytes(request).decode("utf-8"))
        except Exception:
            if len(batch) == 1:
                log(f"Skipping unresolvable iCite PMID {batch[0]}")
                return
            middle = len(batch) // 2
            fetch_one(batch[:middle])
            fetch_one(batch[middle:])
            return
        for record in payload.get("data", []):
            pmid = pmid_string(record.get("pmid"))
            if pmid:
                result[pmid] = record

    for index, batch in enumerate(batches, 1):
        fetch_one(batch)
        if index == 1 or index % 10 == 0 or index == len(batches):
            log(f"iCite {phase}: batch {index:,}/{len(batches):,}; records={len(result):,}")
        time.sleep(0.08)
    return result


def element_text(element: ET.Element | None) -> str:
    if element is None:
        return ""
    return clean_text("".join(element.itertext()))


def fetch_pubmed_abstracts(pmids: Iterable[str], phase: str) -> Dict[str, dict]:
    ids = sorted({pmid_string(p) for p in pmids if pmid_string(p)}, key=int)
    result: Dict[str, dict] = {}
    batches = list(chunks(ids, EFETCH_BATCH))
    log(f"PubMed {phase}: {len(ids):,} PMIDs in {len(batches):,} batches")
    for index, batch in enumerate(batches, 1):
        encoded = urllib.parse.urlencode({
            "db": "pubmed",
            "id": ",".join(batch),
            "retmode": "xml",
        }).encode("utf-8")
        request = urllib.request.Request(
            NCBI_EFETCH,
            data=encoded,
            headers={
                "User-Agent": USER_AGENT,
                "Content-Type": "application/x-www-form-urlencoded",
            },
        )
        try:
            xml_bytes = request_bytes(request)
            root = ET.fromstring(xml_bytes)
        except Exception as exc:
            log(f"PubMed batch {index} failed and was skipped: {exc}")
            time.sleep(0.5)
            continue
        for article in root.findall(".//PubmedArticle"):
            pmid = element_text(article.find(".//MedlineCitation/PMID"))
            if not pmid:
                continue
            title = element_text(article.find(".//Article/ArticleTitle"))
            abstract_parts = []
            for node in article.findall(".//Article/Abstract/AbstractText"):
                label = clean_text(node.attrib.get("Label"))
                text = element_text(node)
                if not text:
                    continue
                abstract_parts.append(f"{label}: {text}" if label else text)
            abstract = " ".join(abstract_parts)
            doi = ""
            for eid in article.findall(".//PubmedData/ArticleIdList/ArticleId"):
                if eid.attrib.get("IdType") == "doi":
                    doi = element_text(eid)
                    break
            result[pmid] = {"title": title, "abstract": abstract, "doi": doi}
        if index == 1 or index % 10 == 0 or index == len(batches):
            log(f"PubMed {phase}: batch {index:,}/{len(batches):,}; records={len(result):,}")
        time.sleep(0.37)
    return result


def load_seeds() -> Tuple[List[str], Dict[str, List[str]]]:
    seeds: List[str] = []
    seed_genes: Dict[str, List[str]] = {}
    with INPUT_TSV.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            pmid = pmid_string(row.get("pmid"))
            if not pmid:
                continue
            genes = [g.strip() for g in (row.get("genes") or "").split("|") if g.strip()]
            seeds.append(pmid)
            seed_genes[pmid] = genes
    if not seeds:
        raise RuntimeError(f"No seeds found in {INPUT_TSV}")
    if len(seeds) != len(set(seeds)):
        raise RuntimeError("Seed PMID list contains duplicates")
    return seeds, seed_genes


def refs_for(record: Mapping[str, object] | None) -> List[str]:
    if not record:
        return []
    refs = record.get("references") or record.get("citedPmids") or []
    output = []
    seen = set()
    for ref in refs if isinstance(refs, list) else []:
        pmid = pmid_string(ref)
        if pmid and pmid not in seen:
            output.append(pmid)
            seen.add(pmid)
    return output


def add_reach(
    reach: MutableMapping[str, Dict[str, str]],
    child: str,
    seed: str,
    parent: str,
) -> None:
    if child not in reach:
        reach[child] = {}
    reach[child].setdefault(seed, parent)


def normalized_gene_aliases(raw_gene: str) -> List[str]:
    candidates = re.split(r"[;|/]", raw_gene)
    aliases = []
    for candidate in candidates:
        gene = candidate.strip()
        if not gene:
            continue
        gene = re.sub(r"\([^)]*\)", "", gene).strip()
        if not gene:
            continue
        if len(gene) < 3 or gene.lower() in GENE_STOPWORDS:
            continue
        aliases.append(gene)
    return aliases


def associated_genes(seed_ids: Iterable[str], seed_genes: Mapping[str, List[str]]) -> List[str]:
    genes: List[str] = []
    seen: Set[str] = set()
    for seed in seed_ids:
        for raw in seed_genes.get(seed, []):
            for gene in normalized_gene_aliases(raw):
                key = gene.casefold()
                if key not in seen:
                    genes.append(gene)
                    seen.add(key)
    return genes


def match_patterns(text: str, regexes: Sequence[re.Pattern[str]], labels: Sequence[str]) -> List[str]:
    matches = []
    for regex, label in zip(regexes, labels):
        if regex.search(text):
            matches.append(label)
    return matches


def gene_matches(text: str, genes: Iterable[str]) -> List[str]:
    matches = []
    for gene in genes:
        escaped = re.escape(gene)
        if re.search(rf"(?<![A-Za-z0-9]){escaped}(?![A-Za-z0-9])", text, flags=re.I):
            matches.append(gene)
    return matches


def score_node(
    pmid: str,
    metadata: Mapping[str, Mapping[str, object]],
    abstracts: Mapping[str, Mapping[str, str]],
    node_seed_ids: Iterable[str],
    seed_genes: Mapping[str, List[str]],
    min_generation: int,
) -> dict:
    meta = metadata.get(pmid, {})
    extra = abstracts.get(pmid, {})
    title = clean_text(extra.get("title") or meta.get("title"))
    abstract = clean_text(extra.get("abstract"))
    full = f"{title} {abstract}".strip()

    title_strong = match_patterns(title, STRONG_RX, TEMP_STRONG_PATTERNS)
    full_strong = match_patterns(full, STRONG_RX, TEMP_STRONG_PATTERNS)
    title_temp = match_patterns(title, TEMP_RX, TEMP_GENERAL_PATTERNS)
    full_temp = match_patterns(full, TEMP_RX, TEMP_GENERAL_PATTERNS)
    title_variant = match_patterns(title, VARIANT_RX, VARIANT_PATTERNS)
    full_variant = match_patterns(full, VARIANT_RX, VARIANT_PATTERNS)
    full_mech = match_patterns(full, MECH_RX, MECHANISM_PATTERNS)
    reviews = match_patterns(full, REVIEW_RX, REVIEW_PATTERNS)

    seed_list = sorted(set(node_seed_ids), key=int)
    genes = associated_genes(seed_list, seed_genes)
    matched_genes = gene_matches(full, genes)
    matched_genes_title = gene_matches(title, genes)

    score = 0.0
    if title_strong:
        score += 7.0
    elif full_strong:
        score += 5.0
    if title_temp:
        score += 3.0
    elif full_temp:
        score += 2.0
    if title_variant:
        score += 2.5
    elif full_variant:
        score += 1.5
    if full_mech:
        score += 1.0
    if matched_genes:
        score += 4.0
    if matched_genes_title:
        score += 1.0
    if (full_strong or full_temp) and full_variant:
        score += 2.0
    if (full_strong or full_temp) and matched_genes:
        score += 2.0
    seed_count = len(seed_list)
    if seed_count >= 2:
        score += min(3.0, math.log2(seed_count))
    if min_generation == 1:
        score += 1.0
    if reviews and not (full_strong or matched_genes):
        score -= 0.5
    score = round(max(score, 0.0), 1)

    if score >= 11:
        priority = "A"
    elif score >= 8:
        priority = "B"
    elif score >= 5.5:
        priority = "C"
    elif seed_count >= 3 and score >= 3:
        priority = "Context"
    else:
        priority = "Low"

    term_matches = []
    for value in title_strong + full_strong + title_temp + full_temp + title_variant + full_variant + full_mech:
        if value not in term_matches:
            term_matches.append(value)
    return {
        "score": score,
        "priority": priority,
        "matched_genes": matched_genes,
        "matched_terms": term_matches,
        "title": title,
        "abstract": abstract,
    }


def broad_abstract_prefilter(
    pmid: str,
    metadata: Mapping[str, Mapping[str, object]],
    seed_ids: Iterable[str],
    seed_genes: Mapping[str, List[str]],
    seed_count: int,
) -> bool:
    title = clean_text(metadata.get(pmid, {}).get("title"))
    if not title:
        return seed_count >= 2
    if any(rx.search(title) for rx in STRONG_RX + TEMP_RX + VARIANT_RX + MECH_RX):
        return True
    if gene_matches(title, associated_genes(seed_ids, seed_genes)):
        return True
    return seed_count >= 2


def merge_metadata(*sources: Mapping[str, Mapping[str, object]]) -> Dict[str, dict]:
    merged: Dict[str, dict] = {}
    for source in sources:
        for pmid, record in source.items():
            if pmid not in merged:
                merged[pmid] = dict(record)
            else:
                for key, value in record.items():
                    if value not in (None, "", [], {}) or key not in merged[pmid]:
                        merged[pmid][key] = value
    return merged


def write_csv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]], gzip_output: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if gzip_output else open
    kwargs = {"mode": "wt", "encoding": "utf-8", "newline": ""}
    with opener(path, **kwargs) as handle:  # type: ignore[arg-type]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def representative_path(
    node: str,
    seed: str,
    g1_reach: Mapping[str, Mapping[str, str]],
    g2_reach: Mapping[str, Mapping[str, str]],
    g3_reach: Mapping[str, Mapping[str, str]],
) -> Tuple[int, List[str]]:
    if seed in g1_reach.get(node, {}):
        return 1, [seed, node]
    if seed in g2_reach.get(node, {}):
        parent = g2_reach[node][seed]
        return 2, [seed, parent, node]
    if seed in g3_reach.get(node, {}):
        parent2 = g3_reach[node][seed]
        parent1 = g2_reach.get(parent2, {}).get(seed)
        if parent1:
            return 3, [seed, parent1, parent2, node]
        return 3, [seed, parent2, node]
    return 99, [seed, node]


def main() -> int:
    start_time = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    seeds, seed_genes = load_seeds()
    seed_set = set(seeds)
    log(f"Loaded {len(seeds):,} unique seeds")

    seed_records = fetch_icite(seeds, include_refs=True, phase="seed generation")
    g1_edges: Set[Tuple[str, str]] = set()
    g1_reach: Dict[str, Dict[str, str]] = {}
    for seed in seeds:
        for child in refs_for(seed_records.get(seed)):
            if child == seed:
                continue
            g1_edges.add((seed, child))
            add_reach(g1_reach, child, seed, seed)
    g1_nodes = set(g1_reach)
    log(f"Generation 1: {len(g1_edges):,} unique edges; {len(g1_nodes):,} unique nodes")

    g1_records = fetch_icite(g1_nodes, include_refs=True, phase="generation 1 expansion")
    g2_edges: Set[Tuple[str, str]] = set()
    g2_reach: Dict[str, Dict[str, str]] = {}
    for parent in sorted(g1_nodes, key=int):
        parent_seeds = g1_reach.get(parent, {})
        for child in refs_for(g1_records.get(parent)):
            if child == parent:
                continue
            g2_edges.add((parent, child))
            for seed in parent_seeds:
                add_reach(g2_reach, child, seed, parent)
    g2_nodes = set(g2_reach)
    log(f"Generation 2: {len(g2_edges):,} unique edges; {len(g2_nodes):,} unique nodes")

    g2_records = fetch_icite(g2_nodes, include_refs=False, phase="generation 2 metadata")
    metadata = merge_metadata(seed_records, g1_records, g2_records)

    pre_node_seeds: Dict[str, Set[str]] = defaultdict(set)
    for node, support in g1_reach.items():
        pre_node_seeds[node].update(support)
    for node, support in g2_reach.items():
        pre_node_seeds[node].update(support)
    pre_min_generation = {
        node: 1 if node in g1_reach else 2
        for node in set(g1_reach) | set(g2_reach)
    }

    abstract_pool = set(g1_nodes) | set(g2_nodes)
    if len(abstract_pool) > MAX_ABSTRACT_FETCH:
        g2_prefiltered = {
            pmid for pmid in g2_nodes
            if broad_abstract_prefilter(
                pmid,
                metadata,
                pre_node_seeds.get(pmid, set()),
                seed_genes,
                len(pre_node_seeds.get(pmid, set())),
            )
        }
        ranked_remaining = sorted(
            g2_nodes - g2_prefiltered,
            key=lambda p: (-len(pre_node_seeds.get(p, set())), int(p)),
        )
        room = max(0, MAX_ABSTRACT_FETCH - len(g1_nodes) - len(g2_prefiltered))
        abstract_pool = set(g1_nodes) | g2_prefiltered | set(ranked_remaining[:room])
        log(f"Abstract pool capped at {len(abstract_pool):,} papers")
    abstracts = fetch_pubmed_abstracts(abstract_pool, phase="G1/G2 abstracts")

    pre_scores: Dict[str, dict] = {}
    for pmid in g2_nodes:
        pre_scores[pmid] = score_node(
            pmid,
            metadata,
            abstracts,
            pre_node_seeds.get(pmid, set()),
            seed_genes,
            pre_min_generation.get(pmid, 2),
        )

    frontier: Set[str] = {
        pmid for pmid, score in pre_scores.items()
        if score["score"] >= 7.0
        or (len(pre_node_seeds.get(pmid, set())) >= 3 and score["score"] >= 4.0)
    }
    for seed in seeds:
        per_seed = [
            pmid for pmid in g2_nodes
            if seed in g2_reach.get(pmid, {}) and pre_scores[pmid]["score"] >= 4.0
        ]
        per_seed.sort(key=lambda p: (-pre_scores[p]["score"], -len(pre_node_seeds.get(p, set())), int(p)))
        frontier.update(per_seed[:8])
    frontier.discard("")
    frontier.discard(None)  # type: ignore[arg-type]
    frontier_sorted = sorted(
        frontier,
        key=lambda p: (-pre_scores[p]["score"], -len(pre_node_seeds.get(p, set())), int(p)),
    )
    if len(frontier_sorted) > MAX_G3_FRONTIER:
        frontier_sorted = frontier_sorted[:MAX_G3_FRONTIER]
    frontier = set(frontier_sorted)
    log(f"Selective generation 3 frontier: {len(frontier):,} G2 papers")

    g2_frontier_records = fetch_icite(frontier, include_refs=True, phase="generation 3 frontier")
    metadata = merge_metadata(metadata, g2_frontier_records)
    g3_edges: Set[Tuple[str, str]] = set()
    g3_reach: Dict[str, Dict[str, str]] = {}
    for parent in frontier_sorted:
        parent_seeds = g2_reach.get(parent, {})
        for child in refs_for(g2_frontier_records.get(parent)):
            if child == parent:
                continue
            g3_edges.add((parent, child))
            for seed in parent_seeds:
                add_reach(g3_reach, child, seed, parent)
    g3_nodes = set(g3_reach)
    log(f"Generation 3 selective: {len(g3_edges):,} unique edges; {len(g3_nodes):,} unique nodes")

    g3_records = fetch_icite(g3_nodes, include_refs=False, phase="generation 3 metadata")
    metadata = merge_metadata(metadata, g3_records)

    node_seeds: Dict[str, Set[str]] = defaultdict(set)
    node_parents: Dict[str, Set[str]] = defaultdict(set)
    node_generations: Dict[str, Set[int]] = defaultdict(set)
    for generation, edges, reach in [
        (1, g1_edges, g1_reach),
        (2, g2_edges, g2_reach),
        (3, g3_edges, g3_reach),
    ]:
        for parent, child in edges:
            node_parents[child].add(parent)
            node_generations[child].add(generation)
        for node, support in reach.items():
            node_seeds[node].update(support)
    min_generation = {node: min(gens) for node, gens in node_generations.items()}

    g3_abstract_pool = set(g3_nodes)
    if len(g3_abstract_pool) > MAX_G3_ABSTRACT_FETCH:
        prefiltered = {
            pmid for pmid in g3_nodes
            if broad_abstract_prefilter(
                pmid,
                metadata,
                node_seeds.get(pmid, set()),
                seed_genes,
                len(node_seeds.get(pmid, set())),
            )
        }
        ranked = sorted(
            g3_nodes - prefiltered,
            key=lambda p: (-len(node_seeds.get(p, set())), int(p)),
        )
        room = max(0, MAX_G3_ABSTRACT_FETCH - len(prefiltered))
        g3_abstract_pool = prefiltered | set(ranked[:room])
        log(f"G3 abstract pool capped at {len(g3_abstract_pool):,} papers")
    g3_abstracts = fetch_pubmed_abstracts(g3_abstract_pool, phase="G3 abstracts")
    abstracts = {**abstracts, **g3_abstracts}

    all_crossref_nodes = set(g1_nodes) | set(g2_nodes) | set(g3_nodes)
    scores: Dict[str, dict] = {}
    for index, pmid in enumerate(sorted(all_crossref_nodes, key=int), 1):
        scores[pmid] = score_node(
            pmid,
            metadata,
            abstracts,
            node_seeds.get(pmid, set()),
            seed_genes,
            min_generation.get(pmid, 99),
        )
        if index % 25_000 == 0:
            log(f"Scored {index:,}/{len(all_crossref_nodes):,} cross-reference nodes")

    candidate_pmids = [
        pmid for pmid in all_crossref_nodes
        if pmid not in seed_set and scores[pmid]["priority"] in {"A", "B", "C", "Context"}
    ]
    candidate_pmids.sort(
        key=lambda p: (
            {"A": 0, "B": 1, "C": 2, "Context": 3}.get(scores[p]["priority"], 9),
            -scores[p]["score"],
            min_generation.get(p, 99),
            -len(node_seeds.get(p, set())),
            int(p),
        )
    )
    rank_by_pmid = {pmid: rank for rank, pmid in enumerate(candidate_pmids, 1)}
    log(f"Shortlisted {len(candidate_pmids):,} unique non-seed candidates")

    node_fields = [
        "pmid", "is_seed", "min_generation", "generations", "title", "year",
        "journal", "doi", "citation_count", "seed_count", "parent_count",
        "score", "priority", "matched_genes", "matched_terms", "has_abstract",
        "pubmed_url", "icite_url",
    ]
    node_rows = []
    all_nodes = seed_set | all_crossref_nodes
    for pmid in sorted(all_nodes, key=int):
        meta = metadata.get(pmid, {})
        score = scores.get(pmid, {"score": "", "priority": "Seed", "matched_genes": [], "matched_terms": []})
        extra = abstracts.get(pmid, {})
        doi = clean_text(extra.get("doi") or meta.get("doi"))
        node_rows.append({
            "pmid": pmid,
            "is_seed": 1 if pmid in seed_set else 0,
            "min_generation": 0 if pmid in seed_set else min_generation.get(pmid, ""),
            "generations": "|".join(str(g) for g in sorted(node_generations.get(pmid, set()))),
            "title": clean_text(extra.get("title") or meta.get("title")),
            "year": meta.get("year", ""),
            "journal": clean_text(meta.get("journal")),
            "doi": doi,
            "citation_count": meta.get("citation_count", ""),
            "seed_count": len(node_seeds.get(pmid, set())),
            "parent_count": len(node_parents.get(pmid, set())),
            "score": score.get("score", ""),
            "priority": score.get("priority", "Seed"),
            "matched_genes": "|".join(score.get("matched_genes", [])),
            "matched_terms": "|".join(score.get("matched_terms", [])),
            "has_abstract": 1 if clean_text(extra.get("abstract")) else 0,
            "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            "icite_url": f"https://icite.od.nih.gov/analysis?pmids={pmid}",
        })
    write_csv(OUT_DIR / "all_nodes.csv.gz", node_fields, node_rows, gzip_output=True)

    edge_fields = ["generation", "parent_pmid", "child_pmid", "parent_is_seed", "child_is_seed", "child_min_generation"]
    edge_rows = []
    for generation, edges in [(1, g1_edges), (2, g2_edges), (3, g3_edges)]:
        for parent, child in sorted(edges, key=lambda pair: (int(pair[0]), int(pair[1]))):
            edge_rows.append({
                "generation": generation,
                "parent_pmid": parent,
                "child_pmid": child,
                "parent_is_seed": 1 if parent in seed_set else 0,
                "child_is_seed": 1 if child in seed_set else 0,
                "child_min_generation": 0 if child in seed_set else min_generation.get(child, ""),
            })
    write_csv(OUT_DIR / "citation_edges.csv.gz", edge_fields, edge_rows, gzip_output=True)

    candidate_fields = [
        "rank", "pmid", "priority", "score", "min_generation", "title", "year",
        "journal", "doi", "citation_count", "seed_count", "parent_count",
        "matched_genes", "matched_terms", "representative_seed_pmid",
        "representative_seed_genes", "representative_path", "abstract_snippet",
        "pubmed_url", "icite_url",
    ]
    candidate_rows = []
    for pmid in candidate_pmids:
        meta = metadata.get(pmid, {})
        extra = abstracts.get(pmid, {})
        score = scores[pmid]
        supports = sorted(node_seeds.get(pmid, set()), key=int)
        representative_seed = supports[0] if supports else ""
        path_generation, path = representative_path(pmid, representative_seed, g1_reach, g2_reach, g3_reach) if representative_seed else (99, [])
        abstract = clean_text(extra.get("abstract"))
        candidate_rows.append({
            "rank": rank_by_pmid[pmid],
            "pmid": pmid,
            "priority": score["priority"],
            "score": score["score"],
            "min_generation": min_generation.get(pmid, path_generation),
            "title": score["title"] or clean_text(meta.get("title")),
            "year": meta.get("year", ""),
            "journal": clean_text(meta.get("journal")),
            "doi": clean_text(extra.get("doi") or meta.get("doi")),
            "citation_count": meta.get("citation_count", ""),
            "seed_count": len(supports),
            "parent_count": len(node_parents.get(pmid, set())),
            "matched_genes": "|".join(score["matched_genes"]),
            "matched_terms": "|".join(score["matched_terms"]),
            "representative_seed_pmid": representative_seed,
            "representative_seed_genes": "|".join(seed_genes.get(representative_seed, [])),
            "representative_path": " > ".join(path),
            "abstract_snippet": abstract[:900],
            "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            "icite_url": f"https://icite.od.nih.gov/analysis?pmids={pmid}",
        })
    write_csv(OUT_DIR / "candidate_papers.csv", candidate_fields, candidate_rows)

    abstract_fields = ["pmid", "title", "abstract"]
    write_csv(
        OUT_DIR / "candidate_abstracts.csv.gz",
        abstract_fields,
        (
            {
                "pmid": pmid,
                "title": scores[pmid]["title"],
                "abstract": scores[pmid]["abstract"],
            }
            for pmid in candidate_pmids
        ),
        gzip_output=True,
    )

    link_fields = ["candidate_pmid", "candidate_rank", "seed_pmid", "seed_genes", "generation", "representative_path"]
    link_rows = []
    for pmid in candidate_pmids:
        for seed in sorted(node_seeds.get(pmid, set()), key=int):
            generation, path = representative_path(pmid, seed, g1_reach, g2_reach, g3_reach)
            link_rows.append({
                "candidate_pmid": pmid,
                "candidate_rank": rank_by_pmid[pmid],
                "seed_pmid": seed,
                "seed_genes": "|".join(seed_genes.get(seed, [])),
                "generation": generation,
                "representative_path": " > ".join(path),
            })
    write_csv(OUT_DIR / "candidate_seed_paths.csv.gz", link_fields, link_rows, gzip_output=True)

    seed_summary_fields = [
        "seed_pmid", "seed_genes", "generation_1_unique", "generation_2_unique",
        "generation_3_unique_selective", "shortlisted_candidates", "priority_A",
        "priority_B", "priority_C", "context", "top_candidate_pmids",
    ]
    seed_summary_rows = []
    for seed in seeds:
        g1_for_seed = {node for node, support in g1_reach.items() if seed in support}
        g2_for_seed = {node for node, support in g2_reach.items() if seed in support}
        g3_for_seed = {node for node, support in g3_reach.items() if seed in support}
        candidates_for_seed = [p for p in candidate_pmids if seed in node_seeds.get(p, set())]
        counts = Counter(scores[p]["priority"] for p in candidates_for_seed)
        seed_summary_rows.append({
            "seed_pmid": seed,
            "seed_genes": "|".join(seed_genes.get(seed, [])),
            "generation_1_unique": len(g1_for_seed),
            "generation_2_unique": len(g2_for_seed),
            "generation_3_unique_selective": len(g3_for_seed),
            "shortlisted_candidates": len(candidates_for_seed),
            "priority_A": counts.get("A", 0),
            "priority_B": counts.get("B", 0),
            "priority_C": counts.get("C", 0),
            "context": counts.get("Context", 0),
            "top_candidate_pmids": "|".join(candidates_for_seed[:10]),
        })
    write_csv(OUT_DIR / "seed_summary.csv", seed_summary_fields, seed_summary_rows)

    priority_counts = Counter(scores[p]["priority"] for p in candidate_pmids)
    unique_crossrefs = all_crossref_nodes - seed_set
    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "seed_papers": len(seeds),
        "seed_records_retrieved_by_icite": len(seed_records),
        "generation_1": {
            "scope": "complete PubMed-indexed backward references of all seeds",
            "unique_edges": len(g1_edges),
            "unique_nodes": len(g1_nodes),
        },
        "generation_2": {
            "scope": "complete PubMed-indexed backward references of all generation-1 papers",
            "unique_edges": len(g2_edges),
            "unique_nodes": len(g2_nodes),
        },
        "generation_3": {
            "scope": "selective backward references of relevance-ranked generation-2 frontier",
            "frontier_papers": len(frontier),
            "unique_edges": len(g3_edges),
            "unique_nodes": len(g3_nodes),
            "frontier_limit": MAX_G3_FRONTIER,
        },
        "unique_non_seed_cross_reference_papers": len(unique_crossrefs),
        "all_unique_nodes_including_seeds": len(all_nodes),
        "shortlisted_non_seed_candidates": len(candidate_pmids),
        "candidate_priority_counts": dict(priority_counts),
        "abstract_records_retrieved": len(abstracts),
        "method": {
            "direction": "backward citation snowballing",
            "deduplication_key": "PMID",
            "reference_source": "NIH iCite / NIH Open Citation Collection",
            "metadata_and_abstract_source": "PubMed / NCBI E-utilities",
            "generation_3_selection": "CalorType keyword, seed-gene, and multi-seed graph support ranking",
            "important_limitations": [
                "Only references resolvable to PubMed IDs are represented; non-PubMed citations are absent.",
                "Generation 3 is deliberately selective rather than exhaustive to control combinatorial expansion.",
                "Relevance ranking is a high-recall triage heuristic and does not replace full-text CalorType extraction.",
            ],
        },
        "runtime_seconds": round(time.time() - start_time, 1),
    }
    with (OUT_DIR / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)

    readme = f"""# CalorType set 2: citation snowball output

Generated: {summary['generated_at_utc']}

- Seeds: {len(seeds):,}
- Generation 1: {len(g1_nodes):,} unique PubMed-indexed references ({len(g1_edges):,} edges)
- Generation 2: {len(g2_nodes):,} unique PubMed-indexed references ({len(g2_edges):,} edges)
- Generation 3: {len(g3_nodes):,} unique references from a selective frontier of {len(frontier):,} generation-2 papers
- Shortlisted candidates: {len(candidate_pmids):,}

## Files

- `candidate_papers.csv`: ranked review table with representative citation path.
- `seed_summary.csv`: per-seed expansion and shortlist counts.
- `all_nodes.csv.gz`: all seed and cross-reference nodes.
- `citation_edges.csv.gz`: normalized graph edges.
- `candidate_seed_paths.csv.gz`: seed-specific paths for shortlisted candidates.
- `candidate_abstracts.csv.gz`: full available abstracts for shortlisted candidates.
- `summary.json`: machine-readable methodology and counts.

## Interpretation

Generations 1 and 2 are exhaustive within the NIH Open Citation Collection's PubMed-to-PubMed links. Generation 3 is selective and follows only the most CalorType-relevant generation-2 papers. Relevance scores are triage signals based on temperature/thermal phenotype terms, mutation/variant terminology, seed-gene mentions, and independent support from multiple seeds.
"""
    (OUT_DIR / "README.md").write_text(readme, encoding="utf-8")

    if BUNDLE.exists():
        BUNDLE.unlink()
    with zipfile.ZipFile(BUNDLE, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path in sorted(OUT_DIR.rglob("*")):
            if path.is_file():
                archive.write(path, arcname=path.relative_to(OUT_DIR.parent))

    log(f"Wrote bundle: {BUNDLE} ({BUNDLE.stat().st_size / 1024 / 1024:.1f} MiB)")
    log(f"Completed in {(time.time() - start_time) / 60:.1f} minutes")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise
    except Exception as exc:
        log(f"FATAL: {exc}")
        raise
