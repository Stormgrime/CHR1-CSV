#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "seed_set2.json"
OUT = ROOT / "output"
OUT.mkdir(parents=True, exist_ok=True)

ICITE = "https://icite.od.nih.gov/api/pubs"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
USER_AGENT = "CalorType-cross-reference-expansion/1.0 (research use)"

THERMAL_PATTERNS = {
    "temperature_sensitive": r"\btemperature[- ]sensitive\b|\bthermosensitiv\w*\b",
    "cold_sensitive": r"\bcold[- ]sensitive\b|\bcold[- ]induc\w*\b",
    "heat_sensitive": r"\bheat[- ]sensitive\b|\bheat[- ]intoler\w*\b",
    "thermolabile": r"\bthermolabil\w*\b|\bthermal instability\b",
    "permissive_temperature": r"\bnon[- ]?permissive temperature\b|\bpermissive temperature\b",
    "fever_febrile": r"\bfebrile\b|\bfever[- ]?(?:induced|triggered|sensitive)?\b",
    "hyperthermia": r"\bhyperthermi\w*\b",
    "hypothermia": r"\bhypothermi\w*\b",
    "heat_shock": r"\bheat shock\b|\bheat[- ]stress\b",
    "temperature_dependence": r"\btemperature[- ]depend\w*\b|\btemperature[- ]induc\w*\b|\btemperature[- ]trigger\w*\b",
    "thermal": r"\bthermal\b|\btemperature\b",
}
MECHANISM_PATTERNS = {
    "mutation": r"\bmutat\w*\b|\bmutant\w*\b",
    "variant": r"\bvariant\w*\b|\bpolymorphism\w*\b",
    "allele": r"\ballel\w*\b",
    "missense": r"\bmissense\b|\bsubstitut\w*\b",
    "folding_stability": r"\bprotein fold\w*\b|\bstabilit\w*\b|\bmisfold\w*\b|\bdenatur\w*\b",
    "trafficking": r"\btraffic\w*\b|\blocali[sz]\w*\b",
    "activity_function": r"\bactivity\b|\bfunction\w*\b|\bkinetic\w*\b",
    "phenotype": r"\bphenotyp\w*\b",
}


def chunks(items: list[str], n: int) -> Iterable[list[str]]:
    for i in range(0, len(items), n):
        yield items[i:i + n]


def get_bytes(url: str, retries: int = 6, pause: float = 1.0) -> bytes:
    last: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "*/*"})
            with urllib.request.urlopen(req, timeout=90) as response:
                return response.read()
        except Exception as exc:
            last = exc
            time.sleep(pause * (2 ** min(attempt, 4)))
    raise RuntimeError(f"Failed after {retries} attempts: {url}") from last


def fetch_icite(pmids: Iterable[str]) -> dict[str, dict[str, Any]]:
    ids = sorted({str(x) for x in pmids if str(x).isdigit()}, key=int)
    out: dict[str, dict[str, Any]] = {}
    for idx, batch in enumerate(chunks(ids, 180), 1):
        query = urllib.parse.urlencode({"pmids": ",".join(batch), "refs": "true"})
        payload = json.loads(get_bytes(f"{ICITE}?{query}").decode("utf-8"))
        for record in payload.get("data", []):
            out[str(record["pmid"])] = record
        print(f"iCite batch {idx}: {len(batch)} requested, {len(out)} cumulative", flush=True)
        time.sleep(0.15)
    return out


def fetch_abstracts(pmids: Iterable[str]) -> dict[str, str]:
    ids = sorted({str(x) for x in pmids if str(x).isdigit()}, key=int)
    out: dict[str, str] = {}
    for idx, batch in enumerate(chunks(ids, 120), 1):
        query = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(batch), "retmode": "xml"})
        root = ET.fromstring(get_bytes(f"{EFETCH}?{query}"))
        for article in root.findall(".//PubmedArticle"):
            pmid_el = article.find(".//MedlineCitation/PMID")
            if pmid_el is None or not pmid_el.text:
                continue
            parts = []
            for element in article.findall(".//Article/Abstract/AbstractText"):
                text = " ".join("".join(element.itertext()).split())
                if not text:
                    continue
                label = element.attrib.get("Label")
                parts.append(f"{label}: {text}" if label else text)
            out[pmid_el.text.strip()] = " ".join(parts)
        print(f"PubMed abstract batch {idx}: {len(batch)} requested, {len(out)} cumulative", flush=True)
        time.sleep(0.36)
    return out


def norm_gene(gene: str) -> str:
    return re.sub(r"[^A-Za-z0-9]", "", gene or "").upper()


def gene_hits(text: str, genes: Iterable[str]) -> list[str]:
    hits = []
    for gene in sorted({g for g in genes if g}):
        token = norm_gene(gene)
        if len(token) < 3:
            continue
        if re.search(rf"(?<![A-Za-z0-9]){re.escape(token)}(?:P)?(?![A-Za-z0-9])", text, re.I):
            hits.append(gene)
    return hits


def pattern_hits(text: str, patterns: dict[str, str]) -> list[str]:
    return [name for name, pattern in patterns.items() if re.search(pattern, text, re.I)]


def safe_refs(record: dict[str, Any] | None) -> list[str]:
    if not record:
        return []
    refs = record.get("references") or record.get("citedPmids") or []
    return [str(x) for x in refs if str(x).isdigit()]


def title_of(record: dict[str, Any] | None) -> str:
    return (record or {}).get("title") or ""


def year_of(record: dict[str, Any] | None) -> Any:
    return (record or {}).get("year") or (record or {}).get("pubYear") or ""


def journal_of(record: dict[str, Any] | None) -> str:
    return (record or {}).get("journal") or (record or {}).get("journalNameIso") or ""


def doi_of(record: dict[str, Any] | None) -> str:
    return (record or {}).get("doi") or ""


def citation_count(record: dict[str, Any] | None) -> int:
    return int((record or {}).get("citation_count") or (record or {}).get("citedByPmidCount") or 0)


def relevance(pmid: str, record: dict[str, Any] | None, abstract: str, roots: set[str], root_genes: dict[str, set[str]], parent_count: int) -> dict[str, Any]:
    text = f"{title_of(record)}\n{abstract}"
    genes = set()
    for root in roots:
        genes.update(root_genes.get(root, set()))
    gh = gene_hits(text, genes)
    th = pattern_hits(text, THERMAL_PATTERNS)
    mh = pattern_hits(text, MECHANISM_PATTERNS)
    strong_thermal = [x for x in th if x != "thermal"]
    root_count = len(roots)
    score = 0.0
    score += min(14.0, 5.0 * len(strong_thermal) + (1.5 if "thermal" in th else 0.0))
    score += min(14.0, 5.0 * len(gh))
    score += min(5.0, 1.0 * len(mh))
    score += min(8.0, 2.0 * math.log2(1 + root_count))
    score += min(4.0, 1.25 * math.log2(1 + parent_count))
    score += min(2.0, 0.35 * math.log2(1 + citation_count(record)))
    if (strong_thermal and gh) or score >= 13:
        tier = "A"
    elif strong_thermal or gh or root_count >= 3 or score >= 7:
        tier = "B"
    else:
        tier = "C"
    return {
        "score": round(score, 3),
        "tier": tier,
        "gene_hits": gh,
        "thermal_hits": th,
        "mechanism_hits": mh,
        "root_count": root_count,
        "parent_count": parent_count,
    }


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            cooked = {}
            for key in fields:
                value = row.get(key, "")
                if isinstance(value, (list, set, tuple)):
                    value = ";".join(str(x) for x in sorted(value))
                elif isinstance(value, dict):
                    value = json.dumps(value, ensure_ascii=False, sort_keys=True)
                cooked[key] = value
            writer.writerow(cooked)


with INPUT.open("r", encoding="utf-8", newline="") as handle:
    seeds = []
    for row in csv.DictReader(handle, delimiter="\t"):
        seeds.append({
            "pmid": str(row["pmid"]).strip(),
            "genes": [g.strip() for g in (row.get("genes") or "").split(";") if g.strip()],
        })
seed_by_pmid = {str(seed["pmid"]): seed for seed in seeds}
seed_pmids = set(seed_by_pmid)
root_genes: dict[str, set[str]] = {pmid: {g for g in seed.get("genes", []) if g} for pmid, seed in seed_by_pmid.items()}

all_meta: dict[str, dict[str, Any]] = {}
all_abstracts: dict[str, str] = {}
node_roots: dict[str, set[str]] = defaultdict(set)
node_parents: dict[str, set[str]] = defaultdict(set)
node_min_generation: dict[str, int] = {}
edges: list[dict[str, Any]] = []
score_by_generation: dict[int, dict[str, dict[str, Any]]] = {}
expansion_rows: list[dict[str, Any]] = []

all_meta.update(fetch_icite(seed_pmids))
for pmid in seed_pmids:
    node_roots[pmid].add(pmid)
    node_min_generation[pmid] = 0


def expand(parents: list[str], generation: int) -> tuple[set[str], dict[str, set[str]]]:
    targets: set[str] = set()
    targets_by_root: dict[str, set[str]] = defaultdict(set)
    for parent in parents:
        roots = set(node_roots.get(parent, set()))
        for child in safe_refs(all_meta.get(parent)):
            targets.add(child)
            node_roots[child].update(roots)
            node_parents[child].add(parent)
            node_min_generation[child] = min(generation, node_min_generation.get(child, generation))
            for root in roots:
                targets_by_root[root].add(child)
            edges.append({"generation": generation, "root_pmids": sorted(roots), "parent_pmid": parent, "child_pmid": child})
    return targets, targets_by_root


def enrich_and_score(nodes: set[str], generation: int) -> dict[str, dict[str, Any]]:
    missing_meta = nodes.difference(all_meta)
    if missing_meta:
        all_meta.update(fetch_icite(missing_meta))
    missing_abstracts = nodes.difference(all_abstracts)
    if missing_abstracts:
        all_abstracts.update(fetch_abstracts(missing_abstracts))
    scores = {}
    for pmid in nodes:
        scores[pmid] = relevance(pmid, all_meta.get(pmid), all_abstracts.get(pmid, ""), node_roots.get(pmid, set()), root_genes, len(node_parents.get(pmid, set())))
    score_by_generation[generation] = scores
    return scores


def select_for_expansion(generation: int, candidates: set[str], scores: dict[str, dict[str, Any]], edges_by_root: dict[str, set[str]], global_cap: int, per_root_top: int) -> tuple[list[str], dict[str, list[str]]]:
    selected: set[str] = set()
    reasons: dict[str, list[str]] = defaultdict(list)
    for pmid in candidates:
        tier = scores.get(pmid, {}).get("tier")
        if tier in {"A", "B"}:
            selected.add(pmid)
            reasons[pmid].append(f"tier_{tier}")
    for root, nodes in edges_by_root.items():
        ranked = sorted(
            (node for node in nodes if node in candidates),
            key=lambda node: (
                scores.get(node, {}).get("score", 0),
                scores.get(node, {}).get("root_count", 0),
                scores.get(node, {}).get("parent_count", 0),
                citation_count(all_meta.get(node)),
            ),
            reverse=True,
        )[:per_root_top]
        for node in ranked:
            selected.add(node)
            reasons[node].append(f"top_{per_root_top}_for_root_{root}")
    ranked_all = sorted(
        selected,
        key=lambda node: (
            scores.get(node, {}).get("score", 0),
            scores.get(node, {}).get("root_count", 0),
            scores.get(node, {}).get("parent_count", 0),
            citation_count(all_meta.get(node)),
        ),
        reverse=True,
    )[:global_cap]
    reasons = {node: reasons[node] for node in ranked_all}
    print(f"Generation {generation}: selected {len(ranked_all)} parents for next expansion", flush=True)
    return ranked_all, reasons


# Generation 1 is exhaustive across every PMID-resolved reference of the 240 seeds.
g1_nodes, g1_by_root = expand(sorted(seed_pmids, key=int), 1)
g1_scores = enrich_and_score(g1_nodes, 1)
g1_parents, g1_reasons = select_for_expansion(1, g1_nodes, g1_scores, g1_by_root, global_cap=1400, per_root_top=12)
for pmid in g1_parents:
    expansion_rows.append({"generation": 1, "pmid": pmid, "reason": g1_reasons.get(pmid, [])})

# Generation 2 follows references of the relevance-guided G1 parent set.
g2_nodes, g2_by_root = expand(g1_parents, 2)
g2_scores = enrich_and_score(g2_nodes, 2)
g2_parents, g2_reasons = select_for_expansion(2, g2_nodes, g2_scores, g2_by_root, global_cap=450, per_root_top=4)
for pmid in g2_parents:
    expansion_rows.append({"generation": 2, "pmid": pmid, "reason": g2_reasons.get(pmid, [])})

# Generation 3 is a terminal frontier; it is scored but not expanded again.
g3_nodes, _ = expand(g2_parents, 3)
g3_scores = enrich_and_score(g3_nodes, 3)

all_nodes = set(node_min_generation)
final_scores: dict[str, dict[str, Any]] = {}
for pmid in all_nodes:
    final_scores[pmid] = relevance(pmid, all_meta.get(pmid), all_abstracts.get(pmid, ""), node_roots.get(pmid, set()), root_genes, len(node_parents.get(pmid, set())))

candidate_rows: list[dict[str, Any]] = []
for pmid in all_nodes:
    generation = node_min_generation.get(pmid, 99)
    record = all_meta.get(pmid)
    score = final_scores.get(pmid) or {"score": 0, "tier": "", "gene_hits": [], "thermal_hits": [], "mechanism_hits": []}
    roots = sorted(node_roots.get(pmid, set()), key=int)
    root_gene_set = set()
    for root in roots:
        root_gene_set.update(root_genes.get(root, set()))
    candidate_rows.append({
        "pmid": pmid,
        "min_generation": generation,
        "is_seed": pmid in seed_pmids,
        "title": title_of(record),
        "year": year_of(record),
        "journal": journal_of(record),
        "doi": doi_of(record),
        "citation_count": citation_count(record),
        "relevance_score": score.get("score", 0),
        "relevance_tier": score.get("tier", ""),
        "root_count": len(roots),
        "parent_count": len(node_parents.get(pmid, set())),
        "root_pmids": roots,
        "root_genes": root_gene_set,
        "matched_genes": score.get("gene_hits", []),
        "thermal_hits": score.get("thermal_hits", []),
        "mechanism_hits": score.get("mechanism_hits", []),
        "abstract": all_abstracts.get(pmid, ""),
        "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
    })

candidate_rows.sort(key=lambda row: (row["is_seed"], row["relevance_tier"] == "A", row["relevance_tier"] == "B", row["relevance_score"], row["root_count"], row["parent_count"], row["citation_count"]), reverse=True)
priority_rows = [row for row in candidate_rows if not row["is_seed"] and (row["relevance_tier"] in {"A", "B"} or row["root_count"] >= 3)]
priority_rows.sort(key=lambda row: (row["relevance_tier"] == "A", row["relevance_tier"] == "B", row["relevance_score"], row["root_count"], row["parent_count"], row["citation_count"]), reverse=True)

edge_rows = []
for edge in edges:
    parent = edge["parent_pmid"]
    child = edge["child_pmid"]
    edge_rows.append({
        "generation": edge["generation"],
        "root_pmids": edge["root_pmids"],
        "parent_pmid": parent,
        "parent_title": title_of(all_meta.get(parent)),
        "child_pmid": child,
        "child_title": title_of(all_meta.get(child)),
        "child_min_generation": node_min_generation.get(child, ""),
        "child_relevance_score": final_scores.get(child, {}).get("score", ""),
    })

seed_rows = []
for seed in seeds:
    pmid = seed["pmid"]
    record = all_meta.get(pmid)
    seed_rows.append({
        "pmid": pmid,
        "title": title_of(record),
        "year": year_of(record),
        "journal": journal_of(record),
        "genes": seed.get("genes", []),
        "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
    })

gene_rollup: dict[str, dict[str, Any]] = {}
for row in priority_rows:
    for gene in row.get("root_genes", []):
        data = gene_rollup.setdefault(gene, {"gene": gene, "priority_candidates": 0, "tier_A": 0, "tier_B": 0, "generation_1": 0, "generation_2": 0, "generation_3": 0, "top_score": 0.0, "top_pmids": []})
        data["priority_candidates"] += 1
        if row["relevance_tier"] == "A":
            data["tier_A"] += 1
        elif row["relevance_tier"] == "B":
            data["tier_B"] += 1
        generation_key = f"generation_{row['min_generation']}"
        if generation_key in data:
            data[generation_key] += 1
        data["top_score"] = max(data["top_score"], row["relevance_score"])
        if len(data["top_pmids"]) < 10:
            data["top_pmids"].append(row["pmid"])
gene_rows = sorted(gene_rollup.values(), key=lambda row: (row["tier_A"], row["priority_candidates"], row["top_score"]), reverse=True)

candidate_fields = ["pmid", "min_generation", "is_seed", "title", "year", "journal", "doi", "citation_count", "relevance_score", "relevance_tier", "root_count", "parent_count", "root_pmids", "root_genes", "matched_genes", "thermal_hits", "mechanism_hits", "abstract", "pubmed_url"]
write_csv(OUT / "seed_papers.csv", seed_rows, ["pmid", "title", "year", "journal", "genes", "pubmed_url"])
write_csv(OUT / "all_candidate_papers.csv", candidate_rows, candidate_fields)
write_csv(OUT / "priority_candidates.csv", priority_rows, candidate_fields)
write_csv(OUT / "citation_edges.csv", edge_rows, ["generation", "root_pmids", "parent_pmid", "parent_title", "child_pmid", "child_title", "child_min_generation", "child_relevance_score"])
write_csv(OUT / "expansion_parents.csv", expansion_rows, ["generation", "pmid", "reason"])
write_csv(OUT / "gene_summary.csv", gene_rows, ["gene", "priority_candidates", "tier_A", "tier_B", "generation_1", "generation_2", "generation_3", "top_score", "top_pmids"])

summary = {
    "method": "Backward citation snowballing through PubMed-indexed references in NIH iCite/OCC.",
    "seed_count": len(seed_pmids),
    "generation_1": {"unique_papers": len(g1_nodes), "edges": sum(1 for edge in edges if edge["generation"] == 1), "expanded_parents": len(g1_parents)},
    "generation_2": {"unique_papers": len(g2_nodes), "edges": sum(1 for edge in edges if edge["generation"] == 2), "expanded_parents": len(g2_parents)},
    "generation_3": {"unique_papers": len(g3_nodes), "edges": sum(1 for edge in edges if edge["generation"] == 3), "expanded_parents": 0},
    "all_unique_nodes_including_seeds": len(all_nodes),
    "priority_candidate_count": len(priority_rows),
    "tier_A_count": sum(row["relevance_tier"] == "A" for row in priority_rows),
    "tier_B_count": sum(row["relevance_tier"] == "B" for row in priority_rows),
    "expansion_caps": {"generation_1_parents": 1400, "generation_2_parents": 450, "per_seed_top_generation_1": 12, "per_seed_top_generation_2": 4},
    "notes": [
        "Generation 1 is exhaustive across all PubMed-resolved references of the 240 seeds.",
        "Generations 2 and 3 are relevance-guided to prevent combinatorial explosion.",
        "Relevance uses seed-gene matches, thermal terms, variant-mechanism terms, shared roots, shared parents, and citation count.",
        "Coverage is limited to references resolved to PMIDs in the NIH Open Citation Collection.",
    ],
}
(OUT / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
print(json.dumps(summary, indent=2), flush=True)
