#!/usr/bin/env python3
from pathlib import Path

source_path = Path(__file__).with_name("expand.py")
source = source_path.read_text(encoding="utf-8")

old_expansion = '''# Generation 1 is exhaustive across every PMID-resolved reference of the 240 seeds.
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
'''

new_expansion = '''# Generation 1 is exhaustive across every PMID-resolved reference of the 240 seeds.
g1_nodes, g1_by_root = expand(sorted(seed_pmids, key=int), 1)
g1_scores = enrich_and_score(g1_nodes, 1)

# Every unique non-seed G1 paper is expanded, making Generation 2 exhaustive.
g1_parents = sorted(g1_nodes.difference(seed_pmids), key=int)
g1_reasons = {pmid: ["exhaustive_generation_2"] for pmid in g1_parents}
for pmid in g1_parents:
    expansion_rows.append({"generation": 1, "pmid": pmid, "reason": g1_reasons.get(pmid, [])})

seen_before_g2 = set(seed_pmids).union(g1_nodes)
g2_targets, g2_by_root = expand(g1_parents, 2)
g2_nodes = g2_targets.difference(seen_before_g2)
g2_scores = enrich_and_score(g2_nodes, 2)

# Generation 3 is relevance-guided. Previously expanded papers are excluded from the parent set.
g2_parents, g2_reasons = select_for_expansion(2, g2_nodes, g2_scores, g2_by_root, global_cap=900, per_root_top=6)
for pmid in g2_parents:
    expansion_rows.append({"generation": 2, "pmid": pmid, "reason": g2_reasons.get(pmid, [])})

# Generation 3 is a terminal frontier; only newly discovered papers are labelled G3.
seen_before_g3 = seen_before_g2.union(g2_nodes)
g3_targets, _ = expand(g2_parents, 3)
g3_nodes = g3_targets.difference(seen_before_g3)
g3_scores = enrich_and_score(g3_nodes, 3)
'''

old_caps = '''    "expansion_caps": {"generation_1_parents": 1400, "generation_2_parents": 450, "per_seed_top_generation_1": 12, "per_seed_top_generation_2": 4},
    "notes": [
        "Generation 1 is exhaustive across all PubMed-resolved references of the 240 seeds.",
        "Generations 2 and 3 are relevance-guided to prevent combinatorial explosion.",
        "Relevance uses seed-gene matches, thermal terms, variant-mechanism terms, shared roots, shared parents, and citation count.",
        "Coverage is limited to references resolved to PMIDs in the NIH Open Citation Collection.",
    ],
'''

new_caps = '''    "expansion_caps": {"generation_1_parents": "all unique non-seed G1 papers", "generation_2_parents": 900, "per_seed_top_generation_2": 6},
    "notes": [
        "Generation 1 is exhaustive across all PubMed-resolved references of the 240 seeds.",
        "Generation 2 is exhaustive across every unique non-seed Generation-1 paper with an iCite record.",
        "Generation 3 is relevance-guided to prevent combinatorial explosion.",
        "Relevance uses seed-gene matches, thermal terms, variant-mechanism terms, shared roots, shared parents, and citation count.",
        "Generation labels are minimum-distance frontiers; previously seen papers are not relabelled in deeper generations.",
        "Coverage is limited to references resolved to PMIDs in the NIH Open Citation Collection.",
    ],
'''

if old_expansion not in source:
    raise RuntimeError("Could not locate expansion block in expand.py")
if old_caps not in source:
    raise RuntimeError("Could not locate summary block in expand.py")

source = source.replace(old_expansion, new_expansion, 1).replace(old_caps, new_caps, 1)
namespace = {"__name__": "__main__", "__file__": str(source_path)}
exec(compile(source, str(source_path), "exec"), namespace)
