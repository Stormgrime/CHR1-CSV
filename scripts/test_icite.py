import json
import urllib.request

url = "https://icite.od.nih.gov/api/pubs?pmids=11038181&refs=true"
with urllib.request.urlopen(url, timeout=60) as response:
    data = json.load(response)
record = data["data"][0]
summary = {
    "top_level_keys": sorted(data.keys()),
    "record_keys": sorted(record.keys()),
    "pmid": record.get("pmid"),
    "reference_count": len(record.get("references") or record.get("citedPmids") or []),
    "references": record.get("references") or record.get("citedPmids") or [],
    "cited_by_count": len(record.get("cited_by") or record.get("citedByPmids") or []),
}
print(json.dumps(summary, indent=2))
with open("test_icite.json", "w", encoding="utf-8") as fh:
    json.dump(summary, fh, indent=2)
