import json
import urllib.request

url = "https://icite.od.nih.gov/api/pubs?pmids=11038181"
with urllib.request.urlopen(url, timeout=60) as response:
    data = json.load(response)
print(json.dumps(data, indent=2))
with open("test_icite.json", "w", encoding="utf-8") as fh:
    json.dump(data, fh, indent=2)
