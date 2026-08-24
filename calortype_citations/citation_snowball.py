#!/usr/bin/env python3
from __future__ import annotations
import csv, json, math, re, time, urllib.parse, urllib.request, urllib.error
import xml.etree.ElementTree as ET
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
INP=ROOT/"calortype_citations"/"seed_records.json"
OUT=ROOT/"results"; OUT.mkdir(exist_ok=True)
ICITE="https://icite.od.nih.gov/api/pubs"
EFETCH="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
UA="CalorType-citation-snowball/1.0"
G1_CAP,G2_CAP=1200,800
G1_THR,G2_THR=4.0,5.0
PRIO=4.0
MAXPATH=3

STRONG=[(r"\btemperature[-\s]?sensitive\b",6,"temperature-sensitive"),
(r"\bthermosensitiv\w*\b",5.5,"thermosensitive"),(r"\bthermolabil\w*\b",5,"thermolabile"),
(r"\bheat[-\s]?sensitiv\w*\b",5,"heat-sensitive"),(r"\bcold[-\s]?sensitiv\w*\b",5,"cold-sensitive"),
(r"\bnon[-\s]?permissive temperature\b",6,"nonpermissive temperature"),
(r"\bpermissive temperature\b",5,"permissive temperature"),
(r"\btemperature[-\s]?(dependent|induced|triggered)\b",4,"temperature-dependent"),
(r"\bfever[-\s]?(induced|triggered|sensitive|associated)\b",5,"fever-triggered"),
(r"\bfebrile\b",4,"febrile"),(r"\bhyperthermi\w*\b",4,"hyperthermia"),
(r"\bhypothermi\w*\b",4,"hypothermia"),(r"\bmalignant hyperthermia\b",5,"malignant hyperthermia")]
TEMP=[(r"\bheat shock\b",2,"heat shock"),(r"\bcold shock\b",2,"cold shock"),
(r"\bthermal stress\b",2,"thermal stress"),(r"\bthermal stability\b",2.5,"thermal stability"),
(r"\btemperature\b",1,"temperature"),(r"\bthermal\b",1,"thermal"),
(r"\bheat\b",.75,"heat"),(r"\bcold\b",.75,"cold"),(r"\bfever\b",1.5,"fever")]
VAR=[(r"\bmutant\w*\b",2,"mutant"),(r"\bmutation\w*\b",2,"mutation"),
(r"\bvariant\w*\b",2,"variant"),(r"\ballele\w*\b",2,"allele"),
(r"\bpolymorphism\w*\b",1.5,"polymorphism"),(r"\bmissense\b",2,"missense"),
(r"\bnonsense\b",2,"nonsense"),(r"\bsubstitution\w*\b",1.5,"substitution"),
(r"\bdeletion\w*\b",1.5,"deletion"),(r"\bconditional allele\b",3,"conditional allele")]
PHENO=[(r"\bepilep\w*\b",1.25,"epilepsy"),(r"\bseizure\w*\b",1.25,"seizure"),
(r"\bchannelopath\w*\b",1.25,"channelopathy"),(r"\barrhythmi\w*\b",1,"arrhythmia"),
(r"\bmyopath\w*\b",1,"myopathy"),(r"\bataxi\w*\b",1,"ataxia"),
(r"\bdystoni\w*\b",1,"dystonia"),(r"\bmisfold\w*\b",1,"misfolding"),
(r"\bprotein stabil\w*\b",1,"protein stability"),(r"\benzyme activit\w*\b",.75,"enzyme activity")]
BADGENE={"HI","CI","WT","DNA","RNA","ATP","GTP","HIS3","NONE","N/A","DCASE","APM1","CCA1"}

def bat(xs,n):
    xs=list(xs)
    for i in range(0,len(xs),n): yield xs[i:i+n]

def get(url,accept="application/json"):
    d=1
    for k in range(6):
        try:
            q=urllib.request.Request(url,headers={"User-Agent":UA,"Accept":accept})
            with urllib.request.urlopen(q,timeout=120) as r: z=r.read()
            time.sleep(.36); return z
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError):
            if k==5: raise
            time.sleep(d); d=min(30,d*2)

def ids(v):
    if not v:return []
    if isinstance(v,str):v=re.split(r"[\s,;]+",v)
    return list(dict.fromkeys(str(x).strip() for x in v if str(x).strip().isdigit()))

def icite(pmids):
    out={}
    pmids=sorted(set(pmids),key=int)
    fields="pmid,title,year,journal,doi,citation_count,relative_citation_ratio,nih_percentile,is_research_article,references"
    for i,b in enumerate(bat(pmids,200),1):
        u=ICITE+"?"+urllib.parse.urlencode({"pmids":",".join(b),"fl":fields})
        p=json.loads(get(u).decode()); rows=p.get("data",p if isinstance(p,list) else [])
        for r in rows:
            p=str(r.get("pmid",""))
            if p.isdigit():
                r["refs"]=ids(r.get("references") or r.get("citedPmids"));out[p]=r
        print("iCite",i,len(b),len(rows),flush=True)
    return out

def txt(e): return " ".join("".join(e.itertext()).split()) if e is not None else ""

def pubmed(pmids):
    out={}; pmids=sorted(set(pmids),key=int)
    for i,b in enumerate(bat(pmids,100),1):
        u=EFETCH+"?"+urllib.parse.urlencode({"db":"pubmed","id":",".join(b),"retmode":"xml",
        "tool":"CalorTypeCitationSnowball","email":"haikg96@gmail.com"})
        root=ET.fromstring(get(u,"application/xml"))
        for x in root.findall(".//PubmedArticle"):
            c=x.find("./MedlineCitation");a=c.find("./Article") if c is not None else None
            p=(c.findtext("./PMID") or "").strip() if c is not None else ""
            if not p.isdigit() or a is None:continue
            ab=" ".join(txt(z) for z in a.findall("./Abstract/AbstractText") if txt(z))
            yr=a.findtext("./Journal/JournalIssue/PubDate/Year") or a.findtext("./ArticleDate/Year") or ""
            if not yr:
                md=a.findtext("./Journal/JournalIssue/PubDate/MedlineDate") or ""
                m=re.search(r"(18|19|20)\d{2}",md);yr=m.group(0) if m else ""
            aid={z.attrib.get("IdType","").lower():txt(z) for z in x.findall("./PubmedData/ArticleIdList/ArticleId")}
            mesh=[txt(z) for z in c.findall("./MeshHeadingList/MeshHeading/DescriptorName") if txt(z)]
            auth=[]
            for z in a.findall("./AuthorList/Author"):
                n=z.findtext("./CollectiveName") or " ".join(filter(None,[z.findtext("./LastName"),z.findtext("./Initials")]))
                if n:auth.append(n)
            out[p]={"title":txt(a.find("./ArticleTitle")),"abstract":ab,"year":yr,
            "journal":txt(a.find("./Journal/Title")),"authors":"; ".join(auth),
            "doi":aid.get("doi",""),"pmcid":aid.get("pmc",""),"mesh":mesh}
        print("PubMed",i,len(b),flush=True)
    return out

def gpat(g):
    if not g or len(g)<3 or g.upper() in BADGENE:return None
    return re.compile(r"(?<![A-Za-z0-9])"+re.escape(g)+r"(?![A-Za-z0-9])",re.I)

def merge(p,I,P):
    i=I.get(p,{});m=P.get(p,{})
    return {"pmid":p,"title":m.get("title") or i.get("title",""),"abstract":m.get("abstract",""),
    "year":m.get("year") or i.get("year",""),"journal":m.get("journal") or i.get("journal",""),
    "authors":m.get("authors",""),"doi":m.get("doi") or i.get("doi",""),"pmcid":m.get("pmcid",""),
    "mesh":m.get("mesh",[]),"citation_count":i.get("citation_count",""),
    "relative_citation_ratio":i.get("relative_citation_ratio",""),"nih_percentile":i.get("nih_percentile",""),
    "is_research_article":i.get("is_research_article",""),"reference_count":len(i.get("refs",[])),
    "icite_found":bool(i),"pubmed_found":bool(m)}

def score(r,ags,gp):
    t=r["title"];a=r["abstract"];s=(t+"\n"+a+"\n"+" ".join(r["mesh"])).lower()
    q=0;terms=[];st=False;vs=False;gm=False;genes=[]
    for pat,n,l in STRONG:
        if re.search(pat,s):q+=n;terms.append(l);st=True
    for pat,n,l in TEMP:
        if re.search(pat,s):q+=n;terms.append(l)
    for pat,n,l in VAR:
        if re.search(pat,s):q+=n;terms.append(l);vs=True
    for pat,n,l in PHENO:
        if re.search(pat,s):q+=n;terms.append(l)
    for g in sorted(ags):
        p=gp.get(g)
        if p and p.search(t):q+=5;genes.append(g);gm=True
        elif p and p.search(a+" "+" ".join(r["mesh"])):q+=3.5;genes.append(g);gm=True
    if not gm:
        for g,p in gp.items():
            if p.search(t):q+=2;genes.append(g);break
            if p.search(a):q+=1;genes.append(g);break
    at=st or any(x in terms for x in ["temperature","thermal","heat","cold","fever","heat shock","cold shock","thermal stress","thermal stability"])
    if at and vs:q+=2;terms.append("temperature × variant")
    if st and gm:q+=2;terms.append("temperature × ancestral gene")
    if gm and vs:q+=1;terms.append("gene × variant")
    gr="High" if q>=10 or(q>=8 and st and gm) else "Medium" if q>=PRIO else "Low" if q>=2 else "Background"
    return round(q,2),gr,sorted(set(terms)),sorted(set(genes)),st,vs,gm

def choose(cands,N,roots,parents,thr,cap,perseed):
    e=set()
    for p in cands:
        r=N[p];q=r["relevance_score"];rc=len(roots[p]);pc=len(parents[p])
        if q>=thr or(r["ancestor_gene_match"] and q>=thr-1.5) or(r["strong_temperature_signal"] and q>=thr-1.5) or(rc>=2 and pc>=2 and q>=max(2,thr-2)):e.add(p)
    ps=defaultdict(list)
    for p in cands:
        r=N[p]
        try:c=float(r["citation_count"] or 0)
        except:c=0
        z=r["relevance_score"]+.6*min(len(roots[p]),5)+.2*math.log1p(c)
        for s in roots[p]:ps[s].append((z,p))
    for L in ps.values():
        for _,p in sorted(L,reverse=True)[:perseed]:
            if N[p]["relevance_score"]>=1.5:e.add(p)
    def key(p):
        r=N[p]
        try:c=float(r["citation_count"] or 0)
        except:c=0
        return r["relevance_score"]+.75*min(len(roots[p]),8)+.35*min(len(parents[p]),8)+.2*math.log1p(c)
    return sorted(e,key=key,reverse=True)[:cap]

def wcsv(path,rows,fields):
    with path.open("w",newline="",encoding="utf-8-sig") as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore");w.writeheader()
        for r in rows:
            w.writerow({k:"; ".join(map(str,r.get(k,[]))) if isinstance(r.get(k), (list,set,tuple)) else ("Yes" if r.get(k) is True else "No" if r.get(k) is False else r.get(k,"")) for k in fields})

def main():
    seeds=json.loads(INP.read_text()); S={x["pmid"]:x for x in seeds}; sp=set(S)
    gps={g:p for g in {g for x in seeds for g in x["genes"]} if (p:=gpat(g))}
    roots=defaultdict(set);genes=defaultdict(set);parents=defaultdict(set);paths=defaultdict(list);ming={p:0 for p in sp}
    for p in sp:roots[p].add(p);genes[p].update(S[p]["genes"]);paths[p]=[(p,)]
    edges=[];seen=set()
    def edge(a,b,g):
        if a==b or(a,b,g) in seen:return
        seen.add((a,b,g));roots[b].update(roots[a]);genes[b].update(genes[a]);parents[b].add(a);ming[b]=min(ming.get(b,g),g)
        for q in paths[a]:
            z=q+(b,)
            if b not in q and z not in paths[b] and len(paths[b])<MAXPATH:paths[b].append(z)
        edges.append({"generation":g,"parent_pmid":a,"child_pmid":b})
    I=icite(sp)
    for a in sp:
        for b in I.get(a,{}).get("refs",[]):edge(a,b,1)
    G1={e["child_pmid"] for e in edges if e["generation"]==1};print("G1",len(G1),flush=True)
    I.update(icite(G1));P=pubmed(G1|sp);N={}
    def sn(xs):
        for p in xs:
            r=merge(p,I,P);q,gr,tm,gs,st,vs,gm=score(r,genes[p],gps)
            r.update({"min_generation":ming.get(p,0 if p in sp else ""),"is_seed":p in sp,
            "seed_confidence":S.get(p,{}).get("confidence",""),"seed_genes":S.get(p,{}).get("genes",[]),
            "root_seed_count":len(roots[p]),"direct_seed_count":len({x for x in parents[p] if x in sp}),
            "parent_count":len(parents[p]),"root_seed_pmids":sorted(roots[p],key=int),"ancestor_genes":sorted(genes[p]),
            "gene_matches":gs,"matched_terms":tm,"relevance_score":q,"relevance_grade":gr,
            "strong_temperature_signal":st,"variant_signal":vs,"ancestor_gene_match":gm,
            "example_path":" → ".join(paths[p][0]) if paths[p] else ""});N[p]=r
    sn(G1|sp)
    X1=choose(G1-sp,N,roots,parents,G1_THR,G1_CAP,3);print("expand G1",len(X1),flush=True)
    for a in X1:
        for b in I.get(a,{}).get("refs",[]):edge(a,b,2)
    G2={e["child_pmid"] for e in edges if e["generation"]==2};new=G2-set(I);I.update(icite(new));P.update(pubmed(G2-set(P)));sn(G2)
    X2=choose(G2-sp,N,roots,parents,G2_THR,G2_CAP,2);print("expand G2",len(X2),flush=True)
    for a in X2:
        for b in I.get(a,{}).get("refs",[]):edge(a,b,3)
    G3={e["child_pmid"] for e in edges if e["generation"]==3};new=G3-set(I);I.update(icite(new));P.update(pubmed(G3-set(P)));sn(G3);sn(set(N))
    for p,r in N.items():
        r["expanded_as_parent"]=p in set(X1)|set(X2);r["expanded_generation"]=2 if p in X1 else 3 if p in X2 else ""
        r["pubmed_url"]=f"https://pubmed.ncbi.nlm.nih.gov/{p}/";r["doi_url"]=f"https://doi.org/{r['doi']}" if r["doi"] else ""
    for e in edges:
        a=N.get(e["parent_pmid"],{});b=N.get(e["child_pmid"],{})
        e.update({"parent_title":a.get("title",""),"child_title":b.get("title",""),"child_year":b.get("year",""),
        "child_journal":b.get("journal",""),"child_relevance_grade":b.get("relevance_grade",""),
        "child_relevance_score":b.get("relevance_score",""),"child_gene_matches":b.get("gene_matches",[]),
        "child_matched_terms":b.get("matched_terms",[]),"root_seed_count":len(roots[e["child_pmid"]]),
        "root_seed_pmids":sorted(roots[e["child_pmid"]],key=int),"root_genes":sorted(genes[e["child_pmid"]]),
        "parent_is_seed":e["parent_pmid"] in sp,"child_is_seed":e["child_pmid"] in sp,
        "child_pubmed_url":f"https://pubmed.ncbi.nlm.nih.gov/{e['child_pmid']}/"})
    rows=sorted(N.values(),key=lambda r:(r["min_generation"],-r["relevance_score"],int(r["pmid"])))
    pr=[r for r in rows if not r["is_seed"] and r["relevance_score"]>=PRIO]
    nf=["pmid","min_generation","is_seed","relevance_grade","relevance_score","expanded_as_parent","expanded_generation",
    "title","abstract","year","journal","authors","doi","pmcid","citation_count","relative_citation_ratio","nih_percentile",
    "is_research_article","reference_count","root_seed_count","direct_seed_count","parent_count","ancestor_genes","gene_matches",
    "matched_terms","strong_temperature_signal","variant_signal","ancestor_gene_match","root_seed_pmids","example_path","mesh",
    "icite_found","pubmed_found","pubmed_url","doi_url"]
    ef=["generation","parent_pmid","parent_title","child_pmid","child_title","child_year","child_journal",
    "child_relevance_grade","child_relevance_score","child_gene_matches","child_matched_terms","root_seed_count",
    "root_seed_pmids","root_genes","parent_is_seed","child_is_seed","child_pubmed_url"]
    wcsv(OUT/"all_nodes.csv",rows,nf);wcsv(OUT/"prioritised_candidates.csv",pr,nf);wcsv(OUT/"citation_edges.csv",edges,ef)
    pa=[]
    for r in pr:
        for j,z in enumerate(paths[r["pmid"]],1):
            s=z[0];pa.append({"candidate_pmid":r["pmid"],"candidate_title":r["title"],"min_generation":r["min_generation"],
            "relevance_grade":r["relevance_grade"],"relevance_score":r["relevance_score"],"seed_pmid":s,
            "seed_title":N.get(s,{}).get("title",""),"seed_genes":S[s]["genes"],"path_number":j,"path_pmids":z,
            "path_length":len(z)-1,"candidate_pubmed_url":r["pubmed_url"]})
    pf=["candidate_pmid","candidate_title","min_generation","relevance_grade","relevance_score","seed_pmid","seed_title",
    "seed_genes","path_number","path_pmids","path_length","candidate_pubmed_url"];wcsv(OUT/"example_paths.csv",pa,pf)
    sm=[]
    for s in sorted(sp,key=int):
        reach={e["child_pmid"] for e in edges if s in roots[e["child_pmid"]]}
        sm.append({"seed_pmid":s,"seed_title":N[s]["title"],"seed_year":N[s]["year"],"seed_journal":N[s]["journal"],
        "seed_confidence":S[s]["confidence"],"seed_genes":S[s]["genes"],
        "g1_reference_edges":sum(1 for e in edges if e["generation"]==1 and e["parent_pmid"]==s),
        "reachable_unique_nodes":len(reach),"high_relevance_nodes":sum(N[p]["relevance_grade"]=="High" for p in reach if p in N),
        "medium_relevance_nodes":sum(N[p]["relevance_grade"]=="Medium" for p in reach if p in N),
        "g1_expansion_parents":sum(s in roots[p] for p in X1),"g2_expansion_parents":sum(s in roots[p] for p in X2),
        "pubmed_url":N[s]["pubmed_url"]})
    sf=["seed_pmid","seed_title","seed_year","seed_journal","seed_confidence","seed_genes","g1_reference_edges",
    "reachable_unique_nodes","high_relevance_nodes","medium_relevance_nodes","g1_expansion_parents","g2_expansion_parents","pubmed_url"]
    wcsv(OUT/"seed_coverage.csv",sm,sf)
    summary={"generated_at_utc":datetime.now(timezone.utc).isoformat(),"seed_count":len(sp),
    "generation_1":{"scope":"Complete iCite/OCC backward references from all seeds","unique_nodes":len(G1),
    "edges":sum(e["generation"]==1 for e in edges),"parents_expanded_to_g2":len(X1),"threshold":G1_THR,"cap":G1_CAP},
    "generation_2":{"scope":"Complete references from relevance-ranked G1 parents","unique_nodes":len(G2),
    "edges":sum(e["generation"]==2 for e in edges),"parents_expanded_to_g3":len(X2),"threshold":G2_THR,"cap":G2_CAP},
    "generation_3":{"scope":"Complete references from relevance-ranked G2 parents; no further expansion","unique_nodes":len(G3),
    "edges":sum(e["generation"]==3 for e in edges)},"network":{"unique_nodes_including_seeds":len(N),
    "citation_edges":len(edges),"prioritised_nonseed_candidates":len(pr),"high_relevance_candidates":sum(r["relevance_grade"]=="High" for r in pr),
    "medium_relevance_candidates":sum(r["relevance_grade"]=="Medium" for r in pr),"example_paths":len(pa)},
    "note":"Relevance scores are discovery triage, not final inclusion decisions.",
    "sources":{"icite":"https://icite.od.nih.gov/api/pubs","pubmed":"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"}}
    (OUT/"summary.json").write_text(json.dumps(summary,indent=2))
    (OUT/"METHOD.md").write_text("# CalorType set 2 citation snowball\n\nBackward references were retrieved from NIH iCite/OCC. Generation 1 is complete for all 240 seeds; generations 2 and 3 expand complete reference lists from relevance-ranked parents. Scores combine temperature/heat/cold/fever language, variant language, exact path-gene matches, phenotype terms and multi-seed convergence. Human review is required.\n")
    print(json.dumps(summary,indent=2),flush=True)
if __name__=="__main__":main()
