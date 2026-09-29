"""
ONE-SHOT PIPELINE (steps 1-4). Run:  python pipeline.py
 1. Look up real GRCh38 coordinates for a panel of MVA / chromosomal-instability /
    growth / tumour / kidney genes.
 2. Pull this patient's PASS, non-reference variants in those genes from the VCF.
 3. Annotate them with Ensembl VEP (consequence + gnomAD frequency). Cached, so
    re-running is fast and safe if the internet drops.
 4. Rank (compound-het pairs, homozygous, single het LoF) and write:
      ranked_candidates.csv   - everything, for you to inspect
      Yobrrrr_vep_filter.csv  - <=10 rows, the Track 1 submission
      Yobrrrr_track1_report.md
Needs:  pip install requests
"""
import gzip, json, os, sys, time, itertools, csv, requests

VCF = r"C:\Users\admin vit\Desktop\Dataset\WGS_EX2312012_HGWCNDSX7.vcf.gz"  # <-- edit if needed
PAD = 2000
MIN_GQ, MIN_DP = 30, 8
MAX_AF = 0.01
CACHE = "vep_cache.json"

# gene -> prior weight (how strongly the gene fits MVA / this phenotype)
PANEL = {
 # MVA / mitotic-checkpoint / centrosome
 "BUB1B":1.0,"CEP57":1.0,"TRIP13":1.0,"CDC20":0.9,"MAD2L1":0.9,"PLK4":0.8,"CENPE":0.7,
 # microcephaly / primordial dwarfism / growth restriction
 "PCNT":0.7,"ORC1":0.7,"ORC4":0.7,"ORC6":0.7,"CDT1":0.7,"CDC6":0.7,"CDC45":0.7,"GMNN":0.7,
 "LIG4":0.7,"NHEJ1":0.7,"XRCC4":0.7,"DNA2":0.7,"CENPJ":0.7,"CEP152":0.7,"CEP63":0.7,
 "CDK5RAP2":0.7,"MCPH1":0.7,"ASPM":0.7,"RNU4ATAC":0.6,"IGF1R":0.6,"IGF2":0.6,
 # DNA repair / chromosome instability / cancer predisposition
 "BLM":0.8,"WRN":0.6,"RECQL4":0.8,"NBN":0.7,"MRE11":0.7,"RAD50":0.7,"ATM":0.7,"ATR":0.7,
 "BRCA2":0.7,"PALB2":0.7,"FANCA":0.6,"FANCD2":0.6,"TP53":0.8,"DICER1":0.8,"NF1":0.7,
 "HRAS":0.7,"CDKN1C":0.7,"ESCO2":0.8,"DDX11":0.7,"NIPBL":0.6,"SMC1A":0.6,"SMC3":0.6,"RAD21":0.6,
 # nephrocalcinosis
 "CLDN16":0.6,"CLDN19":0.6,"CYP24A1":0.6,"SLC34A1":0.6,"SLC34A3":0.6,"ATP6V1B1":0.6,
 "ATP6V0A4":0.6,"SLC4A1":0.6,"CASR":0.6,"KCNJ1":0.6,"SLC12A1":0.6,"BSND":0.6,"CLCNKB":0.6,"HNF1B":0.6,
}
IMPACT_RANK = {"HIGH":3,"MODERATE":2,"LOW":1,"MODIFIER":0}
W = {"HIGH":1.0,"MODERATE":0.5}

# ---------------- 1. gene coordinates ----------------
def _ensembl(sym):
    url=f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{sym}?content-type=application/json"
    for _ in range(3):
        try:
            r=requests.get(url,timeout=30)
            if r.status_code==200:
                j=r.json(); return str(j["seq_region_name"]),j["start"],j["end"]
        except Exception: pass
        time.sleep(2)
def _mygene(sym):
    try:
        r=requests.get(f"https://mygene.info/v3/query?q=symbol:{sym}&species=human&fields=symbol,genomic_pos",timeout=30)
        for h in r.json().get("hits",[]):
            if h.get("symbol")!=sym: continue
            gp=h.get("genomic_pos"); gp=gp if isinstance(gp,list) else [gp]
            for p in gp:
                if p and str(p["chr"]) in [str(i) for i in range(1,23)]+["X"]:
                    return str(p["chr"]),p["start"],p["end"]
    except Exception: pass
def get_genes():
    g={}
    for s in PANEL:
        res=_ensembl(s) or _mygene(s)
        if res is None: print("  !! could not locate",s,"- skipping"); continue
        g[s]=res
    print(f"Located {len(g)}/{len(PANEL)} genes")
    return g

# ---------------- 2. extract ----------------
def extract(genes):
    by_chrom={}
    for s,(c,a,b) in genes.items(): by_chrom.setdefault(c,[]).append((a-PAD,b+PAD,s))
    out=[]
    with gzip.open(VCF,"rt") as f:
        for line in f:
            if line[0]=="#": continue
            t=line.split("\t",9)
            if t[0] not in by_chrom or t[6]!="PASS": continue
            pos=int(t[1])
            for a,b,s in by_chrom[t[0]]:
                if a<=pos<=b:
                    fmt=dict(zip(t[8].split(":"),t[9].split()[0].split(":")))
                    gt=fmt.get("GT","./.")
                    if gt in("0/0","./.","."): continue
                    try: gq=int(fmt.get("GQ",0)); dp=int(fmt.get("DP",0))
                    except ValueError: continue
                    if gq<MIN_GQ or dp<MIN_DP: continue
                    out.append(dict(gene=s,chrom=t[0],pos=pos,ref=t[3],alt=t[4],gt=gt,gq=gq,dp=dp,ad=fmt.get("AD",""),rsid=t[2]))
    print(f"Extracted {len(out)} PASS variants (GQ>={MIN_GQ}, DP>={MIN_DP}) in panel genes")
    return out

# ---------------- 3. VEP ----------------
def vkey(v): return f"{v['chrom']} {v['pos']} . {v['ref']} {v['alt']} . . ."
PARAM_SETS=["canonical=1&mane=1&hgvs=1&af_gnomadg=1&af_gnomade=1","canonical=1&hgvs=1&af=1","canonical=1"]
def vep_batch(batch):
    for ps in PARAM_SETS:
        url="https://rest.ensembl.org/vep/human/region?"+ps
        for _ in range(4):
            try:
                r=requests.post(url,headers={"Content-Type":"application/json","Accept":"application/json"},
                                data=json.dumps({"variants":batch}),timeout=120)
            except Exception: time.sleep(3); continue
            if r.status_code==200: return r.json()
            if r.status_code==429: time.sleep(float(r.headers.get("Retry-After",3))); continue
            if r.status_code==400: break            # bad option -> try simpler param set
            time.sleep(3)
    return None
def annotate(vars_):
    cache=json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    todo=[vkey(v) for v in vars_ if vkey(v) not in cache]
    print(f"VEP: {len(cache)} cached, {len(todo)} to annotate")
    for i in range(0,len(todo),100):
        b=todo[i:i+100]; res=vep_batch(b)
        if res is None:
            print("  batch failed, single-variant fallback")
            res=[]
            for k in b:
                r1=vep_batch([k]); res+= (r1 or [])
                time.sleep(0.3)
        for r in res: cache[r["input"]]=r
        for k in b: cache.setdefault(k,{"input":k,"_missing":True})
        json.dump(cache,open(CACHE,"w")); print(f"  annotated {min(i+100,len(todo))}/{len(todo)}"); time.sleep(0.5)
    return cache
def parse(res,gene):
    best=None
    for tc in res.get("transcript_consequences",[]):
        if tc.get("gene_symbol")!=gene: continue
        key=(IMPACT_RANK.get(tc.get("impact"),0),1 if tc.get("mane_select") else 0,tc.get("canonical",0) or 0)
        if best is None or key>best[0]: best=(key,tc)
    af=0.0
    for cv in res.get("colocated_variants",[]) or []:
        for fr in (cv.get("frequencies") or {}).values():
            for x in fr.values():
                if isinstance(x,(int,float)): af=max(af,x)
    if best is None: return "none",",".join([res.get("most_severe_consequence","?")]),"",af
    tc=best[1]
    return tc.get("impact","MODIFIER"),",".join(tc.get("consequence_terms",[])),tc.get("hgvsp") or tc.get("hgvsc") or "",af

# ---------------- 4. rank ----------------
def rarity(af): return 1.0 if af==0 else 1-min(af/MAX_AF,1)*0.5
def rank(vars_,cache):
    good={}
    for v in vars_:
        r=cache.get(vkey(v),{})
        imp,cons,hg,af=parse(r,v["gene"])
        v.update(impact=imp,cons=cons,hgvs=hg,af=af)
        if imp in W and af<MAX_AF: good.setdefault(v["gene"],[]).append(v)
    cands=[]
    for g,vs in good.items():
        prior=PANEL[g]
        vs=sorted(vs,key=lambda x:-W[x["impact"]]*rarity(x["af"]))[:8]
        for v in vs:
            w=W[v["impact"]]*rarity(v["af"])*prior
            if v["gt"] in("1/1","1|1"):
                cands.append(dict(kind="homozygous",score=min(w,0.99),v1=v,v2=None))
            else:
                cands.append(dict(kind="single_het",score=0.45*w,v1=v,v2=None))
        hets=[v for v in vs if v["gt"] not in("1/1","1|1")]
        for a,b in itertools.combinations(hets,2):
            if a["pos"]==b["pos"]: continue
            s=(W[a["impact"]]*rarity(a["af"])+W[b["impact"]]*rarity(b["af"]))/2*prior
            if "HIGH" in(a["impact"],b["impact"]): s=min(s+0.1,0.99)
            cands.append(dict(kind="compound_het",score=min(s,0.99),v1=min(a,b,key=lambda x:x["pos"]),v2=max(a,b,key=lambda x:x["pos"])))
    cands.sort(key=lambda c:-c["score"]); return cands
def desc(v): return f"{v['gene']} {v['chrom']}:{v['pos']} {v['ref']}>{v['alt']} {v['gt']} {v['impact']} {v['cons']} {v['hgvs']} AF={v['af']:.4g}"

def main():
    genes=get_genes(); vs=extract(genes)
    cache=annotate(vs); cands=rank(vs,cache)
    with open("ranked_candidates.csv","w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["rank","kind","score","variant_1","variant_2"])
        for i,c in enumerate(cands,1): w.writerow([i,c["kind"],round(c["score"],3),desc(c["v1"]),desc(c["v2"]) if c["v2"] else ""])
    print("\n===== TOP 20 CANDIDATES (paste this to Claude) =====")
    for i,c in enumerate(cands[:20],1):
        print(f"{i:2}. [{c['kind']}] score={c['score']:.2f}\n     {desc(c['v1'])}"+(f"\n     {desc(c['v2'])}" if c["v2"] else ""))
    if not cands: print("No qualifying candidates - paste this message to Claude."); return
    cols=["proband_id","chrom_1","pos_1","ref_1","alt_1","chrom_2","pos_2","ref_2","alt_2","epcr","finding_type","notes"]
    if os.path.exists("Yobrrrr_vep_filter.csv"): os.replace("Yobrrrr_vep_filter.csv","Yobrrrr_vep_filter_attempt1.csv")
    with open("Yobrrrr_vep_filter.csv","w",newline="") as fh:
        w=csv.writer(fh); w.writerow(cols)
        for c in cands[:10]:
            a,b=c["v1"],c["v2"]
            w.writerow(["PROBAND01",a["chrom"],a["pos"],a["ref"],a["alt"],
                        b["chrom"] if b else "",b["pos"] if b else "",b["ref"] if b else "",b["alt"] if b else "",
                        round(c["score"],3),"primary",f"{c['kind']} {a['gene']}: {a['cons']}"+(f" + {b['cons']}" if b else "")])
    top=cands[:10]
    with open("Yobrrrr_track1_report.md","w") as fh:
        fh.write("# Track 1 Methods Report\n\n**Proband:** PROBAND01 (MVA)  \n**Input:** WGS VCF (GRCh38), plus clinical phenotype document\n\n")
        fh.write("## Method\n1. Built a gene panel (MVA/mitotic checkpoint: BUB1B, CEP57, TRIP13, CDC20, MAD2L1, PLK4, CENPE; plus genes for microcephaly/primordial dwarfism, DNA-repair/chromosome-instability, rhabdomyosarcoma predisposition and nephrocalcinosis matching the phenotype). Gene coordinates were fetched live from Ensembl (GRCh38), not assumed.\n")
        fh.write(f"2. Extracted variants with FILTER=PASS, non-reference genotype, GQ>={MIN_GQ}, DP>={MIN_DP}, within gene bodies +/-{PAD} bp.\n")
        fh.write("3. Annotated with Ensembl VEP REST (consequence, impact, MANE/canonical transcript, gnomAD frequency).\n")
        fh.write(f"4. Kept HIGH/MODERATE impact variants with population AF < {MAX_AF}. Scored by impact x rarity x gene prior; considered compound-heterozygous pairs, homozygous variants, and single heterozygous variants (partial credit).\n\n## Top candidates\n")
        for i,c in enumerate(top,1):
            fh.write(f"{i}. **{c['kind']}** (score {c['score']:.2f}): {desc(c['v1'])}"+(f" ; {desc(c['v2'])}" if c["v2"] else "")+"\n")
        fh.write("\n## Limitations\nNo parental data (phase unknown); gene-panel approach may miss novel genes; no structural-variant/CNV analysis.\n\n## AI assistant disclosure\n- **Provider:** Anthropic\n- **Model:** Claude (via claude.ai)\n- **Data-handling terms:** [FILL IN from your Claude plan/settings before submitting]\n")
    print("\nWrote: ranked_candidates.csv, Yobrrrr_vep_filter.csv, Yobrrrr_track1_report.md")
if __name__=="__main__": main()
