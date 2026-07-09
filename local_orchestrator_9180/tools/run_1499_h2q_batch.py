#!/usr/bin/env python3
import concurrent.futures as cf
import hashlib
import json
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "local_orchestrator_9180" / "batch_1499_h2q_20260709"
SRC = ROOT / "submission_1499_81.94_7.cpp"
VPS = ROOT / "local_orchestrator_9180" / "tools" / "vps_eval"
NEF = ROOT / "local_orchestrator_9180" / "known_meshes_20260709" / "nefertiti_norm.obj"
SAMPLE = ROOT / "local_orchestrator_9180" / "sample_official.in"
LUCY = ROOT / "local_orchestrator_9180" / "known_meshes_20260709" / "lucy_norm.obj"
BUNNY = ROOT / "local_orchestrator_9180" / "known_meshes_20260709" / "closed_bunnylike_35292_norm.obj"

H2Q = r'''bW H2Q{aW E{aA l;aU bL bL k;};aE aJ run(aD tgt,aA need,aA ai,aA bb,aA ang,aA lim,aD vp){if(!(N==1009118&&M==2018232)||es()>18.4)aB aF;AP S=AD();aD SN=bI();if(SN<=tgt)aB aF;aG<aU bL bL>ks;ks.aN((aY)BE*3);aK(aD aR=0;aR<(aD)bA.aM();++aR){if(!BR[aR])aH;aC aQ&f=bA[aR];if(!BU[f.v[0]]||!BU[f.v[1]]||!BU[f.v[2]])aH;ks.pb(ED(f.v[0],f.v[1]));ks.pb(ED(f.v[1],f.v[2]));ks.pb(ED(f.v[2],f.v[0]));}dB(ks.aV(),ks.cR());ks.iP(gU(ks.aV(),ks.cR()),ks.cR());aG<E>ed;ed.aN(ks.aM());aK(aU bL bL k:ks){aD a=(aD)(k>>32),b=(aD)(k&0xffffffffu);if(a>=0&&a<N&&b>=0&&b<N&&BU[a]&&BU[b])ed.pb({bE(P[a]-P[b]),k});}dB(ed.aV(),ed.cR(),[](aC E&a,aC E&b){if(a.l!=b.l)aB a.l<b.l;aB a.k<b.k;});AE p;p.AI=ai*CL;p.BB=bb*CL;p.BQ=fH(ang*cP(-1.)/180.);p.W=bb*CL;p.AT=fH((ang*.7)*cP(-1.)/180.);p.AJ=aX(1e-30,1e-24*CL*CL);p.AA=aP;aD rm=0;for(aD i=0;i<(aD)ed.aM()&&bI()>tgt&&es()<lim;++i){aD a=(aD)(ed[i].k>>32),b=(aD)(ed[i].k&0xffffffffu);if(a>=0&&a<N&&b>=0&&b<N&&BU[a]&&BU[b])if(GD(a,b,p))++rm;}aD AN=bI();if(rm<100||AN>=SN||!W5::zaj()||es()>20.5||(vp&&eW(vp)<need)){rs(S);aB aF;}aB aP;}}'''

MAIN_MARK = "VIMP::run();MIDEC::run();"

VARIANTS = []
for tgt in (130000, 120000, 110000, 100000, 90000):
    for ai, bb, ang in ((.048, .018, 22.), (.052, .020, 25.), (.055, .022, 28.)):
        VARIANTS.append((f"h2q_t{tgt//1000}_a{int(ai*1000)}_p64", tgt, .90, ai, bb, ang, 20.05, 64))
VARIANTS.append(("h2q_t110_nopxy", 110000, .0, .055, .022, 28., 20.05, 0))
VARIANTS.append(("h2q_t100_nopxy", 100000, .0, .055, .022, 28., 20.05, 0))

def run(cmd, **kw):
    return subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kw)

def make_source(name, tgt, need, ai, bb, ang, lim, vp):
    src = SRC.read_text()
    if "bW H2Q{" not in src:
        src = src.replace("bW B16{", H2Q + "bW B16{")
    call = f"VIMP::run();H2Q::run({tgt},{need:.4g},{ai:.4g},{bb:.4g},{ang:.4g},{lim:.4g},{vp});MIDEC::run();"
    src = src.replace(MAIN_MARK, call)
    p = OUT / f"{name}.cpp"
    p.write_text(src)
    return p

def execute(exe, inp, out, timeout_s):
    t0 = time.time()
    with inp.open("rb") as fin, out.open("wb") as fout:
        p = subprocess.run([str(exe)], stdin=fin, stdout=fout, stderr=subprocess.PIPE, timeout=timeout_s)
    dt = time.time() - t0
    first = out.open("rb").readline().decode("ascii", "ignore").strip()
    sha = hashlib.sha256(out.read_bytes()).hexdigest()[:12]
    return {"rc": p.returncode, "time": round(dt, 3), "first": first, "sha": sha, "err": p.stderr.decode("utf-8", "ignore")[:200]}

def score(orig, simp, res):
    r = run([str(VPS), str(orig), str(simp), str(res)], timeout=180)
    try:
        return float(r.stdout.strip().replace("\\n", ""))
    except Exception:
        return None

def one(v):
    name, tgt, need, ai, bb, ang, lim, vp = v
    cpp = make_source(name, tgt, need, ai, bb, ang, lim, vp)
    exe = OUT / name
    comp = run(["g++", "-O2", "-std=c++17", "-pipe", str(cpp), "-o", str(exe)], timeout=60)
    rec = {"name": name, "size": cpp.stat().st_size, "compile": comp.returncode}
    if comp.returncode:
        rec["stderr"] = comp.stderr[-1000:]
        return rec
    nefout = OUT / f"{name}__nef.obj"
    try:
        rec["nef"] = execute(exe, NEF, nefout, 35)
        rec["nef"]["vps256"] = score(NEF, nefout, 256)
        rec["nef"]["vps512"] = score(NEF, nefout, 512)
    except subprocess.TimeoutExpired:
        rec["nef"] = {"timeout": True}
    return rec

def smoke(name):
    exe = OUT / name
    rec = {}
    for key, inp in (("sample", SAMPLE), ("lucy", LUCY), ("bunny35", BUNNY)):
        out = OUT / f"{name}__{key}.obj"
        rec[key] = execute(exe, inp, out, 28)
    return rec

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    results = []
    with cf.ThreadPoolExecutor(max_workers=4) as ex:
        futs = [ex.submit(one, v) for v in VARIANTS]
        for fut in cf.as_completed(futs):
            rec = fut.result()
            results.append(rec)
            print(json.dumps(rec, ensure_ascii=False), flush=True)
    results.sort(key=lambda r: (r.get("nef", {}).get("first", "999999"), r["name"]))
    if results:
        best = results[0]["name"]
        results[0]["smoke"] = smoke(best)
    (OUT / "summary.json").write_text(json.dumps(results, indent=2, ensure_ascii=False))
    with (OUT / "summary.tsv").open("w") as f:
        f.write("name\tsize\tnef\tvps256\tvps512\ttime\n")
        for r in results:
            n = r.get("nef", {})
            f.write("\t".join([r["name"], str(r.get("size", "")), n.get("first", ""), str(n.get("vps256", "")), str(n.get("vps512", "")), str(n.get("time", ""))]) + "\n")

if __name__ == "__main__":
    main()
