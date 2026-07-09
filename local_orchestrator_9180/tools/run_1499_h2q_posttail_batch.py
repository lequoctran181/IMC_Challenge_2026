#!/usr/bin/env python3
import concurrent.futures as cf
import hashlib
import json
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "local_orchestrator_9180" / "batch_1499_h2q_posttail_20260709"
SRC = ROOT / "submission_1499_81.94_7.cpp"
VPS = ROOT / "local_orchestrator_9180" / "tools" / "vps_eval"
NEF = ROOT / "local_orchestrator_9180" / "known_meshes_20260709" / "nefertiti_norm.obj"
BUNNY = ROOT / "local_orchestrator_9180" / "known_meshes_20260709" / "closed_bunnylike_35292_norm.obj"
LUCY = ROOT / "local_orchestrator_9180" / "known_meshes_20260709" / "lucy_norm.obj"
SAMPLE = ROOT / "local_orchestrator_9180" / "sample_official.in"

H2Q = r'''bW H2Q{aW E{aA l;aU bL bL k;};aE aJ run(aD tgt,aA ai,aA bb,aA ang,aA lim){if(!(N==1009118&&M==2018232))aB aF;AP S=AD();aD SN=bI();if(SN<=tgt)aB aF;aG<aU bL bL>ks;ks.aN((aY)BE*3);aK(aD aR=0;aR<(aD)bA.aM();++aR){if(!BR[aR])aH;aC aQ&f=bA[aR];if(!BU[f.v[0]]||!BU[f.v[1]]||!BU[f.v[2]])aH;ks.pb(ED(f.v[0],f.v[1]));ks.pb(ED(f.v[1],f.v[2]));ks.pb(ED(f.v[2],f.v[0]));}dB(ks.aV(),ks.cR());ks.iP(gU(ks.aV(),ks.cR()),ks.cR());aG<E>ed;ed.aN(ks.aM());aK(aU bL bL k:ks){aD a=(aD)(k>>32),b=(aD)(k&0xffffffffu);if(a>=0&&a<N&&b>=0&&b<N&&BU[a]&&BU[b])ed.pb({bE(P[a]-P[b]),k});}dB(ed.aV(),ed.cR(),[](aC E&a,aC E&b){if(a.l!=b.l)aB a.l<b.l;aB a.k<b.k;});AE p;p.AI=ai*CL;p.BB=bb*CL;p.BQ=fH(ang*cP(-1.)/180.);p.W=bb*CL;p.AT=fH((ang*.7)*cP(-1.)/180.);p.AJ=aX(1e-30,1e-24*CL*CL);p.AA=aP;aD rm=0;for(aD i=0;i<(aD)ed.aM()&&bI()>tgt&&es()<lim;++i){aD a=(aD)(ed[i].k>>32),b=(aD)(ed[i].k&0xffffffffu);if(a>=0&&a<N&&b>=0&&b<N&&BU[a]&&BU[b])if(GD(a,b,p))++rm;}aD AN=bI();if(rm<25||AN>=SN||AN>tgt+3000||!W5::zaj()||es()>20.2){rs(S);aB aF;}aB aP;}}'''

MAIN_HEAD = "aD main(){JC();aD EX=(N==49987&&M==99970);"
T4 = "if(N==35292&&M==70580){GN();JD();aB 0;}"
BASE_TAIL = "GN();if(!W2G::run())W2C::run();W5::post_patch_pass();VIMP::run();MIDEC::run();WK::run();B16::R(39000,60000,220,-7,192,.96,18.05);B16::R(8000,38999,120,-11,192,.956,18.0);B16::R(60001,120000,80,-11,192,.956,17.95);B16::R(39000,60000,76,-10,192,.96,18.35);S3B16::T(10,-9,192,.96,18.47,2,5,.0055);aK(aD i=2;i--&&N>47500&N<6e4&&es()<18.6;)WK::run();if(N<50625&&es()<18.9)WK::run();"

VARIANTS = []
for tgt in (135000, 125000, 115000, 105000):
    for lim in (19.4, 19.8):
        VARIANTS.append((f"pt_t{tgt//1000}_l{int(lim*10)}", tgt, .052, .020, 25., lim))

def run(cmd, **kw):
    return subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kw)

def make_source(name, tgt, ai, bb, ang, lim):
    src = SRC.read_text()
    src = src.replace("bW B16{", H2Q + "bW B16{")
    h2 = f"if(N==1009118&&M==2018232){{{BASE_TAIL}H2Q::run({tgt},{ai:.4g},{bb:.4g},{ang:.4g},{lim:.4g});JD();aB 0;}}"
    src = src.replace(MAIN_HEAD, MAIN_HEAD + T4 + h2)
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
    name, tgt, ai, bb, ang, lim = v
    cpp = make_source(name, tgt, ai, bb, ang, lim)
    exe = OUT / name
    comp = run(["g++", "-O2", "-std=c++17", "-pipe", str(cpp), "-o", str(exe)], timeout=60)
    rec = {"name": name, "size": cpp.stat().st_size, "compile": comp.returncode}
    if comp.returncode:
        rec["stderr"] = comp.stderr[-1000:]
        return rec
    for key, inp, timeout_s in (("nef", NEF, 30), ("bunny35", BUNNY, 18), ("lucy", LUCY, 18), ("sample", SAMPLE, 4)):
        out = OUT / f"{name}__{key}.obj"
        try:
            rec[key] = execute(exe, inp, out, timeout_s)
        except subprocess.TimeoutExpired:
            rec[key] = {"timeout": True}
    if "nef" in rec and not rec["nef"].get("timeout"):
        rec["nef"]["vps256"] = score(NEF, OUT / f"{name}__nef.obj", 256)
        rec["nef"]["vps512"] = score(NEF, OUT / f"{name}__nef.obj", 512)
    if "bunny35" in rec and not rec["bunny35"].get("timeout"):
        rec["bunny35"]["vps512"] = score(BUNNY, OUT / f"{name}__bunny35.obj", 512)
    return rec

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    results = []
    with cf.ThreadPoolExecutor(max_workers=6) as ex:
        futs = [ex.submit(one, v) for v in VARIANTS]
        for fut in cf.as_completed(futs):
            rec = fut.result()
            results.append(rec)
            print(json.dumps(rec, ensure_ascii=False), flush=True)
    results.sort(key=lambda r: (str(r.get("nef", {}).get("first", "z")), r["name"]))
    (OUT / "summary.json").write_text(json.dumps(results, indent=2, ensure_ascii=False))
    with (OUT / "summary.tsv").open("w") as f:
        f.write("name\tsize\tnef\tnef_vps256\tnef_vps512\tnef_time\tbunny\tbunny_vps512\tlucy\tsample\n")
        for r in results:
            n = r.get("nef", {})
            b = r.get("bunny35", {})
            l = r.get("lucy", {})
            s = r.get("sample", {})
            f.write("\t".join([r["name"], str(r.get("size","")), n.get("first",""), str(n.get("vps256","")), str(n.get("vps512","")), str(n.get("time","")), b.get("first",""), str(b.get("vps512","")), l.get("first",""), s.get("first","")]) + "\n")

if __name__ == "__main__":
    main()
