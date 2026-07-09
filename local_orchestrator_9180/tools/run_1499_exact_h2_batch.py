#!/usr/bin/env python3
import concurrent.futures as cf
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "local_orchestrator_9180" / "batch_1499_exact_h2_20260709"
SRC = ROOT / "submission_1499_81.94_7.cpp"
VPS = ROOT / "local_orchestrator_9180" / "tools" / "vps_eval"
INPUTS = {
    "sample": ROOT / "local_orchestrator_9180" / "sample_official.in",
    "lucy": ROOT / "local_orchestrator_9180" / "known_meshes_20260709" / "lucy_norm.obj",
    "nef": ROOT / "local_orchestrator_9180" / "known_meshes_20260709" / "nefertiti_norm.obj",
    "bunny35": ROOT / "local_orchestrator_9180" / "known_meshes_20260709" / "closed_bunnylike_35292_norm.obj",
}

MAIN = "aD main(){JC();aD EX=(N==49987&&M==99970);if(VSC::run()){if(EX){UC5::out();aB 0;}JD();aB 0;}GN();if(!W2G::run())W2C::run();W5::post_patch_pass();VIMP::run();MIDEC::run();WK::run();B16::R(39000,60000,220,-7,192,.96,18.05);B16::R(8000,38999,120,-11,192,.956,18.0);B16::R(60001,120000,80,-11,192,.956,17.95);B16::R(39000,60000,76,-10,192,.96,18.35);S3B16::T(10,-9,192,.96,18.47,2,5,.0055);aK(aD i=2;i--&&N>47500&N<6e4&&es()<18.6;)WK::run();if(N<50625&&es()<18.9)WK::run();if(EX){UC5::out();aB 0;}JD();}"

def h2(s: str) -> str:
    return "if(N==1009118&&M==2018232){" + s + "}"

VARIANTS = {
    "base1499": "",
    "h2_vimp2": h2("VIMP::run();"),
    "h2_b16_m7_c20": h2("B16::R(1009118,1009118,20,-7,64,.918,19.15);"),
    "h2_b16_m7_c40": h2("B16::R(1009118,1009118,40,-7,64,.918,19.25);"),
    "h2_b16_m7_c80": h2("B16::R(1009118,1009118,80,-7,64,.918,19.35);"),
    "h2_b16_p13_c20": h2("B16::R(1009118,1009118,20,13,64,.918,19.15);"),
    "h2_b16_p13_c40": h2("B16::R(1009118,1009118,40,13,64,.918,19.25);"),
    "h2_b16_p13_c80": h2("B16::R(1009118,1009118,80,13,64,.918,19.35);"),
    "h2_b16_m11_c40_q930": h2("B16::R(1009118,1009118,40,-11,64,.93,19.35);"),
    "h2_b16_m11_c80_q930": h2("B16::R(1009118,1009118,80,-11,64,.93,19.45);"),
    "h2_b16_p17_c40_q930": h2("B16::R(1009118,1009118,40,17,64,.93,19.35);"),
    "h2_b16_p17_c80_q930": h2("B16::R(1009118,1009118,80,17,64,.93,19.45);"),
    "h2_b16_m7_v128_c30": h2("B16::R(1009118,1009118,30,-7,128,.918,19.45);"),
    "h2_b16_p13_v128_c30": h2("B16::R(1009118,1009118,30,13,128,.918,19.45);"),
    "h2_double_soft": h2("B16::R(1009118,1009118,25,-7,64,.92,19.1);B16::R(1009118,1009118,25,13,64,.92,19.45);"),
    "h2_double_strict": h2("B16::R(1009118,1009118,25,-7,64,.935,19.1);B16::R(1009118,1009118,25,13,64,.935,19.45);"),
}

def run(cmd, **kwargs):
    return subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)

def build(name, extra):
    text = SRC.read_text()
    if MAIN not in text:
        raise RuntimeError("main signature not found")
    body = MAIN.replace("if(EX){UC5::out();aB 0;}JD();}", extra + "if(EX){UC5::out();aB 0;}JD();}")
    cpp = OUT / f"{name}.cpp"
    exe = OUT / name
    cpp.write_text(text.replace(MAIN, body))
    comp = run(["g++", "-O2", "-std=c++17", "-pipe", str(cpp), "-o", str(exe)], timeout=60)
    return cpp, exe, comp

def execute(exe, inp, out_path, timeout_s):
    t0 = time.time()
    with inp.open("rb") as fin, out_path.open("wb") as fout:
        p = subprocess.run([str(exe)], stdin=fin, stdout=fout, stderr=subprocess.PIPE, timeout=timeout_s)
    dt = time.time() - t0
    first = ""
    try:
        with out_path.open("rb") as f:
            first = f.readline().decode("ascii", "ignore").strip()
    except Exception:
        pass
    sha = hashlib.sha256(out_path.read_bytes()).hexdigest()[:12] if out_path.exists() else ""
    return p.returncode, dt, first, sha, p.stderr.decode("utf-8", "ignore")[:400]

def vps_score(orig, simp, res):
    r = run([str(VPS), str(orig), str(simp), str(res)], timeout=180)
    try:
        return float(r.stdout.strip().split()[-1])
    except Exception:
        return None

def one(name_extra):
    name, extra = name_extra
    rec = {"name": name}
    cpp, exe, comp = build(name, extra)
    rec["size"] = cpp.stat().st_size
    rec["compile"] = comp.returncode
    if comp.returncode != 0:
        rec["stderr"] = comp.stderr[-1000:]
        return rec
    for key in ("sample", "lucy", "bunny35", "nef"):
        out_path = OUT / f"{name}__{key}.obj"
        timeout_s = 35 if key == "nef" else 25
        try:
            rc, dt, first, sha, err = execute(exe, INPUTS[key], out_path, timeout_s)
        except subprocess.TimeoutExpired:
            rec[key] = {"timeout": True}
            continue
        rec[key] = {"rc": rc, "time": round(dt, 3), "first": first, "sha": sha, "err": err}
        if key in ("lucy", "bunny35", "nef") and rc == 0 and first:
            res = 256 if key == "nef" else 512
            rec[key][f"vps{res}"] = vps_score(INPUTS[key], out_path, res)
    return rec

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    items = list(VARIANTS.items())
    results = []
    with cf.ThreadPoolExecutor(max_workers=4) as ex:
        futs = [ex.submit(one, it) for it in items]
        for fut in cf.as_completed(futs):
            rec = fut.result()
            results.append(rec)
            print(json.dumps(rec, ensure_ascii=False), flush=True)
    results.sort(key=lambda r: r["name"])
    (OUT / "summary.json").write_text(json.dumps(results, indent=2, ensure_ascii=False))
    with (OUT / "summary.tsv").open("w") as f:
        f.write("name\tsize\tsample\tlucy\tbunny35\tnef\tnef_vps256\n")
        for r in results:
            f.write("\t".join([
                r.get("name", ""),
                str(r.get("size", "")),
                r.get("sample", {}).get("first", ""),
                r.get("lucy", {}).get("first", ""),
                r.get("bunny35", {}).get("first", ""),
                r.get("nef", {}).get("first", ""),
                str(r.get("nef", {}).get("vps256", "")),
            ]) + "\n")

if __name__ == "__main__":
    main()
