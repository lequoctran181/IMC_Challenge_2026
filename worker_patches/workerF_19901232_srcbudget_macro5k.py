#!/usr/bin/env python3
import re
import sys
from pathlib import Path

if len(sys.argv) != 3:
    raise SystemExit("usage: workerF_19901232_srcbudget_macro5k.py input.cpp output.cpp")

src_path = Path(sys.argv[1])
out_path = Path(sys.argv[2])
s = src_path.read_text()

main_anchor = (
    'int main(){JC();GN();if(!W2G::run())W2C::run();W5::post_patch_pass();'
    'VIMP::run();MIDEC::run();WK::run();B16::R(39000,60000,220,-7,192,.96,18.05);'
    'if(0&&"B16P515A"){}B16::R(39000,60000,76,-10,192,.96,18.35);'
    'for(int i=2;i--&&N>47500&N<6e4&&es()<18.6;)WK::run();'
    'if(N<50625&&es()<18.9)WK::run();JD();}'
)
if s.count(main_anchor) != 1:
    raise SystemExit("abort: exact 19901232/543 main anchor not found once")

for needle in [
    "namespace B16",
    "namespace W2G",
    "namespace W2C",
    "namespace W5",
    "namespace VIMP",
    "namespace MIDEC",
    "namespace WK",
]:
    if s.count(needle) != 1:
        raise SystemExit(f"abort: anchor {needle!r} count is {s.count(needle)}")

for stale in ["W5L473", "76,-13", "HH::", "FANX"]:
    if stale in s:
        raise SystemExit(f"abort: stale/drop marker present: {stale}")

macro_map = [
    ("O0", "double"),
    ("O1", "static"),
    ("O2", "const"),
    ("O3", "inline"),
    ("O4", "vector"),
    ("O5", "unsigned"),
    ("O6", "namespace"),
    ("O7", "struct"),
    ("O8", "template"),
    ("O9", "typename"),
    ("P0", "return false"),
    ("P1", "return true"),
    ("P2", "continue"),
    ("P3", "return"),
    ("P4", "int"),
    ("P5", "bool"),
    ("P6", "void"),
    ("P7", "class"),
]

for name, _ in macro_map:
    if re.search(rf"\b{re.escape(name)}\b", s):
        raise SystemExit(f"abort: macro name already present as token: {name}")

macro_block = "".join(f"#define {name} {value}\n" for name, value in macro_map)
pos = s.find("using namespace std;")
if pos < 0:
    raise SystemExit("abort: using namespace std; not found")
s = s[:pos] + macro_block + s[pos:]

single = {
    "double": "O0",
    "static": "O1",
    "const": "O2",
    "inline": "O3",
    "vector": "O4",
    "unsigned": "O5",
    "namespace": "O6",
    "struct": "O7",
    "template": "O8",
    "typename": "O9",
    "continue": "P2",
    "return": "P3",
    "int": "P4",
    "bool": "P5",
    "void": "P6",
    "class": "P7",
}
single_re = re.compile(r"\b(" + "|".join(re.escape(k) for k in sorted(single, key=len, reverse=True)) + r")\b")

def rewrite_code_piece(piece: str) -> str:
    piece = re.sub(r"\breturn\s+false\b", "P0", piece)
    piece = re.sub(r"\breturn\s+true\b", "P1", piece)
    return single_re.sub(lambda m: single[m.group(1)], piece)

def rewrite_line(line: str) -> str:
    if line.lstrip().startswith("#"):
        return line
    out = []
    code = []
    i = 0
    n = len(line)
    while i < n:
        ch = line[i]
        if ch in ("'", '"'):
            if code:
                out.append(rewrite_code_piece("".join(code)))
                code.clear()
            q = ch
            j = i + 1
            esc = False
            while j < n:
                c = line[j]
                if esc:
                    esc = False
                elif c == "\\":
                    esc = True
                elif c == q:
                    j += 1
                    break
                j += 1
            out.append(line[i:j])
            i = j
        elif ch == "/" and i + 1 < n and line[i + 1] == "/":
            if code:
                out.append(rewrite_code_piece("".join(code)))
                code.clear()
            out.append(line[i:])
            break
        else:
            code.append(ch)
            i += 1
    if code:
        out.append(rewrite_code_piece("".join(code)))
    return "".join(out)

out = "".join(rewrite_line(line) for line in s.splitlines(True))
saved = len(src_path.read_bytes()) - len(out.encode())
if saved < 5000:
    raise SystemExit(f"abort: saved only {saved} bytes")
if len(out.encode()) > 131052:
    raise SystemExit(f"abort: output too large: {len(out.encode())}")

out_path.write_text(out)
print(f"input bytes: {len(src_path.read_bytes())}")
print(f"output bytes: {len(out.encode())}")
print(f"saved bytes: {saved}")
