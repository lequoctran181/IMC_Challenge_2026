#!/usr/bin/env python3
import sys, pathlib

if len(sys.argv) < 3:
    print("usage: python3 gen_workerL_b16_case5.py fetched_sources/19901232.cpp workerL_b16_case5.cpp", file=sys.stderr)
    sys.exit(2)

src_path = pathlib.Path(sys.argv[1])
out_path = pathlib.Path(sys.argv[2])
s = src_path.read_text()

r1_old = 'if(0&&"B16P515A"){}'
r1_new = 'B16::R(47500,50625,700,11,256,.984,18.45);'

r2_old = 'for(int i=2;i--&&N>47500&N<6e4&&es()<18.6;)WK::run();if(N<50625&&es()<18.9)WK::run();JD();}'
r2_new = 'for(int i=2;i--&&N>47500&N<6e4&&es()<18.6;)WK::run();if(N<50625&&es()<18.75)B16::R(47500,50625,900,13,256,.982,19.25);if(N<50625&&es()<18.9)WK::run();JD();}'

if s.count(r1_old) != 1:
    print("fail-closed: anchor r1 not found exactly once", file=sys.stderr)
    sys.exit(1)
if s.count(r2_old) != 1:
    print("fail-closed: anchor r2 not found exactly once", file=sys.stderr)
    sys.exit(1)

t = s.replace(r1_old, r1_new).replace(r2_old, r2_new)

if len(t.encode()) > 131072:
    print("fail-closed: patched source too large:", len(t.encode()), file=sys.stderr)
    sys.exit(1)
if "int main(){JC();GN();" not in t or "JD();}" not in t:
    print("fail-closed: main/tail sanity check failed", file=sys.stderr)
    sys.exit(1)

out_path.write_text(t)
print(len(t.encode()))