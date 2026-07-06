python3 make_q35_case35_generator.py fetched_sources/kattis_19903326_fetched.cpp q35_case35_submit.cpp
# or, from the repo root, omit args and it will try known current-best source names.

g++ -std=c++17 -O2 -pipe -static -s q35_case35_submit.cpp -o q35_case35_submit
# Submit q35_case35_submit.cpp