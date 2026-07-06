#!/usr/bin/env bash
set -euo pipefail
g++ -O2 -std=c++17 -pipe -static -s B92pack_idx_periodic.cpp -o B92pack_idx_periodic
g++ -O2 -std=c++17 -pipe -static -s B92pack_idx_periodic_safe.cpp -o B92pack_idx_periodic_safe
g++ -O2 -std=c++17 -pipe -static -s B92pack_idx_sphere.cpp -o B92pack_idx_sphere
g++ -O2 -std=c++17 -pipe -static -s B92pack_idx_sphere_safe.cpp -o B92pack_idx_sphere_safe
g++ -O2 -std=c++17 -pipe -static -s B92pack_star_strict.cpp -o B92pack_star_strict
g++ -O2 -std=c++17 -pipe -static -s B92pack_star_loose.cpp -o B92pack_star_loose
g++ -O2 -std=c++17 -pipe -static -s B92pack_qnet_dsu.cpp -o B92pack_qnet_dsu
g++ -O2 -std=c++17 -pipe -static -s B92pack_combo_grid_sphere.cpp -o B92pack_combo_grid_sphere
g++ -O2 -std=c++17 -pipe -static -s B92pack_combo_all.cpp -o B92pack_combo_all
