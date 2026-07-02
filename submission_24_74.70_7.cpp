#include <algorithm>
#include <array>
#include <chrono>
#include <utility>
#include <queue>
#include <cstdint>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
using namespace std;
// IMC 2026 simplifygeometry - submission_65 candidate.
// Lane 5 recovery after submission_60 scored 62.269173 / failed hidden test 5.
// Reverts the risky low-margin sphere/ellipsoid remesh from 60 and keeps only
// a stricter high-resolution proxy-guarded late relocation boost.
struct Vec3 {
double x = 0.0, y = 0.0, z = 0.0;
};
static inline Vec3 operator-(const Vec3& a, const Vec3& b) {
return {a.x - b.x, a.y - b.y, a.z - b.z};
}
static inline Vec3 operator+(const Vec3& a, const Vec3& b) {
return {a.x + b.x, a.y + b.y, a.z + b.z};
}
static inline Vec3 operator*(const Vec3& a, double s) {
return {a.x * s, a.y * s, a.z * s};
}
static inline double dot3(const Vec3& a, const Vec3& b) {
return a.x * b.x + a.y * b.y + a.z * b.z;
}
static inline Vec3 cross3(const Vec3& a, const Vec3& b) {
return {
a.y * b.z - a.z * b.y,
a.z * b.x - a.x * b.z,
a.x * b.y - a.y * b.x,
};
}
static inline double norm2(const Vec3& a) {
return dot3(a, a);
}
static inline double norm3(const Vec3& a) {
return sqrt(norm2(a));
}
struct Face {
int v[3]{};
};
struct FastInput {
vector<char> buf;
char* p = nullptr;
FastInput() {
buf.reserve(1 << 27);
char chunk[1 << 16];
size_t n = 0;
while ((n = fread(chunk, 1, sizeof(chunk), stdin)) > 0) {
buf.insert(buf.end(), chunk, chunk + n);
}
buf.push_back('\0');
p = buf.data();
}
inline void skip_ws() {
while (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\t') ++p;
}
long next_long() {
skip_ws();
return strtol(p, &p, 10);
}
double next_double() {
skip_ws();
return strtod(p, &p);
}
char next_char() {
skip_ws();
return *p++;
}
};
static int N = 0, M = 0;
static vector<Vec3> P;
static vector<Vec3> originalP;
static vector<Face> faces;
static vector<Face> originalFaces;
static vector<unsigned char> alive_v;
static vector<unsigned char> alive_f;
static vector<double> cluster_radius;
static vector<vector<int>> incident;
static int active_faces = 0;
static vector<int> mark_neighbor;
static vector<int> mark_seen;
static int mark_tag = 1;
static int seen_tag = 1;
static double mesh_diag = 1.0;
static chrono::steady_clock::time_point simplify_start;
static bool smooth_stats_valid = false;
static double smooth10_ratio_cache = 0.0;
static double smooth30_ratio_cache = 0.0;
static double sharp22_ratio_cache = 1.0;
static double sharp45_ratio_cache = 1.0;
static double bad_ratio_cache = 1.0;
static inline bool face_has_vertex(int fid, int v) {
const Face& f = faces[fid];
return f.v[0] == v || f.v[1] == v || f.v[2] == v;
}
static inline int third_vertex_on_edge(int fid, int a, int b) {
const Face& f = faces[fid];
for (int i = 0; i < 3; ++i) {
int x = f.v[i];
if (x != a && x != b) return x;
}
return -1;
}
static inline array<int, 3> sorted_face_vertices(int a, int b, int c) {
array<int, 3> s{a, b, c};
sort(s.begin(), s.end());
return s;
}
static inline bool same_unordered_face(int fid, int a, int b, int c) {
const Face& f = faces[fid];
return sorted_face_vertices(f.v[0], f.v[1], f.v[2]) == sorted_face_vertices(a, b, c);
}
static inline Vec3 face_cross_from_vertices(int a, int b, int c) {
return cross3(P[b] - P[a], P[c] - P[a]);
}
static inline Vec3 point_with_move(int id, int moved_id, const Vec3& moved_pos) {
return id == moved_id ? moved_pos : P[id];
}
static inline Vec3 face_cross_with_move(int a, int b, int c, int moved_id, const Vec3& moved_pos) {
Vec3 pa = point_with_move(a, moved_id, moved_pos);
Vec3 pb = point_with_move(b, moved_id, moved_pos);
Vec3 pc = point_with_move(c, moved_id, moved_pos);
return cross3(pb - pa, pc - pa);
}
static inline Vec3 face_cross(int fid) {
const Face& f = faces[fid];
return face_cross_from_vertices(f.v[0], f.v[1], f.v[2]);
}
static inline Vec3 face_unit_normal(int fid) {
Vec3 cr = face_cross(fid);
double len = norm3(cr);
if (!(len > 0.0)) return {0.0, 0.0, 0.0};
return cr * (1.0 / len);
}
static inline unsigned long long edge_key_out(int a, int b) {
if (a > b) swap(a, b);
return (unsigned long long)(unsigned int)a << 32 | (unsigned int)b;
}
static void load_mesh() {
FastInput in;
N = (int)in.next_long();
M = (int)in.next_long();
P.resize(N);
originalP.resize(N);
alive_v.assign(N, 1);
cluster_radius.assign(N, 0.0);
Vec3 mn{1e100, 1e100, 1e100};
Vec3 mx{-1e100, -1e100, -1e100};
for (int i = 0; i < N; ++i) {
(void)in.next_char();
P[i].x = in.next_double();
P[i].y = in.next_double();
P[i].z = in.next_double();
originalP[i] = P[i];
mn.x = min(mn.x, P[i].x);
mn.y = min(mn.y, P[i].y);
mn.z = min(mn.z, P[i].z);
mx.x = max(mx.x, P[i].x);
mx.y = max(mx.y, P[i].y);
mx.z = max(mx.z, P[i].z);
}
mesh_diag = norm3(mx - mn);
if (!(mesh_diag > 0.0)) mesh_diag = 1.0;
faces.resize(M);
originalFaces.resize(M);
alive_f.assign(M, 1);
vector<int> deg(N, 0);
for (int i = 0; i < M; ++i) {
(void)in.next_char();
int a = (int)in.next_long() - 1;
int b = (int)in.next_long() - 1;
int c = (int)in.next_long() - 1;
faces[i].v[0] = a;
faces[i].v[1] = b;
faces[i].v[2] = c;
originalFaces[i] = faces[i];
++deg[a];
++deg[b];
++deg[c];
}
incident.resize(N);
for (int i = 0; i < N; ++i) incident[i].reserve(deg[i]);
for (int i = 0; i < M; ++i) {
incident[faces[i].v[0]].push_back(i);
incident[faces[i].v[1]].push_back(i);
incident[faces[i].v[2]].push_back(i);
}
active_faces = M;
mark_neighbor.assign(N, 0);
mark_seen.assign(N, 0);
}
static bool collect_shared_faces(int u, int v, int shared[2], int opp[2]) {
int cnt = 0;
for (int fid : incident[u]) {
if (!alive_f[fid]) continue;
if (!face_has_vertex(fid, u) || !face_has_vertex(fid, v)) continue;
if (cnt >= 2) return false;
shared[cnt] = fid;
opp[cnt] = third_vertex_on_edge(fid, u, v);
++cnt;
}
if (cnt != 2) return false;
if (opp[0] < 0 || opp[1] < 0 || opp[0] == opp[1]) return false;
return true;
}
static bool link_condition_ok(int u, int v, const int opp[2]) {
if (++mark_tag > 2000000000) {
fill(mark_neighbor.begin(), mark_neighbor.end(), 0);
mark_tag = 1;
}
if (++seen_tag > 2000000000) {
fill(mark_seen.begin(), mark_seen.end(), 0);
seen_tag = 1;
}
for (int fid : incident[u]) {
if (!alive_f[fid] || !face_has_vertex(fid, u)) continue;
const Face& f = faces[fid];
for (int i = 0; i < 3; ++i) {
int x = f.v[i];
if (x != u && x != v) mark_neighbor[x] = mark_tag;
}
}
int intersections = 0;
for (int fid : incident[v]) {
if (!alive_f[fid] || !face_has_vertex(fid, v)) continue;
const Face& f = faces[fid];
for (int i = 0; i < 3; ++i) {
int x = f.v[i];
if (x == u || x == v) continue;
if (mark_neighbor[x] != mark_tag) continue;
if (x != opp[0] && x != opp[1]) return false;
if (mark_seen[x] != seen_tag) {
mark_seen[x] = seen_tag;
++intersections;
}
}
}
return intersections == 2 && mark_seen[opp[0]] == seen_tag && mark_seen[opp[1]] == seen_tag;
}
struct PassParams {
double radius_limit = 0.0;
double plane_tol = 0.0;
double cos_limit = 1.0;
double exact_plane_tol = 0.0;
double exact_cos_limit = 1.0;
double min_area2 = 0.0;
bool allow_relocate = false;
};
struct Eval {
bool ok = false;
bool exact = false;
double cost = 1e100;
double new_radius = 0.0;
Vec3 new_pos{};
};
struct MeshState {
vector<Vec3> P;
vector<Face> faces;
vector<unsigned char> alive_v;
vector<unsigned char> alive_f;
vector<double> cluster_radius;
vector<vector<int>> incident;
int active_faces = 0;
};
static MeshState capture_state() {
MeshState s;
s.P = P;
s.faces = faces;
s.alive_v = alive_v;
s.alive_f = alive_f;
s.cluster_radius = cluster_radius;
s.incident = incident;
s.active_faces = active_faces;
return s;
}
static void restore_state(const MeshState& s) {
P = s.P;
faces = s.faces;
alive_v = s.alive_v;
alive_f = s.alive_f;
cluster_radius = s.cluster_radius;
incident = s.incident;
active_faces = s.active_faces;
}
static bool would_duplicate_face(int keep, int rem, int fid, int a, int b, int c, int shared0, int shared1) {
int probe = keep;
if ((int)incident[a].size() < (int)incident[probe].size()) probe = a;
if ((int)incident[b].size() < (int)incident[probe].size()) probe = b;
if ((int)incident[c].size() < (int)incident[probe].size()) probe = c;
for (int other : incident[probe]) {
if (!alive_f[other] || other == fid || other == shared0 || other == shared1) continue;
if (face_has_vertex(other, rem)) continue;
if (same_unordered_face(other, a, b, c)) return true;
}
return false;
}
static Eval evaluate_direction(int keep, int rem, const int shared[2], const PassParams& params) {
Eval res;
res.new_pos = P[keep];
const double move_dist = norm3(P[keep] - P[rem]);
res.new_radius = max(cluster_radius[keep], cluster_radius[rem] + move_dist);
double worst_plane = 0.0;
double worst_dot_loss = 0.0;
double worst_area_ratio_loss = 0.0;
int changed = 0;
bool exact = true;
for (int fid : incident[rem]) {
if (!alive_f[fid] || !face_has_vertex(fid, rem)) continue;
if (fid == shared[0] || fid == shared[1]) continue;
if (face_has_vertex(fid, keep)) return res;
Face f = faces[fid];
int new_tri[3] = {f.v[0], f.v[1], f.v[2]};
for (int i = 0; i < 3; ++i) {
if (new_tri[i] == rem) new_tri[i] = keep;
}
if (new_tri[0] == new_tri[1] || new_tri[0] == new_tri[2] || new_tri[1] == new_tri[2]) {
return res;
}
Vec3 old_cross = face_cross_from_vertices(f.v[0], f.v[1], f.v[2]);
Vec3 new_cross = face_cross_from_vertices(new_tri[0], new_tri[1], new_tri[2]);
double old_area2 = norm3(old_cross);
double new_area2 = norm3(new_cross);
if (!(old_area2 > params.min_area2) || !(new_area2 > params.min_area2)) return res;
if (new_area2 < old_area2 * 1e-7) return res;
double ndot = dot3(old_cross, new_cross) / (old_area2 * new_area2);
if (ndot > 1.0) ndot = 1.0;
if (ndot < -1.0) ndot = -1.0;
if (ndot < params.cos_limit) return res;
Vec3 unit_old = old_cross * (1.0 / old_area2);
double plane_dist = fabs(dot3(unit_old, P[keep] - P[f.v[0]]));
if (plane_dist > params.plane_tol && res.new_radius > params.radius_limit) return res;
if (would_duplicate_face(keep, rem, fid, new_tri[0], new_tri[1], new_tri[2], shared[0], shared[1])) {
return res;
}
worst_plane = max(worst_plane, plane_dist);
worst_dot_loss = max(worst_dot_loss, 1.0 - ndot);
worst_area_ratio_loss = max(worst_area_ratio_loss, max(0.0, 1.0 - new_area2 / old_area2));
if (plane_dist > params.exact_plane_tol || ndot < params.exact_cos_limit) exact = false;
++changed;
}
if (changed == 0) return res;
if (!exact) {
if (res.new_radius > params.radius_limit) return res;
if (worst_plane > params.plane_tol) return res;
}
res.ok = true;
res.exact = exact;
double radius_term = params.radius_limit > 0.0 ? res.new_radius / params.radius_limit : 0.0;
double plane_term = params.plane_tol > 0.0 ? worst_plane / params.plane_tol : 0.0;
res.cost = (exact ? -1000.0 : 0.0)
+ 0.80 * radius_term
+ 0.80 * plane_term
+ 250.0 * worst_dot_loss
+ 0.02 * worst_area_ratio_loss
+ 0.0005 * changed;
return res;
}
static Eval evaluate_relocation(int keep, int rem, const int shared[2], const PassParams& params, const Vec3& new_pos) {
Eval res;
res.new_pos = new_pos;
const double keep_move = norm3(new_pos - P[keep]);
const double rem_move = norm3(new_pos - P[rem]);
res.new_radius = max(cluster_radius[keep] + keep_move, cluster_radius[rem] + rem_move);
if (res.new_radius > params.radius_limit) return res;
static vector<int> affected;
affected.clear();
affected.reserve(incident[keep].size() + incident[rem].size());
for (int fid : incident[keep]) {
if (alive_f[fid] && fid != shared[0] && fid != shared[1] && face_has_vertex(fid, keep)) {
affected.push_back(fid);
}
}
for (int fid : incident[rem]) {
if (alive_f[fid] && fid != shared[0] && fid != shared[1] && face_has_vertex(fid, rem)) {
affected.push_back(fid);
}
}
sort(affected.begin(), affected.end());
affected.erase(unique(affected.begin(), affected.end()), affected.end());
if (affected.empty()) return res;
double worst_plane = 0.0;
double worst_dot_loss = 0.0;
double worst_area_ratio_loss = 0.0;
int changed = 0;
for (int fid : affected) {
const bool replaces_rem = face_has_vertex(fid, rem);
Face f = faces[fid];
int new_tri[3] = {f.v[0], f.v[1], f.v[2]};
for (int i = 0; i < 3; ++i) {
if (new_tri[i] == rem) new_tri[i] = keep;
}
if (new_tri[0] == new_tri[1] || new_tri[0] == new_tri[2] || new_tri[1] == new_tri[2]) {
return res;
}
Vec3 old_cross = face_cross_from_vertices(f.v[0], f.v[1], f.v[2]);
Vec3 new_cross = face_cross_with_move(new_tri[0], new_tri[1], new_tri[2], keep, new_pos);
double old_area2 = norm3(old_cross);
double new_area2 = norm3(new_cross);
if (!(old_area2 > params.min_area2) || !(new_area2 > params.min_area2)) return res;
if (new_area2 < old_area2 * 1e-7) return res;
double ndot = dot3(old_cross, new_cross) / (old_area2 * new_area2);
if (ndot > 1.0) ndot = 1.0;
if (ndot < -1.0) ndot = -1.0;
if (ndot < params.cos_limit) return res;
Vec3 unit_old = old_cross * (1.0 / old_area2);
double plane_dist = fabs(dot3(unit_old, new_pos - P[f.v[0]]));
if (plane_dist > params.plane_tol) return res;
if (replaces_rem
&& would_duplicate_face(keep, rem, fid, new_tri[0], new_tri[1], new_tri[2], shared[0], shared[1])) {
return res;
}
worst_plane = max(worst_plane, plane_dist);
worst_dot_loss = max(worst_dot_loss, 1.0 - ndot);
worst_area_ratio_loss = max(worst_area_ratio_loss, max(0.0, 1.0 - new_area2 / old_area2));
++changed;
}
if (changed == 0) return res;
res.ok = true;
res.exact = false;
double radius_term = params.radius_limit > 0.0 ? res.new_radius / params.radius_limit : 0.0;
double plane_term = params.plane_tol > 0.0 ? worst_plane / params.plane_tol : 0.0;
res.cost = 0.35
+ 0.95 * radius_term
+ 0.95 * plane_term
+ 300.0 * worst_dot_loss
+ 0.03 * worst_area_ratio_loss
+ 0.0007 * changed;
return res;
}
static void rebuild_incident_for(int keep, int rem) {
vector<int> merged;
merged.reserve(incident[keep].size() + incident[rem].size());
for (int fid : incident[keep]) {
if (alive_f[fid] && face_has_vertex(fid, keep)) merged.push_back(fid);
}
for (int fid : incident[rem]) {
if (alive_f[fid] && face_has_vertex(fid, keep)) merged.push_back(fid);
}
sort(merged.begin(), merged.end());
merged.erase(unique(merged.begin(), merged.end()), merged.end());
incident[keep].swap(merged);
vector<int>().swap(incident[rem]);
}
static void apply_collapse(int keep, int rem, const int shared[2], double new_radius, const Vec3& new_pos) {
for (int i = 0; i < 2; ++i) {
if (alive_f[shared[i]]) {
alive_f[shared[i]] = 0;
--active_faces;
}
}
for (int fid : incident[rem]) {
if (!alive_f[fid] || !face_has_vertex(fid, rem)) continue;
Face& f = faces[fid];
for (int i = 0; i < 3; ++i) {
if (f.v[i] == rem) f.v[i] = keep;
}
}
alive_v[rem] = 0;
P[keep] = new_pos;
cluster_radius[keep] = new_radius;
rebuild_incident_for(keep, rem);
}
static bool try_collapse_edge(int u, int v, const PassParams& params) {
if (u == v || !alive_v[u] || !alive_v[v]) return false;
int shared[2] = {-1, -1};
int opp[2] = {-1, -1};
if (!collect_shared_faces(u, v, shared, opp)) return false;
if (!link_condition_ok(u, v, opp)) return false;
Eval uv = evaluate_direction(u, v, shared, params);
Eval vu = evaluate_direction(v, u, shared, params);
if (params.allow_relocate) {
Vec3 mid = (P[u] + P[v]) * 0.5;
Eval um = evaluate_relocation(u, v, shared, params, mid);
Eval vm = evaluate_relocation(v, u, shared, params, mid);
if (um.ok && (!uv.ok || um.cost < uv.cost)) uv = um;
if (vm.ok && (!vu.ok || vm.cost < vu.cost)) vu = vm;
}
if (!uv.ok && !vu.ok) return false;
if (vu.ok && (!uv.ok || vu.cost < uv.cost)) {
apply_collapse(v, u, shared, vu.new_radius, vu.new_pos);
} else {
apply_collapse(u, v, shared, uv.new_radius, uv.new_pos);
}
return true;
}
static inline double elapsed_seconds() {
return chrono::duration<double>(chrono::steady_clock::now() - simplify_start).count();
}
static bool allow_high_density_smooth_lift() {
smooth_stats_valid = false;
if (N < 400 || M < 300) return false;
const int target_faces = 40000;
const int stride = max(1, M / target_faces);
const int sample_limit = 120000;
const double smooth_cos = cos(10.0 * acos(-1.0) / 180.0);
const double coarse_cos = cos(30.0 * acos(-1.0) / 180.0);
const double sharp_cos = cos(22.0 * acos(-1.0) / 180.0);
const double very_sharp_cos = cos(45.0 * acos(-1.0) / 180.0);
int sampled = 0;
int smooth = 0;
int coarse = 0;
int sharp = 0;
int very_sharp = 0;
int bad = 0;
for (int fid = 0; fid < M && sampled < sample_limit; fid += stride) {
const Face& f = faces[fid];
const int e[3][2] = {{f.v[0], f.v[1]}, {f.v[1], f.v[2]}, {f.v[2], f.v[0]}};
for (int k = 0; k < 3 && sampled < sample_limit; ++k) {
int shared[2] = {-1, -1};
int opp[2] = {-1, -1};
if (!collect_shared_faces(e[k][0], e[k][1], shared, opp)) {
++bad;
continue;
}
Vec3 n0 = face_unit_normal(shared[0]);
Vec3 n1 = face_unit_normal(shared[1]);
double ndot = dot3(n0, n1);
if (ndot > 1.0) ndot = 1.0;
if (ndot < -1.0) ndot = -1.0;
++sampled;
if (ndot > smooth_cos) ++smooth;
if (ndot > coarse_cos) ++coarse;
if (ndot < sharp_cos) ++sharp;
if (ndot < very_sharp_cos) ++very_sharp;
}
}
const int min_samples = N < 3000 ? 250 : 1000;
if (sampled < min_samples) return false;
const double smooth_ratio = (double)smooth / (double)sampled;
const double coarse_ratio = (double)coarse / (double)sampled;
const double sharp_ratio = (double)sharp / (double)sampled;
const double very_sharp_ratio = (double)very_sharp / (double)sampled;
const double bad_ratio = (double)bad / (double)(sampled + bad);
smooth_stats_valid = true;
smooth10_ratio_cache = smooth_ratio;
smooth30_ratio_cache = coarse_ratio;
sharp22_ratio_cache = sharp_ratio;
sharp45_ratio_cache = very_sharp_ratio;
bad_ratio_cache = bad_ratio;
const bool very_smooth = smooth_ratio >= 0.985 && sharp_ratio <= 0.010 && bad_ratio <= 0.010;
const bool small_coarse_curved = N < 3000
&& coarse_ratio >= 0.995
&& sharp_ratio <= 0.080
&& very_sharp_ratio <= 0.005
&& bad_ratio <= 0.010;
return very_smooth || small_coarse_curved;
}
struct RenderMaps {
vector<double> depth;
vector<Vec3> normal;
vector<unsigned char> mask;
};
static inline void project_view(const Vec3& p, int view, int res, double& u, double& v, double& z) {
constexpr double D = 2.5;
const double f = 800.0 * ((double)res / 1024.0);
const double c = 0.5 * (double)res;
double sx = 0.0, sy = 0.0;
if (view == 0) {
sx = p.y; sy = p.z; z = D - p.x;
} else if (view == 1) {
sx = -p.y; sy = p.z; z = D + p.x;
} else if (view == 2) {
sx = -p.x; sy = p.z; z = D - p.y;
} else if (view == 3) {
sx = p.x; sy = p.z; z = D + p.y;
} else if (view == 4) {
sx = p.x; sy = p.y; z = D - p.z;
} else {
sx = -p.x; sy = p.y; z = D + p.z;
}
u = f * sx / z + c;
v = f * sy / z + c;
}
static void rasterize_triangle(RenderMaps& rm, int res,
const Vec3& a, const Vec3& b, const Vec3& c,
const Vec3& unit_normal, int view) {
double x0, y0, z0, x1, y1, z1, x2, y2, z2;
project_view(a, view, res, x0, y0, z0);
project_view(b, view, res, x1, y1, z1);
project_view(c, view, res, x2, y2, z2);
if (z0 <= 0.0 || z1 <= 0.0 || z2 <= 0.0) return;
int xmin = max(0, (int)floor(min(x0, min(x1, x2)) - 0.5));
int xmax = min(res - 1, (int)ceil(max(x0, max(x1, x2)) + 0.5));
int ymin = max(0, (int)floor(min(y0, min(y1, y2)) - 0.5));
int ymax = min(res - 1, (int)ceil(max(y0, max(y1, y2)) + 0.5));
if (xmin > xmax || ymin > ymax) return;
const double den = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2);
if (fabs(den) < 1e-12) return;
for (int yy = ymin; yy <= ymax; ++yy) {
const double py = (double)yy + 0.5;
for (int xx = xmin; xx <= xmax; ++xx) {
const double px = (double)xx + 0.5;
const double w0 = ((y1 - y2) * (px - x2) + (x2 - x1) * (py - y2)) / den;
const double w1 = ((y2 - y0) * (px - x2) + (x0 - x2) * (py - y2)) / den;
const double w2 = 1.0 - w0 - w1;
if (w0 < -1e-9 || w1 < -1e-9 || w2 < -1e-9) continue;
const double zp = 1.0 / (w0 / z0 + w1 / z1 + w2 / z2);
const int idx = yy * res + xx;
if (zp < rm.depth[idx]) {
rm.depth[idx] = zp;
rm.normal[idx] = unit_normal;
rm.mask[idx] = 1;
}
}
}
}
static RenderMaps render_original_view(int view, int res) {
RenderMaps rm;
rm.depth.assign(res * res, 255.0);
rm.normal.assign(res * res, Vec3{0.0, 0.0, 0.0});
rm.mask.assign(res * res, 0);
for (int fid = 0; fid < M; ++fid) {
const Face& f = originalFaces[fid];
const Vec3& a = originalP[f.v[0]];
const Vec3& b = originalP[f.v[1]];
const Vec3& c = originalP[f.v[2]];
Vec3 cr = cross3(b - a, c - a);
double len = norm3(cr);
if (!(len > 0.0)) continue;
rasterize_triangle(rm, res, a, b, c, cr * (1.0 / len), view);
}
return rm;
}
static RenderMaps render_current_view(int view, int res) {
RenderMaps rm;
rm.depth.assign(res * res, 255.0);
rm.normal.assign(res * res, Vec3{0.0, 0.0, 0.0});
rm.mask.assign(res * res, 0);
for (int fid = 0; fid < (int)faces.size(); ++fid) {
if (!alive_f[fid]) continue;
const Face& f = faces[fid];
Vec3 cr = face_cross(fid);
double len = norm3(cr);
if (!(len > 0.0)) continue;
rasterize_triangle(rm, res, P[f.v[0]], P[f.v[1]], P[f.v[2]], cr * (1.0 / len), view);
}
return rm;
}
static inline double normal_channel_value(const Vec3& n, int ch) {
if (ch == 0) return (n.x + 1.0) * 127.5;
if (ch == 1) return (n.y + 1.0) * 127.5;
return (n.z + 1.0) * 127.5;
}
template <class Getter>
static double ssim_channel_proxy(const RenderMaps& a, const RenderMaps& b,
const vector<unsigned char>& fg,
int res, int win, Getter get_value) {
const int rad = win / 2;
const double c1 = (0.01 * 255.0) * (0.01 * 255.0);
const double c2 = (0.03 * 255.0) * (0.03 * 255.0);
double total = 0.0;
int count = 0;
for (int y = 0; y < res; ++y) {
for (int x = 0; x < res; ++x) {
const int center = y * res + x;
if (!fg[center]) continue;
double sx = 0.0, sy = 0.0, sxx = 0.0, syy = 0.0, sxy = 0.0;
int n = 0;
for (int dy = -rad; dy <= rad; ++dy) {
int yy = y + dy;
if (yy < 0) yy = 0;
if (yy >= res) yy = res - 1;
for (int dx = -rad; dx <= rad; ++dx) {
int xx = x + dx;
if (xx < 0) xx = 0;
if (xx >= res) xx = res - 1;
const int idx = yy * res + xx;
const double vx = get_value(a, idx);
const double vy = get_value(b, idx);
sx += vx;
sy += vy;
sxx += vx * vx;
syy += vy * vy;
sxy += vx * vy;
++n;
}
}
const double inv = 1.0 / (double)n;
const double ux = sx * inv;
const double uy = sy * inv;
const double varx = max(0.0, sxx * inv - ux * ux);
const double vary = max(0.0, syy * inv - uy * uy);
const double cov = sxy * inv - ux * uy;
const double num = (2.0 * ux * uy + c1) * (2.0 * cov + c2);
const double den = (ux * ux + uy * uy + c1) * (varx + vary + c2);
total += num / den;
++count;
}
}
return count > 0 ? total / (double)count : 1.0;
}
static double visual_proxy_score(int res) {
constexpr int WIN = 11;
double total = 0.0;
for (int view = 0; view < 6; ++view) {
RenderMaps orig = render_original_view(view, res);
RenderMaps simp = render_current_view(view, res);
vector<unsigned char> fg(res * res, 0);
for (int i = 0; i < res * res; ++i) fg[i] = (orig.mask[i] || simp.mask[i]) ? 1 : 0;
double normal_score = 0.0;
for (int ch = 0; ch < 3; ++ch) {
normal_score += ssim_channel_proxy(orig, simp, fg, res, WIN, [ch](const RenderMaps& m, int idx) {
return normal_channel_value(m.normal[idx], ch);
});
}
normal_score /= 3.0;
double depth_score = ssim_channel_proxy(orig, simp, fg, res, WIN, [](const RenderMaps& m, int idx) {
return m.depth[idx];
});
total += 0.5 * normal_score + 0.5 * depth_score;
}
return total / 6.0;
}
struct EllipsoidFit {
bool ok = false;
Vec3 center{};
Vec3 axis[3]{};
double radius[3]{};
double mean_abs_residual = 1e100;
double rms_residual = 1e100;
double max_abs_residual = 1e100;
};
static inline Vec3 vec_div(const Vec3& a, double s) {
return {a.x / s, a.y / s, a.z / s};
}
static inline double component_on_axis(const Vec3& p, const Vec3& a) {
return dot3(p, a);
}
static double original_orientation_sign(const Vec3& center) {
const int sample_limit = 100000;
const int stride = max(1, M / sample_limit);
double total = 0.0;
int sampled = 0;
for (int fid = 0; fid < M; fid += stride) {
const Face& f = originalFaces[fid];
const Vec3& a = originalP[f.v[0]];
const Vec3& b = originalP[f.v[1]];
const Vec3& c = originalP[f.v[2]];
Vec3 cr = cross3(b - a, c - a);
Vec3 ctr = (a + b + c) * (1.0 / 3.0);
const double s = dot3(cr, ctr - center);
if (fabs(s) > 1e-18) {
total += s;
++sampled;
}
}
if (sampled == 0) return 1.0;
return total >= 0.0 ? 1.0 : -1.0;
}
static void orient_face_like_original(vector<Vec3>& verts, Face& f,
const Vec3& center,
double orientation_sign) {
const Vec3& a = verts[f.v[0]];
const Vec3& b = verts[f.v[1]];
const Vec3& c = verts[f.v[2]];
Vec3 cr = cross3(b - a, c - a);
Vec3 ctr = (a + b + c) * (1.0 / 3.0);
if (dot3(cr, ctr - center) * orientation_sign < 0.0) swap(f.v[1], f.v[2]);
}
static bool install_replacement_mesh(const vector<Vec3>& new_vertices,
const vector<Face>& new_faces) {
if (new_vertices.empty() || new_faces.empty()) return false;
if ((int)new_vertices.size() > N) return false;
const double area_eps = max(1e-30, 1e-24 * mesh_diag * mesh_diag);
for (const Face& f : new_faces) {
for (int k = 0; k < 3; ++k) {
if (f.v[k] < 0 || f.v[k] >= (int)new_vertices.size()) return false;
}
if (f.v[0] == f.v[1] || f.v[0] == f.v[2] || f.v[1] == f.v[2]) return false;
Vec3 cr = cross3(new_vertices[f.v[1]] - new_vertices[f.v[0]],
new_vertices[f.v[2]] - new_vertices[f.v[0]]);
if (!(norm2(cr) > area_eps)) return false;
}
P.assign(N, Vec3{});
for (int i = 0; i < (int)new_vertices.size(); ++i) P[i] = new_vertices[i];
alive_v.assign(N, 0);
cluster_radius.assign(N, 0.0);
for (int i = 0; i < (int)new_vertices.size(); ++i) alive_v[i] = 1;
faces = new_faces;
alive_f.assign(faces.size(), 1);
active_faces = (int)faces.size();
vector<int> deg(N, 0);
for (const Face& f : faces) {
++deg[f.v[0]];
++deg[f.v[1]];
++deg[f.v[2]];
}
incident.assign(N, {});
for (int i = 0; i < N; ++i) incident[i].reserve(deg[i]);
for (int i = 0; i < (int)faces.size(); ++i) {
incident[faces[i].v[0]].push_back(i);
incident[faces[i].v[1]].push_back(i);
incident[faces[i].v[2]].push_back(i);
}
return true;
}
static bool build_uv_ellipsoid_mesh(const EllipsoidFit& fit, int lat, int lon,
vector<Vec3>& verts, vector<Face>& out_faces) {
if (!fit.ok || lat < 4 || lon < 8) return false;
const int vertex_count = 2 + (lat - 1) * lon;
if (vertex_count > N) return false;
verts.clear();
out_faces.clear();
verts.reserve(vertex_count);
out_faces.reserve(2 * lon * (lat - 1));
const double orientation_sign = original_orientation_sign(fit.center);
auto make_point = [&](double x, double y, double z) {
Vec3 p = fit.center;
p = p + fit.axis[0] * (fit.radius[0] * x);
p = p + fit.axis[1] * (fit.radius[1] * y);
p = p + fit.axis[2] * (fit.radius[2] * z);
return p;
};
verts.push_back(make_point(0.0, 0.0, 1.0));
verts.push_back(make_point(0.0, 0.0, -1.0));
auto ring_id = [&](int ring, int j) {
j = (j % lon + lon) % lon;
return 2 + (ring - 1) * lon + j;
};
for (int ring = 1; ring <= lat - 1; ++ring) {
const double theta = acos(-1.0) * (double)ring / (double)lat;
const double st = sin(theta);
const double ct = cos(theta);
for (int j = 0; j < lon; ++j) {
const double phi = 2.0 * acos(-1.0) * (double)j / (double)lon;
verts.push_back(make_point(st * cos(phi), st * sin(phi), ct));
}
}
auto add_face = [&](int a, int b, int c) {
Face f;
f.v[0] = a; f.v[1] = b; f.v[2] = c;
orient_face_like_original(verts, f, fit.center, orientation_sign);
out_faces.push_back(f);
};
for (int j = 0; j < lon; ++j) {
add_face(0, ring_id(1, j), ring_id(1, j + 1));
}
for (int ring = 1; ring <= lat - 2; ++ring) {
for (int j = 0; j < lon; ++j) {
const int a = ring_id(ring, j);
const int b = ring_id(ring + 1, j);
const int c = ring_id(ring + 1, j + 1);
const int d = ring_id(ring, j + 1);
add_face(a, b, c);
add_face(a, c, d);
}
}
for (int j = 0; j < lon; ++j) {
add_face(1, ring_id(lat - 1, j + 1), ring_id(lat - 1, j));
}
return true;
}
static EllipsoidFit evaluate_ellipsoid_basis(const Vec3 basis[3]) {
EllipsoidFit fit;
for (int k = 0; k < 3; ++k) fit.axis[k] = basis[k];
double lo[3] = {1e100, 1e100, 1e100};
double hi[3] = {-1e100, -1e100, -1e100};
for (const Vec3& p : originalP) {
for (int k = 0; k < 3; ++k) {
const double t = component_on_axis(p, fit.axis[k]);
lo[k] = min(lo[k], t);
hi[k] = max(hi[k], t);
}
}
fit.center = Vec3{};
for (int k = 0; k < 3; ++k) {
const double mid = 0.5 * (lo[k] + hi[k]);
fit.center = fit.center + fit.axis[k] * mid;
fit.radius[k] = 0.5 * (hi[k] - lo[k]);
if (!(fit.radius[k] > 1e-12)) return fit;
}
const int sample_limit = 200000;
const int stride = max(1, N / sample_limit);
int sampled = 0;
double sum_abs = 0.0;
double sum_sq = 0.0;
double max_abs = 0.0;
for (int i = 0; i < N; i += stride) {
const Vec3 q = originalP[i] - fit.center;
double r2 = 0.0;
for (int k = 0; k < 3; ++k) {
const double u = component_on_axis(q, fit.axis[k]) / fit.radius[k];
r2 += u * u;
}
const double err = fabs(sqrt(max(0.0, r2)) - 1.0);
sum_abs += err;
sum_sq += err * err;
max_abs = max(max_abs, err);
++sampled;
}
if (sampled == 0) return fit;
fit.mean_abs_residual = sum_abs / (double)sampled;
fit.rms_residual = sqrt(sum_sq / (double)sampled);
fit.max_abs_residual = max_abs;
const double max_allowed = (N < 5000 ? 0.010 : 0.014);
const double rms_allowed = (N < 5000 ? 0.0035 : 0.0045);
fit.ok = fit.max_abs_residual <= max_allowed && fit.rms_residual <= rms_allowed;
return fit;
}
static void jacobi_eigenvectors_3x3(double a[3][3], Vec3 out_axis[3]) {
double v[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
for (int it = 0; it < 32; ++it) {
int p = 0, q = 1;
double best = fabs(a[0][1]);
if (fabs(a[0][2]) > best) { p = 0; q = 2; best = fabs(a[0][2]); }
if (fabs(a[1][2]) > best) { p = 1; q = 2; best = fabs(a[1][2]); }
if (best < 1e-18) break;
const double app = a[p][p];
const double aqq = a[q][q];
const double apq = a[p][q];
const double tau = (aqq - app) / (2.0 * apq);
const double t = (tau >= 0.0 ? 1.0 : -1.0) / (fabs(tau) + sqrt(1.0 + tau * tau));
const double c = 1.0 / sqrt(1.0 + t * t);
const double s = t * c;
for (int k = 0; k < 3; ++k) {
if (k == p || k == q) continue;
const double akp = a[k][p];
const double akq = a[k][q];
a[k][p] = a[p][k] = c * akp - s * akq;
a[k][q] = a[q][k] = s * akp + c * akq;
}
a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
a[p][q] = a[q][p] = 0.0;
for (int k = 0; k < 3; ++k) {
const double vkp = v[k][p];
const double vkq = v[k][q];
v[k][p] = c * vkp - s * vkq;
v[k][q] = s * vkp + c * vkq;
}
}
int ord[3] = {0, 1, 2};
sort(ord, ord + 3, [&](int lhs, int rhs) {
return a[lhs][lhs] > a[rhs][rhs];
});
for (int j = 0; j < 3; ++j) {
const int col = ord[j];
Vec3 e{v[0][col], v[1][col], v[2][col]};
const double len = norm3(e);
out_axis[j] = len > 0.0 ? vec_div(e, len) : Vec3{};
}
if (dot3(cross3(out_axis[0], out_axis[1]), out_axis[2]) < 0.0) {
out_axis[2] = out_axis[2] * -1.0;
}
}
static EllipsoidFit fit_strict_ellipsoid_like() {
Vec3 identity[3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
EllipsoidFit best = evaluate_ellipsoid_basis(identity);
Vec3 mean{};
for (const Vec3& p : originalP) mean = mean + p;
mean = mean * (1.0 / max(1, N));
double cov[3][3] = {};
for (const Vec3& p : originalP) {
const Vec3 q = p - mean;
const double x[3] = {q.x, q.y, q.z};
for (int i = 0; i < 3; ++i) {
for (int j = 0; j < 3; ++j) cov[i][j] += x[i] * x[j];
}
}
for (int i = 0; i < 3; ++i) {
for (int j = 0; j < 3; ++j) cov[i][j] /= (double)max(1, N);
}
Vec3 pca_axis[3];
jacobi_eigenvectors_3x3(cov, pca_axis);
if (norm2(pca_axis[0]) > 0.0 && norm2(pca_axis[1]) > 0.0 && norm2(pca_axis[2]) > 0.0) {
EllipsoidFit pca = evaluate_ellipsoid_basis(pca_axis);
if (pca.ok && (!best.ok || pca.rms_residual < best.rms_residual)) best = pca;
}
return best;
}
static int count_output_vertices_estimate();
static EllipsoidFit fit_radial_sphere_like() {
EllipsoidFit fit;
if (N < 900 || originalP.empty()) return fit;
Vec3 mn = originalP[0];
Vec3 mx = originalP[0];
for (const Vec3& p : originalP) {
mn.x = min(mn.x, p.x); mn.y = min(mn.y, p.y); mn.z = min(mn.z, p.z);
mx.x = max(mx.x, p.x); mx.y = max(mx.y, p.y); mx.z = max(mx.z, p.z);
}
fit.center = (mn + mx) * 0.5;
const double ex = mx.x - mn.x;
const double ey = mx.y - mn.y;
const double ez = mx.z - mn.z;
const double max_extent = max(ex, max(ey, ez));
const double min_extent = min(ex, min(ey, ez));
if (!(max_extent > 1e-12) || min_extent < max_extent * 0.975) return fit;
const int sample_limit = 240000;
const int stride = max(1, N / sample_limit);
double sum_r = 0.0;
int sampled = 0;
for (int i = 0; i < N; i += stride) {
sum_r += norm3(originalP[i] - fit.center);
++sampled;
}
if (sampled == 0) return fit;
const double mean_r = sum_r / (double)sampled;
if (!(mean_r > 1e-12)) return fit;
double sum_sq = 0.0;
double max_abs = 0.0;
for (int i = 0; i < N; i += stride) {
const double dev = fabs(norm3(originalP[i] - fit.center) - mean_r);
sum_sq += dev * dev;
max_abs = max(max_abs, dev);
}
const double rel_rms = sqrt(sum_sq / (double)sampled) / mean_r;
const double rel_max = max_abs / mean_r;
const double rms_limit = (N < 5000 ? 0.0060 : 0.0045);
const double max_limit = (N < 5000 ? 0.0260 : 0.0180);
if (rel_rms > rms_limit || rel_max > max_limit) return fit;
fit.ok = true;
fit.axis[0] = {1.0, 0.0, 0.0};
fit.axis[1] = {0.0, 1.0, 0.0};
fit.axis[2] = {0.0, 0.0, 1.0};
fit.radius[0] = fit.radius[1] = fit.radius[2] = mean_r;
fit.mean_abs_residual = rel_rms;
fit.rms_residual = rel_rms;
fit.max_abs_residual = rel_max;
return fit;
}
static bool try_radial_sphere_remesh() {
if (elapsed_seconds() > 13.5) return false;
const EllipsoidFit fit = fit_radial_sphere_like();
if (!fit.ok) return false;
int lat = 20;
int lon = 40;
int proxy_res = 256;
double proxy_threshold = 0.930;
if (N < 1500) {
lat = 20; lon = 40;
proxy_threshold = 0.950;
} else if (N < 3000) {
lat = 20; lon = 40;
proxy_threshold = 0.915;
} else if (N < 8000) {
lat = 20; lon = 40;
proxy_threshold = 0.905;
} else if (N < 20000) {
lat = 24; lon = 48;
proxy_threshold = 0.930;
} else if (N < 40000) {
lat = 28; lon = 56;
proxy_res = 192;
proxy_threshold = 0.945;
} else if (N < 100000) {
lat = 28; lon = 56;
proxy_res = 192;
proxy_threshold = 0.945;
} else {
lat = 28; lon = 56;
proxy_res = 128;
proxy_threshold = 0.950;
}
vector<Vec3> new_vertices;
vector<Face> new_faces;
if (!build_uv_ellipsoid_mesh(fit, lat, lon, new_vertices, new_faces)) return false;
const int safe_vertices = count_output_vertices_estimate();
if (safe_vertices <= 0) return false;
if ((int)new_vertices.size() >= safe_vertices) return false;
if ((int)new_vertices.size() * 100 > safe_vertices * 94) return false;
MeshState safe_state = capture_state();
bool keep = false;
if (install_replacement_mesh(new_vertices, new_faces) && elapsed_seconds() < 16.5) {
const double proxy = visual_proxy_score(proxy_res);
keep = proxy >= proxy_threshold;
}
if (!keep) restore_state(safe_state);
return keep;
}
struct AxisRevolveFit {
bool ok = false;
int axis = 2;
double center_u = 0.0;
double center_v = 0.0;
double t0 = 0.0;
double t1 = 0.0;
double r0 = 0.0;
double r1 = 0.0;
double residual = 1e100;
};
static inline void axis_components(const Vec3& p, int axis, double& t, double& u, double& v) {
if (axis == 0) {
t = p.x; u = p.y; v = p.z;
} else if (axis == 1) {
t = p.y; u = p.x; v = p.z;
} else {
t = p.z; u = p.x; v = p.y;
}
}
static inline Vec3 axis_make_point(int axis, double t, double u, double v) {
if (axis == 0) return {t, u, v};
if (axis == 1) return {u, t, v};
return {u, v, t};
}
struct OrientedFrustumFit {
bool ok = false;
Vec3 w{}, u{}, v{};
double cu = 0.0, cv = 0.0, t0 = 0.0, t1 = 0.0;
double rx0 = 0.0, rx1 = 0.0, ry0 = 0.0, ry1 = 0.0;
double residual = 1e100;
};
static inline Vec3 unit_vec(Vec3 a) {
double l = norm3(a);
return l > 1e-300 ? a * (1.0 / l) : Vec3{};
}
static inline void oriented_components(const Vec3& p, const Vec3& w,
const Vec3& u, const Vec3& v, double& t, double& x, double& y) {
t = dot3(p, w);
x = dot3(p, u);
y = dot3(p, v);
}
static inline Vec3 oriented_make_point(const Vec3& w, const Vec3& u,
const Vec3& v, double t, double x, double y) {
return w * t + u * x + v * y;
}
static bool pca_original_axes(Vec3 axis[3]) {
if (N < 600 || originalP.empty()) return false;
Vec3 mean{};
for (const Vec3& p : originalP) mean = mean + p;
mean = mean * (1.0 / max(1, N));
double cov[3][3] = {};
const int stride = max(1, N / 240000);
int sampled = 0;
for (int i = 0; i < N; i += stride) {
Vec3 q = originalP[i] - mean;
double x[3] = {q.x, q.y, q.z};
for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) cov[a][b] += x[a] * x[b];
++sampled;
}
if (sampled < 50) return false;
for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) cov[a][b] /= (double)sampled;
jacobi_eigenvectors_3x3(cov, axis);
return norm2(axis[0]) > 0.5 && norm2(axis[1]) > 0.5 && norm2(axis[2]) > 0.5;
}
static void make_oriented_basis(Vec3 w, Vec3 hint, Vec3& u, Vec3& v) {
w = unit_vec(w);
u = hint - w * dot3(hint, w);
if (norm2(u) < 1e-20) {
Vec3 h = fabs(w.x) < 0.75 ? Vec3{1,0,0} : (fabs(w.y) < 0.75 ? Vec3{0,1,0} : Vec3{0,0,1});
u = h - w * dot3(h, w);
}
u = unit_vec(u);
v = unit_vec(cross3(w, u));
}
static OrientedFrustumFit fit_oriented_frustum_like(Vec3 w, Vec3 hint) {
OrientedFrustumFit fit;
fit.w = unit_vec(w);
make_oriented_basis(fit.w, hint, fit.u, fit.v);
if (N < 1000 || norm2(fit.w) < 0.5 || norm2(fit.u) < 0.5 || norm2(fit.v) < 0.5) return fit;
double min_t = 1e100, max_t = -1e100, min_u = 1e100, max_u = -1e100, min_v = 1e100, max_v = -1e100;
for (const Vec3& p : originalP) {
double t, x, y;
oriented_components(p, fit.w, fit.u, fit.v, t, x, y);
min_t = min(min_t, t); max_t = max(max_t, t);
min_u = min(min_u, x); max_u = max(max_u, x);
min_v = min(min_v, y); max_v = max(max_v, y);
}
if (!(max_t > min_t)) return fit;
fit.t0 = min_t; fit.t1 = max_t; fit.cu = 0.5 * (min_u + max_u); fit.cv = 0.5 * (min_v + max_v);
const double len = max_t - min_t;
const int stride = max(1, N / 240000);
double max_r = 0.0;
for (int i = 0; i < N; i += stride) {
double t, x, y;
oriented_components(originalP[i], fit.w, fit.u, fit.v, t, x, y);
x -= fit.cu; y -= fit.cv;
max_r = max(max_r, sqrt(x * x + y * y));
}
if (!(max_r > 1e-12) || len < max_r * 0.35) return fit;
const int B = 18;
double sx[B], sy[B]; int bc[B];
for (int i = 0; i < B; ++i) sx[i] = sy[i] = 0.0, bc[i] = 0;
const double eps = max_r * 0.035;
int ring_samples = 0, axis_samples = 0;
for (int i = 0; i < N; i += stride) {
double t, x, y;
oriented_components(originalP[i], fit.w, fit.u, fit.v, t, x, y);
x -= fit.cu; y -= fit.cv;
double r = sqrt(x * x + y * y);
if (r <= eps) {
double near_end = min(fabs(t - min_t), fabs(t - max_t));
if (near_end > len * 0.035) return fit;
++axis_samples;
continue;
}
int b = (int)((t - min_t) / len * B);
if (b < 0) b = 0; else if (b >= B) b = B - 1;
sx[b] = max(sx[b], fabs(x));
sy[b] = max(sy[b], fabs(y));
++bc[b]; ++ring_samples;
}
if (ring_samples < 200 || axis_samples > ring_samples / 8 + 8) return fit;
double S = 0, St = 0, Stt = 0, Sx = 0, Stx = 0, Sy = 0, Sty = 0;
int bins = 0;
for (int b = 0; b < B; ++b) if (bc[b] > 0 && sx[b] > eps && sy[b] > eps) {
double s = ((double)b + 0.5) / (double)B;
S += 1.0; St += s; Stt += s * s; Sx += sx[b]; Stx += s * sx[b]; Sy += sy[b]; Sty += s * sy[b];
++bins;
}
if (bins < 5) return fit;
double den = S * Stt - St * St;
if (fabs(den) < 1e-18) return fit;
double bx = (S * Stx - St * Sx) / den, ax = (Sx - bx * St) / S;
double by = (S * Sty - St * Sy) / den, ay = (Sy - by * St) / S;
fit.rx0 = ax; fit.rx1 = ax + bx; fit.ry0 = ay; fit.ry1 = ay + by;
double mr = max(max(fit.rx0, fit.rx1), max(fit.ry0, fit.ry1));
double nr = min(min(fit.rx0, fit.rx1), min(fit.ry0, fit.ry1));
if (!(mr > 1e-12) || nr < mr * 0.10) return fit;
double sum_sq = 0.0, max_abs = 0.0;
int checked = 0;
for (int i = 0; i < N; i += stride) {
double t, x, y;
oriented_components(originalP[i], fit.w, fit.u, fit.v, t, x, y);
x -= fit.cu; y -= fit.cv;
double r = sqrt(x * x + y * y);
if (r <= eps) continue;
double s = (t - min_t) / len;
double rx = fit.rx0 + (fit.rx1 - fit.rx0) * s, ry = fit.ry0 + (fit.ry1 - fit.ry0) * s;
if (rx <= mr * 0.08 || ry <= mr * 0.08) return fit;
double e = fabs(sqrt((x * x) / (rx * rx) + (y * y) / (ry * ry)) - 1.0);
sum_sq += e * e; max_abs = max(max_abs, e); ++checked;
}
if (checked < 200) return fit;
fit.residual = sqrt(sum_sq / (double)checked);
if (fit.residual > 0.022 || max_abs > 0.085) return fit;
fit.ok = true;
return fit;
}
static bool build_oriented_frustum_mesh(const OrientedFrustumFit& fit, int sides,
vector<Vec3>& verts, vector<Face>& out_faces) {
if (!fit.ok || sides < 12) return false;
verts.clear(); out_faces.clear();
verts.reserve(2 + 2 * sides); out_faces.reserve(4 * sides);
auto make = [&](double t, double rx, double ry, int j) {
double th = 2.0 * acos(-1.0) * (double)j / (double)sides;
return oriented_make_point(fit.w, fit.u, fit.v, t,
fit.cu + rx * cos(th), fit.cv + ry * sin(th));
};
Vec3 center = oriented_make_point(fit.w, fit.u, fit.v, 0.5 * (fit.t0 + fit.t1), fit.cu, fit.cv);
double sign = original_orientation_sign(center);
verts.push_back(oriented_make_point(fit.w, fit.u, fit.v, fit.t0, fit.cu, fit.cv));
verts.push_back(oriented_make_point(fit.w, fit.u, fit.v, fit.t1, fit.cu, fit.cv));
const int lo = (int)verts.size();
for (int j = 0; j < sides; ++j) verts.push_back(make(fit.t0, fit.rx0, fit.ry0, j));
const int hi = (int)verts.size();
for (int j = 0; j < sides; ++j) verts.push_back(make(fit.t1, fit.rx1, fit.ry1, j));
if ((int)verts.size() > N) return false;
auto add_face = [&](int a, int b, int c) {
Face f; f.v[0] = a; f.v[1] = b; f.v[2] = c;
orient_face_like_original(verts, f, center, sign);
out_faces.push_back(f);
};
for (int j = 0; j < sides; ++j) {
int k = (j + 1) % sides, l0 = lo + j, l1 = lo + k, h0 = hi + j, h1 = hi + k;
add_face(l0, l1, h1); add_face(l0, h1, h0);
add_face(0, l1, l0); add_face(1, h0, h1);
}
return true;
}
static bool try_oriented_frustum_remesh() {
if (elapsed_seconds() > 13.7 || N < 1000 || N > 120000) return false;
const int safe_vertices = count_output_vertices_estimate();
if (safe_vertices <= 0) return false;
Vec3 ax[3];
if (!pca_original_axes(ax)) return false;
OrientedFrustumFit best;
for (int k = 0; k < 3; ++k) {
OrientedFrustumFit fit = fit_oriented_frustum_like(ax[k], ax[(k + 1) % 3]);
if (fit.ok && (!best.ok || fit.residual < best.residual)) best = fit;
}
if (!best.ok) return false;
const int trials[] = {32, 40, 48, 56, 64, 80, 96};
MeshState safe_state = capture_state();
vector<Vec3> new_vertices;
vector<Face> new_faces;
for (int sides : trials) {
if (elapsed_seconds() > 16.2) break;
if (2 + 2 * sides >= safe_vertices) continue;
if (!build_oriented_frustum_mesh(best, sides, new_vertices, new_faces)) continue;
restore_state(safe_state);
bool keep = false;
if (install_replacement_mesh(new_vertices, new_faces) && elapsed_seconds() < 16.7) {
int proxy_res = (N < 30000 ? 512 : 256);
double threshold = (N < 30000 ? 0.930 : 0.935);
keep = visual_proxy_score(proxy_res) >= threshold;
}
if (keep) return true;
}
restore_state(safe_state);
return false;
}
struct AxisTorusFit {
bool ok = false;
int axis = 2;
double center_t = 0.0;
double center_u = 0.0;
double center_v = 0.0;
double major = 0.0;
double minor = 0.0;
double rms_rel = 1e100;
double max_rel = 1e100;
};
static AxisTorusFit fit_axis_torus_like(int axis) {
AxisTorusFit fit;
fit.axis = axis;
if (N < 600 || originalP.empty()) return fit;
double min_t = 1e100, max_t = -1e100;
double min_u = 1e100, max_u = -1e100;
double min_v = 1e100, max_v = -1e100;
for (const Vec3& p : originalP) {
double t, u, v;
axis_components(p, axis, t, u, v);
min_t = min(min_t, t); max_t = max(max_t, t);
min_u = min(min_u, u); max_u = max(max_u, u);
min_v = min(min_v, v); max_v = max(max_v, v);
}
if (!(max_t > min_t) || !(max_u > min_u) || !(max_v > min_v)) return fit;
fit.center_t = 0.5 * (min_t + max_t);
fit.center_u = 0.5 * (min_u + max_u);
fit.center_v = 0.5 * (min_v + max_v);
double min_r = 1e100;
double max_r = 0.0;
for (const Vec3& p : originalP) {
double t, u, v;
axis_components(p, axis, t, u, v);
const double du = u - fit.center_u;
const double dv = v - fit.center_v;
const double r = sqrt(du * du + dv * dv);
min_r = min(min_r, r);
max_r = max(max_r, r);
}
if (!(max_r > min_r) || !(min_r > 1e-10)) return fit;
const double major = 0.5 * (max_r + min_r);
const double minor_r = 0.5 * (max_r - min_r);
const double minor_t = 0.5 * (max_t - min_t);
const double minor = 0.5 * (minor_r + minor_t);
if (!(major > 1e-10) || !(minor > 1e-10)) return fit;
if (major < minor * 1.35) return fit;
if (fabs(minor_r - minor_t) > minor * 0.22) return fit;
const int sample_limit = 240000;
const int stride = max(1, N / sample_limit);
double sum_sq = 0.0;
double max_abs = 0.0;
int sampled = 0;
for (int i = 0; i < N; i += stride) {
double t, u, v;
axis_components(originalP[i], axis, t, u, v);
const double du = u - fit.center_u;
const double dv = v - fit.center_v;
const double rho = sqrt(du * du + dv * dv);
const double tube = sqrt((rho - major) * (rho - major)
+ (t - fit.center_t) * (t - fit.center_t));
const double err = fabs(tube - minor);
sum_sq += err * err;
max_abs = max(max_abs, err);
++sampled;
}
if (sampled < 200) return fit;
fit.major = major;
fit.minor = minor;
fit.rms_rel = sqrt(sum_sq / (double)sampled) / minor;
fit.max_rel = max_abs / minor;
const double rms_limit = (N < 3000 ? 0.018 : 0.012);
const double max_limit = (N < 3000 ? 0.080 : 0.055);
fit.ok = fit.rms_rel <= rms_limit && fit.max_rel <= max_limit;
return fit;
}
static bool build_axis_torus_mesh(const AxisTorusFit& fit, int major_steps, int minor_steps,
vector<Vec3>& verts, vector<Face>& out_faces) {
if (!fit.ok || major_steps < 12 || minor_steps < 6) return false;
const int vertex_count = major_steps * minor_steps;
if (vertex_count > N) return false;
verts.clear();
out_faces.clear();
verts.reserve(vertex_count);
out_faces.reserve(vertex_count * 2);
const double pi = acos(-1.0);
auto torus_normal_at = [&](const Vec3& p) {
double t, u, v;
axis_components(p, fit.axis, t, u, v);
const double du = u - fit.center_u;
const double dv = v - fit.center_v;
const double rho = sqrt(du * du + dv * dv);
if (!(rho > 1e-12) || !(fit.minor > 1e-12)) return Vec3{0.0, 0.0, 0.0};
const double cu = du / rho;
const double sv = dv / rho;
double cp = (rho - fit.major) / fit.minor;
double sp = (t - fit.center_t) / fit.minor;
const double len = sqrt(cp * cp + sp * sp);
if (len > 1e-12) {
cp /= len;
sp /= len;
}
return axis_make_point(fit.axis, sp, cp * cu, cp * sv);
};
double orient_sum = 0.0;
int orient_count = 0;
const int sample_limit = 100000;
const int stride = max(1, M / sample_limit);
for (int fid = 0; fid < M; fid += stride) {
const Face& f = originalFaces[fid];
const Vec3& a = originalP[f.v[0]];
const Vec3& b = originalP[f.v[1]];
const Vec3& c = originalP[f.v[2]];
Vec3 cr = cross3(b - a, c - a);
const double clen = norm3(cr);
if (!(clen > 0.0)) continue;
const Vec3 ctr = (a + b + c) * (1.0 / 3.0);
Vec3 pred = torus_normal_at(ctr);
const double plen = norm3(pred);
if (!(plen > 0.0)) continue;
orient_sum += dot3(cr * (1.0 / clen), pred * (1.0 / plen));
++orient_count;
}
const bool flip_orientation = orient_count > 0 && orient_sum < 0.0;
for (int i = 0; i < major_steps; ++i) {
const double theta = 2.0 * pi * (double)i / (double)major_steps;
const double ct = cos(theta);
const double st = sin(theta);
for (int j = 0; j < minor_steps; ++j) {
const double phi = 2.0 * pi * (double)j / (double)minor_steps;
const double cp = cos(phi);
const double sp = sin(phi);
const double radial = fit.major + fit.minor * cp;
verts.push_back(axis_make_point(fit.axis,
fit.center_t + fit.minor * sp,
fit.center_u + radial * ct,
fit.center_v + radial * st));
}
}
auto id = [&](int i, int j) {
i = (i % major_steps + major_steps) % major_steps;
j = (j % minor_steps + minor_steps) % minor_steps;
return i * minor_steps + j;
};
auto add_face = [&](int a, int b, int c) {
Face f;
f.v[0] = a; f.v[1] = b; f.v[2] = c;
if (flip_orientation) swap(f.v[1], f.v[2]);
out_faces.push_back(f);
};
for (int i = 0; i < major_steps; ++i) {
for (int j = 0; j < minor_steps; ++j) {
const int a = id(i, j);
const int b = id(i + 1, j);
const int c = id(i + 1, j + 1);
const int d = id(i, j + 1);
add_face(a, b, c);
add_face(a, c, d);
}
}
return true;
}
static bool try_axis_torus_remesh() {
if (elapsed_seconds() > 13.8 || N < 600 || N > 30000) return false;
const int safe_vertices = count_output_vertices_estimate();
if (safe_vertices <= 0) return false;
AxisTorusFit best;
for (int axis = 0; axis < 3; ++axis) {
const AxisTorusFit fit = fit_axis_torus_like(axis);
if (fit.ok && (!best.ok || fit.rms_rel < best.rms_rel)) best = fit;
}
if (!best.ok) return false;
struct TorusTrial {
int major_steps;
int minor_steps;
int proxy_res;
double threshold;
double keep_ratio;
};
vector<TorusTrial> trials;
if (N < 900) {
trials.push_back({36, 10, 256, 0.965, 0.98});
} else if (N < 1500) {
trials.push_back({48, 14, 512, 0.930, 0.98});
trials.push_back({56, 16, 256, 0.945, 0.98});
} else if (N < 5000) {
trials.push_back({56, 16, 512, 0.925, 0.98});
trials.push_back({64, 16, 512, 0.940, 0.98});
trials.push_back({72, 18, 256, 0.965, 0.98});
trials.push_back({80, 20, 256, 0.970, 0.95});
} else if (N < 15000) {
trials.push_back({80, 20, 256, 0.950, 0.90});
trials.push_back({96, 24, 256, 0.960, 0.98});
} else {
trials.push_back({104, 24, 256, 0.955, 0.90});
trials.push_back({112, 28, 256, 0.965, 0.98});
}
MeshState safe_state = capture_state();
vector<Vec3> new_vertices;
vector<Face> new_faces;
for (const TorusTrial& trial : trials) {
if (elapsed_seconds() > 16.0) break;
const int target_vertices = trial.major_steps * trial.minor_steps;
if (target_vertices <= 0 || target_vertices >= safe_vertices) continue;
if ((double)target_vertices > (double)safe_vertices * trial.keep_ratio) continue;
if (!build_axis_torus_mesh(best, trial.major_steps, trial.minor_steps,
new_vertices, new_faces)) continue;
restore_state(safe_state);
bool keep = false;
if (install_replacement_mesh(new_vertices, new_faces) && elapsed_seconds() < 16.5) {
const double proxy = visual_proxy_score(trial.proxy_res);
keep = proxy >= trial.threshold;
}
if (keep) return true;
}
restore_state(safe_state);
return false;
}
static AxisRevolveFit fit_axis_revolved_like(int axis) {
AxisRevolveFit fit;
fit.axis = axis;
if (N < 1000 || originalP.empty()) return fit;
double min_t = 1e100, max_t = -1e100;
double min_u = 1e100, max_u = -1e100;
double min_v = 1e100, max_v = -1e100;
for (const Vec3& p : originalP) {
double t, u, v;
axis_components(p, axis, t, u, v);
min_t = min(min_t, t); max_t = max(max_t, t);
min_u = min(min_u, u); max_u = max(max_u, u);
min_v = min(min_v, v); max_v = max(max_v, v);
}
if (!(max_t > min_t)) return fit;
fit.center_u = 0.5 * (min_u + max_u);
fit.center_v = 0.5 * (min_v + max_v);
fit.t0 = min_t;
fit.t1 = max_t;
double max_r = 0.0;
for (const Vec3& p : originalP) {
double t, u, v;
axis_components(p, axis, t, u, v);
const double du = u - fit.center_u;
const double dv = v - fit.center_v;
max_r = max(max_r, sqrt(du * du + dv * dv));
}
if (!(max_r > 1e-10)) return fit;
const int sample_limit = 240000;
const int stride = max(1, N / sample_limit);
const double radial_eps = max_r * 0.055;
double s_t = 0.0, s_r = 0.0, s_tt = 0.0, s_tr = 0.0;
int ring_samples = 0;
int axis_samples = 0;
for (int i = 0; i < N; i += stride) {
double t, u, v;
axis_components(originalP[i], axis, t, u, v);
const double du = u - fit.center_u;
const double dv = v - fit.center_v;
const double r = sqrt(du * du + dv * dv);
if (r <= radial_eps) {
const double near_end = min(fabs(t - min_t), fabs(t - max_t));
if (near_end > (max_t - min_t) * 0.030) return fit;
++axis_samples;
continue;
}
s_t += t;
s_r += r;
s_tt += t * t;
s_tr += t * r;
++ring_samples;
}
if (ring_samples < 200) return fit;
const double den = (double)ring_samples * s_tt - s_t * s_t;
if (fabs(den) < 1e-18) return fit;
const double a = ((double)ring_samples * s_tr - s_t * s_r) / den;
const double b = (s_r - a * s_t) / (double)ring_samples;
double r0 = a * min_t + b;
double r1 = a * max_t + b;
if (r0 < -0.035 * max_r || r1 < -0.035 * max_r) return fit;
if (fabs(r0) < 0.035 * max_r) r0 = 0.0;
if (fabs(r1) < 0.035 * max_r) r1 = 0.0;
if (max(r0, r1) < 0.25 * max_r) return fit;
double sum_sq = 0.0;
double max_abs = 0.0;
int checked = 0;
for (int i = 0; i < N; i += stride) {
double t, u, v;
axis_components(originalP[i], axis, t, u, v);
const double du = u - fit.center_u;
const double dv = v - fit.center_v;
const double r = sqrt(du * du + dv * dv);
if (r <= radial_eps) continue;
const double pred = max(0.0, a * t + b);
const double err = fabs(r - pred);
sum_sq += err * err;
max_abs = max(max_abs, err);
++checked;
}
if (checked < 200) return fit;
const double rms = sqrt(sum_sq / (double)checked);
if (rms > max_r * 0.0060 || max_abs > max_r * 0.030) return fit;
if (axis_samples > ring_samples / 3) return fit;
fit.r0 = r0;
fit.r1 = r1;
fit.residual = rms / max_r;
fit.ok = true;
return fit;
}
static bool build_axis_revolved_mesh(const AxisRevolveFit& fit, int sides,
vector<Vec3>& verts, vector<Face>& out_faces) {
if (!fit.ok || sides < 12) return false;
const double eps = max(fit.r0, fit.r1) * 1e-6;
const bool cone0 = fit.r0 <= eps;
const bool cone1 = fit.r1 <= eps;
if (cone0 && cone1) return false;
const Vec3 center = axis_make_point(fit.axis, 0.5 * (fit.t0 + fit.t1), fit.center_u, fit.center_v);
const double sign = original_orientation_sign(center);
const double pi = acos(-1.0);
verts.clear();
out_faces.clear();
auto make = [&](double t, double r, int j) {
const double th = 2.0 * pi * (double)j / (double)sides;
return axis_make_point(fit.axis, t,
fit.center_u + r * cos(th),
fit.center_v + r * sin(th));
};
auto add_face = [&](int a, int b, int c) {
Face f;
f.v[0] = a; f.v[1] = b; f.v[2] = c;
orient_face_like_original(verts, f, center, sign);
out_faces.push_back(f);
};
if (!cone0 && !cone1) {
const int bottom_center = 0;
const int top_center = 1;
verts.push_back(axis_make_point(fit.axis, fit.t0, fit.center_u, fit.center_v));
verts.push_back(axis_make_point(fit.axis, fit.t1, fit.center_u, fit.center_v));
const int bottom_ring = (int)verts.size();
for (int j = 0; j < sides; ++j) verts.push_back(make(fit.t0, fit.r0, j));
const int top_ring = (int)verts.size();
for (int j = 0; j < sides; ++j) verts.push_back(make(fit.t1, fit.r1, j));
if ((int)verts.size() > N) return false;
for (int j = 0; j < sides; ++j) {
const int j2 = (j + 1) % sides;
const int b0 = bottom_ring + j;
const int b1 = bottom_ring + j2;
const int t0 = top_ring + j;
const int t1 = top_ring + j2;
add_face(b0, b1, t0);
add_face(b1, t1, t0);
add_face(bottom_center, b0, b1);
add_face(top_center, t1, t0);
}
} else {
const bool apex_at_bottom = cone0;
const double apex_t = apex_at_bottom ? fit.t0 : fit.t1;
const double base_t = apex_at_bottom ? fit.t1 : fit.t0;
const double base_r = apex_at_bottom ? fit.r1 : fit.r0;
const int apex = 0;
const int base_center = 1;
verts.push_back(axis_make_point(fit.axis, apex_t, fit.center_u, fit.center_v));
verts.push_back(axis_make_point(fit.axis, base_t, fit.center_u, fit.center_v));
const int ring = (int)verts.size();
for (int j = 0; j < sides; ++j) verts.push_back(make(base_t, base_r, j));
if ((int)verts.size() > N) return false;
for (int j = 0; j < sides; ++j) {
const int j2 = (j + 1) % sides;
add_face(apex, ring + j, ring + j2);
add_face(base_center, ring + j2, ring + j);
}
}
return true;
}
struct AxisCapsuleFit {
bool ok = false;
int axis = 2;
double center_t = 0.0;
double center_u = 0.0;
double center_v = 0.0;
double radius = 0.0;
double half_cylinder = 0.0;
double rms_rel = 1e100;
double max_rel = 1e100;
};
static AxisCapsuleFit fit_axis_capsule_like(int axis) {
AxisCapsuleFit fit;
fit.axis = axis;
if (N < 1500 || originalP.empty()) return fit;
double min_t = 1e100, max_t = -1e100;
double min_u = 1e100, max_u = -1e100;
double min_v = 1e100, max_v = -1e100;
for (const Vec3& p : originalP) {
double t, u, v;
axis_components(p, axis, t, u, v);
min_t = min(min_t, t); max_t = max(max_t, t);
min_u = min(min_u, u); max_u = max(max_u, u);
min_v = min(min_v, v); max_v = max(max_v, v);
}
const double eu = max_u - min_u;
const double ev = max_v - min_v;
const double et = max_t - min_t;
if (!(eu > 1e-12) || !(ev > 1e-12) || !(et > 1e-12)) return fit;
if (min(eu, ev) < max(eu, ev) * 0.985) return fit;
fit.center_t = 0.5 * (min_t + max_t);
fit.center_u = 0.5 * (min_u + max_u);
fit.center_v = 0.5 * (min_v + max_v);
fit.radius = 0.25 * (eu + ev);
fit.half_cylinder = 0.5 * et - fit.radius;
if (!(fit.radius > 1e-12) || fit.half_cylinder <= fit.radius * 0.18) return fit;
const int sample_limit = 240000;
const int stride = max(1, N / sample_limit);
const double rel_limit = (N < 20000 ? 0.018 : 0.014);
const double max_limit = (N < 20000 ? 0.060 : 0.045);
int sampled = 0;
int side_samples = 0;
int cap_samples = 0;
int failed = 0;
double sum_sq = 0.0;
double max_abs_rel = 0.0;
for (int i = 0; i < N; i += stride) {
double t, u, v;
axis_components(originalP[i], axis, t, u, v);
const double dt = t - fit.center_t;
const double du = u - fit.center_u;
const double dv = v - fit.center_v;
const double rho = sqrt(du * du + dv * dv);
const double at = fabs(dt);
double rel = 1e100;
if (at <= fit.half_cylinder) {
rel = fabs(rho - fit.radius) / fit.radius;
if (at < fit.half_cylinder * 0.70) ++side_samples;
} else {
const double q = sqrt(rho * rho + (at - fit.half_cylinder) * (at - fit.half_cylinder));
rel = fabs(q - fit.radius) / fit.radius;
if (at > fit.half_cylinder + fit.radius * 0.20) ++cap_samples;
}
if (at > fit.half_cylinder + fit.radius * 1.015) rel = 1e100;
if (!(rel <= max_limit)) {
++failed;
if (failed > max(2, sampled / 60 + 2)) return fit;
}
if (rel < 1e90) {
sum_sq += rel * rel;
max_abs_rel = max(max_abs_rel, rel);
}
++sampled;
}
if (sampled < 200) return fit;
if (side_samples < max(20, sampled / 25)) return fit;
if (cap_samples < max(20, sampled / 25)) return fit;
fit.rms_rel = sqrt(sum_sq / (double)sampled);
fit.max_rel = max_abs_rel;
fit.ok = fit.rms_rel <= rel_limit && fit.max_rel <= max_limit;
return fit;
}
static int capsule_vertex_count(int sides, int hemi_steps, int cyl_segments) {
if (sides < 8 || hemi_steps < 3 || cyl_segments < 1) return 0;
const int ring_count = 2 * hemi_steps + cyl_segments - 1;
return 2 + ring_count * sides;
}
static bool build_axis_capsule_mesh(const AxisCapsuleFit& fit, int sides, int hemi_steps,
int cyl_segments, vector<Vec3>& verts,
vector<Face>& out_faces) {
if (!fit.ok) return false;
const int ring_count = 2 * hemi_steps + cyl_segments - 1;
const int vertex_count = capsule_vertex_count(sides, hemi_steps, cyl_segments);
if (vertex_count <= 0 || vertex_count > N) return false;
verts.clear();
out_faces.clear();
verts.reserve(vertex_count);
out_faces.reserve(2 * ring_count * sides);
const Vec3 center = axis_make_point(fit.axis, fit.center_t, fit.center_u, fit.center_v);
const double orientation_sign = original_orientation_sign(center);
const double pi = acos(-1.0);
auto make = [&](double rr, double w, int j) {
const double th = 2.0 * pi * (double)j / (double)sides;
return axis_make_point(fit.axis,
fit.center_t + w,
fit.center_u + rr * cos(th),
fit.center_v + rr * sin(th));
};
verts.push_back(axis_make_point(fit.axis, fit.center_t + fit.half_cylinder + fit.radius,
fit.center_u, fit.center_v));
vector<int> ring_start;
ring_start.reserve(ring_count);
auto add_ring = [&](double rr, double w) {
ring_start.push_back((int)verts.size());
for (int j = 0; j < sides; ++j) verts.push_back(make(rr, w, j));
};
for (int i = 1; i <= hemi_steps; ++i) {
const double phi = 0.5 * pi * (double)i / (double)hemi_steps;
add_ring(fit.radius * sin(phi), fit.half_cylinder + fit.radius * cos(phi));
}
for (int i = 1; i <= cyl_segments; ++i) {
const double t = (double)i / (double)cyl_segments;
add_ring(fit.radius, fit.half_cylinder * (1.0 - 2.0 * t));
}
for (int i = 1; i < hemi_steps; ++i) {
const double phi = 0.5 * pi + 0.5 * pi * (double)i / (double)hemi_steps;
add_ring(fit.radius * sin(phi), -fit.half_cylinder + fit.radius * cos(phi));
}
const int bottom = (int)verts.size();
verts.push_back(axis_make_point(fit.axis, fit.center_t - fit.half_cylinder - fit.radius,
fit.center_u, fit.center_v));
if ((int)verts.size() != vertex_count || (int)ring_start.size() != ring_count) return false;
auto ring = [&](int r, int j) {
return ring_start[r] + ((j % sides + sides) % sides);
};
auto add_face = [&](int a, int b, int c) {
Face f;
f.v[0] = a; f.v[1] = b; f.v[2] = c;
orient_face_like_original(verts, f, center, orientation_sign);
out_faces.push_back(f);
};
for (int j = 0; j < sides; ++j) add_face(0, ring(0, j + 1), ring(0, j));
for (int r = 0; r + 1 < ring_count; ++r) {
for (int j = 0; j < sides; ++j) {
const int a = ring(r, j);
const int b = ring(r, j + 1);
const int c = ring(r + 1, j);
const int d = ring(r + 1, j + 1);
add_face(a, b, c);
add_face(b, d, c);
}
}
for (int j = 0; j < sides; ++j) add_face(bottom, ring(ring_count - 1, j), ring(ring_count - 1, j + 1));
return true;
}
static bool try_axis_capsule_remesh() {
if (elapsed_seconds() > 13.8 || N < 1500) return false;
const int safe_vertices = count_output_vertices_estimate();
if (safe_vertices <= 0) return false;
AxisCapsuleFit best;
const int axes[3] = {2, 0, 1};
for (int i = 0; i < 3; ++i) {
const AxisCapsuleFit fit = fit_axis_capsule_like(axes[i]);
if (fit.ok && (!best.ok || fit.rms_rel < best.rms_rel)) best = fit;
}
if (!best.ok) return false;
struct CapsuleTrial {
int sides;
int hemi_steps;
int cyl_segments;
int proxy_res;
double threshold;
double keep_ratio;
};
vector<CapsuleTrial> trials;
if (N < 20000) {
trials.push_back({24, 8, 1, 512, 0.925, 0.92});
trials.push_back({32, 10, 1, 512, 0.930, 0.97});
} else if (N < 100000) {
trials.push_back({28, 8, 1, 256, 0.925, 0.92});
trials.push_back({36, 10, 1, 256, 0.932, 0.97});
} else {
trials.push_back({32, 8, 1, 192, 0.928, 0.92});
trials.push_back({40, 10, 1, 192, 0.935, 0.97});
}
MeshState safe_state = capture_state();
vector<Vec3> new_vertices;
vector<Face> new_faces;
for (const CapsuleTrial& trial : trials) {
if (elapsed_seconds() > 16.0) break;
const int target_vertices = capsule_vertex_count(trial.sides, trial.hemi_steps, trial.cyl_segments);
if (target_vertices <= 0 || target_vertices >= safe_vertices) continue;
if ((double)target_vertices > (double)safe_vertices * trial.keep_ratio) continue;
if (!build_axis_capsule_mesh(best, trial.sides, trial.hemi_steps, trial.cyl_segments,
new_vertices, new_faces)) continue;
restore_state(safe_state);
bool keep = false;
if (install_replacement_mesh(new_vertices, new_faces) && elapsed_seconds() < 16.5) {
const double proxy = visual_proxy_score(trial.proxy_res);
keep = proxy >= trial.threshold;
}
if (keep) return true;
}
restore_state(safe_state);
return false;
}
static bool try_axis_revolved_primitive_remesh() {
if (elapsed_seconds() > 13.8) return false;
const int safe_vertices = count_output_vertices_estimate();
if (safe_vertices <= 0) return false;
for (int axis = 0; axis < 3 && elapsed_seconds() < 13.8; ++axis) {
const AxisRevolveFit fit = fit_axis_revolved_like(axis);
if (!fit.ok) continue;
const bool cone_like = min(fit.r0, fit.r1) <= max(fit.r0, fit.r1) * 0.08;
const int sides = cone_like ? 32 : 32;
vector<Vec3> new_vertices;
vector<Face> new_faces;
if (!build_axis_revolved_mesh(fit, sides, new_vertices, new_faces)) continue;
if ((int)new_vertices.size() >= safe_vertices) continue;
if ((int)new_vertices.size() * 100 > safe_vertices * 90) continue;
MeshState safe_state = capture_state();
bool keep = false;
if (install_replacement_mesh(new_vertices, new_faces) && elapsed_seconds() < 16.5) {
const int proxy_res = (N < 20000 ? 512 : 256);
const double proxy = visual_proxy_score(proxy_res);
keep = proxy >= 0.945;
}
if (keep) return true;
restore_state(safe_state);
}
return false;
}
static bool try_strict_ellipsoid_remesh() {
if (N < 900 || elapsed_seconds() > 13.0) return false;
const EllipsoidFit fit = fit_strict_ellipsoid_like();
if (!fit.ok) return false;
int lat = 20;
int lon = 40;
int proxy_res = 512;
double proxy_threshold = 0.930;
if (N < 1500) {
lat = 20;
lon = 40;
proxy_res = 256;
proxy_threshold = 0.950;
} else if (N < 3000) {
lat = 20;
lon = 40;
proxy_res = 256;
proxy_threshold = 0.905;
} else if (N < 5000) {
lat = 20;
lon = 40;
proxy_res = 256;
proxy_threshold = 0.904;
} else if (N < 30000) {
lat = 20;
lon = 40;
proxy_threshold = 0.930;
} else if (N < 100000) {
lat = 22;
lon = 44;
proxy_res = 256;
proxy_threshold = 0.925;
} else {
lat = 22;
lon = 44;
proxy_res = 256;
proxy_threshold = 0.925;
}
vector<Vec3> new_vertices;
vector<Face> new_faces;
if (!build_uv_ellipsoid_mesh(fit, lat, lon, new_vertices, new_faces)) return false;
const int safe_vertices = count_output_vertices_estimate();
if (safe_vertices <= 0) return false;
if ((int)new_vertices.size() >= safe_vertices) return false;
if ((int)new_vertices.size() * 100 > safe_vertices * 92) return false;
MeshState safe_state = capture_state();
bool keep = false;
if (install_replacement_mesh(new_vertices, new_faces) && elapsed_seconds() < 16.5) {
const double proxy = visual_proxy_score(proxy_res);
keep = proxy >= proxy_threshold;
}
if (!keep) restore_state(safe_state);
return keep;
}
static int count_output_vertices_estimate() {
vector<unsigned char> used(N, 0);
int cnt = 0;
for (int fid = 0; fid < (int)faces.size(); ++fid) {
if (!alive_f[fid]) continue;
const Face& f = faces[fid];
for (int k = 0; k < 3; ++k) {
int v = f.v[k];
if (v >= 0 && v < N && alive_v[v] && !used[v]) {
used[v] = 1;
++cnt;
}
}
}
return cnt;
}
struct StrictSphereFit {
bool ok = false;
Vec3 center{};
double radius = 1.0;
};
static StrictSphereFit fit_strict_sphere_like() {
StrictSphereFit s;
if (N < 900 || originalP.empty()) return s;
Vec3 mn = originalP[0], mx = originalP[0];
for (const Vec3& p : originalP) {
mn.x = min(mn.x, p.x); mn.y = min(mn.y, p.y); mn.z = min(mn.z, p.z);
mx.x = max(mx.x, p.x); mx.y = max(mx.y, p.y); mx.z = max(mx.z, p.z);
}
const double ex = mx.x - mn.x;
const double ey = mx.y - mn.y;
const double ez = mx.z - mn.z;
const double max_extent = max(ex, max(ey, ez));
const double min_extent = min(ex, min(ey, ez));
if (!(max_extent > 1e-12) || min_extent < max_extent * 0.982) return s;
s.center = (mn + mx) * 0.5;
const int sample_limit = 200000;
const int stride = max(1, N / sample_limit);
double sum_r = 0.0;
double sum_r2 = 0.0;
int sampled = 0;
for (int i = 0; i < N; i += stride) {
const double r = norm3(originalP[i] - s.center);
sum_r += r;
sum_r2 += r * r;
++sampled;
}
if (sampled == 0) return s;
const double inv = 1.0 / (double)sampled;
const double mean_r = sum_r * inv;
if (!(mean_r > 1e-12)) return s;
const double std_r = sqrt(max(0.0, sum_r2 * inv - mean_r * mean_r));
double max_abs = 0.0;
for (int i = 0; i < N; i += stride) {
max_abs = max(max_abs, fabs(norm3(originalP[i] - s.center) - mean_r));
}
const double rel_std = std_r / mean_r;
const double rel_max = max_abs / mean_r;
const double std_limit = N < 5000 ? 1.2e-4 : 1.8e-4;
const double max_limit = N < 5000 ? 8.5e-4 : 1.2e-3;
if (rel_std > std_limit || rel_max > max_limit) return s;
s.radius = mean_r;
s.ok = true;
return s;
}
static bool build_uv_sphere_mesh(const StrictSphereFit& sphere, int lat, int lon,
vector<Vec3>& verts, vector<Face>& out_faces) {
const int out_vertices = 2 + (lat - 1) * lon;
if (!sphere.ok || lat < 4 || lon < 8 || out_vertices > N) return false;
verts.clear();
out_faces.clear();
verts.reserve(out_vertices);
out_faces.reserve(2 * lat * lon);
const double pi = acos(-1.0);
const Vec3 c = sphere.center;
const double r = sphere.radius;
const double orientation_sign = original_orientation_sign(c);
verts.push_back({c.x, c.y, c.z + r});
for (int i = 1; i < lat; ++i) {
const double phi = pi * (double)i / (double)lat;
const double z = cos(phi);
const double rr = sin(phi);
for (int j = 0; j < lon; ++j) {
const double th = 2.0 * pi * (double)j / (double)lon;
verts.push_back({c.x + r * rr * cos(th), c.y + r * rr * sin(th), c.z + r * z});
}
}
verts.push_back({c.x, c.y, c.z - r});
auto vid = [&](int ring, int j) {
return 1 + (ring - 1) * lon + ((j % lon + lon) % lon);
};
auto add_face = [&](int a, int b, int cc) {
Face f;
f.v[0] = a; f.v[1] = b; f.v[2] = cc;
orient_face_like_original(verts, f, c, orientation_sign);
out_faces.push_back(f);
};
const int bottom = out_vertices - 1;
for (int j = 0; j < lon; ++j) add_face(0, vid(1, j + 1), vid(1, j));
for (int i = 1; i < lat - 1; ++i) {
for (int j = 0; j < lon; ++j) {
const int a = vid(i, j);
const int b = vid(i, j + 1);
const int cc = vid(i + 1, j);
const int d = vid(i + 1, j + 1);
add_face(a, b, cc);
add_face(b, d, cc);
}
}
for (int j = 0; j < lon; ++j) add_face(bottom, vid(lat - 1, j), vid(lat - 1, j + 1));
return true;
}
static bool try_strict_sphere_remesh() {
if (N < 900 || N >= 50000 || elapsed_seconds() > 13.0) return false;
const StrictSphereFit sphere = fit_strict_sphere_like();
if (!sphere.ok) return false;
const int safe_vertices = count_output_vertices_estimate();
if (safe_vertices <= 0) return false;
struct SphereTrial {
int lat;
int lon;
int proxy_res;
double threshold;
double keep_ratio;
};
vector<SphereTrial> trials;
if (N < 1500) {
trials.push_back({12, 24, 1024, 0.922, 0.94});
trials.push_back({14, 28, 768, 0.928, 0.96});
trials.push_back({16, 32, 512, 0.934, 0.98});
} else if (N < 5000) {
trials.push_back({16, 32, 512, 0.925, 0.92});
trials.push_back({18, 36, 512, 0.935, 0.96});
} else if (N < 20000) {
trials.push_back({18, 36, 512, 0.920, 0.92});
trials.push_back({20, 40, 512, 0.935, 0.98});
} else if (N < 50000) {
trials.push_back({22, 44, 256, 0.920, 0.92});
trials.push_back({24, 48, 256, 0.935, 0.96});
trials.push_back({28, 56, 256, 0.945, 0.98});
}
MeshState safe_state = capture_state();
vector<Vec3> new_vertices;
vector<Face> new_faces;
for (const SphereTrial& trial : trials) {
if (elapsed_seconds() > 16.0) break;
const int target_vertices = 2 + (trial.lat - 1) * trial.lon;
if (target_vertices <= 0 || target_vertices >= safe_vertices) continue;
if ((double)target_vertices > (double)safe_vertices * trial.keep_ratio) continue;
if (!build_uv_sphere_mesh(sphere, trial.lat, trial.lon, new_vertices, new_faces)) continue;
restore_state(safe_state);
bool keep = false;
if (install_replacement_mesh(new_vertices, new_faces) && elapsed_seconds() < 16.5) {
keep = visual_proxy_score(trial.proxy_res) >= trial.threshold;
}
if (keep) return true;
}
restore_state(safe_state);
return false;
}
namespace AggressiveHeapCandidate {
struct Vec3 {
double x = 0.0, y = 0.0, z = 0.0;
};
static inline Vec3 operator+(const Vec3& a, const Vec3& b) {
return {a.x + b.x, a.y + b.y, a.z + b.z};
}
static inline Vec3 operator-(const Vec3& a, const Vec3& b) {
return {a.x - b.x, a.y - b.y, a.z - b.z};
}
static inline Vec3 operator*(const Vec3& a, double k) {
return {a.x * k, a.y * k, a.z * k};
}
static inline double dot3(const Vec3& a, const Vec3& b) {
return a.x * b.x + a.y * b.y + a.z * b.z;
}
static inline Vec3 cross3(const Vec3& a, const Vec3& b) {
return {a.y * b.z - a.z * b.y,
a.z * b.x - a.x * b.z,
a.x * b.y - a.y * b.x};
}
static inline double norm2(const Vec3& a) {
return dot3(a, a);
}
static inline double dist2(const Vec3& a, const Vec3& b) {
return norm2(a - b);
}
static inline bool normalize(Vec3& a) {
double n2 = norm2(a);
if (n2 <= 1e-30) return false;
double inv = 1.0 / sqrt(n2);
a = a * inv;
return true;
}
struct Face {
int v[3];
unsigned char active = 1;
};
struct Quadric {
// Symmetric 4x4 matrix stored as:
// xx, xy, xz, xw, yy, yz, yw, zz, zw, ww
double q[10];
Quadric() { memset(q, 0, sizeof(q)); }
void add_plane(double a, double b, double c, double d, double w = 1.0) {
q[0] += w * a * a;
q[1] += w * a * b;
q[2] += w * a * c;
q[3] += w * a * d;
q[4] += w * b * b;
q[5] += w * b * c;
q[6] += w * b * d;
q[7] += w * c * c;
q[8] += w * c * d;
q[9] += w * d * d;
}
void add(const Quadric& o) {
for (int i = 0; i < 10; ++i) q[i] += o.q[i];
}
double eval(const Vec3& p) const {
double x = p.x, y = p.y, z = p.z;
return q[0] * x * x
+ 2.0 * q[1] * x * y
+ 2.0 * q[2] * x * z
+ 2.0 * q[3] * x
+ q[4] * y * y
+ 2.0 * q[5] * y * z
+ 2.0 * q[6] * y
+ q[7] * z * z
+ 2.0 * q[8] * z
+ q[9];
}
};
struct FastInput {
vector<char> buf;
char* p = nullptr;
void read_all() {
char tmp[1 << 16];
size_t n;
while ((n = fread(tmp, 1, sizeof(tmp), stdin)) > 0) {
buf.insert(buf.end(), tmp, tmp + n);
}
buf.push_back('\0');
p = buf.data();
}
inline void skip_ws() {
while (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\t') ++p;
}
int next_int() {
skip_ws();
int sign = 1;
if (*p == '-') {
sign = -1;
++p;
}
int x = 0;
while (*p >= '0' && *p <= '9') {
x = x * 10 + (*p - '0');
++p;
}
return x * sign;
}
double next_double() {
skip_ws();
char* e = nullptr;
double x = strtod(p, &e);
p = e;
return x;
}
char next_char() {
skip_ws();
return *p++;
}
};
struct EdgeFace {
uint64_t key;
int face;
bool operator<(const EdgeFace& o) const {
if (key != o.key) return key < o.key;
return face < o.face;
}
};
struct Node {
double cost;
int u, v;
int vu, vv;
bool operator<(const Node& o) const {
return cost > o.cost;
}
};
static int N, M;
static vector<Vec3> Orig;
static vector<Vec3> P;
static vector<Face> F;
static vector<vector<int>> incident;
static vector<Quadric> Q;
static vector<unsigned char> active_vertex;
static vector<int> head_member, tail_member, next_member, cluster_size;
static vector<int> version_id;
static vector<int> mark_seen;
static int mark_token = 7;
static int active_vertices = 0;
static int active_faces = 0;
static double hausdorff_limit2 = 0.0;
static double diagonal_len = 0.0;
static int target_vertices = 0;
static double target_ratio_override = -1.0;
static priority_queue<Node> pq;
static chrono::steady_clock::time_point start_time;
static double TIME_LIMIT_SECONDS = 19.25;
static const double MIN_NORMAL_COS = 0.55;
static inline uint64_t edge_key(int a, int b) {
if (a > b) swap(a, b);
return (uint64_t)(uint32_t)a << 32 | (uint32_t)b;
}
static inline int key_a(uint64_t key) {
return (int)(key >> 32);
}
static inline int key_b(uint64_t key) {
return (int)(key & 0xffffffffu);
}
static inline bool time_left() {
auto now = chrono::steady_clock::now();
double elapsed = chrono::duration<double>(now - start_time).count();
return elapsed < TIME_LIMIT_SECONDS;
}
static inline bool face_has_vertex(const Face& f, int v) {
return f.v[0] == v || f.v[1] == v || f.v[2] == v;
}
static inline bool face_has_edge(const Face& f, int a, int b) {
return face_has_vertex(f, a) && face_has_vertex(f, b);
}
static inline int third_vertex(const Face& f, int a, int b) {
for (int k = 0; k < 3; ++k) {
int x = f.v[k];
if (x != a && x != b) return x;
}
return -1;
}
static inline Vec3 face_normal_from_indices(int a, int b, int c) {
return cross3(P[b] - P[a], P[c] - P[a]);
}
static bool are_adjacent(int a, int b) {
const vector<int>& ia = incident[a];
const vector<int>& ib = incident[b];
const vector<int>& small = (ia.size() < ib.size()) ? ia : ib;
for (int fid : small) {
if (!F[fid].active) continue;
if (face_has_edge(F[fid], a, b)) return true;
}
return false;
}
static bool collect_edge_opposites(int a, int b, int opp[2]) {
int cnt = 0;
const vector<int>& ia = incident[a];
const vector<int>& ib = incident[b];
const vector<int>& small = (ia.size() < ib.size()) ? ia : ib;
for (int fid : small) {
if (!F[fid].active) continue;
const Face& f = F[fid];
if (!face_has_edge(f, a, b)) continue;
if (cnt >= 2) return false;
int t = third_vertex(f, a, b);
if (t < 0) return false;
opp[cnt++] = t;
}
if (cnt != 2) return false;
return opp[0] != opp[1];
}
static bool link_condition_ok(int a, int b) {
int opp[2];
if (!collect_edge_opposites(a, b, opp)) return false;
if (mark_token > 2000000000) {
fill(mark_seen.begin(), mark_seen.end(), 0);
mark_token = 7;
}
int token_u = mark_token++;
int token_common = mark_token++;
for (int fid : incident[a]) {
if (!F[fid].active) continue;
const Face& f = F[fid];
if (!face_has_vertex(f, a)) continue;
for (int k = 0; k < 3; ++k) {
int x = f.v[k];
if (x != a && x != b) mark_seen[x] = token_u;
}
}
int common_count = 0;
int got0 = 0, got1 = 0;
for (int fid : incident[b]) {
if (!F[fid].active) continue;
const Face& f = F[fid];
if (!face_has_vertex(f, b)) continue;
for (int k = 0; k < 3; ++k) {
int x = f.v[k];
if (x == a || x == b) continue;
if (mark_seen[x] == token_u) {
mark_seen[x] = token_common;
++common_count;
if (x == opp[0]) got0 = 1;
if (x == opp[1]) got1 = 1;
if (common_count > 2) return false;
}
}
}
return common_count == 2 && got0 && got1;
}
static bool cluster_can_move_to(int v, const Vec3& to) {
for (int m = head_member[v]; m != -1; m = next_member[m]) {
if (dist2(Orig[m], to) > hausdorff_limit2) return false;
}
return true;
}
static bool normals_ok_after_collapse(int a, int b, const Vec3& to) {
auto scan = [&](int src) -> bool {
for (int fid : incident[src]) {
if (!F[fid].active) continue;
const Face& f = F[fid];
bool has_a = face_has_vertex(f, a);
bool has_b = face_has_vertex(f, b);
if (!has_a && !has_b) continue;
if (has_a && has_b) continue;
Vec3 old_p[3] = {P[f.v[0]], P[f.v[1]], P[f.v[2]]};
Vec3 new_p[3] = {old_p[0], old_p[1], old_p[2]};
for (int k = 0; k < 3; ++k) {
if (f.v[k] == a || f.v[k] == b) new_p[k] = to;
}
Vec3 old_n = cross3(old_p[1] - old_p[0], old_p[2] - old_p[0]);
Vec3 new_n = cross3(new_p[1] - new_p[0], new_p[2] - new_p[0]);
double old_len2 = norm2(old_n);
double new_len2 = norm2(new_n);
if (old_len2 <= 1e-30 || new_len2 <= 1e-30) return false;
double d = dot3(old_n, new_n);
double limit = MIN_NORMAL_COS * sqrt(old_len2 * new_len2);
if (d <= limit) return false;
}
return true;
};
return scan(a) && scan(b);
}
static bool solve_optimal_position(const Quadric& q, Vec3& out) {
double a00 = q.q[0], a01 = q.q[1], a02 = q.q[2];
double a10 = q.q[1], a11 = q.q[4], a12 = q.q[5];
double a20 = q.q[2], a21 = q.q[5], a22 = q.q[7];
double b0 = -q.q[3], b1 = -q.q[6], b2 = -q.q[8];
double det = a00 * (a11 * a22 - a12 * a21)
- a01 * (a10 * a22 - a12 * a20)
+ a02 * (a10 * a21 - a11 * a20);
if (fabs(det) < 1e-14) return false;
double dx = b0 * (a11 * a22 - a12 * a21)
- a01 * (b1 * a22 - a12 * b2)
+ a02 * (b1 * a21 - a11 * b2);
double dy = a00 * (b1 * a22 - a12 * b2)
- b0 * (a10 * a22 - a12 * a20)
+ a02 * (a10 * b2 - b1 * a20);
double dz = a00 * (a11 * b2 - b1 * a21)
- a01 * (a10 * b2 - b1 * a20)
+ b0 * (a10 * a21 - a11 * a20);
out = {dx / det, dy / det, dz / det};
return isfinite(out.x) && isfinite(out.y) && isfinite(out.z);
}
static bool candidate_ok(int a, int b, const Vec3& pos) {
if (!cluster_can_move_to(a, pos)) return false;
if (!cluster_can_move_to(b, pos)) return false;
return normals_ok_after_collapse(a, b, pos);
}
static bool best_collapse_position(int a, int b, Vec3& best_pos, double& best_cost) {
Quadric q = Q[a];
q.add(Q[b]);
Vec3 opt;
Vec3 cand[6];
int cnt = 0;
if (solve_optimal_position(q, opt)) cand[cnt++] = opt;
cand[cnt++] = (P[a] + P[b]) * 0.5;
cand[cnt++] = P[a];
cand[cnt++] = P[b];
cand[cnt++] = P[a] * 0.75 + P[b] * 0.25;
cand[cnt++] = P[a] * 0.25 + P[b] * 0.75;
best_cost = 1e100;
bool ok = false;
for (int i = 0; i < cnt; ++i) {
const Vec3& pos = cand[i];
if (!candidate_ok(a, b, pos)) continue;
double c = q.eval(pos) + 0.0003 * (dist2(pos, P[a]) + dist2(pos, P[b]));
if (c < best_cost) {
best_cost = c;
best_pos = pos;
ok = true;
}
}
return ok;
}
static double cheap_edge_cost(int a, int b) {
Quadric q = Q[a];
q.add(Q[b]);
Vec3 opt;
double best = 1e100;
if (solve_optimal_position(q, opt)) best = min(best, q.eval(opt));
best = min(best, q.eval((P[a] + P[b]) * 0.5));
best = min(best, q.eval(P[a]));
best = min(best, q.eval(P[b]));
return best + 0.0003 * dist2(P[a], P[b]);
}
static void push_edge(int a, int b) {
if (a == b) return;
if (!active_vertex[a] || !active_vertex[b]) return;
double d2 = dist2(P[a], P[b]);
if (d2 > 4.00001 * hausdorff_limit2) return;
double c = cheap_edge_cost(a, b);
pq.push({c, a, b, version_id[a], version_id[b]});
}
static void compact_incident(int v) {
vector<int>& ids = incident[v];
if (ids.size() < 128) return;
size_t alive = 0;
for (int fid : ids) {
if (F[fid].active && face_has_vertex(F[fid], v)) ++alive;
}
if (alive * 3 + 32 >= ids.size()) return;
vector<int> keep;
keep.reserve(alive + 8);
for (int fid : ids) {
if (F[fid].active && face_has_vertex(F[fid], v)) keep.push_back(fid);
}
ids.swap(keep);
}
static void merge_members(int src, int dst) {
if (head_member[src] == -1) return;
next_member[tail_member[dst]] = head_member[src];
tail_member[dst] = tail_member[src];
cluster_size[dst] += cluster_size[src];
head_member[src] = tail_member[src] = -1;
cluster_size[src] = 0;
}
static void do_collapse(int src, int dst, const Vec3& pos) {
Q[dst].add(Q[src]);
P[dst] = pos;
for (int fid : incident[src]) {
if (!F[fid].active) continue;
Face& f = F[fid];
bool has_src = false, has_dst = false;
for (int k = 0; k < 3; ++k) {
if (f.v[k] == src) has_src = true;
if (f.v[k] == dst) has_dst = true;
}
if (!has_src) continue;
if (has_dst) {
f.active = 0;
--active_faces;
} else {
for (int k = 0; k < 3; ++k) {
if (f.v[k] == src) f.v[k] = dst;
}
incident[dst].push_back(fid);
}
}
active_vertex[src] = 0;
--active_vertices;
++version_id[src];
++version_id[dst];
merge_members(src, dst);
compact_incident(src);
compact_incident(dst);
for (int fid : incident[dst]) {
if (!F[fid].active) continue;
const Face& f = F[fid];
if (!face_has_vertex(f, dst)) continue;
push_edge(f.v[0], f.v[1]);
push_edge(f.v[1], f.v[2]);
push_edge(f.v[2], f.v[0]);
}
}
static bool attempt_collapse(int a, int b) {
if (a == b) return false;
if (!active_vertex[a] || !active_vertex[b]) return false;
if (dist2(P[a], P[b]) > 4.00001 * hausdorff_limit2) return false;
if (!link_condition_ok(a, b)) return false;
Vec3 pos;
double cost = 0.0;
if (!best_collapse_position(a, b, pos, cost)) return false;
int src, dst;
size_t wa = incident[a].size() + (size_t)cluster_size[a] * 2;
size_t wb = incident[b].size() + (size_t)cluster_size[b] * 2;
if (wa <= wb) {
src = a;
dst = b;
} else {
src = b;
dst = a;
}
do_collapse(src, dst, pos);
return true;
}
static void load_mesh() {
FastInput in;
in.read_all();
N = in.next_int();
M = in.next_int();
P.resize(N);
Orig.resize(N);
F.resize(M);
Vec3 mn{1e100, 1e100, 1e100};
Vec3 mx{-1e100, -1e100, -1e100};
for (int i = 0; i < N; ++i) {
(void)in.next_char();
P[i].x = in.next_double();
P[i].y = in.next_double();
P[i].z = in.next_double();
Orig[i] = P[i];
mn.x = min(mn.x, P[i].x);
mn.y = min(mn.y, P[i].y);
mn.z = min(mn.z, P[i].z);
mx.x = max(mx.x, P[i].x);
mx.y = max(mx.y, P[i].y);
mx.z = max(mx.z, P[i].z);
}
vector<int> deg(N, 0);
for (int i = 0; i < M; ++i) {
(void)in.next_char();
int a = in.next_int() - 1;
int b = in.next_int() - 1;
int c = in.next_int() - 1;
F[i].v[0] = a;
F[i].v[1] = b;
F[i].v[2] = c;
++deg[a];
++deg[b];
++deg[c];
}
Vec3 d = mx - mn;
diagonal_len = sqrt(norm2(d));
double hausdorff_limit = 0.05 * diagonal_len * 0.999999;
hausdorff_limit2 = hausdorff_limit * hausdorff_limit;
incident.resize(N);
for (int i = 0; i < N; ++i) incident[i].reserve(deg[i] + 8);
for (int i = 0; i < M; ++i) {
incident[F[i].v[0]].push_back(i);
incident[F[i].v[1]].push_back(i);
incident[F[i].v[2]].push_back(i);
}
Q.assign(N, Quadric());
active_vertex.assign(N, 1);
head_member.resize(N);
tail_member.resize(N);
next_member.assign(N, -1);
cluster_size.assign(N, 1);
version_id.assign(N, 0);
mark_seen.assign(N, 0);
for (int i = 0; i < N; ++i) {
head_member[i] = tail_member[i] = i;
}
active_vertices = N;
active_faces = M;
}
static double initialize_quadrics_and_edges() {
vector<Vec3> face_normals(M);
vector<EdgeFace> edges;
edges.reserve((size_t)M * 3);
for (int i = 0; i < M; ++i) {
Face& f = F[i];
Vec3 n = face_normal_from_indices(f.v[0], f.v[1], f.v[2]);
if (!normalize(n)) n = {0.0, 0.0, 0.0};
face_normals[i] = n;
double dd = -dot3(n, P[f.v[0]]);
Q[f.v[0]].add_plane(n.x, n.y, n.z, dd);
Q[f.v[1]].add_plane(n.x, n.y, n.z, dd);
Q[f.v[2]].add_plane(n.x, n.y, n.z, dd);
edges.push_back({edge_key(f.v[0], f.v[1]), i});
edges.push_back({edge_key(f.v[1], f.v[2]), i});
edges.push_back({edge_key(f.v[2], f.v[0]), i});
}
sort(edges.begin(), edges.end());
long long unique_edges = 0;
long long feature_edges = 0;
const double feature_cos = cos(35.0 * acos(-1.0) / 180.0);
for (size_t i = 0; i < edges.size();) {
size_t j = i + 1;
while (j < edges.size() && edges[j].key == edges[i].key) ++j;
++unique_edges;
if (j - i == 2) {
double d = dot3(face_normals[edges[i].face], face_normals[edges[i + 1].face]);
if (d < feature_cos) ++feature_edges;
}
int a = key_a(edges[i].key);
int b = key_b(edges[i].key);
push_edge(a, b);
i = j;
}
vector<EdgeFace>().swap(edges);
vector<Vec3>().swap(face_normals);
if (unique_edges == 0) return 0.0;
return (double)feature_edges / (double)unique_edges;
}
static void choose_target(double feature_ratio) {
double ratio = 0.089 + 0.035 * min(0.22, feature_ratio);
if (N <= 8000) ratio = max(ratio, 0.095);
if (N <= 1000) ratio = max(ratio, 0.160);
if (N <= 50) ratio = 0.01;
if (target_ratio_override > 0.0) {
ratio = max(0.040, min(0.500, target_ratio_override));
} else {
ratio = max(0.086, min(0.115, ratio));
}
target_vertices = max(4, (int)ceil((double)N * ratio));
}
static void simplify_mesh() {
double feature_ratio = initialize_quadrics_and_edges();
choose_target(feature_ratio);
long long pops = 0;
while (active_vertices > target_vertices && !pq.empty()) {
if ((++pops & 4095) == 0 && !time_left()) break;
Node cur = pq.top();
pq.pop();
int a = cur.u, b = cur.v;
if (a == b) continue;
if (!active_vertex[a] || !active_vertex[b]) continue;
if (cur.vu != version_id[a] || cur.vv != version_id[b]) {
if (are_adjacent(a, b)) push_edge(a, b);
continue;
}
if (!attempt_collapse(a, b)) continue;
}
}
static void save_mesh() {
vector<int> remap(N, -1);
int out_vertices = 0;
for (int i = 0; i < N; ++i) {
if (active_vertex[i]) remap[i] = out_vertices++;
}
int out_faces = 0;
for (int i = 0; i < M; ++i) {
if (!F[i].active) continue;
int a = F[i].v[0], b = F[i].v[1], c = F[i].v[2];
if (a == b || b == c || c == a) continue;
if (remap[a] < 0 || remap[b] < 0 || remap[c] < 0) continue;
++out_faces;
}
static char outbuf[1 << 20];
setvbuf(stdout, outbuf, _IOFBF, sizeof(outbuf));
printf("%d %d\n", out_vertices, out_faces);
for (int i = 0; i < N; ++i) {
if (!active_vertex[i]) continue;
printf("v %.15g %.15g %.15g\n", P[i].x, P[i].y, P[i].z);
}
for (int i = 0; i < M; ++i) {
if (!F[i].active) continue;
int a = F[i].v[0], b = F[i].v[1], c = F[i].v[2];
if (a == b || b == c || c == a) continue;
if (remap[a] < 0 || remap[b] < 0 || remap[c] < 0) continue;
printf("f %d %d %d\n", remap[a] + 1, remap[b] + 1, remap[c] + 1);
}
}
static bool build_from_parent(const vector<::Vec3>& source_vertices,
const vector<::Face>& source_faces,
double parent_diag,
vector<::Vec3>& out_vertices,
vector<::Face>& out_faces,
double seconds_limit,
double target_ratio = -1.0) {
start_time = chrono::steady_clock::now();
TIME_LIMIT_SECONDS = max(0.25, seconds_limit);
target_ratio_override = target_ratio;
N = (int)source_vertices.size();
M = (int)source_faces.size();
if (N < 4 || M < 4) return false;
Orig.assign(N, Vec3{});
P.assign(N, Vec3{});
F.assign(M, Face{});
Vec3 mn{1e100, 1e100, 1e100};
Vec3 mx{-1e100, -1e100, -1e100};
for (int i = 0; i < N; ++i) {
P[i] = {source_vertices[i].x, source_vertices[i].y, source_vertices[i].z};
Orig[i] = P[i];
mn.x = min(mn.x, P[i].x);
mn.y = min(mn.y, P[i].y);
mn.z = min(mn.z, P[i].z);
mx.x = max(mx.x, P[i].x);
mx.y = max(mx.y, P[i].y);
mx.z = max(mx.z, P[i].z);
}
vector<int> deg(N, 0);
for (int i = 0; i < M; ++i) {
F[i].v[0] = source_faces[i].v[0];
F[i].v[1] = source_faces[i].v[1];
F[i].v[2] = source_faces[i].v[2];
F[i].active = 1;
if (F[i].v[0] < 0 || F[i].v[0] >= N || F[i].v[1] < 0 || F[i].v[1] >= N || F[i].v[2] < 0 || F[i].v[2] >= N) return false;
++deg[F[i].v[0]];
++deg[F[i].v[1]];
++deg[F[i].v[2]];
}
Vec3 d = mx - mn;
diagonal_len = parent_diag > 0.0 ? parent_diag : sqrt(norm2(d));
double hausdorff_limit = 0.05 * diagonal_len * 0.999999;
hausdorff_limit2 = hausdorff_limit * hausdorff_limit;
incident.assign(N, {});
for (int i = 0; i < N; ++i) incident[i].reserve(deg[i] + 8);
for (int i = 0; i < M; ++i) {
incident[F[i].v[0]].push_back(i);
incident[F[i].v[1]].push_back(i);
incident[F[i].v[2]].push_back(i);
}
Q.assign(N, Quadric());
active_vertex.assign(N, 1);
head_member.resize(N);
tail_member.resize(N);
next_member.assign(N, -1);
cluster_size.assign(N, 1);
version_id.assign(N, 0);
mark_seen.assign(N, 0);
mark_token = 7;
for (int i = 0; i < N; ++i) head_member[i] = tail_member[i] = i;
active_vertices = N;
active_faces = M;
pq = priority_queue<Node>();
simplify_mesh();
if (active_vertices <= 0 || active_faces <= 0) return false;
vector<int> remap(N, -1);
out_vertices.clear();
out_faces.clear();
out_vertices.reserve(active_vertices);
for (int i = 0; i < N; ++i) {
if (!active_vertex[i]) continue;
remap[i] = (int)out_vertices.size();
out_vertices.push_back(::Vec3{P[i].x, P[i].y, P[i].z});
}
for (int i = 0; i < M; ++i) {
if (!F[i].active) continue;
int a = F[i].v[0], b = F[i].v[1], c = F[i].v[2];
if (a == b || b == c || c == a) continue;
if (a < 0 || b < 0 || c < 0 || a >= N || b >= N || c >= N) continue;
if (remap[a] < 0 || remap[b] < 0 || remap[c] < 0) continue;
::Face f;
f.v[0] = remap[a];
f.v[1] = remap[b];
f.v[2] = remap[c];
out_faces.push_back(f);
}
return !out_vertices.empty() && !out_faces.empty() && out_vertices.size() < source_vertices.size();
}
} // namespace AggressiveHeapCandidate
static void simplify_mesh() {
simplify_start = chrono::steady_clock::now();
const double d = mesh_diag;
const double min_area2 = max(1e-24, 1e-18 * d * d);
const double exact_plane_tol = max(1e-11, 1e-9 * d);
const double exact_cos = 1.0 - 1e-10;
const bool smooth_lift_allowed = allow_high_density_smooth_lift();
const bool medium_detail_profile = smooth_stats_valid
&& N >= 8000 && N < 50000
&& smooth10_ratio_cache >= 0.900
&& smooth30_ratio_cache >= 0.990
&& sharp22_ratio_cache >= 0.003
&& sharp22_ratio_cache <= 0.050
&& sharp45_ratio_cache <= 0.010
&& bad_ratio_cache <= 0.005;
const bool large_detail_profile = smooth_stats_valid
&& N >= 50000
&& smooth30_ratio_cache >= 0.900
&& sharp45_ratio_cache <= 0.050
&& sharp22_ratio_cache >= 0.020
&& sharp22_ratio_cache <= 0.160
&& bad_ratio_cache <= 0.005;
const bool small_coarse_smooth_profile = smooth_stats_valid
&& N >= 400 && N < 8000
&& smooth30_ratio_cache >= 0.995
&& sharp45_ratio_cache <= 0.006
&& bad_ratio_cache <= 0.005
&& (smooth10_ratio_cache <= 0.760
|| sharp22_ratio_cache <= 0.006
|| sharp22_ratio_cache >= 0.050);
const bool small_detail_smooth_profile = smooth_stats_valid
&& N >= 1000 && N < 8000
&& !small_coarse_smooth_profile
&& smooth10_ratio_cache >= 0.800
&& smooth30_ratio_cache >= 0.980
&& sharp22_ratio_cache >= 0.005
&& sharp22_ratio_cache <= 0.060
&& sharp45_ratio_cache <= 0.012
&& bad_ratio_cache <= 0.005;
const bool small_proxy_chase_profile = small_coarse_smooth_profile || small_detail_smooth_profile;
const bool attempt_detail_profile = medium_detail_profile || large_detail_profile;
const bool small_detail_aggressive_profile = small_detail_smooth_profile
&& N >= 1500 && N < 8000
&& sharp22_ratio_cache >= 0.010
&& sharp45_ratio_cache <= 0.010;
const bool tiny_bumpy_chase_profile = smooth_stats_valid
&& N >= 900 && N < 1500
&& !small_coarse_smooth_profile
&& bad_ratio_cache <= 0.006
&& smooth30_ratio_cache >= 0.820
&& sharp22_ratio_cache >= 0.010
&& sharp22_ratio_cache <= 0.160
&& sharp45_ratio_cache <= 0.070;
const bool certified_boost_profile = smooth_stats_valid
&& N >= 60000 && N <= 260000
&& bad_ratio_cache <= 0.006
&& smooth10_ratio_cache <= 0.820
&& smooth30_ratio_cache >= 0.860
&& smooth30_ratio_cache <= 0.980
&& sharp22_ratio_cache >= 0.070
&& sharp22_ratio_cache <= 0.180
&& sharp45_ratio_cache >= 0.015
&& sharp45_ratio_cache <= 0.070;
const bool attempt_aggressive_profile = (N >= 200000)
|| smooth_lift_allowed
|| attempt_detail_profile
|| small_detail_aggressive_profile;
const bool attempt_relocation_profile = smooth_lift_allowed && N >= 3000 && N <= 50000;
MeshState initial_state;
if (attempt_aggressive_profile || small_proxy_chase_profile || tiny_bumpy_chase_profile || certified_boost_profile) initial_state = capture_state();
const bool small_visual_restart_profile = N >= 300 && N < 3000;
MeshState small_visual_restart_state;
if (small_visual_restart_profile) small_visual_restart_state = capture_state();
double final_radius = 0.020 * d;
double final_plane = 0.0058 * d;
double final_angle_deg = 4.4;
if (N >= 200000) {
final_radius = 0.027 * d;
final_plane = 0.0082 * d;
final_angle_deg = 6.1;
} else if (N >= 30000) {
final_radius = 0.024 * d;
final_plane = 0.0072 * d;
final_angle_deg = 5.3;
} else if (N <= 1000) {
final_radius = 0.014 * d;
final_plane = 0.0036 * d;
final_angle_deg = 3.0;
}
vector<PassParams> passes;
auto add_pass = [&](double radius, double plane, double angle_deg, bool allow_relocate = false) {
PassParams p;
p.radius_limit = radius;
p.plane_tol = plane;
p.cos_limit = cos(angle_deg * acos(-1.0) / 180.0);
p.exact_plane_tol = exact_plane_tol;
p.exact_cos_limit = exact_cos;
p.min_area2 = min_area2;
p.allow_relocate = allow_relocate;
passes.push_back(p);
};
add_pass(0.004 * d, 0.0010 * d, 1.0);
add_pass(0.008 * d, 0.0020 * d, 1.8);
add_pass(0.012 * d, 0.0032 * d, 2.8);
add_pass(final_radius, final_plane, final_angle_deg);
add_pass(final_radius, final_plane, final_angle_deg);
constexpr double WORK_SECONDS = 16.5;
auto run_pass_list = [&](const vector<PassParams>& pass_list, bool allow_early_break) -> bool {
for (const PassParams& params : pass_list) {
int reduced = 0;
for (int fid = 0; fid < (int)faces.size(); ++fid) {
if ((fid & 4095) == 0 && elapsed_seconds() > WORK_SECONDS) return false;
if (!alive_f[fid]) continue;
Face f = faces[fid];
if (!alive_v[f.v[0]] || !alive_v[f.v[1]] || !alive_v[f.v[2]]) continue;
struct EdgeCandidate {
double len2;
int a;
int b;
};
EdgeCandidate e[3] = {
{norm2(P[f.v[0]] - P[f.v[1]]), f.v[0], f.v[1]},
{norm2(P[f.v[1]] - P[f.v[2]]), f.v[1], f.v[2]},
{norm2(P[f.v[2]] - P[f.v[0]]), f.v[2], f.v[0]},
};
sort(e, e + 3, [](const EdgeCandidate& lhs, const EdgeCandidate& rhs) {
return lhs.len2 < rhs.len2;
});
for (const EdgeCandidate& cand : e) {
if (try_collapse_edge(cand.a, cand.b, params)) {
++reduced;
break;
}
}
}
if (allow_early_break && reduced == 0 && params.radius_limit >= final_radius * 0.999) break;
}
return true;
};
if (!run_pass_list(passes, true)) return;
if (small_proxy_chase_profile && elapsed_seconds() < 12.4) {
MeshState safe_state = capture_state();
const int safe_vertices = count_output_vertices_estimate();
restore_state(initial_state);
vector<PassParams> chase_passes;
auto add_chase_pass = [&](double radius, double plane, double angle_deg, bool allow_relocate = false) {
PassParams p;
p.radius_limit = radius;
p.plane_tol = plane;
p.cos_limit = cos(angle_deg * acos(-1.0) / 180.0);
p.exact_plane_tol = exact_plane_tol;
p.exact_cos_limit = exact_cos;
p.min_area2 = min_area2;
p.allow_relocate = allow_relocate;
chase_passes.push_back(p);
};
double chase_radius = 0.026 * d;
double chase_plane = 0.0078 * d;
double chase_angle = 5.8;
double keep_ratio_percent = 92.0;
double proxy_threshold = 0.970;
int chase_proxy_res = 512;
if (small_coarse_smooth_profile) {
if (N <= 1000 && sharp22_ratio_cache >= 0.050) {
chase_radius = 0.050 * d;
chase_plane = 0.0150 * d;
chase_angle = 10.0;
keep_ratio_percent = 96.0;
proxy_threshold = 0.970;
} else if (N < 3000 && sharp22_ratio_cache <= 0.006) {
chase_radius = 0.034 * d;
chase_plane = 0.0102 * d;
chase_angle = 7.7;
keep_ratio_percent = 92.0;
proxy_threshold = 0.955;
} else {
chase_radius = 0.030 * d;
chase_plane = 0.0090 * d;
chase_angle = 6.8;
keep_ratio_percent = 92.0;
proxy_threshold = 0.965;
}
} else if (small_detail_smooth_profile) {
if (N < 3000) {
chase_radius = 0.034 * d;
chase_plane = 0.0102 * d;
chase_angle = 7.7;
keep_ratio_percent = 96.0;
proxy_threshold = 0.952;
chase_proxy_res = (N < 1500 ? 768 : 512);
} else {
chase_radius = 0.024 * d;
chase_plane = 0.0072 * d;
chase_angle = 5.3;
keep_ratio_percent = 92.0;
proxy_threshold = 0.965;
}
}
add_chase_pass(0.004 * d, 0.0010 * d, 1.0);
add_chase_pass(0.008 * d, 0.0020 * d, 1.8);
add_chase_pass(0.012 * d, 0.0032 * d, 2.8);
add_chase_pass(chase_radius, chase_plane, chase_angle);
add_chase_pass(chase_radius, chase_plane, chase_angle);
if (small_coarse_smooth_profile && N >= 1500) {
add_chase_pass(chase_radius, chase_plane, chase_angle);
}
bool keep_chase = false;
if (run_pass_list(chase_passes, true) && elapsed_seconds() < 16.2) {
const int chase_vertices = count_output_vertices_estimate();
if (chase_vertices > 0
&& safe_vertices > 0
&& (double)chase_vertices * 100.0 <= (double)safe_vertices * keep_ratio_percent) {
const double proxy = visual_proxy_score(chase_proxy_res);
keep_chase = proxy >= proxy_threshold;
}
}
if (!keep_chase) restore_state(safe_state);
}
if (N >= 300 && N <= 6000 && elapsed_seconds() < 10.8) {
struct SmallVisualTrial {
double radius;
double plane;
double angle_deg;
double proxy_threshold;
double required_percent;
int repeats;
bool allow_relocate;
};
vector<SmallVisualTrial> visual_trials;
const bool visually_smooth = smooth_stats_valid
&& smooth30_ratio_cache >= 0.960
&& sharp45_ratio_cache <= 0.035
&& bad_ratio_cache <= 0.010;
const bool very_safe_smooth = smooth_stats_valid
&& smooth30_ratio_cache >= 0.990
&& sharp45_ratio_cache <= 0.012
&& bad_ratio_cache <= 0.006;
if (N < 1000) {
visual_trials.push_back({0.026 * d, 0.0078 * d, 5.8, 0.975, 99.0, 1, false});
visual_trials.push_back({0.038 * d, 0.0114 * d, 8.4, visually_smooth ? 0.950 : 0.965, 99.0, 2, visually_smooth});
visual_trials.push_back({0.050 * d, 0.0150 * d, 10.8, very_safe_smooth ? 0.935 : 0.955, 98.0, 2, very_safe_smooth});
} else if (N < 3000) {
visual_trials.push_back({0.028 * d, 0.0084 * d, 6.2, 0.970, 99.0, 1, false});
visual_trials.push_back({0.038 * d, 0.0114 * d, 8.4, visually_smooth ? 0.945 : 0.960, 98.8, 2, visually_smooth});
visual_trials.push_back({0.049 * d, 0.0147 * d, 10.8, very_safe_smooth ? 0.930 : 0.950, 98.0, 2, very_safe_smooth});
} else {
visual_trials.push_back({0.030 * d, 0.0090 * d, 6.8, 0.965, 99.0, 1, false});
visual_trials.push_back({0.041 * d, 0.0123 * d, 9.2, visually_smooth ? 0.940 : 0.955, 98.5, 1, visually_smooth});
visual_trials.push_back({0.050 * d, 0.0150 * d, 11.0, very_safe_smooth ? 0.925 : 0.945, 97.8, 1, very_safe_smooth});
}
for (const SmallVisualTrial& trial : visual_trials) {
if (elapsed_seconds() > 14.2) break;
MeshState trial_safe_state = capture_state();
const int before_vertices = count_output_vertices_estimate();
if (before_vertices <= 0) break;
vector<PassParams> trial_passes;
PassParams p;
p.radius_limit = trial.radius;
p.plane_tol = trial.plane;
p.cos_limit = cos(trial.angle_deg * acos(-1.0) / 180.0);
p.exact_plane_tol = exact_plane_tol;
p.exact_cos_limit = exact_cos;
p.min_area2 = min_area2;
p.allow_relocate = trial.allow_relocate;
for (int i = 0; i < trial.repeats; ++i) trial_passes.push_back(p);
bool keep_visual_trial = false;
if (run_pass_list(trial_passes, true) && elapsed_seconds() < 16.2) {
const int after_vertices = count_output_vertices_estimate();
if (after_vertices > 0
&& (double)after_vertices * 100.0 <= (double)before_vertices * trial.required_percent) {
const double proxy = visual_proxy_score(512);
keep_visual_trial = proxy >= trial.proxy_threshold;
}
}
if (!keep_visual_trial) {
restore_state(trial_safe_state);
if (trial.proxy_threshold <= 0.940) break;
}
}
}
if (small_detail_aggressive_profile && elapsed_seconds() < 12.2) {
MeshState safe_state = capture_state();
const int safe_vertices = count_output_vertices_estimate();
struct DetailChaseTrial {
double small_radius;
double small_plane;
double small_angle;
double big_radius;
double big_plane;
double big_angle;
int repeats;
int proxy_res;
double proxy_threshold;
double required_percent;
};
vector<DetailChaseTrial> detail_trials;
detail_trials.push_back({0.030 * d, 0.0090 * d, 6.8,
0.032 * d, 0.0096 * d, 7.2,
3, 1024, 0.900, 98.0});
detail_trials.push_back({0.024 * d, 0.0072 * d, 5.8,
0.026 * d, 0.0078 * d, 6.2,
2, 512, 0.905, 92.0});
bool accepted_detail_chase = false;
for (const DetailChaseTrial& trial : detail_trials) {
if (elapsed_seconds() > 13.4) break;
restore_state(initial_state);
vector<PassParams> detail_chase_passes;
auto add_detail_chase_pass = [&](double radius, double plane, double angle_deg) {
PassParams p;
p.radius_limit = radius;
p.plane_tol = plane;
p.cos_limit = cos(angle_deg * acos(-1.0) / 180.0);
p.exact_plane_tol = exact_plane_tol;
p.exact_cos_limit = exact_cos;
p.min_area2 = min_area2;
detail_chase_passes.push_back(p);
};
double detail_chase_radius = trial.small_radius;
double detail_chase_plane = trial.small_plane;
double detail_chase_angle = trial.small_angle;
if (N >= 3000) {
detail_chase_radius = trial.big_radius;
detail_chase_plane = trial.big_plane;
detail_chase_angle = trial.big_angle;
}
add_detail_chase_pass(0.004 * d, 0.0010 * d, 1.0);
add_detail_chase_pass(0.008 * d, 0.0020 * d, 1.8);
add_detail_chase_pass(0.012 * d, 0.0032 * d, 2.8);
for (int i = 0; i < trial.repeats; ++i) {
add_detail_chase_pass(detail_chase_radius, detail_chase_plane, detail_chase_angle);
}
bool keep_detail_chase = false;
if (safe_vertices > 0
&& run_pass_list(detail_chase_passes, true)
&& elapsed_seconds() < 16.1) {
const int trial_vertices = count_output_vertices_estimate();
if (trial_vertices > 0
&& (double)trial_vertices * 100.0 <= (double)safe_vertices * trial.required_percent) {
keep_detail_chase = visual_proxy_score(trial.proxy_res) >= trial.proxy_threshold;
}
}
if (keep_detail_chase) {
accepted_detail_chase = true;
break;
}
restore_state(safe_state);
}
if (!accepted_detail_chase) restore_state(safe_state);
}
if (small_visual_restart_profile && elapsed_seconds() < 12.0) {
MeshState current_state = capture_state();
const int current_vertices = count_output_vertices_estimate();
restore_state(small_visual_restart_state);
vector<PassParams> restart_passes;
auto add_restart_pass = [&](double radius, double plane, double angle_deg, bool allow_relocate = false) {
PassParams p;
p.radius_limit = radius;
p.plane_tol = plane;
p.cos_limit = cos(angle_deg * acos(-1.0) / 180.0);
p.exact_plane_tol = exact_plane_tol;
p.exact_cos_limit = exact_cos;
p.min_area2 = min_area2;
p.allow_relocate = allow_relocate;
restart_passes.push_back(p);
};
double restart_radius = 0.033 * d;
double restart_plane = 0.0099 * d;
double restart_angle = 7.5;
double restart_proxy_threshold = 0.960;
double restart_required_percent = 98.0;
if (N < 700) {
restart_radius = 0.033 * d;
restart_plane = 0.0099 * d;
restart_angle = 7.5;
restart_proxy_threshold = 0.965;
restart_required_percent = 98.0;
} else if (N < 1500) {
restart_radius = 0.034 * d;
restart_plane = 0.0102 * d;
restart_angle = 7.8;
restart_proxy_threshold = 0.942;
restart_required_percent = 97.0;
} else {
restart_radius = 0.036 * d;
restart_plane = 0.0108 * d;
restart_angle = 8.2;
restart_proxy_threshold = 0.940;
restart_required_percent = 97.0;
}
add_restart_pass(0.004 * d, 0.0010 * d, 1.0);
add_restart_pass(0.008 * d, 0.0020 * d, 1.8);
add_restart_pass(0.012 * d, 0.0032 * d, 2.8);
add_restart_pass(restart_radius, restart_plane, restart_angle);
add_restart_pass(restart_radius, restart_plane, restart_angle);
bool keep_restart = false;
if (current_vertices > 0
&& run_pass_list(restart_passes, true)
&& elapsed_seconds() < 16.2) {
const int restart_vertices = count_output_vertices_estimate();
if (restart_vertices > 0
&& restart_vertices < current_vertices
&& (double)restart_vertices * 100.0 <= (double)current_vertices * restart_required_percent) {
const double proxy = visual_proxy_score(512);
keep_restart = proxy >= restart_proxy_threshold;
}
}
if (!keep_restart) restore_state(current_state);
}
if (tiny_bumpy_chase_profile && elapsed_seconds() < 12.3) {
MeshState current_state = capture_state();
const int current_vertices = count_output_vertices_estimate();
restore_state(initial_state);
vector<PassParams> tiny_passes;
auto add_tiny_pass = [&](double radius, double plane, double angle_deg) {
PassParams p;
p.radius_limit = radius;
p.plane_tol = plane;
p.cos_limit = cos(angle_deg * acos(-1.0) / 180.0);
p.exact_plane_tol = exact_plane_tol;
p.exact_cos_limit = exact_cos;
p.min_area2 = min_area2;
tiny_passes.push_back(p);
};
add_tiny_pass(0.004 * d, 0.0010 * d, 1.0);
add_tiny_pass(0.008 * d, 0.0020 * d, 1.8);
add_tiny_pass(0.012 * d, 0.0032 * d, 2.8);
add_tiny_pass(0.032 * d, 0.0096 * d, 7.8);
add_tiny_pass(0.032 * d, 0.0096 * d, 7.8);
add_tiny_pass(0.032 * d, 0.0096 * d, 7.8);
bool keep_tiny = false;
if (current_vertices > 0
&& run_pass_list(tiny_passes, true)
&& elapsed_seconds() < 16.2) {
const int tiny_vertices = count_output_vertices_estimate();
if (tiny_vertices > 0
&& tiny_vertices < current_vertices
&& (double)tiny_vertices * 100.0 <= (double)current_vertices * 96.0) {
keep_tiny = visual_proxy_score(768) >= 0.935;
}
}
if (!keep_tiny) restore_state(current_state);
}
if (attempt_aggressive_profile && elapsed_seconds() < 11.5) {
MeshState safe_state = capture_state();
const int safe_vertices = count_output_vertices_estimate();
restore_state(initial_state);
vector<PassParams> aggressive_passes;
auto add_aggressive_pass = [&](double radius, double plane, double angle_deg, bool allow_relocate = false) {
PassParams p;
p.radius_limit = radius;
p.plane_tol = plane;
p.cos_limit = cos(angle_deg * acos(-1.0) / 180.0);
p.exact_plane_tol = exact_plane_tol;
p.exact_cos_limit = exact_cos;
p.min_area2 = min_area2;
p.allow_relocate = allow_relocate;
aggressive_passes.push_back(p);
};
double aggressive_radius = 0.0265 * d;
double aggressive_plane = 0.0079 * d;
double aggressive_angle = 5.9;
if (N >= 200000) {
aggressive_radius = 0.033 * d;
aggressive_plane = 0.0100 * d;
aggressive_angle = 7.4;
} else if (N < 10000) {
aggressive_radius = 0.024 * d;
aggressive_plane = 0.0072 * d;
aggressive_angle = 5.3;
} else if (N < 30000) {
aggressive_radius = 0.0265 * d;
aggressive_plane = 0.0079 * d;
aggressive_angle = 5.9;
}
add_aggressive_pass(0.004 * d, 0.0010 * d, 1.0);
add_aggressive_pass(0.008 * d, 0.0020 * d, 1.8);
add_aggressive_pass(0.012 * d, 0.0032 * d, 2.8);
add_aggressive_pass(aggressive_radius, aggressive_plane, aggressive_angle);
add_aggressive_pass(aggressive_radius, aggressive_plane, aggressive_angle);
bool keep_aggressive = false;
if (run_pass_list(aggressive_passes, true) && elapsed_seconds() < 14.0) {
if (small_detail_aggressive_profile && !smooth_lift_allowed && !attempt_detail_profile) {
const int trial_vertices = count_output_vertices_estimate();
if (trial_vertices > 0
&& safe_vertices > 0
&& (double)trial_vertices * 100.0 <= (double)safe_vertices * 92.0
&& elapsed_seconds() < 16.2) {
const double proxy = visual_proxy_score(512);
keep_aggressive = proxy >= 0.915;
}
} else {
const double proxy = visual_proxy_score(64);
keep_aggressive = proxy >= 0.955;
}
}
if (!keep_aggressive) restore_state(safe_state);
}
if (N >= 200000 && elapsed_seconds() < 12.5) {
MeshState current_state = capture_state();
const int current_vertices = count_output_vertices_estimate();
restore_state(initial_state);
vector<PassParams> ultra_passes;
auto add_ultra_pass = [&](double radius, double plane, double angle_deg) {
PassParams p;
p.radius_limit = radius;
p.plane_tol = plane;
p.cos_limit = cos(angle_deg * acos(-1.0) / 180.0);
p.exact_plane_tol = exact_plane_tol;
p.exact_cos_limit = exact_cos;
p.min_area2 = min_area2;
ultra_passes.push_back(p);
};
add_ultra_pass(0.004 * d, 0.0010 * d, 1.0);
add_ultra_pass(0.008 * d, 0.0020 * d, 1.8);
add_ultra_pass(0.012 * d, 0.0032 * d, 2.8);
add_ultra_pass(0.045 * d, 0.0135 * d, 10.2);
add_ultra_pass(0.045 * d, 0.0135 * d, 10.2);
bool keep_ultra = false;
if (run_pass_list(ultra_passes, true) && elapsed_seconds() < 15.8) {
const int ultra_vertices = count_output_vertices_estimate();
if (ultra_vertices > 0 && ultra_vertices * 100 <= current_vertices * 96) {
const double proxy = visual_proxy_score(64);
keep_ultra = proxy >= 0.955;
}
}
if (!keep_ultra) restore_state(current_state);
}
if (attempt_detail_profile && elapsed_seconds() < 12.5) {
MeshState current_state = capture_state();
const int current_vertices = count_output_vertices_estimate();
restore_state(initial_state);
vector<PassParams> detail_passes;
auto add_detail_pass = [&](double radius, double plane, double angle_deg) {
PassParams p;
p.radius_limit = radius;
p.plane_tol = plane;
p.cos_limit = cos(angle_deg * acos(-1.0) / 180.0);
p.exact_plane_tol = exact_plane_tol;
p.exact_cos_limit = exact_cos;
p.min_area2 = min_area2;
detail_passes.push_back(p);
};
double detail_radius = 0.0310 * d;
double detail_plane = 0.00930 * d;
double detail_angle = 6.95;
double detail_threshold = 0.955;
int detail_proxy_res = 128;
int final_repeats = 3;
if (large_detail_profile) {
if (N >= 200000) {
detail_radius = 0.050 * d;
detail_plane = 0.0150 * d;
detail_angle = 11.0;
} else {
detail_radius = 0.047 * d;
detail_plane = 0.0141 * d;
detail_angle = 10.5;
}
detail_threshold = 0.940;
final_repeats = 4;
}
add_detail_pass(0.004 * d, 0.0010 * d, 1.0);
add_detail_pass(0.008 * d, 0.0020 * d, 1.8);
add_detail_pass(0.012 * d, 0.0032 * d, 2.8);
for (int i = 0; i < final_repeats; ++i) {
add_detail_pass(detail_radius, detail_plane, detail_angle);
}
bool keep_detail = false;
if (run_pass_list(detail_passes, true) && elapsed_seconds() < 16.2) {
const int detail_vertices = count_output_vertices_estimate();
if (detail_vertices > 0 && current_vertices > 0
&& detail_vertices * 100 <= current_vertices * 96) {
const double proxy = visual_proxy_score(detail_proxy_res);
keep_detail = proxy >= detail_threshold;
}
}
if (!keep_detail) restore_state(current_state);
}
if (large_detail_profile && N >= 30000 && N < 100000 && elapsed_seconds() < 12.8) {
MeshState current_state = capture_state();
const int current_vertices = count_output_vertices_estimate();
restore_state(initial_state);
vector<PassParams> boost_passes;
auto add_boost_pass = [&](double radius, double plane, double angle_deg) {
PassParams p;
p.radius_limit = radius;
p.plane_tol = plane;
p.cos_limit = cos(angle_deg * acos(-1.0) / 180.0);
p.exact_plane_tol = exact_plane_tol;
p.exact_cos_limit = exact_cos;
p.min_area2 = min_area2;
boost_passes.push_back(p);
};
add_boost_pass(0.004 * d, 0.0010 * d, 1.0);
add_boost_pass(0.008 * d, 0.0020 * d, 1.8);
add_boost_pass(0.012 * d, 0.0032 * d, 2.8);
for (int i = 0; i < 4; ++i) {
add_boost_pass(0.055 * d, 0.0165 * d, 12.2);
}
bool keep_boost = false;
if (run_pass_list(boost_passes, true) && elapsed_seconds() < 16.3) {
const int boost_vertices = count_output_vertices_estimate();
if (boost_vertices > 0 && current_vertices > 0
&& boost_vertices * 100 <= current_vertices * 95) {
const double proxy = visual_proxy_score(128);
keep_boost = proxy >= 0.960;
}
}
if (!keep_boost) restore_state(current_state);
}
if (certified_boost_profile && elapsed_seconds() < 13.2) {
MeshState boost_safe_state = capture_state();
const int boost_start_vertices = count_output_vertices_estimate();
vector<PassParams> certified_passes;
auto add_certified_pass = [&](double radius, double plane, double angle_deg, bool allow_relocate = false) {
PassParams p;
p.radius_limit = radius;
p.plane_tol = plane;
p.cos_limit = cos(angle_deg * acos(-1.0) / 180.0);
p.exact_plane_tol = exact_plane_tol;
p.exact_cos_limit = exact_cos;
p.min_area2 = min_area2;
p.allow_relocate = allow_relocate;
certified_passes.push_back(p);
};
double boost_radius = 0.036 * d;
double boost_plane = 0.0108 * d;
double boost_angle = 8.0;
double keep_ratio_percent = 88.0;
double proxy_threshold = 0.955;
int proxy_res = 384;
bool relocate_on_first = false;
if (N < 8000) {
boost_radius = 0.048 * d;
boost_plane = 0.0144 * d;
boost_angle = 10.5;
keep_ratio_percent = 96.0;
proxy_threshold = 0.935;
proxy_res = 512;
relocate_on_first = smooth10_ratio_cache >= 0.900 && sharp22_ratio_cache <= 0.035;
} else if (N < 30000) {
boost_radius = 0.052 * d;
boost_plane = 0.0156 * d;
boost_angle = 11.5;
keep_ratio_percent = 95.0;
proxy_threshold = 0.930;
proxy_res = 384;
relocate_on_first = smooth10_ratio_cache >= 0.930 && sharp22_ratio_cache <= 0.040;
} else if (N < 80000) {
boost_radius = 0.050 * d;
boost_plane = 0.0150 * d;
boost_angle = 11.1;
keep_ratio_percent = 95.0;
proxy_threshold = 0.910;
proxy_res = 128;
} else if (N < 140000) {
boost_radius = 0.050 * d;
boost_plane = 0.0150 * d;
boost_angle = 11.1;
keep_ratio_percent = 95.0;
proxy_threshold = 0.910;
proxy_res = 128;
} else {
boost_radius = 0.056 * d;
boost_plane = 0.0168 * d;
boost_angle = 12.4;
keep_ratio_percent = 95.0;
proxy_threshold = 0.920;
proxy_res = 96;
}
add_certified_pass(boost_radius, boost_plane, boost_angle, relocate_on_first);
add_certified_pass(boost_radius, boost_plane, boost_angle, false);
bool keep_certified = false;
if (boost_start_vertices > 0
&& run_pass_list(certified_passes, true)
&& elapsed_seconds() < 16.6) {
const int boosted_vertices = count_output_vertices_estimate();
if (boosted_vertices > 0
&& (double)boosted_vertices * 100.0 <= (double)boost_start_vertices * keep_ratio_percent) {
keep_certified = visual_proxy_score(proxy_res) >= proxy_threshold;
}
}
if (!keep_certified) restore_state(boost_safe_state);
}
if (attempt_relocation_profile && elapsed_seconds() < 12.5) {
MeshState endpoint_state = capture_state();
const int endpoint_vertices = count_output_vertices_estimate();
restore_state(initial_state);
vector<PassParams> relocate_passes;
auto add_relocate_pass = [&](double radius, double plane, double angle_deg, bool allow_relocate = false) {
PassParams p;
p.radius_limit = radius;
p.plane_tol = plane;
p.cos_limit = cos(angle_deg * acos(-1.0) / 180.0);
p.exact_plane_tol = exact_plane_tol;
p.exact_cos_limit = exact_cos;
p.min_area2 = min_area2;
p.allow_relocate = allow_relocate;
relocate_passes.push_back(p);
};
double relocate_radius = 0.024 * d;
double relocate_plane = 0.0072 * d;
double relocate_angle = 5.3;
if (N >= 10000) {
relocate_radius = 0.0265 * d;
relocate_plane = 0.0079 * d;
relocate_angle = 5.9;
}
add_relocate_pass(0.004 * d, 0.0010 * d, 1.0);
add_relocate_pass(0.008 * d, 0.0020 * d, 1.8);
add_relocate_pass(0.012 * d, 0.0032 * d, 2.8);
add_relocate_pass(relocate_radius, relocate_plane, relocate_angle, true);
add_relocate_pass(relocate_radius, relocate_plane, relocate_angle, true);
bool keep_relocated = false;
if (run_pass_list(relocate_passes, true) && elapsed_seconds() < 15.5) {
const int relocated_vertices = count_output_vertices_estimate();
if (relocated_vertices > 0 && relocated_vertices * 100 <= endpoint_vertices * 94) {
const double proxy = visual_proxy_score(64);
keep_relocated = proxy >= 0.970;
}
}
if (!keep_relocated) restore_state(endpoint_state);
}
if (elapsed_seconds() < 13.0) {
(void)try_strict_sphere_remesh();
}
if (elapsed_seconds() < 13.0) {
(void)try_strict_ellipsoid_remesh();
}
if (elapsed_seconds() < 13.8) {
(void)try_radial_sphere_remesh();
}
if (elapsed_seconds() < 13.8) {
(void)try_axis_torus_remesh();
}
if (elapsed_seconds() < 13.8) {
(void)try_axis_revolved_primitive_remesh();
}
if (elapsed_seconds() < 13.8) {
(void)try_axis_capsule_remesh();
}
if (smooth_stats_valid
&& N >= 1000 && N < 60000
&& bad_ratio_cache <= 0.010
&& smooth30_ratio_cache >= 0.760
&& sharp45_ratio_cache <= 0.090
&& elapsed_seconds() < 12.6) {
MeshState generic_safe_state = capture_state();
const int generic_start_vertices = count_output_vertices_estimate();
vector<PassParams> generic_passes;
auto add_generic_pass = [&](double radius, double plane, double angle_deg, bool allow_relocate = false) {
PassParams p;
p.radius_limit = radius;
p.plane_tol = plane;
p.cos_limit = cos(angle_deg * acos(-1.0) / 180.0);
p.exact_plane_tol = exact_plane_tol;
p.exact_cos_limit = exact_cos;
p.min_area2 = min_area2;
p.allow_relocate = allow_relocate;
generic_passes.push_back(p);
};
double generic_radius = 0.026 * d;
double generic_plane = 0.0078 * d;
double generic_angle = 6.0;
const double required_percent = 96.0;
double proxy_threshold = 0.925;
int proxy_res = 1024;
if (N >= 6000) {
generic_radius = 0.034 * d;
generic_plane = 0.0102 * d;
generic_angle = 7.8;
proxy_threshold = 0.925;
proxy_res = 1024;
}
if (N >= 20000) {
generic_radius = 0.040 * d;
generic_plane = 0.0120 * d;
generic_angle = 9.0;
proxy_threshold = 0.930;
proxy_res = 512;
}
if (sharp22_ratio_cache > 0.120 || sharp45_ratio_cache > 0.050) {
generic_radius *= 0.82;
generic_plane *= 0.82;
generic_angle *= 0.86;
proxy_threshold += 0.020;
}
add_generic_pass(0.004 * d, 0.0010 * d, 1.0, false);
add_generic_pass(0.008 * d, 0.0020 * d, 1.8, false);
add_generic_pass(0.012 * d, 0.0032 * d, 2.8, false);
add_generic_pass(generic_radius, generic_plane, generic_angle, true);
add_generic_pass(generic_radius, generic_plane, generic_angle, true);
bool keep_generic = false;
if (generic_start_vertices > 0
&& run_pass_list(generic_passes, true)
&& elapsed_seconds() < 16.4) {
const int generic_vertices = count_output_vertices_estimate();
if (generic_vertices > 0
&& generic_vertices < generic_start_vertices
&& (double)generic_vertices * 100.0 <= (double)generic_start_vertices * required_percent) {
keep_generic = visual_proxy_score(proxy_res) >= proxy_threshold;
}
}
if (!keep_generic) restore_state(generic_safe_state);
}
if (N >= 1000 && N <= 70000 && elapsed_seconds() < 13.2) {
MeshState heap_safe_state = capture_state();
const int heap_safe_vertices = count_output_vertices_estimate();
vector<Vec3> heap_vertices;
vector<Face> heap_faces;
const double heap_budget = (N < 30000 ? 4.2 : 5.2);
bool keep_heap_candidate = false;
if (heap_safe_vertices > 0 && N >= 30000) {
if (AggressiveHeapCandidate::build_from_parent(originalP, originalFaces, mesh_diag,
heap_vertices, heap_faces, heap_budget)
&& (int)heap_vertices.size() < heap_safe_vertices
&& elapsed_seconds() < 17.2) {
if (install_replacement_mesh(heap_vertices, heap_faces) && elapsed_seconds() < 17.4) {
const double proxy = visual_proxy_score(1024);
keep_heap_candidate = proxy >= 0.902;
}
}
} else if (heap_safe_vertices > 0) {
const double ratio_trials[] = {0.090, -1.0, 0.115, 0.140, 0.180, 0.250, 0.350, 0.500};
for (double ratio : ratio_trials) {
if (elapsed_seconds() >= 17.0) break;
heap_vertices.clear();
heap_faces.clear();
restore_state(heap_safe_state);
if (!AggressiveHeapCandidate::build_from_parent(originalP, originalFaces, mesh_diag,
heap_vertices, heap_faces, heap_budget,
ratio)) {
continue;
}
if ((int)heap_vertices.size() >= heap_safe_vertices || elapsed_seconds() >= 17.2) continue;
if (!install_replacement_mesh(heap_vertices, heap_faces) || elapsed_seconds() >= 17.4) continue;
const double proxy = visual_proxy_score(1024);
const double proxy_threshold = (ratio < 0.0 ? 0.904 : 0.925);
if (proxy >= proxy_threshold) {
keep_heap_candidate = true;
break;
}
}
}
if (!keep_heap_candidate) restore_state(heap_safe_state);
}
if (N >= 300 && N <= 70000 && elapsed_seconds() < 12.9) {
MeshState certified_base_state = capture_state();
const int certified_base_vertices = count_output_vertices_estimate();
auto run_sorted_edge_pass = [&](const PassParams& params, int repeats) -> bool {
for (int rep = 0; rep < repeats; ++rep) {
if (elapsed_seconds() > 16.2) return false;
vector<unsigned long long> keys;
keys.reserve((size_t)active_faces * 3);
for (int fid = 0; fid < (int)faces.size(); ++fid) {
if ((fid & 8191) == 0 && elapsed_seconds() > 16.2) return false;
if (!alive_f[fid]) continue;
const Face& f = faces[fid];
if (!alive_v[f.v[0]] || !alive_v[f.v[1]] || !alive_v[f.v[2]]) continue;
keys.push_back(edge_key_out(f.v[0], f.v[1]));
keys.push_back(edge_key_out(f.v[1], f.v[2]));
keys.push_back(edge_key_out(f.v[2], f.v[0]));
}
sort(keys.begin(), keys.end());
keys.erase(unique(keys.begin(), keys.end()), keys.end());
struct SortedEdge {
double len2;
unsigned long long key;
};
vector<SortedEdge> edges;
edges.reserve(keys.size());
for (unsigned long long key : keys) {
int a = (int)(key >> 32);
int b = (int)(key & 0xffffffffu);
if (a >= 0 && a < N && b >= 0 && b < N && alive_v[a] && alive_v[b]) {
edges.push_back({norm2(P[a] - P[b]), key});
}
}
sort(edges.begin(), edges.end(), [](const SortedEdge& lhs, const SortedEdge& rhs) {
return lhs.len2 < rhs.len2;
});
int reduced = 0;
for (const SortedEdge& e : edges) {
if ((reduced & 511) == 0 && elapsed_seconds() > 16.2) return false;
int a = (int)(e.key >> 32);
int b = (int)(e.key & 0xffffffffu);
if (a >= 0 && a < N && b >= 0 && b < N && alive_v[a] && alive_v[b]) {
if (try_collapse_edge(a, b, params)) ++reduced;
}
}
if (reduced == 0) break;
}
return true;
};
if (certified_base_vertices > 0) {
struct CertifiedTrial {
double radius;
double plane;
double angle_deg;
int repeats;
int proxy_res;
double proxy_threshold;
double required_percent;
bool relocate;
bool sorted;
};
vector<CertifiedTrial> trials;
const bool smoothish = smooth_stats_valid
&& bad_ratio_cache <= 0.012
&& smooth30_ratio_cache >= 0.760
&& sharp45_ratio_cache <= 0.110;
if (N < 1000) {
trials.push_back({0.036 * d, 0.0140 * d, 10.0, 2, 768, 0.930, 99.0, true, false});
trials.push_back({0.048 * d, 0.0270 * d, 20.0, 4, 768, smoothish ? 0.906 : 0.925, 98.0, true, true});
} else if (N < 3000) {
trials.push_back({0.040 * d, 0.0180 * d, 13.5, 3, 1024, 0.918, 99.0, true, false});
trials.push_back({0.049 * d, 0.0340 * d, 28.0, 6, 1024, smoothish ? 0.903 : 0.922, 98.0, true, true});
} else if (N < 8000) {
trials.push_back({0.044 * d, 0.0220 * d, 17.0, 4, 1024, 0.915, 99.0, true, false});
trials.push_back({0.049 * d, 0.0400 * d, 34.0, 7, 1024, smoothish ? 0.904 : 0.924, 98.0, true, true});
} else if (N < 25000) {
trials.push_back({0.047 * d, 0.0250 * d, 19.0, 4, 1024, 0.916, 99.0, true, false});
trials.push_back({0.049 * d, 0.0410 * d, 34.0, 6, 1024, smoothish ? 0.902 : 0.930, 98.0, true, true});
} else if (N < 60000) {
trials.push_back({0.047 * d, 0.0230 * d, 17.5, 3, 512, 0.932, 99.0, true, false});
trials.push_back({0.049 * d, 0.0350 * d, 28.0, 4, 512, smoothish ? 0.914 : 0.940, 98.0, true, true});
} else if (N < 200000) {
trials.push_back({0.047 * d, 0.0220 * d, 16.0, 3, 256, 0.940, 99.0, true, false});
trials.push_back({0.049 * d, 0.0320 * d, 24.0, 4, 256, smoothish ? 0.930 : 0.948, 98.0, true, N < 120000});
} else {
trials.push_back({0.047 * d, 0.0220 * d, 16.0, 3, 128, 0.950, 99.0, true, false});
trials.push_back({0.049 * d, 0.0300 * d, 22.0, 4, 128, smoothish ? 0.940 : 0.955, 98.0, true, false});
}
MeshState best_state = certified_base_state;
int best_vertices = certified_base_vertices;
for (const CertifiedTrial& trial : trials) {
if (elapsed_seconds() > 14.5) break;
restore_state(certified_base_state);
PassParams p;
p.radius_limit = trial.radius;
p.plane_tol = trial.plane;
p.cos_limit = cos(trial.angle_deg * acos(-1.0) / 180.0);
p.exact_plane_tol = exact_plane_tol;
p.exact_cos_limit = exact_cos;
p.min_area2 = min_area2;
p.allow_relocate = trial.relocate;
bool ran = false;
if (trial.sorted) {
vector<PassParams> warmup;
PassParams w = p;
w.radius_limit = min(p.radius_limit, 0.020 * d);
w.plane_tol = min(p.plane_tol, 0.0060 * d);
w.cos_limit = cos(min(trial.angle_deg, 4.8) * acos(-1.0) / 180.0);
w.allow_relocate = false;
warmup.push_back(w);
ran = run_pass_list(warmup, true) && run_sorted_edge_pass(p, trial.repeats);
} else {
vector<PassParams> pass_list;
for (int i = 0; i < trial.repeats; ++i) pass_list.push_back(p);
ran = run_pass_list(pass_list, true);
}
bool keep_trial = false;
if (ran && elapsed_seconds() < 16.5) {
const int trial_vertices = count_output_vertices_estimate();
if (trial_vertices > 0
&& trial_vertices < best_vertices
&& (double)trial_vertices * 100.0 <=
(double)certified_base_vertices * trial.required_percent) {
const double proxy = visual_proxy_score(trial.proxy_res);
keep_trial = proxy >= trial.proxy_threshold;
}
if (keep_trial) {
best_state = capture_state();
best_vertices = trial_vertices;
}
}
}
restore_state(best_state);
}
}
if (elapsed_seconds() < 16.6) {
(void)try_oriented_frustum_remesh();
}
}
static void append_output_line(string& out, const char* line, int len) {
if (out.size() + (size_t)len > (1 << 20)) {
fwrite(out.data(), 1, out.size(), stdout);
out.clear();
}
out.append(line, line + len);
}
static void save_original_mesh() {
string out;
out.reserve(1 << 20);
char line[128];
int len = snprintf(line, sizeof(line), "%d %d\n", N, M);
append_output_line(out, line, len);
for (int i = 0; i < N; ++i) {
len = snprintf(line, sizeof(line), "v %.10g %.10g %.10g\n",
originalP[i].x, originalP[i].y, originalP[i].z);
append_output_line(out, line, len);
}
for (int fid = 0; fid < M; ++fid) {
const Face& f = originalFaces[fid];
len = snprintf(line, sizeof(line), "f %d %d %d\n",
f.v[0] + 1, f.v[1] + 1, f.v[2] + 1);
append_output_line(out, line, len);
}
if (!out.empty()) fwrite(out.data(), 1, out.size(), stdout);
}
static bool build_compacted_mesh(vector<int>& used_old, vector<Face>& out_faces) {
vector<int> id(N, -1);
used_old.clear();
out_faces.clear();
used_old.reserve(N);
out_faces.reserve(active_faces);
const double area_eps = max(1e-30, 1e-24 * mesh_diag * mesh_diag);
for (int fid = 0; fid < (int)faces.size(); ++fid) {
if (!alive_f[fid]) continue;
const Face& f = faces[fid];
for (int k = 0; k < 3; ++k) {
if (f.v[k] < 0 || f.v[k] >= N || !alive_v[f.v[k]]) return false;
}
if (f.v[0] == f.v[1] || f.v[0] == f.v[2] || f.v[1] == f.v[2]) return false;
Vec3 cr = face_cross(fid);
if (!(norm2(cr) > area_eps)) return false;
Face compact;
for (int k = 0; k < 3; ++k) {
int old = f.v[k];
if (id[old] < 0) {
id[old] = (int)used_old.size();
used_old.push_back(old);
}
compact.v[k] = id[old];
}
out_faces.push_back(compact);
}
if (used_old.empty() || out_faces.empty() || (int)used_old.size() > N) return false;
vector<unsigned long long> edges;
edges.reserve(out_faces.size() * 3);
for (const Face& f : out_faces) {
edges.push_back(edge_key_out(f.v[0], f.v[1]));
edges.push_back(edge_key_out(f.v[1], f.v[2]));
edges.push_back(edge_key_out(f.v[2], f.v[0]));
}
sort(edges.begin(), edges.end());
for (size_t i = 0; i < edges.size();) {
size_t j = i + 1;
while (j < edges.size() && edges[j] == edges[i]) ++j;
if (j - i != 2) return false;
i = j;
}
return true;
}
static void save_mesh() {
vector<int> used_old;
vector<Face> out_faces;
if (!build_compacted_mesh(used_old, out_faces)) {
save_original_mesh();
return;
}
string out;
out.reserve(1 << 20);
char line[128];
int len = snprintf(line, sizeof(line), "%d %d\n",
(int)used_old.size(), (int)out_faces.size());
append_output_line(out, line, len);
const bool high_precision_vertices = (int)used_old.size() * 2 <= N;
for (int old : used_old) {
if (high_precision_vertices) {
len = snprintf(line, sizeof(line), "v %.15g %.15g %.15g\n",
P[old].x, P[old].y, P[old].z);
} else {
len = snprintf(line, sizeof(line), "v %.10g %.10g %.10g\n",
P[old].x, P[old].y, P[old].z);
}
append_output_line(out, line, len);
}
for (const Face& f : out_faces) {
len = snprintf(line, sizeof(line), "f %d %d %d\n",
f.v[0] + 1, f.v[1] + 1, f.v[2] + 1);
append_output_line(out, line, len);
}
if (!out.empty()) fwrite(out.data(), 1, out.size(), stdout);
}
int main() {
load_mesh();
simplify_mesh();
save_mesh();
return 0;
}
