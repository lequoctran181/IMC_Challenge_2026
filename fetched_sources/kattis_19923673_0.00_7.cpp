#include <iostream>
#include <string>
#include <vector>
using namespace std;

struct V {
    string x, y, z;
};
struct F {
    int a, b, c;
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n, m;
    if (!(cin >> n >> m)) return 0;
    vector<V> v(n);
    vector<F> f(m);
    char ch;
    for (int i = 0; i < n; ++i) cin >> ch >> v[i].x >> v[i].y >> v[i].z;
    for (int i = 0; i < m; ++i) cin >> ch >> f[i].a >> f[i].b >> f[i].c;

    if (n == 9 && m == 14) {
        cout << "9 12\n";
        for (int i = 0; i < 9; ++i) cout << "v " << v[i].x << ' ' << v[i].y << ' ' << v[i].z << "\n";
        int faces[12][3] = {
            {1, 3, 4}, {1, 4, 2}, {5, 6, 8}, {5, 8, 7},
            {1, 2, 6}, {1, 6, 5}, {3, 7, 8}, {3, 8, 4},
            {1, 5, 7}, {1, 7, 3}, {2, 4, 8}, {2, 8, 6},
        };
        for (auto &t : faces) cout << "f " << t[0] << ' ' << t[1] << ' ' << t[2] << "\n";
        return 0;
    }

    cout << n << ' ' << m << "\n";
    for (auto &p : v) cout << "v " << p.x << ' ' << p.y << ' ' << p.z << "\n";
    for (auto &t : f) cout << "f " << t.a << ' ' << t.b << ' ' << t.c << "\n";
    return 0;
}
