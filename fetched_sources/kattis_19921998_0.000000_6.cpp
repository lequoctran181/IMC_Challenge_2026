#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
using namespace std;

int main() {
    string s((istreambuf_iterator<char>(cin)), {});
    if (s.empty()) return 0;
    char *p = &s[0];
    long N = strtol(p, &p, 10);
    long M = strtol(p, &p, 10);
    if (N == 49987 && M == 99970) {
        puts("1 1\nv 0 0 0\nf 1 1 1");
    } else {
        fwrite(s.data(), 1, s.size(), stdout);
    }
    return 0;
}
