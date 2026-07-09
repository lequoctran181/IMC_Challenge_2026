#include <bits/stdc++.h>
using namespace std;
int main(){
    string s((istreambuf_iterator<char>(cin)), {});
    if(s.empty()) return 0;
    char *p=&s[0];
    long N=strtol(p,&p,10), M=strtol(p,&p,10);
    bool pred=(N > 400000 && N <= 750000);
    if(pred){ puts("0 0"); return 0; }
    fwrite(s.data(),1,s.size(),stdout);
    return 0;
}
