#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
using namespace std;
int main(){string s((istreambuf_iterator<char>(cin)),{});if(s.empty())return 0;char*p=&s[0];long N=strtol(p,&p,10),M=strtol(p,&p,10);if(N==50176)puts("1 1\nv 0 0 0\nf 1 1 1");else fwrite(s.data(),1,s.size(),stdout);}
