---
layout:     post   				    # 使用的布局（不需要改）
title:      多项式系列算法代码集				# 标题 
subtitle:   最近背背板子吧。。 #副标题
date:       2019-1-30			# 时间
author:     signron					# 作者
header-img: img/bg.jpg 	#这篇文章标题背景图片
catalog: true 						# 是否归档
tags:								#标签
    - FFT
    - NTT
    - MTT
    - Codes
    - Templates
---

## 多项式乘法-FFT

没有模数，~~于是拿NTT水过去了~~

```cpp
#include<cstdio>
#include<iostream>
using namespace std;

const int MAXN=5e6+5,mo=998244353,g=3;

int n,m,limit,len,invn,a[MAXN],b[MAXN],p[MAXN];

inline int read() {
	int x=0,f=1;
	char ch=getchar();
	while(ch<'0'||ch>'9') {
		if(ch=='-')f=-1;
		ch=getchar();
	}
	while(ch>='0'&&ch<='9') {
		x=x*10+ch-'0';
		ch=getchar();
	}
	return x*f;
}

int power(int a,int b) {
	int res=1;
	while (b) {
		if (b&1) res=1ll*res*a%mo;
		a=1ll*a*a%mo;
		b>>=1;
	}
	return res;
}

void NTT(int *A,int type) {
	for (int i=0; i<limit; i++)
		if (i<p[i]) swap(A[i],A[p[i]]);
	for (int l=1; l<limit; l<<=1) {
		int wn=power(g,(mo-1)/(l<<1));
		if (type==-1) wn=power(wn,mo-2);
		for (int i=0; i<limit; i+=(l<<1)) {
			int w=1;
			for (int j=0; j<l; j++,w=1ll*w*wn%mo) {
				int t=1ll*w*A[i+j+l]%mo;
				A[i+j+l]=(A[i+j]-t+mo)%mo;
				A[i+j]=(A[i+j]+t)%mo;
			}
		}
	}
}

int main() {
	n=read();
	m=read();
	for (int i=0; i<=n; i++) a[i]=read();
	for (int i=0; i<=m; i++) b[i]=read();
	limit=1;
	while (limit<n+m+1) limit<<=1,len++;
	for(int i=0; i<limit; i++) p[i]=(p[i>>1]>>1)|((i&1)<<(len-1));
	NTT(a,1);
	NTT(b,1);
	for (int i=0; i<limit; i++) a[i]=1ll*a[i]*b[i]%mo;
	NTT(a,-1);
	invn=power(limit,mo-2);
	for (int i=0; i<=n+m; i++) printf("%lld ",1ll*a[i]*invn%mo);
	printf("\n");
	return 0;
}
```

## 任意模数NTT

模数:1e9+7

```cpp
#include <cstdio>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <cmath>
using namespace std;

typedef long long LL;
const int MAXN = 4e5 + 5;
const LL p[] = {0, 469762049ll, 998244353ll, 1004535809ll}, g = 3, MX = 469762049ll * 998244353ll;

LL N, M, MOD,A[MAXN], B[MAXN], C[4][MAXN], F[MAXN], G[MAXN], ans[MAXN];

LL qpow(LL x, int n, LL mod) {
	LL res = 1;
	while(n) {
		if(n & 1) res = res * x % mod;
		x = x * x % mod;
		n >>= 1;
	}
	return res;
}

void bit_reverse(int len, LL *a) {
	int i, j, k;
	for(i = 0, j = 0; i < len; i++) {
		if(i > j) swap(a[i], a[j]);
		for(k = len >> 1; (j & k); j ^= k, k >>= 1);
		j ^= k;
	}
}

void DNT(int len, LL *a, int type, LL mod) {
	bit_reverse(len, a);
	int i, j, l;
	for(l = 2; l <= len; l <<= 1) {
		int mid = l >> 1;
		LL wn = qpow(g, (mod - 1) / l, mod);
		if(type == -1) wn = qpow(wn, mod - 2, mod);
		for(i = 0; i < len; i += l) {
			LL w = 1;
			for(j = 0; j < mid; j++, w = w * wn % mod) {
				LL x = a[i + j], y = w * a[i + j + mid] % mod;
				a[i + j] = (x + y) % mod;
				a[i + j + mid] = (x - y + mod) % mod;
			}
		}
	}
	if(type == -1) {
		int inv = qpow(len, mod - 2, mod);
		for(i = 0; i <= len; i++) a[i] = a[i] * inv % mod;
	}
}

LL multi(LL a, LL b, LL mod) {
	a %= mod, b %= mod;
	return ((a * b - (LL)((LL)((long double)a / mod * b + 1e-3) * mod)) % mod + mod) % mod;
}

void getC(int I, int len) {
	int i;
	memset(F, 0, sizeof(F)), memset(G, 0, sizeof(G));
	for(i = 0; i <= N; i++) F[i] = A[i];
	for(i = 0; i <= M; i++) G[i] = B[i];
	DNT(len, F, 1, p[I]), DNT(len, G, 1, p[I]);
	for(i = 0; i <= len; i++)
		C[I][i] = F[i] * G[i] % p[I];
	DNT(len, C[I], -1, p[I]);
}

int main() {
	ios::sync_with_stdio(false);
	cin >> N >> M >> MOD;
	int i;
	for(i = 0; i <= N; i++) cin >> A[i];
	for(i = 0; i <= M; i++) cin >> B[i];
	int len = 1;
	while(len <= M + N) len <<= 1;
	getC(1, len), getC(2, len), getC(3, len);
	for(i = 0; i <= N + M; i++) {
		LL x, k1;
		x = (multi(C[1][i] * p[2] % MX, qpow(p[2] % p[1], p[1] - 2, p[1]), MX) +
		     multi(C[2][i] * p[1] % MX, qpow(p[1] % p[2], p[2] - 2, p[2]), MX)) % MX;
		k1 = (multi((C[3][i] % p[3] - x % p[3] + p[3]) % p[3], qpow(MX % p[3], p[3] - 2, p[3]), p[3]));
		ans[i] = ((k1 % MOD) * (MX % MOD) + x % MOD) % MOD;
		cout << ans[i] << " ";
	}
	return 0;
}
```

##多项式求逆-NTT

模数:998244353

```cpp
#include <cstdio>
#include <algorithm>
#include <iostream>
using namespace std;

const int Q=500010,MOD=998244353;

inline int add(int a,int b){a+=b;return a<MOD?a:a-MOD;}
inline int sub(int a,int b){a-=b;return a<0?a+MOD:a;}
inline int mul(int a,int b){return 1LL*a*b%MOD;}

inline int ksm(int a,int b){
    int ans=1;
    while(b){
        if(b&1)ans=mul(ans,a);
        b>>=1,a=mul(a,a);
    }
    return ans;
}


int w[Q];

void NTT(int a[],int n,int f){
    register int i,j,k,m,cnt=1;
    for(i=j=0;i<n;i++){
        if(i<j)swap(a[i],a[j]);
        for(k=n>>1;(j^=k)<k;k>>=1);
    }
    w[0]=1;
    for(m=1;m<n;m<<=1,++cnt){
        int _now=ksm(3,(f*((MOD-1)>>cnt)+MOD-1)%(MOD-1));
        for(i=1;i<m;i++)w[i]=mul(w[i-1],_now);
        for(i=0;i<n;i+=(m<<1))
            for(j=0;j<m;j++){
                int p=a[i+j],q=mul(a[i+j+m],w[j]);
                a[i+j]=add(p,q),a[i+j+m]=sub(p,q);
            }
    }
    if(f<0){
        int n_ni=ksm(n,MOD-2);
        for(i=0;i<n;i++)a[i]=mul(a[i],n_ni);
    }
}

int it[Q];

void GetInv(int a[],int ni[],int toap) { 
    int i,now; 
    ni[0]=ksm(a[0],MOD-2); 
    for(now=2;now<=toap;now<<=1){ 
        fill(ni+(now>>1),ni+(now<<1),0); 
        copy(a,a+now,it); 
        fill(it+now,it+(now<<1),0); 
        NTT(ni,now<<1,1); 
        NTT(it,now<<1,1); 
        for(i=0;i<(now<<1);i++) 
            ni[i]=mul(ni[i],sub(2,mul(it[i],ni[i]))); 
        NTT(ni,now<<1,-1); 
    } 
} 

int b[Q],b_ni[Q];

int main(){
	int p,now;
	scanf("%d",&p);--p;
	for(int i=0;i<=p;++i) scanf("%d",b+i);
	for(now=1;now<=p;now<<=1);
	GetInv(b,b_ni,now);
	for(int i=0;i<=p;++i) printf("%d ",b_ni[i]%MOD);
}
```

## 多项式求逆-MTT

模数:1e9+7，~~居然NTT+卡常能过~~

```cpp
#include<bits/stdc++.h>
#define rep(i, j) for ( int i = 0, i##_end_ = j; i < i##_end_; ++ i)
#define For(i, j ,k) for ( int i = (j), i##_end_ = (k); i <= i##_end_; ++ i)
using namespace std;

typedef long long LL;
const int maxn = 1e5 + 7, MOD = 1e9 + 7;
const long double Pi = acos(-1.);

inline int read() {
    char ch;
    int x=0, f=1;
    for (ch = getchar(); !isdigit(ch); ch = getchar()) if (ch == '-')  f = -1;
    for ( ; isdigit(ch); ch = getchar()) x = (x << 1) + (x << 3) + (ch ^ 48);
    return x * f;
}

int rev[maxn << 2];

inline int ad(int x, int y) {
    if ((x += y) >= MOD) return x - MOD;
    return x;
}

inline LL Pow(LL a, LL b) {
    static LL Ans;
    for (Ans = 1, a %= MOD; b; b >>= 1, (a *= a) %= MOD) if (b & 1) (Ans *= a) %= MOD;
    return Ans;
}

struct Complex {
    long double x, y;
};

inline Complex operator + (const Complex &aa, const Complex &bb) {
    return (Complex) {
        aa.x + bb.x, aa.y + bb.y
    };
}

inline Complex operator - (const Complex &aa, const Complex &bb) {
    return (Complex) {
        aa.x - bb.x, aa.y - bb.y
    };
}

inline Complex operator * (const Complex &aa, const Complex &bb) {
    return (Complex) {
        aa.x * bb.x - aa.y * bb.y, aa.x * bb.y + aa.y * bb.x
    };
}

inline void getr(int n) {
    static int cnt, nn;
    for (cnt = 0, nn = 1; nn < n; nn <<= 1, ++ cnt);
    rep(i, n) rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (cnt - 1));
}

inline void FFT(Complex* a, int n, int ty) {
    rep(i, n) if (i < rev[i]) swap(a[i], a[rev[i]]);
    for ( int i = 2, p = 1; i <= n; p = i, i <<= 1) {
        Complex w0 = Complex {cos(2 * Pi / i), sin(2 * Pi / i) * ty};
        for ( int j = 0; j < n; j += i) {
            Complex w = Complex {1, 0};
            rep(k, p) {
                Complex x = a[j + k], y = a[j + k + p] * w;
                a[j + k] = x + y, a[j + k + p] = x - y;
                w = w * w0;
            }
        }
    }
    if (ty < 0) rep(i, n) a[i].x /= n;
}

inline void multi(int* aa, int* bb, int* cc, int n) {
    static Complex a[maxn << 2], b[maxn << 2], c[maxn << 2], d[maxn << 2];
    rep(i, n) a[i].x = aa[i] >> 15, b[i].x = aa[i] & 32767,c[i].x = bb[i] >> 15, d[i].x = bb[i] & 32767, a[i].y = 0, b[i].y = 0, c[i].y = 0, d[i].y = 0;
    FFT(a, n, 1), FFT(b, n, 1), FFT(c, n, 1), FFT(d, n, 1);
    rep(i, n) {
        Complex ta = a[i], tb = b[i], tc = c[i], td = d[i];
        a[i] = ta * tc, b[i] = tb * tc + ta * td, c[i] = tb * td;
    }
    FFT(a, n, -1), FFT(b, n, -1), FFT(c, n, -1);
    rep(i, n) cc[i] = ad(ad((LL)round(a[i].x) % MOD * ((1ll << 30) % MOD) % MOD, (LL)round(b[i].x) % MOD * (1ll << 15) % MOD), (LL)round(c[i].x) % MOD);
}

int aa[maxn << 2], bb[maxn << 2], cc[maxn << 2];
void Inv(int* a, int* b, int len) {
    if (len == 1) {
        b[0] = Pow(a[0], MOD - 2);
        return;
    }
    Inv(a, b, len >> 1);
    rep(i, len) aa[i] = a[i], bb[i] = b[i];
    For(i, len, (len << 1)) aa[i] = bb[i] = cc[i] = 0;
    getr(len << 1), multi(aa, bb, cc, len << 1), multi(cc, bb, cc, len << 1);
    rep(i, len) b[i] = ad(ad(bb[i], bb[i]), MOD - cc[i]);
}

int n, nn;
int a[maxn << 2], b[maxn << 2];

int main() {
    n = read();
    rep(i, n) a[i] = read();
    for (nn = 1; nn <= n; nn <<= 1);
    Inv(a, b, nn);
    rep(i, n) printf("%d ", b[i]);
    return 0;
}
```

## 分治FFT

模数:998244353，~~于是拿多项式求逆水过去了~~

```cpp
#include <cstdio>
#include <algorithm>
#include <iostream>
using namespace std;

const int Q=500010,MOD=998244353;

inline int add(int a,int b){a+=b;return a<MOD?a:a-MOD;}
inline int sub(int a,int b){a-=b;return a<0?a+MOD:a;}
inline int mul(int a,int b){return 1LL*a*b%MOD;}

inline int ksm(int a,int b){
    int ans=1;
    while(b){
        if(b&1)ans=mul(ans,a);
        b>>=1,a=mul(a,a);
    }
    return ans;
}


int w[Q];

void NTT(int a[],int n,int f){
    register int i,j,k,m,cnt=1;
    for(i=j=0;i<n;i++){
        if(i<j)swap(a[i],a[j]);
        for(k=n>>1;(j^=k)<k;k>>=1);
    }
    w[0]=1;
    for(m=1;m<n;m<<=1,++cnt){
        int _now=ksm(3,(f*((MOD-1)>>cnt)+MOD-1)%(MOD-1));
        for(i=1;i<m;i++)w[i]=mul(w[i-1],_now);
        for(i=0;i<n;i+=(m<<1))
            for(j=0;j<m;j++){
                int p=a[i+j],q=mul(a[i+j+m],w[j]);
                a[i+j]=add(p,q),a[i+j+m]=sub(p,q);
            }
    }
    if(f<0){
        int n_ni=ksm(n,MOD-2);
        for(i=0;i<n;i++)a[i]=mul(a[i],n_ni);
    }
}

int it[Q];

void GetInv(int a[],int ni[],int toap) { 
    int i,now; 
    ni[0]=ksm(a[0],MOD-2); 
    for(now=2;now<=toap;now<<=1){ 
        fill(ni+(now>>1),ni+(now<<1),0); 
        copy(a,a+now,it); 
        fill(it+now,it+(now<<1),0); 
        NTT(ni,now<<1,1); 
        NTT(it,now<<1,1); 
        for(i=0;i<(now<<1);i++) 
            ni[i]=mul(ni[i],sub(2,mul(it[i],ni[i]))); 
        NTT(ni,now<<1,-1); 
    } 
} 

int b[Q],b_ni[Q];

int main(){
    int p,now;
    scanf("%d",&p);
    for(int i=1;i<p;++i) scanf("%d",b+i);
    for(now=1;now<=p;now<<=1);
    for(int i=1;i<p;i++) b[i]=MOD-b[i];(++b[0])%=MOD;
    GetInv(b,b_ni,now);
    for(int i=0;i<p;++i) printf("%d ",b_ni[i]%MOD);
}
```

## 多项式除法-NTT

模数:998244353

```cpp
#include <cstdio>
#include <algorithm>
#include <iostream>
using namespace std;

const int Q=500010,MOD=998244353;

inline int add(int a,int b){a+=b;return a<MOD?a:a-MOD;}
inline int sub(int a,int b){a-=b;return a<0?a+MOD:a;}
inline int mul(int a,int b){return 1LL*a*b%MOD;}

inline int ksm(int a,int b){
    int ans=1;
    while(b){
        if(b&1)ans=mul(ans,a);
        b>>=1,a=mul(a,a);
    }
    return ans;
}


int w[Q];

void NTT(int a[],int n,int f){
    register int i,j,k,m,cnt=1;
    for(i=j=0;i<n;i++){
        if(i<j)swap(a[i],a[j]);
        for(k=n>>1;(j^=k)<k;k>>=1);
    }
    w[0]=1;
    for(m=1;m<n;m<<=1,++cnt){
        int _now=ksm(3,(f*((MOD-1)>>cnt)+MOD-1)%(MOD-1));
        for(i=1;i<m;i++)w[i]=mul(w[i-1],_now);
        for(i=0;i<n;i+=(m<<1))
        for(j=0;j<m;j++){
            int p=a[i+j],q=mul(a[i+j+m],w[j]);
            a[i+j]=add(p,q),a[i+j+m]=sub(p,q);
        }
    }
    if(f<0){
        int n_ni=ksm(n,MOD-2);
        for(i=0;i<n;i++)a[i]=mul(a[i],n_ni);
    }
}

int it[Q];

void GetInv(int a[],int ni[],int toap) { 
    int i,now; 
    ni[0]=ksm(a[0],MOD-2); 
    for(now=2;now<=toap;now<<=1){ 
        fill(ni+(now>>1),ni+(now<<1),0); 
        copy(a,a+now,it); 
        fill(it+now,it+(now<<1),0); 
        NTT(ni,now<<1,1); 
        NTT(it,now<<1,1); 
        for(i=0;i<(now<<1);i++) ni[i]=mul(ni[i],sub(2,mul(it[i],ni[i]))); 
        NTT(ni,now<<1,-1); 
    } 
} 

int a[Q],b[Q],c[Q],r[Q],b_ni[Q];

int main(){
    int p,q,i,now;
    scanf("%d%d",&p,&q);
    for(i=0;i<=p;i++)scanf("%d",&a[i]);
    for(i=0;i<=q;i++)scanf("%d",&b[i]);
    reverse(a,a+p+1),reverse(b,b+q+1);
    for(now=1;now<=p-q+1;now<<=1);
    GetInv(b,b_ni,now);
    fill(b_ni+now,b_ni+(now<<1),0);
    copy(a,a+now,c);
    fill(c+now,c+(now<<1),0);
    NTT(c,now<<1,1),NTT(b_ni,now<<1,1);
    for(int i=0;i<(now<<1);i++)c[i]=mul(c[i],b_ni[i]);
    NTT(c,now<<1,-1);
    fill(c+p-q+1,c+(now<<1),0);
    reverse(a,a+p+1);reverse(b,b+q+1);reverse(c,c+p-q+1);
    for(i=0;i<=p-q;i++)printf("%d ",c[i]);
    putchar('\n');
    for(now=1;now<=p;now<<=1);
    NTT(c,now,1),NTT(b,now,1),NTT(a,now,1);
    for(i=0;i<now;i++)r[i]=sub(a[i],mul(b[i],c[i]));
    NTT(r,now,-1);for(;r[q-1]==0;--q);
    for(i=0;i<q;i++)printf("%d ",r[i]);
    return 0;
}
```

多项式开根-NTT

模数:988244353

```cpp
#include<cstdio>
#include<cstring>
#include<algorithm>
using namespace std;

const int maxn=500010,p=998244353,g=3,inv_2=499122177;

int qpow(int a,int b,int p){
    int ans=1;
    for(;b;b>>=1,a=(long long)a*a%p)if(b&1)ans=(long long)ans*a%p;
    return ans;
}

void NTT(int *aa,int n,int tp){
    for(int i=1,j=0,k;i<n-1;i++){
        k=n;
        do j^=(k>>=1);while(j<k);
        if(i<j)swap(aa[i],aa[j]);
    }
    for(int k=2;k<=n;k<<=1){
        int wn=qpow(g,(tp>0?(p-1)/k:(p-1)/k*(long long)(p-2)%(p-1)),p);
        for(int i=0;i<n;i+=k){
            int w=1;
            for(int j=0;j<(k>>1);j++,w=(long long)w*wn%p){
                int a=aa[i+j],b=(long long)w*aa[i+j+(k>>1)]%p;
                aa[i+j]=(a+b)%p;
                aa[i+j+(k>>1)]=(a-b+p)%p;
            }
        }
    }
    if(tp<0){
        int inv=qpow(n,p-2,p);
        for(int i=0;i<n;i++)aa[i]=(long long)aa[i]*inv%p;
    }
}

void inv(int *aa,int *cc,int n){
    static int bb[maxn];
    fill(cc,cc+(n<<1),0);
    cc[0]=qpow(aa[0],p-2,p);
    for(int k=2;k<=n;k<<=1){
        copy(aa,aa+k,bb);
        fill(bb+k,bb+(k<<1),0);
        NTT(bb,k<<1,1);
        NTT(cc,k<<1,1);
        for(int i=0;i<(k<<1);i++)cc[i]=cc[i]*((2-(long long)bb[i]*cc[i]%p+p)%p)%p;
        NTT(cc,k<<1,-1);
        fill(cc+k,cc+(k<<1),0);
    }
}

void sqt(int *aa,int *cc,int n){
    static int bb[maxn],D[maxn];
    fill(cc,cc+(n<<1),0);
    cc[0]=1;
    for(int k=2;k<=n;k<<=1){
        copy(aa,aa+k,bb);
        fill(bb+k,bb+(k<<1),0);
        inv(cc,D,k);
        NTT(bb,k<<1,1);
        NTT(D,k<<1,1);
        for(int i=0;i<(k<<1);i++)bb[i]=(long long)bb[i]*D[i]%p;
        NTT(bb,k<<1,-1);
        for(int i=0;i<k;i++)cc[i]=(cc[i]+bb[i])%p*(long long)inv_2%p;
    }
}

int n,N,x,aa[maxn],bb[maxn];

int main(){
    scanf("%d",&n);
    for(N=1;N<=n;N<<=1);
    for(int i=0;i<n;++i) scanf("%d",aa+i);
    sqt(aa,bb,N);
    for(int i=0;i<n;i++)printf("%d ",bb[i]);
    return 0;
}
```

## 多项式对数函数

留坑。。

## 多项式指数函数

留坑。。
