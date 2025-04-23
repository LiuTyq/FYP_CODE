#include<iostream>
#include<vector>

using namespace std;


// Fast-power Algorithm, compute b =g^A (mod N)
long long FPA(int g, int A, int N){
    long long a = g, b = 1;
    while(A > 0){
        if (A & 1) b = (b * a) % N;
        a = (a * a) % N;
        A = A / 2;
    }
    return b;
}



// Butterfly transform, no need to use recursion
void butterfly_transform(vector<int>& a) {
    int n = a.size();
    for (int i = 1, j = 0; i < n; i++) {
        int k;
        for (k = n >> 1; j & k; k >>= 1)
            j ^= k;
        j ^= k;
        if (i < j)
            swap(a[i], a[j]);
    }
}


// If "inv" is 1, then execute NTT, or execute INTT
void ntt(vector<int>& a, int inv, int p, int primitive_nth_root, int primitive_nth_root_inverse) {
    int n = a.size();
    butterfly_transform(a);

    for (int len = 2; len <= n; len <<= 1) {
        int wn = FPA(inv == 1 ? primitive_nth_root : primitive_nth_root_inverse, (p - 1) / len, p); // using FPA to generate n-th primitive root
        // int wn = inv == 1 ? primitive_nth_root : primitive_nth_root_inverse;
        for (int i = 0; i < n; i += len) {
            int w = 1;
            for (int j = 0; j < len / 2; j++) {
                int u = a[i + j];
                int v = (w * a[i + j + len / 2]) % p;  
                a[i + j] = (u + v) % p;
                a[i + j + len / 2] = (u - v + p) % p;
                w = (w * wn) % p;
            }
        }
    }

    if (inv == -1) {
        int inv_n = FPA(n, p - 2, p);
        for (int& x : a) x = x * inv_n % p;
    }
}

vector<int> multiply(vector<int>& a, vector<int>& b, int p, int primitive_nth_root, int primitive_nth_root_inverse){
    int n = 1;
    while (n < a.size() + b.size()) n <<= 1;

    a.resize(n), b.resize(n);
    cout << "a: ";
    for(int i = 0; i < n; i++){
        cout << a[i] << " ";
    }
    cout << endl;
    cout << "b: ";
    for(int i = 0; i < n; i++){
        cout << b[i] << " ";
    }
    cout << endl;
    ntt(a, 1, p, primitive_nth_root, primitive_nth_root_inverse);
    ntt(b, 1, p, primitive_nth_root, primitive_nth_root_inverse);

    cout << "ntt_a: ";
    for(int i = 0; i < n; i++){
        cout << a[i] << " ";
    }
    cout << endl;
    cout << "ntt_b: ";
    for(int i = 0; i < n; i++){
        cout << b[i] << " ";
    }

    cout << endl;
    // for(int i = 0; i < n; i++){
    //     a[i] = (a[i] * b[i]) % p;
    //     cout << a[i] << " ";
    // }
    
    ntt(a, -1, p, primitive_nth_root, primitive_nth_root_inverse);
    for(int i = 0; i < n; i++){
        a[i] = (a[i] * b[i]) % p;
        cout << a[i] << " ";
    }
    while (!a.empty() && a.back() == 0) a.pop_back();
    return a;

    
}

int main(){
     cout << FPA(998244352, 2, 998244353) << endl;
    // vector<int> a = {1, 2, 3, 4, 5, 6, 7, 8};
    // butterfly_transform(a);

    // vector<int> A = {0, -1, -1, 0, 1, 0, 1};
    // vector<int> B = {37, 2, 40, 21, 31, 26, 8};

    // vector<int> result = multiply(A, B, 97, 5, 39);


    // cout << endl;
    // for(int x:result){
    //     cout << x << " ";
    // }

    // ntt(A, 1, 7681, 3383, 4298);
    
    // for (int x : A){
    //     cout << x << " ";
    // }

    // cout << endl;

    // ntt(A, -1, 7681, 3383, 4298);

    // for (int x : A){
    //     cout << x << " ";
    // }

}