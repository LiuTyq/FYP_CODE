#include "polynomial.h"

using namespace std;

// Fast-power Algorithm, compute b =g^A (mod N)
int FPA(int g, int A, int N){
    int a = g, b = 1;
    while(A > 0){
        if (A & 1) b = (b * a) % N;
        a = (a * a) % N;
        A = A / 2;
    }
    return b;
}



// Butterfly transform, no need to use recursion
// inline void butterfly_transform(vector<int>& a) {
//     int n = a.size();
//     for (int i = 1, j = 0; i < n; i++) {
//         int k;
//         for (k = n >> 1; j & k; k >>= 1)
//             j ^= k;
//         j ^= k;
//         if (i < j)
//             swap(a[i], a[j]);
//     }
// }

inline void butterfly_transform(vector<int>& a) {
   vector<int> rev(a.size());
    int n = a.size();
    for (int i = 0; i < n; ++i) {
    //   rev[i] = rev[i >> 1] >> 1;
    //   if (i & 1) {  // 如果最后一位是 1，则翻转成 len/2
    //     rev[i] |= n >> 1;
    //   }
      rev[i] = (rev[i >> 1] >> 1) | ((i & 1) ? (n >> 1) : 0);
    }
    for (int i = 0; i < n; ++i) {
      if (i < rev[i]) {  // 保证每对数只翻转一次
        swap(a[i], a[rev[i]]);
      }
    }
    return;
  }

// inline void butterfly_transform(vector<int>& a) {
//     int n = a.size();
//     if (n <= 1) return;
//     int log_n = 0;
//     while ((1 << log_n) < n) ++log_n;
//     vector<int> rev(n);
//     rev[0] = 0;
//     for (int i = 1; i < n; ++i) {
//         rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (log_n - 1));
//         if (i < rev[i]) swap(a[i], a[rev[i]]);
//     }
// }

// inline void butterfly_transform(vector<int>& a) {
//     int n = a.size();
//     for (int i = 1, j = 0; i < n; ++i) {
//         int bit = n >> 1;
//         while (j & bit) bit >>= 1;
//         j ^= bit;
//         if (i < j) swap(a[i], a[j]);
//     }
// }



vector<int> ntt(vector<int> &a, int inv, int p, int primitive_root, int primitive_root_inverse){
    // int n = this->coeffs.size();
    int n = a.size();
    butterfly_transform(a);

    for (int len = 2; len <= n; len <<= 1) {
        int wn = FPA(inv == 1 ? primitive_root : primitive_root_inverse, (p - 1) / len, p); // using FPA to generate n-th primitive root
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

    return a;

}


//--------------------% Ring 类的函数定义 % ----------------------
//--------------------% Ring 类的函数定义 % ----------------------

Ring::Ring(int n, initializer_list<int> c) : RingBase(n, vector<int>(c)) {
    if (coeffs.size() < N) coeffs.resize(N, 0); // 确保长度至少为 N
    reduce(); // 归约
}

// 允许 std::vector<int> 作为输入
Ring::Ring(int n, const vector<int>& c) : RingBase(n, c) {
    if (coeffs.size() < N) coeffs.resize(N, 0);
    reduce();
}


// 加法运算
Ring Ring::operator+(const Ring& other) const {
    vector<int> result(N);
    for (int i = 0; i < N; ++i) {
        result[i] = coeffs[i] + other.coeffs[i];
    }
    return Ring(N, result);
}

// 乘法运算
Ring Ring::operator*(const Ring& other) const { 
    //乘法的代码逻辑是，假设 a(x)=a_0+...+x^{N-1}，b(x)=b_0+...+x^{N-1}，先直接乘出来系数行向量，再利用 x^N=1，将行向量的第 N+1 (x^N) 项到第 2N-1 项 (x^{2N-2}) 换成 x^N
    vector<int> result(2 * N - 1, 0); 
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            result[i + j] += coeffs[i] * other.coeffs[j];
        }
    }
    return Ring(N, result);
}


Ring_p Ring::natural_map(int p) {
    transform(coeffs.begin(), coeffs.end(), coeffs.begin(),
    [p](int x) {return (x % p + p) % p; });
    
    return Ring_p(N, p, coeffs);
}

//--------------------% Ring 类的函数定义 % ----------------------
//--------------------% Ring 类的函数定义 % ----------------------




Ring_p::Ring_p(int n, int p, initializer_list<int> c) : RingBase(n, c), p(p){
    if (coeffs.size() <= N) coeffs.resize(N, 0); // 确保长度至少为 N
    else reduce(); // 归约
    transform(coeffs.begin(), coeffs.end(), coeffs.begin(),[p](int x) {return (x % p + p) % p; });
    if(p == 97){
        primitive_root = 5;
        primitive_root_inverse = 39;
    }

    if(p == 257){
        primitive_root = 3;
        primitive_root_inverse = 86;
    }
    if(p == 7681){
        primitive_root = 17;
        primitive_root_inverse = 2711;
    }
    if(p == 12289){
        primitive_root = 11;
        primitive_root_inverse = 5586;
    }
    if(p == 40961){
        primitive_root = 3;
        primitive_root_inverse = 13654;
    }
    if(p == 998244353){
        primitive_root = 3;
        primitive_root_inverse = 13654;
    }
}

Ring_p::Ring_p(int n, int p, const std::vector<int>& c) : RingBase(n, c), p(p)  {
    if (coeffs.size() <= N) coeffs.resize(N, 0); // 确保长度至少为 N
    else reduce(); // 归约
    transform(coeffs.begin(), coeffs.end(), coeffs.begin(),[p](int x) {return (x % p + p) % p; });
    if(p == 97){
        primitive_root = 5;
        primitive_root_inverse = 39;
    }
    if(p == 257){
        primitive_root = 3;
        primitive_root_inverse = 86;
    }
    if(p == 7681){
        primitive_root = 17;
        primitive_root_inverse = 2711;
    }
    if(p == 12289){
        primitive_root = 11;
        primitive_root_inverse = 5586;
    }
    if(p == 40961){
        primitive_root = 3;
        primitive_root_inverse = 13654;
    }
    if(p == 998244353){
        primitive_root = 3;
        primitive_root_inverse = 332748118;
    }

}



// 加法运算
Ring_p Ring_p::operator+(const Ring_p& other) const {
    vector<int> result(N);
    for (int i = 0; i < N; ++i) {
        result[i] = this->coeffs[i] + other.coeffs[i];
    }   
    transform(result.begin(), result.end(), result.begin(),
              [this](int x) {return (x % p + p) % p; }); // 由于是 Z_p 内的加法，因此需要模去 p！
    
    return Ring_p(N, p, result);
}

Ring_p Ring_p::operator+(const Ring& other) const {
    vector<int> result(N);
    vector<int> other_coeffs = other.getCoeffs();
    for (int i = 0; i < N; ++i) {
        result[i] = this->coeffs[i] + other_coeffs[i];
    }
    transform(result.begin(), result.end(), result.begin(),
              [this](int x) {return (x % p + p) % p; }); // 由于是 Z_p 内的加法，因此需要模去 p！
    
    return Ring_p(N, p, result);
}

Ring_p Ring_p::operator-(const Ring_p& other) const {
    vector<int> result(N);
    for (int i = 0; i < N; ++i) {
        result[i] = this->coeffs[i] - other.coeffs[i];
    }
    transform(result.begin(), result.end(), result.begin(),
              [this](int x) {return (x % p + p) % p; }); // 由于是 Z_p 内的加法，因此需要模去 p！
    
    return Ring_p(N, p, result);
}


Ring_p Ring_p::operator*(const Ring_p& other) const { 
    //乘法的代码逻辑是，假设 a(x)=a_0+...+x^{N-1}，b(x)=b_0+...+x^{N-1}，先直接乘出来系数行向量，再利用 x^N=1，
    //将行向量的第 N+1 (x^N) 项到第 2N-1 项 (x^{2N-2}) 换成 x^N
    vector<int> result(2 * N - 1, 0); 
    vector<int> result_1(2 * N - 1, 0); 
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            result[i + j] += this->coeffs[i] * other.coeffs[j];
        }
    }
    return Ring_p(N, p, result);
}




// Ring_p Ring_p::operator*(const Ring_p& other) const { 
//     // 改良乘法——运用快速傅里叶变换，将复杂度从 O(N^2) 降到 O(NlogN)

//     vector<int> a = this->coeffs;
//     vector<int> b = other.coeffs;
//     int temp = 1;
//     while (temp < a.size() + b.size()) temp <<= 1;
//     a.resize(temp), b.resize(temp);
//     // cout << "a.size(): " << a.size(); cout << ", b.size(): " << b.size() << endl;
//     a = ntt(a, 1, p, primitive_root, primitive_root_inverse);
//     b = ntt(b, 1, p, primitive_root, primitive_root_inverse);
//     // cout << "ntt_a: ";
//     // for(int i = 0;i < temp;i++){
//     //     cout << a[i] << " ";
//     // }
//     // cout << endl;
//     // cout << "ntt_b: ";
//     // for(int i = 0;i < temp;i++){
//     //     cout << b[i] << " ";
//     // }
//     for(int i = 0; i < temp; i++){
//         a[i] = (a[i] * b[i]) % p;
//     }
//     a = ntt(a, -1, p, primitive_root, primitive_root_inverse);
//     return Ring_p(N, p, a);
// }


Ring_p Ring_p::operator*(int scalar) const{
    vector<int> result(N, 0);
    for (int i = 0; i < N; ++i){
        result[i] = (coeffs[i] * scalar) % p;
        if(result[i] < 0) result[i] += p;
    }
    return Ring_p(N, p, result);
}

bool Ring_p::operator!=(const Ring_p& other) const {
    return !(this->coeffs == other.coeffs);
}

pair<Ring_p, Ring_p> Ring_p::divide(const Ring_p& divisor) const{
    vector<int> dividend = this->coeffs;
    vector<int> divisor_coeffs = divisor.coeffs;

    int divisor_deg = divisor.degree();
    int leading = divisor_coeffs[divisor_deg]; //获取除数的最高项系数

    // cout << "leading:" << leading << endl;
    int leading_inv = mod_inverse(leading);
    // cout << "leading_inv:" << leading_inv;

    vector<int> r = dividend; //初始余数为被除数
    Ring_p remainder(N, p, r);

    vector<int> q(N, 0); //初始商为0
    Ring_p quotient(N, p, q); 

    int temp = this->degree();
    
    while (true){
        int r_deg = temp;
        // cout << "r_deg:" << r_deg << endl;
        if (r_deg < divisor_deg) break;

        int delta_deg = r_deg - divisor_deg;
        int t_coeff = (r[r_deg] * leading_inv) % p;  // 商项系数

        // 构造商多项式
        vector<int> term(delta_deg + 1, 0);
        term[delta_deg] = t_coeff;
        Ring_p term_1(N, p, term);

        // 更新商多项式
        quotient = quotient + term_1;


        Ring_p mid = divisor * term_1;
        remainder = remainder - mid;
        r = remainder.coeffs;
        temp = remainder.degree();
    }
    Ring_p Quotient = quotient;
    Ring_p Remainder = remainder;

    return make_pair(Quotient, Remainder);
}


Ring_p Ring_p::transfer_N(int N_1) const{
    vector<int> result(N_1, 0);
    if(N < N_1){
        for(int i = 0; i < N; i++){
            result[i] = this->coeffs[i];
        }
    }
    else {// N >= N_1
        for(int i = 0; i < N - N_1; i++){
            result[i] = (this->coeffs[i] + this->coeffs[N_1 + i]) % p;
        }
        for(int i = N - N_1; i < N_1; i++){
            result[i] = this->coeffs[i];
        }
        
    }
    return Ring_p(N_1, p, result);
}


Ring_p Ring_p::transfer_p(int p_1)const{
    vector<int> result(N, 0);
    if (p <= p_1) return Ring_p(N, p_1, coeffs);
    else {
    for(int i = 0; i < coeffs.size(); i++){
        result[i] = coeffs[i] % p_1;
    }
    return Ring_p(N, p_1, result);
}
}

int Ring_p::mod_inverse(int a) const {
    //a = (a % p + p) % p;
    //cout << "p:" << p << endl;
    for (int i = 1; i < p; ++i) {
        if ((a * i) % p == 1) return i;
    }
    return -1; // 无逆元
}

int Ring_p::output_N() const{
    return N;
}

Ring_p Ring_p::inverse() const{
    vector<int> a_coeffs = this->coeffs;
    a_coeffs.push_back(0);

    vector<int> mod_poly(N + 1, 0); mod_poly[0] = -1; mod_poly[N] = 1; // x^N - 1 的系数表示（模 p 后为 -1 ≡ p-1）


    //---------% 为方便进行带余除法 %----------
    Ring_p g(N + 1, p, mod_poly); Ring_p temp_1 = g; // g \in Z_p[x] / (x^{N+1} + 1)，g = x^N - 1;
    Ring_p y(N + 1, p, a_coeffs); Ring_p temp_2 = y; // y \in Z_p[x] / (x^{N+1} + 1)
    //---------% 为方便进行带余除法 %----------


    auto [q, t] = g.divide(y); // 设当前多项式为 f，计算 x^N - 1 = q*f + y，返回 q 和 y
    q = q.transfer_N(N);
    t = t.transfer_N(N);
    g = g.transfer_N(N); // transfer 后，g = 0 
    y = y.transfer_N(N);

    Ring_p u(N, p, {1});
    Ring_p x(N, p, {0});


    Ring_p s = u - q*x;
    u = x;
    g = y;
    x = s;
    y = t;
    // cout << "g:"; g.print(); 
    // cout << "y:"; y.print(); cout << endl;
    // cout << "-----------" << endl;

    while (y != Ring_p(N, p, {0})) {
        
        auto [q, t] = g.divide(y);  // 执行多项式除法 g = q*y + t

        // cout << "q:"; q.print();
        // cout << "t:"; t.print();
        s = u - q * x;

        // 更新变量
        u = x;
        g = y;
        x = s;
        y = t;
        // cout << "g:"; g.print(); 
        // cout << "y:"; y.print(); cout << endl;
        // cout << "-----------" << endl;
    }

    int zero_count = count(g.coeffs.begin(), g.coeffs.end(), 0); // 统计系数中 0 的个数
    
    
    // 检查 gcd 是否为 1
    if (zero_count != N - 1 || g.coeffs[0] == 0) {
        throw invalid_argument("Inverse does not exist");
    }
    // else{
    //     cout << "exist!" << endl;
    // }


    int new_degree = max(g.degree(), (temp_1*u).degree());

    //---------% 同先前相同，为方便计算，进行 transfer %----------
    g = g.transfer_N(N + new_degree);
    u = u.transfer_N(N + new_degree);
    temp_1 = temp_1.transfer_N(N + new_degree);
    temp_2 = temp_2.transfer_N(N + new_degree);
    auto [m, n] = (g - temp_1 * u).divide(temp_2);
    //---------% 同先前相同，为方便计算，进行 transfer %----------
    
    
    m = m.transfer_N(N);
    m.reduce(); // 应用 x^N -1 归约
    return m * mod_inverse(g.coeffs[0]);
}

Ring_p Ring_p::inverse_ppower()const{ 
    int q = 2;
    int n = 7; // 假设 p = q^n，q 是一个素数
    Ring_p temp = this->transfer_p(q); // 定义 temp 量，表示 temp \in R_q
    Ring_p F_x = temp.inverse(); // 确定 F_x 在 Z_q[x]/(x^N - 1) 中的逆元
    int i = 1;
    while(i < n){
        q = q * 2;
        Ring_p ring_const_2(N, q, {2});
        temp = this->transfer_p(q);
        F_x = F_x.transfer_p(q);
        F_x = F_x * (ring_const_2 - temp * F_x);
        // cout << "q:" << q << endl;
        // cout << "temp: "; temp.print();
        // cout << "F_x: "; F_x.print();
        // cout << "F_x * f(temp) = ";(F_x * temp).print(); // 验证
        i += 1;
        // cout << "----------" << endl;                                                                    
    }
    return F_x;
}



Ring Ring_p::center_lift () const{
        vector<int> lifted_coeffs(coeffs);
        int half_p = (p % 2 == 0) ? (p / 2) : (p / 2 + 1);

        transform(lifted_coeffs.begin(), lifted_coeffs.end(), lifted_coeffs.begin(),
            [=](int coef) {return (coef >= half_p) ? coef - p : coef; });
        //需要注意的是，对于本代码中的Ring_p类，所有系数都是位于[0,p) 中的
        //若系数小于 p/2，那么就不需要操作；若系数大于 p/2，那么就减去 p，使其落在 (-p/2,p/2] 中

        return Ring(N, lifted_coeffs);
    
}
