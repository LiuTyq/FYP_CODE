#include "polynomial.h"

using namespace std;




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
    transform(coeffs.begin(), coeffs.end(), coeffs.begin(),
              [p](int x) {return (x % p + p) % p; });
}

Ring_p::Ring_p(int n, int p, const std::vector<int>& c) : RingBase(n, c), p(p)  {
    if (coeffs.size() <= N) coeffs.resize(N, 0); // 确保长度至少为 N
    else reduce(); // 归约
    transform(coeffs.begin(), coeffs.end(), coeffs.begin(),
              [p](int x) {return (x % p + p) % p; });
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
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            result[i + j] += this->coeffs[i] * other.coeffs[j];
        }
    }
    return Ring_p(N, p, result);
}

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
    int leading = divisor_coeffs[divisor_deg]; //获取除数的首项系数


    int leading_inv = mod_inverse(leading);

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
    a = (a % p + p) % p;
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
    while (y != Ring_p(N, p, {0})) {
        
        
        auto [q, t] = g.divide(y);  // 执行多项式除法 g = q*y + t

        s = u - q * x;

        // 更新变量
        u = x;
        g = y;
        x = s;
        y = t;
    }

    int zero_count = count(g.coeffs.begin(), g.coeffs.end(), 0); // 统计系数中 0 的个数
    
    
    // 检查 gcd 是否为 1
    if (zero_count != N - 1 || g.coeffs[0] == 0) {
        throw invalid_argument("Inverse does not exist");
    }


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




Ring Ring_p::center_lift () const{
        vector<int> lifted_coeffs(coeffs);
        int half_p = (p % 2 == 0) ? (p / 2) : (p / 2 + 1);

        transform(lifted_coeffs.begin(), lifted_coeffs.end(), lifted_coeffs.begin(),
            [=](int coef) {return (coef >= half_p) ? coef - p : coef; });
        //需要注意的是，对于本代码中的Ring_p类，所有系数都是位于[0,p) 中的
        //若系数小于 p/2，那么就不需要操作；若系数大于 p/2，那么就减去 p，使其落在 (-p/2,p/2] 中

        return Ring(N, lifted_coeffs);
    
}
