#pragma once

#include <iostream>
#include <algorithm>
#include <vector>
#include <initializer_list>
#include <numeric>
#include <cmath>
using namespace std;





template<typename T>
class RingBase {
protected:
    int N;
    vector<T> coeffs;
public:
    // 构造函数
    RingBase(int n, const vector<T>& c) : N(n), coeffs(c) {}

    // 打印函数
    void print() const {
        bool first = true; // 判断 x^i 是不是第一项，如果是的话前面就不需要加 "+" 号了！
        if (coeffs[0] != 0){
            cout << coeffs[0] << " ";
            first = false;
        }
        if (coeffs[1] !=0){
            if (coeffs[1] > 0) {
                if (!first) cout << " + ";
                cout << coeffs[1] << "x";
            }else {
                if (!first) cout << " + ";
                cout << "(" << coeffs[1] << ")"<< "x";
            }
            first = false;
        }
        for (int i = 2; i < N; ++i) {
            if (coeffs[i] != 0) {
                if (coeffs[i] > 0) {
                    if (!first) cout << " + ";
                    cout << coeffs[i] << "x^" << i;
                }else {
                    if (!first) cout << " + ";
                    cout << "(" << coeffs[i] << ")"<< "x^" << i;
                }
                first = false;
            }
        }
        if (first) cout << "0"; // 表示多项式为 0
        cout << endl;
    }


    //归约函数，保证函数是落于 R[x] / (x^N - 1) 中的，R = Z 或者 Z_p 或者 Z_q
    void reduce() { 
        if (coeffs.size() > N) {
            for (size_t i = N; i < coeffs.size(); ++i) {
                coeffs[i % N] += coeffs[i]; // x^N = 1，因此合并系数
            }
            coeffs.resize(N);
        }
    }



};
//------------------------------------------------------模板类





//--------------------% 提前声明 % ----------------------
class Ring_p;
//--------------------% 提前声明 % ----------------------




//--------------------% Ring 类 % ----------------------

class Ring : public RingBase<int> {
public:
    // 构造函数
    Ring(int n, initializer_list<int> c);
    Ring(int n, const vector<int>& c);

    // 运算符重载
    Ring operator+(const Ring& other) const;
    Ring operator*(const Ring& other) const;
    Ring operator-(const Ring& other) const;

    // 自然同构至 R_p
    Ring_p natural_map(int p);



    vector<int> getCoeffs() const {return coeffs;}
};

//--------------------% Ring 类 % ----------------------





//--------------------% Ring_p 类 % ----------------------



class Ring_p : public RingBase<int> {
private:
    int p;
    int primitive_root;        // 原根
    int primitive_root_inverse; // 原根的逆元
    int degree() const {        // 获取多项式的度数
        for (int i = coeffs.size() - 1; i >= 0; --i) {
            if (coeffs[i] != 0) return i;
        }
        return -1; // 零多项式
    }
public:
     // 构造函数
     Ring_p(int n, int p, std::initializer_list<int> c);
     Ring_p(int n, int p, const std::vector<int>& c);
 


     //---------------% 运算符重载 %---------------
     Ring_p operator+(const Ring_p& other) const;
     Ring_p operator+(const Ring& other) const;
     Ring_p operator*(const Ring_p& other) const;
     Ring_p operator*(int scalar) const;
     Ring_p operator-(const Ring_p& other) const;
     pair<Ring_p, Ring_p> divide(const Ring_p& divisor) const;
     bool operator!=(const Ring_p& other) const;
     //---------------% 运算符重载 %---------------




     // 计算某个数 a 模 p 的逆元，当 p 为素数时，算法始终成立
     int mod_inverse(int a) const; 
     int output_N() const;

     Ring_p transfer_N(int N_1) const; // 多项式环的转化
     Ring_p transfer_p(int p_1) const;
     Ring_p inverse() const; // 求取元素逆元
     Ring_p inverse_ppower() const;
     // Ring_p generate_f(int d) const; // 生成随机多项式 m，其中 m \in \mathcal{T}(d + 1, d)
     Ring center_lift() const;// Center-Lift (中心提升)


     
};
//--------------------% Ring_p 类 % ----------------------