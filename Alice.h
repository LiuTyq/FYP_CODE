#pragma once
#include "polynomial.cpp"
#include <random>

pair<Ring_p, vector<int>> generate_f(int N, int p, int d){   // 生成随机多项式 f，其中 f \in \mathcal{T}(d + 1, d)
    std::vector<int> v(N, 0);  // 初始化长度为 N，所有元素为 0

    // 填充 d + 1 个 1
    for (int i = 0; i < d + 1; ++i) {
        v[i] = 1;
    }
    
    // 填充 d 个 -1
    for (int i = d + 1; i < 2 * d + 1; ++i) {
        v[i] = -1;
    }

    // 随机打乱
    std::random_device rd;  
    std::mt19937 g(rd());  
    std::shuffle(v.begin(), v.end(), g);

    return {Ring_p(N, p, v), v};
}

Ring_p generate_g(int N, int p, int d){   // 生成随机多项式 g，其中 g \in \mathcal{T}(d, d)
    std::vector<int> v(N, 0);  // 初始化长度为 N，所有元素为 0

    // 填充 d 个 1
    for (int i = 0; i < d; ++i) {
        v[i] = 1;
    }
    
    // 填充 d 个 -1
    for (int i = d; i < 2 * d ; ++i) {
        v[i] = -1;
    }

    // 随机打乱
    std::random_device rd;  
    std::mt19937 g(rd());  
    std::shuffle(v.begin(), v.end(), g);

    return Ring_p(N, p, v);
}
