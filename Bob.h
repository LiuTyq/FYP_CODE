#include "Alice.h"


pair<Ring, Ring_p> generate_m(int N, int p){   // 生成随机多项式 m，其中 r \in \mathcal{T}(d, d)
    std::vector<int> v_1(N);  // 初始化长度为 N
    std::vector<int> v_2(N);
    // 随机打乱
    random_device rd;  
    mt19937 gen(rd()); 
    int half_p = (p % 2 == 0) ? (p / 2) : (p / 2 + 1);
    uniform_real_distribution<> dist(-p / 2, p / 2); // 生成范围 [-p/2, p/2]

    // 填充向量
    for (int i = 0; i < N; ++i) {
        v_1[i] = dist(gen);
        if(v_1[i] >= 0 && v_1[i] <= half_p) v_2[i] = v_1[i];
        else if (v_1[i] < 0 && v_1[i] > -half_p) v_2[i] = v_1[i] + p;
    }
    Ring m_1(N, v_1); // v_1 \in (-p/2, p/2]
    Ring_p m_2(N, p, v_2); // v_2 \in (0, p]
    return make_pair(m_1, m_2);

}

Ring_p generate_r(int N, int p, int d){   // 生成随机多项式 r，其中 r \in \mathcal{T}(d, d)
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