#include "Bob.h"
#include <chrono>
#include <fstream>

int N = 7, p = 3, q = 40961, d = 2; // Parameters chosen by trusted party.
int attemp = 100;
bool found = false;
Ring_p target_f_q(N, q, {0});
Ring_p target_f_p(N, p, {0});
int kk = 0;
int main(){
    vector<double> time(500);
    while(kk < 500){
        auto start = chrono::high_resolution_clock::now();
        Ring_p f1(N, q, {-1, 0, 1, 1, -1, 0, 1});
        Ring_p f2(N, q, {-1, 0, 1, 1, -1, 0, 1});
        //Ring_p candidate_1(N, q, {-1, 0, 1, 1, -1, 0, 1});
        //Ring_p candidate_2(N, p, {-1, 0, 1, 1, -1, 0, 1});
        f1.print();
        f2.print();
        (f1*f2).print();
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double, std::milli> duration = end - start;
        cout << "Running time: " << duration.count() << "ms" << endl;
        time[kk] = duration.count();
        kk += 1;
    }
    ofstream fout("time_data.txt");
    for (int i = 0; i < 500; ++i) {
        fout << time[i] << "\n";
    }
    fout.close();
}
