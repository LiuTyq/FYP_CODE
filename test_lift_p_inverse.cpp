#include "Bob.h"
#include <chrono>
#include <fstream>

int N = 7, p = 3, q = 512, d = 2; // Parameters chosen by trusted party.
int attemp = 100;
bool found = false;
Ring_p target_f_q(N, q, {0});
Ring_p target_f_p(N, p, {0});
int kk = 0;
int main(){
     Ring_p f1(N, q, {-1,0,1,1,-1,0,1});
     Ring_p f2(N, q, {488, 12, 250, 131, 191, 417, 48});
     (f1*f2).print();
    //  f1.print();
    //  f1.inverse_ppower().print();
    // Ring_p f2(N, q, {7,6,5,4,3,2,1});
    // (f1*f2).print();
}