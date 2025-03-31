#include "Bob.h"
#include <chrono>


int N = 7, p = 3, q = 41, d = 2; // Parameters chosen by trusted party.
int attemp = 100;
bool found = false;
Ring_p target_f_q(N, q, {0});
Ring_p target_f_p(N, p, {0});
int main(){
     auto start = chrono::high_resolution_clock::now();
/*-----------------------------------------
% First, Alice generate f and g %
-----------------------------------------*/
     for(int i = 0; i < attemp && !found; i++){
        //auto[candidate_1, generate_coeffs] = generate_f(N, q, d); // candidate_1 express an element in R_q, and generate_coeffs return a vector.
        Ring_p candidate_1(N, q, {-1, 0, 1, 1, -1, 0, 1});
        Ring_p candidate_2(N, p, {-1, 0, 1, 1, -1, 0, 1}); // candidate_2 using generate_coeffs to create an element in R_p 
        try {
            Ring_p inverse_1 = candidate_1.inverse();
            Ring_p inverse_2 = candidate_2.inverse();

            //-----** successfully generate Alice private key f, such that F_p and F_g exists! **-----//
            target_f_q = candidate_1; 
            target_f_p = candidate_2;
            //-----** successfully generate Alice private key f, such that F_p and F_g exists! **-----//
            found = true;

            /* 
            Remark: My logic is, using the same generate_coeffs to generate f, where target_f_q represents f in R_q,
            and target_f_p represents f in R_p. 
            */
           
        } catch (const invalid_argument& e) {
            // if fail, then continue loop
            continue;
        }

        if(!found){
            cerr << "Can not find a polynomial invertible, try " << attemp << " times" << endl;
        }
        else{
            cout << "Successfully generate f such that F_p and F_q exists, run " << i + 1 << " times." << endl;
        }
}

cout << endl;

//Ring_p target_g = generate_g(N, q, d); // g \in R_q
Ring_p target_g(N, q, {0, -1, -1, 0, 1, 0, 1}); // g \in R_q
Ring_p F_q = target_f_q.inverse(); // inverse f^{-1} in R_q
Ring_p F_p = target_f_p.inverse();


//------------- %% Important!! THE generation of PUBLIC KET h(x) %%-------------
Ring_p h = F_q * target_g;
cout << "THE PUBLIC Key h(x) is: "; h.print();
cout << endl;
//------------- %% Important!! THE generation of PUBLIC KET h(x) %%-------------


/*-----------------------------------------
% Bob chooses his message and random polynomial %
-----------------------------------------*/



// auto [m_in_R, m_in_Rp] = generate_m(N, p); // generate message m(x) \in R, such that m is the center lift of a polynomial in R_p (also denote it as m(x))
// Ring_p random = generate_r(N, q, d);

Ring_p mess(N, q, {1, -1, 1, 1, 0, -1});
Ring_p random(N, q, {-1, 1, 0, 0, 0, -1, 1});

Ring_p cipher_e = (random * h) * p + mess;
cout << "Bob generates his message as a polynomial, and the ciphertext e(x) is :"; cipher_e.transfer_p(q).print(); cout << endl;
Ring_p a = cipher_e * target_f_q;
Ring a_1 = a.center_lift();
cout << "The center-lift of a(x):"; a_1.print();
Ring_p a_2 = a_1.natural_map(p);
cout << "Reduce a(x) modulo p:";a_2.print();
Ring_p a_3 = F_p * a_2;
cout << "F_p(x) * a(x) = "; a_3.print();

Ring a_4 = a_3.center_lift();
cout << "Alice encrypts the message as a(x): "; a_4.print();

auto end = chrono::high_resolution_clock::now();
chrono::duration<double, std::milli> duration = end - start;
cout << "Running time: " << duration.count() << " ms" << endl;
}