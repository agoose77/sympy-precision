
#include <mpfr.h>
#include <iostream>
#include <cmath>

#pragma cling load("libmpfr")

int odd_even_sign(int x){
    float sign = pow(-1, x);
    return (sign > 0) ? 1 : -1;
}

void R(mpfr_t result, const mpfr_t alpha, const mpfr_t x, const mpfr_t u, const unsigned int k, mpfr_rnd_t rnd){
    if (k == 0){
        mpfr_set_ui(result, 1, rnd);
    }

    mpfr_ptr* terms = new mpfr_ptr[k];

    // k_fact = factorial(k)
    mpfr_t k_fact;
    mpfr_init(k_fact);
    mpfr_fac_ui(k_fact, k, rnd);

    for (int m=0; m<k; m++){
        // m_fact = factorial(m)
        mpfr_t m_fact;
        mpfr_init(m_fact);
        mpfr_fac_ui(m_fact, m, rnd);

        // k_fact_div_m_fact = k_fact / m_fact
        mpfr_t k_fact_div_m_fact;
        mpfr_init(k_fact_div_m_fact);
        mpfr_div(k_fact_div_m_fact, k_fact, m_fact, rnd);

        // sgn = (-1)**m
        int sgn = odd_even_sign(m);

        // alpha_pow_m = alpha**m
        mpfr_t alpha_pow_m;
        mpfr_init(alpha_pow_m);
        mpfr_pow_ui(alpha_pow_m, alpha, m, rnd);

        // one_sub_alpha_pow_k_sub_m_add_one = (1-alpha)**(k-m+1)
        unsigned int k_sub_m_add_1 = k - m + 1;
        mpfr_t one_sub_alpha_pow_k_sub_m_add_one;
        mpfr_init(one_sub_alpha_pow_k_sub_m_add_one);
        mpfr_ui_sub(one_sub_alpha_pow_k_sub_m_add_one, 1, alpha, rnd);
        mpfr_pow_ui(one_sub_alpha_pow_k_sub_m_add_one, one_sub_alpha_pow_k_sub_m_add_one, k_sub_m_add_1, rnd);

        // k_mul_x = k*x
        mpfr_t k_mul_x;
        mpfr_init(k_mul_x);
        mpfr_mul_ui(k_mul_x, x, k, rnd);

        // u_sub_k_mul_x_pow_m = (u-k*x)**m
        mpfr_t u_sub_k_mul_x_pow_m;
        mpfr_init(u_sub_k_mul_x_pow_m);
        mpfr_sub(u_sub_k_mul_x_pow_m, u, k_mul_x, rnd);
        mpfr_pow_ui(u_sub_k_mul_x_pow_m, u_sub_k_mul_x_pow_m, m, rnd);

        mpfr_ptr term = new mpfr_t;
        mpfr_init(term);
        mpfr_set(term, k_fact_div_m_fact, rnd);
        terms[m] = term;

        mpfr_mul_si(terms[m], terms[m], sgn, rnd);
        mpfr_mul(terms[m], terms[m], alpha_pow_m, rnd);
        mpfr_mul(terms[m], terms[m], one_sub_alpha_pow_k_sub_m_add_one, rnd);
        mpfr_mul(terms[m], terms[m], u_sub_k_mul_x_pow_m, rnd);
    }

    mpfr_sum(result, terms, k, rnd);
    for (int j=0; j < k; j++)
    {
        delete terms[j];
    }
    delete[] terms;

    // alpha_pow_k = alpha**k
    mpfr_t alpha_pow_k;
    mpfr_init(alpha_pow_k);
    mpfr_pow_ui(alpha_pow_k, alpha, k, rnd);

    // one_sub_alpha_pow_k = (1-alpha)**k
    mpfr_t one_sub_alpha_pow_k;
    mpfr_init(one_sub_alpha_pow_k);
    mpfr_ui_sub(one_sub_alpha_pow_k, 1, alpha, rnd);
    mpfr_pow_ui(one_sub_alpha_pow_k, one_sub_alpha_pow_k, k, rnd);

    //denom = (mp.factorial(k)*α**k*(1-α)**k)
    mpfr_t denominator;
    mpfr_init(denominator);
    mpfr_mul(denominator, k_fact, alpha_pow_k, rnd);
    mpfr_mul(denominator, denominator, one_sub_alpha_pow_k, rnd);
    //numerator
    // (-1)**k
    int sgn_k = odd_even_sign(k);
    // (u-k*X)**k
    mpfr_t u_sub_k_mul_x_pow_k;
    mpfr_init(u_sub_k_mul_x_pow_k);
    mpfr_mul_ui(u_sub_k_mul_x_pow_k, x, k, rnd);
    mpfr_sub(u_sub_k_mul_x_pow_k, u, u_sub_k_mul_x_pow_k, rnd);
    mpfr_pow_ui(u_sub_k_mul_x_pow_k, u_sub_k_mul_x_pow_k, k, rnd);
    mpfr_t numerator;
    mpfr_init(numerator);
    mpfr_mul_si(numerator, alpha_pow_k, sgn_k, rnd);
    mpfr_mul(numerator, numerator, u_sub_k_mul_x_pow_k, rnd);
    mpfr_add(result, numerator, result, rnd);

    mpfr_div(result, result, denominator, rnd);

}

void n_average_collisions_numerical(mpfr_t result, const mpfr_t alpha, const mpfr_t x, const mpfr_t u, const unsigned int i, const mpfr_rnd_t rnd)
{
    mpfr_ptr* terms = new mpfr_ptr[i+1];
    for (int k=0; k<i+1; k++){
        // 1/alpha
        mpfr_t one_div_alpha;
        mpfr_init(one_div_alpha);
        mpfr_ui_div(one_div_alpha, 1, alpha, rnd);

        // (R(alpha, x, u, k) * exp((alpha*u-k*x)/(1-alpha)) -
        mpfr_t r;
        mpfr_init(r);
        R(r, alpha, x, u, k, rnd);

        // k*x
        mpfr_t k_mul_x;
        mpfr_init(k_mul_x);
        mpfr_mul_ui(k_mul_x, x, k, rnd);

        mpfr_t one_sub_alpha;
        mpfr_init(one_sub_alpha);
        mpfr_ui_sub(one_sub_alpha, 1, alpha, rnd);

        // exp((alpha*u-k*x)/(1-alpha))
        mpfr_t exp_term;
        mpfr_init(exp_term);

        mpfr_mul(exp_term, alpha, u, rnd);
        mpfr_sub(exp_term, exp_term, k_mul_x, rnd);
        mpfr_div(exp_term, exp_term, one_sub_alpha, rnd);
        mpfr_exp(exp_term, exp_term, rnd);

        //R*exp((alpha*u-k*x)/(1-alpha))
        mpfr_mul(exp_term, r, exp_term, rnd);

        //R*exp((alpha*u-k*x)/(1-alpha)) - 1
        mpfr_sub_ui(exp_term, exp_term, 1, rnd);
        mpfr_mul(exp_term, one_div_alpha, exp_term, rnd);
        mpfr_add_ui(exp_term, exp_term, 1, rnd);

        mpfr_ptr term = new mpfr_t;
        terms[k] = term;
        mpfr_init(term);
        mpfr_set(term, exp_term, rnd);

    }
    mpfr_sum(result, terms, i+1, rnd);
    for (int j=0; j < i+1; j++){
        delete terms[j];
    }
    delete[] terms;
}

{
    unsigned int A = 17;

    // Set rounding mode arg
    mpfr_rnd_t rnd = MPFR_RNDN;

    // Set default precision to 40*3 bits ~ 40dps
    mpfr_set_default_prec ((int)(10000*3.33));

    mpfr_t result, alpha, x, u;
    mpfr_inits(result, alpha, x, u, NULL);

    // alpha = ((A-1)/(A+1))**2
    mpfr_set_ui(alpha, A-1, rnd);
    mpfr_div_ui(alpha, alpha, A+1, rnd);
    mpfr_sqr(alpha, alpha, rnd);

    // x = log(1/alpha)
    mpfr_ui_div(x, 1, alpha, rnd);
    mpfr_log(x, x, rnd);

    // u = log(2E6/1)
    mpfr_t k;
    mpfr_init(k);
    mpfr_set_ui(k, 2000000, rnd);
    mpfr_log(u, k, rnd);
    
    n_average_collisions_numerical(result, alpha, x, u, 43, rnd);
    std::cout << mpfr_get_d(result, rnd) << std::endl;
}


