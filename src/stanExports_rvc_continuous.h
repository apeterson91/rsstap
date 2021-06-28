// Generated by rstantools.  Do not edit by hand.

/*
    rsstap is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    rsstap is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with rsstap.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_rvc_continuous_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_rvc_continuous");
    reader.add_event(34, 34, "include", "data/weights.stan");
    reader.add_event(34, 0, "start", "data/weights.stan");
    reader.add_event(36, 2, "end", "data/weights.stan");
    reader.add_event(36, 35, "restart", "model_rvc_continuous");
    reader.add_event(64, 61, "end", "model_rvc_continuous");
    return reader;
}
template <typename T0__, typename T1__, typename T2__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__>::type, Eigen::Dynamic, 1>
pw_gaussian(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& y,
                const Eigen::Matrix<T1__, Eigen::Dynamic, 1>& eta,
                const T2__& sigma, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        current_statement_begin__ = 13;
        return stan::math::promote_scalar<fun_return_scalar_t__>(subtract((-(0.5) * stan::math::log((6.283185307179586232 * sigma))), multiply(0.5, square(divide(subtract(y, eta), sigma)))));
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct pw_gaussian_functor__ {
    template <typename T0__, typename T1__, typename T2__>
        Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__>::type, Eigen::Dynamic, 1>
    operator()(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& y,
                const Eigen::Matrix<T1__, Eigen::Dynamic, 1>& eta,
                const T2__& sigma, std::ostream* pstream__) const {
        return pw_gaussian(y, eta, sigma, pstream__);
    }
};
#include <stan_meta_header.hpp>
class model_rvc_continuous
  : public stan::model::model_base_crtp<model_rvc_continuous> {
private:
        int N;
        int P;
        int L;
        int M;
        int K;
        int w_num;
        int v_num;
        std::vector<int> v;
        vector_d w;
        std::vector<int> u;
        vector_d y;
        matrix_d X_smooth;
        matrix_d X_fixef;
        matrix_d S1;
        matrix_d S2;
        matrix_d D;
        int has_weights;
        vector_d weights;
public:
    model_rvc_continuous(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_rvc_continuous(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_rvc_continuous_namespace::model_rvc_continuous";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 18;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            check_greater_or_equal(function__, "N", N, 0);
            current_statement_begin__ = 19;
            context__.validate_dims("data initialization", "P", "int", context__.to_vec());
            P = int(0);
            vals_i__ = context__.vals_i("P");
            pos__ = 0;
            P = vals_i__[pos__++];
            check_greater_or_equal(function__, "P", P, 0);
            current_statement_begin__ = 20;
            context__.validate_dims("data initialization", "L", "int", context__.to_vec());
            L = int(0);
            vals_i__ = context__.vals_i("L");
            pos__ = 0;
            L = vals_i__[pos__++];
            check_greater_or_equal(function__, "L", L, 0);
            current_statement_begin__ = 21;
            context__.validate_dims("data initialization", "M", "int", context__.to_vec());
            M = int(0);
            vals_i__ = context__.vals_i("M");
            pos__ = 0;
            M = vals_i__[pos__++];
            check_greater_or_equal(function__, "M", M, 0);
            current_statement_begin__ = 22;
            context__.validate_dims("data initialization", "K", "int", context__.to_vec());
            K = int(0);
            vals_i__ = context__.vals_i("K");
            pos__ = 0;
            K = vals_i__[pos__++];
            check_greater_or_equal(function__, "K", K, 0);
            current_statement_begin__ = 23;
            context__.validate_dims("data initialization", "w_num", "int", context__.to_vec());
            w_num = int(0);
            vals_i__ = context__.vals_i("w_num");
            pos__ = 0;
            w_num = vals_i__[pos__++];
            check_greater_or_equal(function__, "w_num", w_num, 0);
            current_statement_begin__ = 24;
            context__.validate_dims("data initialization", "v_num", "int", context__.to_vec());
            v_num = int(0);
            vals_i__ = context__.vals_i("v_num");
            pos__ = 0;
            v_num = vals_i__[pos__++];
            check_greater_or_equal(function__, "v_num", v_num, 0);
            current_statement_begin__ = 25;
            validate_non_negative_index("v", "v_num", v_num);
            context__.validate_dims("data initialization", "v", "int", context__.to_vec(v_num));
            v = std::vector<int>(v_num, int(0));
            vals_i__ = context__.vals_i("v");
            pos__ = 0;
            size_t v_k_0_max__ = v_num;
            for (size_t k_0__ = 0; k_0__ < v_k_0_max__; ++k_0__) {
                v[k_0__] = vals_i__[pos__++];
            }
            size_t v_i_0_max__ = v_num;
            for (size_t i_0__ = 0; i_0__ < v_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "v[i_0__]", v[i_0__], 0);
            }
            current_statement_begin__ = 26;
            validate_non_negative_index("w", "w_num", w_num);
            context__.validate_dims("data initialization", "w", "vector_d", context__.to_vec(w_num));
            w = Eigen::Matrix<double, Eigen::Dynamic, 1>(w_num);
            vals_r__ = context__.vals_r("w");
            pos__ = 0;
            size_t w_j_1_max__ = w_num;
            for (size_t j_1__ = 0; j_1__ < w_j_1_max__; ++j_1__) {
                w(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 27;
            validate_non_negative_index("u", "(N + 1)", (N + 1));
            context__.validate_dims("data initialization", "u", "int", context__.to_vec((N + 1)));
            u = std::vector<int>((N + 1), int(0));
            vals_i__ = context__.vals_i("u");
            pos__ = 0;
            size_t u_k_0_max__ = (N + 1);
            for (size_t k_0__ = 0; k_0__ < u_k_0_max__; ++k_0__) {
                u[k_0__] = vals_i__[pos__++];
            }
            size_t u_i_0_max__ = (N + 1);
            for (size_t i_0__ = 0; i_0__ < u_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "u[i_0__]", u[i_0__], 0);
            }
            current_statement_begin__ = 28;
            validate_non_negative_index("y", "N", N);
            context__.validate_dims("data initialization", "y", "vector_d", context__.to_vec(N));
            y = Eigen::Matrix<double, Eigen::Dynamic, 1>(N);
            vals_r__ = context__.vals_r("y");
            pos__ = 0;
            size_t y_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < y_j_1_max__; ++j_1__) {
                y(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 29;
            validate_non_negative_index("X_smooth", "M", M);
            validate_non_negative_index("X_smooth", "L", L);
            context__.validate_dims("data initialization", "X_smooth", "matrix_d", context__.to_vec(M,L));
            X_smooth = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(M, L);
            vals_r__ = context__.vals_r("X_smooth");
            pos__ = 0;
            size_t X_smooth_j_2_max__ = L;
            size_t X_smooth_j_1_max__ = M;
            for (size_t j_2__ = 0; j_2__ < X_smooth_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < X_smooth_j_1_max__; ++j_1__) {
                    X_smooth(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 30;
            validate_non_negative_index("X_fixef", "N", N);
            validate_non_negative_index("X_fixef", "P", P);
            context__.validate_dims("data initialization", "X_fixef", "matrix_d", context__.to_vec(N,P));
            X_fixef = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(N, P);
            vals_r__ = context__.vals_r("X_fixef");
            pos__ = 0;
            size_t X_fixef_j_2_max__ = P;
            size_t X_fixef_j_1_max__ = N;
            for (size_t j_2__ = 0; j_2__ < X_fixef_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < X_fixef_j_1_max__; ++j_1__) {
                    X_fixef(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 31;
            validate_non_negative_index("S1", "L", L);
            validate_non_negative_index("S1", "L", L);
            context__.validate_dims("data initialization", "S1", "matrix_d", context__.to_vec(L,L));
            S1 = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(L, L);
            vals_r__ = context__.vals_r("S1");
            pos__ = 0;
            size_t S1_j_2_max__ = L;
            size_t S1_j_1_max__ = L;
            for (size_t j_2__ = 0; j_2__ < S1_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < S1_j_1_max__; ++j_1__) {
                    S1(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 32;
            validate_non_negative_index("S2", "L", L);
            validate_non_negative_index("S2", "L", L);
            context__.validate_dims("data initialization", "S2", "matrix_d", context__.to_vec(L,L));
            S2 = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(L, L);
            vals_r__ = context__.vals_r("S2");
            pos__ = 0;
            size_t S2_j_2_max__ = L;
            size_t S2_j_1_max__ = L;
            for (size_t j_2__ = 0; j_2__ < S2_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < S2_j_1_max__; ++j_1__) {
                    S2(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 33;
            validate_non_negative_index("D", "M", M);
            validate_non_negative_index("D", "K", K);
            context__.validate_dims("data initialization", "D", "matrix_d", context__.to_vec(M,K));
            D = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(M, K);
            vals_r__ = context__.vals_r("D");
            pos__ = 0;
            size_t D_j_2_max__ = K;
            size_t D_j_1_max__ = M;
            for (size_t j_2__ = 0; j_2__ < D_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < D_j_1_max__; ++j_1__) {
                    D(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 35;
            context__.validate_dims("data initialization", "has_weights", "int", context__.to_vec());
            has_weights = int(0);
            vals_i__ = context__.vals_i("has_weights");
            pos__ = 0;
            has_weights = vals_i__[pos__++];
            check_greater_or_equal(function__, "has_weights", has_weights, 0);
            check_less_or_equal(function__, "has_weights", has_weights, 1);
            current_statement_begin__ = 36;
            validate_non_negative_index("weights", "(has_weights ? N : 0 )", (has_weights ? N : 0 ));
            context__.validate_dims("data initialization", "weights", "vector_d", context__.to_vec((has_weights ? N : 0 )));
            weights = Eigen::Matrix<double, Eigen::Dynamic, 1>((has_weights ? N : 0 ));
            vals_r__ = context__.vals_r("weights");
            pos__ = 0;
            size_t weights_j_1_max__ = (has_weights ? N : 0 );
            for (size_t j_1__ = 0; j_1__ < weights_j_1_max__; ++j_1__) {
                weights(j_1__) = vals_r__[pos__++];
            }
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 40;
            validate_non_negative_index("beta_fixef", "P", P);
            num_params_r__ += P;
            current_statement_begin__ = 41;
            validate_non_negative_index("beta_smooth", "L", L);
            num_params_r__ += L;
            current_statement_begin__ = 42;
            validate_non_negative_index("alpha", "K", K);
            num_params_r__ += K;
            current_statement_begin__ = 43;
            num_params_r__ += 1;
            current_statement_begin__ = 44;
            num_params_r__ += 1;
            current_statement_begin__ = 45;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_rvc_continuous() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 40;
        if (!(context__.contains_r("beta_fixef")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta_fixef missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta_fixef");
        pos__ = 0U;
        validate_non_negative_index("beta_fixef", "P", P);
        context__.validate_dims("parameter initialization", "beta_fixef", "vector_d", context__.to_vec(P));
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta_fixef(P);
        size_t beta_fixef_j_1_max__ = P;
        for (size_t j_1__ = 0; j_1__ < beta_fixef_j_1_max__; ++j_1__) {
            beta_fixef(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(beta_fixef);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta_fixef: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 41;
        if (!(context__.contains_r("beta_smooth")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta_smooth missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta_smooth");
        pos__ = 0U;
        validate_non_negative_index("beta_smooth", "L", L);
        context__.validate_dims("parameter initialization", "beta_smooth", "vector_d", context__.to_vec(L));
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta_smooth(L);
        size_t beta_smooth_j_1_max__ = L;
        for (size_t j_1__ = 0; j_1__ < beta_smooth_j_1_max__; ++j_1__) {
            beta_smooth(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(beta_smooth);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta_smooth: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 42;
        if (!(context__.contains_r("alpha")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable alpha missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("alpha");
        pos__ = 0U;
        validate_non_negative_index("alpha", "K", K);
        context__.validate_dims("parameter initialization", "alpha", "vector_d", context__.to_vec(K));
        Eigen::Matrix<double, Eigen::Dynamic, 1> alpha(K);
        size_t alpha_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < alpha_j_1_max__; ++j_1__) {
            alpha(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(alpha);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable alpha: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 43;
        if (!(context__.contains_r("sigma")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable sigma missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("sigma");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "sigma", "double", context__.to_vec());
        double sigma(0);
        sigma = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, sigma);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable sigma: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 44;
        if (!(context__.contains_r("tau1")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable tau1 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("tau1");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "tau1", "double", context__.to_vec());
        double tau1(0);
        tau1 = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, tau1);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable tau1: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 45;
        if (!(context__.contains_r("tau2")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable tau2 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("tau2");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "tau2", "double", context__.to_vec());
        double tau2(0);
        tau2 = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, tau2);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable tau2: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 40;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta_fixef;
            (void) beta_fixef;  // dummy to suppress unused var warning
            if (jacobian__)
                beta_fixef = in__.vector_constrain(P, lp__);
            else
                beta_fixef = in__.vector_constrain(P);
            current_statement_begin__ = 41;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta_smooth;
            (void) beta_smooth;  // dummy to suppress unused var warning
            if (jacobian__)
                beta_smooth = in__.vector_constrain(L, lp__);
            else
                beta_smooth = in__.vector_constrain(L);
            current_statement_begin__ = 42;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> alpha;
            (void) alpha;  // dummy to suppress unused var warning
            if (jacobian__)
                alpha = in__.vector_constrain(K, lp__);
            else
                alpha = in__.vector_constrain(K);
            current_statement_begin__ = 43;
            local_scalar_t__ sigma;
            (void) sigma;  // dummy to suppress unused var warning
            if (jacobian__)
                sigma = in__.scalar_lb_constrain(0, lp__);
            else
                sigma = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 44;
            local_scalar_t__ tau1;
            (void) tau1;  // dummy to suppress unused var warning
            if (jacobian__)
                tau1 = in__.scalar_lb_constrain(0, lp__);
            else
                tau1 = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 45;
            local_scalar_t__ tau2;
            (void) tau2;  // dummy to suppress unused var warning
            if (jacobian__)
                tau2 = in__.scalar_lb_constrain(0, lp__);
            else
                tau2 = in__.scalar_lb_constrain(0);
            // model body
            {
            current_statement_begin__ = 49;
            validate_non_negative_index("eta", "M", M);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> eta(M);
            stan::math::initialize(eta, DUMMY_VAR__);
            stan::math::fill(eta, DUMMY_VAR__);
            current_statement_begin__ = 50;
            validate_non_negative_index("eta_", "N", N);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> eta_(N);
            stan::math::initialize(eta_, DUMMY_VAR__);
            stan::math::fill(eta_, DUMMY_VAR__);
            current_statement_begin__ = 51;
            lp_accum__.add(normal_log<propto__>(alpha, 0, 1));
            current_statement_begin__ = 52;
            lp_accum__.add(gamma_log<propto__>(tau1, 1, 1));
            current_statement_begin__ = 53;
            lp_accum__.add(gamma_log<propto__>(tau2, 1, 1));
            current_statement_begin__ = 54;
            lp_accum__.add(cauchy_log<propto__>(sigma, 0, 5));
            current_statement_begin__ = 55;
            lp_accum__.add(multi_normal_prec_log<propto__>(beta_smooth, rep_vector(0, L), add(multiply(S1, tau1), multiply(S2, tau2))));
            current_statement_begin__ = 56;
            stan::math::assign(eta, elt_multiply(stan::math::exp(multiply(D, alpha)), multiply(X_smooth, beta_smooth)));
            current_statement_begin__ = 57;
            stan::math::assign(eta_, add(csr_matrix_times_vector(N, M, w, v, u, eta), multiply(X_fixef, beta_fixef)));
            current_statement_begin__ = 58;
            if (as_bool(has_weights)) {
                current_statement_begin__ = 59;
                lp_accum__.add(dot_product(weights, pw_gaussian(y, eta_, sigma, pstream__)));
            } else {
                current_statement_begin__ = 61;
                lp_accum__.add(normal_log<propto__>(y, eta_, sigma));
            }
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("beta_fixef");
        names__.push_back("beta_smooth");
        names__.push_back("alpha");
        names__.push_back("sigma");
        names__.push_back("tau1");
        names__.push_back("tau2");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(P);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(L);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(K);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_rvc_continuous_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta_fixef = in__.vector_constrain(P);
        size_t beta_fixef_j_1_max__ = P;
        for (size_t j_1__ = 0; j_1__ < beta_fixef_j_1_max__; ++j_1__) {
            vars__.push_back(beta_fixef(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta_smooth = in__.vector_constrain(L);
        size_t beta_smooth_j_1_max__ = L;
        for (size_t j_1__ = 0; j_1__ < beta_smooth_j_1_max__; ++j_1__) {
            vars__.push_back(beta_smooth(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> alpha = in__.vector_constrain(K);
        size_t alpha_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < alpha_j_1_max__; ++j_1__) {
            vars__.push_back(alpha(j_1__));
        }
        double sigma = in__.scalar_lb_constrain(0);
        vars__.push_back(sigma);
        double tau1 = in__.scalar_lb_constrain(0);
        vars__.push_back(tau1);
        double tau2 = in__.scalar_lb_constrain(0);
        vars__.push_back(tau2);
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            if (!include_gqs__ && !include_tparams__) return;
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_rvc_continuous";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t beta_fixef_j_1_max__ = P;
        for (size_t j_1__ = 0; j_1__ < beta_fixef_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta_fixef" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t beta_smooth_j_1_max__ = L;
        for (size_t j_1__ = 0; j_1__ < beta_smooth_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta_smooth" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t alpha_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < alpha_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "alpha" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "tau1";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "tau2";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t beta_fixef_j_1_max__ = P;
        for (size_t j_1__ = 0; j_1__ < beta_fixef_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta_fixef" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t beta_smooth_j_1_max__ = L;
        for (size_t j_1__ = 0; j_1__ < beta_smooth_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta_smooth" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t alpha_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < alpha_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "alpha" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "tau1";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "tau2";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_rvc_continuous_namespace::model_rvc_continuous stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif