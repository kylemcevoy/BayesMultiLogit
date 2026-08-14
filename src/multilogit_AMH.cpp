#include "RcppArmadillo.h"
#include "numerical_utils.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

namespace {

double log_add_exp(const double a, const double b) {
  const double maximum = std::max(a, b);
  return maximum + std::log1p(std::exp(-std::abs(a - b)));
}

} // namespace

//' Adaptive Metropolis-Hastings Multinomial Logistic Regression
//'
//' Samples a Bayesian multinomial logistic regression posterior using
//' coordinate-wise random-walk Metropolis-Hastings proposals. Proposal scales
//' adapt during burn-in and remain fixed while retained draws are generated.
//' The selected reference category is held at zero.
//'
//' The implementation caches the linear predictors and uses log-sum-exp
//' calculations throughout. Coefficients for each non-reference category have
//' a common multivariate normal prior.
//'
//' @param Y An N by C matrix of non-negative category counts with a positive
//' row total.
//' @param X An N by P numeric design matrix.
//' @param n_sample Positive number of retained draws.
//' @param n_burn Non-negative number of burn-in draws.
//' @param n_sigma_check Positive adaptation-window length.
//' @param step_size Positive initial random-walk standard deviation.
//' @param acceptance_lower Lower target acceptance rate.
//' @param acceptance_upper Upper target acceptance rate.
//' @param step_increase Multiplicative scale increase greater than one.
//' @param step_decrease Multiplicative scale decrease between zero and one.
//' @param prior_mean Length-P normal-prior mean.
//' @param prior_var P by P positive-definite normal-prior covariance matrix.
//' @param reference_cat Zero-based reference-category index. Most users should
//' call `multilogit_AMH()` instead.
//' @param probs If TRUE, return fitted probabilities for retained draws.
//' @param progress If TRUE, print progress every 1,000 iterations.
//'
//' @return A list containing `posterior_coef`, `acceptance_rate`, and
//' `final_step_size`; when requested, it also contains `posterior_prob`.
// [[Rcpp::export]]
List multilogit_AMH_C(
    arma::mat const& Y,
    arma::mat const& X,
    size_t n_sample,
    size_t n_burn,
    size_t n_sigma_check,
    double step_size,
    double acceptance_lower,
    double acceptance_upper,
    double step_increase,
    double step_decrease,
    arma::vec const& prior_mean,
    arma::mat const& prior_var,
    size_t reference_cat,
    bool probs,
    bool progress) {

  const size_t n_sub = X.n_rows;
  const size_t n_pred = X.n_cols;
  const size_t n_cat = Y.n_cols;

  if (Y.n_rows != n_sub || n_sub == 0 || n_pred == 0 || n_cat < 2)
    Rcpp::stop("Invalid X or Y dimensions.");
  if (n_sample == 0 || n_sigma_check == 0)
    Rcpp::stop("n_sample and n_sigma_check must be positive.");
  if (prior_mean.n_elem != n_pred || prior_var.n_rows != n_pred ||
      prior_var.n_cols != n_pred)
    Rcpp::stop("Prior dimensions do not match X.");
  if (reference_cat >= n_cat)
    Rcpp::stop("reference_cat is out of range.");

  arma::mat prior_precision;
  if (!arma::inv_sympd(prior_precision, prior_var))
    Rcpp::stop("prior_var must be positive definite.");

  arma::mat beta(n_pred, n_cat, fill::zeros);
  arma::mat eta(n_sub, n_cat, fill::zeros);
  arma::mat candidate_sigma(n_pred, n_cat, fill::value(step_size));
  candidate_sigma.col(reference_cat).zeros();

  arma::mat window_accept(n_pred, n_cat, fill::zeros);
  arma::mat sample_accept(n_pred, n_cat, fill::zeros);
  const arma::mat YtX = Y.t() * X;
  const arma::vec n_obs = sum(Y, 1);

  arma::cube beta_out(n_pred, n_cat, n_sample, fill::zeros);
  arma::cube prob_out;
  if (probs) prob_out.set_size(n_sub, n_cat, n_sample);

  size_t window_size = 0;

  for (size_t iter = 0; iter < n_burn + n_sample; ++iter) {
    if ((iter & 63U) == 0U) Rcpp::checkUserInterrupt();
    if ((iter & 255U) == 0U) eta = X * beta;

    if (progress && ((iter % 1000U) == 0U || iter + 1 == n_burn + n_sample))
      Rcout << "iteration " << iter + 1 << " of " << n_burn + n_sample << "\n";

    for (size_t category = 0; category < n_cat; ++category) {
      if (category == reference_cat) continue;
      if ((category & 15U) == 0U) Rcpp::checkUserInterrupt();

      // Other categories do not change during this category's coordinate sweep.
      arma::vec log_other(n_sub);
      arma::vec log_denom(n_sub);
      for (size_t row = 0; row < n_sub; ++row)
        log_other(row) = row_log_sum_exp_excluding(eta.row(row), category);
      for (size_t row = 0; row < n_sub; ++row)
        log_denom(row) = log_add_exp(log_other(row), eta(row, category));

      arma::vec centered = beta.col(category) - prior_mean;

      for (size_t predictor = 0; predictor < n_pred; ++predictor) {
        const double current = beta(predictor, category);
        const double proposal = current +
          R::rnorm(0.0, candidate_sigma(predictor, category));
        const double delta = proposal - current;
        if (!std::isfinite(proposal) || !std::isfinite(delta)) continue;

        arma::vec eta_proposal = eta.col(category) + X.col(predictor) * delta;
        if (!eta_proposal.is_finite()) continue;
        arma::vec proposal_denom(n_sub);
        double log_likelihood_ratio = YtX(category, predictor) * delta;

        for (size_t row = 0; row < n_sub; ++row) {
          proposal_denom(row) =
            log_add_exp(log_other(row), eta_proposal(row));
          log_likelihood_ratio += n_obs(row) *
            (log_denom(row) - proposal_denom(row));
        }

        const double precision_cross =
          dot(prior_precision.row(predictor), centered);
        const double log_prior_ratio = -0.5 *
          (2.0 * delta * precision_cross +
           delta * delta * prior_precision(predictor, predictor));

        if (std::log(R::runif(0.0, 1.0)) <
            log_likelihood_ratio + log_prior_ratio) {
          beta(predictor, category) = proposal;
          eta.col(category) = eta_proposal;
          centered(predictor) += delta;
          log_denom = proposal_denom;
          if (iter < n_burn)
            window_accept(predictor, category) += 1.0;
          else
            sample_accept(predictor, category) += 1.0;
        }
      }
    }

    if (iter < n_burn) {
      ++window_size;
      const bool end_window = window_size == n_sigma_check;
      const bool end_burn = iter + 1 == n_burn;
      if (end_window || end_burn) {
        for (size_t category = 0; category < n_cat; ++category) {
          if (category == reference_cat) continue;
          for (size_t predictor = 0; predictor < n_pred; ++predictor) {
            const double rate = window_accept(predictor, category) /
              static_cast<double>(window_size);
            if (rate > acceptance_upper)
              candidate_sigma(predictor, category) *= step_increase;
            else if (rate < acceptance_lower)
              candidate_sigma(predictor, category) *= step_decrease;
            if (!std::isfinite(candidate_sigma(predictor, category)) ||
                candidate_sigma(predictor, category) <= 0.0)
              Rcpp::stop("A proposal scale became non-finite during adaptation.");
          }
        }
        window_accept.zeros();
        window_size = 0;
      }
    } else {
      const size_t draw = iter - n_burn;
      beta_out.slice(draw) = beta;
      if (probs) {
        for (size_t row = 0; row < n_sub; ++row)
          prob_out.slice(draw).row(row) = softmax(eta.row(row));
      }
    }
  }

  arma::mat acceptance_rate = sample_accept /
    static_cast<double>(n_sample);
  acceptance_rate.col(reference_cat).fill(NA_REAL);
  candidate_sigma.col(reference_cat).fill(NA_REAL);

  if (probs) {
    return List::create(
      _["posterior_prob"] = prob_out,
      _["posterior_coef"] = beta_out,
      _["acceptance_rate"] = acceptance_rate,
      _["final_step_size"] = candidate_sigma
    );
  }

  return List::create(
    _["posterior_coef"] = beta_out,
    _["acceptance_rate"] = acceptance_rate,
    _["final_step_size"] = candidate_sigma
  );
}
