#pragma once
#include <cmath>
#include <cstdint>

// Given a target absolute error tolerance `epsilon` and failure probability
// `delta`, returns the minimum number of Monte Carlo samples N such that
//
//     P(| p_hat - p_true | > epsilon) <= delta
//
// via Hoeffding's Inequality for a mean of bounded i.i.d. Bernoulli trials:
//
//     N >= (1 / (2 * epsilon^2)) * ln(2 / delta)
//
// Example: epsilon = 0.01, delta = 0.05 (95% confidence, +/-1% error)
//   -> N ~= 18445 samples, independent of graph size.
inline uint64_t hoeffdingSampleSize(double epsilon, double delta) {
    if (epsilon <= 0.0 || epsilon >= 1.0) {
        throw std::invalid_argument("epsilon must be in (0, 1)");
    }
    if (delta <= 0.0 || delta >= 1.0) {
        throw std::invalid_argument("delta must be in (0, 1)");
    }
    double n = (1.0 / (2.0 * epsilon * epsilon)) * std::log(2.0 / delta);
    return static_cast<uint64_t>(std::ceil(n));
}
