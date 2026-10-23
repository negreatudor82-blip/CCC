#ifndef FACTORY_LINEAR_PROGRAM_HPP
#define FACTORY_LINEAR_PROGRAM_HPP

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <vector>

// A self-contained two-phase simplex solver.  The caller must still validate
// integer plans; a floating-point LP result is not an exact certificate.
namespace factory_lp {

enum class LPStatus { Optimal, Infeasible, Unbounded, TimedOut, NumericalFailure };

struct LPSolution {
    LPStatus status = LPStatus::NumericalFailure;
    double objective = 0;
    std::vector<double> x;
    std::size_t pivots = 0;
};

namespace detail {
using Real = long double;
using Clock = std::chrono::steady_clock;

class Simplex {
    static constexpr Real eps = 1e-11L;
    int m, n;
    std::vector<int> basic, nonbasic;
    std::vector<std::vector<Real>> tab;
    Clock::time_point deadline;
    std::size_t pivots = 0;

    bool expired() const { return Clock::now() >= deadline; }

    bool pivot(int r, int s) {
        if (expired()) return false;
        const Real p = tab[r][s];
        if (!std::isfinite(p) || std::fabs(p) < 1e-25L) return false;
        const Real inverse = 1 / p;
        for (int i = 0; i < m + 2; ++i) {
            if (i == r) continue;
            const Real factor = tab[i][s] * inverse;
            if (factor != 0) {
                for (int j = 0; j < n + 2; ++j)
                    if (j != s) tab[i][j] -= tab[r][j] * factor;
            }
        }
        for (int j = 0; j < n + 2; ++j)
            if (j != s) tab[r][j] *= inverse;
        for (int i = 0; i < m + 2; ++i)
            if (i != r) tab[i][s] *= -inverse;
        tab[r][s] = inverse;
        std::swap(basic[r], nonbasic[s]);
        ++pivots;
        return true;
    }

    LPStatus simplex(int phase) {
        const int objective = phase == 1 ? m + 1 : m;
        while (true) {
            if (expired()) return LPStatus::TimedOut;
            // Bland's rule avoids cycling on degenerate zero-RHS networks.
            int entering = -1;
            for (int j = 0; j <= n; ++j) {
                if (phase == 2 && nonbasic[j] == -1) continue;
                if (tab[objective][j] < -eps &&
                    (entering < 0 || nonbasic[j] < nonbasic[entering]))
                    entering = j;
            }
            if (entering < 0) return LPStatus::Optimal;
            int leaving = -1;
            for (int i = 0; i < m; ++i) {
                if (tab[i][entering] <= eps) continue;
                if (leaving < 0) {
                    leaving = i;
                    continue;
                }
                Real lhs = tab[i][n + 1] / tab[i][entering];
                Real rhs = tab[leaving][n + 1] / tab[leaving][entering];
                Real tolerance = eps * (1 + std::fabs(lhs) + std::fabs(rhs));
                if (lhs < rhs - tolerance ||
                    (std::fabs(lhs - rhs) <= tolerance && basic[i] < basic[leaving]))
                    leaving = i;
            }
            if (leaving < 0) return LPStatus::Unbounded;
            if (!pivot(leaving, entering))
                return expired() ? LPStatus::TimedOut : LPStatus::NumericalFailure;
        }
    }

public:
    Simplex(const std::vector<std::vector<double>>& a,
            const std::vector<double>& b, const std::vector<double>& c,
            Clock::time_point limit)
        : m(static_cast<int>(b.size())), n(static_cast<int>(c.size())),
          basic(m), nonbasic(n + 1),
          tab(m + 2, std::vector<Real>(n + 2)), deadline(limit) {
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) tab[i][j] = a[i][j];
            basic[i] = n + i;
            tab[i][n] = -1;
            tab[i][n + 1] = b[i];
        }
        for (int j = 0; j < n; ++j) {
            nonbasic[j] = j;
            tab[m][j] = -c[j];
        }
        nonbasic[n] = -1;
        tab[m + 1][n] = 1;
    }

    LPSolution solve() {
        LPSolution answer;
        if (expired()) {
            answer.status = LPStatus::TimedOut;
            return answer;
        }
        int row = -1;
        for (int i = 0; i < m; ++i)
            if (row < 0 || tab[i][n + 1] < tab[row][n + 1]) row = i;
        if (row >= 0 && tab[row][n + 1] < -eps) {
            if (!pivot(row, n)) {
                answer.status = expired() ? LPStatus::TimedOut : LPStatus::NumericalFailure;
                return answer;
            }
            LPStatus phase1 = simplex(1);
            if (phase1 != LPStatus::Optimal) {
                answer.status = phase1 == LPStatus::Unbounded
                    ? LPStatus::NumericalFailure : phase1;
                answer.pivots = pivots;
                return answer;
            }
            if (tab[m + 1][n + 1] < -eps) {
                answer.status = LPStatus::Infeasible;
                answer.pivots = pivots;
                return answer;
            }
            if (std::fabs(tab[m + 1][n + 1]) > 1e-7L) {
                answer.status = LPStatus::NumericalFailure;
                answer.pivots = pivots;
                return answer;
            }
            // Remove an artificial basic variable if its row is not redundant.
            for (int i = 0; i < m; ++i) {
                if (basic[i] != -1) continue;
                int column = -1;
                for (int j = 0; j <= n; ++j)
                    if (nonbasic[j] != -1 && std::fabs(tab[i][j]) > eps &&
                        (column < 0 || nonbasic[j] < nonbasic[column])) column = j;
                if (column >= 0 && !pivot(i, column)) {
                    answer.status = expired() ? LPStatus::TimedOut : LPStatus::NumericalFailure;
                    answer.pivots = pivots;
                    return answer;
                }
            }
        }
        answer.status = simplex(2);
        answer.pivots = pivots;
        if (answer.status != LPStatus::Optimal) return answer;
        answer.x.assign(n, 0);
        for (int i = 0; i < m; ++i)
            if (basic[i] >= 0 && basic[i] < n)
                answer.x[basic[i]] = static_cast<double>(tab[i][n + 1]);
        answer.objective = static_cast<double>(tab[m][n + 1]);
        return answer;
    }
};
} // namespace detail

// Maximize c*x subject to A*x <= b and x >= 0.
// TimedOut results intentionally have no partial solution or usable bound.
inline LPSolution solveLP(
    const std::vector<std::vector<double>>& a, const std::vector<double>& b,
    const std::vector<double>& c,
    std::chrono::steady_clock::time_point deadline =
        std::chrono::steady_clock::time_point::max()) {
    LPSolution failure;
    if (a.size() != b.size()) return failure;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].size() != c.size() || !std::isfinite(b[i])) return failure;
        for (double value : a[i]) if (!std::isfinite(value)) return failure;
    }
    for (double value : c) if (!std::isfinite(value)) return failure;
    detail::Simplex solver(a, b, c, deadline);
    LPSolution answer = solver.solve();
    if (answer.status != LPStatus::Optimal) return answer;
    for (double& value : answer.x) {
        if (!std::isfinite(value) || value < -1e-7) {
            answer.status = LPStatus::NumericalFailure;
            return answer;
        }
        if (value < 0) value = 0;
    }
    for (std::size_t i = 0; i < a.size(); ++i) {
        long double lhs = 0, magnitude = std::fabs(b[i]);
        for (std::size_t j = 0; j < c.size(); ++j) {
            long double term = static_cast<long double>(a[i][j]) * answer.x[j];
            lhs += term;
            magnitude += std::fabs(term);
        }
        if (lhs > b[i] + 1e-7L * (1 + magnitude)) {
            answer.status = LPStatus::NumericalFailure;
            return answer;
        }
    }
    long double objective = 0;
    for (std::size_t j = 0; j < c.size(); ++j)
        objective += static_cast<long double>(c[j]) * answer.x[j];
    if (!std::isfinite(objective)) answer.status = LPStatus::NumericalFailure;
    else answer.objective = static_cast<double>(objective);
    return answer;
}

} // namespace factory_lp
#endif
