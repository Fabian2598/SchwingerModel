#ifndef HMC_TUNER_H
#define HMC_TUNER_H

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>

// ---------------------------------------------------------------------------
// Automatic step-size tuning for HMC.
//
// Two independent pieces, use either or both:
//
//   (1) target_dH(p_acc)  /  rescale_eps(...)
//       One-shot estimate based on  <dH> = C * eps^4  (2nd order integrator)
//       and  P_acc = erfc( sqrt(<dH>/2) ).  Converges in ~20-50 trajectories.
//
//   (2) HMCTuner
//       Nesterov dual averaging on log(eps), driven by the Rao-Blackwellised
//       acceptance  a_i = min(1, exp(-dH_i)).  Robust, no scaling-law
//       assumption, this is what Stan uses for its warm-up.
//
// IMPORTANT: freeze() before the measurement phase.  While eps is still
// adapting the chain is not Markov and detailed balance does not hold.
// ---------------------------------------------------------------------------

namespace hmc {

// Inverse complementary error function by bisection (plenty fast, called once).
inline double erfc_inv(double y) {
    if (y <= 0.0) return 10.0;
    if (y >= 2.0) return -10.0;
    double lo = -10.0, hi = 10.0;
    for (int i = 0; i < 200; ++i) {
        double mid = 0.5 * (lo + hi);
        if (std::erfc(mid) > y) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

// <dH> that corresponds to a desired acceptance rate, assuming dH is
// Gaussian with <dH> = Var(dH)/2 (which follows from <exp(-dH)> = 1).
inline double target_dH(double p_acc) {
    double x = erfc_inv(p_acc);
    return 2.0 * x * x;
}

// New step size from a measured <dH> at step size eps_old.
// order = 2 for leapfrog / Omelyan (<dH> ~ eps^4), 4 for a 4th order
// integrator (<dH> ~ eps^8).
inline double rescale_eps(double eps_old, double dH_measured, double p_target,
                          int order = 2, double max_factor = 2.0) {
    double dH_t = target_dH(p_target);
    dH_measured = std::max(dH_measured, 1e-12);
    double f = std::pow(dH_t / dH_measured, 1.0 / (2.0 * order));
    f = std::min(std::max(f, 1.0 / max_factor), max_factor);  // safety clip
    return eps_old * f;
}

// ---------------------------------------------------------------------------

class HMCTuner {
public:
    // tau        : trajectory length to keep fixed (set <= 0 to let tau float
    //              and tune the bare step size instead)
    // n_steps0   : initial guess for the number of MD steps
    // p_target   : desired acceptance (0.65 is cost-optimal in theory,
    //              0.75-0.80 is safer with clover fermions)
    // n_warmup   : number of trajectories to adapt over, then call freeze()
    HMCTuner(double tau, int n_steps0, double p_target = 0.78, int n_warmup = 300)
        : tau_(tau), p_target_(p_target), n_warmup_(n_warmup),
          eps_(tau > 0.0 ? tau / n_steps0 : 1.0 / n_steps0),
          n_steps_(n_steps0) {
        mu_ = std::log(10.0 * eps_);
        log_eps_ = std::log(eps_);
        log_eps_bar_ = log_eps_;
    }

    // Under MPI, call set_verbose(mpi::rank == 0) so only one rank prints.
    void set_verbose(bool v) { verbose_ = v; }

    // ---- what the integrator should use for the next trajectory ----------
    int    n_steps() const { return n_steps_; }
    double eps()     const { return eps_; }
    bool   adapting() const { return !frozen_; }

    // Feed dH = H_new - H_old after every trajectory of the warm-up.
    void record(double dH) {
        if (frozen_) return;

        // Rao-Blackwellised acceptance: uses the actual dH, not the coin flip,
        // so it carries ~10x more information per trajectory.
        double a = std::min(1.0, std::exp(-dH));
        if (!std::isfinite(a)) a = 0.0;

        ++m_;
        const double gamma = 0.05, t0 = 10.0, kappa = 0.75;
        double w = 1.0 / (m_ + t0);
        Hbar_ = (1.0 - w) * Hbar_ + w * (p_target_ - a);
        log_eps_ = mu_ - std::sqrt(static_cast<double>(m_)) / gamma * Hbar_;
        double eta = std::pow(static_cast<double>(m_), -kappa);
        log_eps_bar_ = eta * log_eps_ + (1.0 - eta) * log_eps_bar_;

        set_eps(std::exp(log_eps_));

        if (m_ >= n_warmup_) freeze();
    }

    // Switch to the averaged step size and stop adapting.
    void freeze() {
        set_eps(std::exp(log_eps_bar_));
        frozen_ = true;
        if (verbose_)
            std::printf("# HMCTuner: frozen at eps = %.6f, n_steps = %d, tau = %.4f\n",
                        eps_, n_steps_, eps_ * n_steps_);
    }

    // Optional: randomise the trajectory length each update (breaks
    // resonances with periodic modes and improves ergodicity).  Call this
    // *before* each trajectory and pass the result to the integrator;
    // n must be drawn independently of the current configuration.
    template <class RNG>
    int random_n_steps(RNG& rng, double spread = 0.25) const {
        int lo = std::max(1, static_cast<int>(std::floor(n_steps_ * (1.0 - spread))));
        int hi = std::max(lo, static_cast<int>(std::ceil(n_steps_ * (1.0 + spread))));
        std::uniform_int_distribution<int> d(lo, hi);
        return d(rng);
    }

private:
    void set_eps(double e) {
        e = std::min(std::max(e, 1e-4), 1.0);   // sanity bounds
        if (tau_ > 0.0) {
            // Keep tau fixed => eps is quantised as tau/n.
            n_steps_ = std::max(1, static_cast<int>(std::lround(tau_ / e)));
            eps_ = tau_ / n_steps_;
        } else {
            eps_ = e;
        }
    }

    double tau_, p_target_;
    int    n_warmup_;
    double eps_;
    int    n_steps_;
    double mu_, log_eps_, log_eps_bar_, Hbar_ = 0.0;
    long   m_ = 0;
    bool   frozen_ = false;
    bool   verbose_ = false;
};

}  // namespace hmc

#endif  // HMC_TUNER_H

// ---------------------------------------------------------------------------
// Usage sketch
// ---------------------------------------------------------------------------
//
//   hmc::HMCTuner tuner(/*tau=*/1.0, /*n_steps0=*/10, /*p_target=*/0.78,
//                       /*n_warmup=*/300);
//
//   for (int i = 0; i < n_therm + n_meas; ++i) {
//       double dH = do_trajectory(tuner.n_steps(), tuner.eps());
//       if (tuner.adapting()) tuner.record(dH);
//       if (i == n_therm - 1 && tuner.adapting()) tuner.freeze();
//       if (i >= n_therm) measure();
//   }
//
// Or the cheap one-shot version if you just want a good starting point:
//
//   double eps = tau / n0, sum = 0.0;
//   for (int i = 0; i < 30; ++i) sum += do_trajectory(std::lround(tau/eps), eps);
//   eps = hmc::rescale_eps(eps, sum / 30.0, 0.78);
//   int n = std::lround(tau / eps);  eps = tau / n;
//
// ---------------------------------------------------------------------------