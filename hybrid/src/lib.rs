use pyo3::prelude::*;
use rand::Rng;
use std::cmp::Ordering;
use std::f64::consts::{E, PI};
use std::f64::{consts::LN_2, INFINITY};
const CONST: f64 = 1.01; // Placeholder for `const` from original Python code
const SIGMA_E: f64 = 3.19;

use nrfind::find_root;

extern crate num;

use crate::num::ToPrimitive;
use num::bigint::BigInt;
use num::bigint::ToBigInt;
use num::traits::One;
use num::traits::Zero;

use roots::{find_root_brent, SimpleConvergency};
use rug::Float;
use std::ops::DivAssign;

use gomez::nalgebra as na;
use gomez::{Domain, Problem, SolverDriver, System};
use na::{Dyn, IsContiguous};

fn numerical_logq_hybrid(n: i32, l: i32, h: f64) -> f64 {
    let mut initial_guess = 700;
    let mut res = numerical_logq_hybrid_runoptimize(n, l, 3.19, h, initial_guess);

    let mut bound_trials = 200;
    while res.0.is_empty() && bound_trials > 0 {
        bound_trials -= 1;
        initial_guess = rand::thread_rng().gen_range(30..=2500);
        res = numerical_logq_hybrid_runoptimize(n, l, 3.19, h, initial_guess);
        if !res.0.is_empty() {
            break;
        }
    }

    if res.0.is_empty() {
        println!(
            "numerical_logq_hybrid_runoptimize couldn't find a solution after {} trial(s) for this parameters set",
            200 - bound_trials
        );
        return f64::INFINITY;
    }

    let initial_points: Vec<[f64; 3]> = res.1.iter().map(|el| [el[1], el[2], el[3]]).collect();

    let best_logq = check_candidates_logq(&res.0, h, l as f64, n, 1e2, &initial_points);
    best_logq
}

// === Helper Functions ===

fn log2(x: f64) -> f64 {
    x.log2()
}

fn log_delta(beta: f64) -> f64 {
    if beta <= 1.0 {
        return f64::NEG_INFINITY; // Prevent division by zero or log of negative
    }

    let term1 = (1.0 / (2.0 * (beta - 1.0))) * log2(beta / (2.0 * PI * E));
    let term2 = (1.0 / (2.0 * (beta - 1.0))) * (1.0 / beta) * log2(PI * beta);

    term1 + term2
}

fn approx_binom(n: f64, k: f64) -> f64 {
    if k <= 0.0 || k >= n {
        return 0.0;
    }

    let h = entropy(k / n);
    n * h + 0.5 * log2(n / (8.0 * k * (n - k)))
}

fn entropy(x: f64) -> f64 {
    if x <= 0.0 || x >= 1.0 {
        return 0.0;
    }

    -x * log2(x) - (1.0 - x) * log2(1.0 - x)
}

fn probability_enum(n: i32, h: f64, ng: f64, w: f64) -> f64 {
    let h = h as i32;
    let ng = ng as i32;
    let w = w as i32;

    // Use BigFloat for arbitrary-precision floating-point calculations
    let mut prob = Float::with_val(128, 0.0);

    for i in 0..w {
        let numerator_bigint =
            binomial((n - h) as u64, (ng - i) as u64) * binomial(h as u64, i as u64);
        let denominator_bigint = binomial(n as u64, ng as u64);

        // Convert BigInt to f64 directly
        let numerator = Float::with_val(128, numerator_bigint.to_f64().unwrap_or(0.0));
        let denominator = Float::with_val(128, denominator_bigint.to_f64().unwrap_or(1.0));

        if denominator.is_zero() {
            continue;
        }

        let mut term = numerator;
        term.div_assign(denominator); // term = numerator / denominator
        prob += term;
    }

    if prob <= Float::with_val(128, 0.0) {
        f64::NEG_INFINITY // avoid log2(0)
    } else {
        // Compute log2 of prob
        prob.to_f64().log2()
    }
}

//compute binomial coefficients using bigint to avoid overflow
fn binomial(n: u64, k: u64) -> BigInt {
    if k > n {
        return BigInt::zero();
    }

    let mut num = BigInt::one();
    let mut denom = BigInt::one();

    for i in 0..k {
        num *= (n - i).to_bigint().unwrap();
        denom *= (i + 1).to_bigint().unwrap();
    }

    num / denom
}

fn check_candidates_logq(
    logqs: &[f64],
    h: f64,
    l: f64,
    n: i32,
    best_diff: f64,
    initial_guesses: &[[f64; 3]],
) -> f64 {
    let diff_tolerance = 110.0;

    if best_diff >= diff_tolerance {
        println!("best_diff exceeds < {}. Aborting.", diff_tolerance);
        return f64::INFINITY;
    }

    let sigma_e = 3.19;
    let mut best_logq: Option<f64> = None;
    let mut best_diff_mut = best_diff;

    if logqs.len() == 1 {
        return logqs[0];
    }

    for (i, &logq) in logqs.iter().enumerate() {
        let l_ = numerical_lambda_hybrid_v2(n as f64, logq, sigma_e, h, Some(initial_guesses[i]));

        let diff = (l_ - l).abs();
        if diff < best_diff_mut {
            best_logq = Some(logq);
            best_diff_mut = diff;
        }
    }

    match best_logq {
        Some(val) => val,
        None => {
            let new_diff = best_diff + 10.0;
            println!("No solution found setting best_diff to {}", new_diff);
            check_candidates_logq(logqs, h, l, n, new_diff, initial_guesses)
        }
    }
}

fn _delta(beta: f64) -> f64 {
    if beta <= 1.0 {
        return f64::INFINITY; // prevent division by zero or negative root
    }

    let base = beta / (2.0 * PI * E) * (PI * beta).powf(1.0 / beta);
    base.powf(1.0 / (2.0 * (beta - 1.0)))
}

// Define the system for numerical_lambda_hybrid_v2
struct NumericalLambdaHybridV2 {
    n: f64,
    logq: f64,
    sigma_e: f64,
    xi: f64,
    wg: f64,
}

impl Problem for NumericalLambdaHybridV2 {
    type Field = f64;

    fn domain(&self) -> Domain<Self::Field> {
        Domain::unconstrained(3) // 3 variables: ng_, beta_, d_
    }
}

impl System for NumericalLambdaHybridV2 {
    fn eval<Sx, Srx>(
        &self,
        x: &na::Vector<Self::Field, Dyn, Sx>,
        rx: &mut na::Vector<Self::Field, Dyn, Srx>,
    ) where
        Sx: na::storage::Storage<Self::Field, Dyn> + IsContiguous,
        Srx: na::storage::StorageMut<Self::Field, Dyn>,
    {
        let ng_ = x[0];
        let beta_ = x[1];
        let d_ = x[2];

        let delta = _delta(beta_);
        let lnq = self.logq * f64::ln(2.0);

        rx[0] =
            approx_binom(ng_, self.wg) + self.wg + d_.log2() - (0.292 * beta_ + 16.4 + 3.0) + 2.0;
        rx[1] = d_ - (self.n * lnq / delta.ln()).sqrt().ceil() + ng_ - 1.0;
        rx[2] = (-d_ + 1.0) * delta.log2()
            + ((d_ - self.n + ng_ - 1.0) * self.logq + (self.n - ng_) * self.xi.log2()) / d_
            - (2.0 * self.sigma_e * self.sigma_e).log2();
    }
}

#[pyfunction]
fn numerical_lambda_hybrid_v2(
    n: f64,
    logq: f64,
    sigma_e: f64,
    h: f64,
    mut initial_guess: Option<[f64; 3]>,
) -> f64 {
    let lnq = logq * f64::ln(2.0);
    let sigma_s = (h / n).sqrt();
    let xi = sigma_e / sigma_s;
    let mut rt_min = f64::INFINITY;

    for wg in 2..52 {
        let system = NumericalLambdaHybridV2 {
            n,
            logq,
            sigma_e,
            xi,
            wg: wg as f64,
        };

        let mut solver = SolverDriver::builder(&system)
            .with_initial(initial_guess.unwrap_or([n / 4.0, 40.0, 1377.6]).to_vec())
            .build();

        let tolerance = 1e-4;

        let result = solver
            .find(|state| state.norm() <= tolerance || state.iter() >= 50)
            .map(|state| state.0.to_vec())
            .ok();

        if let Some(res) = result {
            let sol_tol = {
                let [f1, f2, f3] = [
                    approx_binom(res[0], wg as f64) + wg as f64 + res[2].log2()
                        - (0.292 * res[1] + 16.4 + 3.0)
                        + 2.0,
                    res[2] - (n * lnq / _delta(res[1]).ln()).sqrt().ceil() + res[0] - 1.0,
                    (-res[2] + 1.0) * _delta(res[1]).log2()
                        + ((res[2] - n + res[0] - 1.0) * logq + (n - res[0]) * xi.log2()) / res[2]
                        - (2.0 * sigma_e * sigma_e).log2(),
                ];
                f1.abs() + f2.abs() + f3.abs()
            };

            if sol_tol < 7.0 {
                let rt = 0.292 * res[1] + (8.0 * res[2]).log2() + 16.4
                    - probability_enum(n as i32, h, res[0], wg as f64);
                if rt < rt_min {
                    rt_min = rt;
                }
            }
        }
    }

    rt_min
}

// Define the system for numerical_logq_hybrid_runoptimize
struct NumericalLogqHybridRunOptimize {
    n: f64,
    l: f64,
    sigma_e: f64,
    xi: f64,
    wg: f64,
    h: f64,
}

impl Problem for NumericalLogqHybridRunOptimize {
    type Field = f64;

    fn domain(&self) -> Domain<Self::Field> {
        Domain::unconstrained(4) // 4 variables: ng, beta, d, logq
    }
}

impl System for NumericalLogqHybridRunOptimize {
    fn eval<Sx, Srx>(
        &self,
        x: &na::Vector<Self::Field, Dyn, Sx>,
        rx: &mut na::Vector<Self::Field, Dyn, Srx>,
    ) where
        Sx: na::storage::Storage<Self::Field, Dyn> + IsContiguous,
        Srx: na::storage::StorageMut<Self::Field, Dyn>,
    {
        let ng = x[0];
        let beta = x[1];
        let d = x[2];
        let logq = x[3];

        let delta = log_delta(beta);

        rx[0] = self.l - (0.292 * beta + 16.4 + 3.0 + log2(d))
            + approx_binom(self.n - self.h, ng - self.wg)
            + approx_binom(self.h, self.wg)
            - approx_binom(self.n, ng)
            - 1.0;

        rx[1] = self.l - approx_binom(ng, self.wg) - self.wg - 2.0 * log2(d)
            + approx_binom(self.n - self.h, ng - self.wg)
            + approx_binom(self.h, self.wg)
            - approx_binom(self.n, ng)
            - 12.0;

        rx[2] = d - (self.n * logq / delta).sqrt().ceil() + ng - 1.0;

        rx[3] = (-d + 1.0) * delta
            + ((d - self.n + ng - 1.0) * logq + (self.n - ng) * log2(self.xi)) / d
            - log2(2.0 * self.sigma_e * self.sigma_e);
    }
}

fn numerical_logq_hybrid_runoptimize(
    n: i32,
    l: i32,
    sigma_e: f64,
    h: f64,
    initial_guess: i32,
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let sigma_s = (h / n as f64).sqrt();
    let xi = sigma_e / sigma_s;
    let mut sol_qs = vec![];
    let mut sols = vec![];

    let beta_initial_guess = (l as f64 - log2(8.0 * 4.0 * n as f64) - 16.4) / 0.292;
    let logq_initial_guess = initial_guess as f64;

    let d_initial_guess = (n as f64 * logq_initial_guess * LN_2 / log_delta(beta_initial_guess))
        .sqrt()
        .ceil();

    for wg in 2..55 {
        let system = NumericalLogqHybridRunOptimize {
            n: n as f64,
            l: l as f64,
            sigma_e,
            xi,
            wg: wg as f64,
            h,
        };

        let mut solver = SolverDriver::builder(&system)
            .with_initial(vec![
                n as f64 / 16.0,
                beta_initial_guess,
                d_initial_guess,
                logq_initial_guess,
            ])
            .build();

        let tolerance = 1e-4;

        let result = solver
            .find(|state| state.norm() <= tolerance || state.iter() >= 50)
            .map(|state| state.0.to_vec())
            .ok();

        if let Some(res) = result {
            let sol_tolerance = (res[0] - res[0].round()).abs()
                + (res[1] - res[1].round()).abs()
                + (res[2] - res[2].round()).abs()
                + (res[3] - res[3].round()).abs();

            if sol_tolerance < 2.5 {
                let rt = 0.292 * res[1] + log2(8.0 * res[2]) + 16.4
                    - probability_enum(n, h, res[0], wg as f64);

                if (l as f64 - rt).abs() < 2.5 {
                    sol_qs.push(res[3].round());
                    sols.push(vec![
                        wg as f64,
                        res[0].round(),
                        res[1].round(),
                        res[2].round(),
                        res[3].round(),
                    ]);
                }
            }
        }
    }

    (sol_qs, sols)
}

#[pymodule]
fn hybrid(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(numerical_lambda_hybrid_v2, m)?)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_numerical_logq_hybrid() {
        let result = numerical_logq_hybrid(1024, 80, 64.0);
        assert!(result.is_finite());
    }

    #[test]
    fn test_log2() {
        assert_eq!(log2(8.0), 3.0);
        assert_eq!(log2(1.0), 0.0);
    }

    #[test]
    fn test_log_delta() {
        assert!(log_delta(2.0).is_finite());
        assert_eq!(log_delta(1.0), f64::NEG_INFINITY);
    }

    #[test]
    fn test_approx_binom() {
        assert!(approx_binom(10.0, 5.0).is_finite());
        assert_eq!(approx_binom(10.0, 0.0), 0.0);
    }

    #[test]
    fn test_entropy() {
        assert_eq!(entropy(0.5), 1.0);
        assert_eq!(entropy(0.0), 0.0);
    }

    #[test]
    fn test_probability_enum() {
        let result = probability_enum(10, 5.0, 3.0, 2.0);
        assert!(result.is_finite());
    }

    #[test]
    fn test_binomial() {
        assert_eq!(binomial(5, 2), BigInt::from(10));
        assert_eq!(binomial(5, 0), BigInt::from(1));
    }

    #[test]
    fn test_delta() {
        assert!(matches!(_delta(2.0), x if x.is_finite()));
        assert_eq!(_delta(1.0), f64::INFINITY);
    }

    #[test]
    fn test_numerical_lambda_hybrid_v2() {
        let result = numerical_lambda_hybrid_v2(1024.0, 40.0, 3.19, 64.0, None);
        assert!(result.is_finite());
    }

    // #[test]
    // fn test_numerical_logq_hybrid_runoptimize() {
    //     let (sol_qs, sols) = numerical_logq_hybrid_runoptimize(1024, 80, 3.19, 64.0, 700);
    //     assert!(!sol_qs.is_empty());
    //     assert!(!sols.is_empty());
    //}
}
