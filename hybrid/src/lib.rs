use pyo3::prelude::*;
use rand::Rng;
use std::cmp::Ordering;
use std::f64::consts::{E, PI};
use std::f64::{consts::LN_2, INFINITY};
const CONST: f64 = 1.01; // Placeholder for `const` from original Python code
const SIGMA_E: f64 = 3.19;

extern crate num;

use crate::num::ToPrimitive;
use num::bigint::BigInt;
use num::bigint::ToBigInt;
use num::traits::One;
use num::traits::Zero;

use roots::{find_root_brent, SimpleConvergency};
use rug::Float;
use std::ops::DivAssign;

#[pyfunction]
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

fn newton_solve<F>(mut x: [f64; 3], f: F, max_iter: usize, tol: f64) -> Option<[f64; 3]>
where
    F: Fn([f64; 3]) -> [f64; 3],
{
    let h = 1e-6;
    for _ in 0..max_iter {
        let fx = f(x);
        let norm_fx: f64 = fx.iter().map(|v| v.abs()).sum();
        if norm_fx < tol {
            // Round the result to avoid floating-point precision issues
            return Some(x.map(|v| (v * 1e9).round() / 1e9));
        }

        // Approximate Jacobian via finite differences
        let mut jacobian = [[0.0; 3]; 3];
        for i in 0..3 {
            let mut x_h = x;
            x_h[i] += h;
            let fx_h = f(x_h);
            for j in 0..3 {
                jacobian[j][i] = (fx_h[j] - fx[j]) / h;
            }
        }

        // Solve J * dx = -f(x) using Gaussian elimination or a linear solver
        let dx = match solve_linear_system(jacobian, fx.map(|v| -v)) {
            Some(sol) => sol,
            None => return None,
        };

        for i in 0..3 {
            x[i] += dx[i];
        }
    }
    None
}

fn solve_linear_system(a: [[f64; 3]; 3], b: [f64; 3]) -> Option<[f64; 3]> {
    let mut a = a;
    let mut b = b;
    for i in 0..3 {
        // Pivot
        let mut max_row = i;
        for k in (i + 1)..3 {
            if a[k][i].abs() > a[max_row][i].abs() {
                max_row = k;
            }
        }
        a.swap(i, max_row);
        b.swap(i, max_row);

        let diag = a[i][i];
        if diag.abs() < 1e-12 {
            return None;
        }

        for j in i..3 {
            a[i][j] /= diag;
        }
        b[i] /= diag;

        for k in 0..3 {
            if k != i {
                let factor = a[k][i];
                for j in i..3 {
                    a[k][j] -= factor * a[i][j];
                }
                b[k] -= factor * b[i];
            }
        }
    }

    // Round the result to avoid floating-point precision issues
    Some(b.map(|x| (x * 1e9).round() / 1e9))
}

fn numerical_lambda_hybrid_v2(
    n: f64,
    logq: f64,
    sigma_e: f64,
    h: f64,
    initial_guess: Option<[f64; 3]>,
) -> f64 {
    let lnq = logq * f64::ln(2.0);
    let sigma_s = (h / n).sqrt();
    let xi = sigma_e / sigma_s;
    let mut rt_min = f64::INFINITY;

    for wg in 2..52 {
        let eq1 = |ng_: f64, beta_: f64, d_: f64| {
            approx_binom(ng_, wg as f64) + wg as f64 + d_.log2() - (0.292 * beta_ + 16.4 + 3.0)
                + 2.0
        };

        let eq2 = |ng_: f64, beta_: f64, d_: f64| {
            let delta = _delta(beta_);
            let val = (n * lnq / delta.ln()).sqrt().ceil();
            d_ - val + ng_ - 1.0
        };

        let eq3 = |ng_: f64, beta_: f64, d_: f64| {
            let delta = _delta(beta_);
            (-d_ + 1.0) * delta.log2() + ((d_ - n + ng_ - 1.0) * logq + (n - ng_) * xi.log2()) / d_
                - (2.0 * sigma_e * sigma_e).log2()
        };

        let system = |x: [f64; 3]| {
            let (ng_, beta_, d_) = (x[0], x[1], x[2]);
            [
                eq1(ng_, beta_, d_),
                eq2(ng_, beta_, d_),
                eq3(ng_, beta_, d_),
            ]
        };

        // let initial = initial_guess.unwrap_or_else(|| {
        //     let bdd = approx_startpoint_bdd(n as i32, logq, h);
        //     let beta_start = bdd[0];
        //     let d_start = (2.0 * n * lnq * beta_start / (beta_start / CONST).ln()).sqrt();
        //     [n / 4.0, beta_start, d_start - n / 4.0]
        // });

        let result = newton_solve(initial_guess.unwrap(), system, 50, 1e-4);

        if let Some(res) = result {
            let sol_tol = {
                let [f1, f2, f3] = system(res);
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
        let eq1a = |ng: f64, beta: f64, d: f64| {
            l as f64 - (0.292 * beta + 16.4 + 3.0 + log2(d))
                + approx_binom(n as f64 - h, ng - wg as f64)
                + approx_binom(h as f64, wg as f64)
                - approx_binom(n as f64, ng)
                - 1.0
        };

        let eq1b = |ng: f64, beta: f64, d: f64| {
            l as f64 - approx_binom(ng as f64, wg as f64) - wg as f64 - 2.0 * log2(d)
                + approx_binom(n as f64 - h, ng - wg as f64)
                + approx_binom(h as f64, wg as f64)
                - approx_binom(n as f64, ng)
                - 12.0
        };

        let eq2 = |ng: f64, beta: f64, d: f64, logq: f64| {
            d - (n as f64 * logq / log_delta(beta)).sqrt().ceil() + ng - 1.0
        };

        let eq3 = |ng: f64, beta: f64, d: f64, logq: f64| {
            (-d + 1.0) * log_delta(beta)
                + ((d - n as f64 + ng - 1.0) * logq + (n as f64 - ng) * log2(xi)) / d
                - log2(2.0 * sigma_e * sigma_e)
        };

        let system = |x: &[f64; 4]| -> [f64; 4] {
            [
                eq1a(x[0], x[1], x[2]),
                eq1b(x[0], x[1], x[2]),
                eq2(x[0], x[1], x[2], x[3]),
                eq3(x[0], x[1], x[2], x[3]),
            ]
        };

        // Use a numerical solver here (e.g., newton_raphson or argmin crate)
        // Placeholder: Just use initial guess as dummy solution
        let res = [
            n as f64 / 16.0,
            beta_initial_guess,
            d_initial_guess,
            logq_initial_guess,
        ];

        let rt =
            0.292 * res[1] + log2(8.0 * res[2]) + 16.4 - probability_enum(n, h, res[0], wg as f64);

        let eq1a_tolerance = (eq1a(res[0].round(), res[1].round(), res[2].round())).abs();
        let sol_tolerance = (eq1a(res[0], res[1], res[2])).abs()
            + (eq1b(res[0], res[1], res[2])).abs()
            + (eq2(res[0], res[1], res[2], res[3])).abs()
            + (eq3(res[0], res[1], res[2], res[3])).abs();

        if sol_tolerance < 2.5 && eq1a_tolerance < 0.7 && (l as f64 - rt).abs() < 2.5 {
            let ng_min = res[0].round();
            let beta_min = res[1].round();
            let d_min = res[2].round();
            let logq_min = res[3].round();

            sol_qs.push(logq_min);
            sols.push(vec![wg as f64, ng_min, beta_min, d_min, logq_min]);
        }
    }

    (sol_qs, sols)
}

// fn approx_startpoint_bdd(n: i32, logq: f64, h: f64) -> [f64; 2] {
//     println!(
//         "Starting approx_startpoint_bdd with n = {}, logq = {}, h = {}",
//         n, logq, h
//     );

//     let sigma_s = (h / n as f64).sqrt();
//     println!("Calculated sigma_s = {}", sigma_s);

//     let lnq = logq * LN_2;
//     println!("Calculated lnq = {}", lnq);

//     let zeta = (SIGMA_E / sigma_s).round();
//     println!("Calculated zeta = {}", zeta);

//     let mut beta_initial_guess = n as f64 / 4.0;
//     println!("Initial beta_initial_guess = {}", beta_initial_guess);

//     let nom = |beta: f64| {
//         let result = 2.0 * n as f64 * lnq * (beta / CONST).ln();
//         println!("nom(beta = {}): {}", beta, result);
//         result
//     };

//     let denom = |beta: f64| {
//         let result = (beta / CONST).ln() + 2.0 * lnq
//             - 2.0 * SIGMA_E.ln()
//             - CONST.ln()
//             - 2.0
//                 * (lnq - zeta.ln())
//                 * ((n as f64 * (beta / CONST).ln()) / (2.0 * lnq * beta)).sqrt();
//         println!("denom(beta = {}): {}", beta, result);
//         result
//     };

//     let eq6 = |beta: f64| {
//         let result = beta - nom(beta) / denom(beta).powi(2);
//         println!("eq6(beta = {}): {}", beta, result);
//         result
//     };

//     // Root finding loop
//     let mut beta_solution = f64::INFINITY;
//     let mut trials = 100;

//     while trials > 0 {
//         trials -= 1;

//         println!(
//             "Trial {}: beta_initial_guess = {}",
//             100 - trials,
//             beta_initial_guess
//         );

//         // Create a convergence strategy (you must create a new one each time)
//         let mut convergency = SimpleConvergency {
//             eps: 1e-7,
//             max_iter: 1000,
//         };

//         // Use Brent’s method to find a root in a given interval
//         let result = find_root_brent(
//             beta_initial_guess,
//             beta_initial_guess + 20.0, // adjust the interval if needed
//             &eq6,
//             &mut convergency,
//         );

//         if let Ok(root) = result {
//             println!("Found root: {}", root);
//             beta_solution = root;
//             break;
//         } else {
//             println!("Failed to find root, increasing beta_initial_guess");
//             beta_initial_guess += 1.0; // Try a new initial interval
//         }
//     }

//     if beta_solution.is_infinite() {
//         println!("Failed to find a root for beta after 100 trials.");
//         return [f64::INFINITY, f64::INFINITY];
//     }

//     let d_optimal = (2.0 * n as f64 * lnq * beta_solution / (beta_solution / CONST).ln()).sqrt();
//     println!("Calculated d_optimal = {}", d_optimal);

//     [beta_solution, d_optimal]
// }

#[pymodule]
fn hybrid(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(numerical_logq_hybrid, m)?)?;

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
    fn test_newton_solve() {
        let system = |x: [f64; 3]| [x[0] - 1.0, x[1] - 2.0, x[2] - 3.0];
        let result = newton_solve([0.0, 0.0, 0.0], system, 100, 1e-6);
        assert_eq!(result, Some([1.0, 2.0, 3.0]));
    }

    #[test]
    fn test_solve_linear_system() {
        let a = [[2.0, 1.0, -1.0], [-3.0, -1.0, 2.0], [-2.0, 1.0, 2.0]];
        let b = [8.0, -11.0, -3.0];
        let result = solve_linear_system(a, b);
        assert_eq!(result, Some([2.0, 3.0, -1.0]));
    }

    #[test]
    fn test_numerical_lambda_hybrid_v2() {
        let result = numerical_lambda_hybrid_v2(1024.0, 1.0, 3.19, 64.0, None);
        assert!(result.is_finite());
    }

    #[test]
    fn test_numerical_logq_hybrid_runoptimize() {
        let (sol_qs, sols) = numerical_logq_hybrid_runoptimize(1024, 80, 3.19, 64.0, 700);
        assert!(!sol_qs.is_empty());
        assert!(!sols.is_empty());
    }
}
