use pyo3::prelude::*;
use rand::Rng;
use std::cmp::Ordering;
use std::f64::consts::{E, PI};
use std::f64::{consts::LN_2, INFINITY};
const CONST: f64 = 1.01; // Placeholder for `const` from original Python code
const SIGMA_E: f64 = 3.19;

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

    let mut prob = 0.0;

    for i in 0..w {
        let numerator = binomial((n - h) as u64, (ng - i) as u64) * binomial(h as u64, i as u64);
        let denominator = binomial(n as u64, ng as u64);
        prob += numerator as f64 / denominator as f64;
    }

    if prob <= 0.0 {
        f64::NEG_INFINITY // avoid log2(0)
    } else {
        prob.log2()
    }
}

fn binomial(n: u64, k: u64) -> u64 {
    if k > n {
        return 0;
    }
    if k == 0 || k == n {
        return 1;
    }

    let k = k.min(n - k); // Use symmetry
    let mut result = 1u64;

    for i in 0..k {
        result = result * (n - i) / (i + 1);
    }

    result
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
        let l_ = numerical_lambda_hybrid_v2(n, logq, sigma_e, h, Some(initial_guesses[i]));

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

fn numerical_lambda_hybrid_v2(
    n: i32,
    logq: f64,
    sigma_e: f64,
    h: f64,
    initial_guess: Option<[f64; 3]>,
) -> f64 {
    let lnq = logq * LN_2;
    let sigma_s = (h / n as f64).sqrt();
    let xi = sigma_e / sigma_s;
    let mut rt_min = f64::INFINITY;

    for wg in 2..52 {
        let eq1 = |ng: f64, beta: f64, d: f64| {
            approx_binom(ng, wg as f64) + wg as f64 + log2(d) - (0.292 * beta + 16.4 + 3.0) + 2.0
        };

        let eq2 = |ng: f64, beta: f64, d: f64| {
            d - (n as f64 * lnq / log2(_delta(beta))).sqrt().ceil() + ng - 1.0 //TODO: CHANGE log2 to log
        };

        let eq3 = |ng: f64, beta: f64, d: f64| {
            (-d + 1.0) * log2(_delta(beta))
                + ((d - n as f64 + ng - 1.0) * logq + (n as f64 - ng) * log2(xi)) / d
                - log2(2.0 * sigma_e * sigma_e)
        };

        let mut guess = if let Some(g) = initial_guess {
            g
        } else {
            let initial_bdd = approx_startpoint_bdd(n, logq, h);
            let beta_start = initial_bdd[0];
            let d_start = (2.0 * n as f64 * lnq * beta_start / log2(beta_start / CONST)).sqrt(); //TODO: CHANGE log2 to log
            [n as f64 / 4.0, beta_start, d_start - n as f64 / 4.0]
        };

        let mut bound_trials = 250;
        let mut sol_tolerance = 100.0;
        let mut res = [0.0; 3];

        while bound_trials > 0 {
            bound_trials -= 1;

            // === Replace this with a real solver ===
            // res = dummy_fsolve(
            //     &|x: [f64; 3]| {
            //         [
            //             eq1(x[0], x[1], x[2]),
            //             eq2(x[0], x[1], x[2]),
            //             eq3(x[0], x[1], x[2]),
            //         ]
            //     },
            //     guess,
            // );

            // Placeholder: Just use initial guess as dummy solution
            res = guess;
            guess = [res[0] + 0.1, res[1] + 0.1, res[2] + 0.1];

            sol_tolerance = eq1(res[0], res[1], res[2]).abs()
                + eq2(res[0], res[1], res[2]).abs()
                + eq3(res[0], res[1], res[2]).abs();

            if sol_tolerance < 7.0 {
                break;
            }

            let mut rng = rand::thread_rng();
            let beta_start = (guess[1] + rng.gen_range(-1500.0..30.0)).max(40.0);
            let d_start = (2.0 * n as f64 * lnq * beta_start / log2(beta_start / CONST)).sqrt(); //TODO: CHANGE log2 to log
            guess = [n as f64 / 4.0, beta_start, d_start - n as f64 / 4.0];
        }

        let rt =
            0.292 * res[1] + log2(8.0 * res[2]) + 16.4 - probability_enum(n, h, res[0], wg as f64);
        if rt < rt_min && sol_tolerance < 7.0 {
            rt_min = rt;
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

fn approx_startpoint_bdd(n: i32, logq: f64, h: f64) -> [f64; 2] {
    let sigma_s = (h / n as f64).sqrt();
    let lnq = logq * LN_2;
    let zeta = (SIGMA_E / sigma_s).round();
    let beta_initial_guess = n as f64 / 4.0;

    // nominator of the equation
    let nom = |beta: f64| 2.0 * n as f64 * lnq * (beta / CONST).ln();

    // denominator of the equation
    let denom = |beta: f64| {
        (beta / CONST).ln() + 2.0 * lnq
            - 2.0 * SIGMA_E.ln()
            - CONST.ln()
            - 2.0
                * (lnq - (zeta.ln()))
                * ((n as f64 * (beta / CONST).ln() / (2.0 * lnq * beta)).sqrt())
    };

    // Full equation
    let eq6 = |beta: f64| beta - nom(beta) / denom(beta).powi(2);

    // Solve using Newton-Raphson or similar
    let beta_solution = 40 as f64; //TODO CHANGE THIS

    // Compute d
    let d_optimal = (2.0 * n as f64 * lnq * beta_solution / (beta_solution / CONST).ln()).sqrt();

    [beta_solution, d_optimal]
}

#[pymodule]
fn logq_hybrid(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(numerical_logq_hybrid, m)?)?;

    Ok(())
}
