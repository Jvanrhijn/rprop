use std::error;
use std::fmt;

#[derive(Copy, Clone, Debug)]
pub struct MathError;

impl fmt::Display for MathError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Invalid math result")
    }
}

impl error::Error for MathError {
    fn description(&self) -> &str {
        "invalid first item to double"
    }

    fn cause(&self) -> Option<&error::Error> {
        None
    }
}

// Trapezoidal integration tule
pub fn trapz(x: &[f64], y: &[f64]) -> Result<f64, MathError> {
    match x.len() == y.len() {
        true => {
            let mut res = 0.0;
            let num_points = x.len();
            for i in 0..num_points - 2 {
                let dx = x[i + 1] - x[i];
                res += 0.5 * dx * (y[i + 1] + y[i]);
            }
            let dx = x[num_points - 1] - x[num_points - 2];
            res += 0.5 * dx * (y[num_points - 1] + y[num_points - 2]);
            Ok(res)
        },
        false => Err(MathError)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn trapz_test() {
        let x: Vec<f64> = (0..1001).map(|x| x as f64 * 0.001).collect();
        let y: Vec<f64> = x.iter().map(|x| x.powi(2)).collect();
        let integral_exact = 1.0/3.0;
        println!("{}", trapz(&x, &y).unwrap());
        assert!((trapz(&x, &y).unwrap() - integral_exact).abs() < 1e-4);
    }
}
