use std::error;
use std::fmt;
use std::num::ParseIntError;
use std::vec::Vec;

struct MathError;

impl fmt::Display for MathError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Invalid math result")
    }
}

impl error::Error for DoubleError {
    fn description(&self) -> &str {
        "invalid first item to double"
    }

    fn cause(&self) -> Option<&error::Error> {
        None
    }
}

pub struct Airfoil {
    x_values: Vec<f64>,
    y_values: Vec<f64>
}

impl Airfoil {
    pub fn new(xvals: Vec<f64>, yvals: Vec<f64>) -> Self {
        Self{x_values: xvals, y_values: yvals}
    }
    
    pub fn from_file(fpath: &str) -> Self {
        Self{x_values: vec![1.0], y_values: vec![1.0]}
    }

    pub fn area(&self) -> Result<f64, MathError> {
        trapz(self.x_values, self.y_values)
    }
}

// Trapezoidal integration tule
fn trapz(x: &[f64], y: &[f64]) -> Result<f64, MathError> {
    match x.size() == y.size() {
        true => {
            let mut res = 0.0;
            let num_points = x.len();
            for i in 0..num_points - 2 {
                let dx = x[i + 1] - x[i];
                result += 0.5 * dx * (y[i + 1] + y[i]);
            }
            let dx = x[num_points - 1] - x[num_points - 2];
            result += 0.5 * dx * y[num_points - 1] + y[num_points - 2];
            Ok(result)
        },
        false => DoubleError
    }
}

#[cfg(test)]
mod tests {
    use super::Airfoil;
    const EPS: f64 = 1e-16;
    #[test]
    fn area() {
        let airfoil = Airfoil::new(vec![1.0, 0.0, 0.0, 1.0], vec![0.0, 0.01, -0.01, 0.0]);
        assert!((airfoil.area() - 2*1.0*0.01) < EPS);
    }
}
