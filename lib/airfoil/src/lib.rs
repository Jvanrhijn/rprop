use std::vec::Vec;
use math;

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

    pub fn area(&self) -> Result<f64, math::MathError> {
        math::trapz(&self.x_values, &self.y_values)
    }
}

#[cfg(test)]
mod tests {
    use super::Airfoil;
    const EPS: f64 = 1e-16;
    #[test]
    fn area() {
        let base = 1.0;
        let height = 1.0;
        let airfoil = Airfoil::new(vec![base, 0.0, 0.0, base], vec![height, 0.00, 0.00, -height]);
        let area = base*height;
        assert!((airfoil.area().unwrap() - area) < EPS);
    }
}
