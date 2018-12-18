use std::vec::Vec;

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
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
