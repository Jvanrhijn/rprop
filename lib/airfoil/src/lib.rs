use std::vec::Vec;
use std::io;
use std::io::{BufReader, BufRead};
use std::fs::File;
use math;

pub struct Airfoil {
    x_values: Vec<f64>,
    y_values: Vec<f64>
}

impl Airfoil {
    pub fn new(xvals: Vec<f64>, yvals: Vec<f64>) -> Self {
        Self{x_values: xvals, y_values: yvals}
    }
    
    pub fn from_file(fpath: &str) -> Result<Self, io::Error> {
        let mut xvals = Vec::<f64>::new();
        let mut yvals = Vec::<f64>::new();
        let file = File::open(fpath)?;
        for line in BufReader::new(&file).lines() {
            let coordinates = line?.split_whitespace().map(|x| x.parse::<f64>()
                .expect("Failed to parse f64")).collect::<Vec<f64>>();
            xvals.push(coordinates[0]);
            yvals.push(coordinates[1]);
        }
        Ok(Self{x_values: xvals, y_values: yvals})
    }

    pub fn area(&self) -> Result<f64, math::MathError> {
        let (xtop, ytop) = self.top_surface();
        let (xbot, ybot) = self.bottom_surface();
        let ytop_min = ytop.iter().min_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
        let ybot_min = ybot.iter().min_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
        let dy = ytop_min.min(*ybot_min).abs();
        // shift airfoil up so it is entirely above x-axis
        let ytop = ytop.iter().map(|y| y + dy).collect::<Vec<_>>();
        let ybot = ybot.iter().map(|y| y + dy).collect::<Vec<_>>();
        // top area must be negated since x-axis is reversed
        let area_top = -math::trapz(&xtop, &ytop)?;
        let area_bot = math::trapz(&xbot, &ybot)?;
        Ok(area_top - area_bot)
    }

    fn top_surface(&self) -> (&[f64], &[f64]) {
        let num_points = self.x_values.len();
        (&self.x_values[..num_points/2], &self.y_values[..num_points/2])
    }

    fn bottom_surface(&self) -> (&[f64], &[f64]) {
        let num_points = self.x_values.len();
        (&self.x_values[num_points/2..], &self.y_values[num_points/2..])
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
