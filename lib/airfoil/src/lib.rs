use std::vec::Vec;
use std::io::prelude::*;
use std::io::{self, BufReader, BufRead};
use std::fs::File;
use math;

#[derive(Clone)]
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
                .expect("Failed to parse f64")).collect::<Vec<_>>();
            xvals.push(coordinates[0]);
            yvals.push(coordinates[1]);
        }
        Ok(Self{x_values: xvals, y_values: yvals})
    }

    pub fn save_to_file(&self, fpath: &str) -> Result<(), io::Error> {
        let mut file = File::create(fpath)?;
        for (x, y) in self.x_values.iter().zip(self.y_values.iter()) {
            file.write(&format!("{} {}\n", x, y).as_bytes())?;
        }
        Ok(())
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

impl Default for Airfoil {
    fn default() -> Self {
        naca2414()
    }
}

fn naca2414() -> Airfoil {
    Airfoil::new(
        vec![
                1.00000,
                0.99739,
                0.98929,
                0.97587,
                0.95729,
                0.93372,
                0.90542,
                0.87267,
                0.83582,
                0.79527,
                0.75143,
                0.70480,
                0.65586,
                0.60515,
                0.55324,
                0.50069,
                0.44808,
                0.39598,
                0.34454,
                0.29482,
                0.24740,
                0.20285,
                0.16169,
                0.12440,
                0.09141,
                0.06310,
                0.03977,
                0.02165,
                0.00892,
                0.00169,
                0.00000,
                0.00379,
                0.01293,
                0.02730,
                0.04669,
                0.07087,
                0.09957,
                0.13246,
                0.16918,
                0.20937,
                0.25260,
                0.29844,
                0.34644,
                0.39611,
                0.44739,
                0.49931,
                0.55129,
                0.60276,
                0.65316,
                0.70194,
                0.74857,
                0.79252,
                0.83331,
                0.87048,
                0.90360,
                0.93230,
                0.95626,
                0.97518,
                0.98886,
                0.99713,
                1.00000
            ],
        vec![
                0.00147,
                0.00210,
                0.00396,
                0.00700,
                0.01112,
                0.01620,
                0.02207,
                0.02857,
                0.03552,
                0.04274,
                0.05004,
                0.05723,
                0.06412,
                0.07053,
                0.07629,
                0.08120,
                0.08512,
                0.08787,
                0.08913,
                0.08866,
                0.08645,
                0.08255,
                0.07707,
                0.07014,
                0.06198,
                0.05281,
                0.04289,
                0.03245,
                0.02171,
                0.01085,
                0.00000,
                -0.01031,
                -0.01956,
                -0.02770,
                -0.03471,
                -0.04054,
                -0.04516,
                -0.04858,
                -0.05082,
                -0.05195,
                -0.05208,
                -0.05133,
                -0.04987,
                -0.04787,
                -0.04537,
                -0.04232,
                -0.03886,
                -0.03516,
                -0.03132,
                -0.02745,
                -0.02365,
                -0.01998,
                -0.01650,
                -0.01328,
                -0.01035,
                -0.00776,
                -0.00557,
                -0.00381,
                -0.00252,
                -0.00173,
                -0.00147
        ])
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
        assert!((airfoil.area().unwrap() - area).abs() < EPS);
    }

    #[test]
    fn default_foil() {
        let foil = Airfoil::default();
        assert!((foil.area().unwrap() - 0.096).abs() < 0.001);
    }
}
