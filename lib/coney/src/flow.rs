use std::f64::consts::PI;
use propeller::Propeller;
use math;

fn coney_u(y0: f64, y: f64, prop: &Propeller) -> f64 {
    let z = *prop.geometry().num_blades();
    let exp_factor = ((1.0 + y*y).sqrt() - (1.0 + y0*y0).sqrt()).exp();
    let frac = y0*((1.0 + y*y).sqrt() - 1.0)/(y*((1.0 + y0*y0).sqrt() - 1.0));
    (frac*exp_factor).powi(z as i32)
}

fn coney_f1(y0: f64, y: f64, prop: &Propeller) -> f64 {
    let z = *prop.geometry().num_blades();
    let prefactor = -1.0/(2.0*(z as f64)*y0)*((1.0 + y0.powi(2))/(1.0 + y.powi(2))).powf(0.25);
    let u = coney_u(y0, y, prop);
    let ufrac = 1.0/(1.0/u - 1.0);
    let first_frac = (9.0*y0*y0 + 2.0)/(1.0 + y0*y0).powf(1.5);
    let second_frac = (3.0*y*y - 2.0)/(1.0 + y*y).powf(1.5);
    prefactor*(ufrac + 1.0/(24.0*(z as f64))*(first_frac + second_frac)*(1.0 + ufrac).ln())
}

fn coney_f2(y0: f64, y: f64, prop: &Propeller) -> f64 {
    let z = *prop.geometry().num_blades();
    let prefactor = 1.0/(2.0*(z as f64)*y0)*((1.0 + y0.powi(2))/(1.0 + y.powi(2))).powf(0.25);
    let u = coney_u(y0, y, prop);
    let ufrac = 1.0/(u - 1.0);
    let first_frac = (9.0*y0*y0 + 2.0)/(1.0 + y0*y0).powf(1.5);
    let second_frac = (3.0*y*y - 2.0)/(1.0 + y*y).powf(1.5);
    prefactor*(ufrac - 1.0/(24.0*(z as f64))*(first_frac + second_frac)*(1.0 + ufrac).ln())
}

fn c1_jl(xf: f64, rcjm: f64, rvlp: f64) -> f64 {
    let q = 1.0 + (xf*xf + (rcjm - rvlp).powi(2))/(2.0*rcjm*rvlp);
    let s = (xf/(xf*xf + (rcjm - rvlp).powi(2)).sqrt()).asin();
    let t = (4.0*rcjm*rvlp/(xf*xf + (rcjm + rvlp).powi(2))).sqrt();
    let frac = 0.5*xf/(rcjm*rvlp).sqrt() * math::legendre_2mhalf(q);
    let heuman = 0.5*PI*math::heuman(s, t.asin());
    if rcjm > rvlp {
        frac - heuman
    } else {
        PI + frac + heuman
    }
}

fn c2_jl(xf: f64, rcjm: f64, rvlp: f64) -> f64 {
    let q = 1.0 + (xf*xf + (rcjm - rvlp).powi(2))/(2.0*rcjm*rvlp);
    let s = (xf/(xf*xf + (rcjm - rvlp).powi(2)).sqrt()).asin();
    let t = (4.0*rcjm*rvlp/(xf*xf + (rcjm + rvlp).powi(2))).sqrt();
    let frac = 0.5*xf/(rcjm*rvlp).sqrt() * math::legendre_2mhalf(q);
    let heuman = 0.5*PI*math::heuman(s, t.asin());
    if rcjm <= rvlp {
        frac - heuman
    } else {
        PI + frac + heuman
    }
}

struct FlowField {

}