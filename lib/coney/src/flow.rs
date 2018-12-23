use std::f64::consts::PI;
use math;

fn coney_u(y0: f64, y: f64, z: usize) -> f64 {
    let exp_factor = ((1.0 + y*y).sqrt() - (1.0 + y0*y0).sqrt()).exp();
    let frac = y0*((1.0 + y*y).sqrt() - 1.0)/(y*((1.0 + y0*y0).sqrt() - 1.0));
    (frac*exp_factor).powi(z as i32)
}

fn coney_f1(y0: f64, y: f64, z: usize) -> f64 {
    let prefactor = -1.0/(2.0*(z as f64)*y0)*((1.0 + y0.powi(2))/(1.0 + y.powi(2))).powf(0.25);
    let u = coney_u(y0, y, z);
    let ufrac = 1.0/(1.0/u - 1.0);
    let first_frac = (9.0*y0*y0 + 2.0)/(1.0 + y0*y0).powf(1.5);
    let second_frac = (3.0*y*y - 2.0)/(1.0 + y*y).powf(1.5);
    prefactor*(ufrac + 1.0/(24.0*(z as f64))*(first_frac + second_frac)*(1.0 + ufrac).ln())
}

fn coney_f2(y0: f64, y: f64, z: usize) -> f64 {
    let prefactor = 1.0/(2.0*(z as f64)*y0)*((1.0 + y0.powi(2))/(1.0 + y.powi(2))).powf(0.25);
    let u = coney_u(y0, y, z);
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

/// axial component of velocity field induced at r = rc by z helical vortex lines
/// originating at r = rv at pitch angle vp
fn axial_velocity_vortex_line(rc: f64, rv: f64, vp: f64, z: usize) -> f64 {
    let tanvp = vp.tan();
    let y = rc/(rv*tanvp);
    let y0 = 1.0/tanvp;
    if rc <= rv {
        // this is always 31 for some reason: z as f64 / (2.0*rv*tanvp)/(2.0*PI));
        let result = (z as f64 / (2.0*rv*tanvp) - y*z.pow(2) as f64*y0*coney_f1(y0, y, z)/rc)/(2.0*PI);
        result
    } else {
        // this is always 0 for some reason
        let result = (-y*z.pow(2) as f64*y0/rc*coney_f2(y0, y, z))/(2.0*PI);
        result
    }
}

/// tangential component of velocity field induced at r = rc by z helical vortex lines
/// originating at r = rv at pitch angle vp
fn tangential_velocity_vortex_line(rc: f64, rv: f64, vp: f64, z: usize) -> f64 {
    let tanvp = vp.tan();
    let y = rc/(rv*tanvp);
    let y0 = 1.0/tanvp;
    if rc <= rv {
        let result = (z.pow(2) as f64/rc*y0*coney_f1(y0, y, z))/(2.0*PI);
        result
    } else {
        let result = (z as f64/(2.0*rc) + z.pow(2) as f64*y0*coney_f2(y0, y, z)/rc)/(2.0*PI);
        result
    }
}

/// Total axial velocity induced by horseshoe vortex at
/// rv[vi], rvi[vi+1] on the control point at rc[ci]
/// including hub image effects due to hub of radius rh, for
/// a propeller of z blades
pub fn axial_velocity(ci: usize, vi: usize, rc: &[f64], rv: &[f64], vpitch: &[f64], rh: f64, z: usize) -> f64 {
    // vortex line contributions
    //println!("{}", (vpitch[vi].tan()*rc[vi]));
    let uaw1 = axial_velocity_vortex_line(rc[ci], rv[vi+1], (vpitch[vi].tan()*rc[vi]/rv[vi+1]).atan(), z);
    let uaw2 = axial_velocity_vortex_line(rc[ci], rv[vi], (vpitch[vi].tan()*rc[vi]/rv[vi]).atan(), z);
    // image hub contributions
    let uawh1 = axial_velocity_vortex_line(rc[ci], rh*rh/rv[vi+1], (vpitch[vi].tan()*rc[vi]/(rh*rh/rv[vi+1])).atan(), z);
    let uawh2 = axial_velocity_vortex_line(rc[ci], rh*rh/rv[vi], (vpitch[vi].tan()*rc[vi]/(rh*rh/rv[vi])).atan(), z);
    // horsehoe vortex + image hub vortex
    (uaw1 - uawh2) - (uaw2 - uawh2)
}

/// Total axial velocity induced by horseshoe vortex at
/// rv[vi], rvi[vi+1] on the control point at rc[ci]
/// including hub image effects due to hub of radius rh, for
/// a propeller of z blades
pub fn tangential_velocity(ci: usize, vi: usize, rc: &[f64], rv: &[f64], vpitch: &[f64], rh: f64, z: usize) -> f64 {
    // vortex line contributions
    let utw1 = tangential_velocity_vortex_line(rc[ci], rv[vi+1], (vpitch[vi].tan()*rc[vi]/rv[vi+1]).atan(), z);
    let utw2 = tangential_velocity_vortex_line(rc[ci], rv[vi], (vpitch[vi].tan()*rc[vi]/rv[vi]).atan(), z);
    // image hub contributions
    let utwh1 = tangential_velocity_vortex_line(rc[ci], rh*rh/rv[vi+1], (vpitch[vi].tan()*rc[vi]/(rh*rh/rv[vi+1])).atan(), z);
    let utwh2 = tangential_velocity_vortex_line(rc[ci], rh*rh/rv[vi], (vpitch[vi].tan()*rc[vi]/(rh*rh/rv[vi])).atan(), z);
    // horsehoe vortex + image hub vortex
    (utw1 - utwh2) - (utw2 - utwh2)
}
