#[macro_use]
extern crate itertools;
#[macro_use]
extern crate ndarray;
use ndarray::{
    Array1,
    Array2,
};
use ndarray_linalg::{
    Solve
};
use std::vec::Vec;
use propeller::Propeller;
mod flow;

type Result<T> = std::result::Result<T, ndarray_linalg::error::LinalgError>;

/// coney solver trait
/// defines methods that must be implemented for all coney solver structs
/// coney solver should take ownership of propeller, compute and
/// set its optimal hydrodynamic data,
/// and hand ownership back to its owner upon completion.
/// the ConeySolver can be consumed upon completion of the process
pub trait ConeySolver {
    /// Start the optimization process, consume the solver
    fn optimize_propulsor(self, threshold: f64) -> Result<Vec<Propeller>>;
}

pub struct ConeySolverSingle {
    prop: Propeller,
    matrix: Array2<f64>,
    vector: Array1<f64>,
    lagrange_mult: f64,
    vortex_pitch: Vec<f64>
}

impl ConeySolverSingle {

    pub fn new(prop: Propeller) -> Self {
        let n = *prop.specs().num_panels();
        let matrix = Array2::<f64>::zeros((n+1, n+1));
        let vector = Array1::<f64>::zeros(n+1);
        Self{prop, matrix, vector, lagrange_mult: -1.0, vortex_pitch: vec![0.0; n+1]}
    }

    fn fill_matrix(&mut self) {
        // scalar params
        let z = *self.prop.geometry().num_blades();
        let w = *self.prop.specs().rot_speed();
        let n = *self.prop.specs().num_panels();
        let rh = *self.prop.geometry().hub_radius();

        // per-panel data
        let rc = self.prop.control_points();
        let rv = self.prop.vortex_points();
        let dr = self.prop.radial_increment();
        let vt = self.prop.hydro_data().tangential_inflow();
        let ut = self.prop.hydro_data().tangential_vel_ind();

        for i in 0..n {
            self.matrix[[i, n]] = (vt[i] + w*rc[i])*dr[i];
            for j in 0..n {
                let uaij = flow::axial_velocity(i, j, rc, rv, &self.vortex_pitch, rh, z);
                let uaji = flow::axial_velocity(j, i, rc, rv, &self.vortex_pitch, rh, z);
                let utij = flow::tangential_velocity(i, j, rc, rv, &self.vortex_pitch, rh, z);
                let utji = flow::tangential_velocity(j, i, rc, rv, &self.vortex_pitch, rh, z);
                self.matrix[[i, j]] = uaij*rc[i]*dr[i] + uaji*rc[j]*dr[j]
                          + self.lagrange_mult*(utij*dr[i] + utji*dr[j]);
            }
            self.matrix[[n, i]] =  z as f64*(vt[i] + w*rc[i] + ut[i])*dr[i];
        }

    }

    fn fill_vector(&mut self) {
        let cd = self.prop.hydro_data().drag_coeffs();
        let pitch = self.prop.hydro_data().hydro_pitch();
        let rc = self.prop.control_points();
        let dr = self.prop.radial_increment();
        let z = *self.prop.geometry().num_blades();
        let thrust = *self.prop.specs().thrust();
        let w = *self.prop.specs().rot_speed();

        // velocities
        let vt = self.prop.hydro_data().tangential_inflow();
        let va = self.prop.hydro_data().axial_inflow();
        let ut = self.prop.hydro_data().tangential_vel_ind();
        let ua = self.prop.hydro_data().axial_vel_ind();

        // TODO: initialize chords properly
        let c = self.prop.geometry().chords();

        let v = izip!(va, vt, ua, ut, rc)
            .map(|(vai, vti, uai, uti, rci)| ((vai + uai).powi(2) + (vti + uti + rci*w).powi(2)).sqrt())
            .collect::<Vec<_>>();

        let drag: f64 = izip!(cd, v, pitch, dr, c)
            .map(|(cdi, vi, pi, dri, ci)| 0.5*z as f64*cdi*vi*vi*ci*(pi).sin()*dri).sum();

        let mut vector = izip!(va, rc, dr).map(|(vai, rci, dri)| -vai*rci*dri).collect::<Vec<_>>();
        vector.push(thrust + drag);
        self.vector = Array1::<f64>::from_vec(vector);
    }

    // Perform wake alignment, update propeller data (velocities, pitch)
    fn align_wake(&mut self, threshold: f64) -> Result<Vec<f64>> {
        // calculate initial pitch values
        self.update_pitch();


        let circulation = loop {
            let solution = self.optimize_circulation(threshold)?;
            let prev_pitch = self.prop.hydro_data().hydro_pitch().clone();
            self.update_pitch();
            let wake_res = prev_pitch.iter().zip(self.prop.hydro_data().hydro_pitch().iter())
                .map(|(prev, pitch)| (prev - pitch).abs()).collect::<Vec<_>>();

            // check convergence
            if wake_res.iter().filter(|&&x| x > threshold).collect::<Vec<_>>().len() == 0 {
                break solution;
            }
        };

        Ok(circulation.slice(s![..circulation.len()-1]).to_vec())
    }

    // Perform circulation optimization
    // return result of optimization
    fn optimize_circulation(&mut self, threshold: f64) -> Result<Array1<f64>> {
        let n = *self.prop.specs().num_panels();
        let mut sol_prev= Array1::<f64>::ones(n+1);

        loop {
            self.fill_matrix();
            self.fill_vector();

            let solution = self.matrix.solve(&self.vector)?;

            let sol_res = (&solution - &sol_prev).to_vec().into_iter().filter(|res| res.abs() > threshold)
                .collect::<Vec<_>>();

            if sol_res.len() == 0 {
                break Ok(solution)
            }
            sol_prev = solution.clone();

            self.update_velocities(&solution);

            // TODO: log stuff
            self.lagrange_mult = solution[n];
            //println!("{}", self.lagrange_mult);
        }

    }

    fn update_velocities(&mut self, solution: &Array1<f64>) {

        let rc = self.prop.control_points();
        let z = *self.prop.geometry().num_blades();
        let rv = self.prop.vortex_points();
        let rh = *self.prop.geometry().hub_radius();
        let n = *self.prop.specs().num_panels();

        let mut axial_velocity = vec![0.0; n];
        let mut tangential_velocity = vec![0.0; n];

        for i in 0..n {
            for j in 0..n {
                axial_velocity[i] += solution[j]*flow::axial_velocity(i, j, rc, rv, &self.vortex_pitch, rh, z);
                tangential_velocity[i] += solution[j]*flow::tangential_velocity(i, j, rc, rv, &self.vortex_pitch, rh, z);
            }
        }

        *self.prop.hydro_data_mut().axial_vel_ind_mut() = axial_velocity;
        *self.prop.hydro_data_mut().tangential_vel_ind_mut() = tangential_velocity;
    }

    fn update_pitch(&mut self) {
        {
            let rc = self.prop.control_points();
            let w = *self.prop.specs().rot_speed();
            let va = self.prop.hydro_data().axial_inflow();
            let vt = self.prop.hydro_data().tangential_inflow();
            let ua = self.prop.hydro_data().axial_vel_ind(); // this is fucked somehow
            let ut = self.prop.hydro_data().tangential_vel_ind();
            {
                *self.prop.hydro_data_mut().hydro_pitch_mut() = izip!(rc, va, vt, ua, ut)
                    .map(|(rci, vai, vti, uai, uti)| ((vai + uai) / (vti + uti + w*rci)).atan()).collect::<Vec<_>>();
            }
        }
        let rc = self.prop.control_points();
        let n = *self.prop.specs().num_panels();
        let hpitch = self.prop.hydro_data().hydro_pitch();
        let rv = self.prop.vortex_points();

        // interpolate to find pitch at vortex points
        for i in 0..n-1 {
            let slope = (hpitch[i+1] - hpitch[i])/(rc[i+1] - rc[i]);
            self.vortex_pitch[i+1] = hpitch[i] + slope*(rv[i+1] - rc[i]);
        }
        // extrapolate first and last vortex points
        let slope = (hpitch[1] - hpitch[0])/(rc[1] - rc[0]);
        self.vortex_pitch[0] = hpitch[0] + slope*(rv[0] - rc[0]);

        let slope = (hpitch[n-1] - hpitch[n-2])/(rc[n-1] - rc[n-2]);
        self.vortex_pitch[n] = hpitch[n-1] + slope*(rv[n] - rc[n-1]);
    }

}

impl ConeySolver for ConeySolverSingle {
    fn optimize_propulsor(mut self, threshold: f64) -> Result<Vec<Propeller>> {
        // TODO: have align_wake return convergence error if not converged, LinalgError if ndarray failed
        *self.prop.hydro_data_mut().circulation_mut() = self.align_wake(threshold)?;
        Ok(vec![self.prop])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use propeller::PropellerBuilder;

    #[test]
    fn optimize_single_prop() {
        let prop = PropellerBuilder::new(0.1, 0.025, 150.0, 200.0, 10.0, 20)
            .build();
        let coney_solver = ConeySolverSingle::new(prop);
        let prop = coney_solver.optimize_propulsor(1e-6).unwrap().remove(0);
        for &g in prop.hydro_data().circulation() {
            assert!(g > 0.0);
        }

    }
}
