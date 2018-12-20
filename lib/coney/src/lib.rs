#[macro_use]
extern crate itertools;
use std::{error, fmt, vec::Vec};
extern crate rgsl;
use rgsl::{
    types::{
        matrix::MatrixF64,
        vector::VectorF64
    },
    linear_algebra
};
use propeller::Propeller;
mod flow;

#[derive(Copy, Clone, Debug)]
pub struct GslError;

impl fmt::Display for GslError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Error in GSL library")
    }
}

impl error::Error for GslError {
    fn description(&self) -> &str {
        "Error in GSL library"
    }

    fn cause(&self) -> Option<&error::Error> {
        None
    }
}

type Result<T> = std::result::Result<T, GslError>;

/// coney solver trait
/// defines methods that must be implemented for all coney solver structs
/// coney solver should take ownership of propeller, compute and
/// set its optimal hydrodynamic data,
/// and hand ownership back to its owner upon completion.
/// the ConeySolver can be consumed upon completion of the process
pub trait ConeySolver {
    /// Start the optimization process, consume the solver
    fn optimize_propulsor(self, threshold: f64) -> Vec<Propeller>;
}

pub struct ConeySolverSingle {
    prop: Propeller,
    matrix: MatrixF64,
    vector: VectorF64,
    lagrange_mult: f64,
}

impl ConeySolverSingle {
    pub fn new(prop: Propeller) -> Result<Self> {
        let n = *prop.specs().num_panels();
        let matrix = MatrixF64::new(n+1, n+1)
            .ok_or(GslError)?;
        let vector = VectorF64::new(n+1)
            .ok_or(GslError)?;
        Ok(Self{prop, matrix, vector, lagrange_mult: -1.0})
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
        let va = self.prop.hydro_data().axial_inflow();
        // TODO: make sure these are initialized properly
        let ut = self.prop.hydro_data().tangential_vel_ind();
        let pitch = self.prop.hydro_data().hydro_pitch();
        for i in 0..n {
            self.matrix.set(i, n, (vt[i] + w*rc[i])*dr[i]);
            for j in 0..n {
                let uaij = flow::axial_velocity(i, j, rc, rv, pitch, rh, z);
                let uaji = flow::axial_velocity(j, i, rc, rv, pitch, rh, z);
                let utij = flow::tangential_velocity(i, j, rc, rv, pitch, rh, z);
                let utji = flow::tangential_velocity(j, i, rc, rv, pitch, rh, z);
                self.matrix.set(i, j, uaij*rc[i]*dr[i] + uaji*rc[j]*dr[j]
                          + self.lagrange_mult + (utij*dr[i]       + utji*dr[j]));
            }
            self.matrix.set(n, i, z as f64*(vt[i] + w*rc[i] + ut[i])*dr[i]);
        }
    }

    fn fill_vector(&mut self) {

    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
