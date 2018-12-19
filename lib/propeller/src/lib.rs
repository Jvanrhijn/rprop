use std::vec::Vec;
use airfoil::{Airfoil, naca2414};

#[derive(Clone)]
pub struct Propeller {
    geometry: Geometry,
    specs: DesignSpecs,
    hydro_data: HydrodynamicData,
    control_points: Vec<f64>,
    vortex_points: Vec<f64>,
    radial_increment: Vec<f64>,
    dimensional: bool
}

#[derive(Clone)]
struct Geometry {
    radius: f64,
    hub_radius: f64,
    chords: Vec<f64>,
    base_airfoil: Airfoil
}

#[derive(Clone)]
struct DesignSpecs {
    rot_speed: f64,
    ship_speed: f64,
    thrust: f64,
    num_panels: usize,
}

#[derive(Clone)]
struct HydrodynamicData {
    axial_inflow: Vec<f64>,
    tangential_inflow: Vec<f64>,
    axial_vel_ind: Vec<f64>,
    tangential_vel_ind: Vec<f64>,
    hydro_pitch: Vec<f64>,
    drag_coeffs: Vec<f64>,
    circulation: Vec<f64>
}

pub struct PropellerBuilder {
    // required data for all propellers:
    radius: f64,
    hub_radius: f64,
    thrust: f64,
    rot_speed: f64,
    ship_speed: f64,
    num_panels: usize,
    // optional data with default values
    // geometric:
    chords: Option<Vec<f64>>,
    airfoil: Option<Airfoil>,
    // hydrodynamic:
    axial_inflow: Option<Vec<f64>>, // default filled with ship_speed
    tangential_inflow: Option<Vec<f64>>, // default zero
    drag_coeffs: Option<Vec<f64>>, // default zero
    // configurational data
    dim: bool // default false

}

impl PropellerBuilder {
    // required data for building propeller
    pub fn new(radius: f64, hub_radius: f64, thrust: f64, rot_speed: f64, ship_speed: f64, num_panels: usize) -> Self {
        Self {
            radius,
            hub_radius,
            thrust,
            rot_speed,
            ship_speed,
            num_panels,
            chords: None,
            airfoil: None,
            axial_inflow: None,
            tangential_inflow: None,
            drag_coeffs: None,
            dim: false
        }
    }

    pub fn dimensional(mut self, dim: bool) -> Self {
        self.dim = dim;
        self
    }

    pub fn axial_inflow(mut self, axial_inflow: Vec<f64>) -> Self {
        self.axial_inflow = Some(axial_inflow);
        self
    }

    pub fn tangential_inflow(mut self, tangential_inflow: Vec<f64>) -> Self {
        self.tangential_inflow = Some(tangential_inflow);
        self
    }

    pub fn build(self) -> Propeller {
        // TODO: apply nondimensionalization
        let geometry = Geometry{
            radius: self.radius,
            hub_radius: self.hub_radius,
            chords: self.chords.unwrap_or(vec![0.0; self.num_panels]),
            base_airfoil: self.airfoil.unwrap_or(naca2414()),
        };
        let hydro_data = HydrodynamicData{
            axial_inflow: self.axial_inflow.unwrap_or(vec![self.ship_speed; self.num_panels]),
            tangential_inflow: self.tangential_inflow.unwrap_or(vec![0.0; self.num_panels]),
            // TODO: provide good default values for these
            axial_vel_ind: vec![0.0; self.num_panels],
            tangential_vel_ind: vec![0.0; self.num_panels],
            hydro_pitch: vec![0.0; self.num_panels],
            drag_coeffs: self.drag_coeffs.unwrap_or(vec![0.0; self.num_panels]),
            circulation: vec![0.0; self.num_panels],
        };
        let specs = DesignSpecs{
            rot_speed: self.rot_speed,
            ship_speed: self.ship_speed,
            thrust: self.thrust,
            num_panels: self.num_panels
        };
        // generate vortex/control points
        let dr = (self.radius - self.hub_radius)/(self.num_panels as f64 + 0.25);
        let mut vortex_points = vec![self.hub_radius];
        let mut control_points = vec![self.hub_radius + 0.5*dr];
        for i in 1..self.num_panels {
            vortex_points.push(vortex_points[i-1] + dr);
            control_points.push(control_points[i-1] + dr);

            vortex_points.push(self.radius - 0.25*dr); // tip inset
        }
        let radial_increment = vortex_points.iter().zip(vortex_points[1..].iter())
            .map(|(rv, rvp1)| rvp1 - rv).collect::<Vec<_>>();
        Propeller{
            geometry,
            specs,
            hydro_data,
            control_points,
            vortex_points,
            radial_increment,
            dimensional: self.dim,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn build_propeller_simple() {
        let propeller = PropellerBuilder::new(1.0, 0.25, 1000.0, 200.0, 10.0, 20)
            .tangential_inflow(vec![0.1; 20])
            .axial_inflow(vec![11.0; 20])
            .build();
    }
}
