#[macro_use]
extern crate getset;
#[macro_use]
extern crate itertools;
// std
use std::vec::Vec;
// first party
use airfoil::Airfoil;

const WATER_DENSITY: f64 = 1000.0;


#[derive(Clone, Getters, MutGetters)]
pub struct Propeller {
    #[get = "pub"] #[get_mut = "pub"]
    geometry: Geometry,
    #[get = "pub"] #[get_mut = "pub"]
    specs: DesignSpecs,
    #[get = "pub"] #[get_mut = "pub"]
    hydro_data: HydrodynamicData,
    #[get = "pub"] #[get_mut = "pub"]
    control_points: Vec<f64>,
    #[get = "pub"] #[get_mut = "pub"]
    vortex_points: Vec<f64>,
    #[get = "pub"] #[get_mut = "pub"]
    radial_increment: Vec<f64>,
    dimensional: bool
}

#[derive(Clone, Getters, MutGetters)]
pub struct Geometry {
    #[get = "pub"] #[get_mut = "pub"]
    radius: f64,
    #[get = "pub"] #[get_mut = "pub"]
    hub_radius: f64,
    #[get = "pub"] #[get_mut = "pub"]
    num_blades: usize,
    #[get = "pub"] #[get_mut = "pub"]
    chords: Vec<f64>,
    base_airfoil: Airfoil
}

#[derive(Clone, Getters, MutGetters)]
pub struct DesignSpecs {
    #[get = "pub"] #[get_mut = "pub"]
    rot_speed: f64,
    #[get = "pub"] #[get_mut = "pub"]
    ship_speed: f64,
    #[get = "pub"] #[get_mut = "pub"]
    thrust: f64,
    #[get = "pub"] #[get_mut = "pub"]
    num_panels: usize,
}

#[derive(Clone, Getters, MutGetters)]
pub struct HydrodynamicData {
    #[get = "pub"] #[get_mut = "pub"]
    axial_inflow: Vec<f64>,
    #[get = "pub"] #[get_mut = "pub"]
    tangential_inflow: Vec<f64>,
    #[get = "pub"] #[get_mut = "pub"]
    axial_vel_ind: Vec<f64>,
    #[get = "pub"] #[get_mut = "pub"]
    tangential_vel_ind: Vec<f64>,
    #[get = "pub"] #[get_mut = "pub"]
    hydro_pitch: Vec<f64>,
    #[get = "pub"] #[get_mut = "pub"]
    drag_coeffs: Vec<f64>,
    #[get = "pub"] #[get_mut = "pub"]
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
    num_blades: usize, // default 2
    // hydrodynamic:
    axial_inflow: Option<Vec<f64>>, // default filled with ship_speed
    tangential_inflow: Option<Vec<f64>>, // default zero
    drag_coeffs: Option<Vec<f64>>, // default zero
    // configurational data
    dim: bool // default false

}

struct DimensionScale {
    pub length: f64,
    pub speed: f64,
    pub rotation_speed: f64,
    pub force: f64,
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
            chords: Some(vec![0.0; num_panels]),
            airfoil: Some(Airfoil::default()),
            num_blades: 2,
            axial_inflow: None,
            tangential_inflow: None,
            drag_coeffs: None,
            dim: false
        }
    }

    pub fn num_blades(mut self, num: usize) -> Self {
        self.num_blades = num;
        self
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

    pub fn chords(mut self, chords: Vec<f64>) -> Self {
        self.chords = Some(chords);
        self
    }

    pub fn build(mut self) -> Propeller {
        // TODO: more idiomatic solution for this
        let dim = DimensionScale{
            length: if self.dim {1.0} else {self.radius},
            speed: if self.dim {1.0} else {self.ship_speed},
            rotation_speed: if self.dim{1.0} else {self.ship_speed/self.radius},
            force: if self.dim {1.0} else { WATER_DENSITY*self.radius.powi(2)*self.ship_speed.powi(2)}
        };
        let geometry = Geometry{
            radius: self.radius/dim.length,
            hub_radius: self.hub_radius/dim.length,
            num_blades: self.num_blades,
            chords: self.chords.take().unwrap().into_iter().map(|c| c/dim.length).collect(),
            base_airfoil: self.airfoil.take().unwrap()
        };
        let specs = DesignSpecs{
            rot_speed: self.rot_speed/dim.rotation_speed,
            ship_speed: self.ship_speed/dim.speed,
            thrust: self.thrust/dim.force,
            num_panels: self.num_panels
        };
        // generate vortex/control points
        let dr = (self.radius - self.hub_radius)/(self.num_panels as f64 + 0.25)/dim.length;
        let mut vortex_points = vec![self.hub_radius/dim.length];
        let mut control_points = vec![self.hub_radius/dim.length + 0.5*dr];
        for i in 1..self.num_panels {
            vortex_points.push(vortex_points[i-1] + dr);
            control_points.push(control_points[i-1] + dr);
        }
        vortex_points.push(self.radius/dim.length - 0.25*dr); // tip inset
        let radial_increment = vortex_points.iter().zip(vortex_points[1..].iter())
            .map(|(rv, rvp1)| rvp1 - rv).collect::<Vec<_>>();
        let axial_inflow = self.axial_inflow.get_or_insert(vec![self.ship_speed; self.num_panels])
            .iter().map(|v| v/dim.speed).collect::<Vec<_>>();
        let tangential_inflow = self.tangential_inflow.get_or_insert(vec![0.0; self.num_panels])
            .iter().map(|v| v/dim.speed).collect::<Vec<_>>();
        let hydro_pitch = izip!(axial_inflow.iter(), tangential_inflow.iter(), control_points.iter())
            .map(|(va, vt, rc)| (va/(vt + self.rot_speed*rc)).atan()).collect::<Vec<_>>();
        let hydro_data = HydrodynamicData{
            axial_inflow,
            tangential_inflow,
            // TODO: provide good default values for these
            axial_vel_ind: vec![0.0; self.num_panels]
                .into_iter().map(|v| v/dim.speed).collect(),
            tangential_vel_ind: vec![0.0; self.num_panels]
                .into_iter().map(|v| v/dim.speed).collect(),
            hydro_pitch,
            drag_coeffs: self.drag_coeffs.unwrap_or(vec![0.0; self.num_panels]),
            circulation: vec![0.0; self.num_panels]
                .into_iter().map(|c| c/(dim.length*dim.speed)).collect(),
        };
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
            .dimensional(true)
            .build();
        propeller.hydro_data().axial_inflow().iter().for_each(|v| assert_eq!(*v, 11.0));
        propeller.hydro_data().tangential_inflow().iter().for_each(|v| assert_eq!(*v, 0.1));
    }

    #[test]
    fn build_propeller_nondimensional() {
        let propeller = PropellerBuilder::new(1.0, 0.25, 1000.0, 200.0, 10.0, 20)
            .tangential_inflow(vec![0.1; 20])
            .axial_inflow(vec![11.0; 20])
            .dimensional(false)
            .build();
        propeller.hydro_data().axial_inflow().iter().for_each(|v| assert_eq!(*v, 1.1));
        propeller.hydro_data().tangential_inflow().iter().for_each(|v| assert_eq!(*v, 0.01));
    }
}
