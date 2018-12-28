extern crate rs_xfoil;
use std::vec::Vec;
use propeller::Propeller;

pub struct GeometryOptimizer {
    prop: Propeller,
    inflow_vel: Vec<f64>
}

impl GeometryOptimizer {
    pub fn new(prop: Propeller) -> Self {
        let rc = prop.control_points();
        let w = prop.specs().rot_speed();
        // velocities
        let vt = prop.hydro_data().tangential_inflow();
        let va = prop.hydro_data().axial_inflow();
        let ut = prop.hydro_data().tangential_vel_ind();
        let ua = prop.hydro_data().axial_vel_ind();

        // TODO: initialize chords properly
        let c = prop.geometry().chords();

        let v = izip!(va, vt, ua, ut, rc)
            .map(|(vai, vti, uai, uti, rci)| ((vai + uai).powi(2) + (vti + uti + rci*w).powi(2)).sqrt())
            .collect::<Vec<_>>();

        Self{prop, inflow_vel: v}
    }

    pub fn xfoil_iteration(&mut self) {
        const VISCOSITY: f64 = 8.917e-7;
        let lift_coefficients = izip!(
            self.prop.hydro_data().circulation(), self.prop.geometry().chords(), self.inflow_vel.iter()
        ).map(|(g, c, v)| g/(v*c)).collect::<Vec<_>>();

        // result vectors
        let mut angles = Vec::<f64>::new();
        let mut drag_coeffs = Vec::<f64>::new();

        for (cl, c, v) in izip!(lift_coefficients, self.prop.geometry().chords(), self.inflow_vel.iter()) {
            // TODO: make viscosity configurable
            let reynolds = (v*c*self.prop.specs().ship_speed()*self.prop.geometry().radius()/VISCOSITY) as usize;

            // TODO: proper random name generation
            let fname = "foo";
            self.prop.geometry().base_airfoil().save_to_file(fname);

            // TODO: properly scale foil thickness
            //TODO: check for convergence
            let xfoil_result = rs_xfoil::Config::new("/usr/local/bin/xfoil")
                .airfoil_polar_file(fname)
                .reynolds(reynolds)
                .lift_coefficient(cl)
                .get_runner().unwrap()
                .dispatch()
                .unwrap();

            angles.push(xfoil_result["CL"][0]);
            drag_coeffs.push(xfoil_result["CD"][0]);
        }
        // TODO: extrapolate/interpolate cd/alpha at stations where xfoil failed
    }
}