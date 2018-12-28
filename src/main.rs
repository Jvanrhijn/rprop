#[macro_use]
extern crate itertools;
use coney::{
    ConeySolver,
    ConeySolverSingle
};
use propeller::{
    PropellerBuilder,
};

mod geometry;


fn main() {
    let prop = PropellerBuilder::new(0.1, 0.025, 100.0, 200.0, 10.0, 50)
        .build();
    let coney_solver = ConeySolverSingle::new(prop);
    let prop = coney_solver.optimize_propulsor(1e-6).unwrap().remove(0);
    println!("{:?}", prop.hydro_data().circulation());
}
