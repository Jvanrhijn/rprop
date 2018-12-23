use coney::{
    ConeySolver,
    ConeySolverSingle
};
use propeller::{
    PropellerBuilder,
};


fn main() {
    let prop = PropellerBuilder::new(0.1, 0.025, 100.0, 200.0, 10.0, 10)
        .build();
    let coney_solver = ConeySolverSingle::new(prop).unwrap();
    let _prop = coney_solver.optimize_propulsor(1e-6).unwrap().remove(0);
}
