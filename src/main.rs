use coney::{
    ConeySolver,
    ConeySolverSingle
};
use propeller::{
    PropellerBuilder,
};


fn main() {
    let prop = PropellerBuilder::new(1.0, 0.25, 1000.0, 200.0, 10.0, 20)
        .build();
    let coney_solver = ConeySolverSingle::new(prop).unwrap();
    let _prop = coney_solver.optimize_propulsor(1e-6).unwrap().remove(0);
}
