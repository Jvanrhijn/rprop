use std::vec::Vec;

pub struct Propeller {
    
}

struct Geometry {
    radius: f64,
    hub_radius: f64,
    chords: Vec<f64>,
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
