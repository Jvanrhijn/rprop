use airfoil::Airfoil;

fn main() {
    let foil = Airfoil::from_file("NACA2414.dat").unwrap();
    println!("{}", foil.area().unwrap());
}
