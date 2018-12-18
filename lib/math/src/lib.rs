use std::ops::Mul;

pub fn square<T: Mul<Output=T> + Copy>(x: T) -> T {
    x*x
}

#[cfg(test)]
mod tests {
    use super::square;
    #[test]
    fn it_works() {
        assert_eq!(square(1), 1);
        assert_eq!(square(2), 4);
    }
}
