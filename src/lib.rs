#![allow(non_snake_case)]
pub mod complex;

#[cfg(test)]
mod tests {
    use crate::complex::Complex;
    use crate::complex::Complexf;

    #[test]
    fn it_works() {
        let testa : Complexf = Complex::new( 1.0, 3.0 );
        let testb : Complexf = Complex::new( 4.0, 5.0 );
        
        let mul = testa * testb;

        let div = testa / testb;

        let add = testa + testb;

        let res = testa - testb;
        
        //let sqr = Complex::<f32>::sqrd(&testa);
        let sqr = testa.sqrt();

        println!("{:?}", mul);
        println!("{:?}", div);
        println!("{:?}", add);
        println!("{:?}", res);
        println!("{:?}", sqr);
    }
}
