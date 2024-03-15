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
        
        let sqr = testa.sqrt();
        let mag = testa.magnitude();

        let pow = testa.pow(2.0);

        let pol = testa.to_polar();

        println!("mul {:?}", mul);
        println!("div {:?}", div);
        println!("add {:?}", add);
        println!("res {:?}", res);
        println!("sqr {:?}", sqr);
        println!("mag {:?}", mag);
        println!("pow {:?}", pow);
        println!("pol {:?}", pol);
    }
}
