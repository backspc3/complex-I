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

        let lnc = testa.ln();
        let log = testa.log(2.0);
        let cos = testa.cos();
        let sin = testa.sin();

        // Right now I have to figure out how to do a test
        // correctly in rust, so I am basically just hand checing
        // the answers with this calculator:
        // <https://www.calculators-math.com/complex-number-calculator/>
        println!("mul {:?}", mul);
        println!("div {:?}", div);
        println!("add {:?}", add);
        println!("res {:?}", res);
        println!("sqr {:?}", sqr);
        println!("mag {:?}", mag);
        println!("pow {:?}", pow);
        println!("pol {:?}", pol);
        println!("pol {:?}", lnc);
        println!("pol {:?}", log);
        println!("pol {:?}", cos);
        println!("pol {:?}", sin);
    }
}
