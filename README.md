# complex-I
A complex number library in rust.

## Overview

This library provides a type complex for floating point numbers
which offers various common mathematical operations. 

see docs for all available operations.

## Example

```rust
    use complex::Complex;
    use complex::Complexf;

    fn main() {

        // Create two complex numbers.
        let complexa : Complexf = Complex::new( 1.0, 3.0 ); 
        let complexb : Complexf = Complex::new( 2.0, 1.0 );

        // Can use these types with various operators
        // and transformations.
        let add = complexa + complexb;
        let mul = complexa * complexb;
    
        complexa /= complexb;

        let cos = complexa.cos();
        // Check docs for all available operations and transformations.

        println!("Final: {:?}", cos);
    }
```

## Resources:
- Wikipedia:\
     [Complex numbers](https://en.wikipedia.org/wiki/Complex_number)\
     [De Moivre's formula](https://en.wikipedia.org/wiki/De_Moivre%27s_formula)\
     [Trigonometric functions](https://en.wikipedia.org/wiki/Trigonometric_functions)\
     [Hyperbolic functions](https://en.wikipedia.org/wiki/Hyperbolic_functions)\
     [Inverse trig functions](https://en.wikipedia.org/wiki/Inverse_trigonometric_functions)\
     [Inverse hyperbolic functions](https://en.wikipedia.org/wiki/Inverse_hyperbolic_functions)
- Proofwiki:\
     [proofs/definitions](https://proofwiki.org/wiki/Main_Page)
