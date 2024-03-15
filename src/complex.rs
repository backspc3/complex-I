// Complex number implementation in pure rust
// https://en.wikipedia.org/wiki/Complex_number

#![allow(dead_code)]

extern crate num_traits;

use num_traits::Float;
use std::ops::{Add, Div, Mul, Sub, Neg, AddAssign, SubAssign, MulAssign, DivAssign};

#[derive(Debug, Clone, Copy)]
pub struct Complex<T> 
where T : Float {
    a: T,
    b: T,
}

pub type Complexf = Complex<f32>;
pub type Complexd = Complex<f64>;

impl<T> Neg for Complex<T> where T : Float {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self {
            a: -self.a,
            b: -self.b
        }
    }
}

impl<T> Add for Complex<T> where T : Float {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            a: self.a + rhs.a,
            b: self.b + rhs.b
        }
    }
}

impl<T> AddAssign for Complex<T> where T : Float {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.a = self.a + rhs.a;
        self.b = self.b + rhs.b;
    }
}

impl<T> SubAssign for Complex<T> where T : Float {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.a = self.a - rhs.a;
        self.b = self.b - rhs.b;
    }
}

impl<T> Sub for Complex<T> where T : Float {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            a: self.a - rhs.a,
            b: self.b - rhs.b
        }
    }
}

impl<T> Mul for Complex<T> where T : Float {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            a: self.a * rhs.a - self.b * rhs.b,
            b: self.a * rhs.b - rhs.a * self.b,
        }
    }
}

impl<T> MulAssign for Complex<T> where T : Float {

    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        self.a = self.a * rhs.a - self.b * rhs.b;
        self.b = self.a * rhs.b - rhs.a * self.b;
    }
}

impl<T> Div for Complex<T> where T : Float {
    type Output = Self;

    #[inline]
    fn div(self, rhs: Self) -> Self::Output {
        Self {
            a: (self.a * rhs.a + self.b * rhs.b) / ((rhs.a * rhs.a) + (rhs.b * rhs.b)),
            b: (self.b * rhs.a - self.a * rhs.b) / ((rhs.a * rhs.a) + (rhs.b * rhs.b))
        }
    }
}

impl<T> DivAssign for Complex<T> where T : Float {

    #[inline]
    fn div_assign(&mut self, rhs: Self) {
        self.a = (self.a * rhs.a + self.b * rhs.b) / ((rhs.a * rhs.a) + (rhs.b * rhs.b));
        self.b = (self.b * rhs.a - self.a * rhs.b) / ((rhs.a * rhs.a) + (rhs.b * rhs.b));
    }
}

impl<T> Complex<T> where T : Float {

    #[inline]
    // Returns a zeroed complex number meaning
    // 0 + 0 * i
    pub fn zero() -> Self {
        Self { a: T::zero(), b: T::zero() }
    }

    #[inline]
    // Functions as a 'constructor' type to construct
    // a complex number given the real and imaginary parts.
    pub fn new( a : T, b : T ) -> Self {
        Self { a: a, b: b }
    }

    #[inline]
    // Returns a complex nubmer using Cis notation
    // -> exp( b * phase )
    // https://es.wikipedia.org/wiki/Cis_(matem%C3%A1ticas)
    pub fn cis( phase : T ) -> Self {
        Self::new( phase.cos() , phase.sin() )
    }

    #[inline]
    // Returns a valid complex number from polar representation
    pub fn from_polar( r : T, theta : T ) -> Self {
        Self {
            a: r * theta.cos(),
            b: r * theta.sin()
        }
    }

    #[inline]
    // Scales the complex number by a given scalar.
    pub fn scale(&self, scale : T) -> Self {
        Self { 
            a: self.a * scale , 
            b: self.b * scale 
        }
    }

    #[inline]
    // Returns the angle of the vector using atan2
    pub fn arg(&self) -> T {
        self.a.atan2(self.b)
    }

    #[inline]
    // Computes the square root of the complex number
    pub fn sqrt( &self ) -> Self {
        //                                  DOT PRODUCT
        let mut m =( self.a * self.a + self.b * self.b ).sqrt();
        let mut angle = self.b.atan2(self.a);
        m = m.sqrt();
        // THIS IS CURSED... NO GOOD I THINK
        angle = angle / T::from(2.0).unwrap();
        Self { 
            a: m * angle.cos(), 
            b: m * angle.sin() 
        }
    }

    #[inline]
    /// Returns the conjugate of the complex
    /// num, such as ret = ( a, -b ).
    pub fn conjugate(&self) -> Self {
        Self { 
            a: self.a, 
            b: -self.b 
        }
    }

    #[inline]
    /// Returns the magnitude or |abs| value
    /// of the complex number, also known
    /// as length.
    pub fn magnitude(&self) -> T {
        (self.a * self.a + self.b * self.b).sqrt()
    }

    #[inline]
    /// Returns a tuple consisting of the 
    /// complex numbers mag and angle
    /// required to use polar form.
    pub fn polar(&self) -> (T , T) {
        ( self.magnitude(), self.arg() )
    }

    #[inline]
    /// Returns the squared magnitude or 
    /// |abs| value of the complex number, 
    /// also known as length.
    pub fn squared_magnitude(&self) -> T {
        self.a * self.a + self.b * self.b
    }
    
}