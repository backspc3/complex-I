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

    fn neg(self) -> Self::Output {
        Self {
            a: -self.a,
            b: -self.b
        }
    }
}

impl<T> Add for Complex<T> where T : Float {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            a: self.a + rhs.a,
            b: self.b + rhs.b
        }
    }
}

impl<T> AddAssign for Complex<T> where T : Float {
    fn add_assign(&mut self, rhs: Self) {
        self.a = self.a + rhs.a;
        self.b = self.b + rhs.b;
    }
}

impl<T> SubAssign for Complex<T> where T : Float {
    fn sub_assign(&mut self, rhs: Self) {
        self.a = self.a - rhs.a;
        self.b = self.b - rhs.b;
    }
}

impl<T> Sub for Complex<T> where T : Float {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            a: self.a - rhs.a,
            b: self.b - rhs.b
        }
    }
}

impl<T> Mul for Complex<T> where T : Float {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            a: self.a * rhs.a - self.b * rhs.b,
            b: self.a * rhs.b - rhs.a * self.b,
        }
    }
}

impl<T> MulAssign for Complex<T> where T : Float {
    fn mul_assign(&mut self, rhs: Self) {
        self.a = self.a * rhs.a - self.b * rhs.b;
        self.b = self.a * rhs.b - rhs.a * self.b;
    }
}

impl<T> Div for Complex<T> where T : Float {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        Self {
            a: (self.a * rhs.a + self.b * rhs.b) / ((rhs.a * rhs.a) + (rhs.b * rhs.b)),
            b: (self.b * rhs.a - self.a * rhs.b) / ((rhs.a * rhs.a) + (rhs.b * rhs.b))
        }
    }
}

impl<T> DivAssign for Complex<T> where T : Float {
    fn div_assign(&mut self, rhs: Self) {
        self.a = (self.a * rhs.a + self.b * rhs.b) / ((rhs.a * rhs.a) + (rhs.b * rhs.b));
        self.b = (self.b * rhs.a - self.a * rhs.b) / ((rhs.a * rhs.a) + (rhs.b * rhs.b));
    }
}

impl<T> Complex<T> where T : Float {
    pub fn zero() -> Self {
        Self { a: T::zero(), b: T::zero() }
    }

    pub fn new( a : T, b : T ) -> Self {
        Self { a: a, b: b }
    }

    pub fn scale(&self, scale : T) -> Self {
        Self { 
            a: self.a * scale , 
            b: self.b * scale 
        }
    }

    pub fn sqrt( &self ) -> Self {
        //                                  DOT PRODUCT
        let mut m = Float::sqrt( self.a * self.a + self.b * self.b );
        let mut angle = Float::atan2(self.b, self.a);
        m = Float::sqrt(m);
        // THIS IS CURSED... NO GOOD
        angle = angle / T::from(2.0).unwrap();
        Self { 
            a: m * Float::cos(angle), 
            b: m * Float::sin(angle) 
        }
    }
}

impl Complexf  {
    pub fn sqrd(&self) -> f32 {
        (self.a * self.a) - (self.b * self.b) + 2.0 * (self.a * self.b)
    }
}

impl Complexd {
    pub fn sqrd(&self) -> f64 {
        (self.a * self.a) - (self.b * self.b) + 2.0 * (self.a * self.b)
    }
}