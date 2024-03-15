// Complex number implementation in pure rust

// -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// This library defines a type Complex<T>, where T has to be have the trait Float
// for floating point numbers.

// Inspired by Daniel Shiffmans's Pi Day video:
// - https://www.youtube.com/watch?v=6UlGLB_jiCs

// An exercise in both complex numbers and rust.

// -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

// Doc:
// - https://en.wikipedia.org/wiki/Complex_number
// - https://en.wikipedia.org/wiki/De_Moivre%27s_formula

#![allow(dead_code)]

extern crate num_traits;

use num_traits::Float;
use std::ops::{Add, Div, Mul, Sub, Neg, AddAssign, SubAssign, MulAssign, DivAssign};

#[derive(PartialEq, Eq, Copy, Clone, Hash, Debug, Default)]
pub struct Complex<T> 
where T : Float {
    a: T,
    b: T,
}

pub type Complexf = Complex<f32>;
pub type Complexd = Complex<f64>;

// -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// - List of traits to allow for arithmetic using operators:

// complexa + complexb
// complexa += complexb
// complexa * complexb
// complexa *= compexb

// complexa - complexb
// complexa /= compexb
// complexa / complexb
// complexa /= compexb

// -complexa
// -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

impl<T> Neg for Complex<T> where T : Float {
    type Output = Self;

    #[inline]
    /// Negates self:
    /// 
    /// a = x + yi
    /// 
    /// -a = -a - yi
    /// 
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
    /// Adds two complex numbers
    /// given by:
    /// 
    /// a = x + yi
    /// b = u + vi
    /// 
    /// z = a + b
    /// 
    /// z = (x + yi) + (u + vi)
    /// 
    /// z = (x + u) + (y + v)i
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            a: self.a + rhs.a,
            b: self.b + rhs.b
        }
    }
}

impl<T> AddAssign for Complex<T> where T : Float {
    #[inline]
    /// See definition of add.
    fn add_assign(&mut self, rhs: Self) {
        self.a = self.a + rhs.a;
        self.b = self.b + rhs.b;
    }
}

impl<T> SubAssign for Complex<T> where T : Float {
    #[inline]
    /// See definition of sub.
    fn sub_assign(&mut self, rhs: Self) {
        self.a = self.a - rhs.a;
        self.b = self.b - rhs.b;
    }
}

impl<T> Sub for Complex<T> where T : Float {
    type Output = Self;

    #[inline]
    /// Subtracts two complex numbers
    /// given by:
    /// 
    /// a = x + yi
    /// b = u + vi
    /// 
    /// z = a - b
    /// 
    /// z = (x + yi) - (u + vi)
    /// 
    /// z = (x - u) + (y - v)i
    /// 
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
    /// Returns the multiplication of two complex numbers
    /// given by:
    /// 
    /// a = x + yi
    /// b = u + vi
    /// 
    /// z = a * b
    /// 
    /// z = (x + yi) * (u + vi)
    /// 
    /// expanded to:
    /// 
    /// z = xu + xvi + yui + yvi^2
    /// 
    /// then:
    /// 
    /// z = ( xu - yv ) + ( xv + yu )i
    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            a: self.a * rhs.a - self.b * rhs.b,
            b: self.a * rhs.b + self.b * rhs.a,
        }
    }
}

impl<T> MulAssign for Complex<T> where T : Float {

    #[inline]
    /// See definition of mul.
    fn mul_assign(&mut self, rhs: Self) {
        self.a = self.a * rhs.a - self.b * rhs.b;
        self.b = self.a * rhs.b + self.b * rhs.a;
    }
}

impl<T> Div for Complex<T> where T : Float {
    type Output = Self;

    #[inline]
    /// Returns the divison of two complex numbers
    /// given by:
    /// 
    /// a = x + yi
    /// b = u + vi
    /// 
    /// z = a / b
    /// z = (x + yi) / (u + vi)
    /// 
    /// to simplify divide both numerator and denominator
    /// by the conjugate of the denominator
    /// 
    /// z = (x + yi) * (u - vi) / (u + vi) * (u - vi)
    /// 
    /// z = ( xu + yv ) + ( yu - xv )i / (u^2 + v^2)
    /// 
    /// which means
    /// 
    ///             REAL PART                       IMAGINARY PART
    /// z = ( ( xu + yv ) / (u^2 + v^2) )  + ( ( yu - xv ) / (u^2 + v^2) ) i; 
    /// 
    fn div(self, rhs: Self) -> Self::Output {
        Self {
            a: (self.a * rhs.a + self.b * rhs.b) / ((rhs.a * rhs.a) + (rhs.b * rhs.b)),
            b: (self.b * rhs.a - self.a * rhs.b) / ((rhs.a * rhs.a) + (rhs.b * rhs.b))
        }
    }
}

impl<T> DivAssign for Complex<T> where T : Float {

    #[inline]
    /// See definition of div.
    fn div_assign(&mut self, rhs: Self) {
        self.a = (self.a * rhs.a + self.b * rhs.b) / ((rhs.a * rhs.a) + (rhs.b * rhs.b));
        self.b = (self.b * rhs.a - self.a * rhs.b) / ((rhs.a * rhs.a) + (rhs.b * rhs.b));
    }
}

// -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

// Implementation body for methods and transformations that can be applied to
// a complex number.

impl<T> Complex<T> where T : Float {

    #[inline]
    /// Functions as a 'constructor' type to construct
    /// a complex number given the real and imaginary parts.
    pub fn new( a : T, b : T ) -> Self {
        Self { a: a, b: b }
    }

    #[inline]
    /// Returns a zeroed complex number meaning
    /// 0 + 0 * i
    pub fn zero() -> Self {
        Self::new( T::zero(), T::zero() )
    }

    #[inline]
    /// Returns a true when both the imaginary and real
    /// parts of a complex number are zero... false when they
    /// are different than zero.
    pub fn is_zero(&self) -> bool {
        self.a == T::zero() && self.b == T::zero()
    }


    #[inline]
    /// Returns a complex number by inputting a phase value 
    /// using Cis notation:
    /// cis(phase) = cos(phase) + i * sin(phase)
    /// https://en.wikipedia.org/wiki/Cis_(mathematics)
    pub fn cis( phase : T ) -> Self {
        Self::new( phase.cos() , phase.sin() )
    }

    #[inline]
    /// Returns a tuple consisting of the 
    /// complex numbers polar form:
    /// - mag.
    /// - angle.
    /// 
    /// In that order
    pub fn to_polar(&self) -> (T , T) {
        ( self.magnitude(), self.angle() )
    }

    #[inline]
    /// Returns a valid complex number from polar representation
    pub fn from_polar( r : T, theta : T ) -> Self {
        Self::new(r * theta.cos(), r * theta.sin())
    }

    #[inline]
    /// Scales the complex number by a given scalar.
    /// meaning:
    /// scale * re + i * scale * im
    pub fn scale(&self, scale : T) -> Self {
        Self::new(self.a * scale, self.b * scale)
    }

    #[inline]
    /// Unscaled the complex number by a given scalar.
    /// meaning:
    /// re / scale + (i * im) / scale
    pub fn unscale(&self, scale : T) -> Self {
        Self::new(self.a / scale, self.b / scale)
    }

    #[inline]
    /// Returns the angle between the cartesian representation
    /// of the complex number and the x plane or real plane
    pub fn angle(&self) -> T {
        self.b.atan2(self.a)
    }

    #[inline]
    /// Computes the square root of the complex number
    /// 
    /// Using De Moivre's theorem/formula.
    pub fn sqrt( &self ) -> Self {
        // If the imaginary number is zero, skip computations,
        // and simply return a zeroed complex number.
        if self.is_zero() {
            return Self::zero();
        }

        // if not zero, compute the polar form of the complex number
        let (mut m, mut angle) = self.to_polar();
        // and using De Moivre's theorem/formula
        // compute the sqrt

        // we need the square root of the numbers's
        // magnitude sqrt(mag).
        m = m.sqrt();
        // and half of the numbers angle with the x plane
        // angle / 2
        // THIS IS CURSED... NO GOOD I THINK
        angle = angle / T::from(2).unwrap();
        // then convert back to complex notation
        Self::from_polar(m, angle)
    }

    #[inline]
    /// Returns self^n as a complex number
    /// 
    /// Using De Moivre's theorem
    /// 
    /// z^n = r^n [ cos(On) + i sin(On) ]
    /// 
    /// EXPANDED
    /// 
    /// z^n = r^n * cos(O * n) + i * r^n * sin(O * n)
    pub fn pow(&self, exponent : T) -> Self {
        let polar = self.to_polar();
        Self::from_polar( polar.0.powf(exponent), polar.1 * exponent)
    }

    #[inline]
    /// Returns the conjugate of the complex
    /// num, such as ret = ( a, -b ).
    pub fn conjugate(&self) -> Self {
        Self::new(self.a, -self.b)
    }

    #[inline]
    /// Returns the magnitude or |abs| value
    /// of the complex number, also known
    /// as length.
    pub fn magnitude(&self) -> T {
        (self.a * self.a + self.b * self.b).sqrt()
    }

    #[inline]
    /// Returns the squared magnitude or 
    /// |abs| value of the complex number, 
    /// also known as squared length.
    pub fn squared_magnitude(&self) -> T {
        self.a * self.a + self.b * self.b
    }

    #[inline]
    /// Returns the imaginary component 
    pub const fn im(&self) -> T {
        self.b
    }

    #[inline]
    // Returns the real component
    pub const fn real(&self) -> T {
        self.a
    } 
    
}