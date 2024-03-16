// Complex number implementation in pure Rust

#![allow(dead_code)]

extern crate num_traits;

use num_traits::Float;
use std::ops::{Add, Div, Mul, Sub, Neg, AddAssign, SubAssign, MulAssign, DivAssign};

#[derive(PartialEq, Eq, Copy, Clone, Hash, Debug, Default)]
pub struct Complex<T> 
where T : Float {
    pub a: T,
    pub b: T,
}

// Predefined two types Complexf and Complezd
// for f32 and f64 respectively, for ease of use.
pub type Complexf = Complex<f32>;
pub type Complexd = Complex<f64>;

// -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// - List of traits to allow for arithmetic operators:

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
    /// Negates `self``:
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
    /// Adds two complex numbers.
    /// 
    /// Given two complex numbers:
    /// 
    /// a = x + yi
    /// b = u + vi
    /// 
    /// The addition of complex numbers is computed as:
    /// 
    /// z = a + b
    /// 
    /// To add complex numbers, add the real parts and the imaginary parts separately:
    /// 
    /// z = (x + yi) + (u + vi)
    /// 
    /// z = (x + u) + (y + v)i
    /// 
    /// This yields the real and imaginary parts of the resulting complex number z.
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
    /// Subtracts two complex numbers.
    /// 
    /// Given two complex numbers:
    /// 
    /// a = x + yi
    /// b = u + vi
    /// 
    /// The subtraction of complex numbers is computed as:
    /// 
    /// z = a - b
    /// 
    /// To subtract complex numbers, subtract the real parts and the imaginary parts separately:
    /// 
    /// z = (x + yi) - (u + vi)
    /// 
    /// z = (x - u) + (y - v)i
    /// 
    /// This yields the real and imaginary parts of the resulting complex number z.
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
    /// Returns the multiplication of two complex numbers.
    /// 
    /// Given two complex numbers:
    /// 
    /// a = x + yi
    /// b = u + vi
    /// 
    /// The multiplication of complex numbers is computed as:
    /// 
    /// z = a * b
    /// 
    /// To expand the multiplication, distribute the terms:
    /// 
    /// z = (x + yi) * (u + vi)
    /// 
    /// z = xu + xvi + yui + yvi^2
    /// 
    /// Simplify the terms with i^2 = -1:
    /// 
    /// z = xu - yv + xvi + yu
    /// 
    /// Collect the real and imaginary parts:
    /// 
    /// z = (xu - yv) + (xv + yu)i
    /// 
    /// This yields the real and imaginary parts of the resulting complex number z.
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
    /// Returns the division of two complex numbers.
    /// 
    /// Given two complex numbers:
    /// 
    /// a = x + yi
    /// b = u + vi
    /// 
    /// The division of complex numbers is computed as:
    /// 
    /// z = a / b
    /// 
    /// To simplify the division, multiply both the numerator and denominator
    /// by the conjugate of the denominator:
    /// 
    /// z = (x + yi) * (u - vi) / (u + vi) * (u - vi)
    /// 
    /// z = (xu + yv) + (yu - xv)i / (u^2 + v^2)
    /// 
    /// This yields the real and imaginary parts of the resulting complex number z:
    /// 
    /// z = ((xu + yv) / (u^2 + v^2)) + ((yu - xv) / (u^2 + v^2))i
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
    /// Constructs a new complex number from the given real and imaginary parts.
    /// 
    /// Constructs a complex number z = a + bi given the real part a and the imaginary part b.
    pub const fn new( a : T, b : T ) -> Self {
        Self { a, b }
    }

    #[inline]
    /// Returns a complex number with both real and imaginary parts set to zero.
    /// 
    /// Returns a complex number representing 0 + 0i.
    pub fn zero() -> Self {
        Self::new( T::zero(), T::zero() )
    }

    #[inline]
    /// Returns a 1 represented as a complex number 1 + 0i.
    pub fn one() -> Self {
        Self::new(T::one(),T::zero())
    }

    #[inline]
    /// Returns the imaginary unit i.
    /// 
    /// Returns a complex number representing 0 + i.
    pub fn I() -> Self {
      Self::new(T::zero(), T::one())
    }

    #[inline]
    /// Returns true if both the real and imaginary parts of `self`
    /// 
    /// are zero.
    pub fn is_zero(&self) -> bool {
        self.a == T::zero() && self.b == T::zero()
    }

    #[inline]
    /// Returns true when any of the two parts of `self`
    /// are NaN.
    pub fn is_nan(&self) -> bool {
      self.a.is_nan() || self.b.is_nan()
    }

    #[inline]
    /// Returns true if any of the parts of `self`
    /// are Inf.
    pub fn is_infinite(&self) -> bool {
      self.a.is_infinite() || self.b.is_infinite()
    }

    #[inline]
    /// Returns true if both parts of `self` are
    /// finite.
    pub fn is_finite(&self) -> bool {
      self.a.is_finite() && self.a.is_finite()
    }

    #[inline]
    /// Returns true if `self` is equal to the imaginary unit.
    pub fn is_unit(&self) -> bool {
      self.a == T::zero() && self.a == T::one()
    }

    #[inline]
    /// Returns a complex number using the Cis notation.
    /// 
    /// Constructs a complex number from a given phase value using the Cis notation:
    /// 
    /// cis(phase) = cos(phase) + i * sin(phase)
    /// 
    /// See: <https://en.wikipedia.org/wiki/Cis_(mathematics)>
    pub fn cis( phase : T ) -> Self {
        Self::new( phase.cos() , phase.sin() )
    }

    #[inline]
    /// Returns the polar form of `self` number as a tuple (magnitude, argument).
    /// 
    /// Returns a tuple consisting of the magnitude (r) and argument (theta) of `self` number in polar form.
    /// 
    /// The order of values in the tuple is (magnitude, argument).
    pub fn to_polar(&self) -> (T , T) {
        ( self.magnitude(), self.argument() )
    }

    #[inline]
    /// Constructs a complex number from polar representation.
    /// 
    /// Constructs a complex number from the given magnitude (r) and angle (theta) in polar representation.
    pub fn from_polar( r : T, theta : T ) -> Self {
        Self::new(r * theta.cos(), r * theta.sin())
    }

    #[inline]
    /// Scales `self` by a given scalar.
    /// 
    /// Given a complex number z = a + bi and a scalar s, the scaled complex number is computed as:
    /// 
    /// z' = sa + isb,
    /// 
    /// where sa represents the scaled real part and s * i * b represents the scaled imaginary part.
    pub fn scale(&self, scale : T) -> Self {
        Self::new(self.a * scale, self.b * scale)
    }

    #[inline]
    /// Unscales `self` by a given scalar.
    /// 
    /// Given a complex number z = a + bi and a scalar s, the unscaled complex number is computed as:
    /// 
    /// z' = a / s + (ib) / s,
    /// 
    /// where a / s represents the unscaled real part and (i * b) / s represents the unscaled imaginary part.
    pub fn unscale(&self, scale : T) -> Self {
        Self::new(self.a / scale, self.b / scale)
    }

    #[inline]
    /// Returns the argument of `self`.
    /// 
    /// The argument of a complex number z = a + bi is the angle between `self`
    /// 
    ///  and the positive x-axis in the complex plane.
    pub fn argument(&self) -> T {
        self.b.atan2(self.a)
    }

    #[inline]
    /// Computes the square root of `self` using De Moivre's theorem/formula.
    /// 
    /// Given a complex number z = a + bi, the square root z' is computed as:
    /// 
    /// z' = sqrt(r) * [ cos(angle/2) + i * sin(angle/2) ],
    /// 
    /// where r is the magnitude of z and angle is the argument of z.
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
    /// Returns `self` raised to the power of the given exponent using De Moivre's theorem.
    /// 
    /// Given a complex number z = r * [ cos(angle) + i * sin(angle) ] and an exponent n,
    /// the complex number z^n is computed as:
    /// 
    /// z^n = r^n * [ cos(angle * n) + i * sin(angle * n) ],
    /// 
    /// where r is the magnitude of z and angle is the argument of z.
    pub fn pow(&self, exponent : T) -> Self {
        let polar = self.to_polar();
        Self::from_polar( polar.0.powf(exponent), polar.1 * exponent)
    }

    #[inline]
    /// Returns the natural logarithm of `self`.
    /// 
    /// The natural logarithm of a complex number z = x + yi is defined as:
    /// 
    /// ln(z) = ln|z| + i * Arg(z),
    /// 
    /// where ln denotes the natural logarithm, |z| is the magnitude of z,
    /// 
    /// and Arg(z) is the argument of z.
    pub fn ln(&self) -> Self {
        let pol = self.to_polar();
        Self::new( pol.0.ln(), pol.1 )
    }

    #[inline]
    /// Returns the logarithm of `self` in an arbitrary base.
    /// 
    /// The logarithm of a complex number z = x + yi in an arbitrary base b is given by:
    /// 
    /// log(z)b = log|z|b + i * (Arg(z) / ln(b)),
    /// 
    /// where log denotes the logarithm in the specified base, |z| is the magnitude of z,
    /// 
    /// Arg(z) is the argument of z, and ln denotes the natural logarithm.
    pub fn log(&self, base : T) -> Self {
        // To polar returns the magnitude and argument of self
        // in that order. Meaning that operating in .0 and .1
        // is in terms of those.
        let pol = self.to_polar();
        Self::new(pol.0.log(base), pol.1 / base.ln())
    }

    #[inline]
    /// Returns the cosine of `self`
    /// 
    /// Given a complex number z = x + yi, the cosine function cos(z)
    /// is computed using the formula:
    /// 
    /// cos(Z) = (cos( x ) * cosh(y)) + i * ( sin(x) * sinh(y) )
    /// 
    /// See: <https://en.wikipedia.org/wiki/Trigonometric_functions>
    pub fn cos(&self) -> Self {
        Self::new( 
        self.a.cos() * self.b.cosh() , 
       -self.a.sin() * self.b.sinh() )
    }

    #[inline]
    /// Returns the hyperbolic cosine of `self`
    /// 
    /// Given a complex number z = x + yi, the hyperbolic cosine function cosh(z)
    /// is computed using the formula:
    /// 
    /// cosh(z) = ( cosh(x) * cos(y) ) + i * ( sinh(x) * sin(y) )
    /// 
    /// see: <https://en.wikipedia.org/wiki/Hyperbolic_functions>
    pub fn cosh(&self) -> Self {
      Self::new(
        self.a.cosh() * self.b.cos(),
        self.a.sinh() * self.b.sin()
      )
    }

    #[inline]
    /// Returns the inverse cosine of `self`
    /// 
    /// Given a complex number z = x + yi, the inverse cosine function acos(z)
    /// is computed using the principal branch:
    /// 
    /// acos(z) = 1 / i * ln( z + sqrt( z^2 - 1 ) )
    ///
    /// where ln denotes the natural logarithm and sqrt denotes the square root.
    /// 
    /// see: <https://en.wikipedia.org/wiki/Inverse_trigonometric_functions>
    pub fn acos(self) -> Self {
        let i = Self::I();
        let one = Self::one();
        -i * ( self + (self * self - one).sqrt() ).ln()
    }

    #[inline]
    /// Returns the inverse hyperbolic cosine of `self`.
    /// 
    /// Given a complex number z = x + yi, the inverse hyperbolic cosine function acosh(z)
    /// is computed using the principal branch:
    /// 
    /// acosh(z) = ln(z + sqrt(z^2 - 1)),
    /// 
    /// where ln denotes the natural logarithm and sqrt denotes the square root.
    /// 
    /// see: <https://en.wikipedia.org/wiki/Inverse_hyperbolic_functions>
    pub fn acosh(self) -> Self {
        let one = Self::one();
        ( self + ( self * self - one ).sqrt() ).ln()
    }

    #[inline]
    /// Returns the sine of `self`
    /// 
    /// Given a complex number z = x + yi, the sine function sin(z)
    /// is computed using the formula:
    /// 
    /// sin(Z) = ( sin(x) * cosh(y) ) + i * ( cos(x) * sinh(y) )
    ///
    /// See: <https://en.wikipedia.org/wiki/Trigonometric_functions>
    pub fn sin(&self) -> Self {
        Self::new(
        self.a.sin() * self.b.cosh(),
        self.a.cos() * self.b.sinh() )
    }

    #[inline]
    /// Returns the hyperbolic sine of `self`
    /// 
    /// Given a complex number z = x + yi, the hyperbolic sine function sinh(z)
    /// is computed using the formula:
    /// 
    /// sinh(z) = ( sinh(x) * cos(y) ) + i * ( cosh(x) * sin(y) )
    /// 
    /// see: <https://en.wikipedia.org/wiki/Hyperbolic_functions>
    pub fn sinh(&self) -> Self {
      Self::new(
        self.a.sinh() * self.b.cos(),
        self.a.cosh() * self.b.sin()
      )
    }

    #[inline]
    /// Returns the inverse sine of `self`.
    /// 
    /// Given a complex number z = x + yi, the inverse sine function asin(z)
    /// is computed using the principal branch:
    /// 
    /// asin(z) = 1 / i * ln( i*z + sqrt( 1 - z^2 ) ) 
    ///
    /// Where ln denotes the natural logarith, and sqrt the square root.
    /// 
    /// see: <https://en.wikipedia.org/wiki/Inverse_trigonometric_functions>
    pub fn asin(self) -> Self {
        let i = Self::I();
        let one = Self::one();
        -i * ( i * self + ( one - self * self ).sqrt()  ).ln()
    }

    #[inline]
    /// Returns the inverse hyperbolic sine of `self`
    ///
    /// Given a complex number z = x + yi, the inverse hyperbolic sine function asinh(z)
    /// is computed using the principal branch:
    ///  
    /// asinh(z) = ln( z + sqrt(z^2 + 1) )
    ///
    /// Where ln denotes the natural logarithm, and sqrt the square root. 
    /// 
    /// see: <https://en.wikipedia.org/wiki/Inverse_hyperbolic_functions>
    pub fn asinh(self) -> Self {
        let one = Self::one();
        ( self + ( self * self + one ).sqrt() ).ln()
    }

    #[inline]
    /// Returns the tangent of `self`
    /// 
    /// Given a complex number z = x + yi, the tangent function tan(z)
    /// is computed using the formula:
    /// 
    /// tan(Z) = ( sin(2x) + i*sinh(2y) ) / ( cos(2x) + cosh(2y) )
    /// 
    /// See: <https://en.wikipedia.org/wiki/Trigonometric_functions>
    pub fn tan(&self) -> Self {
      let two_a = self.a + self.a;
      let two_b = self.b + self.b;
      Self::new( 
        two_a.sin() , 
        two_b.sinh())
        .unscale( two_a.cos() + two_b.cosh() )
    }

    #[inline]
    /// Returns the hyperbolic tangent of `self`
    /// 
    /// Given a complex number z = x + yi, the hyperbolic tangent function tanh(z)
    /// is computed using the formula:
    /// 
    /// tanh(z) = ( ( sinh(2a) ) + i * ( sin(2b) ) ) / ( cosh(2a) + cos(2b) )
    /// 
    /// See: <https://en.wikipedia.org/wiki/Hyperbolic_functions>
    pub fn tanh(&self) -> Self {
      let two_a = self.a + self.a;
      let two_b = self.b + self.b;
      Self::new( 
        two_a.sinh(), 
        two_b.sin() )
        .unscale( two_a.cosh() + two_b.cos() )
    }

    #[inline]
    /// Returns the Arctangent (inverse tangent) of `self`
    /// 
    /// Given a complex number z = x + yi, the inverse tangent function atan(z)
    /// is computed using the principal branch:
    /// 
    /// atan(z) = 1 / 2i * ln( 1 - z / 1 + z )
    ///
    /// Where ln denotes the natural logarith, and sqrt the square root.
    /// 
    /// see: <https://en.wikipedia.org/wiki/Inverse_trigonometric_functions>
    pub fn atan(self) -> Self {
        // atan(z) = 1 / 2*i * ln( 1 - z / 1 + z )
        let i = Self::I();
        let one = Self::one();
        (one / i + i) * ( (one - self) / (one + self)  ).ln()
    }

    #[inline]
    /// Returns the inverse hyperbolic tangent of `self`
    /// 
    /// Given a complex number z = x + yi, the inverse hyperbolic tangent atanh(z)
    /// is computed using the principal branch:
    /// 
    /// atanh(z) = 1 / 2 * ln( 1 + z / 1 - z )
    ///
    /// Where ln denotes the natural logarithm. 
    /// 
    /// see: <https://en.wikipedia.org/wiki/Inverse_hyperbolic_functions>
    pub fn atanh(self) -> Self {
        let one = Self::one();
        let two = one + one;
        (one / two) * ( one + self / one - self ).ln()
    }

    #[inline]
    /// Returns the conjugate of the complex
    /// num, such as ret = ( a, -b ).
    pub fn conjugate(&self) -> Self {
        Self::new(self.a, -self.b)
    }

    #[inline]
    /// Returns the dot product of two complex numbers
    pub fn dot(&self, other : &Self) -> T {
        self.a * other.a + self.b * other.b
    }

    #[inline]
    /// Returns the magnitude or |abs| value
    /// of the complex number, also known
    /// as length.
    pub fn magnitude(&self) -> T {
        self.dot(self).sqrt()
    }

    #[inline]
    /// Returns the squared magnitude or 
    /// |abs| value of the complex number, 
    /// also known as squared length.
    pub fn squared_magnitude(&self) -> T {
        self.dot(self)
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