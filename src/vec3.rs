use std::ops::*;
use std::f64;

#[derive(Debug, Copy, Clone)]
pub struct Vec3 {
    pub x: f64, 
    pub y: f64, 
    pub z: f64 
}


impl Vec3 {
    pub fn new(x: f64, y: f64, z: f64) -> Self { 
        Vec3{x, y, z} 
    }

    pub fn zero() -> Self { 
        Vec3{x: 0., y: 0., z : 0.} 
    }

    pub fn length_sqr(&self) -> f64 {
        self.x*self.x + self.y*self.y + self.z*self.z
    }

    pub fn length(&self) -> f64 {
        f64::sqrt(self.length_sqr())
    }

    pub fn normalize(self) -> Self {
        self / self.length()
    }

    pub fn abs(self) -> Self {
        Vec3::new(f64::abs(self.x), f64::abs(self.y), f64::abs(self.z))
    }

    pub fn dot(self, other : Vec3) -> f64 {
        self.x*other.x + self.y*other.y + self.z*other.z
    }

    pub fn cross(self, other : Vec3) -> Self {
        Vec3::new(self.y*other.z - self.z*other.y,
                  self.z*other.x - self.x*other.z,
                  self.x*other.y - self.y*other.x)
    }

    pub fn rotate_to(self, normal : &Vec3) -> Self {
        // Given a direction vector self sampled around the z-axis of a
        // local coordinate system, this function applies the same
        // rotation to self as is needed to rotate the z-axis to the
        // actual direction n that v should have been sampled around
        // [Frisvad, Journal of Graphics Tools 16, 2012;
        //  Duff et al., Journal of Computer Graphics Techniques 6, 2017].
        let sign = if normal[2].is_sign_positive() { 1. } else { -1.};
        let a = -1./(1. + f64::abs(normal[2]));
        let b = normal[0]*normal[1]*a;
        Vec3::new(1. + normal[0]*normal[0]*a, b, 
                  -sign*normal[0])*self[0] + Vec3::new(sign*b, sign*(1. + normal[1]*normal[1]*a), 
                  -normal[1])
                      *self[1] + *normal*self[2]
    }

    //Create a coodinate system from a single vector.
    pub fn create_tangent_vectors(&self) -> (Vec3, Vec3) {
        let v2 = if f64::abs(self.x) > f64::abs(self.y) {
            Vec3::new(-self.z, 0., self.x) /
                f64::sqrt(self.x * self.x + self.z * self.z)
        } else {
            Vec3::new(0., self.z, -self.y) /
                f64::sqrt(self.y * self.y + self.z * self.z)
        };
        let v3 = self.cross(v2);

        (v2, v3)
    }
}

impl From<[f32; 3]> for Vec3 {
    fn from(array : [f32; 3]) -> Self {
        Vec3::new(array[0].into(), array[1].into(), array[2].into()) 
    }
}

impl Sub for Vec3 {
    type Output = Vec3;
    fn sub(self, other: Vec3) -> Self {
        Vec3::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

impl Neg for Vec3 {
    type Output = Vec3;

    fn neg(self) -> Self {
        Vec3::new(-self.x, -self.y, -self.z)
    }
}

impl Add for Vec3 {
    type Output = Vec3;
    fn add(self, other: Vec3) -> Self {
        Vec3::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

impl Mul for Vec3 {
    type Output = Vec3;
    fn mul(self, other: Vec3) -> Self {
        Vec3::new(self.x * other.x, self.y * other.y, self.z * other.z)
    }
}

impl Mul<f64> for Vec3 {
    type Output = Vec3;
    fn mul(self, other: f64) -> Self {
        Vec3::new(self.x * other, self.y * other, self.z * other)
    }
}

impl Mul<Vec3> for f64 {
    type Output = Vec3;
    fn mul(self, other: Vec3) -> Vec3 {
        other * self
    }
}

impl Div<f64> for Vec3 {
    type Output = Vec3;
    fn div(self, other: f64) -> Self {
        Vec3::new(self.x / other, self.y / other, self.z / other)
    }
}

impl MulAssign<f64> for Vec3 {
    fn mul_assign(&mut self, other: f64) {
        self.x *= other;
        self.y *= other;
        self.z *= other
    }
}

impl AddAssign for Vec3 {
    fn add_assign(&mut self, other: Vec3) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z
    }
}

impl Index<usize> for Vec3 {
    type Output = f64;

    fn index(&self, idx : usize) -> &Self::Output {
        match idx {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Vec3 out of bounds")
        }
    }
}

#[derive(Debug, Copy, Clone)]
pub struct Vec2 {
    pub x: f64, 
    pub y: f64, 
}


impl Vec2 {
    pub fn new(x: f64, y: f64) -> Self { 
        Vec2{x, y} 
    }

    pub fn zero() -> Self { 
        Vec2{x: 0., y: 0.} 
    }

    pub fn length_sqr(&self) -> f64 {
        self.x*self.x + self.y*self.y
    }

    pub fn length(&self) -> f64 {
        f64::sqrt(self.length_sqr())
    }

    pub fn normalize(self) -> Self {
        self / self.length()
    }

    pub fn abs(self) -> Self {
        Vec2::new(f64::abs(self.x), f64::abs(self.y))
    }

    pub fn dot(self, other : Vec2) -> f64 {
        self.x*other.x + self.y*other.y
    }

}

impl From<[f32; 2]> for Vec2 {
    fn from(array : [f32; 2]) -> Self {
        Vec2::new(array[0].into(), array[1].into()) 
    }
}

impl Sub for Vec2 {
    type Output = Vec2;
    fn sub(self, other: Vec2) -> Self {
        Vec2::new(self.x - other.x, self.y - other.y)
    }
}

impl Neg for Vec2 {
    type Output = Vec2;

    fn neg(self) -> Self {
        Vec2::new(-self.x, -self.y)
    }
}

impl Add for Vec2 {
    type Output = Vec2;
    fn add(self, other: Vec2) -> Self {
        Vec2::new(self.x + other.x, self.y + other.y)
    }
}

impl Mul for Vec2 {
    type Output = Vec2;
    fn mul(self, other: Vec2) -> Self {
        Vec2::new(self.x * other.x, self.y * other.y)
    }
}

impl Mul<f64> for Vec2 {
    type Output = Vec2;
    fn mul(self, other: f64) -> Self {
        Vec2::new(self.x * other, self.y * other)
    }
}

impl Mul<Vec2> for f64 {
    type Output = Vec2;
    fn mul(self, other: Vec2) -> Vec2 {
        other * self
    }
}

impl Div<f64> for Vec2 {
    type Output = Vec2;
    fn div(self, other: f64) -> Self {
        Vec2::new(self.x / other, self.y / other)
    }
}

impl MulAssign<f64> for Vec2 {
    fn mul_assign(&mut self, other: f64) {
        self.x *= other;
        self.y *= other;
    }
}

impl AddAssign for Vec2 {
    fn add_assign(&mut self, other: Vec2) {
        self.x += other.x;
        self.y += other.y;
    }
}

impl Index<usize> for Vec2 {
    type Output = f64;

    fn index(&self, idx : usize) -> &Self::Output {
        match idx {
            0 => &self.x,
            1 => &self.y,
            _ => panic!("Vec2 out of bounds")
        }
    }
}
