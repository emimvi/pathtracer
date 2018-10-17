pub use extern_bvh::{ray::RayT, Boundable, BBox, vector::Vector, BVH as extern_BVH};

use geometry::{Intersectable, Surface, Ray};
use vec3::*;

impl<'a> RayT for &'a mut Ray {
    #[inline]
    fn direction(&self) -> Vector { self.direction.into() } 

    #[inline]
    fn origin(&self) -> Vector { self.origin.into() }

    #[inline]
    fn t_max(&self) -> f32 { self.t_max as f32 }

    #[inline]
    fn t_min(&self) -> f32 { self.t_min as f32}
}

pub trait Object : Intersectable + Boundable {}
impl<T> Object for T where T: Intersectable + Boundable {}

pub struct Geometry {
    pub object : Box<dyn Object>
}

impl Geometry {
    pub fn new(obj : impl Object + 'static) -> Geometry {
        Geometry {
            object: Box::new(obj)
        }
    }
}

impl Intersectable for Geometry {
    #[inline]
    fn intersect(&self, ray : &mut Ray) -> Option<Surface> {
        self.object.intersect(ray)
    }
}

impl Boundable for Geometry {
    #[inline]
    fn bounds(&self, a : f32, b : f32) -> BBox {
        self.object.bounds(a, b)
    }
}

pub type BVH<T> = extern_BVH<T>;

impl From<Vec3> for Vector {
    fn from(vec: Vec3) -> Self {
        Self::new(vec.x as f32, vec.y as f32, vec.z as f32)
    }
}
