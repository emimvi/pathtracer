
use extern_bvh::{ray::RayT, vector::Vector};

pub use extern_bvh::{Boundable, BBox, BVH as e_BVH};

use geometry::{Intersectable, Surface, Ray};
use vec3::*;
use material::*;
use std::sync::Arc;


pub struct BVH<T> where T : Shape {
    geometry : e_BVH<T>
}

impl<T> BVH<T> where T : Shape {
    pub fn new(geometry : Vec<T>) -> Self {
        BVH {
            geometry : e_BVH::unanimated(8, geometry)
        }
    }

    pub fn empty() -> Self {
        BVH {
            geometry : e_BVH::empty()
        }
    }
}

impl<T> Boundable for BVH<T> where T : Shape {
    #[inline]
    fn bounds(&self, a : f32, b : f32) -> BBox {
        self.geometry.bounds(a, b)
    }
}


impl<T> Intersectable for BVH<T> where T : Shape {
    #[inline]
    fn intersect(&self, mut ray : &mut Ray) -> Option<Surface> {
        self.geometry.intersect(&mut ray, |r, i| i.intersect(r))
    }
}

impl<'a> RayT for &'a mut Ray {
    #[inline]
    fn direction(&self) -> Vector { self.direction.into() } 

    #[inline]
    fn origin(&self) -> Vector { self.origin.into() }

    #[inline]
    fn t_max(&self) -> f32 { self.t_max as f32 }

    #[inline]
    fn t_min(&self) -> f32 { self.t_min as f32 }
}

pub trait Shape : Intersectable + Boundable + 'static {}
impl<T> Shape for T where T: Intersectable + Boundable + 'static {}

pub struct Geometry {
    pub shape : Box<dyn Shape>,
    pub material : Box<dyn _Material>
}

impl Geometry {
    pub fn new(obj : impl Shape, material : Box<dyn _Material>) -> Geometry {
        Geometry {
            shape: Box::new(obj),
            material : material
        }
    }
}

impl Intersectable for Geometry {
    #[inline]
    fn intersect(&self, ray : &mut Ray) -> Option<Surface> {
        if let Some(mut surface) = self.shape.intersect(ray) {
            surface.material = Some(&*self.material);
            return Some(surface)
        } 
        None
    }
}

impl Boundable for Geometry {
    #[inline]
    fn bounds(&self, a : f32, b : f32) -> BBox {
        self.shape.bounds(a, b)
    }
}

impl From<Vec3> for Vector {
    fn from(vec: Vec3) -> Self {
        Self::new(vec.x as f32, vec.y as f32, vec.z as f32)
    }
}
