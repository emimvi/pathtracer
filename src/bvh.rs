use extern_bvh::{ray::RayT, vector::Vector};

pub use extern_bvh::{BBox, Boundable, BVH as e_BVH};

use algebra::vec3::*;
use geometry::{Intersectable, Ray, Surface, SurfaceHit};
use material::*;
use std::sync::Arc;

pub struct BVH
{
    geometry: e_BVH<Geometry>,
}

impl BVH
{
    pub fn new(geometry: Vec<Geometry>) -> Self {
        BVH {
            geometry: e_BVH::unanimated(8, geometry),
        }
    }

    pub fn empty() -> Self {
        BVH {
            geometry: e_BVH::empty(),
        }
    }

    pub fn intersect(&self, ray : &mut Ray) -> Option<SurfaceHit> {
        self.geometry.intersect(ray, |r, i| i.intersect_geometry(r))
    }
}

impl Boundable for BVH
{
    #[inline]
    fn bounds(&self, a: f32, b: f32) -> BBox {
        self.geometry.bounds(a, b)
    }
}


impl RayT for Ray {
    #[inline]
    fn direction(&self) -> Vector {
        self.direction.into()
    }

    #[inline]
    fn origin(&self) -> Vector {
        self.origin.into()
    }

    #[inline]
    fn t_max(&self) -> f32 {
        self.t_max as f32
    }

    #[inline]
    fn t_min(&self) -> f32 {
        self.t_min as f32
    }
}

pub trait Shape: Intersectable + Boundable  {}
impl<T> Shape for T where T: Intersectable + Boundable  {}

pub struct Geometry {
    pub shape: Box<dyn Shape>,
    pub material: Box<dyn _Material>,
}

impl Geometry {
    pub fn new(obj: Box<dyn Shape>, material: Box<dyn _Material>) -> Geometry {
        Geometry {
            shape: obj,
            material,
        }
    }

    #[inline]
    fn intersect_geometry(&self, mut ray: &mut Ray) -> Option<SurfaceHit> {
        if let Some(surface) = self.shape.intersect(&mut ray) {
            Some(SurfaceHit { material : self.material.as_ref(), surface })
        } else {
            None
        }
    }
}

impl Boundable for Geometry {
    #[inline]
    fn bounds(&self, a: f32, b: f32) -> BBox {
        self.shape.bounds(a, b)
    }
}

impl From<Vec3> for Vector {
    fn from(vec: Vec3) -> Self {
        Self::new(vec.x as f32, vec.y as f32, vec.z as f32)
    }
}
