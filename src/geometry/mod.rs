pub mod ray;
pub use self::ray::*;

pub mod surface;
pub use self::surface::Surface;

use material::*;
use std::sync::*;
use std::f64;
use vec3::*;

pub trait Intersectable : Sync {
    fn intersect(&self, ray : &mut Ray) -> Option<Surface>;
}

#[derive(Debug, Clone)]
pub struct Sphere {
    pub center : Vec3, 
    pub radius : f64, 
    pub material : Arc<Material>  
}
#[derive(Debug, Clone)]
pub struct Plane  {
    pub origin : Vec3, 
    pub normal : Vec3, 
    pub material : Arc<Material> 
}

impl Intersectable for Plane {

    fn intersect(&self, ray : &mut Ray) -> Option<Surface> {

        let denom = self.normal.dot(ray.direction);
        if f64::abs(denom) < 1e-6 {
            return None
        }

        let t = (self.origin - ray.origin).dot(self.normal) / denom;

        if t > ray.t_max || t < ray.t_min { return None }

        ray.t_max = t;
        let surface = Surface::new(ray.origin + ray.direction*t, 
                                self.normal, 
                                Arc::clone(&self.material));
        Some(surface)
    }
}

impl Intersectable for Sphere {

    //Solve sphere equation, substituting in ray parametrisation
    fn intersect(&self, ray : &mut Ray) -> Option<Surface> {
        let p = ray.origin - self.center;
        let b = 2. * ray.direction.dot(p);
        let c = p.dot(p) - self.radius*self.radius;
        let d = b*b - 4.*c;
        if d < 0. { return None }

        let t1 = (-b - f64::sqrt(d))/2.;
        let t2 = (-b + f64::sqrt(d))/2.;

        let t = if t1 > ray.t_min  && t1 < ray.t_max {
            t1
        } else if t2 > ray.t_min && t2 < ray.t_max {
            t2
        } else {
            return None
        };
        ray.t_max = t;
        let hit_position = ray.origin + ray.direction*t;
        let normal = (hit_position - self.center).normalize();
        Some(Surface::new(hit_position, normal, Arc::clone(&self.material)))
    }
}

