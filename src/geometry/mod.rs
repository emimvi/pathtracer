pub mod ray;
pub use self::ray::*;

pub mod surface;
pub use self::surface::{Surface, SurfaceHit};

pub mod triangle_mesh;
pub use self::triangle_mesh::*;

use algebra::vec3::*;
use material::*;
use std::f64;

use bvh::*;

pub trait Intersectable : Sync {
    fn intersect(&self, ray: &mut Ray) -> Option<Surface>;
}

#[derive(Debug, Clone)]
pub struct Sphere {
    pub center: Vec3,
    pub radius: f64,
}
#[derive(Debug, Clone)]
pub struct Plane {
    pub origin: Vec3,
    pub normal: Vec3,
}

impl Intersectable for Plane {
    fn intersect(&self, ray: &mut Ray) -> Option<Surface> {
        let denom = self.normal.dot(ray.direction);
        if f64::abs(denom) < 1e-6 {
            return None;
        }

        let t = (self.origin - ray.origin).dot(self.normal) / denom;

        if t > ray.t_max || t < ray.t_min {
            return None;
        }

        ray.t_max = t;
        let surface = Surface::new(ray.origin + ray.direction * t, self.normal);
        Some(surface)
    }
}

impl Intersectable for Sphere {
    //Solve sphere equation, substituting in ray parametrisation
    #[allow(clippy::many_single_char_names)]
    fn intersect(&self, ray: &mut Ray) -> Option<Surface> {
        let p = ray.origin - self.center;
        let b = ray.direction.dot(p);
        let c = p.dot(p) - self.radius * self.radius;
        if b * b < 0. {
            return None;
        }

        let d = f64::sqrt(b * b - c);
        let t1 = (-b - d);
        let t2 = (-b + d);

        let t = if t1 > ray.t_min && t1 < ray.t_max {
            t1
        } else if t2 > ray.t_min && t2 < ray.t_max {
            t2
        } else {
            return None;
        };
        ray.t_max = t;
        let hit_position = ray.origin + ray.direction * t;
        let normal = (hit_position - self.center).normalize();

        //Differentials
        let hit_local = hit_position - self.center;
        let phi = {
            let phi = f64::atan2(hit_local.y, hit_local.x);
            if phi < 0. {
                phi + 2. * f64::consts::PI
            } else {
                phi
            }
        };
        let mut inner = hit_local.z / self.radius;
        if inner < -1. || inner >= 1. {
            inner = 1.
        }
        let mut theta = f64::acos(inner);

        let zRadius = f64::sqrt(hit_local.x * hit_local.x + hit_local.y * hit_local.y);
        let invZRadius = 1. / zRadius;
        let cos_phi = hit_local.x * invZRadius;
        let sin_phi = hit_local.y * invZRadius;

        let dpdu = Vec3::new(
            -2. * f64::consts::PI * hit_local.y,
            2. * f64::consts::PI * hit_local.x,
            0.,
        );
        let dpdv = Vec3::new(
            hit_local.z * cos_phi,
            hit_local.z * sin_phi,
            -self.radius * f64::sin(theta),
        ) * -f64::consts::PI;

        // Compute sphere $\dndu$ and $\dndv$
        let d2Pduu =
            -2. * f64::consts::PI * 2. * f64::consts::PI * Vec3::new(hit_local.x, hit_local.y, 0.);
        let d2Pduv = (-f64::consts::PI)
            * hit_local.z
            * 2.
            * f64::consts::PI
            * Vec3::new(-sin_phi, cos_phi, 0.);
        let d2Pdvv = -(-f64::consts::PI)
            * (-f64::consts::PI)
            * Vec3::new(hit_local.x, hit_local.y, hit_local.z);

        // Compute coefficients for fundamental forms
        let E = Vec3::dot(dpdu, dpdu);
        let F = Vec3::dot(dpdu, dpdv);
        let G = Vec3::dot(dpdv, dpdv);
        let N = Vec3::normalize(Vec3::cross(dpdu, dpdv));
        let e = Vec3::dot(N, d2Pduu);
        let f = Vec3::dot(N, d2Pduv);
        let g = Vec3::dot(N, d2Pdvv);

        // Compute $\dndu$ and $\dndv$ from fundamental form coefficients
        let invEGF2 = 1. / (E * G - F * F);
        let dndu = ((f * F - e * G) * invEGF2 * dpdu + (e * F - f * E) * invEGF2 * dpdv);
        let dndv = ((g * F - f * G) * invEGF2 * dpdu + (f * F - g * E) * invEGF2 * dpdv);

        //println!("{} {} {} {}", hit_local, d2Pduu, d2Pduv, d2Pdvv );

        let mut surface = Surface {
            position: hit_position,
            normal,
            dpdu,
            dpdv,
            dndu,
            dndv,
            ..Surface::default()
        };
        surface.calculate_differentials(&ray);

        Some(surface)
    }
}

impl Boundable for Sphere {
    fn bounds(&self, _: f32, _: f32) -> BBox {
        BBox::span(
            self.center + Vec3::new(-self.radius, -self.radius, -self.radius),
            self.center + Vec3::new(self.radius, self.radius, self.radius),
        )
    }
}
