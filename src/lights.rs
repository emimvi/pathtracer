use algebra::vec3::*;
use geometry::*;
use obj::*;
use std::f64;
use Scene;

pub fn sample_sphere() -> Vec3 {
    let rnd = rand::random::<f64>();
    let rnd2 = rand::random::<f64>();

    let z = 1. - 2. * rnd;
    let r = f64::sqrt(f64::max(0., 1. - z * z));
    let phi = 2. * f64::consts::PI * rnd2;
    return Vec3::new(r * f64::cos(phi), r * f64::sin(phi), z);
}

//Direction towards the light, and light intensity, and pdf.
type LightSample = Option<(Vec3, Vec3, f64)>;

pub trait Light: Sync {
    fn sample(&self, surface: &Surface, scene: &'static Scene) -> LightSample;
}

#[derive(Debug, Clone)]
pub struct DirectionalLight {
    pub direction: Vec3,
    pub radiance: f64,
}
//All lights are assumed white.
#[derive(Debug, Clone)]
pub struct PointLight {
    pub position: Vec3,
    pub intensity: f64,
}

//All lights are assumed white.
pub struct AreaLight {
    pub mesh: &'static Sphere,
    pub intensity: Vec3
}

impl Light for AreaLight {
    fn sample(&self, surface: &Surface, scene: &'static Scene) -> LightSample {
        //Sample point on light
        let (light_pos, normal, pdf) = self.mesh.sample_surface();
        let diff = light_pos - surface.position;
        let light_direction = diff.normalize();
        if normal.dot(-light_direction) < 0. {
            return None;
        }
        let mut shadow_ray = Ray::new(surface.position, light_direction);
        shadow_ray.t_max = diff.length();

        // Check if surface is shadowed by another object
        if scene.trace_any(&mut shadow_ray).is_none() {
            let irradiance = self.intensity / diff.length_sqr();
            Some((light_direction, irradiance, pdf))
        } else {
            None
        }
    }
}

pub trait SurfaceSample {

    //point, normal, pdf
    fn sample_surface(&self) -> (Vec3, Vec3, f64);
}

impl SurfaceSample for Sphere {
    fn sample_surface(&self) -> (Vec3, Vec3, f64) {
        let unit_sample = sample_sphere();
        let normal = unit_sample;
        let surface_point = unit_sample * self.radius + self.center;

        (surface_point, normal, 0.25*f64::consts::FRAC_1_PI)
    }
}

//impl Light for PointLight {
//    fn sample(&self, surface: &Surface, scene: &Scene) -> LightSample {
//        let diff = self.position - surface.position;
//        let light_direction = diff.normalize();
//        let mut shadow_ray = Ray::new(surface.position, light_direction);
//        shadow_ray.t_max = diff.length();
//
//        // Check if surface is shadowed by another object
//        if scene.trace_any(&mut shadow_ray).is_none() {
//            let irradiance = self.intensity / diff.length_sqr();
//            Some((light_direction, irradiance))
//        } else {
//            None
//        }
//    }
//}
//
//impl Light for DirectionalLight {
//    fn sample(&self, surface: &Surface, scene: &Scene) -> LightSample {
//        let mut shadow_ray = Ray::new(surface.position, -self.direction);
//
//        // Check if surface is shadowed by another object
//        if scene.trace_any(&mut shadow_ray).is_none() {
//            Some((-self.direction, self.radiance))
//        } else {
//            None
//        }
//    }
//}

fn EstimateDirect(it : &SurfaceHit, /*uScattering, */
                       light : &dyn Light, /* uLight, */
                        scene : &'static Scene, /* Sampler &sampler, */
                        ray_direction : Vec3,
                        specular : bool) -> Vec3 {
    let bsdfFlags = specular;

    // Sample light source with multiple importance sampling
    let light_sample = light.sample(&it.surface, scene);
    if light_sample.is_none() {
        return Vec3::zero()
    }
    let (wi, light_intensity, lightPdf) = light_sample.unwrap();

    if lightPdf < 1e-10 {
        return Vec3::zero();
    }

    // Compute BSDF or phase function's value for light sample
    // Evaluate BSDF for light sampling strategy
    //
    let isect = it;
    let f : Vec3 = isect.material.brdf(Vec3::zero(), Vec3::zero(), Vec3::zero()) * f64::abs(wi.dot(isect.surface.normal));
    let scatteringPdf = isect.material.pdf(ray_direction, wi, isect.surface.normal);
    // Compute effect of visibility for light source sample
    let Ld = f * light_intensity / lightPdf;
    Ld
}
