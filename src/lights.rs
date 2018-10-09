use vec3::*;
use geometry::*;
use std::f64;
use Scene;

//Direction towards the light, and light intensity.
type LightSample = Option<(Vec3, f64)>;


pub trait Light : Sync  {
    fn sample(&self, surface : &Surface, scene : &Scene) -> LightSample;
}

#[derive(Debug, Copy, Clone)]
pub struct DirectionalLight {
    pub direction : Vec3,
    pub radiance : f64
}
//All lights are assumed white.
#[derive(Debug, Copy, Clone)]
pub struct PointLight {
    pub position : Vec3,
    pub intensity : f64
}

impl Light for PointLight {
    fn sample(&self, surface : &Surface, scene : &Scene) -> LightSample {

        let diff = self.position - surface.position;
        let light_direction = diff.normalize();
        let mut shadow_ray = Ray::new(surface.position, light_direction);
        shadow_ray.t_max = diff.length();

        // Check if surface is shadowed by another object
        if scene.trace_any(&mut shadow_ray).is_none() 
        {
            let irradiance = self.intensity/(f64::consts::PI*4.*diff.length_sqr());
            Some( (light_direction, irradiance) )
        } else {
            None
        }
    }
}

impl Light for DirectionalLight {

    fn sample(&self, surface : &Surface, scene : &Scene) -> LightSample {
        let mut shadow_ray = Ray::new(surface.position, -self.direction);

        // Check if surface is shadowed by another object
        if scene.trace_any(&mut shadow_ray).is_none() 
        {
            Some( (-self.direction, self.radiance) )
        } else {
            None
        }
    }
}
