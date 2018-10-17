#![allow(unused_imports)]
extern crate rand;
extern crate rayon;
extern crate tobj;
extern crate image;
extern crate bvh as extern_bvh;
mod vec3;
mod lights;
mod scene;
mod geometry;
mod obj;
mod mipmap;
mod bvh;
use geometry::*;
use vec3::*;
use scene::*;
use std::f64;
use rayon::prelude::*;
use mipmap::Image;

fn sample_cosine_weighted_hemisphere(normal : &Vec3) -> Vec3
{
    let rnd = rand::random::<f64>();
    let cos_theta = f64::sqrt(rnd);
    let sin_theta = f64::sqrt(1.0 - rnd);
    let phi = 2.0*f64::consts::PI*rand::random::<f64>();

    // Calculate new direction as if the z-axis were the normal
    let spherical_direction = Vec3::new(sin_theta*f64::cos(phi), sin_theta*f64::sin(phi), cos_theta);

    spherical_direction.rotate_to(normal)
}

pub fn shade(surface : &Surface, incident_ray : &Ray, scene : &Scene) -> Vec3 {
    if incident_ray.trace_depth > 1 { 
        return Vec3::zero() 
    };
    match surface.material.illumination_model {
        m @ IlluminationModel::Lambertian | m @ IlluminationModel::Glossy  => {
            let color = surface.material.get_diffuse(&surface);
            let direct_illumination = if let Some((light_dir, radiance)) = scene.lights.sample(surface, scene)
            {
                let angle = f64::max(0., light_dir.dot(surface.normal));
                color*angle*radiance/f64::consts::PI
            } else {
                Vec3::new(0., 0., 0.)
            };

            let _indirect_illumination = {
                //Use luminance as the probability of a ray getting absorped. 
                let luminance = (color[0]+color[1]+color[2])/3.; //TODO Importance sample luminance
                if luminance < rand::random::<f64>() {
                    Vec3::new(0., 0., 0.)
                } else {
                    let indirect_samples = 1;
                    let mut illumination = Vec3::zero();
                    for _ in 0..indirect_samples {
                        let incident_direction = sample_cosine_weighted_hemisphere(&surface.normal);
                        let mut trace_ray = incident_ray.clone();
                        trace_ray.origin = surface.position;
                        trace_ray.direction = incident_direction;
                        trace_ray.trace_depth += 1;
                        let use_differentials = true;
                        if use_differentials {
                            trace_ray.rx_origin = surface.position + surface.dpdx;
                            trace_ray.ry_origin = surface.position + surface.dpdy;
                            let (v0, v1) = incident_direction.create_tangent_vectors();
                            trace_ray.rx_direction = (incident_direction + 0.2*v0).normalize();
                            trace_ray.ry_direction = (incident_direction + 0.2*v1).normalize();
                        } else {
                            trace_ray.rx_origin = Vec3::zero();
                            trace_ray.ry_origin = Vec3::zero();
                            trace_ray.rx_direction = Vec3::zero();
                            trace_ray.ry_direction = Vec3::zero();
                        }

                        let mirr_obj = scene.trace_closest(&mut trace_ray);
                        let color_bleed = mirr_obj.map_or(scene.background, |next_surface| {
                            shade(&next_surface, &mut trace_ray, &scene)
                        });


                        illumination += color_bleed * color / luminance
                    }
                    illumination/indirect_samples as f64
                }
            };

            match m {
                IlluminationModel::Lambertian => direct_illumination,
                IlluminationModel::Glossy  =>  _indirect_illumination,
                _ => unimplemented!{}
            }
        },
        IlluminationModel::Mirror => {
            let reflected_dir = (incident_ray.direction - 2.*incident_ray.direction.dot(surface.normal)*surface.normal).normalize();
            let mut trace_ray = incident_ray.clone();
            trace_ray.origin = surface.position;
            trace_ray.direction = reflected_dir;
            trace_ray.trace_depth += 1;

            let mirr_obj = scene.trace_closest(&mut trace_ray);
            mirr_obj.map_or(scene.background, |next_surface| {
                shade(&next_surface, &mut trace_ray, &scene)
            })
        },
        _ => {
            panic!(format!("Uknown illumination model: {:?}\n", surface.material.illumination_model))
        },
    }
}


pub fn render(x:f64, y:f64, pixel_delta_x : f64, pixel_delta_y : f64, scene : &Scene) -> Vec3 {
    let mut cam_ray = scene.camera.ray(x,y, pixel_delta_x, pixel_delta_y);

    let closest_object = scene.trace_closest(&mut cam_ray);
    closest_object.map_or(scene.background, |surface| {
        shade(&surface, &mut cam_ray, scene)
    })
}

fn main() {
    let h = 512;
    let w = 512;
    let pixel_size_x = 1./(w as f64);
    let pixel_size_y = 1./(h as f64);
    let samples_per_pixel = 1;

    let scene = Scene::glossy_planes();

    //Generate jitter subsamples in range (0,1), and scale it to the pixelsize.
    let jitters : Vec<(f64, f64)> = generate_jitter(samples_per_pixel)
            .iter().map(|&(dx,dy)| (dx*pixel_size_x, dy*pixel_size_y)).collect();

    let pixels : Vec<Vec3> = (0..h*w).into_par_iter().map(|i| { 
        let x = i % w; // x = 0..width
        let y = i / h; // y = 0..height
        let pixel_ndc_x = ( (x as f64) + 0.5 ) / w as f64; //0..1, starting at top left.
        let pixel_ndc_y = ( (y as f64) + 0.5 ) / h as f64; //0..1, starting at top left.
        let pixel_ss_x = 2. * pixel_ndc_x - 1.; //  -1..1, starting at bottom left.
        let pixel_ss_y = 1. - 2. * pixel_ndc_y; //  -1..1, starting at bottom left.

        let pixel_result = jitters.iter().fold(Vec3::zero(), 
                           |sum, &(dx,dy)| sum + render(pixel_ss_x + dx, pixel_ss_y + dy, 2.*pixel_size_x/samples_per_pixel as f64, 2.*pixel_size_y/samples_per_pixel as f64,  &scene));

        pixel_result/jitters.len() as f64
    }).collect();

    let image = Image {
        height : h,
        width: w,
        data : pixels
    };
    
    if let Err(e) = image.save_ppm("default.ppm") {
        println!("{}", e);
    };

}

fn generate_jitter(num_samples : u32) -> Vec<(f64, f64)> {
    let step = 1./num_samples as f64;
    let mut jitters = Vec::new();
    for i in 0..num_samples {
        for j in 0..num_samples {
            jitters.push( ( step * (rand::random::<f64>() + i as f64) - 1./2., 
                            step * (rand::random::<f64>() + j as f64) - 1./2.) );
        }
    }
    jitters
}
