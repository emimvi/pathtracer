#![allow(dead_code)]
extern crate rand;
extern crate rayon;
extern crate tobj;
extern crate image;
mod vec3;
mod lights;
mod scene;
mod geometry;
mod obj;
mod texture;
use vec3::*;
use scene::*;
use geometry::*;
use std::f64;
use rayon::prelude::*;

use std::fs::File;
use std::io;
use std::io::{Write, BufWriter};

fn sample_cosine_weighted_hemisphere() -> Vec3
{
	let rnd = rand::random::<f64>();
	let cos_theta = f64::sqrt(rnd);
	let sin_theta = f64::sqrt(1.0 - rnd);
	let phi = 2.0*f64::consts::PI*rand::random::<f64>();

	// Calculate new direction as if the z-axis were the normal
	let spherical_direction = Vec3::new(sin_theta*f64::cos(phi), sin_theta*f64::sin(phi), cos_theta);

	spherical_direction
}

pub fn shade(surface : &Surface, incident_ray : &Ray, scene : &Scene) -> Vec3 {
    if incident_ray.trace_depth > 2 { 
        return Vec3::zero() 
    };
    match surface.material.illumination_model {
        Some(IlluminationModel::Lambertian)  => {
            let color = surface.material.get_diffuse(&surface.uv);
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
                        let new_direction = sample_cosine_weighted_hemisphere().rotate_to(surface.normal);
                        let mut trace_ray = incident_ray.clone();
                        trace_ray.origin = surface.position;
                        trace_ray.direction = new_direction;
                        trace_ray.trace_depth += 1;
                        let mirr_obj = scene.trace_closest(&mut trace_ray);
                        illumination += mirr_obj.map_or(scene.background, |next_surface| {
                            shade(&next_surface, &mut trace_ray, &scene)
                        }) * color / luminance
                    }
                    illumination/indirect_samples as f64
                }
            };

            direct_illumination //+ _indirect_illumination
        },
        Some(IlluminationModel::Mirror) => {
            //if incident_ray.trace_depth > 64 
            //{
            //    Vec3::zero();
            //}
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
        None => Vec3::new(0.8, 0., 1.),
        _ => {
            panic!(format!("Uknown illumination model: {:?}\n", surface.material.illumination_model))
        },
    }
}


pub fn render(x:f64, y:f64, pixel_delta_x : f64, pixel_delta_y : f64, scene : &Scene) -> Vec3 {
    //println!("{} {}", x, y);
    let eye = Vec3{ x : 250., 
                    y : 250., 
                    z : -500.};
    let direction = Vec3::new(x, y, 2.).normalize();
    let mut cam_ray = Ray::new(eye, direction);
    cam_ray.rx_origin = eye;
    cam_ray.ry_origin = eye;
    cam_ray.rx_direction = Vec3::new(x + pixel_delta_x, y, 2.).normalize();
    cam_ray.ry_direction = Vec3::new(x, y + pixel_delta_y, 2.).normalize();

    //unsafe {
    //    static mut COUNT : u32 = 0;
    //    if COUNT < 10 {
    //        println!("{:#?}", cam_ray);
    //        COUNT += 1;
    //    }
    //}

    let closest_object = scene.trace_closest(&mut cam_ray);
    closest_object.map_or(scene.background, |surface| {
        shade(&surface, &mut cam_ray, scene)
    })
}

fn main() {
    let h = 1024;
    let w = 1024;
    let pixel_size_x = 1./(w as f64);
    let pixel_size_y = 1./(h as f64);
    let samples_per_pixel = 2; 

    let scene = Scene::cornell_box();

    //Generate jitter subsamples in range (0,1), and scale it to the pixelsize.
    let jitters : Vec<(f64, f64)> = generate_jitter(samples_per_pixel)
            .iter().map(|&(dx,dy)| (dx*pixel_size_x, dy*pixel_size_y)).collect();

    let image : Vec<Vec3> = (0..h*w).into_par_iter().map(|i| { 
        let x = i % w;
        let y = i / h;
        let pixel_ndc_x = ( (x as f64) + 0.5 ) / w as f64;
        let pixel_ndc_y = ( (y as f64) + 0.5 ) / h as f64;
        let pixel_ss_x = 2. * pixel_ndc_x - 1.;
        let pixel_ss_y = 1. - 2. * pixel_ndc_y;

        let pixel_result = jitters.iter().fold(Vec3::zero(), 
                           |sum, &(dx,dy)| sum + render(pixel_ss_x + dx, pixel_ss_y + dy, pixel_size_x/samples_per_pixel as f64, pixel_size_y/samples_per_pixel as f64,  &scene));

        pixel_result/jitters.len() as f64
    }).collect();
    
    
    if let Err(e) = save_ppm(h, w, &image) {
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


fn clamp(x : f64) -> f64 {
    if x < 0. { 0. } 
    else if x > 1. { 1. }
    else { x }
}

fn save_ppm(height:usize, width:usize, image_data: &[Vec3]) -> Result<(), io::Error> {
    let to_ppm = |x:f64| (clamp(x)*255f64 + 0.5) as u8;

    let file = File::create("default.ppm")?;
    let mut buffer = BufWriter::new(file);

    write!(buffer, "P3\n{} {}\n255\n", width, height)?;
    for pixel in image_data.iter() {
        write!(buffer, "{} {} {}\n", to_ppm(pixel.x), to_ppm(pixel.y), to_ppm(pixel.z))?;
    }
    write!(buffer, "")?;
    Ok(())
}
