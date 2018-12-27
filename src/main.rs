#![allow(dead_code)]
#![allow(non_snake_case)]
#![allow(unused)]
extern crate rand;
extern crate rayon;
extern crate tobj;
extern crate image;
extern crate bvh as extern_bvh;
#[macro_use]
extern crate clap;
#[macro_use]
extern crate lazy_static;
mod vec3;
mod lights;
mod scene;
mod geometry;
mod obj;
mod mipmap;
mod bvh;
mod material;
mod microfacet;
use microfacet::*;
use material::*;
use geometry::*;
use vec3::*;
use scene::*;
use std::f64;
use rayon::prelude::*;
use mipmap::Image;
use clap::App;

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

pub fn shade(surface : &Surface, incident_ray : &Ray, scene : &Scene, trace_depth : usize) -> Vec3 {
    if trace_depth > CONFIG.bounces { 
        return Vec3::zero() 
    };
    let trace_depth = trace_depth + 1;

    let material = if let Some(m) = surface.material.as_ref() {
        m
    } else {
        return Vec3::new(1.,0.7,0.85) //Pink!
    };
    match material.illumination_model {
        IlluminationModel::Constant => { 
            material.get_diffuse(&surface) 
        },
        IlluminationModel::Lambertian => {
            let color = material.get_diffuse(&surface);
            let f = color*f64::consts::FRAC_1_PI;
            let mut direct_illumination : Vec3 = scene.lights.iter().fold(Vec3::zero(), |prev, light| {
                if let Some((light_dir, radiance)) = light.sample(surface, scene)
                {
                    let angle = f64::max(0., light_dir.dot(surface.normal));
                    prev + f*angle*radiance
                } else {
                    prev
                }
            });
            direct_illumination += material.ambient*color;

            let _indirect_diffuse = {
                //Use luminance as the probability of a ray getting absorped. 
                let luminance = (color[0]+color[1]+color[2])/3.; //TODO Importance sample luminance
                if luminance < rand::random::<f64>() {
                    Vec3::new(0., 0., 0.)
                } else {
                    let mut illumination = Vec3::zero();
                    let sampled_direction = sample_cosine_weighted_hemisphere(&surface.normal);
                    let angle = surface.normal.dot(sampled_direction);
                    let pdf = angle*f64::consts::FRAC_1_PI;
                    if pdf <= 0. {
                        Vec3::zero()
                    } else {
                        //let mut trace_ray = Ray::new(surface.position, sampled_direction);
                        let mut trace_ray = reflect_diff(&surface, incident_ray, sampled_direction, pdf);
                        let mirr_obj = scene.trace_closest(&mut trace_ray);
                        let color_bleed = mirr_obj.map_or(scene.background, |next_surface| {
                            shade(&next_surface, &mut trace_ray, &scene, trace_depth)
                        });

                        illumination += color_bleed * f * angle/ pdf;
                        illumination / luminance
                    }
                }
            };
            direct_illumination + _indirect_diffuse
        },
        IlluminationModel::Transparent => {
            let fresnel = material.fresnel.as_ref().unwrap();
            let cos_theta = (-incident_ray.direction).dot(surface.normal);
            let eta = fresnel.eta(cos_theta);
            let fresnel_r= fresnel.dielectric(cos_theta);
            let normal = if cos_theta < 0. {
                -surface.normal
            } else {
                surface.normal
            };

            let reflected_dir = incident_ray.direction.reflect(normal);
            let mut trace_ray = reflect_diff(&surface, incident_ray, reflected_dir, 1.);
            let mirr_obj = scene.trace_closest(&mut trace_ray);
            let reflected_radiance = mirr_obj.map_or(scene.background, |next_surface| {
                shade(&next_surface, &mut trace_ray, &scene, trace_depth)
            });

            let refracted_dir = incident_ray.direction.refract(normal, eta);
            let mut refracted_ray = refract_diff(&surface, incident_ray, refracted_dir, 1.);
            let mirr_obj = scene.trace_closest(&mut refracted_ray);
            let refracted_radiance = mirr_obj.map_or(scene.background, |next_surface| {
                shade(&next_surface, &mut refracted_ray, &scene, trace_depth)
            });
            fresnel_r*reflected_radiance + (1.-fresnel_r)*refracted_radiance
        },
         IlluminationModel::Glossy | IlluminationModel::Micro  => {
            // let (color, _, pdf) = material.micro.as_ref().unwrap().shade(incident_ray, surface); 
            // if pdf == 0. { return Vec3::zero() }

            let mat = material.micro.as_ref().unwrap();
            let direct_illum = scene.lights.iter().fold(Vec3::zero(), |prev, light| {
                            if let Some((light_dir, radiance)) = light.sample(surface, scene)
                            {
                                let color = mat.f(-incident_ray.direction, light_dir, surface.normal); 
                                let angle = f64::max(0.,light_dir.dot(surface.normal));
                                let r = color*radiance*angle;
                                prev + r
                            } else {
                                prev
                            }
                        });

            let indirect_illum = {
                let (f, direction, pdf) = mat.shade(incident_ray, surface);
                let cos_theta = direction.dot(surface.normal);
                if pdf == 0. || cos_theta <= 0. {
                     Vec3::zero()
                } else {
                    let mut trace_ray = reflect_diff(&surface, incident_ray, direction, pdf);

                    let mirr_obj = scene.trace_closest(&mut trace_ray);
                    let refl_col = mirr_obj.map_or(scene.background, |next_surface| {
                      shade(&next_surface, &mut trace_ray, &scene, trace_depth)
                    });
                    refl_col * f * cos_theta / pdf
                }
            };

            direct_illum + indirect_illum
        },
        IlluminationModel::Mirror => {
            let reflected_dir = incident_ray.direction.reflect(surface.normal);
            let mut trace_ray = reflect_diff(&surface, incident_ray, reflected_dir, 1.);
            let mirr_obj = scene.trace_closest(&mut trace_ray);
            mirr_obj.map_or(scene.background, |next_surface| {
                shade(&next_surface, &mut trace_ray, &scene, trace_depth)
            })
        },
        _ => {
            panic!(format!("Uknown illumination model: {:?}\n", material.illumination_model))
        },
    }
}

fn refract_diff(isect : &Surface, inc_ray : &Ray, new_direction : Vec3, pdf : f64) -> Ray {
    let mut new_ray = Ray::new(isect.position, new_direction);
    if !CONFIG.use_diffs {
        return new_ray;
    }

    let wo = -inc_ray.direction;
    let wi = new_ray.direction;
    let cos_theta = wo.dot(isect.normal);
    let eta = isect.material.as_ref().unwrap_or_else(|| panic!("No material")).fresnel.as_ref().unwrap_or_else(|| panic!("No fresnel")).eta(cos_theta);
    let normal = if cos_theta < 0. {
        -isect.normal
    } else {
        isect.normal
    };

    new_ray.rx_origin = isect.position + isect.dpdx;
    new_ray.ry_origin = isect.position + isect.dpdy;
    let dndx = isect.dndu * isect.dudx + isect.dndv * isect.dvdx;
    let dndy = isect.dndu * isect.dudy + isect.dndv * isect.dvdy;
    let dwodx = -inc_ray.rx_direction - wo;
    let dwody = -inc_ray.ry_direction - wo;
    let dDNdx = Vec3::dot(dwodx, normal) + Vec3::dot(wo, dndx);
    let dDNdy = Vec3::dot(dwody, normal) + Vec3::dot(wo, dndy);
    let mu = eta * Vec3::dot(-wo, normal) - Vec3::dot(wi, normal);
    let dmudx = (eta - (eta * eta * Vec3::dot(-wo, normal)) / Vec3::dot(wi, normal)) * dDNdx;
    let dmudy = (eta - (eta * eta * Vec3::dot(-wo, normal)) / Vec3::dot(wi, normal)) * dDNdy;

    new_ray.rx_direction = wi + eta * dwodx - (mu * dndx + dmudx * normal);
    new_ray.ry_direction = wi + eta * dwody - (mu * dndy + dmudy * normal);
    new_ray
}


fn reflect_diff(isect : &Surface, inc_ray : &Ray, new_direction : Vec3, pdf : f64) -> Ray {
    let mut new_ray = Ray::new(isect.position, new_direction);
    if !CONFIG.use_diffs {
        return new_ray;
    }

    //Set differentials for new ray depending on material.
    use IlluminationModel::*;
    let material_type = isect.material.as_ref().unwrap().illumination_model;
    let x = match material_type {
        Constant | Lambertian | Glossy => { 
            let diff_size = CONFIG.diff_size;
            new_ray.rx_origin = isect.position + isect.dpdx;
            new_ray.ry_origin = isect.position + isect.dpdy;
            let (v0, v1) = new_direction.create_tangent_vectors();
            new_ray.rx_direction = (new_direction + diff_size*v0).normalize();
            new_ray.ry_direction = (new_direction + diff_size*v1).normalize();
        },
        Mirror | Transparent => { 
            let wo = -inc_ray.direction;
            let wi = new_ray.direction;
            let cos_theta = wo.dot(isect.normal);
            let normal = if cos_theta < 0. {
                -isect.normal
            } else {
                isect.normal
            };
            new_ray.rx_origin = isect.position + isect.dpdx;
            new_ray.ry_origin = isect.position + isect.dpdy;
            let dndx = isect.dndu * isect.dudx + isect.dndv * isect.dvdx;
            let dndy = isect.dndu * isect.dudy + isect.dndv * isect.dvdy;
            let dwodx = -inc_ray.rx_direction - wo;
            let dwody = -inc_ray.ry_direction - wo;
            let dDNdx = Vec3::dot(dwodx, normal) + Vec3::dot(wo, dndx);
            let dDNdy = Vec3::dot(dwody, normal) + Vec3::dot(wo, dndy);
            new_ray.rx_direction = wi - dwodx +
                2. * (Vec3::dot(wo, normal) * dndx + dDNdx * normal);
            new_ray.ry_direction = wi - dwody +
                2. * (Vec3::dot(wo, normal) * dndy + dDNdy * normal);
        },
        _ => { panic!("Unsupported illumination model") }
    };

    new_ray
}


pub fn render(x:f64, y:f64, pixel_delta_x : f64, pixel_delta_y : f64, scene : &Scene) -> Vec3 {
    let mut cam_ray = scene.camera.ray(x,y, pixel_delta_x, pixel_delta_y);

    let closest_object = scene.trace_closest(&mut cam_ray);
    closest_object.map_or(scene.background, |surface| {
        shade(&surface, &mut cam_ray, scene, 0)
    })
}



pub struct Config {
    dim : usize,
    samples_pr_pixel : usize,
    diff_size : f64,
    output : String,
    bounces : usize,
    use_diffs : bool
}


lazy_static! {
    static ref CONFIG: Config = {
        let matches = App::new("myapp")
            .args_from_usage("--dim [size]        'Sets dimension'
                             --samples [num]      'Samples pr. pixel'
                             --diff [float]       'Differential size'
                             -b [num]            'max Bounces'
                             -o [str]             'Output file'
                             --no-diffs            'Disable ray differentials'
                             ").get_matches();

        let dim = value_t!(matches, "dim", usize).unwrap_or(512);
        let samples_pr_pixel = value_t!(matches, "samples", usize).unwrap_or(1);
        let bounces = value_t!(matches, "b", usize).unwrap_or(5);
        let diff_size = value_t!(matches, "diff", f64).unwrap_or(0.0);
        let diff_string = diff_size.to_string().replace(".", "");
        let use_diffs = !matches.is_present("no-diffs");

        let output = matches.value_of("o").unwrap_or("res");
        let fname = format!("{}_{}_{}spp_df{}_td{}.ppm", 
                            output,
                            dim, 
                            samples_pr_pixel, 
                            if use_diffs {diff_string} else {"nodiffs".to_string()}, 
                            bounces);
        
        Config {
            dim,
            samples_pr_pixel,
            diff_size,
            output : String::from(fname),
            bounces,
            use_diffs
        }
    };
}

fn main() {
    let h = CONFIG.dim;
    let w = CONFIG.dim;
    let pixel_size_x = 1./(w as f64);
    let pixel_size_y = 1./(h as f64);
    let samples_per_pixel = CONFIG.samples_pr_pixel; //This is really n^2!!

    let scene = Scene::cornell_box();
    //let scene = Scene::glossy_planes();
    //let scene = Scene::quad();

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
                           |sum, &(dx,dy)| sum + render(pixel_ss_x + dx, pixel_ss_y + dy, pixel_size_x as f64, pixel_size_y as f64,  &scene));

        pixel_result/jitters.len() as f64
    }).collect();

    let image = Image {
        height : h,
        width: w,
        data : pixels
    };
    
    if let Err(e) = image.save_ppm(&CONFIG.output) {
        println!("{}", e);
    };

}

fn generate_jitter(num_samples : usize) -> Vec<(f64, f64)> {
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

