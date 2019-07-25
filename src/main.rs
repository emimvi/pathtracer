#![allow(dead_code)]
#![allow(non_snake_case)]
#![allow(unused)]
extern crate bvh as extern_bvh;
extern crate image;
extern crate minifb;
extern crate rand;
extern crate rayon;
extern crate tobj;
#[macro_use]
extern crate clap;
#[macro_use]
extern crate lazy_static;
mod algebra;
mod bvh;
mod geometry;
mod lights;
mod material;
mod microfacet;
mod mipmap;
mod obj;
mod scene;
use algebra::vec3::*;
use clap::App;
use geometry::*;
use material::*;
use microfacet::*;
use minifb::{Key, Window, WindowOptions};
use mipmap::Image;
use rayon::prelude::*;
use scene::*;
use std::f64;

fn sample_cosine_weighted_hemisphere(normal: &Vec3) -> Vec3 {
    let rnd = rand::random::<f64>();
    let cos_theta = f64::sqrt(rnd);
    let sin_theta = f64::sqrt(1.0 - rnd);
    let phi = 2.0 * f64::consts::PI * rand::random::<f64>();

    // Calculate new direction as if the z-axis were the normal
    let spherical_direction = Vec3::new(
        sin_theta * f64::cos(phi),
        sin_theta * f64::sin(phi),
        cos_theta,
    );

    spherical_direction.rotate_to(normal)
}
pub fn shade(surface: &Surface, incident_ray: &Ray, scene: &Scene, trace_depth: usize) -> Vec3 {
    Vec3::zero()
}

//pub fn shade(surface : &Surface, incident_ray : &Ray, scene : &Scene, trace_depth : usize) -> Vec3 {
//    if trace_depth > CONFIG.bounces {
//        return Vec3::zero()
//    };
//    let trace_depth = trace_depth + 1;
//
//    let material = if let Some(m) = surface.material.as_ref() {
//        m
//    } else {
//        return Vec3::new(1.,0.7,0.85) //Pink!
//    };
//    match material.illumination_model {
//        IlluminationModel::Constant => {
//            material.get_diffuse(&surface)
//        },
//        IlluminationModel::Lambertian => {
//            let color = material.get_diffuse(&surface);
//            let f = color*f64::consts::FRAC_1_PI;
//            let mut direct_illumination : Vec3 = scene.lights.iter().fold(Vec3::zero(), |prev, light| {
//                if let Some((light_dir, radiance)) = light.sample(surface, scene)
//                {
//                    let angle = f64::max(0., light_dir.dot(surface.normal));
//                    prev + f*angle*radiance
//                } else {
//                    prev
//                }
//            });
//            direct_illumination += material.ambient*color;
//
//
//            let _indirect_diffuse = {
//                //Use luminance as the probability of a ray getting absorped.
//                let luminance = if (trace_depth < 5) {
//                    1.
//                } else {
//                    (color[0]+color[1]+color[2])/3. //TODO Importance sample luminance
//                };
//                if luminance < rand::random::<f64>() {
//                    Vec3::new(0., 0., 0.)
//                } else {
//                    let mut illumination = Vec3::zero();
//                    let sampled_direction = sample_cosine_weighted_hemisphere(&surface.normal);
//                    let angle = surface.normal.dot(sampled_direction);
//                    let pdf = angle*f64::consts::FRAC_1_PI;
//                    if pdf <= 0. {
//                        Vec3::zero()
//                    } else {
//                        //let mut trace_ray = Ray::new(surface.position, sampled_direction);
//                        let mut trace_ray = reflect_diff(&surface, incident_ray, sampled_direction, pdf);
//                        let mirr_obj = scene.trace_closest(&mut trace_ray);
//                        let color_bleed = mirr_obj.map_or(scene.background, |next_surface| {
//                            shade(&next_surface, &mut trace_ray, &scene, trace_depth)
//                        });
//
//                        illumination += color_bleed * f * angle/ pdf;
//                        illumination / luminance
//                    }
//                }
//            };
//
//            direct_illumination + _indirect_diffuse
//        },
//        IlluminationModel::Transparent => {
//            let fresnel = material.fresnel.as_ref().unwrap();
//            let cos_theta = (-incident_ray.direction).dot(surface.normal);
//            let eta = fresnel.eta(cos_theta);
//            let fresnel_r= fresnel.dielectric(cos_theta);
//            let normal = if cos_theta < 0. {
//                -surface.normal
//            } else {
//                surface.normal
//            };
//
//            let reflected_dir = incident_ray.direction.reflect(normal);
//            let mut trace_ray = reflect_diff(&surface, incident_ray, reflected_dir, 1.);
//            let mirr_obj = scene.trace_closest(&mut trace_ray);
//            let reflected_radiance = mirr_obj.map_or(scene.background, |next_surface| {
//                shade(&next_surface, &mut trace_ray, &scene, trace_depth)
//            });
//
//            let refracted_dir = incident_ray.direction.refract(normal, eta);
//            let mut refracted_ray = refract_diff(&surface, incident_ray, refracted_dir, 1.);
//            let mirr_obj = scene.trace_closest(&mut refracted_ray);
//            let refracted_radiance = mirr_obj.map_or(scene.background, |next_surface| {
//                shade(&next_surface, &mut refracted_ray, &scene, trace_depth)
//            });
//            fresnel_r*reflected_radiance + (1.-fresnel_r)*refracted_radiance
//        },
//         IlluminationModel::Glossy  => {
//            // let (color, _, pdf) = material.micro.as_ref().unwrap().shade(incident_ray, surface);
//            // if pdf == 0. { return Vec3::zero() }
//
//            let mat = material.micro.as_ref().unwrap();
//            let direct_illum = scene.lights.iter().fold(Vec3::zero(), |prev, light| {
//                            if let Some((light_dir, radiance)) = light.sample(surface, scene)
//                            {
//                                let color = mat.f(-incident_ray.direction, light_dir, surface.normal);
//                                let angle = f64::max(0.,light_dir.dot(surface.normal));
//                                let r = color*radiance*angle;
//                                prev + r
//                            } else {
//                                prev
//                            }
//                        });
//
//            let indirect_illum = {
//                let (f, direction, pdf) = mat.shade(incident_ray, surface);
//                let cos_theta = direction.dot(surface.normal);
//                if pdf == 0. || cos_theta <= 0. {
//                     Vec3::zero()
//                } else {
//                    let mut trace_ray = reflect_diff(&surface, incident_ray, direction, pdf);
//
//                    let mirr_obj = scene.trace_closest(&mut trace_ray);
//                    let refl_col = mirr_obj.map_or(scene.background, |next_surface| {
//                      shade(&next_surface, &mut trace_ray, &scene, trace_depth)
//                    });
//                    refl_col * f * cos_theta / pdf
//                }
//            };
//
//            direct_illum + indirect_illum
//        },
//        IlluminationModel::Mirror => {
//            let reflected_dir = incident_ray.direction.reflect(surface.normal);
//            let mut trace_ray = reflect_diff(&surface, incident_ray, reflected_dir, 1.);
//            let mirr_obj = scene.trace_closest(&mut trace_ray);
//            mirr_obj.map_or(scene.background, |next_surface| {
//                shade(&next_surface, &mut trace_ray, &scene, trace_depth)
//            })
//        },
//        _ => {
//            panic!(format!("Uknown illumination model: {:?}\n", material.illumination_model))
//        },
//    }
//}

//fn refract_diff(surface : &Surface, inc_ray : &Ray, new_direction : Vec3, pdf : f64) -> Ray {
//    let mut new_ray = Ray::new(surface.position, new_direction);
//    if !CONFIG.use_diffs {
//        return new_ray;
//    }
//
//    let wo = -inc_ray.direction;
//    let wi = new_ray.direction;
//    let cos_theta = wo.dot(surface.normal);
//    let eta = surface.material.as_ref().unwrap_or_else(|| panic!("No material")).fresnel.as_ref().unwrap_or_else(|| panic!("No fresnel")).eta(cos_theta);
//    let normal = if cos_theta < 0. {
//        -surface.normal
//    } else {
//        surface.normal
//    };
//
//    new_ray.rx_origin = surface.position + surface.dpdx;
//    new_ray.ry_origin = surface.position + surface.dpdy;
//    let dndx = surface.dndu * surface.dudx + surface.dndv * surface.dvdx;
//    let dndy = surface.dndu * surface.dudy + surface.dndv * surface.dvdy;
//    let dwodx = -inc_ray.rx_direction - wo;
//    let dwody = -inc_ray.ry_direction - wo;
//    let dDNdx = Vec3::dot(dwodx, normal) + Vec3::dot(wo, dndx);
//    let dDNdy = Vec3::dot(dwody, normal) + Vec3::dot(wo, dndy);
//    let mu = eta * Vec3::dot(-wo, normal) - Vec3::dot(wi, normal);
//    let dmudx = (eta - (eta * eta * Vec3::dot(-wo, normal)) / Vec3::dot(wi, normal)) * dDNdx;
//    let dmudy = (eta - (eta * eta * Vec3::dot(-wo, normal)) / Vec3::dot(wi, normal)) * dDNdy;
//
//    new_ray.rx_direction = wi + eta * dwodx - (mu * dndx + dmudx * normal);
//    new_ray.ry_direction = wi + eta * dwody - (mu * dndy + dmudy * normal);
//    new_ray
//}
//
//
//fn reflect_diff(surface : &Surface, inc_ray : &Ray, new_direction : Vec3, pdf : f64) -> Ray {
//    let mut new_ray = Ray::new(surface.position, new_direction);
//    if !CONFIG.use_diffs {
//        return new_ray;
//    }
//    let wo = -inc_ray.direction;
//    let wi = new_ray.direction;
//
//    new_ray.rx_origin = surface.position + surface.dpdx;
//    new_ray.ry_origin = surface.position + surface.dpdy;
//    use IlluminationModel::*;
//    let material = surface.material.as_ref().unwrap();
//    let material_type = material.illumination_model;
//
//    //Set differentials for new ray depending on material.
//    match material_type {
//        Lambertian => {
//            let diff_size = CONFIG.diff_size;
//            let (v0, v1) = new_direction.create_tangent_vectors();
//            new_ray.rx_direction = (new_direction + diff_size*v0).normalize();
//            new_ray.ry_direction = (new_direction + diff_size*v1).normalize();
//        },
//        Glossy => {
//            let diff_size = if !CONFIG.roughdiff {
//                CONFIG.diff_size
//            } else {
//                CONFIG.diff_size*material.roughness / 0.8
//            };
//            let (v0, v1) = new_direction.create_tangent_vectors();
//            new_ray.rx_direction = (new_direction + diff_size*v0).normalize();
//            new_ray.ry_direction = (new_direction + diff_size*v1).normalize();
//        },
//        Mirror | Transparent => {
//            //Flip normal if we're inside an object.
//            let cos_theta = wo.dot(surface.normal);
//            let normal = if cos_theta < 0. {
//                -surface.normal
//            } else {
//                surface.normal
//            };
//            let dwodx = -inc_ray.rx_direction - wo;
//            let dwody = -inc_ray.ry_direction - wo;
//            let dndx = surface.dndu * surface.dudx + surface.dndv * surface.dvdx;
//            let dndy = surface.dndu * surface.dudy + surface.dndv * surface.dvdy;
//            let dDNdx = Vec3::dot(dwodx, normal) + Vec3::dot(wo, dndx);
//            let dDNdy = Vec3::dot(dwody, normal) + Vec3::dot(wo, dndy);
//            new_ray.rx_direction = wi - dwodx +
//                2. * (Vec3::dot(wo, normal) * dndx + dDNdx * normal);
//            new_ray.ry_direction = wi - dwody +
//                2. * (Vec3::dot(wo, normal) * dndy + dDNdy * normal);
//        },
//        _ => { panic!("Unsupported illumination model") }
//    };
//
//    new_ray
//}

pub fn render(x: f64, y: f64, pixel_delta_x: f64, pixel_delta_y: f64, scene: &Scene) -> Vec3 {
    let cam_ray = scene.camera.ray(x, y, pixel_delta_x, pixel_delta_y);

    let mut trace_ray = cam_ray;
    let mut color_result = Vec3::zero();
    let mut bounces = 0;

    let mut accumulated_pdf = Vec3::one();

    while bounces <= CONFIG.bounces {
        bounces += 1;

        let surface = if let Some(surf) = scene.trace_closest(&mut trace_ray) {
            surf
        } else {
            //Add background light and return
            color_result += accumulated_pdf * scene.background;
            return color_result;
        };

        let material = surface.material.unwrap();

        //Add emitted light
        color_result += accumulated_pdf * material.get_emission();

        //Sample brdf and update accumulated color
        let (color, new_direction, pdf) = material.sample_brdf(trace_ray, surface.normal);
        accumulated_pdf *= color * f64::abs(new_direction.dot(surface.normal)) / pdf;

        //Spawn new ray
        trace_ray = Ray::new(surface.position, new_direction);
    }

    color_result
}

pub struct Config {
    dim: usize,
    samples_pr_pixel: usize,
    diff_size: f64,
    output: String,
    bounces: usize,
    use_diffs: bool,
    roughdiff: bool,
}

lazy_static! {
    static ref CONFIG: Config = {
        let matches = App::new("myapp")
            .args_from_usage(
                "--dim [size]        'Sets dimension'
                             --samples [num]      'Samples pr. pixel'
                             --diff [float]       'Differential size'
                             -b [num]            'max Bounces'
                             -o [str]             'Output file'
                             -r                   'Use roughness based heuristic'
                             --no-diffs            'Disable ray differentials'
                             ",
            )
            .get_matches();

        let dim = value_t!(matches, "dim", usize).unwrap_or(512);
        let samples_pr_pixel = value_t!(matches, "samples", usize).unwrap_or(8);
        let bounces = value_t!(matches, "b", usize).unwrap_or(30);
        let diff_size = value_t!(matches, "diff", f64).unwrap_or(0.0);
        let diff_string = diff_size.to_string().replace(".", "");
        let use_diffs = !matches.is_present("no-diffs");
        let roughdiff = matches.is_present("p");
        let roughdiff_str = if roughdiff { "_p" } else { "" };

        let output = matches.value_of("o").unwrap_or("res");
        let fname = format!(
            "{}_{}_{}spp_df{}_td{}{}.ppm",
            output,
            dim,
            samples_pr_pixel,
            if use_diffs {
                diff_string
            } else {
                "nodiffs".to_string()
            },
            bounces,
            roughdiff_str
        );

        Config {
            dim,
            samples_pr_pixel,
            diff_size,
            output: fname,
            bounces,
            use_diffs,
            roughdiff,
        }
    };
}

fn main() {
    let h = CONFIG.dim;
    let w = CONFIG.dim;
    let pixel_size_x = 1. / (w as f64);
    let pixel_size_y = 1. / (h as f64);
    let samples_per_pixel = 1; //CONFIG.samples_pr_pixel; //This is really n^2!!

    let scene = Scene::cornell_box();
    //let scene = Scene::glossy_planes();
    //let scene = Scene::quad();

    let mut window = Window::new("Test - ESC to exit", w, h, WindowOptions::default())
        .unwrap_or_else(|e| {
            panic!("{}", e);
        });

    let mut pixels: Vec<Vec3> = vec![Vec3::zero(); h * w];
    let mut iteration = 1;
    while window.is_open() && !window.is_key_down(Key::Escape) {
        //Generate jitter subsamples in range (0,1), and scale it to the pixelsize.
        let jitters: Vec<(f64, f64)> = generate_jitter(samples_per_pixel)
            .iter()
            .map(|&(dx, dy)| (dx * pixel_size_x, dy * pixel_size_y))
            .collect();

        let next: Vec<Vec3> = (0..h * w)
            .into_par_iter()
            .map(|i| {
                let x = i % w; // x = 0..width
                let y = i / h; // y = 0..height
                let pixel_ndc_x = ((x as f64) + 0.5) / w as f64; //0..1, starting at top left.
                let pixel_ndc_y = ((y as f64) + 0.5) / h as f64; //0..1, starting at top left.
                let pixel_ss_x = 2. * pixel_ndc_x - 1.; //  -1..1, starting at bottom left.
                let pixel_ss_y = 1. - 2. * pixel_ndc_y; //  -1..1, starting at bottom left.

                let pixel_result = jitters.iter().fold(Vec3::zero(), |sum, &(dx, dy)| {
                    sum + render(
                        pixel_ss_x + dx,
                        pixel_ss_y + dy,
                        2. * pixel_size_x,
                        2. * pixel_size_y,
                        &scene,
                    )
                });

                pixel_result / jitters.len() as f64
            })
            .collect();

        for i in 0..pixels.len() {
            pixels[i] += next[i];
        }
        // let image = Image {
        //     height : h,
        //     width: w,
        //     data : pixels
        // };

        // if let Err(e) = image.save_ppm(&CONFIG.output) {
        //     println!("{}", e);
        // };

        // We unwrap here as we want this code to exit if it fails. Real applications may want to handle this in a different way
        window
            .update_with_buffer(&get_raw_pixels_u32(&pixels, iteration))
            .unwrap();
        iteration += 1;
    }
}

pub fn get_raw_pixels_u32(pixels: &[Vec3], iteration: u32) -> Vec<u32> {
    let denom = 1. / f64::from(iteration);
    pixels
        .iter()
        .map(|vec3| {
            let r = u32::from(Image::f64_to_8bit_color(vec3.x * denom));
            let g = u32::from(Image::f64_to_8bit_color(vec3.y * denom));
            let b = u32::from(Image::f64_to_8bit_color(vec3.z * denom));
            (r << 16 | g << 8 | b)
        })
        .collect()
}

fn generate_jitter(num_samples: usize) -> Vec<(f64, f64)> {
    let step = 1. / num_samples as f64;
    let mut jitters = Vec::new();
    for i in 0..num_samples {
        for j in 0..num_samples {
            jitters.push((
                step * (rand::random::<f64>() + i as f64) - 1. / 2.,
                step * (rand::random::<f64>() + j as f64) - 1. / 2.,
            ));
        }
    }
    jitters
}
