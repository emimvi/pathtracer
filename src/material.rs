use algebra::vec3::Vec3;
use f64;
use geometry::Ray;
use geometry::Surface;
use microfacet::*;
use mipmap::MipMap;
use std::fmt::Debug;

fn sample_cosine_weighted_hemisphere(normal: &Vec3) -> Vec3 {
    let rnd = ::rand::random::<f64>();
    let cos_theta = f64::sqrt(rnd);
    let sin_theta = f64::sqrt(1.0 - rnd);
    let phi = 2.0 * f64::consts::PI * ::rand::random::<f64>();

    // Calculate new direction as if the z-axis were the normal
    let spherical_direction = Vec3::new(
        sin_theta * f64::cos(phi),
        sin_theta * f64::sin(phi),
        cos_theta,
    );

    spherical_direction.rotate_to(normal)
}

pub struct Glass {
    refractive_index: f64,
}

pub type BrdfSample = (Vec3, Vec3, f64);

#[derive(Debug)]
pub struct Lambertian {
    pub color: Vec3,
    pub emission: Vec3,
}

impl _Material for Lambertian {
    fn brdf(&self, _: Vec3, _: Vec3, _: Vec3) -> Vec3 {
        self.color * f64::consts::FRAC_1_PI
    }

    fn sample_brdf(&self, ray: Ray, normal: Vec3) -> BrdfSample {
        let sampled_direction = sample_cosine_weighted_hemisphere(&normal);
        let mut brdf_sample: BrdfSample;
        let color = self.brdf(ray.direction, sampled_direction, normal);
        let pdf = self.pdf(sampled_direction, ray.direction, normal);
        (color, sampled_direction, pdf)
    }

    fn pdf(&self, out_direction: Vec3, _: Vec3, normal: Vec3) -> f64 {
        let angle = normal.dot(out_direction);
        angle * f64::consts::FRAC_1_PI
    }

    fn get_emission(&self) -> Vec3 {
        self.emission
    }
}

pub trait _Material: Sync + Debug + 'static {
    fn brdf(&self, wo: Vec3, wi: Vec3, normal: Vec3) -> Vec3;
    fn sample_brdf(&self, ray: Ray, normal: Vec3) -> BrdfSample;
    fn pdf(&self, in_direction: Vec3, out_direction: Vec3, normal: Vec3) -> f64;

    fn get_emission(&self) -> Vec3;
}

#[derive(Debug, Clone)]
pub struct Material {
    pub name: String,
    pub ambient: Vec3,
    pub diffuse: Vec3,
    pub specular: Vec3,
    pub illumination_model: IlluminationModel,
    pub texture: Option<MipMap>,
    pub roughness: f64,
    pub micro: Option<Glossy>,
    pub fresnel: Option<FresnelDielectric>,
}

impl Material {
    pub fn diffuse(color: Vec3) -> Material {
        Material {
            name: String::from("diffuse"),
            diffuse: color,
            illumination_model: IlluminationModel::Lambertian,
            ..Material::default()
        }
    }

    pub fn constant(color: Vec3) -> Material {
        Material {
            name: String::from("diffuse"),
            diffuse: color,
            illumination_model: IlluminationModel::Constant,
            ..Material::default()
        }
    }

    pub fn get_diffuse(&self, surface: &Surface) -> Vec3 {
        if let Some(tex) = &self.texture {
            let st = surface.uv;
            //compute texture differentials
            let duvdx = [surface.dudx, surface.dvdx];
            let duvdy = [surface.dudy, surface.dvdy];

            tex.sample_mipmap(st[0], st[1], duvdx, duvdy)
        } else {
            self.diffuse
        }
    }

    pub fn is_emissive(&self) -> bool {
        self.ambient[0] > 0. || self.ambient[1] > 0. || self.ambient[2] > 0.
    }

    pub fn glass() -> Material {
        let fresnel = FresnelDielectric::new(1.5);
        Material {
            name: String::from("glass"),
            illumination_model: IlluminationModel::Transparent,
            fresnel: Some(fresnel),
            ..Material::default()
        }
    }

    pub fn mirror() -> Material {
        Material {
            name: String::from("mirror"),
            illumination_model: IlluminationModel::Mirror,
            ..Material::default()
        }
    }

    pub fn micro(color: Vec3, roughness: f64, ior: f64) -> Material {
        let micro_facet = Glossy::new(roughness, ior, color);
        Material {
            name: String::from("glossy"),
            diffuse: color,
            illumination_model: IlluminationModel::Glossy,
            micro: Some(micro_facet),
            ..Material::default()
        }
    }
}
impl Default for Material {
    fn default() -> Material {
        Material {
            name: String::from("none"),
            diffuse: Vec3::zero(),
            specular: Vec3::zero(),
            ambient: Vec3::zero(),
            illumination_model: IlluminationModel::Constant,
            texture: None,
            roughness: 0.,
            micro: None,
            fresnel: None,
        }
    }
}

#[derive(Debug, Copy, Clone)]
pub enum IlluminationModel {
    Constant = 0,
    Lambertian = 1,
    Glossy = 2,
    Mirror = 3,
    Transparent = 4,
    Micro = 20,
}
impl From<u8> for IlluminationModel {
    fn from(int: u8) -> Self {
        match int {
            0 => IlluminationModel::Constant,
            1 => IlluminationModel::Lambertian,
            2 => IlluminationModel::Glossy,
            3 => IlluminationModel::Mirror,
            4 => IlluminationModel::Transparent,
            20 => IlluminationModel::Micro,
            _ => panic!(format!("Uknown illumination model: {}\n", int)),
        }
    }
}
