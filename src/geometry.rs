
use std::f64;
use vec3::*;
use texture::*;

#[derive(Debug, Copy, Clone)]
pub struct Ray {pub origin : Vec3,
                pub direction : Vec3,
                pub t_min : f64,
                pub t_max : f64,
                pub trace_depth : u8,
                pub rx_origin: Vec3,
                pub ry_origin: Vec3,
                pub rx_direction:  Vec3,
                pub ry_direction:  Vec3,
                }

impl Ray {
    pub fn new(origin : Vec3, direction : Vec3) -> Ray {
        Ray {
            origin : origin,
            direction : direction,
            t_min : 1e-8,
            t_max : f64::MAX,
            trace_depth : 0,
            rx_origin: Vec3::zero(),
            ry_origin: Vec3::zero(),
            rx_direction:  Vec3::zero(),
            ry_direction:  Vec3::zero(),
        }
    }

    pub fn scale_differentials(&mut self, scale_factor : f64) {
            self.rx_origin    =  self.origin + (self.rx_origin - self.origin) * scale_factor; 
            self.ry_origin    =  self.origin + (self.ry_origin - self.origin) * scale_factor; 
            self.rx_direction =  self.direction + (self.rx_direction - self.direction) * scale_factor;
            self.ry_direction =  self.direction + (self.ry_direction - self.direction) * scale_factor;
    }
}

#[derive(Debug, Clone)]
pub struct Surface<'a> { pub position : Vec3,
                     pub normal: Vec3,
                     pub material :  &'a Material,
                     pub uv : (f64, f64) }

#[derive(Debug, Clone)]
pub struct Sphere {pub center : Vec3, pub radius : f64, pub material : Material  }
#[derive(Debug, Clone)]
pub struct Plane  {pub origin : Vec3, pub normal : Vec3, pub material : Material }

#[derive(Debug, Clone)]
pub struct Material {
    pub name : String,
    pub ambient: Vec3,
    pub diffuse: Vec3,
    pub specular: Vec3,
    pub illumination_model: Option<IlluminationModel>,
    pub texture : Option<Texture>
}

impl Material {
    pub fn create_diffuse(color : Vec3) -> Material {
        Material {
            name : String::from("diffuse"),
            diffuse: color,
            illumination_model: Some(IlluminationModel::Lambertian),
            ..Material::default()
        }
    }

    pub fn get_diffuse(&self, uv : &(f64, f64)) -> Vec3 {
        if let Some(tex) = &self.texture {
            return tex.sample_nearest(uv.0, uv.1)
        }
        self.diffuse
    }

    pub fn create_mirror() -> Material {
        Material {
            name : String::from("mirror"),
            illumination_model: Some(IlluminationModel::Mirror),
            ..Material::default()
        }
    }
}
impl Default for Material {
    fn default() -> Material {
        Material {
            name : String::from("none"),
            diffuse: Vec3::zero(),
            specular: Vec3::zero(),
            ambient: Vec3::zero(),
            illumination_model: None,
            texture : None,
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
}
impl From<u8> for IlluminationModel {
    fn from(int : u8) -> Self {
        match int {
            0 => IlluminationModel::Constant,
            1 => IlluminationModel::Lambertian,
            2 => IlluminationModel::Glossy,
            3 => IlluminationModel::Mirror,
            4 => IlluminationModel::Transparent,
            _ => panic!(format!("Uknown illumination model: {}\n", int))
        }
    }
}

pub trait Intersectable : Sync {
    fn intersect(&self, ray : &mut Ray) -> Option<Surface>;
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
        let surface = Surface { position: ray.origin + ray.direction*t, 
                                normal : self.normal, 
                                material : &self.material,
                                uv : (0., 0.)
        } ;
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
        Some(Surface{position : hit_position,
                     normal: normal,
                     material : &self.material,
                     uv : (0., 0.)
        })
    }
}
