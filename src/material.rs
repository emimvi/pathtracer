use vec3::Vec3;
use mipmap::MipMap;
use geometry::Surface;

#[derive(Debug, Clone)]
pub struct Material {
    pub name : String,
    pub ambient: Vec3,
    pub diffuse: Vec3,
    pub specular: Vec3,
    pub illumination_model: IlluminationModel,
    pub texture : Option<MipMap>
}

impl Material {
    pub fn diffuse(color : Vec3) -> Material {
        Material {
            name : String::from("diffuse"),
            diffuse: color,
            illumination_model: IlluminationModel::Lambertian,
            ..Material::default()
        }
    }

    pub fn get_diffuse(&self, surface : &Surface) -> Vec3 {
        if let Some(tex) = &self.texture {
            let st = surface.uv;
            //compute texture differentials
            let duvdx = [surface.dudx, surface.dvdx];
            let duvdy = [surface.dudy, surface.dvdy];

            let width = f64::max(f64::max(f64::abs(duvdx[0]),
            f64::abs(duvdx[1])),
            f64::max(f64::abs(duvdy[0]),
            f64::abs(duvdy[1])));
            tex.sample_mipmap(st[0], st[1], width)
        } else {
            self.diffuse
        }
    }

    pub fn mirror() -> Material {
        Material {
            name : String::from("mirror"),
            illumination_model: IlluminationModel::Mirror,
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
            illumination_model: IlluminationModel::Constant,
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

