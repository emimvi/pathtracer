
use std::f64;
use vec3::*;
use mipmap::*;

#[derive(Debug, Copy, Clone)]
pub struct Ray {pub origin : Vec3,
                pub direction : Vec3,
                pub t_min : f64,
                pub t_max : f64,
                pub trace_depth : u8,
                pub rx_origin: Vec3, //Differentials
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

    //pub fn scale_differentials(&mut self, scale_factor : f64) {
    //        self.rx_origin    =  self.origin + (self.rx_origin - self.origin) * scale_factor; 
    //        self.ry_origin    =  self.origin + (self.ry_origin - self.origin) * scale_factor; 
    //        self.rx_direction =  self.direction + (self.rx_direction - self.direction) * scale_factor;
    //        self.ry_direction =  self.direction + (self.ry_direction - self.direction) * scale_factor;
    //}
}

#[derive(Debug, Clone)]
pub struct Surface<'a> { pub position : Vec3,
                         pub normal: Vec3,
                         pub material :  &'a Material,
                         pub uv : Vec2,
                         pub dpdu : Vec3,
                         pub dpdv : Vec3,
                         pub dpdx : Vec3,
                         pub dpdy : Vec3,
                         pub dudx : f64,
                         pub dudy : f64,
                         pub dvdx : f64,
                         pub dvdy : f64,
                        }
impl<'a> Surface<'a> {
    pub fn new(position : Vec3, normal : Vec3, material : &'a Material) -> Surface<'a> {
        Surface { position, 
                  normal, 
                  material,
                  uv : Vec2::from([0., 0.]),
                  dpdu : Vec3::zero(),
                  dpdv : Vec3::zero(),
                  dpdx : Vec3::zero(),
                  dpdy : Vec3::zero(),
                  dudx : 0.,
                  dudy : 0.,
                  dvdx : 0.,
                  dvdy : 0.
        }
    }

    pub fn print_diff(&self, ray : &Ray) {
       if self.position.x > 498. && self.position.x < 500. &&
           self.position.y > 498. && self.position.y < 500. {
            println!("Ray: {:#?}", ray);
            println!("Surface:" );
            println!("position: {:?}", self.position); 
            println!("normal  : {:?}", self.normal  ); 
            println!("uv      : {:?}", self.uv      );
            println!("dpdu    : {:?}", self.dpdu    ) ;
            println!("dpdv    : {:?}", self.dpdv    ) ;
            println!("dpdx    : {:?}", self.dpdx    ) ;
            println!("dpdy    : {:?}", self.dpdy    ) ;
            println!("dudx    : {}", self.dudx    ) ;
            println!("dudy    : {}", self.dudy    ) ;
            println!("dvdx    : {}", self.dvdx    ) ;
            println!("dvdy    : {}", self.dvdy    ) ;
        }
    }

   pub fn calculate_differentials(&mut self, ray : &Ray) {
        //pbrt. Chapter 10, p601
        let d = self.normal.dot(self.position);

        //Distance of the dx ray to the surface
        let tx = -(self.normal.dot(ray.rx_origin) - d) / self.normal.dot(ray.rx_direction);
        let px = ray.rx_origin + tx * ray.rx_direction;

        let ty = -(self.normal.dot(ray.ry_origin) - d) / self.normal.dot(ray.ry_direction);
        let py = ray.ry_origin + ty * ray.ry_direction;

        self.dpdx = px - self.position;
        self.dpdy = py - self.position;

        let dim = if f64::abs(self.normal.x) > f64::abs(self.normal.y) 
                  && f64::abs(self.normal.x) > f64::abs(self.normal.z) {
            (1,2)
        } else if f64::abs(self.normal.y) > f64::abs(self.normal.z) {
            (0,2)
        } else {
            (0,1)
        };

        let mat_a = [ [ self.dpdu[dim.0], self.dpdv[dim.0] ],
                  [ self.dpdu[dim.1], self.dpdv[dim.1] ] ];
        let mat_bx = [ px[dim.0] - self.position[dim.0], px[dim.1] - self.position[dim.1] ];
        let mat_by = [ py[dim.0] - self.position[dim.0], py[dim.1] - self.position[dim.1] ]; 
        let (dudx, dvdx) = solve_linear_system_2x2(mat_a, mat_bx);
        self.dudx = dudx;
        self.dvdx = dvdx;
        let (dudy, dvdy) = solve_linear_system_2x2(mat_a, mat_by);
        self.dudy = dudy;
        self.dvdy = dvdy;

   }
}

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
    pub illumination_model: IlluminationModel,
    pub texture : Option<MipMap>
}

impl Material {
    pub fn create_diffuse(color : Vec3) -> Material {
        Material {
            name : String::from("diffuse"),
            diffuse: color,
            illumination_model: IlluminationModel::Lambertian,
            ..Material::default()
        }
    }

    pub fn get_diffuse(&self, surface : &Surface) -> Vec3 {
        let st = surface.uv;
        //compute texture differentials
        let dstdx = [surface.dudx, surface.dvdx];
        let dstdy = [surface.dudy, surface.dvdy];


        let width = f64::max(f64::max(f64::abs(dstdx[0]),
                                      f64::abs(dstdx[1])),
                             f64::max(f64::abs(dstdy[0]),
                                      f64::abs(dstdy[1])));

        if let Some(tex) = &self.texture {
            tex.sample_mipmap(st[0], st[1], width)
        } else {
            self.diffuse
        }
    }

    pub fn create_mirror() -> Material {
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
        let surface = Surface::new(ray.origin + ray.direction*t, 
                                self.normal, 
                                &self.material);
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
        Some(Surface::new(hit_position, normal, &self.material))
    }
}


fn solve_linear_system_2x2(a : [[f64;2];2], b : [f64;2]) -> (f64,f64)  {
       let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
       if f64::abs(det) < 1e-10 {
           return (0.,0.);
       }
       let x0 = (a[1][1] * b[0] - a[0][1] * b[1]) / det;
       let x1 = (a[0][0] * b[1] - a[1][0] * b[0]) / det;
       if f64::is_nan(x0) || f64::is_nan(x1) {
           return (0.,0.);
       }
       (x0, x1)
}
