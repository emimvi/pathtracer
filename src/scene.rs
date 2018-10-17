use vec3::*;
use geometry::*;
use lights::*;
use obj;
use std::path::Path;

use bvh::{BVH, Geometry};

pub struct Scene {
    pub objects : BVH<Geometry>,
    pub lights : Box<Light>,
    pub background : Vec3,
    pub camera : Camera 
}

pub struct Camera {
    pub eye : Vec3,
    pub view_direction : Vec3,
    pub up_direction : Vec3,
}
impl Camera {
    pub fn new(eye : Vec3, view_direction : Vec3, up_direction : Vec3) -> Camera {
        Camera {
            eye,
            view_direction,
            up_direction
        }
    }

    pub fn ray(&self, x : f64, y: f64, pixel_delta_x : f64, pixel_delta_y: f64) -> Ray {
        let camera_const = 1.;

        let view_dir = self.view_direction * camera_const;

        let ip_axes0 = (self.view_direction.cross(self.up_direction)).normalize();
        let ip_axes1 = (ip_axes0.cross(self.view_direction)).normalize();

        //normalize(ip_normal + ip_axes[0]*coords[0] + ip_axes[1]*coords[1]);
        let ray_dir = (view_dir + ip_axes0*x + ip_axes1*y).normalize();

        let mut ray = Ray::new(self.eye, ray_dir);

        ray.rx_origin = self.eye;
        ray.ry_origin = self.eye;
        ray.rx_direction = (Vec3::new(pixel_delta_x, 0., 0.) + view_dir + ip_axes0*x + ip_axes1*y).normalize();
        ray.ry_direction = (Vec3::new(0., pixel_delta_y, 0.) + view_dir + ip_axes0*x + ip_axes1*y).normalize();
        ray
    }
}

impl Scene {
    pub fn trace_closest(&self, mut ray : &mut Ray) -> Option<Surface> {
        self.objects.intersect(&mut ray, |r, i| {
            i.intersect(r)
        })
    }

    //TODO: Return on first intersection, instead of checking all possible intersections.
    pub fn trace_any(&self, ray : &mut Ray) -> Option<Surface> {
        self.trace_closest(ray)
    }

    pub fn glossy_planes() -> Scene {
        let eye = Vec3::new(0.,4.,1.);
        let direction = Vec3::new(0., -1., 1.).normalize();
        let camera = Camera::new(eye, direction, Vec3::new(0., 1., 0.));

        let objects : Vec<Geometry> = obj::load_obj(&Path::new("/Users/Imbert/Documents/DTU/thesis/tracer/models/glossy_planes.obj"));
        let bvh = BVH::unanimated(16, objects);        

        let background = Vec3::zero();
        let lights = Box::new( DirectionalLight { direction : Vec3::new(0., -1., 1.).normalize(), radiance : 5. } );
        Scene { objects : bvh, lights, background , camera }
    }

    //pub fn cornell_box() -> Scene {
    //    let mut objects : Vec<Box<Intersectable>> = obj::load_obj(&Path::new("./models/Test.obj"));
    //    let sphere = Sphere { center : Vec3::new(350., 100., 350.), 
    //                          radius : 100., 
    //                          material : Material::create_mirror() };
    //    objects.push(Box::new(sphere));
    //    //let background = Vec3::new(0.3, 0.3, 0.7);
    //    let background = Vec3::zero();
    //    let lights = Box::new( PointLight { position : Vec3::new(250., 490., 250.), intensity : 2000000. } );
    //    Scene { objects, lights, background }
    //}

    //pub fn default() -> Scene {
    //    let red  = Vec3::new(1.,0.,0.);
    //    let _green = Vec3::new(0.,1.,0.);
    //    let blue  = Vec3::new(0.,0.,1.);
    //    let gray  = Vec3::new(0.5,0.5,0.5);
    //    //let white = Vec3::new(1f64,1f64,1f64);
    //    let background = Vec3::new(0.3, 0.3, 0.7);
    //    //Circle
    //    let z = 5.;
    //    let circle = Sphere { center : Vec3::new(0., 0., z), 
    //                          radius : 0.5, 
    //                          material : Material::create_mirror() };
    //    let circle2 = Sphere { center : Vec3::new(0.5, 0.0, z-1.), 
    //                          radius : 0.25, 
    //                          material : Material::create_diffuse( _green) };
    //    let floor = Plane { origin : Vec3::new(0., -1., 0.), 
    //                        normal : Vec3::new(0., 1., 0.).normalize(), 
    //                        material : Material::create_diffuse( gray ) };
    //    let roof = Plane { origin : Vec3::new(0., 2., 0.), normal : Vec3::new(0., -1., 0.).normalize(), material : Material::create_diffuse( gray ) };
    //    let back = Plane { origin : Vec3::new(0., 0., z+3.), normal : Vec3::new(0., 0., -1.).normalize(), material : Material::create_diffuse( gray ) };
    //    let right = Plane { origin : Vec3::new(2., 0., 0.), 
    //                        normal : Vec3::new(-1., 0., 0.).normalize(), 
    //                        material : Material::create_diffuse( blue )};
    //    let left = Plane { origin : Vec3::new(-2., 0., 0.), 
    //                       normal : Vec3::new(1., 0., 0.).normalize(), 
    //                       material : Material::create_diffuse( red )};

    //    let objects : Vec<Box<Intersectable>> = vec!(Box::new(floor), Box::new(circle), Box::new(circle2), Box::new(roof), Box::new(back), Box::new(right), Box::new(left));
    //    //let light = DirectionalLight { direction : Vec3::new(0., -1.0, 0.), radiance : 1. };
    //    let lights = Box::new(PointLight { position : Vec3::new(0., 1.0, z), intensity : 100. });

    //    Scene { objects, lights, background }
    //}
}
