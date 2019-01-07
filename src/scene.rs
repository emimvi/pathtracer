use vec3::*;
use geometry::*;
use lights::*;
use mipmap::*;
use obj;
use std::path::Path;
use std::sync::Arc;

use material::*;

use bvh::{BVH, Geometry};

pub struct Scene {
    objects : BVH<Geometry>,
    pub lights : Vec<Box<Light>>,
    pub background : Vec3,
    pub camera : Camera 
}

pub struct Camera {
    eye : Vec3,
    view_direction : Vec3,
    up_direction : Vec3,
    camera_const : f64, //FOV
}

impl Camera {
    pub fn new(eye : Vec3, view_direction : Vec3, up_direction : Vec3, camera_const : f64) -> Camera {
        Camera {
            eye,
            view_direction,
            up_direction,
            camera_const
        }
    }

    pub fn ray(&self, x : f64, y: f64, pixel_delta_x : f64, pixel_delta_y: f64) -> Ray {

        let ip_axes0 = -(self.view_direction.cross(self.up_direction)).normalize();
        let ip_axes1 = -(ip_axes0.cross(self.view_direction)).normalize();


        //normalize(ip_normal + ip_axes[0]*coords[0] + ip_axes[1]*coords[1]);
        let view_dir = self.view_direction * self.camera_const;
        let ray_dir = (view_dir + ip_axes0*x + ip_axes1*y).normalize();

        let mut ray = Ray::new(self.eye, ray_dir);

        if ::CONFIG.use_diffs {
            ray.rx_origin = self.eye;
            ray.ry_origin = self.eye;
            ray.rx_direction = (view_dir + ip_axes0*(x+pixel_delta_x) + ip_axes1*y).normalize();
            ray.ry_direction = (view_dir + ip_axes0*x + ip_axes1*(y+pixel_delta_y)).normalize();
        }

        ray
    }
}

impl Scene {
    pub fn trace_closest(&self, ray : &mut Ray) -> Option<Surface> {
        self.objects.intersect(ray)
    }

    //TODO: Return on first intersection, instead of checking all possible intersections.
    pub fn trace_any(&self, ray : &mut Ray) -> Option<Surface> {
        self.trace_closest(ray)
    }

    pub fn glossy_planes() -> Scene {
        let eye = Vec3::new(0.,3.,2.);
        let direction = Vec3::new(0., -1., 1.).normalize();
        let camera = Camera::new(eye, direction, Vec3::new(0., 1., 0.), 1.);
        
        //let eye = Vec3::new(0.,3.,-7.);
        //let direction = Vec3::new(0., -0.2, 1.).normalize();
        //let camera = Camera::new(eye, direction, Vec3::new(0., 1., 0.), 4.);

        let (mut objects, _) = obj::load_obj(&Path::new("/Users/Imbert/Documents/DTU/thesis/tracer/models/glossy_planes.obj"));

        //let mut material = Material::constant(Vec3::new(0.5, 0.5, 0.5));
        //let rainbow = MipMap::rainbow();
        //material.texture = Some(rainbow);
        //let arc = Arc::new(material);
        //for mut geometry in &mut objects.iter_mut() {
        //    if geometry.material.texture.is_some() {
        //        geometry.material = Arc::clone(&arc);
        //    }
        //}
        
        let sphere = Sphere { center : Vec3::new(0., 0.5, 3.), 
                              radius : 0.5, 
                              };
        //objects.push(Geometry::new(sphere      , Material::mirror()));

        let bvh = BVH::new(objects);
        //let bvh = BVH::unanimated(16, objects);        
        //

        let background = Vec3::zero();
        //let lights : Vec<Box<Light>> = vec!(Box::new( DirectionalLight { direction : Vec3::new(0., -1., 1.).normalize(), radiance : 5. } ));
        let lights : Vec<Box<Light>> = Vec::new();
        Scene { objects : bvh, lights, background , camera }
    }

    pub fn quad() -> Scene {
        let eye = Vec3::new(0.0,0.0,0.);
        let direction = Vec3::new(0., 0., 1.).normalize();
        let camera = Camera::new(eye, direction, Vec3::new(0., 1., 0.), 1.);

        let (mut objects, _) = obj::load_obj(&Path::new("/Users/Imbert/Documents/DTU/thesis/tracer/models/quad.obj"));

        let sphere = Sphere { center : Vec3::new(-0.5, 0.0, 0.), 
                              radius : 0.5, 
                              };
        let sphere_g = Sphere { center : Vec3::new(0.5, 0.0, 0.), 
                              radius : 0.5, 
                              };
        //objects.push(Geometry::new(sphere      , Material::mirror()));
        //objects.push(Geometry::new(sphere_g      , Material::glass()));

        let mut material = Material::constant(Vec3::new(0.5, 0.5, 0.5));
        let rainbow = MipMap::rainbow();
        material.texture = Some(rainbow);
        let arc = Arc::new(material);
        for mut geometry in &mut objects.iter_mut() {
            geometry.material = Arc::clone(&arc);
        }


        let bvh = BVH::new(objects);
        //let bvh = BVH::unanimated(16, objects);        

        let background = Vec3::one()*0.33;
        let lights : Vec<Box<Light>> = vec!(Box::new( DirectionalLight { direction : Vec3::new(0., -1., 1.).normalize(), radiance : 3. } ));
        //let lights : Vec<Box<Light>> = Vec::new();
        //let lights : Box<Light> = Box::new(DirectionalLight { direction : Vec3::new(0., 0., 1.), radiance : 3. });
        Scene { objects : bvh, lights, background , camera }
    }

    pub fn cornell_box() -> Scene {
        let (mut objects, _) = obj::load_obj(&Path::new("./models/Test.obj"));

        // let mut material = Material::diffuse(Vec3::new(0.5, 0.5, 0.5));
        // let rainbow = MipMap::rainbow();
        // material.texture = Some(rainbow);
        // let arc = Arc::new(material);
        // for mut geometry in &mut objects.iter_mut() {
        //     geometry.material = Arc::clone(&arc);
        // }

        let sphere = Sphere { center : Vec3::new(350., 100., 350.), 
                              radius : 100., 
                              };
        let sphere_diff = Sphere { center : Vec3::new(150., 100., 350.), 
                              radius : 100., 
                              };

        let sphere_glass = Sphere { center : Vec3::new(200., 100., 175.), 
                              radius : 100., 
                              };
        let sphere_gloss = Sphere { center : Vec3::new(400., 375., 250.), 
                              radius : 100., 
                              };


        //objects.push(Geometry::new(sphere      , Material::micro(Vec3::one(), 0.07, 1.5)));
        //objects.push(Geometry::new(sphere      , Material::mirror()));
        objects.push(Geometry::new(sphere_glass, Material::glass()));
    
        //objects.push(Geometry::new(sphere, Material::diffuse(GREEN)));
        let bvh = BVH::new(objects);

        //let background = Vec3::new(0.3, 0.3, 0.7);
        let background = Vec3::zero();
        let lights : Box<Light> = Box::new(PointLight { position : Vec3::new(255., 300., 55.), intensity : 250000. });
        //let lights : Box<Light> = Box::new(DirectionalLight { direction : Vec3::new(0., 0., 1.), radiance : 3. });

        let eye = Vec3::new(275.,275.,-600.);
        let direction = Vec3::new(0., 0., 1.).normalize();
        let camera = Camera::new(eye, direction, Vec3::new(0., 1., 0.), 2.);
        Scene { objects : bvh, lights : vec!(), background, camera }
    }

    pub fn default() -> Scene {
        let _red  = Vec3::new(1.,0.,0.);
        let _green = Vec3::new(0.,1.,0.);
        let blue  = Vec3::new(0.,0.,1.);
        let _gray  = Vec3::new(0.5,0.5,0.5);
        let gold = Vec3::new(1.00,0.71,0.29);
        let white = Vec3::new(1f64,1f64,1f64);
        let background = Vec3::new(0.3, 0.3, 0.7);
        //Circle
        let z = 1.5;
        let geo = {
            let circle = Sphere { center : Vec3::new(0.0, 0., z), radius : 0.75 };
            let roughness = 0.1;
            let ior = 1.5;
            let mat_diff = Material::micro(white, roughness, ior);
            Geometry::new(circle, mat_diff)
        };
        let geo2 = {
            let circle = Sphere { center : Vec3::new(1.5, 0., z-1.), radius : 0.75 };
            let mat_diff = Material::diffuse(blue);
            Geometry::new(circle, mat_diff)
        };
        let geo3 = {
            let below = 1000.;
            let circle = Sphere { center : Vec3::new(0., -below-0.75, z), radius : below };
            let mat_diff = Material::diffuse(_red);
            Geometry::new(circle, mat_diff)
        };


        let objects : Vec<Geometry> = vec!(geo, geo2, geo3);
        let bvh = BVH::new(objects);        

        let eye = Vec3::new(0.,0.,0.);
        let direction = Vec3::new(0., 0., 1.).normalize();
        let camera = Camera::new(eye, direction, Vec3::new(0., 1., 0.), 2.);

        let light : Box<Light> = Box::new(DirectionalLight { direction : Vec3::new(0., -1.0, 0.), radiance : ::f64::consts::PI });
        //let light : Box<Light> = Box::new(PointLight { position : Vec3::new(0., 0., 0.), intensity : 200.});
        let lights = vec!(light);
        Scene { objects : bvh, lights, background, camera }
    }
}
