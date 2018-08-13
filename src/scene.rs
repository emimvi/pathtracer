use vec3::*;
use geometry::*;
use lights::*;
use obj;
use std::path::Path;

pub struct Scene {
    pub objects : Vec<Box<Intersectable>>,
    pub lights : PointLight,
    pub background : Vec3
}


impl Scene {
    pub fn trace_closest(&self, mut ray : &mut Ray) -> Option<Surface> {
        self.objects.iter()
                .fold(None, |prev, o| o.intersect(&mut ray).or(prev))
    }

    pub fn trace_any(&self, mut ray : &mut Ray) -> Option<Surface> {
        for o in &self.objects {
            let surface = o.intersect(&mut ray);
            if surface.is_some() {
                return surface;
            }
        }
        None
    }

    pub fn cornell_box() -> Scene {
        let mut objects : Vec<Box<Intersectable>> = obj::load_obj(&Path::new("/Users/Imbert/Documents/Playground/rust/tracer/models/Test.obj"));
        let sphere = Sphere { center : Vec3::new(350., 100., 350.), 
                              radius : 100., 
                              material : Material::create_mirror() };
        objects.push(Box::new(sphere));
        let background = Vec3::new(0.3, 0.3, 0.7);
        let lights = PointLight { position : Vec3::new(250., 490., 250.), intensity : 2000000. };
        Scene { objects, lights, background }
    }

    pub fn default() -> Scene {
        let red  = Vec3::new(1.,0.,0.);
        let _green = Vec3::new(0.,1.,0.);
        let blue  = Vec3::new(0.,0.,1.);
        let gray  = Vec3::new(0.5,0.5,0.5);
        //let white = Vec3::new(1f64,1f64,1f64);
        let background = Vec3::new(0.3, 0.3, 0.7);
        //Circle
        let z = 5.;
        let circle = Sphere { center : Vec3::new(0., 0., z), 
                              radius : 0.5, 
                              material : Material::create_mirror() };
        let circle2 = Sphere { center : Vec3::new(0.5, 0.0, z-1.), 
                              radius : 0.25, 
                              material : Material::create_diffuse( _green) };
        let floor = Plane { origin : Vec3::new(0., -1., 0.), 
                            normal : Vec3::new(0., 1., 0.).normalize(), 
                            material : Material::create_diffuse( gray ) };
        let roof = Plane { origin : Vec3::new(0., 2., 0.), normal : Vec3::new(0., -1., 0.).normalize(), material : Material::create_diffuse( gray ) };
        let back = Plane { origin : Vec3::new(0., 0., z+3.), normal : Vec3::new(0., 0., -1.).normalize(), material : Material::create_diffuse( gray ) };
        let right = Plane { origin : Vec3::new(2., 0., 0.), 
                            normal : Vec3::new(-1., 0., 0.).normalize(), 
                            material : Material::create_diffuse( blue )};
        let left = Plane { origin : Vec3::new(-2., 0., 0.), 
                           normal : Vec3::new(1., 0., 0.).normalize(), 
                           material : Material::create_diffuse( red )};

        let objects : Vec<Box<Intersectable>> = vec!(Box::new(floor), Box::new(circle), Box::new(circle2), Box::new(roof), Box::new(back), Box::new(right), Box::new(left));
        //let light = DirectionalLight { direction : Vec3::new(0., -1.0, 0.), radiance : 1. };
        let lights = PointLight { position : Vec3::new(0., 1.0, z), intensity : 100. };

        Scene { objects, lights, background }
    }
}
