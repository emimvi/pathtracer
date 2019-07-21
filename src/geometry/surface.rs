use std::sync::Arc;
use std::rc::Rc;
use geometry::{Material, Ray};
use material::*;
use vec3::*;

#[derive(Debug, Clone)]
pub struct Surface<'a> { 
    pub position : Vec3,
    pub normal: Vec3,
    pub material :  Option<&'a dyn _Material>,
    pub uv : Vec2,
    pub dpdu : Vec3, pub dpdv : Vec3, // world space / texture space
    pub dndu : Vec3, pub dndv : Vec3, // normal derivative / texture space
    pub dpdx : Vec3, pub dpdy : Vec3, // world space / screen space
    pub dudx : f64 , pub dvdx : f64,  // texture space / screen space x
    pub dudy : f64 , pub dvdy : f64,  // texture space / screen space y
}

impl Default for Surface<'static> {
    fn default() -> Surface<'static> {
        Surface { 
            position : Vec3::zero(), 
            normal : Vec3::new(0.,1.,0.), 
            uv : Vec2::from([0., 0.]),
            material: None,
            dpdu : Vec3::zero(), //How much does the world position change pr. uv-coodinate
            dpdv : Vec3::zero(),
            dpdx : Vec3::zero(), //How much does the world position AT AN INTERSECTION change wrt. a given camera ray.
            dpdy : Vec3::zero(),
            dndu : Vec3::zero(),
            dndv : Vec3::zero(),
            dudx : 0.,
            dudy : 0.,
            dvdx : 0.,
            dvdy : 0.
        }
    }
}

impl Surface<'_> {
    pub fn new(position : Vec3, normal : Vec3) -> Surface<'static> {
        Surface { position, 
                  normal, 
                  ..Surface::default()
        }
    }

   pub fn calculate_differentials(&mut self, ray : &Ray) {
       if !::CONFIG.use_diffs {
           return;
       }

       //pbrt. p601
       let d = self.normal.dot(self.position);
       //Distance of the dx and dy ray to the tangent plane of the surface.
       let tx = -(self.normal.dot(ray.rx_origin) - d) /
                self.normal.dot(ray.rx_direction);
       let ty = -(self.normal.dot(ray.ry_origin) - d) /
                self.normal.dot(ray.ry_direction);

       let px = ray.rx_origin + tx * ray.rx_direction;
       let py = ray.ry_origin + ty * ray.ry_direction;
       if f64::is_infinite(tx) || f64::is_nan(tx) || f64::is_infinite(ty) || f64::is_nan(ty){
           self.dudx = 0.;
           self.dvdx = 0.;
           self.dudy = 0.;
           self.dvdy = 0.;
           self.dpdx = Vec3::zero();
           self.dpdy = Vec3::zero();
           return;
       }
       let dpdx = px - self.position;
       let dpdy = py - self.position;
       self.dpdx = dpdx;
       self.dpdy = dpdy;

       //pbrt. p603
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
       let bx = [ self.dpdx[dim.0], self.dpdx[dim.1] ]; 
       let by = [ self.dpdy[dim.0], self.dpdy[dim.1] ]; 
       let (dudx, dvdx) = solve_linear_system_2x2(mat_a, bx);
       let (dudy, dvdy) = solve_linear_system_2x2(mat_a, by);

       self.dudx = dudx;
       self.dvdx = dvdx;
       self.dudy = dudy;
       self.dvdy = dvdy;

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
