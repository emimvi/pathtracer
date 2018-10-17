use std::sync::Arc;
use geometry::{Material, Ray};
use vec3::*;

#[derive(Debug, Clone)]
pub struct Surface { pub position : Vec3,
                         pub normal: Vec3,
                         pub material :  Arc<Material>,
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
impl Surface {
    pub fn new(position : Vec3, normal : Vec3, material : Arc<Material>) -> Surface {
        Surface { position, 
                  normal, 
                  material : Arc::clone(&material),
                  uv : Vec2::from([0., 0.]),
                  dpdu : Vec3::zero(), //How much does the world position change pr. uv-coodinate
                  dpdv : Vec3::zero(),
                  dpdx : Vec3::zero(), //How much does the world position AT AN INTERSECTION change wrt. a given camera ray.
                  dpdy : Vec3::zero(),
                  dudx : 0.,
                  dudy : 0.,
                  dvdx : 0.,
                  dvdy : 0.
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
