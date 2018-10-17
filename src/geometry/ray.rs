use vec3::Vec3;
use std::f64;

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
}

