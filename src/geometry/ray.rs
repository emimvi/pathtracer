use algebra::vec3::*;
use std::f64;

#[derive(Debug, Clone)]
pub struct Ray {
    pub origin: Vec3,
    pub direction: Vec3,
    pub t_min: f64,
    pub t_max: f64,
    pub rx_origin: Vec3, //Differential footprint
    pub ry_origin: Vec3,
    pub rx_direction: Vec3, //Differential direction
    pub ry_direction: Vec3,
}

impl Ray {
    pub fn new(origin: Vec3, direction: Vec3) -> Ray {
        Ray {
            origin,
            direction,
            t_min: 1e-5,
            t_max: f64::MAX,
            rx_origin: Vec3::zero(),
            ry_origin: Vec3::zero(),
            rx_direction: Vec3::zero(),
            ry_direction: Vec3::zero(),
        }
    }
}
