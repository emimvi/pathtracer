use algebra::vec3::*;
use geometry::*;
use lights::*;
use material::*;
use microfacet::*;
use mipmap::*;
use rand;
use std::path::Path;

use bvh::{BBox, Boundable, Geometry, BVH};

pub struct TriMesh {
    vertices: Vec<Vec3>,
    faces: BVH<Triangle>,
    texture_coordinates: Option<Vec<Vec2>>,
}

unsafe impl Send for Triangle {}
unsafe impl Sync for Triangle {}
pub struct Triangle {
    parent: *const TriMesh,
    face: [usize; 3],
}

impl Intersectable for Box<TriMesh> {
    #[inline]
    fn intersect(&self, ray: &mut Ray) -> Option<Surface> {
        self.faces.intersect(ray)
    }
}

impl TriMesh {
    pub fn new(
        vertices: Vec<Vec3>,
        faces: Vec<[usize; 3]>,
        texture_coordinates: Option<Vec<Vec2>>,
    ) -> Box<TriMesh> {
        let mut mesh = Box::new(TriMesh {
            vertices,
            faces: BVH::empty(),
            texture_coordinates,
        });
        {
            let bvh = BVH::new({
                faces
                    .iter()
                    .map(|face| Triangle {
                        face: *face,
                        parent: &*mesh,
                    })
                    .collect()
            });

            mesh.faces = bvh;
        }
        mesh
    }
}

impl Boundable for Box<TriMesh> {
    fn bounds(&self, _a: f32, _b: f32) -> BBox {
        self.faces.bounds(_a, _b)
    }
}

impl Boundable for Triangle {
    fn bounds(&self, _: f32, _: f32) -> BBox {
        let (v0, v1, v2) = self.get_vertices();
        BBox::singular(v0).point_union(v1).point_union(v2)
    }
}

impl Triangle {
    pub fn get_vertices(&self) -> (Vec3, Vec3, Vec3) {
        let mesh = unsafe { &*self.parent };
        (
            mesh.vertices[self.face[0]],
            mesh.vertices[self.face[1]],
            mesh.vertices[self.face[2]],
        )
    }

    fn get_texture_coordinates(&self) -> [Vec2; 3] {
        let mesh = unsafe { &*self.parent };

        let uvs: [Vec2; 3] = if let Some(tex_coords) = &mesh.texture_coordinates {
            [
                tex_coords[self.face[0]],
                tex_coords[self.face[1]],
                tex_coords[self.face[2]],
            ]
        } else {
            [Vec2::new(0., 0.), Vec2::new(1.0, 0.), Vec2::new(1., 1.)]
        };

        uvs
    }
}

impl Intersectable for Triangle {
    fn intersect(&self, ray: &mut Ray) -> Option<Surface> {
        let (v0, v1, v2) = self.get_vertices();
        let e0 = v1 - v0;
        let e1 = v0 - v2;
        let normal = e0.cross(e1);

        if f64::abs(ray.direction.dot(normal)) < 1.0e-12 {
            return None;
        }
        let t = (v0 - ray.origin).dot(normal) / ray.direction.dot(normal);

        // Check distance to intersection
        if t < ray.t_min || t > ray.t_max {
            return None;
        }

        // Find barycentric coordinates
        let tmp = (v0 - ray.origin).cross(ray.direction);
        let denom = 1. / ray.direction.dot(normal);
        let v = tmp.dot(e1) * denom;
        let w = tmp.dot(e0) * denom;
        let u = 1. - v - w;
        if v < 0.0 || w < 0.0 || v + w > 1.0 {
            return None;
        }
        ray.t_max = t;

        let uvs: [Vec2; 3] = self.get_texture_coordinates();

        //Partial derivatives: pbrt p. 163 triangle.cpp
        let (dp02, dp12): (Vec3, Vec3) = (v0 - v2, v1 - v2);
        let (duv02, duv12): (Vec2, Vec2) = (uvs[0] - uvs[2], uvs[1] - uvs[2]);
        let determinant: f64 = duv02[0] * duv12[1] - duv02[1] * duv12[0];
        let inverse_det: f64 = 1. / determinant;
        let dpdu: Vec3 = (duv12[1] * dp02 - duv02[1] * dp12) * inverse_det;
        let dpdv: Vec3 = (-duv12[0] * dp02 + duv02[0] * dp12) * inverse_det;

        let uv = uvs[0] * u + uvs[1] * v + uvs[2] * w;
        let mut surface = Surface {
            position: ray.origin + ray.direction * t,
            normal: -normal.normalize(), // Let the counterclockwise wound side face forward
            uv,
            dpdu,
            dpdv,
            ..Surface::default()
        };

        surface.calculate_differentials(&ray);
        Some(surface)
    }
}
