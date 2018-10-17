#![allow(unused)]
use std::path::Path;
use tobj;
use vec3::*;
use geometry::*;
use mipmap::*;

use std::sync::{Weak, RwLock, Arc, LockResult, RwLockReadGuard};
use bvh::{Geometry, BBox, Boundable, BVH};

struct TriangleMesh(Arc<RwLock<TriMesh>>);
impl TriangleMesh {
    fn new(vertices : Vec<Vec3> , faces : Vec<[usize ; 3]>, material : Material, texture_coordinates : Option<Vec<Vec2>>) -> TriangleMesh {
        let mesh = TriangleMesh (
            Arc::new(RwLock::new(TriMesh::new(vertices, vec!(), material, texture_coordinates)))
        );
        {
            let bvh = BVH::unanimated(16, {
                faces.iter().map(|face| {
                    Triangle {
                        face : *face,
                        parent : Arc::downgrade(&(mesh.0))
                    }
                }).collect()});

            let mut p = mesh.0.write().unwrap();
            (*p).faces = bvh;
        }
        mesh
    }

}

pub struct TriMesh {
    vertices : Vec<Vec3>,
    faces : BVH<Triangle>,
    texture_coordinates : Option<Vec<Vec2>>,
    material : Arc<Material>
}

struct Triangle {
    parent : Weak<RwLock<TriMesh>>,
    face : [usize ; 3]
}


impl Intersectable for TriangleMesh {

    #[inline]
    fn intersect(&self, ray : &mut Ray) -> Option<Surface> {
        let mesh = self.0.read().unwrap();
        mesh.intersect(ray)
    }
}

impl Intersectable for TriMesh {

    #[inline]
    fn intersect(&self, ray : &mut Ray) -> Option<Surface> {
        self.faces.iter().fold(None, |prev : Option<Surface>, triangle| {
            triangle.intersect(ray).or(prev)
        })
    }
}

impl TriMesh {
    fn new(vertices : Vec<Vec3> , _faces : Vec<[usize ; 3]>, material : Material, texture_coordinates : Option<Vec<Vec2>>) -> TriMesh {
        TriMesh{
            vertices,
            faces : BVH::empty(),
            material: Arc::new(material),
            texture_coordinates
        }
    }
}

impl Boundable for TriangleMesh {
    fn bounds(&self, _a : f32, _b : f32) -> BBox {
        self.0.read().unwrap().bounds(_a,_b)
    }
}

impl Boundable for TriMesh {
    fn bounds(&self, _a : f32, _b : f32) -> BBox {
        self.faces.bounds(_a, _b)
    }
}

impl Boundable for Triangle {
    fn bounds(&self, _ : f32, _ : f32) -> BBox {
        let (v0, v1, v2) = self.get_vertices();
        BBox::singular(v0.into())
            .point_union(&v1.into())
            .point_union(&v2.into())

    }
}

impl Triangle {
    
    //Utility function to avoid having to handle the deep structure of TriangleMesh.
    //Pass in a function that takes a &TriMesh as argument, and returns whatever you want from the
    //mesh. This is neccesary as we cannot return the reference directly due to the RwLockGuard
    //going out of scope when the function exits.
    fn parent<F, R>(&self, func : F) -> R 
        where F : Fn(&TriMesh) -> R  {
        let tri_mesh = &*self.parent.upgrade().unwrap();
        let mesh = tri_mesh.read().unwrap();
        func(&mesh)
    }

    fn get_vertices(&self) -> (Vec3, Vec3, Vec3) {
        self.parent(|mesh|{
            (mesh.vertices[self.face[0]],
             mesh.vertices[self.face[1]],
             mesh.vertices[self.face[2]])
        })
    }

    fn get_texture_coordinates(&self) -> [Vec2 ; 3] {
        self.parent(|mesh|{
            let uvs : [Vec2 ; 3] = if let Some(tex_coords) = &mesh.texture_coordinates {
                [tex_coords[self.face[0]], tex_coords[self.face[1]], tex_coords[self.face[2]]]
            } else { 
                [ Vec2::new(0.,0.), Vec2::new(1.0, 0.), Vec2::new(1., 1.) ] 
            };
            uvs
        })
    }

    fn get_material(&self) -> Arc<Material> {
        self.parent(|mesh|{
            Arc::clone(&mesh.material)
        })
    }
}

impl Intersectable for Triangle {
    fn intersect(&self, ray : &mut Ray) -> Option<Surface> {
        let (v0, v1, v2) = self.get_vertices();
        let e0 = v1 - v0;
        let e1 = v0 - v2;
        let normal = e0.cross(e1);

        if f64::abs(ray.direction.dot(normal)) < 1.0e-12 {
            return None;
        }
        let t = (v0 - ray.origin).dot(normal)/ray.direction.dot(normal);

        // Check distance to intersection
        if t < ray.t_min || t > ray.t_max {
            return None;
        }

        // Find barycentric coordinates
        let tmp = (v0 - ray.origin).cross(ray.direction);
        let v = tmp.dot(e1)/ray.direction.dot(normal);
        let w = tmp.dot(e0)/ray.direction.dot(normal);
        let u = 1. - v - w;
        if v < 0.0 || w < 0.0 || v + w > 1.0 {
            return None;
        }
        ray.t_max = t;

        let uvs = self.get_texture_coordinates();
        let uv = uvs[0]*u + uvs[1]*v + uvs[2]*w;

        //Partial derivatives: pbrt triangle.cpp
        let duv02 = uvs[0] - uvs[2];
        let duv12 = uvs[1] - uvs[2];
        let (dp02,dp12) = (v0 - v2, v1 - v2);
        let determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
        let inverse_det = 1. / determinant;
        let dpdu = ( duv12[1] * dp02 - duv02[1] * dp12) * inverse_det;
        let dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * inverse_det;

        let mut surface = Surface{ 
            position: ray.origin + ray.direction*t, 
            normal : -normal.normalize(),// Let the counterclockwise wound side face forward
            material : Arc::clone(&self.get_material()),
            uv,
            dpdu : dpdu ,
            dpdv : dpdv ,
            dpdx : Vec3::zero(),
            dpdy : Vec3::zero(),
            dudx : 0.0,
            dudy : 0.0,
            dvdx : 0.0,
            dvdy : 0.0 };

        surface.calculate_differentials(&ray);
        Some(surface)
    }

}

pub fn load_obj(obj_name : &Path) -> Vec<Geometry> {
    let obj = tobj::load_obj(obj_name);
    let (models, materials) = obj.expect(&format!("Failed to load {}: ", obj_name.display()));

    let mut meshes : Vec<Geometry> = Vec::new();
    models.iter().for_each(|model| {

        //Vertices: Collect flattened vertices into Vec3's
        let positions = &model.mesh.positions;
        let mut vertices : Vec<Vec3> = Vec::new();
        for xyz in positions.chunks(3) {
            let (v0,v1,v2) = (xyz[0] as f64, xyz[1] as f64, xyz[2] as f64);
            vertices.push(Vec3::new(v0,v1,v2));
        };

        //Material:
        let material = if let Some(id) = model.mesh.material_id {
            if let Some(illum) = &materials[id].illumination_model {
            let tobj_material = &materials[id];
            Material {
                name : tobj_material.name.clone(),
                diffuse: Vec3::from(tobj_material.diffuse),
                specular: Vec3::from(tobj_material.specular),
                ambient: Vec3::from(tobj_material.ambient),
                illumination_model: IlluminationModel::from(*illum),
                texture : if !tobj_material.diffuse_texture.is_empty() {
                    let mip_map = match MipMap::load(Path::new(&tobj_material.diffuse_texture)) {
                        Ok(map) => map,
                        Err(err) => panic!("Failed to load MipMap: {}", err)
                    };
                    Some(mip_map) 
                } else { None }
            }
            } else {Material::default() }
        } else {
            Material::default()
        };

        //Texture 
        let tex_coords = &model.mesh.texcoords;
        let texture_coordinates : Option<Vec<Vec2>> = if tex_coords.len() > 0 {
            let mut uvs : Vec<Vec2> = Vec::new();
            for uv in tex_coords.chunks(2) {
                let coord = Vec2::new(uv[0] as f64, uv[1] as f64);
                uvs.push(coord);
            }
            Some(uvs)
        } else { None };



        //Faces: Collect flattened indices into triplets
        let indices = &model.mesh.indices;
        let mut faces : Vec<[usize ; 3]> = Vec::with_capacity(indices.len()/3);
        for xyz in indices.chunks(3) {
            let face = [xyz[0] as usize, xyz[1] as usize, xyz[2] as usize];
            faces.push(face);
        }

        let tri_mesh = TriangleMesh::new(vertices, faces, material, texture_coordinates);       


        let geo : Geometry = Geometry::new(tri_mesh);

        meshes.push(geo);
    });

    meshes
}



