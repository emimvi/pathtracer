use algebra::vec3::*;
use geometry::*;
use lights::*;
use material::*;
use microfacet::*;
use mipmap::*;
use rand;
use std::path::Path;
use tobj;

use bvh::{BBox, Boundable, Geometry, BVH};

pub fn load_obj(obj_name: &Path) -> (Vec<Geometry>, Vec<Box<dyn Light>>) {
    let obj = tobj::load_obj(obj_name);
    let (models, materials) =
        obj.unwrap_or_else(|_| panic!("Failed to load {}: ", obj_name.display()));

    let mut meshes: Vec<Geometry> = Vec::new();
    let mut lights: Vec<Box<dyn Light>> = Vec::new();
    models.iter().for_each(|model| {
        //Vertices: Collect flattened vertices into Vec3's
        let positions = &model.mesh.positions;
        let mut vertices: Vec<Vec3> = Vec::new();
        for xyz in positions.chunks(3) {
            let (v0, v1, v2) = (f64::from(xyz[0]), f64::from(xyz[1]), f64::from(xyz[2]));
            vertices.push(Vec3::new(v0, v1, v2));
        }

        //Material:
        //let material = if let Some(id) = model.mesh.material_id {
        //    if let Some(illum) = &materials[id].illumination_model {
        //    let tobj_material = &materials[id];
        //    Material {
        //        name : tobj_material.name.clone(),
        //        diffuse: Vec3::from(tobj_material.diffuse),
        //        specular: Vec3::from(tobj_material.specular),
        //        ambient: Vec3::from(tobj_material.ambient),
        //        illumination_model: IlluminationModel::from(*illum),
        //        roughness: tobj_material.shininess as f64,
        //        micro : Some(Glossy::new(tobj_material.shininess as f64, 1.5, Vec3::from(tobj_material.diffuse))),
        //        texture : if !tobj_material.diffuse_texture.is_empty() {
        //            let mip_map = match MipMap::load(Path::new(&tobj_material.diffuse_texture)) {
        //                Ok(map) => map,
        //                Err(err) => panic!("Failed to load MipMap: {}", err)
        //            };
        //            Some(mip_map)
        //        } else { None }
        //        , ..Material::default()
        //    }
        //    } else {Material::default() }
        //} else {
        //    Material::default()
        //};
        //
        let default = Box::new(Lambertian {
            color: Vec3::zero(),
            emission: Vec3::zero(),
        });
        let material: Box<dyn _Material> = if let Some(id) = model.mesh.material_id {
            if let Some(illum) = &materials[id].illumination_model {
                let tobj_material = &materials[id];
                // if (tobj_material.ambient[0] == 0.) {
                //     Box::new(Glossy::new(0.5, 0.5, Vec3::from(tobj_material.diffuse)))
                // } else
                {
                    Box::new(Lambertian {
                        color: Vec3::from(tobj_material.diffuse),
                        emission: Vec3::from(tobj_material.ambient),
                    })
                }
            } else {
                default
            }
        } else {
            default
        };

        //Texture
        let tex_coords = &model.mesh.texcoords;
        let texture_coordinates: Option<Vec<Vec2>> = if !tex_coords.is_empty() {
            let mut uvs: Vec<Vec2> = Vec::new();
            for uv in tex_coords.chunks(2) {
                let coord = Vec2::new(f64::from(uv[0]), f64::from(uv[1]));
                uvs.push(coord);
            }
            Some(uvs)
        } else {
            None
        };

        //Faces: Collect flattened indices into triplets
        let indices = &model.mesh.indices;
        let mut faces: Vec<[usize; 3]> = Vec::with_capacity(indices.len() / 3);
        for xyz in indices.chunks(3) {
            let face = [xyz[0] as usize, xyz[1] as usize, xyz[2] as usize];
            faces.push(face);
        }

        let tri_mesh = TriMesh::new(vertices, faces, texture_coordinates);

        let geo: Geometry = Geometry::new(tri_mesh, material);

        meshes.push(geo);
    });

    (meshes, lights)
}
