use std::path::Path;
use tobj;
use vec3::*;
use geometry::*;
use mipmap::*;

pub fn load_obj(obj_name : &Path) -> Vec<Box<dyn Intersectable>> {
    let obj = tobj::load_obj(obj_name);
    let (models, materials) = obj.expect(&format!("Failed to load {}: ", obj_name.display()));

    let mut meshes : Vec<Box<Intersectable>> = Vec::new();
    models.iter().for_each(|model| {

        //Vertices: Collect flattened vertices into Vec3's
        let positions = &model.mesh.positions;
        let mut vertices : Vec<Vec3> = Vec::new();
        for xyz in positions.chunks(3) {
            let (v0,v1,v2) = (xyz[0] as f64, xyz[1] as f64, xyz[2] as f64);
            vertices.push(Vec3::new(v0,v1,v2));
        };

        //Faces: Collect flattened indices into triplets
        let indices = &model.mesh.indices;
        let mut faces : Vec<[usize ; 3]> = Vec::with_capacity(indices.len()/3);
        for xyz in indices.chunks(3) {
            let face = [xyz[0] as usize, xyz[1] as usize, xyz[2] as usize];
            faces.push(face);
        }

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


        let tri_mesh = TriangleMesh{
            vertices,
            faces,
            material,
            texture_coordinates
        };
        //if let Some(_) = tri_mesh.texture_coordinates {
        //    println!("{:?}", tri_mesh.vertices);
        //    println!("{:?}", tri_mesh.faces);
        //    println!("{:?}", tri_mesh.texture_coordinates);
        //}
        meshes.push(Box::new(tri_mesh))
    });

    meshes
}

#[derive(Debug)]
struct TriangleMesh {
    //TODO: Add mesh name
    vertices : Vec<Vec3>,
    faces : Vec<[usize; 3]>,
    texture_coordinates : Option<Vec<Vec2>>,
    material : Material
}

impl Intersectable for TriangleMesh {
    fn intersect(&self, ray : &mut Ray) -> Option<Surface> {
        self.faces.iter().fold(None, |prev : Option<Surface>, &face_id| {
            let v0 = &self.vertices[face_id[0]];
            let v1 = &self.vertices[face_id[1]];
            let v2 = &self.vertices[face_id[2]];
            let e0 = *v1 - *v0;
            let e1 = *v0 - *v2;
            let normal = e0.cross(e1);

            if f64::abs(ray.direction.dot(normal)) < 1.0e-12 {
                return prev;
            }
            let t = (*v0 - ray.origin).dot(normal)/ray.direction.dot(normal);

            // Check distance to intersection
            if t < ray.t_min || t > ray.t_max {
                return prev;
            }

            // Find barycentric coordinates
            let tmp = (*v0 - ray.origin).cross(ray.direction);
            let v = tmp.dot(e1)/ray.direction.dot(normal);
            let w = tmp.dot(e0)/ray.direction.dot(normal);
            let u = 1. - v - w;
            if v < 0.0 || w < 0.0 || v + w > 1.0 {
                return prev;
            }
            ray.t_max = t;

            let uvs : [Vec2 ; 3] = if let Some(tex_coords) = &self.texture_coordinates {
                [tex_coords[face_id[0]], tex_coords[face_id[1]], tex_coords[face_id[2]]]
            } else { 
                [ Vec2::new(0.,0.), Vec2::new(1.0, 0.), Vec2::new(1., 1.) ] 
            };
            let uv = uvs[0]*u + uvs[1]*v + uvs[2]*w;

            //Partial derivatives: pbrt triangle.cpp
            let duv02 = uvs[0] - uvs[2];
            let duv12 = uvs[1] - uvs[2];
            let (dp02,dp12) = (*v0 - *v2, *v1 - *v2);
            let determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
            let inverse_det = 1. / determinant;
            let dpdu = ( duv12[1] * dp02 - duv02[1] * dp12) * inverse_det;
            let dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * inverse_det;

            let mut surface = Surface{ 
                position: ray.origin + ray.direction*t, 
                normal : -normal.normalize(),// Let the counterclockwise wound side face forward
                material : &self.material,
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
        })
    }
}
