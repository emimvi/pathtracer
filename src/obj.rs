use std::path::Path;
use tobj;
use vec3::*;
use geometry::*;
use texture::*;

pub fn load_obj(obj_name : &Path) -> Vec<Box<dyn Intersectable>> {
    let obj = tobj::load_obj(obj_name);
    let (models, materials) = obj.expect(&format!("Failed to load {}: ", obj_name.display()));

    //tobj::print_model_info(&models, &materials);

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
                illumination_model: Some(IlluminationModel::from(*illum)),
                texture : if !tobj_material.diffuse_texture.is_empty() {
                    Some(Texture::load(Path::new(&tobj_material.diffuse_texture))) 
                } else { None }
            }
            } else {Material::default() }
        } else {
            Material::default()
        };

        //Texture 
        let tex_coords = &model.mesh.texcoords;
        let texture_coordinates : Option<Vec<[f64 ; 2]>> = if tex_coords.len() > 0 {
            let mut uvs : Vec<[f64 ; 2]> = Vec::new();
            for uv in tex_coords.chunks(2) {
                let coord = [uv[0] as f64, uv[1] as f64];
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
        if let Some(_) = tri_mesh.texture_coordinates {
            println!("{:?}", tri_mesh.vertices);
            println!("{:?}", tri_mesh.faces);
            println!("{:?}", tri_mesh.texture_coordinates);
        }
        meshes.push(Box::new(tri_mesh))
    });

    meshes
}

#[derive(Debug)]
struct TriangleMesh {
    //TODO: Add mesh name
    vertices : Vec<Vec3>,
    faces : Vec<[usize; 3]>,
    texture_coordinates : Option<Vec<[f64; 2]>>,
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

            let uv = if let Some(tex_coords) = &self.texture_coordinates {
                let uvs = [tex_coords[face_id[0]], tex_coords[face_id[1]], tex_coords[face_id[2]]];
                let res = (uvs[0][0]*u + uvs[1][0]*v + uvs[2][0]*w,  uvs[0][1]*u + uvs[1][1]*v + uvs[2][1]*w);
                res
            } else { (0.,0.) };

            // Let the counterclockwise wound side face forward
            Some(Surface{ 
                position: ray.origin + ray.direction*t, 
                normal : -normal.normalize(),
                material : &self.material,
                uv
            })
        })
    }
}

pub fn test_load() {
    use std::env;

    // We assume that we are in a valid directory.
    let path = env::current_dir().unwrap();
    println!("The current directory is {}", path.display());
    let cornell_box = tobj::load_obj(&Path::new("./models/CornellBox.obj"));
    assert!(cornell_box.is_ok());
    let (models, materials) = cornell_box.unwrap();

    println!("# of models: {}", models.len());
    println!("# of materials: {}", materials.len());
    for (i, m) in models.iter().enumerate() {
        let mesh = &m.mesh;
        println!("model[{}].name = \'{}\'", i, m.name);
        println!("model[{}].mesh.material_id = {:?}", i, mesh.material_id);

        println!("Size of model[{}].indices: {}", i, mesh.indices.len());
        // Normals and texture coordinates are also loaded, but not printed in this example
        println!("model[{}].vertices: {}", i, mesh.positions.len() / 3);
        assert!(mesh.positions.len() % 3 == 0);
    }

    for (i, m) in materials.iter().enumerate() {
        println!("material[{}].name = \'{}\'", i, m.name);
        println!("    material.Ka = ({}, {}, {})", m.ambient[0], m.ambient[1], m.ambient[2]);
        println!("    material.Kd = ({}, {}, {})", m.diffuse[0], m.diffuse[1], m.diffuse[2]);
        println!("    material.Ks = ({}, {}, {})", m.specular[0], m.specular[1], m.specular[2]);
        println!("    material.Ns = {}", m.shininess);
        println!("    material.d = {}", m.dissolve);
        println!("    material.map_Ka = {}", m.ambient_texture);
        println!("    material.map_Kd = {}", m.diffuse_texture);
        println!("    material.map_Ks = {}", m.specular_texture);
        println!("    material.map_Ns = {}", m.normal_texture);
        println!("    material.map_d = {}", m.dissolve_texture);
        for (k, v) in &m.unknown_param {
            println!("    material.{} = {}", k, v);
        }
    }
}
