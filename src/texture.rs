use image;
use image::GenericImage;
use std::path::Path;
use vec3::Vec3;

#[derive(Debug, Clone)]
pub struct Texture {
    data : Vec<[f64 ; 4]>,
    height : usize,
    width : usize,
}

impl Texture {

    pub fn load(file_name : &Path) -> Texture {
        let img = image::open(file_name).expect(&format!("Couldn't find file: {}", file_name.display()));
        let img_data : Vec<f64> = img.raw_pixels().iter().map(|&pixel| pixel as f64/255.).collect();
        let mut vec4 = Vec::new();
        for chunk in img_data.chunks(4) {
            vec4.push([chunk[0], chunk[1], chunk[2], chunk[3]]); 
        }
        Texture {
            data : vec4,
            height : img.height() as usize,
            width : img.width() as usize,
        }
    }

    pub fn sample_nearest(&self, u: f64, v: f64) -> Vec3 {
        let u = u - f64::floor(u);
        let v = v - f64::floor(v);

        let mut u_int = (u*self.width  as f64 + 0.5)  as usize;
        let mut v_int = (v*self.height as f64 + 0.5) as usize;
        if u_int == self.width {
            u_int =  0;
        }
        if v_int == self.height {
            v_int =  0;
        }
        let idx = u_int + (self.height - v_int - 1)*self.width;
        let color = self.data[idx];

        Vec3::new(color[0],color[1],color[2])
    }
}

pub fn test() {
    // Use the open function to load an image from a Path.
    // ```open``` returns a `DynamicImage` on success.
    let _img1 = Texture::load(&Path::new("./models/lorem.png"));

    // The dimensions method returns the images width and height.
   // println!("dimensions {:?}", img.dimensions());

   // // The color method returns the image's `ColorType`.

   // println!("{:?}", img.get_pixel(1000, 1000));
}
