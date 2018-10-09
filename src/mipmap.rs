use image;
use std::path::Path;
use std::ops::*;
use vec3::*;
use std::fs::File;
use std::io::*;


#[derive(Debug, Clone)]
pub struct MipMap {
    pyramid : Vec<Image>,
}

#[derive(Debug, Clone)]
pub struct Image {
   pub data : Vec<Vec3>,
   pub height : usize,
   pub width : usize
}

impl Image {
    pub fn save_ppm(&self, filename : &str) -> Result<()> {
        let to_ppm = |x:f64| (clamp(x)*255f64 + 0.5) as u8;

        let file = File::create(filename).expect(&format!("Unable to create file: {}", filename));
        let mut file = BufWriter::new(file);

        write!(file, "P3\n{} {}\n255\n", self.width, self.height)?;
        for pixel in self.data.iter() {
            write!(file, "{} {} {}\n", to_ppm(pixel.x) as u8, to_ppm(pixel.y) as u8, to_ppm(pixel.z) as u8)?;
        }
        write!(file, "")?;
        Ok(())
    }

    fn get_pixel(&self, s : usize, t : usize) -> Vec3 {
        if s >= self.width || t >= self.height {
            Vec3::zero()
        } else {  
            self[t][s]
        }
    }
}

impl Index<usize> for Image {
    type Output = [Vec3];

    fn index(&self, row: usize) -> &[Vec3] {
        let start = (self.height - row - 1)*self.width;
        &self.data[start .. start + self.width]
    }
}

impl MipMap  {
    pub fn load(file_name : &Path) -> Result<Self> {
        let img = image::open(file_name).expect(&format!("Couldn't find file: {}", file_name.display()));
        let img = img.to_rgb();
        let h = img.height() as usize;
        let w = img.width() as usize;

        let mut texels = Vec::with_capacity(img.len() / 4 );

        for chunk in img.chunks(3) {
            texels.push(Vec3::new(chunk[0] as f64/255., chunk[1] as f64/255., chunk[2] as f64/255.)); 
        }

        let image = Image {
            data : texels,
            height : h,
            width : w
        };

        Ok(MipMap::new(image))
    }

    pub fn new(image : Image) -> Self {

        let is_power_of_2 = |x: usize| { (x != 0) && ((x & (x - 1)) == 0) };
        if !is_power_of_2(image.height) || !is_power_of_2(image.width) {
            panic!(format!("Image needs to have power of two dimensions"));
        }
        let mut pyramid = Vec::new();
        pyramid.push(image);
        while pyramid[pyramid.len()-1].height > 1 {
            let next = MipMap::half_resolution(&pyramid[pyramid.len()-1]);
            pyramid.push(next);
        }

        MipMap { pyramid }
    }

    fn triangle(&self, level : usize, st : [f64;2])  -> Vec3 {
        //level = Clamp(level, 0, Levels() - 1);
        let s = st[0] * (self.pyramid[level].width as f64) - 0.5;
        let t = st[1] * (self.pyramid[level].height as f64) - 0.5;

        //Floor s & t
        let s0 = s as usize;
        let t0 = t as usize;

        let ds = s - s0 as f64;
        let dt = t - t0 as f64;
        (1. - ds) * (1. - dt) * self.pyramid[level].get_pixel(s0  ,t0  ) +
        (1. - ds) * dt        * self.pyramid[level].get_pixel(s0  ,t0+1) +
        ds       * (1. - dt)  * self.pyramid[level].get_pixel(s0+1,t0  ) +
        ds       * dt         * self.pyramid[level].get_pixel(s0+1,t0+1)
    }

    pub fn sample_mipmap(&self, u: f64, v: f64, width : f64) -> Vec3 {
        let n_levels = self.pyramid.len() as f64;

        let mut level = n_levels - 1. + f64::log2(f64::max(width, 1e-8));
        if level < 0. {
            level = 0.; 
        } else if level >= n_levels - 1. {
            level = n_levels - 1.;
        } else {
            //let iLevel = std::floor(level);
            //Float delta = level - iLevel;
            //return Lerp(delta, triangle(iLevel, st), triangle(iLevel + 1, st));
            level = f64::floor(level);
        }
        let level = level as usize;

        return self.triangle(level, [u, v]);
    }

    pub fn sample_nearest(&self, u: f64, v: f64) -> Vec3 {
        let level = 0;
        return self.triangle(level, [u, v]);
    }

    pub fn save_all(&self) -> Result<()> {
        for (i, image) in self.pyramid.iter().enumerate() {
            image.save_ppm(&format!("tex{}.ppm", i))?;
        }
        Ok(())
    }

    fn half_resolution(image : &Image) -> Image {
        let height = image.height/2;
        let width = image.width/2;
        let mut half_resolution = Image { data : Vec::with_capacity(height*width),
                               height,
                               width };
        for i in ( 0..image.height-1 ).rev().step_by(2) {
            for j in ( 0..image.width ).step_by(2) {
                half_resolution.data.push( (image[i][j]+image[i+1][j]+image[i][j+1]+image[i+1][j+1] ) / 4.);
            }
        }
        half_resolution
    }

}

//fn log2_int(num : u32) -> u32 {
//    let mut num = num;
//    let mut targetlevel = 0;
//    while num >= 1 {
//        targetlevel += 1;
//        num = num >> 1;
//    };
//    targetlevel
//}

fn clamp(x : f64) -> f64 {
    if x < 0. { 0. } 
    else if x > 1. { 1. }
    else { x }
}
