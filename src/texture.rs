use image;
use image::GenericImage;
use std::path::Path;
use vec3::Vec3;
use std::ops::*;
use mipmap::*;

pub enum ImageWrap {
    Repeat, Black, Clamp
}

#[derive(Debug, Clone)]
pub struct Texture {
    mipmap : MipMap
}

impl Texture {



}

fn log2_int(num : usize) -> usize {
    let mut num = num;
    let mut targetlevel = 0;
    while (num >= 1) {
        targetlevel += 1;
        num = num >> 1;
    };
    targetlevel
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
