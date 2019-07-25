use vec3::*;
use std::ops::*;

#[derive(Debug)]
struct Transform {
    matrix: [[f32; 4] ; 4]
}


impl Transform {

    pub fn identity() -> Transform {
        Transform { matrix :
        [[ 1., 0., 0., 0.],
         [ 0., 1., 0., 0.],
         [ 0., 0., 1., 0.],
         [ 0., 0., 0., 1.]] }
    }

    pub fn translate(v : &Vec3) -> Transform {
        Transform { matrix :
        [[ 1., 0., 0., v.x],
         [ 0., 1., 0., v.y],
         [ 0., 0., 1., v.z],
         [ 0., 0., 0., 1.]] }
    }

    pub fn scale(v : &Vec3) -> Transform {
        Transform { matrix :
        [[ v.x, 0., 0., 0.],
         [ 0., v.y, 0., 0.],
         [ 0., 0., v.z, 0.],
         [ 0., 0., 0., 1.]] }
    }

    pub fn rotate_x(radians : f32) -> Transform {
        let sin = f32::sin(radians);
        let cos = f32::cos(radians);
        Transform { matrix :
        [[ 1. , 0.  , 0.   , 0.],
         [ 0. , cos , -sin , 0.],
         [ 0. , sin , cos  , 0.],
         [ 0. , 0.  , 0.   , 1.]] }
    }

    pub fn rotate_y(radians : f32) -> Transform {
        let sin = f32::sin(radians);
        let cos = f32::cos(radians);
        Transform { matrix :
        [[ cos  , 0. , sin , 0.],
         [ 0.   , 1. , 0.  , 0.],
         [ -sin , 0. , cos , 0.],
         [ 0.   , 0. , 0.  , 1.]] }
    }

    pub fn rotate_z(radians : f32) -> Transform {
        let sin = f32::sin(radians);
        let cos = f32::cos(radians);
        Transform { matrix :
        [[ cos , -sin , 0. , 0.],
         [ sin , cos  , 0. , 0.],
         [ 0.  , 0.   , 1. , 0.],
         [ 0.  , 0.   , 0. , 1.]] }
    }

    pub fn look_at(stand_at : Vec3, look_at : Vec3, up_vector : Vec3) -> Transform {
        let view_direction = (look_at - stand_at).normalize();
        let side_vector = view_direction.cross(up_vector).normalize();
        let orthorgonal_up = side_vector.cross(view_direction).normalize();

        let mut m : [[f32; 4] ; 4] = [[0. ; 4]; 4]; //Zero initialize
        m[0][3] = stand_at.x;
        m[1][3] = stand_at.y;
        m[2][3] = stand_at.z;
        m[3][3] = 1.;

        m[0][0] = side_vector.x;
        m[1][0] = side_vector.y;
        m[2][0] = side_vector.z;
        m[3][0] = 0.;
        m[0][1] = orthorgonal_up.x;
        m[1][1] = orthorgonal_up.y;
        m[2][1] = orthorgonal_up.z;
        m[3][1] = 0.;
        m[0][2] = view_direction.x;
        m[1][2] = view_direction.y;
        m[2][2] = view_direction.z;
        m[3][2] = 0.;

        Transform { matrix : m }
    }

}

impl Mul for Transform {
    type Output = Transform;
    fn mul(self, other: Transform) -> Transform {
        let mut m : [[f32; 4] ; 4] = [[0. ; 4]; 4]; //Zero initialize
        for i in 0..4 {
            for j in 0..4 {
                let mut sum = 0.;
                for k in 0..4 {
                    sum += self[i][k]*other[k][j];
                }
                m[i][j] = sum;
            }
        }
        Transform { matrix : m }
    }
}

impl Mul<Vec3> for Transform {
    type Output = Vec3;
    fn mul(self, other: Vec3) -> Vec3 {
        let mut v = [0.; 3];
        for i in 0..3 {
            let mut sum = 0.;
            for j in 0..3 {
                sum += self[i][j]*other[j];
            }
            v[i] = sum;
        }
        Vec3::new(v[0], v[1], v[2])
    }
}

impl Index<usize> for Transform {
    type Output = [f32];

    fn index(&self, idx : usize) -> &Self::Output {
        &self.matrix[idx]
    }
}

#[test]
fn mat_mult() {
    let t1 = Transform { matrix : 
        [[ 1., 1., 1., 1.],
         [ 2., 2., 2., 2.],
         [ 4., 4., 4., 4.],
         [ 5., 5., 5., 5.]] };
    let t2 = Transform { matrix : 
        [[ 1., 2., 3., 4.],
         [ 1., 2., 3., 4.],
         [ 1., 2., 3., 4.],
         [ 1., 2., 3., 4.]] };
    let res = t1*t2;
    let res_cmp = Transform { matrix :
        [[ 1., 8., 12., 16.],
         [ 8., 16., 24., 32.],
         [ 16., 32., 46., 64.],
         [ 20., 40., 60., 60.]] };
    assert_eq!(res[2][1], res_cmp[2][1]);
    assert_eq!(res[2][3], res_cmp[2][3]);
}
