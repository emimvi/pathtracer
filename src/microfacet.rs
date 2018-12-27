use vec3::*;
use std::mem;
use rand;
use ray::Ray;
use surface::Surface;
use std::f64;

#[derive(Debug, Clone)]
pub struct Glossy<T = Beckmann> 
where T : MicrofacetDistribution
{
    distribution : T,
    fresnel : FresnelConductor,
    reflectance_color : Vec3
}

impl Glossy {
    
    pub fn new(roughness : f64, ior : f64, color : Vec3) -> Self {
        let _gold_ior = Vec3::new(4.0624, 2.2039, 1.8534);
        let _gold_k = Vec3::new(0.13100, 0.55758, 1.31);

        Glossy {
            distribution : Beckmann::new(roughness),
            fresnel : FresnelConductor::new(Vec3::one()*32., Vec3::one()*0.1),
            reflectance_color : color,
        }
    }

    ///Sample the material, given a view direction. Returns a color, sampled direction and
    ///probability/pdf of the sampled direction.
    pub fn shade(&self, ray : &Ray, surface : &Surface) -> (Vec3, Vec3, f64) {
        let normal = surface.normal;
        let wo = -ray.direction;

        let half_vector = self.distribution.sample_microfacet_normal(wo, normal);
        let sampled_wi = (-wo).reflect(half_vector);

        let pdf = self.distribution.pdf(half_vector, normal) / (4. * wo.dot(half_vector));
        let color = self.f(wo, sampled_wi, normal);


        (color, sampled_wi, pdf)
    }

    pub fn f(&self, wo : Vec3, wi : Vec3, normal : Vec3) -> Vec3 {

        let wh = (wo+wi).normalize();
        
        let cos_theta_o = f64::abs(wo.dot(normal)); //abs_cosTheta(wo);
        let cos_theta_i = f64::abs(wi.dot(normal)); //abs_cosTheta(wi);

        if cos_theta_i == 0. || cos_theta_o == 0. { 
            return Vec3::zero();
        }
        if wh[0] == 0. && wh[1] == 0. && wh[2] == 0. {
            return Vec3::zero();
        }

        let fresnel = self.fresnel.evaluate(wo.dot(wh));
        let d = self.distribution.D(wh, normal);
        let g = {
            let t = self.distribution.G(wo, wi, normal);
            if t < 0. { 0. } else if t > 1. { 1. } else { t }
        };
        
        assert!(g <= 1.0+1e-3 && g >= 0., "G outside of 0..1: {}", g);

        let r = self.reflectance_color * d * g * fresnel /
               (4. * cos_theta_i * cos_theta_o);
        r
    }
}

pub trait MicrofacetDistribution {

    fn sample_microfacet_normal(&self, out_direction : Vec3, normal : Vec3) -> Vec3;
    #[allow(non_snake_case)]
    fn D(&self, half_vector : Vec3, macro_normal : Vec3) -> f64;
    //#[allow(non_snake_case)]
    //fn G1(&self, w : Vec3, vh : Vec3, n : Vec3) -> f64;

    #[allow(non_snake_case)]
    fn G(&self, wo : Vec3, wi : Vec3, n :Vec3) -> f64 {

        let hv = (wi+wo).normalize();
        self.G1(wo, hv, n)*self.G1(wi, hv, n)
        //1. / (1. + self.lambda(wo, n) + self.lambda(wi, n) )
    }

    //fn lambda(&self, w : Vec3, n : Vec3) -> f64;
    #[allow(non_snake_case)]
    fn G1(&self, wo : Vec3, wi : Vec3, n : Vec3) -> f64;

    fn pdf(&self, half_vector : Vec3, macro_normal : Vec3) -> f64 {
        self.D(half_vector, macro_normal)*abs_cos_theta(&half_vector, &macro_normal)
    }
}

#[derive(Debug, Clone)]
pub struct Beckmann {
    width : f64
}
impl Beckmann {
    pub fn new(roughness : f64) -> Self {
        let roughness = f64::max(roughness, 1e-5);
        //let mut width = f64::ln(roughness);
        //width = 1.62142 + 0.819955 * width + 0.1734 * width * width +
        //       0.0171201 * width * width * width + 0.000640711 * width * width * width * width;
        Beckmann { width : roughness}
    }

}
impl MicrofacetDistribution for Beckmann {

    //fn lambda(&self, w : Vec3, n : Vec3) -> f64 {
    //    let abs_tan_theta = f64::abs(tan_theta(&w,&n));
    //    if f64::is_infinite(abs_tan_theta) {
    //        return 0.;
    //    }
    //    let a = 1. / (self.width * abs_tan_theta);
    //    if a >= 1.6 {
    //        return 0.;
    //    }
    //    (1. - 1.259 * a + 0.396 * a * a) /
    //        (3.535 * a + 2.181 * a * a)
    //}

    ///Samples a half vector for use in microfacet calculations.
    fn sample_microfacet_normal(&self, _ : Vec3, normal : Vec3) -> Vec3 {
        // Get random numbers
        let tan_theta_sqr = -self.width*self.width*f64::ln(1. - rand::random::<f64>());
        let phi = 2.0*f64::consts::PI*rand::random::<f64>();

        // Calculate sampled half-angle vector as if the z-axis were the normal
        let cos_theta = 1. / f64::sqrt(1. + tan_theta_sqr);
        let sin_theta = f64::sqrt(f64::max(0., 1. - cos_theta * cos_theta));
        let hv = spherical_direction(sin_theta, cos_theta, phi);

        // Rotate from z-axis to actual normal
        (hv.rotate_to(&normal)).normalize()
    }

    fn D(&self, m_norm : Vec3, n_norm : Vec3) -> f64 {
	let cos_theta_m = Vec3::dot(m_norm, n_norm);

	if cos_theta_m <= 0. {
		return 0.;
        }
	let cos4theta_m = f64::powi(cos_theta_m, 4);
	let theta_m = f64::acos(cos_theta_m);
	let tan2theta_m = f64::tan(theta_m) * f64::tan(theta_m);
	let alpha2 = self.width*self.width;

	f64::exp(-tan2theta_m / alpha2) / (f64::consts::PI * alpha2*cos4theta_m)
    }

    fn G1(&self, w : Vec3, m_norm : Vec3, n_norm : Vec3) -> f64 {
        let cos_theta_v = Vec3::dot(w, n_norm);

        if Vec3::dot(w, m_norm) / cos_theta_v <= 0.0 {
        	return 0.;
        }

        let theta_v = f64::acos(cos_theta_v);
        let a = 1. / (self.width*f64::tan(theta_v) );
        if a >= 1.6 {
        	return 1.;
        }

        (3.535*a + 2.181*a*a) /
         (1. + 2.276*a + 2.577*a*a)
    }

}

#[derive(Debug, Clone)]
struct GGX { 
    roughness : f64
}
impl GGX {
    pub fn new(roughness : f64) -> Self {
        GGX { roughness }
    }
}
impl MicrofacetDistribution for GGX {

    ///Samples a half vector for use in microfacet calculations.
    fn sample_microfacet_normal(&self, _ : Vec3, normal : Vec3) -> Vec3 {

        let xi1 = rand::random::<f64>();
        let tan_theta_sqr = self.roughness*self.roughness*xi1/(1. - xi1);
        let phi = 2.0*rand::random::<f64>()*f64::consts::PI;

        // Calculate sampled half-angle vector as if the z-axis were the normal
        let cos_theta_sqr = 1./(1. + tan_theta_sqr);
        let cos_theta = f64::sqrt(cos_theta_sqr);
        let sin_theta = f64::sqrt(1. - cos_theta_sqr);
        let hv = spherical_direction(sin_theta, cos_theta, phi);

        // Rotate from z-axis to actual normal
        hv.rotate_to(&normal).normalize()
    }

    fn D(&self, half_vector : Vec3, macro_normal : Vec3) -> f64 {
        let cos_theta = half_vector.dot(macro_normal);
        if cos_theta <= 0.0 { 
            return 0.
        }

        let alpha_g_sqr = self.roughness*self.roughness;

        let tan_theta_sqr = tan2_theta(&half_vector, &macro_normal);
        if tan_theta_sqr.is_infinite() {
            return 0.;
        }
        
        let alpha_tan2 = alpha_g_sqr + tan_theta_sqr;

        alpha_g_sqr /
            ( f64::consts::PI * f64::powi(cos_theta, 4) * (alpha_tan2*alpha_tan2) )
    }

    fn G1(&self, w : Vec3, vh : Vec3, _:Vec3) -> f64 {
        let cos_theta = w.dot(vh);
        let alpha_g_sqr = self.roughness*self.roughness;
        let cos_theta_sqr = cos_theta*cos_theta;
        let tan_theta_sqr = (1. - cos_theta_sqr)/cos_theta_sqr;
        2./(1. + f64::sqrt(1. + alpha_g_sqr*tan_theta_sqr))
    }

}
trait Fresnel {
    fn evaluate(&self, cos_theta_i : f64) -> Vec3;
}

#[derive(Debug, Clone)]
struct FresnelConductor {
    ior_i : Vec3, //Index of refraction for the object the fresnel belongs to
    ior_t : Vec3, //Index of refraction for the object the fresnel is placed in (eg. ~1 for air)
    absorption : Vec3
}

impl FresnelConductor {
    pub fn new(ior : Vec3, absorption : Vec3) -> Self {
        FresnelConductor {
            ior_t : Vec3::one(),
            ior_i : ior,
            absorption
        }
    }

    pub fn conductor(cos_theta_i : f64, ior_i : Vec3, ior_t : Vec3, absorption : Vec3) -> Vec3 {
        let eta = ior_t / ior_i;
        let etak = absorption / ior_i;

        let cosThetaI2 = cos_theta_i * cos_theta_i;
        let sinThetaI2 = Vec3::one()*(1. - cosThetaI2);
        let eta2 = eta * eta;
        let etak2 = etak * etak;

        let t0 = eta2 - etak2 - sinThetaI2;
        let a2plusb2 = (t0 * t0 + 4. * eta2 * etak2).map(&f64::sqrt);
        let t1 = a2plusb2 + Vec3::one()*cosThetaI2;
        let a = (0.5 * (a2plusb2 + t0)).map(&f64::sqrt);
        let t2 = 2. * cos_theta_i * a;
        let Rs = (t1 - t2) / (t1 + t2);

        let t3 = cosThetaI2 * a2plusb2 + sinThetaI2 * sinThetaI2;
        let t4 = t2 * sinThetaI2;
        let Rp = Rs * (t3 - t4) / (t3 + t4);

        return 0.5 * (Rp + Rs);
    }

}

impl Fresnel for FresnelConductor {
    fn evaluate(&self, cos_theta_i : f64) -> Vec3 {
        FresnelConductor::conductor(cos_theta_i, self.ior_i, self.ior_t, self.absorption)
    }
}

#[derive(Debug, Clone)]
pub struct FresnelDielectric {
    ior_i : f64, //Index of refraction for the object the fresnel belongs to
    ior_t : f64, //Index of refraction for the object the fresnel is placed in (eg. ~1 for air)
}

impl Fresnel for FresnelDielectric {
    fn evaluate(&self, cos_theta_i : f64) -> Vec3 {
        Vec3::one()*self.dielectric(cos_theta_i)
    }
}

impl FresnelDielectric {

    ///Create a fresnel for an object placed in vacuum/air
    pub fn new(ior : f64) -> FresnelDielectric {
        FresnelDielectric{
            ior_i : 1.,
            ior_t : ior
        }
    }

    pub fn eta(&self, cos_theta : f64) -> f64 {
        let inside = cos_theta < 0.;

        if !inside {
            self.ior_i/self.ior_t
        } else {
            self.ior_t/self.ior_i
        }
    }

    pub fn dielectric(&self, cos_theta_i : f64) -> f64 {
        //If cos_theta_i is in [0..1.), we are on the outside. If [-1..0) on the inside.
        let inside = cos_theta_i < 0.;

        //If we are inside, swap values so we treat the inside material as incident direction.
        let mut ior_i = self.ior_i;
        let mut ior_t = self.ior_t;
        let mut cos_theta_i = cos_theta_i;
        if inside {
            mem::swap(&mut ior_i, &mut ior_t);
            cos_theta_i = f64::abs(cos_theta_i);
        }
        let eta = ior_i/ior_t;;

        //Compute cos_theta_t using Snell’s law
        let sin_theta_i = f64::sqrt(f64::max(0., 1. - cos_theta_i * cos_theta_i));
        let sin_theta_t = ior_i / ior_t * sin_theta_i;

        let sin_theta_t_sqr = eta*eta*(1. - cos_theta_i * cos_theta_i);

        //Handle total internal reflection
        if sin_theta_t >= 1. {
            return 1.;
        }

        let cos_theta_t = f64::sqrt(f64::max(0., 1. - sin_theta_t_sqr));
        let fresnel_par =  ((ior_t * cos_theta_i) - (ior_i * cos_theta_t)) /
                           ((ior_t * cos_theta_i) + (ior_i * cos_theta_t));
        let fresnel_perp = ((ior_i * cos_theta_i) - (ior_t * cos_theta_t)) / 
                           ((ior_i * cos_theta_i) + (ior_t * cos_theta_t));  

        let res = (fresnel_par * fresnel_par + fresnel_perp * fresnel_perp) / 2.;
        res
    }
}

fn cos_theta(w1 : &Vec3, w2 : &Vec3) -> f64 {
    w1.dot(*w2)
}
fn cos2_theta(w1 : &Vec3, w2 : &Vec3) -> f64 {
    let cos_theta = cos_theta(w1, w2);
    cos_theta * cos_theta 
}
fn abs_cos_theta(w1 : &Vec3, w2 : &Vec3) -> f64 {
    f64::abs(cos_theta(w1, w2)) 
}
fn sin2_theta(w1 : &Vec3, w2 : &Vec3) -> f64 {
    f64::max(0., 1. - cos2_theta(w1, w2))
}
#[allow(unused)]
fn sin_theta(w1 : &Vec3, w2 : &Vec3) -> f64 {
    f64::sqrt(sin2_theta(w1, w2)) 
}
#[allow(unused)]
fn tan_theta(w1 : &Vec3, w2 : &Vec3) -> f64 {
    sin_theta(w1, w2) / cos_theta(w1, w2) 
}
fn tan2_theta(w1 : &Vec3, w2 : &Vec3) -> f64 {
    sin2_theta(w1, w2) / cos2_theta(w1, w2)
}

// Given spherical coordinates, where theta is the 
// polar angle and phi is the azimuthal angle, this
// function returns the corresponding direction vector
fn spherical_direction(sin_theta : f64, cos_theta : f64, phi : f64) -> Vec3
{
  Vec3::new(sin_theta*f64::cos(phi), sin_theta*f64::sin(phi), cos_theta)
}


#[test]
fn test_fresnel() {
    let f = FresnelDielectric::new(1.5);

    let r = f.dielectric(0.7071067811865476);
    let exp = 0.0502399;

    assert!(f64::abs(r - exp) < 0.0001, "Fresnel error: Got {}, expected {}", r, exp);
}
