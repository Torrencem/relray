#![feature(euclidean_division)]
#![allow(non_snake_case)]
#![allow(unused_variables)]
#![allow(dead_code)]

extern crate nalgebra as na;
extern crate image;

use image::ImageBuffer;

use std::f64;

use na::Vector3;
use na::Point3;
use na::Rotation3;
use na::Unit;

// Units: meters kilograms seconds
static CONST_C: f64 = 299792458f64;
static CONST_G: f64 = 6.67408e-11f64;

static SOLAR_MASS: f64 = 2e30f64;

// For a real scale black hole
static EX_MD_BH_MASS: f64 = 100.0 * SOLAR_MASS;

struct GeoSim {
    // Initial psi coordinate
    psi_init: f64,
    // Current psi coordinate
    psi: f64,
    // Current r coordinate
    r: f64,
    // Schwartzchild radius of black hole
    rs: f64,
    // b (see paper)
    b: f64,
    // time step for psi
    dpsi: f64,
    // running in reverse (see below)
    reversed: bool,
}

// Given some initial conditions, create a GeoSim,
// which is a simulator for the path of a photon in
// r,psi 2D radial coordinates
fn geosim(R: f64, psi_init: f64, mut alpha: f64, rs: f64, mut dpsi: f64) -> GeoSim {
    // reversed is to deal with paths that approach the
    // black hole before departing. These correspond to
    // angles greater than pi halves. These are simulated
    // the exact same, except backwards (since light works
    // that way).
    let mut reversed = false;
    if alpha >= 3.1415926f64 / 2.0 {
        alpha = 3.1415926f64 - alpha;
        dpsi *= -1.0;
        reversed = true;
    }
    // b is the parameter which completely encodes the
    // emission angle of the light
    let b = R * alpha.sin() / (1.0 - rs / R).sqrt();
    GeoSim {
        psi_init: psi_init,
        psi: psi_init,
        r: R,
        rs: rs,
        b: b,
        dpsi: dpsi,
        reversed: reversed,
    }
}

impl Iterator for GeoSim {
    // psi,r pairs
    type Item = (f64, f64);
    
    fn next(&mut self) -> Option<(f64, f64)> {
        // The actual photon simulation:
        // takes place in r-psi radial 2D coordinates

        // Catch odd positions and errors
        if !self.r.is_finite() {
            return None;
        }
        // We want to yield our current position,
        // then update it
        let toyield = if !self.reversed {
            (self.psi, self.r)
        } else {
            // If we're reversed, we still want to
            // travel in the increasing psi direction
            (2.0 * self.psi_init - self.psi, self.r)
        };
        
        // Adjusted from formula in paper:
        // compute the change in r with respect
        // to a change in psi
        let mut dr_dpsi =  self.r * self.r * 
            (1.0 / self.b / self.b - 1.0 / self.r / self.r * (1.0 - self.rs / self.r)).sqrt();
        
        // Numerical weirdness: if we've been computing backwards,
        // and now our formula gives NaN, it's because
        // we're at our periastron. Swap the direction
        // to be going "forward" and nudge the angle a
        // tiny bit to keep from staying in NaN land.
        if self.reversed && (1.0 / self.b / self.b - 1.0 / self.r / self.r * (1.0 - self.rs / self.r)) < 0.0 {
            // turn into a forward iterator
            self.b = self.r / (1.0 - self.rs / self.r).sqrt() - 0.003;
            // for this iteration, our radius doesn't change
            dr_dpsi = 0.0;
            self.reversed = false;
            self.dpsi *= -1.0;
            // fix psi
            self.psi = 2.0 * self.psi_init - self.psi;
        }

        // Increase psi by dpsi
        self.psi += self.dpsi;
        // Increase our radius appropriately
        self.r += dr_dpsi * self.dpsi;
        // Give the value we computed last step
        Some(toyield)
    }
}

// Convert a mass of a blackhole into a schwarzchild radius
fn mass_to_rs(mass: f64) -> f64 {
    2.0 * CONST_G * mass / CONST_C / CONST_C
}

// Radial to cartesian coordinates
fn rad_to_cart(psi: f64, r: f64) -> (f64, f64) {
    (r * psi.cos(), r * psi.sin())
}

// Cartesian to radial coordinates
fn cart_to_rad(x: f64, y: f64) -> (f64, f64) {
    (y.atan2(x), (x * x + y * y).sqrt())
}

fn main() {
    // A couple of notes: formula (2) and (3) from this paper:
    // https://www.aanda.org/articles/aa/pdf/2016/11/aa29075-16.pdf
    // is what's used in this program. If you're running this,
    // decrease the resolution a bit and use 'cargo run --release'
    // 
    // Note that light rays are calculated as being emitted from the
    // camera and hitting objects. Obviously that's not how photons
    // actually work, but it's identical in reverse due to optics.

    // For debugging: Trace out a single photon's orbit
    // in 2D radial space for plotting in sage
    // let rs = 1.0;
    // let dpsi = 0.001f64;
    // for (psi, r) in geosim(rs * 2.0, 3.14159 / 2.0, 3.14159 / 1.5, rs, dpsi) {
    //     println!("{}, {}", psi, r);
    // }
    
    // Camera's 3D location
    let camera_pos = Point3::new(6.0, 6.0, 4.0);
    // Point in 3D space the camera is pointing at
    let camera_focus = Point3::new(0.0, 0.0, 1.0);
    
    // Image dimensions
    let img_width = 256 * 8;
    let img_height = 256 * 8;
    // Width of the camera "lens"
    let camera_width = 2.0;
    // Height of the camera "lens"
    let camera_height: f64 = (img_height as f64) / (img_width as f64) * camera_width;
    // Distance between camera "lens" and location
    let camera_ap = 1.5;
    
    // Camera's facing direction
    let camera_direction = Unit::new_normalize(camera_focus - camera_pos);
    // Camera's direction in the XY plane
    let flat_camera = Vector3::new(camera_direction.x, camera_direction.y, 0.0);
    // The camera's relative left and up vectors
    let camera_x_axis = Rotation3::from_axis_angle(&Vector3::z_axis(), -f64::consts::PI / 2.0) * flat_camera;
    let camera_y_axis = Rotation3::from_axis_angle(&Unit::new_normalize(camera_x_axis), f64::consts::PI / 2.0) * *camera_direction;
    
    // Height of the checkered floor in the z-direction
    let floor_height = -2.0;
    // dpsi: angular distance per step
    let march_len = 0.001f64;
    
    // r_s of the black hole
    let bh_r = 0.3;

    let img = ImageBuffer::from_fn(img_width, img_height, |x, y| {
        let y = img_height - y;
        // Calculate the point on the camera "lens"
        // corresponding to this pixel's ray
        let ap_point: Point3<f64> = camera_pos + camera_ap * camera_direction.as_ref() + 
                                 ((x as f64) - (img_width as f64) / 2.0) / (img_width as f64) * camera_width * camera_x_axis +
                                 ((y as f64) - (img_height as f64) / 2.0) / (img_height as f64) * camera_height * camera_y_axis;
        // The travelling direction of the photon
        // is given by the vector lenspoint - camerapoint
        let init_dir = ap_point - camera_pos;
        
        // For the normal plane (no black hole):
        // let mut curr_point = camera_pos.clone();
        // loop {
        //     if curr_point.z <= floor_height {
        //         if (curr_point.x.ceil() as i32 + curr_point.y.ceil() as i32).rem_euclid(2) == 0 {
        //             return image::Rgb([0xF0u8, 0x64u8, 0x49u8]);
        //         } else {
        //             return image::Rgb([0xEDu8, 0xE6u8, 0xE3u8]);
        //         }
        //     }
        //     if curr_point.x.abs() > 10.0 || curr_point.y.abs() > 10.0 || curr_point.z.abs() > 10.0 {
        //         return image::Rgb([0x36u8, 0x38u8, 0x2Eu8]);
        //     }
        //     curr_point += init_dir * march_len;
        // }

        // A magic rotation that maps the plane of the
        // orbit of our photon into the YZ plane
        let transform = Rotation3::look_at_rh(
            &(camera_pos - Point3::origin()),
            &init_dir,
        );
        
        // We start at the camera position
        let init_pos = transform * camera_pos;
        let init_pos = cart_to_rad(init_pos.y, init_pos.z);
        // Alpha is the angle of emission for use in geosim
        let alpha = (transform * camera_pos - Point3::origin()).angle(
            &(transform * init_dir));
        
        // Number of iterations (for limiting super long loops)
        let mut itr = 0;
        
        // Debug: Draw a dot for tracing a single photon's orbit
        // if (x as i64 - 55).abs() < 3 && (y as i64 - 128).abs() < 3 {
            // return image::Rgb([0x00u8, 0xFFu8, 0x00u8]);
            // dbg!(init_pos.1);
        // }
        
        for (psi, r) in geosim(init_pos.1, init_pos.0, alpha, bh_r, march_len) {
            // Our particle is inside the event horizon
            if r <= bh_r {
                return image::Rgb([0u8, 0u8, 0u8]);
            }
            // Debug: trace a single photons orbit
            // for graphing in sage
            // if x == 55 && y == 128 {
            //     println!("{}, {}", psi, r);
            // }
            let (x, y) = rad_to_cart(psi, r);
            let pos = Vector3::new(0.0, x, y);
            // Undo the transformation (get the proper
            // 3D position of the photon)
            let curr_point = transform.inverse() * pos;

            // Check for floor collision
            if curr_point.z <= floor_height {
                if (curr_point.x.ceil() as i32 + curr_point.y.ceil() as i32).rem_euclid(2) == 0 {
                    return image::Rgb([0xF0u8, 0x64u8, 0x49u8]);
                } else {
                    return image::Rgb([0xEDu8, 0xE6u8, 0xE3u8]);
                }
            }
            // Check for sky collision
            if curr_point.x.abs() > 20.0 || curr_point.y.abs() > 20.0 || curr_point.z.abs() > 30.0 {
                return image::Rgb([0x36u8, 0x38u8, 0x2Eu8]);
            }
            // Check for orbits which take way
            // too long (useful for debugging)
            if itr > 50000 {
                return image::Rgb([0u8, 0u8, 0u8]);
            }
            itr += 1;
        }
        // Shouldn't normally happen, call it
        // sky collision
        image::Rgb([0x36u8, 0x38u8, 0x2Eu8])
    });

    img.save("output.png").unwrap();
}
