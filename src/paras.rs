/// Parameters file:
/// This contains all values required to instantiate the simulation, these are:
///     - Upper and lower state bounds
///     - Upper and lower input bounds
///     - Upper and lower input rate bounds
///     - Initial conditions
///     - time step
///     - simulation time

use ndarray::prelude::*;
use ndarray::ArrayBase;

pub struct Lim {
    pub x_ub: Vec<f32>,
    pub x_lb: Vec<f32>,
    pub u_ub: Vec<f32>,
    pub u_lb: Vec<f32>,
    pub udot_ub: Vec<f32>,
    pub udot_lb: Vec<f32>,
}

pub fn lim() -> Lim {

    let m2f: f32 = 3.28084;
    let f2m: f32 = 1.0/m2f;

    //----------------------Limits---------------------------//
    let npos_min: f32        = -f32::INFINITY;       // (m)
    let epos_min: f32        = -f32::INFINITY;       // (m)
    let h_min: f32           = 0.0;             // (m)
    let phi_min: f32         = -f32::INFINITY;       // (deg)
    let theta_min: f32       = -f32::INFINITY;       // (deg)
    let psi_min: f32         = -f32::INFINITY;       // (deg)
    let vt_min: f32           = 0.0;             // (m/s)
    let alpha_min: f32       = -20.;          // (deg)
    let beta_min: f32        = -30.;          // (deg)
    let p_min: f32           = -30.0;           // (deg/s)
    let q_min: f32           = -10.0;           // (deg/s)
    let r_min: f32           = -5.0;            // (deg/s)

    let thrust_min: f32           = 1000.0;          // (lbs)
    let dh_min: f32          = -25.0;           // (deg)
    let da_min: f32          = -21.5;         // (deg)
    let dr_min: f32          = -30.0;          // (deg)
    let lef_min: f32         = 0.0;            // (deg)

    let npos_max: f32        = f32::INFINITY;        // (m)
    let epos_max: f32        = f32::INFINITY;        // (m)
    let h_max: f32           = 100000.0;         // (m)
    let phi_max: f32         = f32::INFINITY;        // (deg)
    let theta_max: f32       = f32::INFINITY;        // (deg)
    let psi_max: f32         = f32::INFINITY;        // (deg)
    let vt_max: f32           = 900.0;           // (m/s)
    let alpha_max: f32       = 90.0;            // (deg)
    let beta_max: f32        = 30.0;            // (deg)
    let p_max: f32           = 30.0;            // (deg/s)
    let q_max: f32           = 10.0;            // (deg/s)
    let r_max: f32           = 5.0;             // (deg/s)

    let thrust_max: f32           = 19000.0;         // (lbs)
    let dh_max: f32          = 25.0;            // (deg)
    let da_max: f32          = 21.5;          // (deg)
    let dr_max: f32          = 30.0;            // (deg)
    let lef_max: f32         = 25.0;            // (deg)

    let x_ub = vec![npos_max, epos_max, h_max, phi_max, theta_max, psi_max, vt_max, alpha_max, beta_max, p_max, q_max, r_max, thrust_max, dh_max, da_max, dr_max, lef_max, f32::INFINITY];
    let x_lb = vec![npos_min, epos_min, h_min, phi_min, theta_min, psi_min, vt_min, alpha_min, beta_min, p_min, q_min, r_min, thrust_min, dh_min, da_min, dr_min, lef_min, -f32::INFINITY];
    
    let u_ub = vec![thrust_max, dh_max, da_max, dr_max];
    let u_lb = vec![thrust_min, dh_min, da_min, dr_min];

    let udot_ub = vec![10000.0, 60.0, 80.0, 120.0];
    let udot_lb = vec![-10000.0, -60.0, -80.0, -120.0];

    let lim = Lim{
        x_ub: x_ub,
        x_lb: x_lb,
        u_ub: u_ub,
        u_lb: u_lb,
        udot_ub: udot_ub,
        udot_lb: udot_lb,
    };

    return lim
}

pub fn ic() -> ndarray::Array {

    let pi: f32 = 3.1415926536897932384626;

    //-----------------Initial Conditions--------------------//
    let npos: f32        = 0.;                // m
    let epos: f32        = 0.;                // m
    let h: f32           = 3048.;             // m
    let phi: f32         = 0.;                // rad
    let theta: f32       = 0.;                // rad
    let psi: f32         = 0.;                // rad

    let vt: f32          = 213.36;            // m/s
    let alpha: f32       = 1.0721 * pi/180.0; // rad
    let beta: f32        = 0.;                // rad
    let p: f32           = 0.;                // rad/s
    let q: f32           = 0.;                // rad/s
    let r: f32           = 0.;                // rad/s

    let thrust: f32      = 2886.6468;         // lbs
    let dh: f32          = -2.0385;           // deg
    let da: f32          = -0.087577;         // deg
    let dr: f32          = -0.03877;          // deg
    let lef: f32         = 0.3986;            // deg

    let x0 = array![npos, epos, h, phi, theta, psi, vt, alpha, beta, p, q, r, thrust, dh, da, dr, lef].t();

    return x0
}

pub fn sim_paras() -> (f32, f32) {
    let dt: f32 = 0.001;
    let sim_time: f32 = 10.0;
    return (dt, sim_time)
}