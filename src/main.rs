#![allow(unused_imports)]
#![allow(unused_variables)]
#![allow(dead_code)]

mod paras;

use ndarray::prelude::*;
use ndarray::Array;

use typenum::{U16, U1024};

fn main() {

    // set flags
    // fi_flag = 1 -> high fidelity model (full Nguyen)
    // fi_flag = 1 -> low fidelity model (Stevens Lewis reduced)
    let fi_flag: i32 = 1;

    // stability_flag only functional for high fidelity model currently!
    // stability_flag = 1 -> unstable xcg 35% model
    // stability_flag = 0 -> stable xcg 25% model
    let stab_flag: i32 = 0;

    // get limits from paras
    let lim = paras::lim();

    // print out some debugging
    // println!{"{:?}",lim.x_ub}
    // println!{"{:?}",lim.x_lb}

    // get simulation parameters from paras
    let (dt, sim_time) = paras::sim_paras();

    let a = array![
            [1.,2.,3.], 
            [4.,5.,6.],
        ]; 
    //assert_eq!(a.ndim(), 2);         // get the number of dimensions of array a
    //assert_eq!(a.len(), 6);          // get the number of elements in array a
    //assert_eq!(a.shape(), [2, 3]);   // get the shape of array a
    //assert_eq!(a.is_empty(), false); // check if the array has zero elements

    //let a = Array::<f64, _>::linspace(0.,5.,11);
    //let b = Array::range(0., 4., 1.);

    println!("{:?}", &a);

    println!("{:?}", a.dot(&a.t()));

}

