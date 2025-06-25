use crate::environment::K;
use crate::fit::polynomial;
use crate::solvent::Solvent;
use crate::water::{density_water, water, water_viscosity};
//Transient cavity dynamics and divergence from the Stokes–Einstein equation in organic aerosol
fn sucrose_water_diffusion(activity:f64) ->f64{
    10.0f64.powf(-15.613 + 1.262*activity - 3.476*activity.powi(2) + 46.468*activity.powi(3) - 60.030*activity.powi(4) + 22.691*activity.powi(5))
}

fn sucrose_water_diffusion_zobrist(activity:f64, temperature:f64) ->f64{
    let a = 7.0 + 0.175*(1.0-46.46*(1.0-activity));
    let b = 262.867*(1.0+10.53*(1.0-activity)-0.3*(1.0-activity).powi(2));
    let T_0 = 127.9*(1.0+0.4514*(1.0-activity)-0.5*(1.0-activity).powf(1.7));
    10.0f64.powf(-(a+b/(temperature-T_0)))
}

fn sucrose_vad_mfs(concentration:f64)->f64{
    if concentration <= 0.0{
        0.0
    } else {
        1.0/((1.0/concentration - 1.0/1580.5)*1000.0+1.0)
    }
}

fn sucrose_activity(mfs:f64,temperature:f64)->f64{
    let a = -1.0;
    let b = -0.99721;
    let c = 0.13599;
    let d = 0.001688;
    let e = -0.005151;
    let f = 0.009607;
    let g = -0.006142;
    let T = 298.15;
    (1.0+a*mfs)/(1.0+b*mfs+c*mfs.powi(2)) + (temperature-T)*(d*mfs+e*mfs.powi(2)+f*mfs.powi(3)+g*mfs.powi(4))
}

fn sucrose_vad_density(mfs:f64)->f64{
    1.0/((1.0-mfs)/1000.0 + mfs/1580.5)
}

fn sucrose_viscosity(rh:f64)->f64{
    10.0f64.powf(15.92 - 0.276*rh + 8.68e4*rh.powi(2))
}
pub fn sucrose()->Binary{
    Binary{
        volatile_diffusion: |mfw,temperature| sucrose_water_diffusion(sucrose_activity(1.0-mfw,temperature)),
        non_volatile_diffusion: |mfs,temperature| sucrose_water_diffusion(sucrose_activity(mfs,temperature)),
        activity: |mfs,temperature|sucrose_activity(mfs,temperature),
        density:|mfs| sucrose_vad_density(mfs),
        viscosity: |mfs, temperature| sucrose_viscosity(sucrose_activity(mfs, temperature)*100.0),
        volatile: water(),
        mfs: sucrose_vad_mfs
    }
}

pub fn sucrose_zobrist()->Binary{
    Binary{
        volatile_diffusion: |mfw,temperature| sucrose_water_diffusion_zobrist(sucrose_activity(1.0-mfw,temperature),temperature),
        non_volatile_diffusion: |mfs,temperature| sucrose_water_diffusion_zobrist(sucrose_activity(mfs,temperature),temperature),
        activity: |mfs,temperature|sucrose_activity(mfs,temperature),
        density:|mfs| sucrose_vad_density(mfs),
        viscosity: |mfs, temperature| sucrose_viscosity(sucrose_activity(mfs, temperature)*100.0),
        volatile: water(),
        mfs: sucrose_vad_mfs
    }
}

pub struct Binary {
    pub volatile_diffusion: fn(f64,f64) -> f64,
    pub non_volatile_diffusion: fn(f64,f64) -> f64,
    pub activity: fn(f64,f64)->f64,
    pub density: fn(f64)->f64,
    pub viscosity: fn(f64,f64) -> f64,
    pub mfs: fn(f64)->f64,
    pub volatile: Solvent
}

impl Binary{
    pub fn concentration(&self, mfs: f64) -> f64 {
        (self.density)(mfs) * mfs
    }
    pub fn get(name:String)->Binary{
        match name.as_str(){
            "sucrose"=>sucrose(),
            "sucrose_zobrist"=>sucrose_zobrist(),
            "water"=>pure_water(),
            "nacl"=>nacl(),
            bad_solution => {panic!("{} IS NOT A KNOWN BINARY SOLUTION",bad_solution)},
        }
    }
}

fn aqueous_NaCl_diffusion(mfs:f64 , temperature:f64)->f64{
    let D_0=1e-9*(1.955 - 20.42*mfs + 141.7*mfs.powi(2) - 539.8*mfs.powi(3) + 995.6*mfs.powi(4) - 698.7*mfs.powi(5));
    D_0*temperature/293.0*water_viscosity(293.0)/water_viscosity(temperature)
}

fn nacl()->Binary{
    Binary{
        volatile_diffusion: |mfw,temperature| aqueous_NaCl_diffusion(1.0-mfw,temperature),
        non_volatile_diffusion: |mfs,temperature| aqueous_NaCl_diffusion(mfs,temperature),
        activity: |mfs,temperature| polynomial(&[48.5226539, -158.04388699, 186.59427048, -93.88696437, 19.28939256,
            -2.99894206, -0.47652352, 1.],mfs),
        density: |mfs|polynomial(&[-940.62808,2895.88613, -2131.05669, 1326.69542 , -55.33776, 998.2],mfs),
        viscosity: |mfs, temperature| water_viscosity(temperature)*aqueous_NaCl_diffusion(0.0, temperature)/aqueous_NaCl_diffusion(mfs, temperature),
        mfs: |concentration| polynomial(&[-3.67338330e-21,  3.36689881e-17, -1.37012771e-13,  3.36008061e-10,
        -5.85656285e-07,  9.89047989e-04,  2.45656466e-04],concentration),
        volatile: water(),
    }
}

fn pure_water() -> Binary {
    Binary{
        volatile_diffusion: |mfw,temperature| 2e-9,
        non_volatile_diffusion: |mfs,temperature| 2e-9,
        activity: |mfs,temperature| 1.0,
        density: |mfs|density_water(293.0),
        viscosity: |mfs,temperature|water_viscosity(temperature),
        mfs: |concentration| 0.0,
        volatile: water(),
    }
}