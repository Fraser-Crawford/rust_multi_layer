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

fn sucrose_density(mfs:f64)->f64{
    1000.0*(0.9989+0.3615*mfs+0.2964*mfs.powi(2)-0.3186*mfs.powi(3)+0.24191*mfs.powi(4))
}

fn sucrose_viscosity(rh:f64)->f64{
    10.0f64.powf(15.92 - 0.276*rh + 8.68e4*rh.powi(2))
}
pub fn sucrose()->Binary{
    Binary{
        volatile_diffusion: |mfw,temperature| sucrose_water_diffusion(sucrose_activity(1.0-mfw,temperature)),
        non_volatile_diffusion: |mfs,temperature| sucrose_water_diffusion(sucrose_activity(mfs,temperature)),
        activity: |mfs,temperature|sucrose_activity(mfs,temperature),
        density:|mfs| sucrose_density(mfs),
        viscosity: |mfs, temperature| sucrose_viscosity(sucrose_activity(mfs, temperature)*100.0),
        volatile: water(),
        mfs: |concentration| polynomial(&[ 4.72908344e-11, -2.88752522e-07,  9.70951175e-04,  2.13331134e-03],concentration),
        solute_vapour_pressure: |mfs, temperature| 0.0,
        solute_latent_heat: |temperature| 0.0,
        solute_molar_mass: 342.3,
        solute_vapour_binary_diffusion_coefficient: |temperature|0.0,
    }
}

pub fn sucrose_zobrist()->Binary{
    Binary{
        volatile_diffusion: |mfw,temperature| sucrose_water_diffusion_zobrist(sucrose_activity(1.0-mfw,temperature),temperature),
        non_volatile_diffusion: |mfs,temperature| sucrose_water_diffusion_zobrist(sucrose_activity(mfs,temperature),temperature),
        activity: |mfs,temperature|sucrose_activity(mfs,temperature),
        density:|mfs| sucrose_density(mfs),
        viscosity: |mfs, temperature| sucrose_viscosity(sucrose_activity(mfs, temperature)*100.0),
        volatile: water(),
        mfs: |concentration| polynomial(&[ 4.72908344e-11, -2.88752522e-07,  9.70951175e-04,  2.13331134e-03],concentration),
        solute_vapour_pressure: |mfs, temperature| 0.0,
        solute_latent_heat: |temperature| 0.0,
        solute_molar_mass: 342.3,
        solute_vapour_binary_diffusion_coefficient: |temperature|0.0,
    }
}

pub struct Binary {
    pub volatile_diffusion: fn(f64,f64) -> f64,
    pub non_volatile_diffusion: fn(f64,f64) -> f64,
    pub activity: fn(f64,f64)->f64,
    pub density: fn(f64)->f64,
    pub viscosity: fn(f64,f64) -> f64,
    pub mfs: fn(f64)->f64,
    pub volatile: Solvent,
    pub solute_vapour_pressure: fn(f64,f64)->f64,
    pub solute_latent_heat: fn(f64)->f64,
    pub solute_molar_mass: f64,
    pub solute_vapour_binary_diffusion_coefficient: fn(f64)->f64,
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
            "acetone_water"=>acetone_water(),
            "ethanol_water"=>ethanol_water(),
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
        density: |mfs|polynomial(&[-940.62808,2895.88613, -2131.05669, 1326.69542 , -55.33776, 998.2],mfs.sqrt()),
        viscosity: |mfs, temperature| water_viscosity(temperature)*aqueous_NaCl_diffusion(0.0, temperature)/aqueous_NaCl_diffusion(mfs, temperature),
        mfs: |concentration| polynomial(&[-3.67338330e-21,  3.36689881e-17, -1.37012771e-13,  3.36008061e-10,
        -5.85656285e-07,  9.89047989e-04,  2.45656466e-04],concentration),
        volatile: water(),
        solute_vapour_pressure: |mfs, temperature| 0.0,
        solute_latent_heat: |temperature| 0.0,
        solute_molar_mass: 58.443,
        solute_vapour_binary_diffusion_coefficient: |temperature|0.0,
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
        solute_vapour_pressure: |mfs, temperature| 0.0,
        solute_latent_heat: |temperature| 0.0,
        solute_molar_mass: 0.0,
        solute_vapour_binary_diffusion_coefficient: |temperature|0.0,
    }
}

fn acetone_water() -> Binary {
    Binary{
        volatile_diffusion: |mfw,temperature| 2e-9,
        non_volatile_diffusion: |mfs,temperature| 2e-9,
        activity: |mfs,temperature| polynomial(&[ -4.20909806,  11.34113248, -10.86495888,   4.68627885,   0.02769658],1.0-mfs).sqrt(),
        density: |mfs| polynomial(&[-11.48319621,   67.8898522,  -268.83710607,  996.97012553],mfs),
        viscosity: |mfs,temperature| water_viscosity(temperature),
        mfs: |concentration| polynomial(&[ 1.26105862e-10,  2.38874563e-07,  1.01002766e-03, -3.77828436e-04],concentration),
        volatile: water(),
        solute_vapour_pressure: |mfs, temperature| 133.32*10.0f64.powf(7.02447-(1167.235/(temperature+224.844-273.15)))*
            polynomial(&[-1.63526684,  5.78856728, -6.33670548,  3.13542919,  0.01447529],mfs),
        solute_latent_heat: |temperature| 31e3/58.08e-3,
        solute_molar_mass: 58.08,
        solute_vapour_binary_diffusion_coefficient: |temperature| 1.22e-5,
    }
}

fn ethanol_water()-> Binary {
    Binary{
        volatile_diffusion:|mfw,temperature| 2e-9,
        non_volatile_diffusion:|mfs,temperature| 2e-9,
        activity: |mfs,temperature| polynomial(&[ 1.06320401, -0.96404147, -1.31940143,  2.27407497, -0.04457325],1.0-mfs).sqrt(),
        density: |mfs| polynomial(&[ -92.42424242, -108.39393939,  998.54545455],mfs),
        viscosity: |mfs,temperature| water_viscosity(temperature),
        mfs: |concentration| polynomial(&[ 5.47018806e-10, -1.95592318e-07,  1.05808057e-03, -1.24664184e-03],concentration),
        volatile: water(),
        solute_vapour_pressure: |mfs, temperature| 133.32*10.0f64.powf(8.04494-(1554.3/(temperature+222.65-273.15)))*
            polynomial(&[-0.59874137,  2.83710217, -3.01591469,  1.74666051,  0.01336109],mfs),
        solute_latent_heat: |temperature| 42.32e3/46.069e-3,
        solute_molar_mass: 46.069,
        solute_vapour_binary_diffusion_coefficient: |temperature| 1.27e-5,
    }
}