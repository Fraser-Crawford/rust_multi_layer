pub fn polynomial(coefficients: &[f64], x: f64) -> f64 {
    let n = coefficients.len();
    coefficients
        .iter()
        .enumerate()
        .fold(0.0, |acc, (i, &coefficient)| {
            acc + coefficient * x.powi((n - i - 1) as i32)
        })
}
pub fn old_convection_coefficients(radius:f64)->[f64; 4]{
    let radius_um = radius*1e6;
    let a = polynomial(&[ 2.37956833e+05, -8.80948549e+06,  1.27526684e+08, -8.02394969e+08,
                       1.80753529e+09,],radius_um);
    let b = polynomial(&[-3.24266198e-06,  2.52456713e-04, -8.64260775e-03,  2.26523227e-01,
                       1.76115671e+00],radius_um);
    let c = polynomial(&[-1.20462112e-05,  9.51099894e-04, -3.34584167e-02,  8.78979791e-01,
                       4.91237655e+00],radius_um);
    let d = polynomial(&[ 3.39923692e-02, -1.94358487e+00, -6.28015228e+01,  4.93658814e+02,
                       -2.45023768e+03],radius_um);
    [a,b,c,d]
}
pub fn convection_coefficients(radius:f64)->[f64; 3]{
    let radius_um = radius*1e6;
    let a = polynomial(&[-1.46116518e-02,  4.48801718e+01,  1.40847081e+02, -6.66552070e+02],radius_um);
    let b = polynomial(&[-3.91066347e-06,  2.94571781e-04, -8.29825186e-03,  1.50240218e-01],radius_um);
    let c = polynomial(&[ 6.39361187e-03,  8.52247064e+01, -6.32772673e+02,  2.84695376e+03],radius_um);
    [a,b,c]
}

pub fn reverse_planck(normalised_radius:f64,coefficients:[f64;4])->f64{
    let X = normalised_radius-1.0;
    if X == 0.0{
        0.0
    } else {
        coefficients[0]*(-X).powf(coefficients[1])/((-coefficients[2]*X).exp()-1.0)+coefficients[3]*X
    }
}

pub fn asymmetric_gaussian(x2:f64, coefficients:[f64;3]) ->f64{
    coefficients[0]/coefficients[1]*(-(x2-0.5).powi(2)/coefficients[1]).exp()+coefficients[2]*(1.0-x2)
}