/* Original implementation:
 *   https://github.com/Signalsmith-Audio/hilbert-iir
 *   Signalsmith's Hilbert IIR: A single-file, dependency-free Hilbert filter.
 *   Copyright (c) 2024 Geraint Luff / Signalsmith Audio Ltd.
 *   Released under 0BSD: BSD Zero Clause License
 */

use nalgebra::Complex;

// Approximate Hilbert analytic transform with an IIR filter.
// The filter uses a 12th-order IIR design with pre-calculated coefficients and
// poles. To increase throughput, this translation is designed to work with
// multiple channels simultaneously.

const ORDER: usize = 12;

/// Pre-calculated complex coefficients for the IIR filter.
const COEFFS: [Complex<f32>; ORDER] = [
    Complex::<f32>::new(-0.000224352093802, 0.00543499018201),
    Complex::<f32>::new(0.010750055781500, -0.01738906856810),
    Complex::<f32>::new(-0.045679587391700, 0.02291669314290),
    Complex::<f32>::new(0.112825005820000, 0.00278413661237),
    Complex::<f32>::new(-0.208067578452000, -0.10462895867500),
    Complex::<f32>::new(0.287178375010000, 0.33619239719000),
    Complex::<f32>::new(-0.254675294431000, -0.68303389965500),
    Complex::<f32>::new(0.048108183502600, 0.95406158937400),
    Complex::<f32>::new(0.227861357867000, -0.89127357456900),
    Complex::<f32>::new(-0.365411839137000, 0.52508831727100),
    Complex::<f32>::new(0.280729061131000, -0.15513120660600),
    Complex::<f32>::new(-0.093506178772800, 0.00512245855404),
];

/// Pre-calculated complex poles for the IIR filter.
const POLES: [Complex<f32>; ORDER] = [
    Complex::<f32>::new(-0.00495335976478, 0.0092579876872),
    Complex::<f32>::new(-0.01785949130200, 0.0273493725543),
    Complex::<f32>::new(-0.04137143731550, 0.0744756910287),
    Complex::<f32>::new(-0.08821484088850, 0.1783496774570),
    Complex::<f32>::new(-0.17922965812000, 0.3960134022300),
    Complex::<f32>::new(-0.33826180075300, 0.8292295333540),
    Complex::<f32>::new(-0.55768869973200, 1.6129853832800),
    Complex::<f32>::new(-0.73515773614800, 2.7998739868200),
    Complex::<f32>::new(-0.71905738117200, 4.1639616612800),
    Complex::<f32>::new(-0.51787102520900, 5.2972482680400),
    Complex::<f32>::new(-0.28019746947100, 5.9959860238800),
    Complex::<f32>::new(-0.08527513545310, 6.3048492377000),
];

/// Direct term for the IIR filter.
const DIRECT: f32 = 0.000262057212648;

// Performance note:
// Internal structures separate the real and imaginary components using two
// float arrays rather than an array of Complex<f32> as this offer more
// opportunity for auto-vectorization.

/// Structure representing the IIR filter with its coefficients and poles.
/// The coefficients are computed once and shared across all channels.
#[derive(Clone, Copy)]
pub struct Filter {
    coeffs_r: [f32; ORDER],
    coeffs_i: [f32; ORDER],
    poles_r: [f32; ORDER],
    poles_i: [f32; ORDER],
    direct: f32,
}

/// Structure representing the internal state of the filter for one channel.
#[derive(Clone, Copy)]
pub struct State {
    real: [f32; ORDER],
    imag: [f32; ORDER],
}

impl State {
    /// Initial state for a channel.
    pub const INITIAL: State = State {
        real: [-1.0; ORDER],
        imag: [0.0; ORDER],
    };
}

impl Filter {
    /// Initialize a filter with the given sample rate and passband gain.
    /// Recommended settings: `Filter::new(48000.0, 2.0)`
    pub fn new(sample_rate: f32, passband_gain: f32) -> Filter {
        let freq_factor = f32::min(0.46, 20000.0 / sample_rate);
        let mut result = Filter {
            coeffs_r: [0.0; ORDER],
            coeffs_i: [0.0; ORDER],
            poles_r: [0.0; ORDER],
            poles_i: [0.0; ORDER],
            direct: 0.0,
        };
        result.direct = DIRECT * 2.0 * passband_gain * freq_factor;
        for i in 0..ORDER {
            let coeff = COEFFS[i] * freq_factor * passband_gain;
            result.coeffs_r[i] = coeff.re;
            result.coeffs_i[i] = coeff.im;
            let pole = (POLES[i] * freq_factor).exp();
            result.poles_r[i] = pole.re;
            result.poles_i[i] = pole.im;
        }
        result
    }

    /// Process the input signal and produce the output signal for N channels.
    pub fn process<const N: usize>(
        &self,
        state: &mut [State; N],
        input: &[&[f32]; N],
        output: &mut [&mut [Complex<f32>]; N],
    ) {
        if N == 0 {
            return;
        };
        let samples = input[0].len();
        for i in 0..samples {
            let mut ra: [f32; N] = input.map(|x| x[i] * self.direct);
            let mut ia: [f32; N] = [0f32; N];
            for j in 0..ORDER {
                let mut rv = [0f32; N];
                let mut iv = [0f32; N];
                for k in 0..N {
                    rv[k] = state[k].real[j] * self.poles_r[j]
                        - state[k].imag[j] * self.poles_i[j]
                        + input[k][i] * self.coeffs_r[j];
                    iv[k] = state[k].real[j] * self.poles_i[j]
                        + state[k].imag[j] * self.poles_r[j]
                        + input[k][i] * self.coeffs_i[j];
                }
                for k in 0..N {
                    ra[k] += rv[k];
                    ia[k] += iv[k];
                    state[k].real[j] = rv[k];
                    state[k].imag[j] = iv[k];
                }
            }
            for k in 0..N {
                output[k][i] = Complex::new(ra[k], ia[k]);
            }
        }
    }

    /// Process the input signal and produce separate real and imaginary output
    /// signals for N channels.
    #[allow(dead_code)]
    pub fn process_split<const N: usize>(
        &self,
        state: &mut [State; N],
        input: &[&[f32]; N],
        real: &mut [&mut [f32]; N],
        imag: &mut [&mut [f32]; N],
    ) {
        if N == 0 {
            return;
        };
        let samples = input[0].len();
        for i in 0..samples {
            let mut ra: [f32; N] = input.map(|x| x[i] * self.direct);
            let mut ia: [f32; N] = [0f32; N];
            for j in 0..ORDER {
                for k in 0..N {
                    let r = state[k].real[j] * self.poles_r[j]
                        - state[k].imag[j] * self.poles_i[j]
                        + input[k][i] * self.coeffs_r[j];
                    let i = state[k].real[j] * self.poles_i[j]
                        + state[k].imag[j] * self.poles_r[j]
                        + input[k][i] * self.coeffs_i[j];
                    ra[k] += r;
                    ia[k] += i;
                    state[k].real[j] = r;
                    state[k].imag[j] = i;
                }
            }
            for k in 0..N {
                real[k][i] = ra[k];
                imag[k][i] = ia[k];
            }
        }
    }
}
