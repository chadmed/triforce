// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * An attempt at an MVDR beamformer for equilateral triangle mic arrays
 * as found on newer Apple Silicon Macs.
 *
 * Currently mono, but could probably be extended to be stereo somewhat
 * easily. I think.
 *
 * Copyright (C) 2024 James Calligeros <jcalligeros99@gmail.com>
 */

use std::f32::consts::PI;

use lv2::prelude::*;

use nalgebra::{linalg::SVD, Matrix3, Vector3};

use rustfft::{num_complex::Complex, num_traits::Zero, FftPlanner};

use std::sync::Mutex;

const C: f32 = 343.00; /* m*s^-1 */

/// The distance of a given element in the array from the zeroth
/// element
#[derive(Copy, Clone)]
struct ElemDistance {
    x: f32,
    y: f32,
}

/// Perform a Hilbert transform on a slice of f32s to give us the analytic signal of
/// our input sample buffer. This is necessary to extract phase information from the
/// signal, and to make matrix operations a bit easier.
fn analytic_signal(planner: &mut FftPlanner<f32>, signal: &[f32]) -> Vec<Complex<f32>> {
    let len: usize = signal.len();

    // Convert each real sample into a complex sample
    let mut complex_signal: Vec<Complex<f32>> =
        signal.iter().map(|&x| Complex::new(x, 0.0)).collect();

    // Set up the fft and inverse fft
    let fft = planner.plan_fft_forward(len);
    let ifft = planner.plan_fft_inverse(len);

    // Mutate the output buffer into the forward FFT
    fft.process(&mut complex_signal);

    // Perform the Hilbert transform on the FFT. To do this, we multiply every
    // positive sample under the Nyquist limit by 2+0j, and destroy every sample
    // above it.
    for i in 0..len {
        if i > 0 && i < len / 2 {
            complex_signal[i] *= Complex::new(2.0, 0.0);
        } else if i >= len / 2 {
            complex_signal[i] = Complex::zero();
        }
    }

    // Turn the original complex buffer into the inverse FFT and then normalise
    ifft.process(&mut complex_signal);
    complex_signal.iter_mut().for_each(|x| *x /= len as f32);

    complex_signal
}

/// The steering vector is a representation of the phase delays at each microphone.
/// It is calculated by taking the dot product of the array geometry matrix and the
/// unit vector of the direction of arrival.
fn steering_vec(theta: f32, phi: f32, f: f32, elems: [ElemDistance; 3]) -> Vector3<Complex<f32>> {
    // Mic positions are relative to Left/Top to preserve x/y axis semantics
    let mic_positions: Vec<Vector3<f32>> =
        elems.iter().map(|e| Vector3::new(e.x, e.y, 0f32)).collect();

    // Calculate angular repetency (2pi/lambda)
    let repetency = (2f32 * PI) / (f / C);

    // Compute the unit vector of the DOA
    let u_dir = Vector3::new(phi.sin() * theta.cos(), phi.sin() * theta.sin(), phi.cos());

    // Calculate the steering vector by taking the array geometry, speed of sound,
    // and DOA unit vector
    let mut steering_vector = Vector3::from_element(Complex::new(0f32, 0f32));

    for (i, mic_pos) in mic_positions.iter().enumerate() {
        let delay = mic_pos.dot(&u_dir) / C;
        let phase = -repetency * delay;
        steering_vector[i] = Complex::new(phase.cos(), phase.sin());
    }

    steering_vector
}

/// There's nothing special about this, it's just a covariance matrix. It is always
/// square.
fn covariance(signals: &Vec<Vec<Complex<f32>>>) -> Matrix3<Complex<f32>> {
    let n_samples = signals[0].len();

    let mut covar = Matrix3::zeros();

    for t in 0..n_samples {
        let discrete: Vector3<Complex<f32>> =
            Vector3::from_iterator(signals.iter().map(|s| s[t]));
        covar += &discrete * discrete.adjoint();
    }

    // Our samples are shit, so we can't get a very nice covariance matrix.
    // Regularise the shit covariance matrix by introducing a constant value
    // across the identity
    let reg = Matrix3::identity()
        .map(|x: f32| Complex::new(x * 1e-4f32, 0f32));

    covar /= Complex::new(n_samples as f32, 0f32);
    covar + reg
}

/// To calculate the weighting vector for the beamformer, we need to invert
/// the covariance matrix, multiply it by the steering vector, then divide that
/// by itself multiplied by the Hermitian transpose of the steering vector
/// w = (cov^-1 * sv) / (sv.adjoint() * (cov^-1 * sv)). Note that the denominator
/// is the same as the conjugate-linear dot product of the steering vector and
/// the numerator.
fn mvdr_weights(cov: &Matrix3<Complex<f32>>, sv: &Vector3<Complex<f32>>) -> Vector3<Complex<f32>> {
    // Since we have a numerically unstable covariance matrix, we can't take the
    // true inverse of it. Let's instead decompose it and take the pseudoinverse.
    let svd = SVD::new(cov.to_owned(), true, true);
    let r_inv = svd.pseudo_inverse(1e-4f32).unwrap();

    let num = r_inv * sv;
    let den = sv.dotc(&num); // Conjugate-linear dot product

    num / den
}

/*
 * Input and output ports used by the plugin
 *
 *
 * Ports:
 *      in_1: channel 1 input (left/top)
 *      in_2: channel 2 input (right/bottom)
 *      in_3: channel 3 input (vertex)
 *      out: output
 *      h_angle: horizontal steering angle in degrees (relative to input 1)
 *      v_angle: vertical steering angle in degrees
 *      opt_freq: frequency to optimise for
 *      t_win: covariance matrix time window
 */
#[derive(PortCollection)]
pub struct Ports {
    in_1: InputPort<Audio>,
    in_2: InputPort<Audio>,
    in_3: InputPort<Audio>,
    out: OutputPort<Audio>,
    h_angle: InputPort<Control>,
    v_angle: InputPort<Control>,
    opt_freq: InputPort<Control>,
    t_win: InputPort<Control>,
    mic2_x: InputPort<Control>,
    mic2_y: InputPort<Control>,
    mic3_x: InputPort<Control>,
    mic3_y: InputPort<Control>,
}

/*
 * Plugin state
 */
#[uri("https://chadmed.au/triforce")]
pub struct Triforce {
    hangle_curr: f32,
    vangle_curr: f32,
    freq_curr: f32,
    sample_rate: f32,
    samples_since_last_update: u32,
    covar_window: Vec<Vec<Complex<f32>>>,
    steering_vector: Vector3<Complex<f32>>,
    covar: Matrix3<Complex<f32>>,
    array_geom: [ElemDistance; 3],
    fft_planner: Mutex<FftPlanner<f32>>,
    weights: Vector3<Complex<f32>>,
}

trait Beamformer: Plugin {
    fn update_params(&mut self, ports: &mut Ports);
}

impl Triforce {

    pub fn with_sample_rate(sample_rate: f32) -> Self {
        Self {
            hangle_curr: 0f32,
            vangle_curr: 0f32,
            freq_curr: 1000f32,
            samples_since_last_update: u32::max_value(),
            sample_rate,
            covar_window: vec![
                vec![Complex::new(0f32, 0f32); 256],
                vec![Complex::new(0f32, 0f32); 256],
                vec![Complex::new(0f32, 0f32); 256],
            ],
            array_geom: [ElemDistance { x: 0f32, y: 0f32 }; 3],
            steering_vector: steering_vec(
                90f32.to_radians(),
                45f32.to_radians(),
                1000f32,
                [ElemDistance { x: 0f32, y: 0f32 }; 3],
            ),
            covar: Matrix3::zeros(),
            fft_planner: Mutex::new(FftPlanner::new()),
            weights: Vector3::zeros(),
        }
    }

    pub fn process_slice(&mut self, mic1: &[f32], mic2: &[f32], mic3: &[f32], output: &mut [f32], t_win: f32) {

        // All three sample buffers will have the same number of samples
        let num_samples = mic1.len();

        // Steering vector is relative to Left/Top mic
        let inputs = {
            let mut planner = self.fft_planner.lock().unwrap();
            vec![
                analytic_signal(&mut planner, mic1),
                analytic_signal(&mut planner, mic2),
                analytic_signal(&mut planner, mic3),
            ]
        };

        // Update the covariance matrix. We use an overlapping window to smooth over
        // the transitions.
        if self.samples_since_last_update as f32 >= (t_win / 1000f32) * self.sample_rate {
            self.samples_since_last_update = 0;
            // We want a 1/3 overlap
            let i = num_samples / 3;
            self.covar_window[0].extend_from_slice(&inputs[0][0..i]);
            self.covar_window[1].extend_from_slice(&inputs[1][0..i]);
            self.covar_window[2].extend_from_slice(&inputs[2][0..i]);
            self.covar = covariance(&self.covar_window);
            self.covar_window[0] = inputs[0][i + 1..num_samples].to_vec();
            self.covar_window[1] = inputs[1][i + 1..num_samples].to_vec();
            self.covar_window[2] = inputs[2][i + 1..num_samples].to_vec();
            self.weights = mvdr_weights(&self.covar, &self.steering_vector);
        }
        else {
            self.samples_since_last_update += num_samples as u32;
        }

        for t in 0..num_samples {
            let discrete: Vector3<Complex<f32>> = Vector3::from_iterator(
                inputs.iter().map(|s| s[t]),
            );

            let out =
                // Conjugate-linear dot product
                self.weights.dotc(&discrete)
                // // Now we need to revert the Hilbert transform and output the signal
                .re;

            // Do all of our NFP and clamping here
            output[t] =
                if out.is_finite() && !out.is_nan() {
                    out.clamp(-10f32, 10f32)
                } else {
                    0f32
                };
        }
    }
}


impl Plugin for Triforce {
    type Ports = Ports;

    type InitFeatures = ();
    type AudioFeatures = ();

    fn new(info: &PluginInfo, _features: &mut ()) -> Option<Self> {
        Some(Self::with_sample_rate(info.sample_rate() as f32))
    }

    fn run(&mut self, ports: &mut Ports, _features: &mut (), samples: u32) {
        Beamformer::update_params(self, ports);
        if samples < 1024 {
            return;
        }
        self.process_slice(&ports.in_1, &ports.in_2, &ports.in_3, &mut ports.out, *ports.t_win);
    }
}

impl Beamformer for Triforce {
    fn update_params(&mut self, ports: &mut Ports) {
        if self.hangle_curr != *ports.h_angle
            || self.freq_curr != *ports.opt_freq
            || self.vangle_curr != *ports.v_angle
            || self.array_geom[1].x != *ports.mic2_x
            || self.array_geom[1].y != *ports.mic2_y
            || self.array_geom[2].x != *ports.mic3_x
            || self.array_geom[2].y != *ports.mic3_y
        {
            self.hangle_curr = *ports.h_angle;
            self.vangle_curr = *ports.v_angle;
            self.freq_curr = *ports.opt_freq;
            self.array_geom = [
                ElemDistance { x: 0f32, y: 0f32 },
                ElemDistance {
                    x: *ports.mic2_x,
                    y: *ports.mic2_y,
                },
                ElemDistance {
                    x: *ports.mic3_x,
                    y: *ports.mic3_y,
                },
            ];
            self.steering_vector = steering_vec(
                self.hangle_curr.to_radians(),
                self.vangle_curr.to_radians(),
                self.freq_curr,
                self.array_geom,
            );

            // The steering vector has changed
            self.weights = mvdr_weights(&self.covar, &self.steering_vector);
        }
    }
}

lv2_descriptors!(Triforce);
