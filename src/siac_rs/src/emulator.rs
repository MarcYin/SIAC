#![allow(non_local_definitions)]
//! Neural Network Emulator
//!
//! Fast forward pass for the two-hidden-layer neural network emulators
//! used to replace computationally expensive 6S radiative transfer calculations.

use ndarray::{s, Array1, Array2, Array3, ArrayView2};
use numpy::{IntoPyArray, PyArray2, PyArray3, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::prelude::*;

/// Two-hidden-layer neural network emulator
///
/// Implements the forward pass (and optional Jacobian) for the SIAC emulators
/// that approximate 6S radiative transfer model outputs.
#[pyclass]
pub struct TwoLayerNN {
    /// First layer weights (input_dim, hidden1)
    w1: Array2<f32>,
    /// First layer biases (hidden1,)
    b1: Array1<f32>,
    /// Second layer weights (hidden1, hidden2)
    w2: Array2<f32>,
    /// Second layer biases (hidden2,)
    b2: Array1<f32>,
    /// Output layer weights (hidden2, output_dim)
    w3: Array2<f32>,
    /// Output layer biases (output_dim,)
    b3: Array1<f32>,
}

#[pymethods]
impl TwoLayerNN {
    /// Create a new neural network emulator from weight arrays
    #[new]
    fn new(
        w1: PyReadonlyArray2<f32>,
        b1: PyReadonlyArray1<f32>,
        w2: PyReadonlyArray2<f32>,
        b2: PyReadonlyArray1<f32>,
        w3: PyReadonlyArray2<f32>,
        b3: PyReadonlyArray1<f32>,
    ) -> Self {
        Self {
            w1: w1.as_array().to_owned(),
            b1: b1.as_array().to_owned(),
            w2: w2.as_array().to_owned(),
            b2: b2.as_array().to_owned(),
            w3: w3.as_array().to_owned(),
            b3: b3.as_array().to_owned(),
        }
    }

    /// Forward pass with optional Jacobian computation
    ///
    /// # Arguments
    /// * `x` - Input array of shape (n_samples, n_inputs)
    /// * `compute_jacobian` - Whether to compute the Jacobian
    ///
    /// # Returns
    /// Tuple of (output, jacobian) where jacobian is None if not requested
    #[pyo3(signature = (x, compute_jacobian=false))]
    fn predict<'py>(
        &self,
        py: Python<'py>,
        x: PyReadonlyArray2<f32>,
        compute_jacobian: bool,
    ) -> PyResult<(&'py PyArray2<f32>, Option<&'py PyArray3<f32>>)> {
        let x = x.as_array();
        let (output, jacobian) = self.forward(x, compute_jacobian);

        let output_py = output.into_pyarray(py);
        let jacobian_py = jacobian.map(|j| j.into_pyarray(py));

        Ok((output_py, jacobian_py))
    }

    /// Get network architecture info
    fn get_architecture(&self) -> (usize, usize, usize, usize) {
        (
            self.w1.shape()[0], // input_dim
            self.w1.shape()[1], // hidden1
            self.w2.shape()[1], // hidden2
            self.w3.shape()[1], // output_dim
        )
    }
}

impl TwoLayerNN {
    /// Forward pass implementation
    fn forward(
        &self,
        x: ArrayView2<f32>,
        compute_jacobian: bool,
    ) -> (Array2<f32>, Option<Array3<f32>>) {
        let n_samples = x.shape()[0];

        // Layer 1: Linear + ReLU
        let a1 = x.dot(&self.w1) + &self.b1;
        let h1 = a1.mapv(|v| v.max(0.0));

        // Layer 2: Linear + ReLU
        let a2 = h1.dot(&self.w2) + &self.b2;
        let h2 = a2.mapv(|v| v.max(0.0));

        // Output layer: Linear (no activation)
        let output = h2.dot(&self.w3) + &self.b3;

        if compute_jacobian {
            // Backpropagation for Jacobian
            // dout/dx = w3^T @ diag(relu'(a2)) @ w2^T @ diag(relu'(a1)) @ w1^T
            let jacobian = self.compute_jacobian(&a1, &a2, n_samples);
            (output, Some(jacobian))
        } else {
            (output, None)
        }
    }

    /// Compute Jacobian via backpropagation
    fn compute_jacobian(
        &self,
        a1: &Array2<f32>,
        a2: &Array2<f32>,
        n_samples: usize,
    ) -> Array3<f32> {
        let n_inputs = self.w1.shape()[0];
        let n_outputs = self.w3.shape()[1];
        let mut jacobian = Array3::zeros((n_samples, n_outputs, n_inputs));

        // Compute Jacobian for each sample
        for i in 0..n_samples {
            // ReLU derivatives
            let d_relu1: Array1<f32> = a1.row(i).mapv(|v| if v > 0.0 { 1.0 } else { 0.0 });
            let d_relu2: Array1<f32> = a2.row(i).mapv(|v| if v > 0.0 { 1.0 } else { 0.0 });

            for output_idx in 0..n_outputs {
                let grad_h2 = self.w3.column(output_idx).to_owned(); // (hidden2,)
                let grad_a2 = &grad_h2 * &d_relu2; // (hidden2,)
                let grad_h1 = self.w2.dot(&grad_a2); // (hidden1,)
                let grad_a1 = &grad_h1 * &d_relu1; // (hidden1,)
                let grad_x = self.w1.dot(&grad_a1); // (n_inputs,)

                jacobian.slice_mut(s![i, output_idx, ..]).assign(&grad_x);
            }
        }

        jacobian
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn ones_network(input: usize, h1: usize, h2: usize, output: usize) -> TwoLayerNN {
        TwoLayerNN {
            w1: Array2::ones((input, h1)),
            b1: Array1::zeros(h1),
            w2: Array2::ones((h1, h2)),
            b2: Array1::zeros(h2),
            w3: Array2::ones((h2, output)),
            b3: Array1::zeros(output),
        }
    }

    #[test]
    fn forward_shape() {
        let nn = ones_network(7, 8, 8, 1);
        let x = Array2::ones((10, 7));
        let (output, jac) = nn.forward(x.view(), false);
        assert_eq!(output.shape(), &[10, 1]);
        assert!(jac.is_none());
    }

    #[test]
    fn jacobian_shape() {
        let nn = ones_network(7, 8, 8, 1);
        let x = Array2::ones((10, 7));
        let (_, jacobian) = nn.forward(x.view(), true);
        let j = jacobian.expect("jacobian requested");
        assert_eq!(j.shape(), &[10, 1, 7]);
    }

    #[test]
    fn architecture_matches_weights() {
        let nn = ones_network(5, 12, 9, 3);
        assert_eq!(nn.get_architecture(), (5, 12, 9, 3));
    }

    #[test]
    fn all_ones_network_known_answer() {
        // Network with all-ones weights and zero biases, ReLU activations,
        // and a positive input of all ones: every activation is positive so
        // ReLU is a no-op. Then the forward pass collapses to repeated
        // matmuls of ones, giving output = input_dim * h1 * h2.
        let (input_dim, h1, h2) = (4, 5, 6);
        let nn = ones_network(input_dim, h1, h2, 1);
        let x = Array2::ones((1, input_dim));
        let (out, _) = nn.forward(x.view(), false);
        let expected = (input_dim as f32) * (h1 as f32) * (h2 as f32);
        assert_relative_eq!(out[[0, 0]], expected, epsilon = 1e-5);
    }

    #[test]
    fn negative_input_silenced_by_relu() {
        // Purely negative pre-activation (through non-positive biases) should
        // produce zero output via ReLU-ReLU-linear path.
        let mut nn = ones_network(3, 4, 4, 2);
        nn.b1.fill(-1000.0);
        nn.b2.fill(-1000.0);
        let x = Array2::ones((2, 3));
        let (out, _) = nn.forward(x.view(), false);
        for v in out.iter() {
            assert_relative_eq!(*v, 0.0, epsilon = 1e-6);
        }
    }

    #[test]
    fn jacobian_zero_when_relu_dead() {
        // If the network is fully saturated at zero, gradients must be zero.
        let mut nn = ones_network(3, 4, 4, 2);
        nn.b1.fill(-1000.0);
        nn.b2.fill(-1000.0);
        let x = Array2::ones((2, 3));
        let (_, jacobian) = nn.forward(x.view(), true);
        let j = jacobian.expect("jacobian requested");
        for v in j.iter() {
            assert_relative_eq!(*v, 0.0, epsilon = 1e-12);
        }
    }
}
