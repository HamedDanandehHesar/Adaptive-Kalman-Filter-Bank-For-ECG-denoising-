<img width="560" height="420" alt="untitled1" src="https://github.com/user-attachments/assets/4a885254-2490-463b-95e3-c060329bcc3f" />
<img width="560" height="420" alt="Untitled" src="https://github.com/user-attachments/assets/91543456-d0a7-46b7-ae5e-c920ec2fcf69" />
<img width="560" height="420" alt="untitled2" src="https://github.com/user-attachments/assets/eefebd2d-3f62-4ac7-85b6-a5adb29bad2e" />


# An Adaptive Kalman Filter Bank for ECG Denoising

This repository contains the MATLAB implementation of the method proposed in our paper:

**Hamed Danandeh Hesar and Maryam Mohebbi**  
**An Adaptive Kalman Filter Bank for ECG Denoising**  
*IEEE Journal of Biomedical and Health Informatics*, 2020  
https://doi.org/10.1109/JBHI.2020.2982935

The proposed approach uses an **adaptive Kalman filter bank** to denoise ECG signals without requiring a predefined ECG morphology model.  
The method separates the ECG into two components:

- **QRS complex**: high-frequency part
- **P and T waves**: low-frequency part

Then, two adaptive Kalman filters are applied independently and their parameters are iteratively updated using the **Expectation-Maximization (EM)** algorithm.

---

# Overview

ECG denoising is challenging because of:

- nonstationary noise
- muscle artifacts
- baseline disturbances
- morphology variations caused by arrhythmias

To address these issues, this method combines:

- **R-peak detection**
- **signal decomposition into QRS and P-T parts**
- **adaptive Kalman filtering**
- **EM-based parameter estimation**
- **Bryson and Henrikson’s technique** for improved robustness to colored and nonstationary noise

---

# Method Summary

The processing pipeline is as follows:

1. Load ECG signal from a `.mat` file
2. Add noise to the signal
3. Detect R-peaks using the Pan–Tompkins algorithm
4. Separate the ECG into:
   - QRS-related samples
   - P-T related samples
5. Construct two independent signals
6. Initialize two adaptive Kalman filter models
7. Iteratively estimate and update the model parameters using EM
8. Smooth both ECG components
9. Reconstruct the denoised ECG signal

---

# Input Data

The code expects a MATLAB `.mat` file containing:

- `x` : ECG signal matrix
- `fs` : sampling frequency in Hz

Example:

```matlab
x  % ECG signal
fs % sampling frequency
```

The first channel of `x` is used as the ECG signal.

---

# Noise Model

To evaluate performance, white Gaussian noise is added to the ECG signal:

```matlab
SNR = 6;
y = awgn(x, SNR, 'measured');
```

The noisy ECG is then used as the input to the denoising pipeline.

---

# QRS Detection

R-peaks are detected using the Pan–Tompkins algorithm:

```matlab
qrs_positions = pantompkins_qrs(x_noisy, fs);
```

The detected R-peaks are used to identify the QRS region and the remaining samples corresponding to the P and T waves.

---

# Signal Decomposition

The ECG is split into two parts:

- **QRS interval**
- **P-T interval**

For each detected R-peak, the QRS region is approximated by samples around the peak:

```matlab
ind_qrs = [ind_Rpeak(i)-round(0.03*fs) : ind_Rpeak(i)+round(0.07*fs)];
```

The remaining samples are assigned to the P-T part.

---

# Adaptive Kalman Filter Bank

The method uses **two separate adaptive Kalman filters**:

- one for the **P-T waveform**
- one for the **QRS waveform**

Each filter has its own state-space model and covariance matrices.

---

## State-Space Form

The general Kalman model can be written as:

```math
x_k = A x_{k-1} + w_k
```

```math
z_k = H x_k + v_k
```

where:

- `x_k` is the state vector
- `z_k` is the observation
- `A` is the state transition matrix
- `H` is the measurement matrix
- `w_k` is process noise
- `v_k` is measurement noise

---

# EM-Based Parameter Adaptation

The filter parameters are updated iteratively using the **EM algorithm**.

At each iteration, the following parameters are adapted:

- `A` : state transition matrix
- `Q` : process noise covariance
- `R` : measurement noise covariance
- `H` : measurement matrix

The code repeats this process for a fixed number of iterations:

```matlab
max_number_of_iterations = 15;
```

This adaptation allows the filter to better fit different ECG morphologies and noise conditions.

---

# Square-Root Implementation

The code uses a **square-root Kalman filtering/smoothing formulation** to improve numerical stability.

This approach is beneficial because it:

- reduces numerical instability
- improves matrix factorization robustness
- avoids covariance degeneration in iterative updates

---

# Reconstruction

After filtering, the two components are reassembled into the original ECG timeline:

- denoised QRS samples
- denoised P-T samples

Finally, the complete ECG signal is reconstructed and plotted against the noisy and original signals.

---

# Output

The script generates the following figures:

## Figure 1
Noisy ECG signal with detected R-peaks.

## Figure 2
Iteration-wise reconstruction during EM updates.

## Figure 3
Comparison of:

- original ECG
- denoised ECG
- noisy ECG

---

# Main Parameters

Important parameters used in the script:

- `SNR = 6`  
  Noise level for simulation

- `max_number_of_iterations = 15`  
  Number of EM iterations

- `Psi_pt = 0.001`  
  Modification parameter for P-T part

- `Psi_qrs = 0.01`  
  Modification parameter for QRS part

- `R_pt = 15`, `R_qrs = 1.5`  
  Initial measurement noise covariance values

- `Q_pt = eye(14)`, `Q_qrs = eye(14)`  
  Initial process noise covariance values

---

# Related Publication

If you use this repository in your research, please cite the following paper:

**Hamed Danandeh Hesar and Maryam Mohebbi**  
**An Adaptive Kalman Filter Bank for ECG Denoising**  
*IEEE Journal of Biomedical and Health Informatics*, 2020  
https://doi.org/10.1109/JBHI.2020.2982935

---

# BibTeX

```bibtex
@article{hesar2020adaptive,
  title={An Adaptive Kalman Filter Bank for ECG Denoising},
  author={Hesar, Hamed Danandeh and Mohebbi, Maryam},
  journal={IEEE Journal of Biomedical and Health Informatics},
  year={2020},
  publisher={IEEE},
  doi={10.1109/JBHI.2020.2982935}
}
```

---

# Requirements

- MATLAB
- Signal Processing Toolbox
- A MATLAB implementation of Pan–Tompkins QRS detection
- Access to `awgn` function for noise generation

---

# Repository Structure

```text
ECG-Kalman-Denoising/
├── main_script.m
├── pantompkins_qrs.m
├── data/
│   ├── 1-Mit-normal.mat
│   ├── 2-Mit-normal.mat
│   └── ...
└── README.md
```

---

# Notes

- This implementation is designed for **research and educational use**
- The method is especially effective for ECG signals corrupted by **stationary and nonstationary noise**
- It is robust to morphology variations caused by **abnormal heart rhythms**

---
