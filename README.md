
# multimedia-transmission
# Unified Digital Communication System for Multimedia Transmission

## Overview
This project implements an end-to-end digital communication system for transmitting multimedia data (image and audio) over noisy wireless channels. It combines hardware implementation and simulation to evaluate performance under practical and theoretical conditions.

## Key Contributions
- End-to-end pipeline: Input → Encoding → Modulation → Transmission → Reception → Reconstruction  
- Hardware implementation of BPSK using SDR/LabVIEW  
- Simulation of QPSK, 8-QAM, and 16-QAM in MATLAB  
- Pulse shaping using Rectangular and Root Raised Cosine (RRC) filters  
- BER analysis under AWGN channel  

## System Highlights
- Supports both image and audio transmission  
- Frame-based transmission with synchronization and metadata  
- Comparison of modulation schemes and filter performance  
- Visualization using eye diagrams and constellation plots  

## Performance Insight
- RRC filters significantly improve performance over rectangular filters  
- Higher roll-off factor reduces ISI and lowers BER  
- BPSK provides robust performance in hardware  
- Higher-order modulations improve data rate but are more noise-sensitive  

## Current Status
- Completed:
  - Full transmitter and receiver design  
  - Hardware testing with BPSK  
  - Simulation of QPSK, 8-QAM, and 16-QAM  
  - BER performance evaluation  

- Limitations:
  - Hardware restricted to BPSK implementation  
  - Synchronization and noise challenges in real-time transmission  

## Tech Stack
- MATLAB (signal processing, simulation, modulation)
- LabVIEW + SDR (hardware implementation)

## Conclusion
The project demonstrates a reliable multimedia communication system and highlights the trade-off between spectral efficiency and noise robustness. Proper filter design plays a critical role in improving system performance.
>>>>>>> 5405a37 (Initial commit with project structure)
