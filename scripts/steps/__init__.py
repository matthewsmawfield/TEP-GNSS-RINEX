"""
TEP-GNSS-RINEX Analysis Pipeline - Step Modules
================================================

This package contains the analysis steps for the Temporal Equivalence Principle
(TEP) validation using raw RINEX GNSS observations.

Step Modules:
    step_1_0_data_acquisition: Download and process RINEX/NAV data from CDDIS
    step_1_1_generate_dynamic_50_metadata: Generate optimized station selection
    step_2_0_raw_spp_analysis: Core spatial correlation analysis
    step_2_1_control_tests: Control and validation tests
    step_2_1_elevation_analysis: Satellite elevation mask analysis
    step_2_2_anisotropy_analysis: Directional anisotropy characterization
    step_2_3_kp_stratification: Geomagnetic storm stratification
    step_2_3_temporal_analysis: Multi-year temporal stability
    step_2_4_diurnal_analysis: Local solar time patterns
    step_2_4_null_tests: Null hypothesis tests
    step_2_5_orbital_coupling: Earth orbital velocity correlation
    step_2_6_planetary_events: Planetary alignment event analysis
    step_2_7_cmb_frame_analysis: Cosmic microwave background frame alignment

Author: Matthew Lukin Smawfield
Date: 9 December 2025
License: CC-BY-4.0
"""
