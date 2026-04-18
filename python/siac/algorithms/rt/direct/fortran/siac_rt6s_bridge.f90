module siac_rt6s_bridge_mod
  use omp_lib
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer, parameter :: SIXS_OUTPUT_COUNT = 115

contains

  subroutine sixs_run_case_direct( &
      month, day, atmospheric_mode, radiosonde_altitude_km, radiosonde_pressure_mb, &
      radiosonde_temperature_k, radiosonde_water_g_m3, radiosonde_ozone_g_m3, &
      aerosol_mode, aerosol_mixture, aerosol_distribution_rmin, aerosol_distribution_rmax, &
      aerosol_distribution_component_count, aerosol_distribution_x1, aerosol_distribution_x2, &
      aerosol_distribution_x3, aerosol_distribution_cij, aerosol_distribution_rn, &
      aerosol_distribution_ri, aerosol_sun_count, aerosol_sun_radius, aerosol_sun_dvlogr, &
      aerosol_layer_count, aerosol_layer_height, aerosol_layer_aot, aerosol_layer_type, &
      reference_reflectance, spectral_wlinf, spectral_wlsup, spectral_response, &
      aerosol_model_path, surface_inhomo, surface_idirec, surface_target_mode, &
      surface_target_constant, surface_target_spectrum, surface_env_mode, &
      surface_env_constant, surface_env_spectrum, surface_radius_km, surface_brdf_model, &
      surface_brdf_params, surface_brdf_options, surface_brdf_struct, surface_brdf_optics, &
      surface_brdf_table_solar, surface_brdf_table_view, surface_brdf_spherical_albedo, &
      surface_brdf_directional_reflectance, atmospheric_correction_mode, &
      atmospheric_correction_value, sza_deg, saa_deg, vza_deg, vaa_deg, aot550, tcwv_cm, &
      tco3_atmcm, elevation_km, output_values, status_code)
    implicit none

    integer, intent(in) :: month
    integer, intent(in) :: day
    integer, intent(in) :: atmospheric_mode
    integer, intent(in) :: aerosol_mode
    integer, intent(in) :: aerosol_distribution_component_count
    integer, intent(in) :: aerosol_sun_count
    integer, intent(in) :: aerosol_layer_count
    integer, intent(in) :: aerosol_layer_type(50)
    integer, intent(in) :: surface_inhomo
    integer, intent(in) :: surface_idirec
    integer, intent(in) :: surface_target_mode
    integer, intent(in) :: surface_env_mode
    integer, intent(in) :: surface_brdf_model
    integer, intent(in) :: surface_brdf_options(5)
    integer, intent(in) :: atmospheric_correction_mode
    double precision, intent(in) :: radiosonde_altitude_km(34)
    double precision, intent(in) :: radiosonde_pressure_mb(34)
    double precision, intent(in) :: radiosonde_temperature_k(34)
    double precision, intent(in) :: radiosonde_water_g_m3(34)
    double precision, intent(in) :: radiosonde_ozone_g_m3(34)
    double precision, intent(in) :: aerosol_mixture(4)
    double precision, intent(in) :: aerosol_distribution_rmin
    double precision, intent(in) :: aerosol_distribution_rmax
    double precision, intent(in) :: aerosol_distribution_x1(4)
    double precision, intent(in) :: aerosol_distribution_x2(4)
    double precision, intent(in) :: aerosol_distribution_x3(4)
    double precision, intent(in) :: aerosol_distribution_cij(4)
    double precision, intent(in) :: aerosol_distribution_rn(20, 4)
    double precision, intent(in) :: aerosol_distribution_ri(20, 4)
    double precision, intent(in) :: aerosol_sun_radius(50)
    double precision, intent(in) :: aerosol_sun_dvlogr(50)
    double precision, intent(in) :: aerosol_layer_height(50)
    double precision, intent(in) :: aerosol_layer_aot(50)
    double precision, intent(in) :: reference_reflectance
    double precision, intent(in) :: spectral_wlinf
    double precision, intent(in) :: spectral_wlsup
    double precision, intent(in) :: spectral_response(1501)
    character(len=*), intent(in) :: aerosol_model_path
    double precision, intent(in) :: surface_target_constant
    double precision, intent(in) :: surface_target_spectrum(1501)
    double precision, intent(in) :: surface_env_constant
    double precision, intent(in) :: surface_env_spectrum(1501)
    double precision, intent(in) :: surface_radius_km
    double precision, intent(in) :: surface_brdf_params(12)
    double precision, intent(in) :: surface_brdf_struct(4)
    double precision, intent(in) :: surface_brdf_optics(3)
    double precision, intent(in) :: surface_brdf_table_solar(10, 13)
    double precision, intent(in) :: surface_brdf_table_view(10, 13)
    double precision, intent(in) :: surface_brdf_spherical_albedo
    double precision, intent(in) :: surface_brdf_directional_reflectance
    double precision, intent(in) :: atmospheric_correction_value
    double precision, intent(in) :: sza_deg
    double precision, intent(in) :: saa_deg
    double precision, intent(in) :: vza_deg
    double precision, intent(in) :: vaa_deg
    double precision, intent(in) :: aot550
    double precision, intent(in) :: tcwv_cm
    double precision, intent(in) :: tco3_atmcm
    double precision, intent(in) :: elevation_km
    double precision, intent(out) :: output_values(SIXS_OUTPUT_COUNT)
    integer, intent(out) :: status_code

    external :: sixs_case_core

    integer :: output_unit
    integer :: io_status
    integer :: i
    logical :: ier_case
    real :: radiosonde_altitude_local(34)
    real :: radiosonde_pressure_local(34)
    real :: radiosonde_temperature_local(34)
    real :: radiosonde_water_local(34)
    real :: radiosonde_ozone_local(34)
    real :: aerosol_mix_local(4)
    real :: aerosol_dist_x1_local(4)
    real :: aerosol_dist_x2_local(4)
    real :: aerosol_dist_x3_local(4)
    real :: aerosol_dist_cij_local(4)
    real :: aerosol_dist_rn_local(20, 4)
    real :: aerosol_dist_ri_local(20, 4)
    real :: aerosol_sun_radius_local(50)
    real :: aerosol_sun_dvlogr_local(50)
    real :: aerosol_layer_height_local(50)
    real :: aerosol_layer_aot_local(50)
    integer :: aerosol_layer_type_local(50)
    real :: spectral_response_local(1501)
    real :: surface_target_spectrum_local(1501)
    real :: surface_env_spectrum_local(1501)
    real :: surface_brdf_params_local(12)
    real :: surface_brdf_struct_local(4)
    real :: surface_brdf_optics_local(3)
    real :: surface_brdf_table_solar_local(10, 13)
    real :: surface_brdf_table_view_local(10, 13)
    real :: output_values_case(SIXS_OUTPUT_COUNT)
    character(len=1024) :: aerosol_file_local
    double precision :: nan64
    real :: nan32

    nan64 = ieee_value(0.0d0, ieee_quiet_nan)
    nan32 = ieee_value(0.0, ieee_quiet_nan)
    status_code = 0
    output_values = nan64

    aerosol_file_local = ' '
    if (len_trim(aerosol_model_path) > 0) then
      aerosol_file_local = aerosol_model_path(1:min(len_trim(aerosol_model_path), len(aerosol_file_local)))
    end if

    do i = 1, 34
      radiosonde_altitude_local(i) = real(radiosonde_altitude_km(i))
      radiosonde_pressure_local(i) = real(radiosonde_pressure_mb(i))
      radiosonde_temperature_local(i) = real(radiosonde_temperature_k(i))
      radiosonde_water_local(i) = real(radiosonde_water_g_m3(i))
      radiosonde_ozone_local(i) = real(radiosonde_ozone_g_m3(i))
    end do
    do i = 1, 4
      aerosol_mix_local(i) = real(max(aerosol_mixture(i), 0.0d0))
      aerosol_dist_x1_local(i) = real(max(aerosol_distribution_x1(i), 0.0d0))
      aerosol_dist_x2_local(i) = real(max(aerosol_distribution_x2(i), 0.0d0))
      aerosol_dist_x3_local(i) = real(max(aerosol_distribution_x3(i), 0.0d0))
      aerosol_dist_cij_local(i) = real(max(aerosol_distribution_cij(i), 0.0d0))
    end do
    do i = 1, 50
      aerosol_sun_radius_local(i) = real(max(aerosol_sun_radius(i), 0.0d0))
      aerosol_sun_dvlogr_local(i) = real(max(aerosol_sun_dvlogr(i), 0.0d0))
      aerosol_layer_height_local(i) = real(max(aerosol_layer_height(i), 0.0d0))
      aerosol_layer_aot_local(i) = real(max(aerosol_layer_aot(i), 0.0d0))
      aerosol_layer_type_local(i) = aerosol_layer_type(i)
    end do
    do i = 1, 1501
      spectral_response_local(i) = real(max(spectral_response(i), 0.0d0))
      surface_target_spectrum_local(i) = real(max(surface_target_spectrum(i), 0.0d0))
      surface_env_spectrum_local(i) = real(max(surface_env_spectrum(i), 0.0d0))
    end do
    do i = 1, 12
      surface_brdf_params_local(i) = real(surface_brdf_params(i))
    end do
    do i = 1, 4
      surface_brdf_struct_local(i) = real(surface_brdf_struct(i))
    end do
    do i = 1, 3
      surface_brdf_optics_local(i) = real(surface_brdf_optics(i))
    end do
    do i = 1, 20
      aerosol_dist_rn_local(i, :) = real(aerosol_distribution_rn(i, :))
      aerosol_dist_ri_local(i, :) = real(aerosol_distribution_ri(i, :))
    end do
    do i = 1, 13
      surface_brdf_table_solar_local(:, i) = real(surface_brdf_table_solar(:, i))
      surface_brdf_table_view_local(:, i) = real(surface_brdf_table_view(:, i))
    end do

    open(newunit=output_unit, status='scratch', action='readwrite', form='formatted', iostat=io_status)
    if (io_status /= 0) then
      status_code = -10
      return
    end if

    output_values_case = nan32
    call sixs_case_core( &
        output_unit, &
        real(max(sza_deg, 0.0d0)), &
        real(modulo(saa_deg, 360.0d0)), &
        real(max(vza_deg, 0.0d0)), &
        real(modulo(vaa_deg, 360.0d0)), &
        month, &
        day, &
        atmospheric_mode, &
        real(max(tcwv_cm, 0.0d0)), &
        real(max(tco3_atmcm, 0.0d0)), &
        radiosonde_altitude_local, &
        radiosonde_pressure_local, &
        radiosonde_temperature_local, &
        radiosonde_water_local, &
        radiosonde_ozone_local, &
        aerosol_mode, &
        aerosol_mix_local, &
        real(max(aerosol_distribution_rmin, 0.0d0)), &
        real(max(aerosol_distribution_rmax, 0.0d0)), &
        aerosol_distribution_component_count, &
        aerosol_dist_x1_local, &
        aerosol_dist_x2_local, &
        aerosol_dist_x3_local, &
        aerosol_dist_cij_local, &
        aerosol_dist_rn_local, &
        aerosol_dist_ri_local, &
        aerosol_sun_count, &
        aerosol_sun_radius_local, &
        aerosol_sun_dvlogr_local, &
        aerosol_layer_count, &
        aerosol_layer_height_local, &
        aerosol_layer_aot_local, &
        aerosol_layer_type_local, &
        aerosol_file_local, &
        real(max(aot550, 0.0d0)), &
        real(max(elevation_km, 0.0d0)), &
        real(spectral_wlinf), &
        real(spectral_wlsup), &
        spectral_response_local, &
        surface_inhomo, &
        surface_idirec, &
        surface_target_mode, &
        real(surface_target_constant), &
        surface_target_spectrum_local, &
        surface_env_mode, &
        real(surface_env_constant), &
        surface_env_spectrum_local, &
        real(max(surface_radius_km, 0.0d0)), &
        surface_brdf_model, &
        surface_brdf_params_local, &
        surface_brdf_options, &
        surface_brdf_struct_local, &
        surface_brdf_optics_local, &
        surface_brdf_table_solar_local, &
        surface_brdf_table_view_local, &
        real(surface_brdf_spherical_albedo), &
        real(surface_brdf_directional_reflectance), &
        real(reference_reflectance), &
        atmospheric_correction_mode, &
        real(atmospheric_correction_value), &
        ier_case, &
        output_values_case &
    )

    close(output_unit)

    if (ier_case) then
      status_code = 1
      return
    end if

    output_values = dble(output_values_case)
  end subroutine sixs_run_case_direct

  subroutine sixs_run_batch_impl( &
      month, day, atmospheric_mode, radiosonde_altitude_km, radiosonde_pressure_mb, &
      radiosonde_temperature_k, radiosonde_water_g_m3, radiosonde_ozone_g_m3, &
      aerosol_mode, aerosol_mixture, aerosol_distribution_rmin, aerosol_distribution_rmax, &
      aerosol_distribution_component_count, aerosol_distribution_x1, aerosol_distribution_x2, &
      aerosol_distribution_x3, aerosol_distribution_cij, aerosol_distribution_rn, &
      aerosol_distribution_ri, aerosol_sun_count, aerosol_sun_radius, aerosol_sun_dvlogr, &
      aerosol_layer_count, aerosol_layer_height, aerosol_layer_aot, aerosol_layer_type, &
      reference_reflectance, spectral_wlinf, spectral_wlsup, spectral_response, &
      aerosol_model_path, surface_inhomo, surface_idirec, surface_target_mode, &
      surface_target_constant, surface_target_spectrum, surface_env_mode, &
      surface_env_constant, surface_env_spectrum, surface_radius_km, surface_brdf_model, &
      surface_brdf_params, surface_brdf_options, surface_brdf_struct, surface_brdf_optics, &
      surface_brdf_table_solar, surface_brdf_table_view, surface_brdf_spherical_albedo, &
      surface_brdf_directional_reflectance, atmospheric_correction_mode, &
      atmospheric_correction_value, n_cases, sza_deg, saa_deg, vza_deg, vaa_deg, aot550, tcwv_cm, &
      tco3_atmcm, elevation_km, n_threads, output_values, status_code)
    implicit none

    integer, intent(in) :: month
    integer, intent(in) :: day
    integer, intent(in) :: atmospheric_mode
    integer, intent(in) :: aerosol_mode
    integer, intent(in) :: aerosol_distribution_component_count
    integer, intent(in) :: aerosol_sun_count
    integer, intent(in) :: aerosol_layer_count
    integer, intent(in) :: aerosol_layer_type(50)
    integer, intent(in) :: surface_inhomo
    integer, intent(in) :: surface_idirec
    integer, intent(in) :: surface_target_mode
    integer, intent(in) :: surface_env_mode
    integer, intent(in) :: surface_brdf_model
    integer, intent(in) :: surface_brdf_options(5)
    integer, intent(in) :: atmospheric_correction_mode
    integer, intent(in) :: n_cases
    integer, intent(in) :: n_threads
    double precision, intent(in) :: radiosonde_altitude_km(34)
    double precision, intent(in) :: radiosonde_pressure_mb(34)
    double precision, intent(in) :: radiosonde_temperature_k(34)
    double precision, intent(in) :: radiosonde_water_g_m3(34)
    double precision, intent(in) :: radiosonde_ozone_g_m3(34)
    double precision, intent(in) :: aerosol_mixture(4)
    double precision, intent(in) :: aerosol_distribution_rmin
    double precision, intent(in) :: aerosol_distribution_rmax
    double precision, intent(in) :: aerosol_distribution_x1(4)
    double precision, intent(in) :: aerosol_distribution_x2(4)
    double precision, intent(in) :: aerosol_distribution_x3(4)
    double precision, intent(in) :: aerosol_distribution_cij(4)
    double precision, intent(in) :: aerosol_distribution_rn(20, 4)
    double precision, intent(in) :: aerosol_distribution_ri(20, 4)
    double precision, intent(in) :: aerosol_sun_radius(50)
    double precision, intent(in) :: aerosol_sun_dvlogr(50)
    double precision, intent(in) :: aerosol_layer_height(50)
    double precision, intent(in) :: aerosol_layer_aot(50)
    double precision, intent(in) :: reference_reflectance
    double precision, intent(in) :: spectral_wlinf
    double precision, intent(in) :: spectral_wlsup
    double precision, intent(in) :: spectral_response(1501)
    character(len=*), intent(in) :: aerosol_model_path
    double precision, intent(in) :: surface_target_constant
    double precision, intent(in) :: surface_target_spectrum(1501)
    double precision, intent(in) :: surface_env_constant
    double precision, intent(in) :: surface_env_spectrum(1501)
    double precision, intent(in) :: surface_radius_km
    double precision, intent(in) :: surface_brdf_params(12)
    double precision, intent(in) :: surface_brdf_struct(4)
    double precision, intent(in) :: surface_brdf_optics(3)
    double precision, intent(in) :: surface_brdf_table_solar(10, 13)
    double precision, intent(in) :: surface_brdf_table_view(10, 13)
    double precision, intent(in) :: surface_brdf_spherical_albedo
    double precision, intent(in) :: surface_brdf_directional_reflectance
    double precision, intent(in) :: atmospheric_correction_value
    double precision, intent(in) :: sza_deg(n_cases)
    double precision, intent(in) :: saa_deg(n_cases)
    double precision, intent(in) :: vza_deg(n_cases)
    double precision, intent(in) :: vaa_deg(n_cases)
    double precision, intent(in) :: aot550(n_cases)
    double precision, intent(in) :: tcwv_cm(n_cases)
    double precision, intent(in) :: tco3_atmcm(n_cases)
    double precision, intent(in) :: elevation_km(n_cases)
    double precision, intent(out) :: output_values(SIXS_OUTPUT_COUNT, n_cases)
    integer, intent(out) :: status_code(n_cases)

    integer :: i

    if (n_threads > 0) then
      call omp_set_num_threads(n_threads)
    end if

    !$omp parallel do schedule(static)
    do i = 1, n_cases
      call sixs_run_case_direct( &
          month, day, atmospheric_mode, radiosonde_altitude_km, radiosonde_pressure_mb, &
          radiosonde_temperature_k, radiosonde_water_g_m3, radiosonde_ozone_g_m3, &
          aerosol_mode, aerosol_mixture, aerosol_distribution_rmin, aerosol_distribution_rmax, &
          aerosol_distribution_component_count, aerosol_distribution_x1, aerosol_distribution_x2, &
          aerosol_distribution_x3, aerosol_distribution_cij, aerosol_distribution_rn, &
          aerosol_distribution_ri, aerosol_sun_count, aerosol_sun_radius, aerosol_sun_dvlogr, &
          aerosol_layer_count, aerosol_layer_height, aerosol_layer_aot, aerosol_layer_type, &
          reference_reflectance, spectral_wlinf, spectral_wlsup, spectral_response, &
          aerosol_model_path, surface_inhomo, surface_idirec, surface_target_mode, &
          surface_target_constant, surface_target_spectrum, surface_env_mode, surface_env_constant, &
          surface_env_spectrum, surface_radius_km, surface_brdf_model, surface_brdf_params, &
          surface_brdf_options, surface_brdf_struct, surface_brdf_optics, surface_brdf_table_solar, &
          surface_brdf_table_view, surface_brdf_spherical_albedo, &
          surface_brdf_directional_reflectance, atmospheric_correction_mode, &
          atmospheric_correction_value, sza_deg(i), saa_deg(i), vza_deg(i), vaa_deg(i), aot550(i), tcwv_cm(i), &
          tco3_atmcm(i), elevation_km(i), output_values(:, i), status_code(i))
    end do
    !$omp end parallel do
  end subroutine sixs_run_batch_impl

end module siac_rt6s_bridge_mod


subroutine sixs_f2py_run_batch( &
    month, day, atmospheric_mode, radiosonde_altitude_km, radiosonde_pressure_mb, &
    radiosonde_temperature_k, radiosonde_water_g_m3, radiosonde_ozone_g_m3, &
    aerosol_mode, aerosol_mixture, aerosol_distribution_rmin, aerosol_distribution_rmax, &
    aerosol_distribution_component_count, aerosol_distribution_x1, aerosol_distribution_x2, &
    aerosol_distribution_x3, aerosol_distribution_cij, aerosol_distribution_rn, &
    aerosol_distribution_ri, aerosol_sun_count, aerosol_sun_radius, aerosol_sun_dvlogr, &
    aerosol_layer_count, aerosol_layer_height, aerosol_layer_aot, aerosol_layer_type, &
    reference_reflectance, spectral_wlinf, spectral_wlsup, spectral_response, aerosol_model_path, &
    surface_inhomo, surface_idirec, surface_target_mode, surface_target_constant, &
    surface_target_spectrum, surface_env_mode, surface_env_constant, surface_env_spectrum, &
    surface_radius_km, surface_brdf_model, surface_brdf_params, surface_brdf_options, &
    surface_brdf_struct, surface_brdf_optics, surface_brdf_table_solar, surface_brdf_table_view, &
    surface_brdf_spherical_albedo, surface_brdf_directional_reflectance, &
    atmospheric_correction_mode, atmospheric_correction_value, n_cases, sza_deg, saa_deg, &
    vza_deg, vaa_deg, aot550, tcwv_cm, tco3_atmcm, elevation_km, n_threads, output_values, status_code)
  use siac_rt6s_bridge_mod, only: SIXS_OUTPUT_COUNT, sixs_run_batch_impl
  implicit none

  integer, intent(in) :: month
  integer, intent(in) :: day
  integer, intent(in) :: atmospheric_mode
  integer, intent(in) :: aerosol_mode
  integer, intent(in) :: aerosol_distribution_component_count
  integer, intent(in) :: aerosol_sun_count
  integer, intent(in) :: aerosol_layer_count
  integer, intent(in) :: aerosol_layer_type(50)
  integer, intent(in) :: surface_inhomo
  integer, intent(in) :: surface_idirec
  integer, intent(in) :: surface_target_mode
  integer, intent(in) :: surface_env_mode
  integer, intent(in) :: surface_brdf_model
  integer, intent(in) :: surface_brdf_options(5)
  integer, intent(in) :: atmospheric_correction_mode
  integer, intent(in) :: n_cases
  integer, intent(in) :: n_threads
  double precision, intent(in) :: radiosonde_altitude_km(34)
  double precision, intent(in) :: radiosonde_pressure_mb(34)
  double precision, intent(in) :: radiosonde_temperature_k(34)
  double precision, intent(in) :: radiosonde_water_g_m3(34)
  double precision, intent(in) :: radiosonde_ozone_g_m3(34)
  double precision, intent(in) :: aerosol_mixture(4)
  double precision, intent(in) :: aerosol_distribution_rmin
  double precision, intent(in) :: aerosol_distribution_rmax
  double precision, intent(in) :: aerosol_distribution_x1(4)
  double precision, intent(in) :: aerosol_distribution_x2(4)
  double precision, intent(in) :: aerosol_distribution_x3(4)
  double precision, intent(in) :: aerosol_distribution_cij(4)
  double precision, intent(in) :: aerosol_distribution_rn(20, 4)
  double precision, intent(in) :: aerosol_distribution_ri(20, 4)
  double precision, intent(in) :: aerosol_sun_radius(50)
  double precision, intent(in) :: aerosol_sun_dvlogr(50)
  double precision, intent(in) :: aerosol_layer_height(50)
  double precision, intent(in) :: aerosol_layer_aot(50)
  double precision, intent(in) :: reference_reflectance
  double precision, intent(in) :: spectral_wlinf
  double precision, intent(in) :: spectral_wlsup
  double precision, intent(in) :: spectral_response(1501)
  character(len=*), intent(in) :: aerosol_model_path
  double precision, intent(in) :: surface_target_constant
  double precision, intent(in) :: surface_target_spectrum(1501)
  double precision, intent(in) :: surface_env_constant
  double precision, intent(in) :: surface_env_spectrum(1501)
  double precision, intent(in) :: surface_radius_km
  double precision, intent(in) :: surface_brdf_params(12)
  double precision, intent(in) :: surface_brdf_struct(4)
  double precision, intent(in) :: surface_brdf_optics(3)
  double precision, intent(in) :: surface_brdf_table_solar(10, 13)
  double precision, intent(in) :: surface_brdf_table_view(10, 13)
  double precision, intent(in) :: surface_brdf_spherical_albedo
  double precision, intent(in) :: surface_brdf_directional_reflectance
  double precision, intent(in) :: atmospheric_correction_value
  double precision, intent(in) :: sza_deg(n_cases)
  double precision, intent(in) :: saa_deg(n_cases)
  double precision, intent(in) :: vza_deg(n_cases)
  double precision, intent(in) :: vaa_deg(n_cases)
  double precision, intent(in) :: aot550(n_cases)
  double precision, intent(in) :: tcwv_cm(n_cases)
  double precision, intent(in) :: tco3_atmcm(n_cases)
  double precision, intent(in) :: elevation_km(n_cases)
  double precision, intent(inout) :: output_values(SIXS_OUTPUT_COUNT, n_cases)
  integer, intent(inout) :: status_code(n_cases)

  call sixs_run_batch_impl( &
      month, day, atmospheric_mode, radiosonde_altitude_km, radiosonde_pressure_mb, &
      radiosonde_temperature_k, radiosonde_water_g_m3, radiosonde_ozone_g_m3, aerosol_mode, &
      aerosol_mixture, aerosol_distribution_rmin, aerosol_distribution_rmax, &
      aerosol_distribution_component_count, aerosol_distribution_x1, aerosol_distribution_x2, &
      aerosol_distribution_x3, aerosol_distribution_cij, aerosol_distribution_rn, &
      aerosol_distribution_ri, aerosol_sun_count, aerosol_sun_radius, aerosol_sun_dvlogr, &
      aerosol_layer_count, aerosol_layer_height, aerosol_layer_aot, aerosol_layer_type, &
      reference_reflectance, spectral_wlinf, spectral_wlsup, spectral_response, aerosol_model_path, &
      surface_inhomo, surface_idirec, surface_target_mode, surface_target_constant, &
      surface_target_spectrum, surface_env_mode, surface_env_constant, surface_env_spectrum, &
      surface_radius_km, surface_brdf_model, surface_brdf_params, surface_brdf_options, &
      surface_brdf_struct, surface_brdf_optics, surface_brdf_table_solar, surface_brdf_table_view, &
      surface_brdf_spherical_albedo, surface_brdf_directional_reflectance, &
      atmospheric_correction_mode, atmospheric_correction_value, n_cases, sza_deg, saa_deg, &
      vza_deg, vaa_deg, aot550, tcwv_cm, tco3_atmcm, elevation_km, n_threads, output_values, status_code)
end subroutine sixs_f2py_run_batch
