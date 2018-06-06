module radar_simulator

contains
   subroutine simulate_radar( &
      errorstatus, & !in
      wavelength, & !in
      particle_spectrum, & !in
      PIA, & !in
      spectral_broadening, & !in
      n_heights, & !in
      radar_Pnoise, & !in
      radar_max_V, & !in
      radar_min_V, & !in
      radar_nfft, & !in
      radar_nfft_aliased, & !in
      radar_no_Ave, & !in
      radar_aliasing_nyquist_interv, & !in
      radar_K2, & !in
      seed, & !in
      noise_turb_spectra & !out
      )
      ! This routine takes the backscattering spectrum depending on Doppler velocity,
      ! adds noise and turbulence and simulates temporal averaging
      !
      ! based on Spectra_simulator by P. Kollias
      ! converted from Matlab to Fortran by M. Maahn (2012)
      !

      use kinds
      use constants
      use report_module
      implicit none

      real(kind=dbl), intent(in) ::  wavelength !heigth of layer in m
      real(kind=dbl), dimension(n_heights, radar_nfft_aliased), intent(in):: particle_spectrum !backscattering particle spectrum per Doppler velocity [mm⁶/m³/(m/s)] NON-SI
      real(kind=dbl), dimension(n_heights), intent(in) ::  PIA !path inetgrated attenuation in dB
      real(kind=dbl), dimension(n_heights), intent(in) ::  spectral_broadening
      real(kind=dbl), dimension(n_heights), intent(in) ::  radar_Pnoise !noise in linear units
      integer, intent(in) :: n_heights
      real(kind=dbl), intent(in) ::  radar_max_V
      real(kind=dbl), intent(in) ::  radar_min_V
      integer, intent(in) ::  radar_nfft
      integer, intent(in) ::  radar_nfft_aliased
      integer, intent(in) ::  radar_aliasing_nyquist_interv
      integer, intent(in) ::  radar_no_Ave
      real(kind=dbl), intent(in) :: radar_K2
      integer, intent(in) ::  seed
      real(kind=dbl), dimension(n_heights, radar_nfft), intent(out):: noise_turb_spectra

      integer :: hh

      integer(kind=long), intent(out) :: errorstatus
      integer(kind=long) :: err
      character(len=80) :: msg
      character(len=15) :: nameOfRoutine = 'radar_simulator'

      if (verbose >= 2) call report(info, 'Start of ', nameOfRoutine)
      err = 0

      noise_turb_spectra(:, :) = -9999.d0

      do hh = 1, n_heights

         if (ANY(ISNAN(particle_spectrum(hh, :))) .or. ISNAN(spectral_broadening(hh)) .or. ISNAN(PIA(hh))) then
            if (verbose >= 2) print *, 'skipping due to NAN', hh
            continue
         end if

         call simulate_radar_one( &
            errorstatus, &
            wavelength, &
            particle_spectrum(hh, :), &
            PIA(hh), &
            spectral_broadening(hh), &
            radar_Pnoise(hh), &
            radar_max_V, &
            radar_min_V, &
            radar_nfft, &
            radar_nfft_aliased, &
            radar_no_Ave, &
            radar_aliasing_nyquist_interv, &
            radar_K2, &
            seed, &
            noise_turb_spectra(hh, :) &
            )

         if (err /= 0) then
            msg = 'error in estimate_spectralBroadening_one!'
            call report(err, msg, nameOfRoutine)
            errorstatus = err
            return
         end if

      end do
      if (verbose >= 2) call report(info, 'End of ', nameOfRoutine)

   end subroutine simulate_radar

   subroutine simulate_radar_one( &
      errorstatus, & !in
      wavelength, & !in
      particle_spectrum, & !in
      PIA, & !in
      spectral_broadening, & !in
      radar_Pnoise, & !in
      radar_max_V, & !in
      radar_min_V, & !in
      radar_nfft, & !in
      radar_nfft_aliased, & !in
      radar_no_Ave, & !in
      radar_aliasing_nyquist_interv, & !in
      radar_K2, & !in
      seed, & !in
      noise_turb_spectra & !out
      )
      ! This routine takes the backscattering spectrum depending on Doppler velocity,
      ! adds noise and turbulence and simulates temporal averaging
      !
      ! based on Spectra_simulator by P. Kollias
      ! converted from Matlab to Fortran by M. Maahn (2012)
      !
      ! out is saved directly to vars_output module

      use kinds
      use constants
      use report_module
      use random_module, only:get_random

      implicit none

      real(kind=dbl), intent(in) ::  wavelength !heigth of layer in m
      real(kind=dbl), dimension(radar_nfft_aliased), intent(in):: particle_spectrum !backscattering particle spectrum per Doppler velocity [mm⁶/m³/(m/s)] NON-SI
      real(kind=dbl), intent(in) ::  PIA !path inetgrated attenuation in dB
      real(kind=dbl), intent(in) ::  spectral_broadening
      real(kind=dbl), intent(in) ::  radar_Pnoise !noise in linear units
      real(kind=dbl), intent(in) ::  radar_max_V
      real(kind=dbl), intent(in) ::  radar_min_V
      integer, intent(in) ::  radar_nfft
      integer, intent(in) ::  radar_nfft_aliased
      integer, intent(in) ::  radar_aliasing_nyquist_interv
      integer, intent(in) ::  radar_no_Ave
      real(kind=dbl), intent(in) :: radar_K2
      integer, intent(in) ::  seed
      real(kind=dbl), dimension(radar_nfft), intent(out):: noise_turb_spectra

      real(kind=dbl) ::  back !volumetric backscattering crossection in m²/m³
      real(kind=dbl), dimension(radar_nfft_aliased) :: particle_spectrum_att
      real(kind=dbl), dimension(radar_nfft_aliased) :: spectra_velo_aliased
      real(kind=dbl), dimension(radar_nfft_aliased):: turb
      real(kind=dbl), dimension(radar_nfft*radar_no_Ave):: x_noise
      real(kind=dbl), dimension(radar_no_Ave, radar_nfft):: noise_turb_spectra_tmp
      real(kind=dbl), dimension(radar_nfft):: snr_turb_spectra, &
                                              spectra_velo, turb_spectra_aliased
      integer::quailty_aliasing
      real(kind=dbl), dimension(2*radar_nfft_aliased - 1):: turb_spectra
      logical, parameter :: use_fft = .true.
      real(kind=dbl):: SNR, del_v, K2, Ze_back, K, &
                       min_V_aliased, max_V_aliased
      integer(kind=long) :: ii, tt, ts_imin, ts_imax, startI, stopI
      integer(kind=long), intent(out) :: errorstatus
      integer(kind=long) :: err = 0
      character(len=80) :: msg
      character(len=19) :: nameOfRoutine = 'radar_simulator_one'

      interface
         subroutine convolution(errorstatus, X, M, A, N, use_fft, Y)
            use kinds
            implicit none
            integer(kind=long), intent(out) :: errorstatus
            INTEGER, intent(in) :: M ! Size of input vector X
            INTEGER, intent(in) :: N ! Size of convolution filter A
            REAL(kind=dbl), intent(in), DIMENSION(M) :: X
            REAL(kind=dbl), intent(in), DIMENSION(N) :: A
            logical, intent(in) :: use_fft
            REAL(kind=dbl), intent(out), DIMENSION(M + N - 1) :: Y
         end subroutine convolution

         subroutine radar_hildebrand_sekhon(errorstatus, spectrum, n_ave, n_ffts, &
                                            noise_mean, noise_max)
            use kinds
            implicit none
            integer(kind=long), intent(out) :: errorstatus
            integer, intent(in) :: n_ave, n_ffts
            real(kind=dbl), dimension(n_ffts), intent(in) :: spectrum
            real(kind=dbl), intent(out) :: noise_mean
            real(kind=dbl), intent(out) :: noise_max
         end subroutine radar_hildebrand_sekhon
      end interface
      if (verbose >= 2) call report(info, 'Start of ', nameOfRoutine)
      err = 0
      back = SUM(particle_spectrum)

      call assert_false(err, (ANY(ISNAN(particle_spectrum))), &
                        "got nan in values in backscattering spectrum")
      call assert_false(err, ISNAN(back) .or. (back < 0.d0), &
                        "got nan or negative value in linear Ze")
      call assert_false(err, (SUM(particle_spectrum) < 0.d0), &
                        "sum particle_spectrum < 0")
      if (err > 0) then
         errorstatus = fatal
         msg = "assertation error"
         call report(errorstatus, msg, nameOfRoutine)
         return
      end if

      ! get |K|**2 and lambda
      K2 = radar_K2

      !transform backscattering in linear reflectivity units, 10*log10(back) would be in dBz
      Ze_back = 1.d18*(1.d0/(K2*pi**5))*back*(wavelength)**4 ![mm⁶/m³]

      !take care of path integrated attenuation
      Ze_back = Ze_back/10d0**(0.1d0*PIA)
      particle_spectrum_att = particle_spectrum(:)/10d0**(0.1d0*PIA)

      if (verbose >= 2) print *, "particle_spectrum_att"
      if (verbose >= 2) print *, particle_spectrum_att

      !get delta velocity
      del_v = (radar_max_V - radar_min_V)/radar_nfft ![m/s]
      !create array from min_v to max_v iwth del_v spacing -> velocity spectrum of radar
      spectra_velo = (/(((ii*del_v) + radar_min_V), ii=0, radar_nfft - 1)/) ! [m/s]

      !same for the extended spectrum
      min_V_aliased = radar_min_V - radar_aliasing_nyquist_interv*(radar_max_V - radar_min_V)
      max_V_aliased = radar_max_V + radar_aliasing_nyquist_interv*(radar_max_V - radar_min_V)
      spectra_velo_aliased = (/(((ii*del_v) + min_V_aliased), ii=0, radar_nfft_aliased - 1)/) ! [m/s]

      ! The convolution function will return a vector of length 2*radar_nfft_aliased-1.
      ! Get the indices to cut out the required part.
      ! Error here result in a shifted spectrum
      ts_imin = (radar_nfft_aliased/2) + 1
      ts_imax = 2*(radar_nfft_aliased/2)

      !get turbulence (no turbulence in clear sky...)
      turb(:) = 0.d0
      if ((spectral_broadening > 0.d0) .and. (back > 0)) then
         do tt = 1, radar_nfft_aliased
            ! gaussian function with same length as radar spectrum, centered around zero
            turb(tt) = exp(-(spectra_velo_aliased(tt) - 0)**2.d0/(2.d0*spectral_broadening**2.d0))
         end do
         turb(:) = turb/sum(turb) !normalize to unity area
      end if

      ! skip in case ss was so little that sum(turb) is zero)
      if ((SUM(turb) > 0.d0) .and. (back > 0)) then
         call assert_false(err, (ANY(ISNAN(turb)) .or. (SUM(turb) <= 0.d0)), &
                           "got nan or negative value in linear turb")
         call assert_true(err, ((ABS(turb(1)) < almostZero) .and. (ABS(turb(radar_nfft_aliased)) < almostZero)), &
                          "increase radar_aliasing_nyquist_interv, turbulence too large")
         if (err > 0) then
            errorstatus = fatal
            msg = "assertation error"
            call report(errorstatus, msg, nameOfRoutine)
            return
         end if !(err > 0)

         !convolute spectrum and noise
         call convolution(err, particle_spectrum_att, radar_nfft_aliased, turb, radar_nfft_aliased, use_fft, turb_spectra)
         if (err /= 0) then
            msg = 'error in convolution!'
            call report(err, msg, nameOfRoutine)
            errorstatus = err
            return
         end if !(err > 0)

         !I don't like Nans values here
         where (ISNAN(turb_spectra)) turb_spectra = 0.d0
         ! negative number are resulting from numerical effects
         where (turb_spectra < 0) turb_spectra = 0.d0
      else !no turb
         turb_spectra(:) = 0.d0
         turb_spectra(ts_imin:ts_imax) = particle_spectrum_att
      end if

      if (verbose >= 10) print *, "turb_spectra"
      if (verbose >= 10) print *, SHAPE(turb_spectra)
      if (verbose >= 10) print *, turb_spectra

      quailty_aliasing = 0
      !lets look for aliasing effects. if we calculated radar_aliasing_nyquist_interv for a broader spectrum than necessary, fold it again:
      if (radar_aliasing_nyquist_interv > 0) then
         turb_spectra_aliased(:) = 0.d0
         do, ii = 1, 1 + 2*radar_aliasing_nyquist_interv
         !get indices
         startI = ts_imin + (ii - 1)*radar_nfft
         stopI = ts_imax - (1 + 2*radar_aliasing_nyquist_interv - ii)*radar_nfft
         if ((ii .ne. radar_aliasing_nyquist_interv + 1) &
             .and. (SUM(turb_spectra(startI:stopI)) .gt. almostZero)) then
            quailty_aliasing = 1
         end if
         !appy aliasing
         turb_spectra_aliased = turb_spectra_aliased + turb_spectra(startI:stopI)
      end do
   else ! aliasing effects not considered
      turb_spectra_aliased = turb_spectra(ts_imin:ts_imax)
   end if

   if (verbose > 2) then
      if (quailty_aliasing .ne. 0) then
         print *, "radar quality: aliasing found"
      else
         print *, "radar quality: NO aliasing found"
      end if
   end if

   call assert_false(err, (ANY(ISNAN(turb_spectra_aliased)) .or. ANY(turb_spectra_aliased < 0.d0)), &
                     "got nan or negative value in linear turb_spectra_aliased")
   ! call assert_false(err,(ALL(turb_spectra_aliased==0)),&
   !     "all values of turb_spectra_aliased == 0")
   if (err > 0) then
      errorstatus = fatal
      msg = "assertation error"
      call report(errorstatus, msg, nameOfRoutine)
      return
   end if

   if (verbose == 666) then
      print *, "##########################################"
      print *, "velocity (m/s)"
      print *, spectra_velo
      print *, "##########################################"
      print *, "particle_spec with turbulence (v) [mm⁶/m³/(m/s)] ", MAXVAL(turb_spectra_aliased)
      print *, SHAPE(turb_spectra_aliased)
      print *, turb_spectra_aliased
      print *, "##########################################"
   end if

   !get the SNR
   SNR = 10.d0*log10(Ze_back/radar_Pnoise)
   !this here is for scaling, if we have now a wrong Ze due to all the turbulence, rescaling etc.
   !Can happen e.g. due to numeric issues when applying very large or very small turbulence. Skip
   !if Ze_back ==0.
   if (Ze_back > 0) then
      K = (Ze_back/SUM(turb_spectra_aliased*del_v))
   else
      K = 1.d0
   end if

   if (verbose >= 4) print *, "first K", K
   snr_turb_spectra = (K*turb_spectra_aliased + radar_Pnoise/(radar_nfft*del_v))
   !   snr_turb_spectra =turb_spectra_aliased + radar_Pnoise/(radar_nfft*del_v)

   call assert_false(err, (ANY(ISNAN(snr_turb_spectra)) .or. ANY(snr_turb_spectra < 0.d0)), &
                     "got nan or negative value in linear snr_turb_spectra")
   call assert_false(err, (ISNAN(K)) .or. (K >= HUGE(K)), &
                     "K is nan or infinitive")
   if (err > 0) then
      errorstatus = fatal
      msg = "assertation error"
      call report(errorstatus, msg, nameOfRoutine)
      return
   end if

   if (radar_no_Ave .eq. 0) then !0 means infinity-> no noise
      noise_turb_spectra = snr_turb_spectra
   else
      !get noise.
      if (verbose > 2) print *, "get noise"
      call get_random(err, radar_no_Ave*radar_nfft, seed, x_noise)
      if (err /= 0) then
         msg = 'error in random!'
         call report(err, msg, nameOfRoutine)
         errorstatus = err
         return
      end if
      do tt = 1, radar_no_Ave
         noise_turb_spectra_tmp(tt, :) = -log(x_noise((tt - 1)*radar_nfft + 1:tt*radar_nfft))*snr_turb_spectra
      end do

      if (radar_no_Ave .eq. 1) then
         noise_turb_spectra = noise_turb_spectra_tmp(1, :)
      else
         noise_turb_spectra = SUM(noise_turb_spectra_tmp, DIM=1)/radar_no_Ave
      end if
   end if

   !spetial output for testing the radar simulator
   if (verbose == 666) then
      print *, "##########################################"
      print *, "particle_spec with turbulence and noise (v) [mm⁶/m³/(m/s)] ", MAXVAL(noise_turb_spectra)
      print *, noise_turb_spectra
      print *, "##########################################"
   end if

   !apply spectral resolution
   noise_turb_spectra = noise_turb_spectra*del_v !now [mm⁶/m³]

   if (verbose >= 4) then
      print *, "first K", K
      print *, "TOTAL", " Ze back", 10*log10(Ze_back)
      print *, "TOTAL", " Ze SUM(particle_spectrum)*del_v", 10*log10(SUM(particle_spectrum)*del_v)
      print *, "TOTAL", " Ze SUM(particle_spectrum_att)*del_v", 10*log10(SUM(particle_spectrum_att)*del_v)
      print *, "TOTAL", " Ze SUM(turb_spectra)*del_v", 10*log10(SUM(turb_spectra)*del_v)
      print *, "TOTAL", " Ze SUM(turb_spectra_aliased)*del_v", 10*log10(SUM(turb_spectra_aliased)*del_v)
      print *, "TOTAL", " Ze SUM(snr_turb_spectra)*del_v", 10*log10(SUM(snr_turb_spectra)*del_v)
      print *, "TOTAL", " Ze SUM(noise_turb_spectra)*del_v", 10*log10(SUM(noise_turb_spectra))
      print *, "TOTAL", " Ze SUM(snr_turb_spectra)*del_v-radar_Pnoise", 10*log10(SUM(snr_turb_spectra)*del_v - radar_Pnoise)
      print *, "TOTAL", " Ze SUM(noise_turb_spectra)*del_v-radar_Pnoise", 10*log10(SUM(noise_turb_spectra) - radar_Pnoise)
   end if

   if (verbose >= 5) then
      print *, " linear spectrum befor receiver uncertainty"
      print *, noise_turb_spectra
      print *, "#####################"
   end if

   if (verbose >= 5) then
      print *, "final linear spectrum"
      print *, noise_turb_spectra
      print *, "#####################"
   end if

   errorstatus = err
   if (verbose >= 2) call report(info, 'End of ', nameOfRoutine)
   return
end subroutine simulate_radar_one

end module radar_simulator

