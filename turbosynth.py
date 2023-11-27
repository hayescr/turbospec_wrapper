import os
import abund_utils as au
from atmos_wrapper.atmos_manager import AtmosManager
from turbospec_wrapper.turbospec_manager import TurbospecManager
from turbospec_wrapper.turbospec_tools import ConvolutionManager


class TurboSynth:
    def __init__(self, path='.', turbopath='.', faltbon_loc='Utilities'):
        '''
        '''
        self.path = path
        self.turbopath = turbopath
        self.faltbon_path = f"{turbopath}/{faltbon_loc}"

        self.turbodaemon = TurbospecManager(
            inpath=path, turbopath=self.turbopath)
        self.convoldaemon = ConvolutionManager(
            inpath=path, faltbon_path=self.faltbon_path)

    def init_params(self, teff, logg, feh, vmicro=None, cfe=0.0, alphafe=0.0):
        '''
        Initializes the parameters
        '''
        self.teff = teff
        self.logg = logg
        self.feh = feh
        self.vmicro = vmicro
        self.cfe = cfe
        self.alphafe = alphafe

    def init_abunds(self, abunds=None, solar_reference='Asplund2009'):
        '''
        Initializes the abundances for Turbospectrum and sets a default
            C12/C13 ratio (15 for giants and 90 for dwarfs) and N14/N15 ratio
        '''
        self.turbodaemon.set_abund(metals=self.feh, alphas=self.alphafe,
                                   abundances=abunds,
                                   solar_reference=solar_reference)

        if self.logg < 3.8:
            c12c13 = 15.
        else:
            c12c13 = 90.

        n14n15 = 330.

        self.isotopes = {6.012: iso_frac2('major', c12c13),
                         6.013: iso_frac2('minor', c12c13),
                         7.014: iso_frac2('major', n14n15),
                         7.015: iso_frac2('minor', n14n15)}

    def set_linelists(self, linelists):
        '''
        Sets the desired libraries
        '''
        self.linelists = linelists

    def set_abund(self, element, abund):
        '''
        Sets the abundance of a single element by giving the element symbol and
            the log epsilon abundance.

        element : string, element symbol
        abund : float, new log epsilon abundance value
        '''
        if self.turbodaemon._abund_flag:
            self.turbodaemon.abundances[au.atomic_sym_to_num(element)] = abund

    def get_abund(self, element):
        '''
        Retrieves a given abundance.
        '''
        if self.turbodaemon._abund_flag:
            return self.turbodaemon.abundances[au.atomic_sym_to_num(element)]
        else:
            return

    def get_allabund(self):
        '''
        return all abundances in a dictionary
        '''
        if self.turbodaemon._abund_flag:
            return {au.atomic_num_to_sym(
                elem): abund for elem,
                abund in self.turbodaemon.abundances.items()}
        else:
            return

    def make_atmosphere(self, star=None):
        '''
        Makes the atmos manager and runs the interpolator
        '''

        cwd = os.getcwd()
        abspath = os.path.join(cwd, self.path)
        atmos = AtmosManager(self.teff, self.logg, self.feh,
                             vmicro=self.vmicro, cfe=self.cfe,
                             alphafe=self.alphafe)
        self.atmosname = atmos.interp_atmos(star=star, file_format='ts',
                                            path=abspath)

    def set_isotope(self, iso_dict):
        '''
        Sets a the isotope ratios given in the input dictionary.  iso_dict
            should have keys that give the atomic number and isotope mass,
            and the fractions should sum up to 1 for a given element, e.g.,
            {6.012 : 0.9375, 6.013 : 0.0625} to set 12C and 13C respectivly.

        iso_dict : dictionary of element.isotope_mass as keys and relative
            fractions as the values.
        '''

        for key, frac in iso_dict.items():
            self.isotopes[key] = frac

    def synth(self, wave_range, delta_lambda=0.01, synth_fname=None,
              verbose=False):
        '''
        Makes a synthetic spectrum.

        wave_range : len 2 list, wavelenght range to synthesize in angstroms
        delta_lambda : float, wavelength step to synthesize in angstroms
        synth_fname : string, output name for the synthesized spectrum
        '''
        self.turbodaemon.set_wave(lambda_range=wave_range,
                                  delta_lambda=delta_lambda)
        self.opac_filename = self.turbodaemon.run_babsma(model=self.atmosname,
                                                         vmicro=self.vmicro)

        if self.logg < 3.5:
            sph_flag = True
        else:
            sph_flag = False

        self.synth_fname = self.turbodaemon.run_bsyn(
            sph_flag=sph_flag, opac_filename=self.opac_filename,
            linelists=self.linelists, isotopes=self.isotopes,
            result_filename=synth_fname, verbose=verbose)

    def eqwidth(self, wave_range, delta_lambda=0.01, abund_fname=None,
                verbose=False):
        '''
        Makes a synthetic spectrum.

        wave_range : len 2 list, wavelenght range to synthesize in angstroms
        delta_lambda : float, wavelength step to synthesize in angstroms
        synth_fname : string, output name for the synthesized spectrum
        '''
        self.turbodaemon.set_wave(lambda_range=wave_range,
                                  delta_lambda=delta_lambda)
        self.opac_filename = self.turbodaemon.run_babsma(model=self.atmosname,
                                                         vmicro=self.vmicro)

        if self.logg < 3.5:
            sph_flag = True
        else:
            sph_flag = False

        self.abund_fname = self.turbodaemon.run_eqwidt(
            sph_flag=sph_flag, opac_filename=self.opac_filename,
            linelists=self.linelists, isotopes=self.isotopes,
            result_filename=abund_fname, verbose=verbose)

    def convol(self, profile=None, fwhm=None, vel=None, synth_fname=None,
               convol_fname=None):
        '''
        Convolves the synthetic spectrum with a chosen profile.

        profile : integer (1-4), indicates what broadening profile to use
            (1=exp, 2=gauss, 3=rad-tan, 4=rot)
        fwhm : float, convolutional broadening FWHM in miliAngstroms
        vel : float, convolutional broadening velocity in km/s
        synth_fname : string, input synthetic spectrum
        convol_fname : string, output convolved spectrum
        '''
        if synth_fname is None:
            synth_fname = self.synth_fname

        output_name = self.convoldaemon.run_faltbon(
            synth_fname, profile=profile, fwhm=fwhm, vel=vel,
            result=convol_fname)

        return output_name


def iso_frac2(key, ratio):
    if key == 'major':
        frac = ratio / (ratio + 1.)
    elif key == 'minor':
        frac = 1 / (ratio + 1.)
    return frac
