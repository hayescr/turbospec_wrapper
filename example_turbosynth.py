from turbosynth import TurboSynth

teff = 4000
logg = 0.5
feh = -1.5
vmicro = 2.5
cfe = 0.0
alphafe = 0.0

teff = 4286
logg = 1.64
feh = -0.53
vmicro = 1.25
cfe = 0.16
alphafe = 0.49
res = 43.0
vmacro = 4.06
vrot = 3.80

profiles = [2, 2, 4]
broads = [{'fwhm': res}, {'vel': vmacro}, {'vel': vrot}]

star = 'arcturus_new_fill0_mol_notivzr'

abundances = {'C': 8.03, 'N': 7.62, 'O': 8.63, 'Na': 5.76, 'Mg': 7.38, 'Al': 6.19, 'Si': 7.15, 'Ca': 5.9, 'Sc': 2.68, 'Ti': 4.59, 'V': 3.68, 'Cr': 5.07, 'Mn': 4.66, 'Fe': 6.93, 'Co': 4.49, 'Ni': 5.77}

rel_path = '../Turbospectrum2019/COM-v19.1/'
linelists = ['DATA/Hlinedata',
'bacchus_linelists/atoms_4200-9200.list',
'bacchus_linelists/12C12C.list',
'bacchus_linelists/12C13C.list',
'bacchus_linelists/13C13C.list',
'bacchus_linelists/12CH.list',
'bacchus_linelists/13CH.list',
'bacchus_linelists/MgH.list',
'bacchus_linelists/12CN.list',
'bacchus_linelists/13CN.list']

linelists = ['DATA/Hlinedata',
'new_linelists/ges_atomic_linelist.turbo.txt',
'bacchus_linelists/12C12C.list',
'bacchus_linelists/12C13C.list',
'bacchus_linelists/13C13C.list',
'bacchus_linelists/12CH.list',
'bacchus_linelists/13CH.list',
'bacchus_linelists/MgH.list',
'bacchus_linelists/12CN.list',
'bacchus_linelists/13CN.list']

linelists = ['DATA/Hlinedata',
'new_linelists/ges_atomic_linelist.turbo.txt',
'new_linelists/12C1H_ges_linelist.turbo.txt',
'new_linelists/13C1H_ges_linelist.turbo.txt',
'new_linelists/12C12C_ges_linelist.turbo.txt',
'new_linelists/12C13C_ges_linelist.turbo.txt',
'new_linelists/13C13C_ges_linelist.turbo.txt',
'new_linelists/12C14N_ges_linelist.turbo.txt',
'new_linelists/12C15N_ges_linelist.turbo.txt',
'new_linelists/13C14N_ges_linelist.turbo.txt',
'new_linelists/14N1H_ges_linelist.turbo.txt',
'new_linelists/16O1H_ges_linelist.turbo.txt',
'new_linelists/24Mg1H_ges_linelist.turbo.txt',
'new_linelists/25Mg1H_ges_linelist.turbo.txt',
'new_linelists/26Mg1H_ges_linelist.turbo.txt',
'new_linelists/28Si1H_ges_linelist.turbo.txt',
'new_linelists/40Ca1H_ges_linelist.turbo.txt',
'new_linelists/56Fe1H_ges_linelist.turbo.txt',
'new_linelists/46Ti16O_ges_linelist.turbo.txt','new_linelists/47Ti16O_ges_linelist.turbo.txt','new_linelists/48Ti16O_ges_linelist.turbo.txt','new_linelists/49Ti16O_ges_linelist.turbo.txt','new_linelists/50Ti16O_ges_linelist.turbo.txt','new_linelists/51V16O_ges_linelist.turbo.txt',
'new_linelists/90Zr16O_ges_linelist.turbo.txt','new_linelists/91Zr16O_ges_linelist.turbo.txt','new_linelists/92Zr16O_ges_linelist.turbo.txt','new_linelists/94Zr16O_ges_linelist.turbo.txt','new_linelists/96Zr16O_ges_linelist.turbo.txt',]

linelists = ['DATA/Hlinedata',
'new_linelists/ges_atomic_linelist.turbo.txt',
'new_linelists/12C1H_ges_linelist.turbo.txt',
'new_linelists/13C1H_ges_linelist.turbo.txt',
'new_linelists/12C12C_ges_linelist.turbo.txt',
'new_linelists/12C13C_ges_linelist.turbo.txt',
'new_linelists/13C13C_ges_linelist.turbo.txt',
'new_linelists/12C14N_ges_linelist.turbo.txt',
'new_linelists/12C15N_ges_linelist.turbo.txt',
'new_linelists/13C14N_ges_linelist.turbo.txt',
'new_linelists/14N1H_ges_linelist.turbo.txt',
'new_linelists/16O1H_ges_linelist.turbo.txt',
'new_linelists/24Mg1H_ges_linelist.turbo.txt',
'new_linelists/25Mg1H_ges_linelist.turbo.txt',
'new_linelists/26Mg1H_ges_linelist.turbo.txt',
'new_linelists/28Si1H_ges_linelist.turbo.txt',
'new_linelists/40Ca1H_ges_linelist.turbo.txt',
'new_linelists/56Fe1H_ges_linelist.turbo.txt']


#linelists = ['DATA/Hlinedata',
#'bacchus_linelists/atoms_4200-9200.list']

linelists = [f'{rel_path}{linelist}' for linelist in linelists]

wave_range = [5000, 5500]


test = TurboSynth(path='turbosynth_test', turbopath='../Turbospectrum2019')
test.init_params(teff, logg, feh, vmicro=vmicro, cfe=cfe, alphafe=alphafe)
test.init_abunds(abunds=abundances)

#c_iso = {6.012: 0.80, 6.013: 0.20}
#test.set_isotope(c_iso)
test.set_linelists(linelists)

print(f'The abundance of C should be {test.get_abund("C")}')
print(f'The abundance dict should be {test.get_allabund()}')

test.make_atmosphere(star=star)
test.synth(wave_range)

for i, (profile, kwargs) in enumerate(zip(profiles, broads)):
    if i == 0:
        fname=test.synth_fname
    else:
        fname=convol
    convol=test.convol(synth_fname=fname, profile=profile, **kwargs)