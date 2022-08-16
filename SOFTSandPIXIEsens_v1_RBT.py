import numpy as np
import matplotlib.pyplot as plt
import fisher, sys
import foregrounds as fg
import spectral_distortions as sd
import run_fisher as rs
#sys.path.append("/Users/rbt/Documents/GitHub/submmIFU/codes/")
sys.path.append("/mnt/c/Users/ritob/Documents/GitHub/submmIFU/codes/")

spath = '/mnt/c/Users/ritob/Documents/JPL_stuff/SDLIMwork/'
from matplotlib.pyplot import cm
import Mather_photonNEP12a as NEP
import scipy.constants as scons
ndp = np.float64

c= scons.c
h = scons.h
kB = scons.Boltzmann

## PIXIE AS IS
sdata = np.loadtxt('templates/Sensitivities.dat', dtype=ndp)
fs_pixie = sdata[:, 0] * 1e9
sens_pixie = sdata[:, 1]

def round_to_multiple(number, multiple):
    return multiple * round(number / multiple)


## SOFTS with Nse as variable
def sigB_giveNse(band_details, Time):

    nu_min = band_details['nu_minGHz']*1E9
    nu_max = nu_min*band_details['BWrat']
    BWR = band_details['BWrat']
    #nu_res = band_details['nu_resGHz']*1E9
    Npx = band_details['N_pixels']
    Nse = band_details['Nse']
    #Nse = int(np.ceil((nu_max-nu_min)/nu_res))
    nu_res = (nu_max-nu_min)/Nse

    NEP_phot1 = NEP.photonNEPdifflim(nu_min, nu_max, 2.73, aef=1.0) #This is CMB Tnoise
    NEP_phot2 = NEP.photonNEPdifflim(nu_min, nu_max, 18.1, aef=1E-5) #This is CIB Tnoise
    NEP_det = 10E-18 # 10 ATTO WATTS per square-root(hz)
    NEP_phot = np.sqrt(NEP_phot1**2 +NEP_phot2**2 + NEP_det**2) #Dont include atmosphere for now
    # in making nu_vec we must be aware of resolution
    nu_vec = np.linspace(nu_min, nu_max, Nse) #in Hz
    AOnu = (c/nu_vec)**2

    if (NEP_phot<1E-18):
        print('Warning photon NEP is below 1 aW/rt.Hz!')
    elif (NEP_phot>1E-15):
        print('Warning photon NEP is above 1 fW/rt.Hz!')

    delP = (BWR/np.pi)*NEP_phot/np.sqrt(Time*Npx)
    sigma_B = delP/(AOnu)/nu_res
    sigma_B_beamweighted = sigma_B*(nu_vec/nu_max)

    return nu_vec, sigma_B_beamweighted, Nse

#### main stuff

config = 'config_1'

BWR = 1.52
numinset, numaxset = np.zeros(7), np.zeros(7)
numinset[0] = 60
numinset[1:] = numinset[0]*(BWR**np.linspace(1,6,6))
numaxset[:-1] = numinset[0]*(BWR**np.linspace(1,6,6))
numaxset[-1] = np.max(numinset[0]*(BWR**np.linspace(1,6,6)))*BWR
print('sub-bands min GHz', numinset)
print('sub-bands max GHz',numaxset)

## use this with sigB_giveNse
Nsemin = 6
Nse1 = Nsemin*1
Nse2 = Nsemin*1
Nse3 = Nsemin*1
Nse4 = Nsemin*1
Nse5 = Nsemin*1
Nse6 = Nsemin*1
Nse7 = Nsemin*1

pix1 = round_to_multiple(1E2*0.3,5)
pix2 = round_to_multiple(3E2*0.25,5)
pix3 = round_to_multiple(5E2*0.3,5)
pix4 = round_to_multiple(6E2*0.3,5)
pix5 = round_to_multiple(7E2*0.3,5)
pix6 = round_to_multiple(1.1E3*0.3,5)
pix7 = round_to_multiple(3E3*0.3,5)

tot_pix = pix1+pix2+pix3+pix4+pix5+pix6+pix7

# implictly scale by f_sky = 0.7
Tobs = 3.942E7 # 2.2706e+8 = 86.4 months, 3.942e+7 = 15 months

print('Total pix', tot_pix)

print('total pixel count SOFTS=',tot_pix)

## use this with sigB_giveNse
Bands_list = [{'name':'Band 1','nu_minGHz':numinset[0],'BWrat':BWR,'Nse':Nse1,'N_pixels':pix1},\
      {'name':'Band 2','nu_minGHz':numinset[1],'BWrat':BWR,'Nse':Nse2,'N_pixels':pix2},\
      {'name':'Band 3','nu_minGHz':numinset[2],'BWrat':BWR,'Nse':Nse3,'N_pixels':pix3},\
      {'name':'Band 4','nu_minGHz':numinset[3],'BWrat':BWR,'Nse':Nse4,'N_pixels':pix4},\
      {'name':'Band 5','nu_minGHz':numinset[4],'BWrat':BWR,'Nse':Nse5,'N_pixels':pix5},\
      {'name':'Band 6','nu_minGHz':numinset[5],'BWrat':BWR,'Nse':Nse6,'N_pixels':pix6},\
      {'name':'Band 7','nu_minGHz':numinset[6],'BWrat':BWR,'Nse':Nse7,'N_pixels':pix7},\
      ]



sensitivity_fname= config+"_SOFTS_sens.dat"
f=open("templates/"+sensitivity_fname,"w")

ncolor = len(Bands_list)
color=iter(cm.rainbow(np.linspace(1,0,ncolor+2)))

Nse_set = np.empty(len(Bands_list))

plt.figure(figsize=(8,5))
for bb in range(ncolor):
    b = Bands_list[bb]
    nu_vec_b, sigma_B_b, Nse = sigB_giveNse(b, Tobs)
    Nse_set[bb] = Nse
    for jj in range(len(nu_vec_b)):
        temp1 = nu_vec_b[jj]*1E-9
        temp2 = sigma_B_b[jj]
        tempall = str(temp1)+'\t'+str(temp2)
        f.write("%s\n" % tempall)

    colr=next(color)
    #plt.plot(nu_vec_b_beam*1E-9, sigma_B_b_beam, alpha=0.8, color=colr, label=b['name'] +' Npx='+str(b['N_pixels']))
    #plt.plot(nu_vec_b_beam*1E-9, sigma_B_b_beam, marker='o', alpha=0.8, color=colr)
    plt.plot(nu_vec_b*1E-9, sigma_B_b, marker='*', alpha=0.8, color=colr, label=b['name'] +' Npx='+str(b['N_pixels']))
plt.loglog(fs_pixie*1E-9,sens_pixie,'--k', label='PIXIE 15 months TOTAL sensitivty')
plt.ylabel('W/m^2/Sr/Hz')
plt.xlabel('GHz')
#plt.legend(title='SOFTS' , bbox_to_anchor=(1,0.5))
#plt.ylim(1E-26, 1E-21)
plt.xlim(20, 2500)
plt.tight_layout()
plt.savefig(spath+'NESB_wbeamcorr_'+config+'.pdf', bbox_inches='tight')
plt.show()
f.close()
#
print(config+' Nse_set', Nse_set)

rs.rbtfooling2(sensitivity_fname)
