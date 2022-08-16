import numpy as np
import matplotlib.pyplot as plt
import fisher, sys
import foregrounds as fg
import spectral_distortions as sd
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

# #let's make sure that our NEP caclulator is correct
# # RJ approximation is
# def NEP_RJfunc(nu0, Delta_nu, Ts):
#     x0 = h*nu0/(kB*Ts)
#     if x0 > 0.5:
#         print('RJ approx is inaccurate... going ahead anyway')
#     Popt = 2.0*kB*Ts*Delta_nu
#     NEPsquared = 2*h*nu0*Popt + (Popt**2)/Delta_nu
#     return np.sqrt(NEPsquared)
#
#
# for ii in range(30,87,2):
#     nu0_try = ii*1E9
#     Delta_nu_try1 = 10E9
#     Delta_nu_try2 = 10E9
#     T_star1 = 1000.1#K
#     T_star2 = 300.1
#     # for a 300K Blackbody the peak is at 30 THz.
#
#     NEP_rbt_1=NEP.photonNEPdifflim(nu0_try-0.5*Delta_nu_try1 , nu0_try+0.5*Delta_nu_try1 , T_star1)
#     NEP_RJ_1 = NEP_RJfunc(nu0_try , Delta_nu_try1, T_star1)
#
#     print(NEP_rbt_1, NEP_RJ_1)
#
#     NEP_rbt_2=NEP.photonNEPdifflim(nu0_try-0.5*Delta_nu_try2 , nu0_try+0.5*Delta_nu_try2 , T_star2)
#     NEP_RJ_2 = NEP_RJfunc(nu0_try , Delta_nu_try2, T_star2)
#     plt.scatter(ii, NEP_rbt_1/NEP_RJ_1, color= 'r',marker='d')
#     plt.scatter(ii, NEP_rbt_2/NEP_RJ_2, color= 'b',marker='*')
#
# plt.ylabel('RBT/RJ')
# plt.show()

# nu = np.linspace(1E1, 3E3, 99)*1E9
# bbspec = (1E-5)* 2*h*(nu**3)/(c**2)/(np.exp(h*nu/(kB*18.8))-1)
# plt.loglog(nu*1E-9, bbspec*1E20)
# plt.xlabel('GHz')
# plt.ylabel('MJy/sr')
# plt.title('CIB as BB - crude approx.')
# plt.show()

## PIXIE AS IS
sdata = np.loadtxt('templates/Sensitivities.dat', dtype=ndp)
fs_pixie = sdata[:, 0] * 1e9
sens_pixie = sdata[:, 1]

def round_to_multiple(number, multiple):
    return multiple * round(number / multiple)

## SOFTS with resolution as variable
def sigB(band_details, Time):

    nu_min = band_details['nu_minGHz']*1E9
    nu_max = nu_min*band_details['BWrat']
    nu_res = band_details['nu_resGHz']*1E9
    Npx = band_details['N_pixels']
    Nse = int(np.ceil((nu_max-nu_min)/nu_res))

    NEP_phot1 = NEP.photonNEPdifflim(nu_min, nu_max, 2.8, aef=1.0) #This is CMB Tnoise
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

    delP = 2*Nse*NEP_phot/np.sqrt(Time*Npx)
    sigma_B = delP/(AOnu)/nu_res

    return nu_vec, sigma_B, Nse

## SOFTS with Nse as variable
def sigB_giveNse(band_details, Time):

    nu_min = band_details['nu_minGHz']*1E9
    nu_max = nu_min*band_details['BWrat']
    #nu_res = band_details['nu_resGHz']*1E9
    Npx = band_details['N_pixels']
    Nse = band_details['Nse']
    #Nse = int(np.ceil((nu_max-nu_min)/nu_res))
    nu_res = (nu_max-nu_min)/Nse

    NEP_phot1 = NEP.photonNEPdifflim(nu_min, nu_max, 2.8, aef=1.0) #This is CMB Tnoise
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

    delP = 2*Nse*NEP_phot/np.sqrt(Time*Npx)
    sigma_B = delP/(AOnu)/nu_res

    return nu_vec, sigma_B, Nse

#### main stuff
## use this with sigB
res = 5.0 # GHz
res1=res*1
res2=res*1.5
res3=res*2
res4=res*3
res5=res*4
res6=res*5
res7=res*6
## use this with sigB_giveNse
Nsemin = 2
Nse1 = Nsemin*1
Nse2 = Nsemin*1
Nse3 = Nsemin*1
Nse4 = Nsemin*1
Nse5 = Nsemin*1
Nse6 = Nsemin*1
Nse7 = Nsemin*1


print(80*(1.5**np.linspace(0,6,7)))

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
Omega_PIXIE = 4*np.pi*0.7
theta_b_SOFTS = 0.5* np.pi/180.0 # 0.5 degree
Omega_SOFTS_beam = 4*np.pi*(np.sin(0.5*theta_b_SOFTS)**2)
N_beams = Omega_PIXIE/Omega_SOFTS_beam
Tobs_beam = Tobs/N_beams

print('N_beams', N_beams)
print('Total pix', tot_pix)

print('total pixel count SOFTS=',tot_pix)

# ## use this with sigB
# Bands_list = [{'name':'Band 1','nu_minGHz':80,'BWrat':1.5,'nu_resGHz':res1,'N_pixels':pix1},\
#       {'name':'Band 2','nu_minGHz':120,'BWrat':1.5,'nu_resGHz':res2,'N_pixels':pix2},\
#       {'name':'Band 3','nu_minGHz':180,'BWrat':1.5,'nu_resGHz':res3,'N_pixels':pix3},\
#       {'name':'Band 4','nu_minGHz':270,'BWrat':1.5,'nu_resGHz':res4,'N_pixels':pix4},\
#       {'name':'Band 5','nu_minGHz':405,'BWrat':1.5,'nu_resGHz':res5,'N_pixels':pix5},\
#       {'name':'Band 6','nu_minGHz':607.5,'BWrat':1.5,'nu_resGHz':res6,'N_pixels':pix6},\
#       {'name':'Band 7','nu_minGHz':911.25,'BWrat':1.5,'nu_resGHz':res7,'N_pixels':pix7},\
#       ]
## use this with sigB_giveNse
Bands_list = [{'name':'Band 1','nu_minGHz':80,'BWrat':1.5,'Nse':Nse1,'N_pixels':pix1},\
      {'name':'Band 2','nu_minGHz':120,'BWrat':1.5,'Nse':Nse2,'N_pixels':pix2},\
      {'name':'Band 3','nu_minGHz':180,'BWrat':1.5,'Nse':Nse3,'N_pixels':pix3},\
      {'name':'Band 4','nu_minGHz':270,'BWrat':1.5,'Nse':Nse4,'N_pixels':pix4},\
      {'name':'Band 5','nu_minGHz':405,'BWrat':1.5,'Nse':Nse5,'N_pixels':pix5},\
      {'name':'Band 6','nu_minGHz':607.5,'BWrat':1.5,'Nse':Nse6,'N_pixels':pix6},\
      {'name':'Band 7','nu_minGHz':911.25,'BWrat':1.5,'Nse':Nse7,'N_pixels':pix7},\
      ]


ncolor = len(Bands_list)
color=iter(cm.rainbow(np.linspace(1,0,ncolor+2)))

Nse_set = np.empty(len(Bands_list))

f=open("templates/config2_SOFTS_sens.dat","w")

plt.figure(figsize=(8,5))
for bb in range(ncolor):
    b = Bands_list[bb]
    nu_vec_b_beam, sigma_B_b_beam, Nse = sigB_giveNse(b, Tobs_beam)
    nu_vec_b, sigma_B_b, Nse = sigB_giveNse(b, Tobs_beam*N_beams)
    Nse_set[bb] = Nse
    for jj in range(len(nu_vec_b)):
        temp1 = nu_vec_b[jj]*1E-9
        temp2 = sigma_B_b[jj]
        tempall = str(temp1)+'\t'+str(temp2)
        f.write("%s\n" % tempall)

    colr=next(color)
    plt.plot(nu_vec_b_beam*1E-9, sigma_B_b_beam, alpha=0.8, color=colr, label=b['name'] +' Npx='+str(b['N_pixels']))
    plt.plot(nu_vec_b_beam*1E-9, sigma_B_b_beam, marker='o', alpha=0.8, color=colr)
    plt.plot(nu_vec_b*1E-9, sigma_B_b, marker='*', alpha=0.8, color=colr)
plt.loglog(fs_pixie*1E-9,sens_pixie,'--k', label='PIXIE 15 months TOTAL sensitivty')
plt.ylabel('W/m^2/Sr/Hz')
plt.xlabel('GHz')
plt.legend(title='SOFTS (colored markers):\n single-beam (o), all beams (*)' , bbox_to_anchor=(1,0.5))
plt.ylim(1E-26, 1E-21)
plt.xlim(50, 2500)
plt.tight_layout()
plt.savefig(spath+'NESB_config2.pdf', bbox_inches='tight')
plt.show()
f.close()
#
# plt.figure()
# plt.scatter(range(len(Bands_list)), Nse_set)
# plt.show()
print('Nse_set', Nse_set)
