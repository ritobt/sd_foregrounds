import numpy as np
import matplotlib.pyplot as plt
import fisher, fisher_SOFTS, sys
import foregrounds as fg
import spectral_distortions as sd
sys.path.append("/mnt/c/Users/ritob/Documents/GitHub/submmIFU/codes/")
import myellp as ellp


def sens_vs_dnu(fmax=61, fstep=0.1, sens=[1., 0.1, 0.01, 0.001]):
    dnu = np.arange(1, fmax, fstep) * 1.e9
    fish = fisher.FisherEstimation()
    args = fish.args
    N = len(args)
    data = {}
    for sen in sens:
        print "on sens ", sen
        data[sen] = {}
        for arg in args:
            data[sen][arg] = []
        for nu in dnu:
            scale = sen * (15.e9 / nu)
            fish = fisher.FisherEstimation(fstep=nu, mult=scale)
            fish.run_fisher_calculation()
            for arg in args:
                data[sen][arg].append(fish.errors[arg])
    x = {}
    x['sens'] = sens
    x['dnu'] = dnu
    x['data'] = data
    np.save('datatest', x)
    return

def sens_vs_nbins(sens=[1., 0.1, 0.01, 0.001]):
    nbins = np.arange(50, 601)[::-1]
    dnu = 3.e12 / nbins
    fish = fisher.FisherEstimation()
    args = fish.args
    N = len(args)
    data = {}
    for sen in sens:
        print "on sens ", sen
        data[sen] = {}
        for arg in args:
            data[sen][arg] = []
        for nu in dnu:
            scale = sen * (15.e9 / nu)
            #ps = {'As':0.1, 'alps': 0.1}
            #ps = {'As':0.01, 'alps': 0.01}
            ps = {}
            fish = fisher.FisherEstimation(fstep=nu, mult=scale, priors=ps)
            fish.run_fisher_calculation()
            for arg in args:
                data[sen][arg].append(fish.errors[arg])
    x = {}
    x['sens'] = sens
    x['dnu'] = dnu
    x['data'] = data
    np.save('senscalc_nbins_0p', x)
    return


def drop_vs_dnu(fmax=61, fstep=0.1, drops=[0, 1, 2]):
    dnu = np.arange(1, fmax, fstep) * 1.e9
    fish = fisher.FisherEstimation()
    args = fish.args
    N = len(args)
    data = {}
    for drop in drops:
        print "on drop ", drop
        data[drop] = {}
        for arg in args:
            data[drop][arg] = []
        for k, nu in enumerate(dnu):
            if k % 100 == 0:
                print "on k ", k
            scale = 15.e9 / nu
            ps = {'As':0.01, 'alps': 0.01}
            #ps = {'As':0.1, 'alps': 0.1}
            fish = fisher.FisherEstimation(fstep=nu, mult=scale, drop=drop, priors=ps)
            fish.run_fisher_calculation()
            for arg in args:
                data[drop][arg].append(fish.errors[arg])
    x = {}
    x['drops'] = drops
    x['dnu'] = dnu
    x['data'] = data
    np.save('fullcalc_1p_drop012_coarse', x)
    return

def drop_vs_dnu_nomu(fmax=61, fstep=0.1, drops=[0, 1, 2]):
    dnu = np.arange(1, fmax, fstep) * 1.e9
    fish = fisher.FisherEstimation()
    fish.set_signals(fncs=[sd.DeltaI_reltSZ_2param_yweight, sd.DeltaI_DeltaT,
                           fg.thermal_dust_rad, fg.cib_rad, fg.jens_freefree_rad,
                           fg.jens_synch_rad, fg.spinning_dust, fg.co_rad])
    args = fish.args
    N = len(args)
    data = {}
    for drop in drops:
        print "on drop ", drop
        data[drop] = {}
        for arg in args:
            data[drop][arg] = []
        for k, nu in enumerate(dnu):
            if k % 100 == 0:
                print "on k ", k
            scale = 15.e9 / nu
            fish = fisher.FisherEstimation(fstep=nu, mult=scale, drop=drop)
            fish.set_signals(fncs=[sd.DeltaI_reltSZ_2param_yweight, sd.DeltaI_DeltaT,
                                   fg.thermal_dust_rad, fg.cib_rad, fg.jens_freefree_rad,
                                   fg.jens_synch_rad, fg.spinning_dust, fg.co_rad])
            fish.run_fisher_calculation()
            for arg in args:
                data[drop][arg].append(fish.errors[arg])
    x = {}
    x['drops'] = drops
    x['dnu'] = dnu
    x['data'] = data
    np.save('fullcalc_10p_drop012_nomu', x)
    return

def drop_vs_nbin(drops=[0, 1, 2]):
    nbins = np.arange(50, 601)[::-1]
    dnu = 3.e12 / nbins
    fish = fisher.FisherEstimation()
    args = fish.args
    N = len(args)
    data = {}
    for drop in drops:
        print "on drop ", drop
        data[drop] = {}
        for arg in args:
            data[drop][arg] = []
        for k, nu in enumerate(dnu):
            if k % 100 == 0:
                print "on k ", k
            scale = 15.e9 / nu
            #ps = {'As':0.01, 'alps': 0.01}
            ps = {'As':0.1, 'alps': 0.1}
            fish = fisher.FisherEstimation(fstep=nu, mult=scale, drop=drop, priors=ps)
            fish.run_fisher_calculation()
            for arg in args:
                data[drop][arg].append(fish.errors[arg])
    x = {}
    x['drops'] = drops
    x['dnu'] = dnu
    x['data'] = data
    np.save('fullcalc_10p_drop012_nbins_onefifthsynch', x)
    return

def drop_vs_nbin_nomu(drops=[0, 1, 2]):
    nbins = np.arange(50, 601)[::-1]
    dnu = 3.e12 / nbins
    fish = fisher.FisherEstimation()
    fish.set_signals(fncs=[sd.DeltaI_reltSZ_2param_yweight, sd.DeltaI_DeltaT,
                           fg.thermal_dust_rad, fg.cib_rad, fg.jens_freefree_rad,
                           fg.jens_synch_rad, fg.spinning_dust, fg.co_rad])
    args = fish.args
    N = len(args)
    data = {}
    for drop in drops:
        print "on drop ", drop
        data[drop] = {}
        for arg in args:
            data[drop][arg] = []
        for k, nu in enumerate(dnu):
            if k % 100 == 0:
                print "on k ", k
            scale = 15.e9 / nu
            fish = fisher.FisherEstimation(fstep=nu, mult=scale, drop=drop)
            fish.set_signals(fncs=[sd.DeltaI_reltSZ_2param_yweight, sd.DeltaI_DeltaT,
                                   fg.thermal_dust_rad, fg.cib_rad, fg.jens_freefree_rad,
                                   fg.jens_synch_rad, fg.spinning_dust, fg.co_rad])
            fish.run_fisher_calculation()
            for arg in args:
                data[drop][arg].append(fish.errors[arg])
    x = {}
    x['drops'] = drops
    x['dnu'] = dnu
    x['data'] = data
    np.save('fullcalc_10p_drop012_nbins_nomu', x)
    return

def rbtfooling1():
    # dnu = 9E9
    fish = fisher.FisherEstimation()
    #fish = fisher_SOFTS.FisherEstimation()
    args = fish.args
    N = len(args)
    data = {}
    for arg in args:
        data[arg] = []
    fish = fisher_SOFTS.FisherEstimation(bandpass=False,doCO=True,doSOFTS=False)
    fish.run_fisher_calculation()
    print(np.shape(fish.cov))
    #print(np.shape(fish.args))
    #print(fish.args)

###
    Matrix_Cov = fish.cov
    N_templates=len(fish.cov)
    labs_templates = fish.args
    plt.figure()
    plt.bar(range(N_templates), np.sqrt(np.diagonal(Matrix_Cov)), align='center')
    plt.yscale('log')
    plt.ylabel('sigma = sqrt(diag(Matrix_Cov))')
    plt.xticks(range(N_templates), labs_templates, size='small')
    plt.show()
    count = 0
    plt.figure(figsize=(3*N_templates,3*N_templates))
    for mm in range(N_templates):
        for nn in range(N_templates):
            count+=1
            xin, yin = 0, 0
            xc, yc = 0, 0
            C = np.array([[0., 0.],[0., 0.]])
            C[0,0] = Matrix_Cov[mm, mm]
            C[0,1] = Matrix_Cov[mm, nn]
            C[1,0] = Matrix_Cov[nn, mm]
            C[1,1] = Matrix_Cov[nn, nn]
            #print(C)

            if mm > nn:
                plt.subplot(N_templates, N_templates, count)
                stuff= ellp.myellp(x_cent=xc, y_cent=yc, cov=C, mass_level=0.95)
                plt.fill(stuff[0],stuff[1],color='c', alpha = 0.2)
                stuff= ellp.myellp(x_cent=xc, y_cent=yc, cov=C, mass_level=0.68)
                plt.fill(stuff[0],stuff[1],color='b', alpha = 0.2)
                plt.plot(xc,yc, marker='+', ms = 12, color='purple')
                plt.plot(xin,yin, marker='*', ms = 10, color='orange')
                plt.grid('True', alpha = 0.5)
                plt.xlabel(labs_templates[mm])
                plt.ylabel(labs_templates[nn])
    plt.tight_layout()
    plt.close
    # plt.savefig('rbtfooling1_FisherEllipses.png')
    # plt.show()

###


    for arg in args:
        data[arg].append(fish.errors[arg])
    x = {}
    x['data'] = data
    np.save('datatest_rbtfooling1', x)
    return

def rbtfooling2(softs_sens_file):
    spath = '/mnt/c/Users/ritob/Documents/JPL_stuff/SDLIMwork/'

## regular PIXIE
    fish = fisher.FisherEstimation()
    #fish = fisher_SOFTS.FisherEstimation()
    args = fish.args
    N = len(args)
    data = {}
    for arg in args:
        data[arg] = []
    fish = fisher_SOFTS.FisherEstimation(bandpass=False,doCO=True,doSOFTS=False)
    fish.run_fisher_calculation()
    print(np.shape(fish.cov))
    #print(np.shape(fish.args))
    #print(fish.args)

    Matrix_Cov = fish.cov
    print(np.sqrt(np.diagonal(Matrix_Cov)))
    N_templates=len(fish.cov)
    labs_templates = fish.args
    plt.figure(figsize=(15,5))
    plt.bar(range(N_templates), np.sqrt(np.diagonal(Matrix_Cov)), alpha =0.5, align='center', label='regular PIXIE')

# for SOFTS
    fish = fisher_SOFTS.FisherEstimation()
    args = fish.args
    N = len(args)
    data = {}
    for arg in args:
        data[arg] = []
    fish = fisher_SOFTS.FisherEstimation(bandpass=False,doCO=True,doSOFTS=True, SOFTS_sensitivity_file = softs_sens_file)
    fish.run_fisher_calculation()
    print(np.shape(fish.cov))
    #print(np.shape(fish.args))
    #print(fish.args)

    Matrix_Cov_SOFTS = fish.cov
    print(np.sqrt(np.diagonal(Matrix_Cov_SOFTS)))
    N_templates=len(fish.cov)
    labs_templates = fish.args

    tf = softs_sens_file.split('_SOFTS_sens.dat')
    plt.bar(range(N_templates), np.sqrt(np.diagonal(Matrix_Cov_SOFTS)), alpha =0.5, align='center', label='SOFTS '+tf[0])
    plt.yscale('log')
    plt.ylabel('sigma = sqrt(diag(Matrix_Cov))')
    plt.xticks(range(N_templates), labs_templates, size='small')
    plt.legend()

    plt.savefig(spath+'fisher_Comparion_'+tf[0]+'.pdf')
    plt.show()

    print('Ratio of stds: PIXIE/SOFTS ...\n')
    print(np.sqrt(np.diagonal(Matrix_Cov))/np.sqrt(np.diagonal(Matrix_Cov_SOFTS)))

    return

def rbtfooling3():
    spath = '/mnt/c/Users/ritob/Documents/JPL_stuff/SDLIMwork/'

## regular PIXIE testing
    fish = fisher.FisherEstimation()
    #fish = fisher_SOFTS.FisherEstimation()
    args = fish.args
    N = len(args)
    data = {}
    for arg in args:
        data[arg] = []
    fish = fisher_SOFTS.FisherEstimation(bandpass=False,doCO=True,doSOFTS=False)
    fish.run_fisher_calculation()
    print(np.shape(fish.cov))
    #print(np.shape(fish.args))
    #print(fish.args)

    Matrix_Cov = fish.cov
    print(np.sqrt(np.diagonal(Matrix_Cov)))
    N_templates=len(fish.cov)
    labs_templates = fish.args
    plt.figure(figsize=(12,5))
    plt.bar(range(N_templates), np.sqrt(np.diagonal(Matrix_Cov)), alpha =0.5, align='center', label='regular PIXIE no BP')

## regular PIXIE testing
    fish = fisher.FisherEstimation()
    #fish = fisher_SOFTS.FisherEstimation()
    args = fish.args
    N = len(args)
    data = {}
    for arg in args:
        data[arg] = []
    fish = fisher_SOFTS.FisherEstimation(bandpass=True,doCO=True,doSOFTS=False)
    fish.run_fisher_calculation()
    print(np.shape(fish.cov))
    #print(np.shape(fish.args))
    #print(fish.args)

    Matrix_Cov = fish.cov
    print(np.sqrt(np.diagonal(Matrix_Cov)))
    N_templates=len(fish.cov)
    labs_templates = fish.args

    plt.bar(range(N_templates), np.sqrt(np.diagonal(Matrix_Cov)), alpha =0.5, align='center', label='regular PIXIE yes BP')


    plt.yscale('log')
    plt.ylabel('sigma = sqrt(diag(Matrix_Cov))')
    plt.xticks(range(N_templates), labs_templates, size='small')
    plt.legend()

    plt.show()
    ## OK seems like yes / no bnadpass makes very little difference
    return


#drop_vs_nbin_nomu()
#drop_vs_nbin()
#drop_vs_dnu(fmax=61, fstep=0.1, drops=[0, 1, 2])
#drop_vs_dnu_nomu(fmax=61, fstep=0.1, drops=[0, 1, 2])
#sens_vs_nbins()
#rbtfooling2()
