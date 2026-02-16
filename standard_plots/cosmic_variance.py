import numpy as np

vol = 28.17131365194189 #cMpc

vol_for_test = 310/11.0

muv, xg, yg, zg= np.loadtxt('uv_z10_sample_lagos19.txt', unpack = True, usecols =[0,2,3,4])

n_per_side = 11

mlow = -25 + 5.0 * np.log10(0.677)
mupp = -10 + 5.0 * np.log10(0.677)
dm = 0.5
mbins = np.arange(mlow,mupp,dm)
xlf   = mbins + dm/2.0


uv_lf = np.zeros(shape = (4**3+2, len(xlf)))

print(min(xg), max(xg))
uv_lf = np.zeros(shape = (4**3+2, len(xlf)))
p = 2
H, bins_edges = np.histogram(muv,bins=np.append(mbins,mupp))
uv_lf[0,:] = xlf
uv_lf[1,:] = H

print(np.arange(0, n_per_side, 3))
for x in np.arange(0, n_per_side, 3):
    for y in np.arange(0, n_per_side, 3):
        for z in np.arange(0, n_per_side, 3):
            ind = np.where((xg >= vol_for_test * x) & (xg < vol_for_test * (x+1)) & (yg >= vol_for_test * y) & (yg < vol_for_test * (y+1)) & (zg >= vol_for_test * z) & (zg < vol_for_test * (z+1)))
            vol_patch = vol_for_test**3
            H, bins_edges = np.histogram(muv[ind],bins=np.append(mbins,mupp))
            uv_lf[p,:] = H #/vol_patch/dm
            p = p + 1
            print(x,y,z,p)


np.savetxt('UV_LF_Counts_Lagos19.txt', np.transpose(uv_lf))
uv_lf_final = np.zeros(shape = (2, len(xlf)))

for m in range(0,len(xlf)):
     
    ind = np.where(uv_lf[:,m] != 0)
    uv_lf_final[0,m] = np.mean(uv_lf[ind,m])
    uv_lf_final[1,m] = np.sqrt((np.std(uv_lf[ind,m]))**2 + np.mean(uv_lf[ind,m])) #summing errors in quadrature

print(uv_lf_final[0,:])
uv_lf_final = uv_lf_final / vol_patch/dm

for a,b,c in zip(xlf, uv_lf_final[0,:], uv_lf_final[1,:]):
    print(a,b,c)

