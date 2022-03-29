from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
import astropy.units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from photutils import DAOStarFinder, CircularAperture, CircularAnnulus, aperture_photometry
import numpy as np
import matplotlib.pyplot as plt

def get_inst_mag(pos, data, r_c=7, r_in=8, r_out=10):
    star = CircularAperture(pos,r=r_c)
    bkgd = CirculurAnnulus(pos, r_in=r_in, r_out=r_out)
    apers = [star,bkgd]
    phot_table = aperture_photometry(data,apers)
    bkgd_flux = phot_table['aperture_sum_1']/bkgd.area*star.area
    flux = phot_table['aperture_sum_0']- bkgd_flux
    mag = -2.5*np.log(flux)
    return mag

### get standard positions

hdu2 = fits.open("wcs_ngc6633_data.fits")
wcs_head = hdu2[0].header
wcs = WCS(wcs_head)

ra,dec,V,B = np.loadtxt("apass_ngc6633.csv",unpack=True,skiprows=1,delimiter=',')

x_pix, y_pix = WCS.all_world2pix(ra*u.deg, dec*u.deg, 1)
positions = np.transpose([x_pix, y_pix])

###get V mags

hdu = fits.open("master_V.fit")
Vdata = np.array(hdu[0].data,dtype=np.int32)

inst_v_mag = get_inst_mag(positions,Vdata)

###calculate b mags

hdu=fits.open("master_B.fit")
Bdata=np.array(hdu[0].data, dtype=np.int32)

inst_b_mag = get_inst_mag(positions, Bdata)

###trim out  bad data for fitting

inst_bv=inst_b_mag-inst_v_mag
app_bv=B-V
blist=B-inst_b_mag

mask=np.isnan(inst_bv)
inst_bv=inst_bv[~mask]
app_bv=app_bv[~mask]
blist=blist[~mask]

mask=(np.abs(inst_bv)<=2)
inst_bv=inst_bv[mask]
app_bv=app_bv[mask]
blist=blist[mask]

mask=(inst_bv>=-0.5)
inst_bv=inst_bv[mask]
app_bv=app_bv[mask]
blist=blist[mask]

###calculate Tbv and Cbv

fit, cov=np.polyfit(inst_bv, app_bv, 1, cov=True)
Tbv,Cbv=fit[0],fit[1]
print(Tbv,Cbv,np.sqrt(np.diag(cov)))


###plot Tbv and Cbv fit to check

plt.scatter(inst_bv,app_bv)
xlist=np.linspace(np.min(inst_bv)*0.8,np.max(inst_bv)*1.1,100)
ylist=Tbv*xlist+Cbv
plt.plot(xlist,ylist,'k', lw=4)
plt.ylim(plt.xlim())
plt.xlabel("B-V")
plt.ylabel("V")
plt.show




#calculate Tv and Cv

fit, cov=np.polyfit(app_bv,blist,1,cov=True)
Tb,Cb=fit[0],fit[1]
print(Tb,Cb,np.sqrt(np.diag(cov)))

#plot Tb and Cb to check
plt.scatter(app_bv, blist)
ylist=Tb*xlist+Cb
plt.plot(xlist,ylist,'k',lw=4)
plt.xlabel("B-V")
plt.ylabel("V")
plt.show()






###get all stars

mean,median,std= sigma_clipped_stats(Vdata,sigma=3.0)
daofind=DAOStarFinder(fwhm=3.0, threshold=3.*std)
Vsources = daofind(Vdata-median)

mean,median,std= sigma_clipped_stats(Bdata,sigma=3.0)
daofind=DAOStarFinder(fwhm=3.0, threshold=3.*std)
Bsources = daofind(Bdata-median)

sources=Table(names=("x","y"))
dist_sq=0.5
for vstar in Vsources:
    for bstar in Bsources:
        if((bstar["xcentroid"]-vstar['xcentroid'])**2+(bstar["ycentroid"]-vstar['ycentroid'])**2)<=dist_sq:
            sources.add_row([bstar['xcentroid'],bstar['ycentroid']])
print(len(sources))

positions=np.transpose([sources['x'],sources['y']])

all_inst_b_mag = get_inst_mag(positions,Bdata)
all_inst_v_mag = get_inst_mag(positions, Vdata)

all_BV = Tbv*(all_inst_b_mag - all_inst_v_mag) +Cbv
all_B = Tb*all_BV + Cb + all_inst_v_mag
all_V=all_B-all_BV

###find distance modulus

isochrone=np.loadtxt("#isochrone file")
iso=isochrone[isochrone[:,0]==4]
plt.plot(iso[:,7], iso[:,5],label = 4)
plt.plot(all_BV,all_V, ',')
plt.xlim(-0.5,2)
plt.ylim(reversed(plt.ylim()))
plt.show()
