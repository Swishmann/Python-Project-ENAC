#!/usr/bin/env python
# -*-coding:Latin-1 -*
#ppack_gnu_435
#export PYTHONPATH='/projects/MTS-mdvs_develop/to28780/matplotlib-1.2.0/build/lib.linux-x86_64-2.6':$PYTHONPATH
#export LD_LIBRARY_PATH=/projects/MTS-mdvs_develop/to28780/matplotlib-1.2.0/build/lib.linux-x86_64-2.6/matplotlib:$LD_LIBRARY_PATH
#export PYTHONPATH='/projects/MTS-PERFOS/DONNEES_BGV/to22716/06_MODULES/xlutils-1.7.1':$PYTHONPATH
#export PYTHONPATH='/projects/MTS-PERFOS/DONNEES_BGV/to22716/06_MODULES/xlrd-0.9.3':$PYTHONPATH
#. init_siremi test
import os
import glob
import sys
sys.path.append("/projects/MTS-PERFOS/DONNEES_BGV/to22716/06_MODULES")
sys.path.insert(1,'/projects/MASHARE/TOOLS/EXPLOIT/UTIL/modules')
from  reader_perfo_files import read_file
from vmu_plot import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties
import matplotlib.colors as colors
from itertools import *
from perfo_bga import *
from numpy.polynomial import Polynomial as P
import numpy.polynomial.polynomial as poly
from filters import *
import Aladyn
import AladynPackage
import Adn.Objects.FtObject as Adn

#--------------------------------------------------------------
fic_adn = sys.argv[1]
teta_target = float(sys.argv[2])

#--------------------------------------------------------------
# liste des tps caracteristique a extraire du fichier adn
#--------------------------------------------------------------
liste_tps_carac_adn = ['Instant de Lift Off','Instant de mise en rotation','Instant de passage des 35ft','Instant de mise en poussee']

#--------------------------------------------------------------
# conversion au format CSV du fichier ADN-EDAS
#--------------------------------------------------------------
adn = adn_edas(path=fic_adn)    
#replace_string_file(fic_adn,outfile='None',replacements={'nan':'0.0'})
#os.system('sed \'s/nan/0.0/g\' '+ fic_adn+' > toto')
#os.system('mv toto '+fic_adn)
adn.save(path=fic_adn+'_tmp.csv', fmt='csv') 
dt = adn.get_parameter('DT', datatype="[Donnees EV]")





#--------------------------------------------------------------
# analyse rotation
#--------------------------------------------------------------
t_vr_eva = adn.get_parameter('Instant de mise en rotation', datatype="[Caracteristiques de l'essai]")[0] 
adn.add_parameter('t_vr_eva', [t_vr_eva]*len(dt), '', datatype="[Donnees EV]")
i_vr_eva = [i for i in range(len(dt)) if abs(dt[i] - t_vr_eva)<0.001 ][0] 
adn.add_parameter('t_vr_eva_visu', [-9999 if i < i_vr_eva else 9999 for i in range(len(dt))], '', datatype="[Donnees EV]")
print ('\n'+'#'*50)*3
print 'instant debut rotation (detection EVA) => ',t_vr_eva

try:
   TIME = adn.get_parameter('DT', datatype="[Donnees EV]")
   TETAL = adn.get_parameter('TETAL', datatype="[Donnees EV]")
   DTETALSDT = derivative(array(TIME),array(TETAL))
   adn.add_parameter('DTETALSDT', DTETALSDT, '', datatype="[Donnees EV]")
   VC = adn.get_parameter('VCADC1L', datatype="[Donnees EV]")
   RCD = adn.get_parameter('RCD', datatype="[Donnees EV]")
   RCG = adn.get_parameter('RCG', datatype="[Donnees EV]")
   DDMPIL = adn.get_parameter('DDM', datatype="[Donnees EV]")
   #DDMPILf=fsmoo(DDMPIL, 0.1, fc=0.5, fmin=0.0, fmax=30.0, Ts=1.0)[0]
   #adn.add_parameter('DDMPILf', DDMPILf, '', datatype="[Donnees EV]")
   DDMCOP = adn.get_parameter('DDMCOP', datatype="[Donnees EV]")
   DDM = [ min(DDMCOP[i],DDMPIL[i]) for i in range(len(DDMPIL))]
   DDDMSDT = derivative(array(TIME),array(DDM))
   adn.add_parameter('DDDMSDT', DDDMSDT, '', datatype="[Donnees EV]")
   adn.add_parameter('DDMBGA', DDM, '', datatype="[Donnees EV]")
   
   tmp = [i for i in range(len(VC))  if DTETALSDT[i]>0.0 and VC[i]>80.0 and RCD[i]>100.0 and (DDM[i]<-3.0 or DDMCOP[i]<-3.0) and DDDMSDT[i]<-3.0 ]   
   pts_select = [ tmp[i] for i in range(len(tmp)) if (tmp[i]==tmp[min(i+1,len(tmp)-1)]-1 and not i==len(tmp)-1)  or (tmp[i]==tmp[max(i-1,0)]+1 and not i==0)]
   
   ddm_extract = [ DDM[i] for i in pts_select]
   adn.add_parameter('ddm_extract', [DDM[i]  if i in pts_select  else 0.0 for i in range(len(DDM))  ], '', datatype="[Donnees EV]")
   dt_extract = [ TIME[i] for i in pts_select]
   poly_coeff_1 = poly.polyfit(dt_extract, ddm_extract, 1)
   poly_1 = P(poly_coeff_1)
   tmp = poly.polyroots(poly_coeff_1)[0] #regression lineaire de l ordre manche
   
   
   i_vr_egv = [i for i in tarange(len(dt)) if (dt[i] > tmp  and DDM[i]<-0.5) or DDM[i]<-4.0 ][0] 
   t_vr_egv_visu = TIME[i_vr_egv]
   adn.add_parameter('t_vr_egv_visu', [-9999 if i < i_vr_egv else 9999 for i in range(len(TIME))], '', datatype="[Donnees EV]")
   print 'instant debut rotation (detection EGV) => ',t_vr_egv_visu
except:
   print 'pb critere debut rotation'
   pass

#Petite bite
#--------------------------------------------------------------
# analyse liftoff
#--------------------------------------------------------------
t_lo_eva = adn.get_parameter('Instant de Lift Off', datatype="[Caracteristiques de l'essai]")[0] 
adn.add_parameter('t_lo_eva', [t_lo_eva]*len(dt), '', datatype="[Donnees EV]")
print 'instant liftoff  (detection EVA) => ',t_lo_eva
i_lo_eva = [i for i in range(len(dt)) if abs(dt[i] - t_lo_eva)<0.001 ][0] 
adn.add_parameter('t_lo_eva_visu', [-9999 if i <= i_lo_eva else 9999 for i in range(len(dt))], '', datatype="[Donnees EV]")

try:
   TIME = adn.get_parameter('DT', datatype="[Donnees EV]")
   VC = adn.get_parameter('VCADC1L', datatype="[Donnees EV]")
   HCINEL = adn.get_parameter('HCINEL', datatype="[Donnees EV]")
   ZD = adn.get_parameter('ZD', datatype="[Donnees EV]")
   ZG = adn.get_parameter('ZG', datatype="[Donnees EV]")
   # valeur moyenne charge train en vol
   mean_load_ground = mean(array([ 0.5 * (ZD[i] + ZG[i]) for i in range(len(VC))  if (VC[i]>20.0 and VC[i]<50) ]))
   if mean_load_ground < 0.0:
      ZD = -ZD
      ZG = -ZG
      mean_load_ground = -mean_load_ground
   mean_load_flight = mean(array([ 0.5 * (ZD[i] + ZG[i]) for i in range(len(HCINEL))  if HCINEL[i]>10.0 and HCINEL[i]<50 ]))
   #print 'Charge moyenne sur les train au sol (20kts<VC<50kts)  ==> ', mean_load_ground
   #print 'Charge moyenne sur les train en vol (HCINEL>10ft)  ==> ', mean_load_flight
   crit_load = mean_load_flight + 0.05 * (mean_load_ground - mean_load_flight)
   i_lo_egv = [i for i in range(len(ZD))  if ( ZD[i]<crit_load and  ZG[i]<crit_load and i>i_vr_egv and HCINEL[i]<20) ][0]
   t_lo_egv_visu = TIME[i_lo_egv]
   adn.add_parameter('t_lo_egv_visu', [-9999 if i < i_lo_egv else 9999 for i in range(len(TIME))], '', datatype="[Donnees EV]")
   print 'instant liftoff  (detection EGV) => ',t_lo_egv_visu
except:
   print 'pb critere liftoff'
   pass




#--------------------------------------------------------------
# analyse passage 35ft
#--------------------------------------------------------------
t_35_eva = adn.get_parameter('Instant de passage des 35ft', datatype="[Caracteristiques de l'essai]")[0] 
adn.add_parameter('t_35_eva', [t_35_eva]*len(dt), '', datatype="[Donnees EV]")
print 'instant passage 35ft  (detection EVA) => ',t_35_eva
try:
   i_35_eva = [i for i in range(len(dt)) if abs(dt[i] - t_35_eva)<0.001 ][0] 
except:
   i_35_eva = len(dt)-1

adn.add_parameter('t_35_eva_visu', [-9999 if i <= i_35_eva else 9999 for i in range(len(dt))], '', datatype="[Donnees EV]")

try:
   TIME = adn.get_parameter('DT', datatype="[Donnees EV]")
   HCINEL = adn.get_parameter('HCINEL', datatype="[Donnees EV]")
   i_35_egv = [i for i in range(len(HCINEL))  if HCINEL[i]>35.0 ][0]
   t_35_egv_visu = TIME[i_35_egv]
   adn.add_parameter('t_35_egv_visu', [-9999 if i <= i_35_egv else 9999 for i in range(len(TIME))], '', datatype="[Donnees EV]")
   print 'instant passage 35ft  (detection EGV) => ',t_35_egv_visu
except:
   print 'pb critere passage 35ft'
   pass

print ('#'*50+'\n')*3



choix_trait = raw_input('utilisation temps caracteristiques EVA ou EGV ou manuel? (EVA / EGV / MAN) : ')
if choix_trait == 'EGV':
   i_vr = i_vr_egv
   i_lo = i_lo_egv
   i_35 = i_35_egv
if choix_trait == 'EVA':
   i_vr = i_vr_eva
   i_lo = i_lo_eva
   i_35 = i_35_eva
if choix_trait == 'MAN':
   t_vr = float(raw_input('Temps debut rotation : '))
   t_lo = float(raw_input('Temps liftoff : '))
   t_35 = float(raw_input('Temps passage 35ft : '))
   i_vr = [i for i in range(len(dt)) if dt[i] > t_vr ][0] 
   i_lo = [i for i in range(len(dt)) if dt[i] > t_lo ][0] 
   i_35 = [i for i in range(len(dt)) if dt[i] > t_35 ][0] 
t_vr = TIME[i_vr]
t_lo = TIME[i_lo]
t_35 = TIME[i_35]

adn.add_parameter('t_vr_visu', [-9999 if i < i_vr else 9999 for i in range(len(dt))], '', datatype="[Donnees EV]")
adn.add_parameter('t_lo_visu', [-9999 if i < i_lo else 9999 for i in range(len(dt))], '', datatype="[Donnees EV]")
adn.add_parameter('t_35_visu', [-9999 if i < i_35 else 9999 for i in range(len(dt))], '', datatype="[Donnees EV]")






vr_adn = adn.get_parameter('VCADC1L', datatype="[Donnees EV]")[i_vr]  
adn.add_parameter('vr_adn', [vr_adn]*len(dt), '', datatype="[Donnees EV]")
print ('#'*50+'\n')*3
print 'CAS debut rotation',vr_adn
teta_vr_adn = adn.get_parameter('TETAL', datatype="[Donnees EV]")[i_vr]  
adn.add_parameter('teta_vr_adn', [teta_vr_adn]*len(dt), '', datatype="[Donnees EV]")
print 'TETA debut rotation',teta_vr_adn
adn.add_parameter('dt_vr', array(dt)-t_vr, '', datatype="[Donnees EV]")

vlo_adn = adn.get_parameter('VCADC1L', datatype="[Donnees EV]")[i_lo]  
adn.add_parameter('vlo_adn', [vlo_adn]*len(dt), '', datatype="[Donnees EV]")
print 'CAS liftoff ',vlo_adn
teta_vlo_adn = adn.get_parameter('TETAL', datatype="[Donnees EV]")[i_lo]  
adn.add_parameter('teta_vlo_adn', [teta_vlo_adn]*len(dt), '', datatype="[Donnees EV]")
print 'TETA liftoff',teta_vlo_adn

v35_adn = adn.get_parameter('VCADC1L', datatype="[Donnees EV]")[i_35]  
adn.add_parameter('v35_adn', [v35_adn]*len(dt), '', datatype="[Donnees EV]")
print 'CAS passage 35ft ',v35_adn
teta_v35_adn = adn.get_parameter('TETAL', datatype="[Donnees EV]")[i_35]  
adn.add_parameter('teta_v35_adn', [teta_v35_adn]*len(dt), '', datatype="[Donnees EV]")
print 'TETA passage 35ft',teta_v35_adn
TETAL = adn.get_parameter('TETAL', datatype="[Donnees EV]")
adn.add_parameter('dTETAL', array(TETAL) - teta_vr_adn, '', datatype="[Donnees EV]")

adn.add_parameter('dq_9deg', [-9.0]*len(dt), '', datatype="[Donnees EV]")
adn.add_parameter('dq_12deg', [-12.0]*len(dt), '', datatype="[Donnees EV]")
adn.add_parameter('dq_13.5deg', [-13.5]*len(dt), '', datatype="[Donnees EV]")


# detection du moment ou teta = teta target
adn.add_parameter('teta_target', [teta_target]*len(dt), '', datatype="[Donnees EV]")
i_teta_target = [ i for i in range(len(TETAL)) if TETAL[i]-teta_target>0.0 ]
if len(i_teta_target) > 0:
   i_teta_target = i_teta_target[0]
else:
   i_teta_target = 0
   print 'teta target non atteint ....'
adn.add_parameter('t_teta_target_visu', [-9999 if i <= i_teta_target else 9999 for i in range(len(dt))], '', datatype="[Donnees EV]")
adn.add_parameter('teta_target_m_tetavr', [teta_target - teta_vr_adn]*len(dt), '', datatype="[Donnees EV]")



#beta target
beta_target=dict()
beta_target['1814']=([120,140,180,220,260,280],[2.2,1.85,1.45,1.1,0.95,0.95])
beta_target['1820']=([120,140,180,220,260,280],[2.2,1.85,1.45,1.1,0.95,0.95])
beta_target['1826']=([120,140,180,220,260,280],[2.2,1.85,1.45,1.1,0.95,0.95])
Configuration = str(adn.get_parameter('Configuration', datatype="[Caracteristiques de l'essai]")[0])[0:4]
print 'Configuration',Configuration
VCADC1L = adn.get_parameter('VCADC1L', datatype="[Donnees EV]")
beta_target_ft = [ interp(VCADC1L[i],beta_target[Configuration][0],beta_target[Configuration][1]) for i in range(len(VCADC1L))]
adn.add_parameter('beta_target_ft', beta_target_ft, '', datatype="[Donnees EV]")
adn.add_parameter('beta_target_ft_neg', -array(beta_target_ft), '', datatype="[Donnees EV]")


# detection panne moteur
n11 = adn.get_parameter('N11L', datatype="[Donnees EV]")
n12 = adn.get_parameter('N12L', datatype="[Donnees EV]")
h = adn.get_parameter('HCINEL', datatype="[Donnees EV]")
dn = mean(abs(array([ n11[i] - n12[i] for i in range(len(n11)) if (h[i]>5 and  h[i]<35.0) ])))
adn.save(path=fic_adn+'.csv', fmt='csv') 
print ('#'*50+'\n')*3
if dn > 10.0:
   os.system('/opt/soft/bin/visu_mdvs --version 14.8p5  -ef template  -fs KPTOATT_OEI.dat -fs  '+fic_adn+'.csv -pfc '+fic_adn+'.pdf  -a')
   os.system('/opt/soft/bin/visu_mdvs --version 14.8p5  -ef template  -fs KPTOATT_OEI.dat -fs  '+fic_adn+'.csv')
if dn <= 5.0:
   os.system('/opt/soft/bin/visu_mdvs --version 14.8p5  -ef template  -fs KPTOATT_AEO.dat -fs  '+fic_adn+'.csv -pfc '+fic_adn+'.pdf  -a')
   os.system('/opt/soft/bin/visu_mdvs --version 14.8p5  -ef template  -fs KPTOATT_AEO.dat -fs  '+fic_adn+'.csv')

  
