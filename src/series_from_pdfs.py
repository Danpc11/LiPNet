"""Group- and patient-level data transcribed from the 13 uploaded PDFs (no values from memory).
z_F = PVF(mL/min/100 g graft) / donor reference (study's own donor PVF/100 g when reported, else 90 = healthy live donor,
      Sainz-Barriga 2010 / Troisi 2003 Table 2 / Chan 2011 Table 2 give 90, 91, 81).
z_P = (PVP - CVP)/5 ; CVP measured when reported, otherwise assumed 5 mmHg and flagged."""
import numpy as np, pandas as pd
R = []  # study, group, n, source, GRWR, PVF100, donor_ref, PVP, CVP, cvp_meas, outcome, events, z_F, z_P, notes
def add(study, group, n, src, grwr, pvf, ref, pvp, cvp, cvpm, outcome, ev, notes=''):
    zF = pvf/ref if (pvf is not None and ref) else np.nan
    zP = (pvp-cvp)/5 if pvp is not None else np.nan
    R.append(dict(estudio=study, grupo=group, n=n, fuente=src, GRWR=grwr, PVF_100g=pvf, ref_donante=ref, PVP=pvp, CVP=cvp,
                  CVP_medida=cvpm, desenlace=outcome, eventos=ev, pct=(100*ev/n if (n and ev is not None) else np.nan), z_F=zF, z_P=zP, notas=notes))
# --- Troisi 2003, Liver Transpl 9:S36 (Tables 2,3,4,5)
add('Troisi 2003','G1 sin GIM',11,'Tab 3 & 5, p.S38-39',1.12,268,91,None,None,False,'SFSS',3,'PVF/GW medio 268±230; donante 91±16 (Tab 2)')
add('Troisi 2003','G2 tras ligadura esplénica',13,'Tab 4 & 5',1.13,240,91,None,None,False,'SFSS',0,'antes GIM 360±143 -> después 240±91')
# --- Troisi 2005, AJT 5:1397 (Tables 3,4; Results)
add('Troisi 2005','G1 sin GIM (GRWR<=0.8)',5,'Tab 4, p.1399',0.73,401,118,None,None,False,'SFSS',3,'donante 117-120 (Tab 3); superv injerto 1a 20%')
add('Troisi 2005','G2 shunt hemiportocava',8,'Tab 4, p.1399',0.71,190,118,None,None,False,'SFSS',0,'PVF1 (shunt clampeado) 537 -> PVF2 190; superv injerto 75%')
# --- Ou 2010, Transplant Proc 42:876 (Table 2, individual)
add('Ou 2010','rPVF>250 sin modulación',2,'Tab 2 & 3, p.878',0.93,268,90,None,None,False,'SFSS',2,'pacientes 1-2: 251, 285 mL/min/100 g')
add('Ou 2010','rPVF>250 con SAL/esplenectomía',6,'Tab 2 & 3, p.878',1.06,186,90,None,None,False,'SFSS',1,'pre 296,416,554,257,413,495 -> post 186,147,236,87,165,296 (media 186)')
# --- Vasavada 2014, Int J Surg 12:177 (Table 1; sect 4.7)
add('Vasavada 2014','PVF>190 (día 0)',27,'sect 4.7, p.179',None,None,90,None,None,False,'disfunción temprana (def. Kyushu)',10,'z_F >=2.1 (umbral); medias: disfunción 193±54 vs no 139±55 (Tab 1)')
R[-1]['z_F']=190/90
add('Vasavada 2014','PVF<=190',107,'sect 4.7, p.179',None,None,90,None,None,False,'disfunción temprana',9,'z_F <2.1; media de grupo no reportada')
R[-1]['z_F']=139/90
# --- Alim 2016, Liver Transpl 22:1643 (Results; Tab 1,2)
add('Alim 2016','GRWR<0.8 (SAL si PVF>250)',43,'Results p.1645; Tab 2',0.76,220,90,None,None,False,'SFSS',1,'PVF mediana 220 (195-250); 5 SAL 336->240; los 6 con desenlace adverso tenían PVF 210-250 (Tab 1)')
# --- Chan 2011, Liver Transpl 17:115 (Table 2,3)
add('Chan 2011','lóbulo derecho con VHM, sin modulación',46,'Tab 2, p.116; Tab 3',None,318,81,14,8,True,'muerte hospitalaria',2,'G/SLV 47%; gradiente PVP-CVP medido 6 (-1..11) -> CVP implícita 8; complicaciones no asociadas a flujo/presión')
R[-1]['z_P']=6/5
# --- Yagi 2006, Transplantation 81:373 (Results; Tab 1)
add('Yagi 2006','PVP modulada <20 (n=27)',27,'Results p.374; Tab 1',1.06,295,90,14.2,5.2,True,'muerte 1 año',2,'PVF 1828±592 mL/min, peso injerto mediana 620 g -> ~295/100 g (cálculo propio); superv 1a 92.3%')
# --- Wang 2014, Surg Today (Table 1; Fig 3c; Tab 2)
add('Wang 2014','esplenectomía',154,'Tab 1; Fig 3c',None,371,90,16.1,5,False,'disfunción primaria del injerto',15,'9.7%; PVP 24.9 -> 16.1 al cierre')
add('Wang 2014','sin esplenectomía',122,'Tab 1; Fig 3c',None,337,90,18.4,5,False,'disfunción primaria del injerto',24,'19.7%')
add('Wang 2014','PVP cierre >=20',54,'Tab 2',None,None,None,21,5,False,'pérdida del injerto <6 m',9,'9/54 vs 21/238 (Tab 2)')
add('Wang 2014','PVP cierre <20',238,'Tab 2',None,None,None,16,5,False,'pérdida del injerto <6 m',21,'PVP medio asumido 16')
# --- Osman 2017, Hepatol Res 47:293 (Table 1)
add('Osman 2017','A: PVP final <15',39,'Tab 1, p.298',1.06,None,None,12.15,5,False,'SFSS',1,'mortalidad 90 d 3/39')
add('Osman 2017','B: PVP final 15-19',37,'Tab 1, p.298',1.00,None,None,17.14,5,False,'SFSS',6,'mortalidad 90 d 9/37')
# --- Ogura 2010, Liver Transpl 16:718 (Results: portocaval gradient)
add('Ogura 2010','PVP final <15',86,'Results p.723',None,None,None,12.2,6.0,True,'muerte 2 años',6,'gradiente medido 6.2±0.3; superv 2a 93.0% -> ~6 muertes')
R[-1]['z_P']=6.2/5
add('Ogura 2010','PVP final >=15',43,'Results p.723',None,None,None,16.7,6.8,True,'muerte 2 años',14,'gradiente medido 9.9±0.4; superv 2a 66.3% -> ~14 muertes')
R[-1]['z_P']=9.9/5
# --- Yagi 2005, Liver Transpl 11:68 (Results)
add('Yagi 2005','L: PVP <20 (POD 1-3)',15,'Results p.70',None,None,None,15.4,5,False,'muerte 1 año',1,'PVP medio 15.4 (13-17.6); superv 92.9%')
add('Yagi 2005','H: PVP >=20 (POD 1-3)',17,'Results p.70',None,None,None,23.7,5,False,'muerte 1 año',7,'PVP medio 23.7 (20-31); superv 58.8%')
# --- Yamada 2008, AJT 8:847 (Table 4 individual PVP/CVP before closure)
F=[(23,17),(19,9),(12,9),(16,9),(20,10),(20,7),(19,8),(16,10),(20,10),(20,12)]
g=np.mean([p-c for p,c in F])
add('Yamada 2008','HPCS por PVP>20 (GRWR 0.58-0.79)',10,'Tab 4 col F, p.851',0.70,None,None,np.mean([p for p,_ in F]),np.mean([c for _,c in F]),True,'SFSS',0,f'gradiente individual medio {g:.1f} mmHg (rango {min(p-c for p,c in F)}-{max(p-c for p,c in F)})')
R[-1]['z_P']=g/5
df=pd.DataFrame(R); df.to_csv('results/series_pdfs_z.tsv',sep='\t',index=False,float_format='%.2f')
print(df[['estudio','grupo','n','PVF_100g','ref_donante','PVP','CVP','CVP_medida','z_F','z_P','desenlace','eventos','pct']].round(2).to_string(index=False))
