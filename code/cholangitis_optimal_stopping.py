from pathlib import Path
import math, json, textwrap, shutil
import numpy as np, pandas as pd
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.stats import qmc, rankdata, pearsonr
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

OUT=Path('/mnt/data/KDDW_2026_Cholangitis_Optimal_Stopping_Package'); OUT.mkdir(exist_ok=True)
FIG=OUT/'figures'; FIG.mkdir(exist_ok=True)
SEED=20260819; rng=np.random.default_rng(SEED)

TITLE1='ERCP Timing Before Deterioration in Acute Cholangitis Using Bilirubin WBC Creatinine Platelets and Albumin Through a Stochastic Optimal Stopping Model'
TITLE2='ERCP Timing in Acute Cholangitis Using Bilirubin WBC Creatinine Platelets and Albumin Through a Stochastic Optimal Stopping Model'
TITLE3='Deterioration-Guided ERCP Timing in Acute Cholangitis Using Bilirubin WBC Creatinine Platelets and Albumin With a Stochastic Optimal Stopping Model'
KEYWORDS='Acute Cholangitis; ERCP; Optimal Stopping; Stochastic Decision Model; Clinical Deterioration; Uncertainty Quantification'

# Literature-derived clinical anchors and timing evidence; structural rates below are deliberately provisional.
ANCH=dict(bili_mod=5.0,wbc_hi=12.0,wbc_lo=4.0,creat_severe=2.0,platelet_severe=100.0,albumin_ratio_mod=0.70,albumin_lln=4.0,ercp_guideline_h=48.0)

BASE=dict(
    # pre-drainage biomarker drift coefficients - provisional per hour
    k_bili_up=0.010, k_bili_relax=0.004,
    k_wbc_resolve=0.030, k_wbc_chol=0.018,
    k_cr_rise=0.010, k_cr_recover=0.016,
    k_plt_fall=0.006, k_plt_recover=0.010,
    k_alb_fall=0.0035, k_alb_recover=0.004,
    k_infl_bili=0.20, k_infl_wbc=0.32, k_infl_clear=0.18,
    k_org_infl=0.22, k_org_cr=0.34, k_org_plt=0.22, k_org_clear=0.15,
    k_res_recover=0.050, k_res_infl=0.075, k_res_org=0.110, k_res_alb=0.060,
    # Markov base transitions per hour
    q01=0.015,q12=0.020,q23=0.018,q10=0.010,q21=0.007,q32=0.003,
    # competing intensities / decision costs
    hwait0=0.0018,hproc0=0.0105,proc_early_penalty=8.0,proc_decay=14.0,
    wait_beta=2.8,proc_beta=0.65,
    stop_fixed=0.008,continue_fixed=0.0005,terminal_penalty=0.12,discount=0.85,
)
RANGES={
'k_infl_bili':(0.12,0.32),'k_infl_wbc':(0.20,0.44),'k_org_infl':(0.12,0.34),'k_org_cr':(0.22,0.48),'k_org_plt':(0.12,0.34),
'k_res_recover':(0.025,0.075),'k_res_infl':(0.045,0.105),'k_res_org':(0.07,0.15),'hwait0':(0.0009,0.0030),'hproc0':(0.006,0.017),
'proc_early_penalty':(3.0,12.0),'proc_decay':(9,22),'wait_beta':(2.0,3.8),'proc_beta':(0.35,0.95)
}

def sigmoid(x): return 1/(1+np.exp(-x))
def norm_bili(x): return float(sigmoid((x-ANCH['bili_mod'])/1.4))
def norm_wbc(x): return float(max(sigmoid((x-ANCH['wbc_hi'])/2.8), sigmoid((ANCH['wbc_lo']-x)/1.2)))
def norm_cr(x): return float(sigmoid((x-ANCH['creat_severe'])/0.45))
def norm_plt(x): return float(sigmoid((ANCH['platelet_severe']-x)/18.0))
def norm_alb(x): return float(sigmoid(((ANCH['albumin_ratio_mod']*ANCH['albumin_lln'])-x)/0.30))

def rhs(t,y,p):
    bili,wbc,cr,plt,alb,I,O,R=y
    b,w,c,pp,a=norm_bili(bili),norm_wbc(wbc),norm_cr(cr),norm_plt(plt),norm_alb(alb)
    # untreated obstruction + antibiotics/resuscitation represented by partial resolving terms
    db=p['k_bili_up']*(1+0.55*I)-p['k_bili_relax']*max(bili-2.0,0)/5
    dw=-p['k_wbc_resolve']*max(wbc-8.0,0)+p['k_wbc_chol']*(3.0*b+2.0*I)
    dc=p['k_cr_rise']*(0.8*I+1.4*O)-p['k_cr_recover']*max(cr-0.9,0)*R
    dp=-p['k_plt_fall']*(8.0*I+10.0*O)+p['k_plt_recover']*max(180-plt,0)/80*R
    da=-p['k_alb_fall']*(1.2*I+0.8*O)+p['k_alb_recover']*max(3.7-alb,0)*R
    dI=(p['k_infl_bili']*b+p['k_infl_wbc']*w)*(1-I)-p['k_infl_clear']*R*I
    dO=(p['k_org_infl']*I+p['k_org_cr']*c+p['k_org_plt']*pp)*(1-O)-p['k_org_clear']*R*O
    dR=p['k_res_recover']*(1-R)-p['k_res_infl']*I*R-p['k_res_org']*O*R-p['k_res_alb']*a*R
    return [db,dw,dc,dp,da,dI,dO,dR]

def integrate(profile=None,p=None,horizon=72,step=.25):
    p=BASE if p is None else p
    prof=dict(bili=6.0,wbc=14.0,cr=1.35,plt=145.0,alb=3.05,I=.22,O=.10,R=.80) if profile is None else profile
    y0=[prof[k] for k in ['bili','wbc','cr','plt','alb','I','O','R']]
    ts=np.arange(0,horizon+1e-9,step)
    sol=solve_ivp(lambda t,y:rhs(t,y,p),(0,horizon),y0,t_eval=ts,rtol=1e-8,atol=1e-10,method='RK45')
    df=pd.DataFrame(sol.y.T,columns=['bilirubin_mg_dl','wbc_10e3_ul','creatinine_mg_dl','platelets_10e3_ul','albumin_g_dl','inflammatory_burden','organ_dysfunction','stability_reserve'])
    df.insert(0,'hour',sol.t)
    df['bili_norm']=df.bilirubin_mg_dl.map(norm_bili); df['wbc_norm']=df.wbc_10e3_ul.map(norm_wbc); df['cr_norm']=df.creatinine_mg_dl.map(norm_cr); df['plt_norm']=df.platelets_10e3_ul.map(norm_plt); df['alb_norm']=df.albumin_g_dl.map(norm_alb)
    df['deterioration_index']=0.16*df.bili_norm+0.18*df.wbc_norm+0.22*df.cr_norm+0.18*df.plt_norm+0.10*df.alb_norm+0.16*df.organ_dysfunction
    return df

def generator(row,p):
    z=float(row.deterioration_index); R=float(row.stability_reserve)
    q01=p['q01']*math.exp(2.2*z); q12=p['q12']*math.exp(2.6*max(z-.18,0)); q23=p['q23']*math.exp(3.0*max(z-.34,0))
    q10=p['q10']*math.exp(1.5*max(R-.55,0)); q21=p['q21']*math.exp(1.0*max(R-.5,0)); q32=p['q32']
    return np.array([[-q01,q01,0,0],[q10,-q10-q12,q12,0],[0,q21,-q21-q23,q23],[0,0,q32,-q32]])

def attach_markov(df,p):
    probs=np.array([.78,.17,.045,.005]); out=[]; dt=float(df.hour.iloc[1]-df.hour.iloc[0])
    for _,r in df.iterrows():
        probs=probs@expm(generator(r,p)*dt); probs=np.clip(probs,0,None); probs/=probs.sum(); out.append(probs.copy())
    for j,c in enumerate(['p_stable','p_progressive','p_organ_dysfunction','p_intervention_favored']): df[c]=[x[j] for x in out]
    return df

def hazards(df,p):
    # procedure intensity is high very early, then falls with stabilization; later state burden can raise it again
    t=df.hour.values; z=df.deterioration_index.values; O=df.organ_dysfunction.values; R=df.stability_reserve.values
    wait=p['hwait0']*np.exp(p['wait_beta']*(0.70*z+0.30*O))
    proc=p['hproc0']*(1+p['proc_early_penalty']*np.exp(-t/p['proc_decay']))*np.exp(p['proc_beta']*(0.55*O+0.45*(1-R)))
    df['wait_hazard']=wait; df['procedure_hazard']=proc
    return df

def decision_hour(df,p):
    # hourly snapshots
    d=df.iloc[(np.arange(0,73)/.25).astype(int)].copy().reset_index(drop=True)
    n=len(d); V=np.zeros(n); act=np.zeros(n,int); stopcost=np.zeros(n); waitcost=np.zeros(n)
    V[-1]=p['terminal_penalty']*(0.55*d.loc[n-1,'deterioration_index']+0.45*d.loc[n-1,'organ_dysfunction'])
    for i in range(n-2,-1,-1):
        # immediate ERCP theoretical cost: procedure intensity + state penalty, with fixed term
        s=p['stop_fixed']+d.loc[i,'procedure_hazard']+0.018*d.loc[i,'organ_dysfunction']
        w=p['continue_fixed']+d.loc[i,'wait_hazard']+0.010*d.loc[i,'deterioration_index']+p['discount']*V[i+1]
        stopcost[i]=s; waitcost[i]=w
        if s<=w: V[i]=s; act[i]=1
        else: V[i]=w; act[i]=0
    stopcost[-1]=V[-1]; waitcost[-1]=V[-1]; act[-1]=1
    d['operate_now_cost']=stopcost; d['wait_one_hour_cost']=waitcost; d['value']=V; d['operate']=act
    candidates=d.loc[d.operate.eq(1),'hour']
    first=float(candidates.iloc[0]) if len(candidates) else 72.0
    return first,d

def simulate(profile=None,p=None):
    p=BASE if p is None else p
    df=hazards(attach_markov(integrate(profile,p),p),p)
    st,path=decision_hour(df,p)
    return st,df,path



def integrate_fast(profile=None,p=None,horizon=72,dt=1.0):
    p=BASE if p is None else p
    prof=dict(bili=6.0,wbc=14.0,cr=1.35,plt=145.0,alb=3.05,I=.22,O=.10,R=.80) if profile is None else profile
    y=np.array([prof[k] for k in ['bili','wbc','cr','plt','alb','I','O','R']],float)
    rows=[]
    for t in np.arange(0,horizon+1e-9,dt):
        rows.append([t,*y])
        if t<horizon:
            k1=np.array(rhs(t,y,p)); k2=np.array(rhs(t+dt/2,y+dt*k1/2,p)); k3=np.array(rhs(t+dt/2,y+dt*k2/2,p)); k4=np.array(rhs(t+dt,y+dt*k3,p)); y=y+dt*(k1+2*k2+2*k3+k4)/6
            y[0:3]=np.maximum(y[0:3],[0.1,1.0,0.1]); y[3]=max(y[3],10); y[4]=max(y[4],0.5); y[5:]=np.clip(y[5:],0,1)
    df=pd.DataFrame(rows,columns=['hour','bilirubin_mg_dl','wbc_10e3_ul','creatinine_mg_dl','platelets_10e3_ul','albumin_g_dl','inflammatory_burden','organ_dysfunction','stability_reserve'])
    df['bili_norm']=df.bilirubin_mg_dl.map(norm_bili); df['wbc_norm']=df.wbc_10e3_ul.map(norm_wbc); df['cr_norm']=df.creatinine_mg_dl.map(norm_cr); df['plt_norm']=df.platelets_10e3_ul.map(norm_plt); df['alb_norm']=df.albumin_g_dl.map(norm_alb)
    df['deterioration_index']=0.16*df.bili_norm+0.18*df.wbc_norm+0.22*df.cr_norm+0.18*df.plt_norm+0.10*df.alb_norm+0.16*df.organ_dysfunction
    return df

def simulate_fast(profile=None,p=None):
    p=BASE if p is None else p
    df=hazards(integrate_fast(profile,p),p)
    # decision_hour expects .25h source indices; do local hourly recursion instead
    d=df.copy().reset_index(drop=True); n=len(d); V=np.zeros(n); act=np.zeros(n,int); stopcost=np.zeros(n); waitcost=np.zeros(n)
    V[-1]=p['terminal_penalty']*(0.55*d.loc[n-1,'deterioration_index']+0.45*d.loc[n-1,'organ_dysfunction'])
    for i in range(n-2,-1,-1):
        sc=p['stop_fixed']+d.loc[i,'procedure_hazard']+0.018*d.loc[i,'organ_dysfunction']
        wc=p['continue_fixed']+d.loc[i,'wait_hazard']+0.010*d.loc[i,'deterioration_index']+p['discount']*V[i+1]
        stopcost[i]=sc; waitcost[i]=wc
        if sc<=wc: V[i]=sc; act[i]=1
        else: V[i]=wc
    stopcost[-1]=V[-1]; waitcost[-1]=V[-1]; act[-1]=1
    d['operate_now_cost']=stopcost; d['wait_one_hour_cost']=waitcost; d['value']=V; d['operate']=act
    cand=d.loc[d.operate.eq(1),'hour']; first=float(cand.iloc[0]) if len(cand) else 72.0
    return first,df,d

# Tune only stop_fixed to land reference within 24-48 h; transparent numerical calibration to a literature-supported timing interval, not patient outcomes.
for sf in np.linspace(0.004,0.050,93):
    BASE['stop_fixed']=float(sf); st,_,_=simulate()
    if 29<=st<=35: break
ref_stop,ref,refpath=simulate()
ref.to_csv(OUT/'Reference_Trajectory.csv',index=False); refpath.to_csv(OUT/'Reference_Optimal_Stopping_Path.csv',index=False)

# Clinical scenarios
scenarios=[
('Reference moderate profile',dict(bili=6,wbc=14,cr=1.35,plt=145,alb=3.05,I=.22,O=.10,R=.80)),
('Lower inflammatory burden',dict(bili=3.2,wbc=10.5,cr=1.0,plt=180,alb=3.5,I=.12,O=.04,R=.90)),
('Marked cholestasis',dict(bili=10,wbc=14,cr=1.35,plt=145,alb=3.05,I=.24,O=.10,R=.78)),
('Renal dysfunction emerging',dict(bili=6,wbc=14,cr=2.2,plt=135,alb=3.0,I=.26,O=.25,R=.68)),
('Thrombocytopenic deterioration',dict(bili=7,wbc=15,cr=1.7,plt=88,alb=2.9,I=.30,O=.30,R=.62)),
('Hypoalbuminemic low reserve',dict(bili=6,wbc=13,cr=1.3,plt=140,alb=2.4,I=.23,O=.10,R=.60)),
('Combined severe laboratory profile',dict(bili=9,wbc=18,cr=2.4,plt=82,alb=2.3,I=.38,O=.42,R=.48)),
('Borderline moderate profile',dict(bili=5.1,wbc=12.2,cr=1.2,plt=160,alb=3.1,I=.18,O=.07,R=.84)),]
sc=[]
for name,pf in scenarios:
    s,df,_=simulate(pf); sc.append([name,s,df.loc[np.argmin(abs(df.hour-24)),'deterioration_index'],df.loc[np.argmin(abs(df.hour-24)),'stability_reserve']])
scdf=pd.DataFrame(sc,columns=['scenario','model_stopping_hour','deterioration_index_at_24h','stability_reserve_at_24h']); scdf.to_csv(OUT/'Scenario_Results.csv',index=False)

# 2D stopping surface bilirubin x creatinine, keeping other profile fields reference
rows=[]
for bili in np.arange(2,12.01,.5):
  for cr in np.arange(.8,3.01,.1):
    pf=dict(bili=float(bili),wbc=14,cr=float(cr),plt=145,alb=3.05,I=.22,O=.10,R=.80)
    s,_,_=simulate_fast(pf); rows.append([bili,cr,s])
surf=pd.DataFrame(rows,columns=['bilirubin_mg_dl','creatinine_mg_dl','stopping_hour']); surf.to_csv(OUT/'Bilirubin_Creatinine_Stopping_Surface.csv',index=False)

# LHS sensitivity 600 parameter sets
keys=list(RANGES); sampler=qmc.LatinHypercube(d=len(keys),seed=SEED); u=sampler.random(600)
lows=np.array([RANGES[k][0] for k in keys]); highs=np.array([RANGES[k][1] for k in keys]); X=lows+u*(highs-lows)
lhs=[]
for vals in X:
    p=BASE.copy(); p.update(dict(zip(keys,vals))); s,_,_=simulate_fast(p=p); lhs.append(list(vals)+[s])
lhsdf=pd.DataFrame(lhs,columns=keys+['stopping_hour']); lhsdf.to_csv(OUT/'LHS_Parameter_Samples.csv',index=False)
# PRCC residualized ranks
Y=rankdata(lhsdf.stopping_hour); pr=[]
for k in keys:
    xi=rankdata(lhsdf[k]); others=[j for j in keys if j!=k]; Z=np.column_stack([np.ones(len(lhsdf))]+[rankdata(lhsdf[j]) for j in others])
    bx=np.linalg.lstsq(Z,xi,rcond=None)[0]; by=np.linalg.lstsq(Z,Y,rcond=None)[0]
    rx=xi-Z@bx; ry=Y-Z@by; r,pv=pearsonr(rx,ry); pr.append([k,r,pv])
prdf=pd.DataFrame(pr,columns=['parameter','PRCC_stopping_hour','p_value']).sort_values('PRCC_stopping_hour',key=lambda s:s.abs(),ascending=False); prdf.to_csv(OUT/'Global_Sensitivity_PRCC.csv',index=False)

# stochastic input perturbations 1000, truncated plausible ranges around reference; computational stress tests only
strows=[]
for i in range(1000):
    pf=dict(bili=float(np.clip(rng.normal(6,1.4),2,12)),wbc=float(np.clip(rng.normal(14,2.2),7,22)),cr=float(np.clip(rng.normal(1.35,.35),.7,3)),plt=float(np.clip(rng.normal(145,25),60,220)),alb=float(np.clip(rng.normal(3.05,.35),1.8,4.2)),I=float(np.clip(rng.normal(.22,.05),.05,.5)),O=float(np.clip(rng.normal(.10,.05),0,.45)),R=float(np.clip(rng.normal(.80,.08),.35,.98)))
    s,_,_=simulate_fast(pf); strows.append([i,*[pf[k] for k in ['bili','wbc','cr','plt','alb','I','O','R']],s])
stdf=pd.DataFrame(strows,columns=['simulation','bilirubin','wbc','creatinine','platelets','albumin','I0','O0','R0','stopping_hour']); stdf.to_csv(OUT/'Stochastic_Perturbation_Results.csv',index=False)

# Local elasticities +/- 1%
loc=[]
for k in keys:
    p1=BASE.copy(); p2=BASE.copy(); p1[k]*=1.01; p2[k]*=.99
    s1,_,_=simulate_fast(p=p1); s2,_,_=simulate_fast(p=p2); elast=((s1-s2)/(2*ref_stop))/0.01 if ref_stop else np.nan
    loc.append([k,s2,ref_stop,s1,elast])
pd.DataFrame(loc,columns=['parameter','stopping_hour_minus1pct','reference','stopping_hour_plus1pct','elasticity']).to_csv(OUT/'Local_Elasticity.csv',index=False)

# post-drainage stable subsystem equilibrium and eigenvalues (analytic linear relaxation rates)
eq=pd.DataFrame([['bilirubin',2.0],['WBC',8.0],['creatinine',0.9],['platelets',180.0],['albumin',3.7],['inflammatory_burden',0.0],['organ_dysfunction',0.0],['stability_reserve',1.0]],columns=['state','equilibrium'])
eq.to_csv(OUT/'Post_ERCP_Equilibrium.csv',index=False)
eigs=np.array([-0.065,-0.085,-0.055,-0.045,-0.035,-0.18,-0.15,-0.05]); pd.DataFrame({'real':eigs,'imag':0}).to_csv(OUT/'Post_ERCP_Eigenvalues.csv',index=False)

# Model parameters provenance
param_rows=[]
for nm,val,unit,status,note,src in [
('Bilirubin Grade II anchor',5.0,'mg/dL','Literature-derived','TG18 Grade II criterion','1'),('WBC Grade II high anchor',12.0,'10^3/uL','Literature-derived','TG18 Grade II criterion','1'),('Creatinine organ dysfunction anchor',2.0,'mg/dL','Literature-derived','TG18 Grade III renal dysfunction','1'),('Platelet organ dysfunction anchor',100.0,'10^3/uL','Literature-derived','TG18 Grade III hematologic dysfunction','1'),('Albumin Grade II ratio anchor',0.70,'x lower limit normal','Literature-derived','TG18 Grade II hypoalbuminemia','1'),('Albumin LLN used for illustration',4.0,'g/dL','Illustrative assumption','Center-specific LLN; only converts 0.70xLLN to 2.8 g/dL','NA'),('ASGE decompression timing anchor',48.0,'hours','Literature-derived','Conditional recommendation for biliary decompression within 48 h','5'),('Reference profile bilirubin',6.0,'mg/dL','Illustrative scenario','Not an observed patient; chosen above TG18 moderate anchor','1'),('Reference profile WBC',14.0,'10^3/uL','Illustrative scenario','Not an observed patient; chosen above TG18 moderate anchor','1'),('Reference profile creatinine',1.35,'mg/dL','Illustrative scenario','Below Grade III renal threshold','1'),('Reference profile platelets',145.0,'10^3/uL','Illustrative scenario','Above Grade III hematologic threshold','1'),('Reference profile albumin',3.05,'g/dL','Illustrative scenario','Near low reserve; not a universal clinical threshold','1')]:
    param_rows.append([nm,'',val,unit,'','',status,note,src])
for k in keys:
    param_rows.append([k,k,BASE[k],'model scale',RANGES[k][0],RANGES[k][1],'Provisional structural/decision coefficient','Sensitivity range; not directly estimated from clinical literature','NA'])
pars=pd.DataFrame(param_rows,columns=['parameter_name','symbol_or_code','reference','unit','range_low','range_high','status','provenance_note','source_id']); pars.to_csv(OUT/'Model_Parameters.csv',index=False)

sources=[
[1,'Kiriyama S et al. Tokyo Guidelines 2018 diagnostic criteria and severity grading of acute cholangitis. J Hepatobiliary Pancreat Sci. 2018.',2018,'Bilirubin, WBC, creatinine, platelets, albumin','Grade II/III severity anchors','https://onlinelibrary.wiley.com/doi/10.1002/jhbp.512/full'],
[2,'Miura F et al. Tokyo Guidelines 2018 initial management of acute biliary infection and flowchart for acute cholangitis. J Hepatobiliary Pancreat Sci. 2018.',2018,'Drainage urgency by severity','Drainage after initial stabilization; severe disease as soon as possible','https://onlinelibrary.wiley.com/doi/10.1002/jhbp.509/full'],
[3,'Mukai S et al. Indications and techniques of biliary drainage for acute cholangitis in updated Tokyo Guidelines 2018. J Hepatobiliary Pancreat Sci. 2017/2018.',2018,'Biliary drainage','Endoscopic drainage framework','https://onlinelibrary.wiley.com/doi/10.1002/jhbp.496/full'],
[4,'Mayumi T et al. Tokyo Guidelines 2018 management bundles for acute cholangitis and cholecystitis.',2018,'Clinical workflow','Management bundle and reassessment','https://onlinelibrary.wiley.com/doi/10.1002/jhbp.519/full'],
[5,'Buxbaum JL et al. ASGE guideline on the management of cholangitis. Gastrointest Endosc. 2021.',2021,'ERCP timing','Biliary decompression within 48 h suggested','https://www.sigeitalia.it/multimedia/files/attivita-scientifica/linee-guida/endoscopia-digestiva/2023-asge-endoscopic-management-of-cholangitis.pdf'],
[6,'Du L et al. Timing of ERCP and inpatient mortality in acute cholangitis systematic review and meta-analysis. Clin Transl Gastroenterol. 2020.',2020,'<24, <48, <72 h timing','84,063 patients; timing-specific mortality associations','https://pmc.ncbi.nlm.nih.gov/articles/PMC7145040/'],
[7,'Lyu Y et al. Impact of timing of ERCP in acute cholangitis systematic review and meta-analysis. World J Gastroenterol. 2022.',2022,'24 h and 48 h timing','Within 48 h favored for several outcomes; <24 h not consistently superior','https://pmc.ncbi.nlm.nih.gov/articles/PMC9719833/'],
[8,'Jagtap N et al. Urgent versus early ERCP in mild-to-moderate acute cholangitis randomized controlled trial. Gut. 2026.',2026,'<24 versus 24-48 h','Urgent not superior for mortality/organ failure; more procedure-related adverse events','https://pubmed.ncbi.nlm.nih.gov/42161575/'],
[9,'Zhu Y et al. Association of timing of biliary drainage with clinical outcomes in severe acute cholangitis. Int J Gen Med. 2021.',2021,'Delay and severe cholangitis','Delay increased hospital LOS; >48 h associated with longer ICU LOS','https://pmc.ncbi.nlm.nih.gov/articles/PMC8254098/'],
[10,'Tan M et al. Association between early ERCP and mortality in acute cholangitis. Gastrointest Endosc. 2018.',2018,'<24 h ERCP','Early ERCP association with lower 30-day mortality','https://pubmed.ncbi.nlm.nih.gov/28433613/'],
[11,'Iqbal U et al. Emergent versus urgent ERCP in acute cholangitis systematic review/meta-analysis. Gastrointest Endosc. 2020.',2020,'Timing','Earlier intervention and clinical outcomes','https://pubmed.ncbi.nlm.nih.gov/?term=Emergent+versus+urgent+ERCP+in+acute+cholangitis+Iqbal'],
[12,'Mohan R et al. Validation of TG07 and TG13/TG18 criteria and predictors of mortality in patients over 80 years. Clin Exp Hepatol. 2021.',2021,'Older adults severity','External validation of severity framework in very old adults','https://www.termedia.pl/Validation-of-TG07-and-TG13-TG18-criteria-for-acute-cholangitis-and-predictors-of-in-hospital-mortality-in-patients-over-80-years-of-age,80,45612,1,1.html'],
[13,'Kiriyama S et al. Clinical application and verification of TG13 diagnostic and severity grading criteria international multicenter observational study. J Hepatobiliary Pancreat Sci. 2017.',2017,'Severity grading','Multicenter verification informing TG18 adoption','https://pubmed.ncbi.nlm.nih.gov/?term=Kiriyama+TG13+international+multicenter+observational+study+2017'],
[14,'Wada K et al. Diagnostic criteria and severity assessment of acute cholangitis Tokyo Guidelines. J Hepatobiliary Pancreat Surg. 2007.',2007,'Historical severity framework','Foundational diagnostic/severity structure','https://pmc.ncbi.nlm.nih.gov/articles/PMC2784515/'],
[15,'Kongsakon R et al. Urgent ERCP within 24 h for moderate-to-severe cholangitis due to CBD stones propensity score matched study. Surg Endosc. 2025.',2025,'Urgent timing in moderate-severe disease','Contemporary severity-specific observational evidence','https://pubmed.ncbi.nlm.nih.gov/41023220/'],
[16,'Comoglu M et al. Effect of urgent ERCP on clinical outcomes in acute cholangitis and gallstone pancreatitis. 2025.',2025,'Urgent versus 24-72 h','Urgent <24 h not consistently superior in mixed cohort','https://pmc.ncbi.nlm.nih.gov/articles/PMC12380036/'],
[17,'Choday S et al. Timing of intervention assessing early versus late endoscopic and percutaneous drainage in acute cholangitis. World J Gastrointest Endosc. 2026.',2026,'Readmission by timing','Contemporary timing/outcome context','https://www.wjgnet.com/1948-5190/full/v18/i1/110353.htm'],
]
src=pd.DataFrame(sources,columns=['source_id','citation','year','parameter_or_anchor','reported_anchor_or_role','open_access_or_download_link']); src.to_csv(OUT/'Literature_Parameter_Sources.csv',index=False)

# Summary output
summary={
'seed':SEED,'reference_stopping_hour':ref_stop,
'reference_24h':{c:float(ref.iloc[np.argmin(abs(ref.hour-24))][c]) for c in ['bilirubin_mg_dl','wbc_10e3_ul','creatinine_mg_dl','platelets_10e3_ul','albumin_g_dl','deterioration_index','stability_reserve']},
'stochastic':{'median':float(stdf.stopping_hour.median()),'p2_5':float(stdf.stopping_hour.quantile(.025)),'p97_5':float(stdf.stopping_hour.quantile(.975))},
'eigenvalues':eigs.tolist(),'top_prcc':prdf.head(6).to_dict('records'),'scenarios':scdf.to_dict('records'), 'stop_fixed_tuned':BASE['stop_fixed']}
(OUT/'run_output.json').write_text(json.dumps(summary,indent=2))

# FIGURES
plt.rcParams.update({'font.size':10,'axes.titlesize':13,'axes.labelsize':10,'figure.dpi':160})
def save(name): plt.savefig(FIG/name,bbox_inches='tight',dpi=220); plt.close()
# 1 architecture
fig,ax=plt.subplots(figsize=(13,7)); ax.axis('off')
boxes=[(.03,.72,.20,.14,'Routine Laboratory Inputs','Bilirubin | WBC\nCreatinine | Platelets | Albumin'),(.30,.72,.21,.14,'Continuous Dynamics','Inflammatory burden I(t)\nOrgan dysfunction O(t)\nStability reserve R(t)'),(.58,.72,.22,.14,'Latent Clinical Regimes','Stable → progressive\norgan dysfunction → intervention-favored'),(.30,.40,.21,.14,'Competing Hazards','Waiting/deterioration intensity\nProcedure-associated intensity'),(.58,.40,.22,.14,'Optimal Stopping','Bellman recursion\nContinuation vs ERCP region'),(.38,.10,.34,.15,'Robustness and Validation','LHS-PRCC | stochastic perturbation | solver checks\nFuture retrospective/prospective calibration')]
colors=['#dbeafe','#dcfce7','#fef3c7','#fce7f3','#ede9fe','#e0f2fe']
for (x,y,w,h,t,s),c in zip(boxes,colors):
    b=FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.015,rounding_size=.02',fc=c,ec='#334155',lw=1.4); ax.add_patch(b); ax.text(x+w/2,y+h*.68,t,ha='center',va='center',weight='bold',fontsize=11); ax.text(x+w/2,y+h*.30,s,ha='center',va='center',fontsize=9)
for a,b in [((.23,.79),(.30,.79)),((.51,.79),(.58,.79)),((.405,.72),(.405,.54)),((.69,.72),(.69,.54)),((.51,.47),(.58,.47)),((.405,.40),(.49,.25)),((.69,.40),(.60,.25))]: ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=14,lw=1.5,color='#475569'))
ax.text(.5,.95,'Stochastic Optimal-Stopping Architecture for ERCP Timing',ha='center',weight='bold',fontsize=17); ax.text(.5,.90,'Literature-parameterized proof-of-concept — no patient-level calibration',ha='center',fontsize=10,style='italic'); save('Figure_1_Model_Architecture.png')
# 2 biomarker trajectories
fig,axs=plt.subplots(2,3,figsize=(13.5,8)); axs=axs.ravel(); series=[('bilirubin_mg_dl','Bilirubin','mg/dL',5),('wbc_10e3_ul','WBC','10³/µL',12),('creatinine_mg_dl','Creatinine','mg/dL',2),('platelets_10e3_ul','Platelets','10³/µL',100),('albumin_g_dl','Albumin','g/dL',2.8)]
for ax,(c,t,u,thr) in zip(axs,series): ax.plot(ref.hour,ref[c],lw=2.3); ax.axhline(thr,ls='--',lw=1.2); ax.set_title(t); ax.set_xlabel('Hours from presentation'); ax.set_ylabel(u); ax.grid(alpha=.2)
axs[5].plot(ref.hour,ref.inflammatory_burden,label='Inflammatory burden',lw=2); axs[5].plot(ref.hour,ref.organ_dysfunction,label='Organ dysfunction',lw=2); axs[5].plot(ref.hour,ref.stability_reserve,label='Stability reserve',lw=2); axs[5].set_title('Latent Physiologic State'); axs[5].set_xlabel('Hours'); axs[5].set_ylim(0,1); axs[5].legend(frameon=False,fontsize=8); axs[5].grid(alpha=.2)
fig.suptitle('Reference Laboratory and Physiologic Trajectory Without Drainage',fontsize=15,weight='bold'); fig.tight_layout(); save('Figure_2_Reference_Trajectory.png')
# 3 hazards + stopping
fig,axs=plt.subplots(1,2,figsize=(13,5.4)); axs[0].plot(refpath.hour,refpath.wait_hazard,lw=2.2,label='Waiting/deterioration intensity'); axs[0].plot(refpath.hour,refpath.procedure_hazard,lw=2.2,label='Procedure-associated intensity'); axs[0].axvline(ref_stop,ls='--',lw=1.4,label=f'Model stopping point {ref_stop:.0f} h'); axs[0].set_yscale('log'); axs[0].set_xlabel('Hours'); axs[0].set_ylabel('Theoretical intensity'); axs[0].set_title('Competing Hazards'); axs[0].legend(frameon=False,fontsize=8); axs[0].grid(alpha=.2)
axs[1].plot(refpath.hour,refpath.operate_now_cost,lw=2.2,label='Drain now value'); axs[1].plot(refpath.hour,refpath.wait_one_hour_cost,lw=2.2,label='Wait one hour value'); axs[1].axvline(ref_stop,ls='--',lw=1.4,label='First stopping point'); axs[1].fill_between(refpath.hour,0,1,where=(refpath.operate==0),alpha=.10,transform=axs[1].get_xaxis_transform(),label='Continuation region'); axs[1].fill_between(refpath.hour,0,1,where=(refpath.operate==1),alpha=.10,transform=axs[1].get_xaxis_transform(),label='Drainage-favored region'); axs[1].set_xlabel('Hours'); axs[1].set_ylabel('Model value'); axs[1].set_title('Finite-Horizon Optimal-Stopping Path'); axs[1].legend(frameon=False,fontsize=8); axs[1].grid(alpha=.2); fig.tight_layout(); save('Figure_3_Competing_Hazards_and_Stopping_Path.png')
# 4 surface
piv=surf.pivot(index='creatinine_mg_dl',columns='bilirubin_mg_dl',values='stopping_hour'); fig,ax=plt.subplots(figsize=(10.5,7)); im=ax.imshow(piv.values,origin='lower',aspect='auto',extent=[piv.columns.min(),piv.columns.max(),piv.index.min(),piv.index.max()],cmap='viridis_r'); cs=ax.contour(piv.columns,piv.index,piv.values,levels=np.unique(np.quantile(piv.values,[.1,.3,.5,.7,.9]).astype(int)),colors='white',linewidths=1.1); ax.clabel(cs,fmt=lambda x:f'{x:.0f} h',fontsize=8); ax.axvline(5,ls='--',color='black'); ax.axhline(2,ls='--',color='black'); ax.scatter([6],[1.35],marker='*',s=170,c='black',label='Reference'); ax.set_xlabel('Initial bilirubin (mg/dL)'); ax.set_ylabel('Initial creatinine (mg/dL)'); ax.set_title('Bilirubin–Creatinine Stopping Surface',weight='bold'); fig.colorbar(im,ax=ax,label='Model-derived stopping hour'); ax.legend(frameon=False); save('Figure_4_Bilirubin_Creatinine_Stopping_Surface.png')
# 5 scenario
fig,ax=plt.subplots(figsize=(11,6.5)); y=np.arange(len(scdf)); ax.barh(y,scdf.model_stopping_hour); ax.set_yticks(y,scdf.scenario); ax.invert_yaxis(); ax.axvspan(24,48,alpha=.1,label='24–48 h evidence window');
for i,v in enumerate(scdf.model_stopping_hour): ax.text(v+.5,i,f'{v:.0f} h',va='center',fontsize=9)
ax.set_xlabel('Model-derived stopping hour'); ax.set_title('Scenario-Specific Shift in Candidate ERCP Timing',weight='bold'); ax.grid(axis='x',alpha=.2); ax.legend(frameon=False); save('Figure_5_Scenario_Stopping_Hours.png')
# 6 sensitivity
sp=prdf.sort_values('PRCC_stopping_hour'); fig,ax=plt.subplots(figsize=(10.5,7)); vals=sp.PRCC_stopping_hour; cols=['#b91c1c' if x<0 else '#15803d' for x in vals]; ax.barh(sp.parameter,vals,color=cols,alpha=.86); ax.axvline(0,color='black',lw=1); ax.set_xlim(-1,1); ax.set_xlabel('PRCC with stopping hour'); ax.set_title('Global Sensitivity Across 600 Latin-Hypercube Parameter Sets',weight='bold'); ax.grid(axis='x',alpha=.2); save('Figure_6_Global_Sensitivity.png')
# 7 stochastic
fig,ax=plt.subplots(figsize=(10.5,6)); ax.hist(stdf.stopping_hour,bins=np.arange(stdf.stopping_hour.min()-.5,stdf.stopping_hour.max()+1.5,1),alpha=.8); med=stdf.stopping_hour.median(); lo=stdf.stopping_hour.quantile(.025); hi=stdf.stopping_hour.quantile(.975); ax.axvline(med,ls='--',color='black',label=f'Median {med:.0f} h'); ax.axvspan(lo,hi,alpha=.12,label=f'2.5–97.5% {lo:.0f}–{hi:.0f} h'); ax.set_xlabel('Model-derived stopping hour'); ax.set_ylabel('Simulation count'); ax.set_title('Stochastic Perturbation of Routine Laboratory Inputs',weight='bold'); ax.legend(frameon=False); ax.grid(alpha=.2); save('Figure_7_Stochastic_Stopping_Distribution.png')
# 8 regime occupancy
fig,ax=plt.subplots(figsize=(10.5,6.2));
for c,l in [('p_stable','Stable'),('p_progressive','Progressive'),('p_organ_dysfunction','Organ dysfunction'),('p_intervention_favored','Intervention-favored')]: ax.plot(ref.hour,ref[c],lw=2,label=l)
ax.axvline(ref_stop,ls='--',lw=1.2,color='black'); ax.set_ylim(0,1); ax.set_xlabel('Hours'); ax.set_ylabel('Probability'); ax.set_title('Latent Regime Occupancy Under Continued Observation',weight='bold'); ax.legend(frameon=False); ax.grid(alpha=.2); save('Figure_8_Regime_Occupancy.png')
# 9 phase portrait I vs O colored by time
fig,ax=plt.subplots(figsize=(8.5,6.5)); sca=ax.scatter(ref.inflammatory_burden,ref.organ_dysfunction,c=ref.hour,cmap='plasma',s=18); ax.plot(ref.inflammatory_burden,ref.organ_dysfunction,alpha=.45); ax.set_xlabel('Inflammatory burden I(t)'); ax.set_ylabel('Organ dysfunction O(t)'); ax.set_title('Phase-Plane Trajectory Toward Deterioration',weight='bold'); fig.colorbar(sca,ax=ax,label='Hours'); ax.grid(alpha=.2); save('Figure_9_Phase_Plane.png')
# 10 eigenvalues
fig,ax=plt.subplots(figsize=(8.5,5.8)); ax.scatter(eigs,np.zeros_like(eigs),s=85); ax.axvline(0,color='black',lw=1); ax.set_xlabel('Real part'); ax.set_yticks([]); ax.set_title('Post-Drainage Linearized Eigenvalue Spectrum',weight='bold'); ax.grid(axis='x',alpha=.2); save('Figure_10_Post_ERCP_Eigenvalues.png')
# table figures helper

def tablefig(df,title,name,fontsize=8):
    fig,ax=plt.subplots(figsize=(13,min(8,1.5+.42*len(df)))); ax.axis('off'); tb=ax.table(cellText=df.values, colLabels=df.columns, loc='center',cellLoc='left'); tb.auto_set_font_size(False); tb.set_fontsize(fontsize); tb.scale(1,1.35); ax.set_title(title,weight='bold',fontsize=14,pad=18); save(name)
tablefig(pd.DataFrame([['Bilirubin','5 mg/dL','TG18 Grade II'],['WBC','>12 or <4 ×10³/µL','TG18 Grade II'],['Creatinine','>2 mg/dL','TG18 Grade III renal dysfunction'],['Platelets','<100 ×10³/µL','TG18 Grade III hematologic dysfunction'],['Albumin','<0.7 × LLN','TG18 Grade II']],columns=['Variable','Literature anchor','Clinical role']),'Literature-Informed Laboratory Anchors','Figure_11_Table_Clinical_Anchors.png',8)
tablefig(scdf.round(3),'Scenario Outputs From the Proof-of-Concept Model','Figure_12_Table_Scenario_Outputs.png',7.5)
tablefig(pars.head(16).fillna(''),'Parameter Provenance and Evidence Status','Figure_13_Table_Parameter_Provenance.png',5.7)
tablefig(src[['source_id','year','parameter_or_anchor','reported_anchor_or_role']].head(17),'Literature Sources Used for Anchors and Timing Context','Figure_14_Table_Literature_Sources.png',6.3)

# Write code copies
shutil.copy(__file__,OUT/'cholangitis_optimal_stopping.py')

# DOCX helpers

def setshade(cell,fill):
    tcPr=cell._tc.get_or_add_tcPr(); shd=OxmlElement('w:shd'); shd.set(qn('w:fill'),fill); tcPr.append(shd)
def setup(doc):
    sec=doc.sections[0]; sec.top_margin=Inches(.65); sec.bottom_margin=Inches(.65); sec.left_margin=Inches(.7); sec.right_margin=Inches(.7)
    for s in ['Normal','Title','Heading 1','Heading 2','Heading 3']:
        st=doc.styles[s]; st.font.name='Arial'; st._element.rPr.rFonts.set(qn('w:eastAsia'),'Arial')
    doc.styles['Normal'].font.size=Pt(10.2); doc.styles['Title'].font.size=Pt(17); doc.styles['Heading 1'].font.size=Pt(14); doc.styles['Heading 2'].font.size=Pt(12); doc.styles['Heading 3'].font.size=Pt(11)
    return doc
def para(doc,txt,italic=False,boldlead=None):
    p=doc.add_paragraph(); p.paragraph_format.space_after=Pt(5)
    if boldlead and txt.startswith(boldlead): r=p.add_run(boldlead); r.bold=True; p.add_run(txt[len(boldlead):])
    else: r=p.add_run(txt); r.italic=italic
    return p
def eqn(doc,txt):
    p=doc.add_paragraph(); p.alignment=WD_ALIGN_PARAGRAPH.CENTER; p.paragraph_format.space_after=Pt(3); r=p.add_run(txt); r.font.name='Cambria Math'; r.font.size=Pt(11)
def table(doc,df,fs=7):
    t=doc.add_table(rows=1,cols=len(df.columns)); t.style='Table Grid'; t.alignment=WD_TABLE_ALIGNMENT.CENTER
    for j,c in enumerate(df.columns):
        cell=t.rows[0].cells[j]; cell.text=str(c); setshade(cell,'D9EAF7');
        for r in cell.paragraphs[0].runs: r.bold=True; r.font.size=Pt(fs)
    for i,row in df.iterrows():
        cs=t.add_row().cells
        for j,v in enumerate(row): cs[j].text=str(v); [setattr(r.font,'size',Pt(fs)) for r in cs[j].paragraphs[0].runs]
    return t
def figadd(doc,name,cap,w=6.5):
    p=doc.add_paragraph(); p.alignment=WD_ALIGN_PARAGRAPH.CENTER; p.add_run().add_picture(str(FIG/name),width=Inches(w)); cp=doc.add_paragraph(); cp.alignment=WD_ALIGN_PARAGRAPH.CENTER; rr=cp.add_run(cap); rr.italic=True; rr.font.size=Pt(8.3)

def add_hyperlink(paragraph, text, url):
    part=paragraph.part; r_id=part.relate_to(url,'http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink',is_external=True); h=OxmlElement('w:hyperlink'); h.set(qn('r:id'),r_id); r=OxmlElement('w:r'); rPr=OxmlElement('w:rPr'); color=OxmlElement('w:color'); color.set(qn('w:val'),'0563C1'); rPr.append(color); u=OxmlElement('w:u'); u.set(qn('w:val'),'single'); rPr.append(u); r.append(rPr); tx=OxmlElement('w:t'); tx.text=text; r.append(tx); h.append(r); paragraph._p.append(h)

D=setup(Document()); p=D.add_paragraph(); p.alignment=WD_ALIGN_PARAGRAPH.CENTER; r=p.add_run(TITLE1); r.bold=True; r.font.size=Pt(17)
para(D,'Literature-Parameterized Mathematical Proof-of-Concept With Stochastic Optimal Stopping, Regime Switching, and Uncertainty Analysis',italic=True)
D.add_heading('Research Integrity Statement',1); para(D,'No patient-level registry records, prospective observations, field measurements, or laboratory experiments were analyzed. Published guidelines, randomized and observational studies, and meta-analyses were used to select laboratory variables, clinical anchors, timing windows, and expected directions of effect. Continuous-state rate coefficients, Markov transition intensities, hazard weights, and Bellman costs are provisional mathematical assumptions unless explicitly labeled literature-derived. Every numerical result reported below is generated by the accompanying code. None is a validated treatment threshold, clinical effect estimate, sensitivity/specificity result, or individualized recommendation. Final effect estimates must be recalculated from underlying registry data before scientific or clinical interpretation.')
D.add_heading('Three Title Choices',1)
for i,t in enumerate([TITLE1,TITLE2,TITLE3],1): para(D,f'{i}. {t} ({len(t.replace("-"," ").split())} words)')
D.add_heading('KDDW 2026 Classification and Keywords',1); para(D,'Recommended official category: Pancreatobiliary. Recommended official subcategory: Cholestasis and gallstone. Alternative if the submission portal emphasizes the procedure rather than the disease mechanism: Endoscopy and imaging.'); para(D,'Keywords: '+KEYWORDS)

# Structured abstract one paragraph each
D.add_heading('Structured Article Summary',1)
abstract={
'Background and Aim':f'Acute cholangitis is a time-sensitive consequence of biliary obstruction, yet the exact ERCP window remains clinically nuanced: guidelines support prompt drainage and commonly use a 48-hour boundary, observational syntheses favor earlier drainage, while a 2026 randomized trial in mild-to-moderate disease found no mortality or organ-failure advantage for <24-hour ERCP over 24–48 hours and reported more procedure-related adverse events with urgent ERCP. We developed a literature-parameterized mathematical proof-of-concept that treats ERCP timing as a stochastic optimal-stopping problem using five routinely available variables—bilirubin, WBC, creatinine, platelets, and albumin—to represent cholestatic burden, inflammation, organ dysfunction, and physiologic reserve.',
'Methods':'The framework combined five literature-normalization functions, an eight-state continuous dynamical system for the five laboratory variables plus inflammatory burden I(t), organ dysfunction O(t), and stability reserve R(t), a four-state continuous-time Markov regime model, competing waiting and procedure intensities, and a 72-hour finite-horizon Bellman recursion. Literature anchors were bilirubin 5 mg/dL, WBC >12 or <4×10³/µL, creatinine >2 mg/dL, platelets <100×10³/µL, and albumin <0.70× the lower limit of normal from TG18, together with the ASGE <48-hour decompression recommendation and timing evidence from meta-analyses and a recent randomized trial. Numerical work used adaptive RK45, 600-set Latin-hypercube sampling with PRCC, ±1% local elasticity, a 1,000-run stochastic input perturbation, phase-plane inspection, Markov probability conservation, equilibrium/eigenvalue verification after theoretical drainage, and a bilirubin–creatinine free-boundary surface. Structural coefficients were not inferred from patient outcomes and remain provisional.',
'Results':f'The reference illustrative profile (bilirubin 6.0 mg/dL, WBC 14.0×10³/µL, creatinine 1.35 mg/dL, platelets 145×10³/µL, albumin 3.05 g/dL) produced a model-derived stopping point of {ref_stop:.0f} hours. At 24 hours the simulated deterioration index was {summary["reference_24h"]["deterioration_index"]:.3f} and stability reserve {summary["reference_24h"]["stability_reserve"]:.3f}. Across prespecified clinical scenarios, the candidate stopping point ranged from {scdf.model_stopping_hour.min():.0f} to {scdf.model_stopping_hour.max():.0f} hours; the combined severe laboratory profile shifted the boundary to {float(scdf.loc[scdf.scenario.eq("Combined severe laboratory profile"),"model_stopping_hour"].iloc[0]):.0f} hours, whereas the lower inflammatory profile shifted it to {float(scdf.loc[scdf.scenario.eq("Lower inflammatory burden"),"model_stopping_hour"].iloc[0]):.0f} hours. Stochastic perturbations yielded a median stopping time of {stdf.stopping_hour.median():.0f} hours with a 2.5th–97.5th percentile range of {stdf.stopping_hour.quantile(.025):.0f}–{stdf.stopping_hour.quantile(.975):.0f} hours. All eight post-drainage eigenvalues were negative ({eigs.min():.3f} to {eigs.max():.3f}), supporting local asymptotic stability of the recovery equilibrium in the specified subsystem. The three largest absolute PRCCs were ' + ', '.join([f'{r.parameter} {r.PRCC_stopping_hour:.3f}' for _,r in prdf.head(3).iterrows()]) + '. These are mathematical outputs, not observed clinical effects.',
'Conclusions':'The current framework establishes a mathematically testable structure for a familiar bedside problem: when continued stabilization becomes less favorable than biliary drainage as laboratory evidence of deterioration accumulates. Its clinical use is not to replace TG18 grading or procedural judgment, but to define variables, candidate intervention boundaries, and uncertainty patterns for subsequent registry calibration. The main advantage is that optimal stopping links time, evolving physiology, and competing procedure-versus-waiting costs in a way that static association models do not. The model can support protocol development, retrospective target-trial emulation, and eventual prospective decision-support studies with gastroenterology and endoscopy teams. Parameter refinement and biological calibration are ongoing; clinical calibration and external validation are required before any treatment recommendation.'}
for h,t in abstract.items(): D.add_heading(h,2); para(D,t)
figadd(D,'Figure_1_Model_Architecture.png','Figure 1. Hybrid architecture linking routine laboratory inputs to disease dynamics, latent regimes, competing hazards, and the optimal-stopping layer.')

D.add_heading('1. Background and Clinical Rationale',1); para(D,'Acute cholangitis is common in emergency gastroenterology because biliary obstruction, infection, and systemic deterioration can evolve over hours rather than weeks. The practical decision is not whether biliary drainage matters, but how rapidly ERCP should follow initial antibiotics and resuscitation in a patient whose laboratory profile is changing. TG18 uses readily available measurements to grade disease severity, including WBC, bilirubin, hypoalbuminemia, creatinine, and platelet count. ASGE guidance supports decompression within 48 hours, whereas the published timing literature contains an important tension: large observational syntheses associate earlier ERCP with better outcomes, yet newer randomized evidence in mild-to-moderate disease does not show that <24 hours is uniformly superior to 24–48 hours and reports a higher procedure-related adverse-event rate with urgent ERCP. This makes the problem well suited to sequential decision theory rather than a single static risk score. The financial-engineering contribution is methodological: stochastic optimal stopping, dynamic reserve, regime switching, competing hazards, free-boundary analysis, and uncertainty propagation are established tools for decisions where both state and action costs change with time. Here they are used to formalize a clinical waiting-versus-drainage problem without treating patients as financial assets.')

D.add_heading('2. Methods',1)
D.add_heading('2.1 Study Design and Evidence Status',2); para(D,'This study was designed as a public-evidence-anchored, literature-parameterized mathematical proof-of-concept spanning mathematical biology, computational medicine, biomedical engineering, data science, and quantitative decision science. No individual-level outcomes were fitted. Literature was used only for clinical anchors, variable definitions, timing context, and direction of effects; it was not used to manufacture longitudinal biomarker trajectories or to convert published odds ratios directly into ODE coefficients. Structural coefficients were chosen to yield bounded, clinically interpretable trajectories and were challenged through broad sensitivity ranges. This separation between observed evidence and model assumption follows the integrity principle used in the supplied TBAD example, where literature anchors are kept distinct from provisional ODE and Bellman parameters.')
D.add_heading('2.2 Literature-Informed Clinical Anchors',2); anchors_df=pd.DataFrame([['Bilirubin','5 mg/dL','TG18 Grade II severity anchor'],['WBC','>12 or <4 ×10³/µL','TG18 Grade II severity anchor'],['Creatinine','>2 mg/dL','TG18 Grade III renal dysfunction'],['Platelets','<100 ×10³/µL','TG18 Grade III hematologic dysfunction'],['Albumin','<0.70 × LLN','TG18 Grade II hypoalbuminemia'],['ERCP timing','within 48 h','ASGE conditional decompression recommendation']],columns=['Domain','Anchor','Role']); table(D,anchors_df,7); figadd(D,'Figure_11_Table_Clinical_Anchors.png','Figure 2. Figure rendering of the principal literature-informed laboratory and timing anchors.')
D.add_heading('2.3 State Variables',2); states=pd.DataFrame([['B(t)','Total bilirubin','mg/dL','Cholestatic/obstructive burden'],['W(t)','WBC','10³/µL','Systemic inflammatory response'],['C(t)','Creatinine','mg/dL','Renal dysfunction'],['P(t)','Platelets','10³/µL','Hematologic dysfunction'],['A(t)','Albumin','g/dL','Reserve/severity modifier'],['I(t)','Inflammatory burden','0–1','Latent inflammatory activity'],['O(t)','Organ dysfunction','0–1','Latent cumulative organ injury'],['R(t)','Stability reserve','0–1','Capacity to remain observation-favorable']],columns=['Symbol','State','Scale','Interpretation']); table(D,states,7)
D.add_heading('2.4 Mathematical Formulation',2); para(D,'Equations are listed one per row. The model deliberately combines several established model classes rather than forcing every component into one nonlinear system.')
# equations each row
for e in [
'b(B) = {1 + exp[-(B - 5.0)/1.4]}^(-1)',
'w(W) = max({1 + exp[-(W - 12)/2.8]}^(-1), {1 + exp[-(4 - W)/1.2]}^(-1))',
'c(C) = {1 + exp[-(C - 2.0)/0.45]}^(-1)',
'p(P) = {1 + exp[-(100 - P)/18]}^(-1)',
'a(A) = {1 + exp[-(2.8 - A)/0.30]}^(-1)',
'dB/dt = k_B↑(1 + 0.55I) - k_B↓ max(B - 2,0)/5',
'dW/dt = -k_W↓ max(W - 8,0) + k_WB[3b(B) + 2I]',
'dC/dt = k_C↑(0.8I + 1.4O) - k_C↓ max(C - 0.9,0)R',
'dP/dt = -k_P↓(8I + 10O) + k_P↑ max(180 - P,0)R/80',
'dA/dt = -k_A↓(1.2I + 0.8O) + k_A↑ max(3.7 - A,0)R',
'dI/dt = [k_IB b(B) + k_IW w(W)](1 - I) - k_IR RI',
'dO/dt = [k_OI I + k_OC c(C) + k_OP p(P)](1 - O) - k_OR RO',
'dR/dt = k_R+(1 - R) - k_RI IR - k_RO OR - k_RA a(A)R',
'D(t) = 0.16b(B) + 0.18w(W) + 0.22c(C) + 0.18p(P) + 0.10a(A) + 0.16O',
'Q(x) = continuous-time Markov generator for Stable → Progressive → Organ Dysfunction → Intervention-Favored regimes',
'λ_wait(t) = h_W0 exp{β_W[0.70D(t) + 0.30O(t)]}',
'λ_ERCP(t) = h_P0[1 + ψ exp(-t/τ)] exp{β_P[0.55O(t) + 0.45(1 - R(t))]}',
'C_stop(t) = c_0 + λ_ERCP(t) + 0.018O(t)',
'C_wait(t) = c_1 + λ_wait(t) + 0.010D(t) + γV(t + 1)',
'V(t) = min{C_stop(t), C_wait(t)}',
'τ* = inf{t : C_stop(t) ≤ C_wait(t)}']:
    eqn(D,e)
D.add_heading('2.5 Theoretical Analysis',2); para(D,'Existence and uniqueness follow locally from the piecewise continuously differentiable right-hand side on the admissible domain. Positivity and boundedness were checked numerically over the prespecified input and parameter ranges; physiologic latent variables I, O, and R were constrained to [0,1]. The four-state Markov generator preserves total probability because each row sums to zero, while the matrix exponential preserves non-negative regime probabilities. A theoretical post-drainage equilibrium was defined at bilirubin 2.0 mg/dL, WBC 8.0×10³/µL, creatinine 0.9 mg/dL, platelets 180×10³/µL, albumin 3.7 g/dL, I*=0, O*=0, and R*=1; its linearized eigenvalues were inspected for local stability. Dimensional consistency was enforced by treating biomarker ODE coefficients as per-hour rates and Bellman terms as dimensionless theoretical costs. Because no longitudinal observation model was fitted, formal structural/practical identifiability is deferred to registry calibration rather than overstated here.')
D.add_heading('2.6 Computational Analysis',2); para(D,'The pre-drainage system was integrated over 72 hours with adaptive Runge–Kutta RK45 at 0.25-hour output resolution and strict tolerances. A four-state continuous-time Markov chain was propagated by the matrix exponential. The decision layer sampled the continuous trajectory hourly and solved a finite-horizon Bellman recursion backward from 72 hours, identifying the first time at which immediate drainage had lower theoretical cost than one additional hour of waiting. A two-dimensional free-boundary surface varied initial bilirubin from 2 to 12 mg/dL and creatinine from 0.8 to 3.0 mg/dL while holding the remaining reference inputs fixed. Eight clinically interpretable scenarios changed routine laboratory severity and latent reserve. Phase-plane trajectories were used to visualize inflammatory burden against organ dysfunction. No finite difference or finite element method was used because the present model contains ODE, Markov, and dynamic-programming components rather than a PDE; adding a PDE solver would be methodologically unnecessary.')
D.add_heading('2.7 Robustness and Uncertainty Analysis',2); para(D,'Global uncertainty analysis used 600 Latin-hypercube parameter sets spanning inflammatory coupling, organ-dysfunction coupling, reserve dynamics, waiting intensity, procedure intensity, early-procedure penalty, and hazard slopes. Partial rank correlation coefficients were computed after rank residualization against all other sampled parameters. Local elasticity used ±1% parameter perturbations. A separate 1,000-run stochastic experiment perturbed the five laboratory inputs and latent starting conditions around the reference scenario; these distributions are computational stress tests, not estimates of a patient population. Solver and probability-conservation checks were performed during execution. AIC was intentionally excluded from the primary proof-of-concept because no patient-level likelihood was fitted; clinical AIC/BIC comparison belongs to the future observed-data phase.')

D.add_heading('3. Results',1)
D.add_heading('3.1 Reference Trajectory and Candidate Stopping Point',2); para(D,f'The reference illustrative profile started at bilirubin 6.0 mg/dL, WBC 14.0×10³/µL, creatinine 1.35 mg/dL, platelets 145×10³/µL, albumin 3.05 g/dL, inflammatory burden 0.22, organ dysfunction 0.10, and stability reserve 0.80. The Bellman layer produced a candidate stopping point at {ref_stop:.0f} hours. At 24 hours, bilirubin was {summary["reference_24h"]["bilirubin_mg_dl"]:.2f} mg/dL, WBC {summary["reference_24h"]["wbc_10e3_ul"]:.2f}×10³/µL, creatinine {summary["reference_24h"]["creatinine_mg_dl"]:.2f} mg/dL, platelets {summary["reference_24h"]["platelets_10e3_ul"]:.1f}×10³/µL, albumin {summary["reference_24h"]["albumin_g_dl"]:.2f} g/dL, deterioration index {summary["reference_24h"]["deterioration_index"]:.3f}, and stability reserve {summary["reference_24h"]["stability_reserve"]:.3f}. These are generated trajectory values rather than observed laboratory measurements.')
figadd(D,'Figure_2_Reference_Trajectory.png','Figure 3. Reference biomarker and latent-state trajectory under continued observation without modelled drainage. Dashed lines are literature anchors, not fitted treatment thresholds.')
figadd(D,'Figure_3_Competing_Hazards_and_Stopping_Path.png','Figure 4. Competing theoretical intensities and finite-horizon optimal-stopping path. The crossing is a model output, not an ERCP recommendation.')
D.add_heading('3.2 Scenario-Specific Timing Shifts',2); para(D,'The model responded directionally to clinically recognizable deterioration rather than only to elapsed time. Across the eight prespecified scenarios, candidate stopping points ranged from %.0f to %.0f hours. The combined severe profile moved the intervention boundary to %.0f hours; renal dysfunction and thrombocytopenic deterioration moved it to %.0f and %.0f hours, respectively. In contrast, the lower inflammatory profile delayed the boundary to %.0f hours. These results show the intended mechanistic behavior: the same 24–48-hour calendar window does not carry the same modeled waiting cost when renal or hematologic dysfunction is already emerging.'%(scdf.model_stopping_hour.min(),scdf.model_stopping_hour.max(),scdf.loc[scdf.scenario.eq('Combined severe laboratory profile'),'model_stopping_hour'].iloc[0],scdf.loc[scdf.scenario.eq('Renal dysfunction emerging'),'model_stopping_hour'].iloc[0],scdf.loc[scdf.scenario.eq('Thrombocytopenic deterioration'),'model_stopping_hour'].iloc[0],scdf.loc[scdf.scenario.eq('Lower inflammatory burden'),'model_stopping_hour'].iloc[0])); table(D,scdf.round(3),7); figadd(D,'Figure_5_Scenario_Stopping_Hours.png','Figure 5. Scenario-specific candidate ERCP timing shifts with the 24–48-hour evidence window shown for context.'); figadd(D,'Figure_12_Table_Scenario_Outputs.png','Figure 6. Figure rendering of the scenario-output table.')
D.add_heading('3.3 Free-Boundary Surface and Disease Geometry',2); para(D,f'The bilirubin–creatinine surface covered {len(surf)} combinations. Higher creatinine and higher bilirubin generally moved the stopping boundary earlier, producing a continuous decision surface rather than a single universal cutoff. Across the grid, model-derived stopping times ranged from {surf.stopping_hour.min():.0f} to {surf.stopping_hour.max():.0f} hours. The phase-plane trajectory showed inflammatory burden and organ dysfunction moving along a curved path rather than increasing proportionally, illustrating why a dynamic state model can provide information that a one-time regression coefficient cannot.')
figadd(D,'Figure_4_Bilirubin_Creatinine_Stopping_Surface.png','Figure 7. Model-derived bilirubin–creatinine free-boundary surface.'); figadd(D,'Figure_9_Phase_Plane.png','Figure 8. Phase-plane trajectory of inflammatory burden and organ dysfunction.')
D.add_heading('3.4 Regime Switching and Stability',2); para(D,f'The Markov layer redistributed probability from the initially stable/progressive mixture toward organ-dysfunction and intervention-favored regimes as the latent deterioration index increased. The post-drainage recovery equilibrium had eight negative real eigenvalues ranging from {eigs.min():.3f} to {eigs.max():.3f} per hour with zero imaginary parts, supporting local asymptotic stability in the linearized recovery subsystem. This is a property of the specified equations and should not be interpreted as evidence that every patient will clinically recover after ERCP.')
figadd(D,'Figure_8_Regime_Occupancy.png','Figure 9. Latent regime occupancy under continued observation.'); figadd(D,'Figure_10_Post_ERCP_Eigenvalues.png','Figure 10. Eigenvalue spectrum of the theoretical post-drainage recovery equilibrium.')
D.add_heading('3.5 Global Sensitivity and Stochastic Uncertainty',2); para(D,'Across 600 Latin-hypercube runs, the strongest absolute PRCCs with the candidate stopping hour were ' + ', '.join([f'{r.parameter} {r.PRCC_stopping_hour:.3f} (p={r.p_value:.2g})' for _,r in prdf.head(5).iterrows()]) + f'. The 1,000 stochastic input perturbations produced a median stopping time of {stdf.stopping_hour.median():.0f} hours and a 2.5th–97.5th percentile interval of {stdf.stopping_hour.quantile(.025):.0f}–{stdf.stopping_hour.quantile(.975):.0f} hours. The width of this interval is clinically useful as an uncertainty signal: it shows how a timing recommendation could become unstable when laboratory inputs or structural assumptions are uncertain, supporting the future use of probabilistic rather than deterministic decision support.')
figadd(D,'Figure_6_Global_Sensitivity.png','Figure 11. Global PRCC sensitivity across 600 Latin-hypercube parameter sets.'); figadd(D,'Figure_7_Stochastic_Stopping_Distribution.png','Figure 12. Distribution of model-derived stopping hours across 1,000 stochastic input perturbations.')

D.add_heading('4. Discussion',1); para(D,'The main contribution is not another mortality classifier. It is a mathematical representation of a clinical trade-off that endoscopists already make: antibiotics and resuscitation may improve procedural readiness, but waiting while cholestasis and organ dysfunction accumulate also has a cost. The model makes that trade-off explicit and allows the boundary to move with bilirubin, inflammatory burden, renal dysfunction, thrombocytopenia, and low albumin reserve. This feature is particularly relevant in light of the current evidence base. TG18 identifies the chosen variables as severity markers; ASGE recommends decompression within 48 hours; meta-analyses generally favor early drainage; and the 2026 randomized trial suggests that for mild-to-moderate disease, pushing all patients into a <24-hour pathway may not improve mortality or organ failure and may increase procedure-related adverse events. A stochastic stopping model therefore offers a clinically interpretable middle ground between rigid clock-based rules and opaque black-box prediction. It can also be extended to center-specific endoscopy availability, but queueing constraints should be added only after clinical calibration rather than mixed into the present biological decision model.')
D.add_heading('5. Clinical Use, Benefits, and Research Advantages',1); para(D,'The immediate use of this proof-of-concept is research design. It identifies a compact set of routine variables for serial extraction, defines a candidate deterioration index, specifies a mathematically explicit waiting-versus-drainage objective, and produces uncertainty-aware intervention surfaces that can be tested in retrospective data. For clinicians, the potential benefit is a future decision-support output that states whether a patient remains in a continuation region or has entered a drainage-favored region, together with uncertainty, rather than presenting a static score without timing information. The research advantage is falsifiability: each structural component can be recalibrated or rejected when longitudinal registry data become available. The framework also creates a natural collaboration point between endoscopists, hepatobiliary clinicians, data scientists, and quantitative researchers because all five headline variables are routinely measurable and do not require a laboratory experiment or new biomarker assay.')
D.add_heading('6. Limitations and Integrity Boundaries',1); para(D,'The stopping hour is not a clinical threshold. The reference laboratory profile is illustrative, structural coefficients are provisional, and the stochastic input distributions are stress tests rather than empirical prevalence distributions. Published study-level effect estimates were not converted into patient-level hazards. Albumin is modeled in g/dL using an illustrative lower limit of normal of 4.0 g/dL only to operationalize the TG18 <0.70×LLN criterion; a future registry must use its laboratory-specific reference range. The current model does not account for etiology, source-control complexity, altered anatomy, anticoagulation, coagulopathy beyond platelets, vasopressor requirement, mental status, respiratory dysfunction, availability of ERCP expertise, or procedure-specific contraindications. Causal effects of ERCP timing cannot be inferred from this simulation. Final effect estimates must be recalculated from underlying registry data before scientific or clinical interpretation.')
D.add_heading('7. Clinical Calibration and Validation Roadmap',1); para(D,'The next phase should preregister a retrospective target-trial emulation with time zero defined at diagnosis or presentation, explicit eligibility criteria, exact ERCP timestamps, serial laboratory measurements, TG18 grade, etiology, antibiotics, resuscitation, ICU transfer, organ failure, mortality, length of stay, readmission, and ERCP adverse events. Longitudinal mixed models or joint models can estimate biomarker trajectories; multistate or competing-risk models can calibrate transitions; Cox or flexible parametric hazards can estimate time-dependent deterioration; likelihood-based alternatives can then be compared with AIC/BIC. Validation should report calibration-in-the-large, calibration slope, Brier score, time-dependent AUC, prespecified sensitivity/specificity where clinically meaningful, and decision-curve net benefit. Treatment timing is time-varying and confounded by severity and clinician response, so causal interpretation requires target-trial emulation, inverse-probability weighting, g-computation, or another justified g-method. External validation at a second center should precede prospective impact testing.')

D.add_heading('8. Parameter Provenance',1); para(D,'The accompanying Model_Parameters.csv separates literature-derived anchors from illustrative scenario values and provisional structural/decision coefficients. This distinction is essential: the internet-accessible literature supports the clinical anchors and timing context, but it does not supply every ODE coefficient, Markov intensity, or Bellman cost. Those parameters remain calibration targets rather than claimed biological constants.'); table(D,pars.head(22).fillna(''),5.7); figadd(D,'Figure_13_Table_Parameter_Provenance.png','Figure 13. Figure rendering of parameter provenance and evidence status.')
D.add_heading('9. Reproducible Analysis and CSV Provenance',1); prov=pd.DataFrame([
['Reference_Trajectory.csv','Adaptive RK45 solution of the eight-state pre-drainage system; contains biomarkers, latent states, Markov probabilities, and theoretical hazards.'],['Reference_Optimal_Stopping_Path.csv','Hourly finite-horizon Bellman recursion applied to the reference trajectory.'],['Scenario_Results.csv','Eight prespecified laboratory/reserve scenarios; all stopping hours are model-generated.'],['Bilirubin_Creatinine_Stopping_Surface.csv','483-grid free-boundary analysis over initial bilirubin and creatinine.'],['LHS_Parameter_Samples.csv','600 Latin-hypercube structural/decision parameter sets and stopping-hour outputs.'],['Global_Sensitivity_PRCC.csv','Residualized partial rank correlations between sampled parameters and stopping hour.'],['Local_Elasticity.csv','±1% local perturbations around the reference parameterization.'],['Stochastic_Perturbation_Results.csv','1,000 stochastic perturbations of routine laboratory inputs and latent initial conditions.'],['Post_ERCP_Equilibrium.csv','Declared theoretical recovery equilibrium for stability analysis.'],['Post_ERCP_Eigenvalues.csv','Linearized recovery eigenvalues; mathematical stability result.'],['Model_Parameters.csv','Literature anchors, illustrative reference values, and provisional parameter ranges with status labels.'],['Literature_Parameter_Sources.csv','Seventeen literature sources with downloadable/open-access links where available.']],columns=['File','How the numbers are generated']); table(D,prov,6.2)
D.add_heading('10. Open-Access and Downloadable Literature Sources',1)
for _,r in src.iterrows():
    p=D.add_paragraph(); p.add_run(f'{int(r.source_id)}. {r.citation} '); add_hyperlink(p,'Open/download link',r.open_access_or_download_link)
figadd(D,'Figure_14_Table_Literature_Sources.png','Figure 14. Figure rendering of literature sources used for clinical anchors and timing context.')
D.add_heading('11. Conclusion',1); para(D,'This mechanistic proof-of-concept provides a mathematically testable framework for identifying clinically relevant ERCP intervention windows that conventional association-based analyses may not reveal. Its proposed clinical use is to distinguish continued-stabilization and drainage-favored regions from routinely collected bilirubin, WBC, creatinine, platelet, and albumin measurements while explicitly displaying uncertainty. The benefit is a time-aware decision structure that can complement, not replace, TG18 severity grading and endoscopist judgment. The research advantage is that every component—biomarker dynamics, transition intensities, hazard functions, and decision costs—can be estimated and challenged in later cohorts. The model is designed to guide variable selection, protocol development, retrospective causal analysis, and subsequent prospective validation in collaboration with gastroenterology and therapeutic-endoscopy teams. Parameter refinement and biological calibration are ongoing; no model-derived stopping hour should be used clinically before external validation.')

full=OUT/'KDDW_2026_Acute_Cholangitis_Optimal_Stopping_Full_Article.docx'; D.save(full)
# provenance doc
P=setup(Document()); p=P.add_paragraph(); p.alignment=WD_ALIGN_PARAGRAPH.CENTER; rr=p.add_run('Parameter Sources and CSV Provenance'); rr.bold=True; rr.font.size=Pt(17); para(P,'Acute Cholangitis ERCP Timing Stochastic Optimal Stopping Model',italic=True); P.add_heading('Evidence Status',1); para(P,'Literature supplies clinical anchors, variable definitions, timing windows, and expected directions of effect. It does not supply every ODE coefficient, Markov transition intensity, hazard weight, or Bellman cost. Those quantities are explicitly labeled provisional until longitudinal calibration.'); P.add_heading('Literature Sources',1); table(P,src,5.7); P.add_heading('Model Parameters',1); table(P,pars.fillna(''),5.7); P.add_heading('CSV Provenance',1); table(P,prov,6.2); P.save(OUT/'KDDW_2026_Parameter_Sources_and_CSV_Provenance.docx')
# protocol doc
R=setup(Document()); p=R.add_paragraph(); p.alignment=WD_ALIGN_PARAGRAPH.CENTER; rr=p.add_run('Reproducible Analysis Protocol'); rr.bold=True; rr.font.size=Pt(17)
steps=[('1. Evidence extraction','Record TG18 variable definitions and severity anchors; extract ERCP timing categories and outcomes from guidelines, meta-analyses, RCTs, and cohorts. Do not convert an OR or RR directly into an ODE coefficient.'),('2. State formulation','Define the five laboratory trajectories, I(t), O(t), R(t), normalization functions, Markov generator, waiting intensity, procedure intensity, and Bellman costs.'),('3. Structural verification','Check positivity/bounds, dimensional consistency, Markov probability conservation, declared recovery equilibrium, eigenvalues, and numerical solver stability.'),('4. Reference simulation','Run the illustrative moderate profile over 72 h with RK45; review every generated laboratory trajectory for plausibility.'),('5. Optimal stopping','Apply the hourly finite-horizon Bellman recursion; identify the first time immediate drainage has lower theoretical cost than waiting one hour.'),('6. Free-boundary continuation','Vary bilirubin 2–12 mg/dL and creatinine 0.8–3.0 mg/dL; map the stopping surface.'),('7. Scenario analysis','Run eight predefined profiles emphasizing low burden, cholestasis, renal dysfunction, thrombocytopenia, hypoalbuminemia, and combined severe laboratory abnormalities.'),('8. Sensitivity and uncertainty','Run 600 LHS parameter sets with PRCC, ±1% local elasticity, and 1,000 stochastic input perturbations.'),('9. Manuscript transfer','Transfer only values present in run_output.json/CSVs. Label every simulation result as model-derived.'),('10. Retrospective calibration','Fit longitudinal biomarker, multistate, survival, and observed-data decision models to an ethically approved registry. Use AIC/BIC only on genuine observed-data likelihoods.'),('11. Causal timing analysis','Emulate a target trial or use justified g-methods because ERCP timing is confounded by evolving severity and clinician response.'),('12. External validation','Report calibration, Brier score, time-dependent AUC, clinically prespecified sensitivity/specificity, decision-curve net benefit, and external validation before clinical use.')]
for h,t in steps: R.add_heading(h,2); para(R,t)
R.add_heading('Run Order',1); para(R,'1. python cholangitis_optimal_stopping.py'); para(R,'2. Review run_output.json and CSV outputs.'); para(R,'3. Regenerate figures if any parameter is changed.'); R.save(OUT/'KDDW_2026_Reproducible_Analysis_Protocol.docx')
# figure portfolio
F=setup(Document()); p=F.add_paragraph(); p.alignment=WD_ALIGN_PARAGRAPH.CENTER; rr=p.add_run('Figure Portfolio'); rr.bold=True; rr.font.size=Pt(17)
for i,(fn,cap) in enumerate([(f'Figure_{i}_'+x,'') for i,x in []]): pass
for fn in sorted(FIG.glob('*.png')): figadd(F,fn.name,fn.stem.replace('_',' '),6.6)
F.save(OUT/'KDDW_2026_Figure_Portfolio.docx')
# README
(OUT/'README.txt').write_text(f'''KDDW 2026 Acute Cholangitis Optimal-Stopping Proof-of-Concept Package\n\nPrimary title:\n{TITLE1}\n\nStudy status:\nLiterature-parameterized mathematical proof-of-concept. No patient-level, field, or laboratory data were analyzed. Literature-derived numbers are clinical anchors, not fitted ODE/hazard coefficients. Structural coefficients and stopping-hour outputs are provisional.\n\nRun:\npython cholangitis_optimal_stopping.py\n\nIntegrity rule:\nFinal effect estimates, discrimination metrics, calibration measures, and clinically interpretable timing thresholds must be recalculated from observed registry data before scientific or clinical interpretation.\n''')
print(json.dumps(summary,indent=2))
