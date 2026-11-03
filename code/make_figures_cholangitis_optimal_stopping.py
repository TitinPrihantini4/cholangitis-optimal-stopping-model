from pathlib import Path
import textwrap, zipfile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Patch

# ============================================================
# REVISED FIGURE-ONLY SCRIPT
# Recreates Figures 1–14 from existing numerical CSV outputs.
# It does NOT rerun the ODE, LHS/PRCC, or stochastic simulations,
# so the numerical essence/results are unchanged.
# ============================================================

ROOT = Path('/mnt/data/KDDW_2026_Cholangitis_Optimal_Stopping_Package')
FIG = ROOT / 'figures_revised'
FIG.mkdir(parents=True, exist_ok=True)

ref = pd.read_csv(ROOT / 'Reference_Trajectory.csv')
refpath = pd.read_csv(ROOT / 'Reference_Optimal_Stopping_Path.csv')
surf = pd.read_csv(ROOT / 'Bilirubin_Creatinine_Stopping_Surface.csv')
scdf = pd.read_csv(ROOT / 'Scenario_Results.csv')
prdf = pd.read_csv(ROOT / 'Global_Sensitivity_PRCC.csv')
stdf = pd.read_csv(ROOT / 'Stochastic_Perturbation_Results.csv')
eigdf = pd.read_csv(ROOT / 'Post_ERCP_Eigenvalues.csv')
eqdf = pd.read_csv(ROOT / 'Post_ERCP_Equilibrium.csv')
pars = pd.read_csv(ROOT / 'Model_Parameters.csv')
src = pd.read_csv(ROOT / 'Literature_Parameter_Sources.csv')

ref_stop = float(scdf.loc[scdf['scenario'].eq('Reference moderate profile'),'model_stopping_hour'].iloc[0])

# Muted, publication-style palette: deliberately avoids bright primary RGB colors.
INK = '#30343b'; GRID = '#9aa0a6'; THRESH = '#6f7782'
MUTED = ['#7b6d8d','#8b7d6b','#667c75','#8a6f6f','#6f7488','#897c91']
PURPLE = '#76658a'; TAUPE = '#8a7868'; SLATE = '#687987'; MAUVE = '#8b6f82'
SAND = '#a18d72'; DUST = '#8b7786'; POS = '#6f7896'; NEG = '#9a7a70'

plt.rcParams.update({
    'font.size': 10, 'axes.titlesize': 13, 'axes.labelsize': 10,
    'figure.dpi': 160, 'axes.edgecolor': INK, 'axes.linewidth': 0.9,
    'font.family': 'DejaVu Sans'
})

def save(name):
    plt.savefig(FIG / name, bbox_inches='tight', dpi=240, facecolor='white')
    plt.close()

def style_ax(ax):
    ax.grid(alpha=.16, color=GRID, linewidth=.8)
    ax.set_axisbelow(True)
    for s in ax.spines.values(): s.set_color('#555b63')

def jagged(ax, x, y, color, label=None, marker='o', markevery=12, lw=2.0, ls='--'):
    # Existing data are untouched; dash + sparse markers create the requested "bergerigi" visual style.
    return ax.plot(x, y, color=color, alpha=.82, lw=lw, ls=ls,
                   marker=marker, markersize=3.2, markevery=markevery,
                   markerfacecolor='white', markeredgewidth=.8, label=label)

# ---------------- Figure 1: snake architecture ----------------
fig, ax = plt.subplots(figsize=(14, 8.2)); ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')
ax.text(.5,.965,'Stochastic Optimal-Stopping Architecture for ERCP Timing',ha='center',va='top',weight='bold',fontsize=18)
ax.text(.5,.905,'Literature-parameterized proof-of-concept: no patient level calibration',ha='center',va='top',fontsize=11,style='italic',color='#444')
# Snake: row 1 L->R, row 2 R->L, row 3 centered
boxes = [
 (.045,.68,.25,.145,'Routine Laboratory Inputs','Bilirubin | WBC\nCreatinine | Platelets | Albumin','#e6e1eb'),
 (.375,.68,.25,.145,'Continuous Dynamics','Inflammatory burden I(t)\nOrgan dysfunction O(t)\nStability reserve R(t)','#e5e9e6'),
 (.705,.68,.25,.145,'Latent Clinical Regimes','Stable → progressive\nOrgan dysfunction →\nintervention-favored','#eee8dd'),
 (.705,.39,.25,.145,'Optimal Stopping','Bellman recursion\nContinuation vs ERCP region','#e7e4ec'),
 (.375,.39,.25,.145,'Competing Hazards','Waiting/deterioration intensity\nProcedure-associated intensity','#eee3e8'),
 (.255,.105,.49,.15,'Robustness and Validation','LHS-PRCC | stochastic perturbation | solver checks\nFuture retrospective/prospective calibration','#e2e8eb'),
]
for x,y,w,h,t,s,c in boxes:
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.012,rounding_size=.018',fc=c,ec='#48505a',lw=1.45))
    ax.text(x+w/2,y+h*.68,t,ha='center',va='center',weight='bold',fontsize=11.5)
    ax.text(x+w/2,y+h*.32,s,ha='center',va='center',fontsize=9.3,linespacing=1.25)
arrow_kw=dict(arrowstyle='-|>',mutation_scale=16,lw=1.6,color='#59616b',connectionstyle='arc3,rad=0')
for a,b in [((.295,.752),(.375,.752)),((.625,.752),(.705,.752)),((.83,.68),(.83,.535)),((.705,.462),(.625,.462)),((.50,.39),(.50,.255))]:
    ax.add_patch(FancyArrowPatch(a,b,**arrow_kw))
save('Figure_1_Model_Architecture.png')

# ---------------- Figure 2 ----------------
fig, axs = plt.subplots(2,3,figsize=(13.8,8.3)); axs=axs.ravel()
series=[('bilirubin_mg_dl','Bilirubin','mg/dL',5),('wbc_10e3_ul','WBC','10³/µL',12),('creatinine_mg_dl','Creatinine','mg/dL',2),('platelets_10e3_ul','Platelets','10³/µL',100),('albumin_g_dl','Albumin','g/dL',2.8)]
for i,(ax,(c,t,u,thr)) in enumerate(zip(axs,series)):
    jagged(ax,ref.hour,ref[c],MUTED[i],marker=['o','s','D','^','v'][i],markevery=16,ls=(0,(4,2)))
    ax.axhline(thr,ls=(0,(3,2)),lw=1.15,color=THRESH,alpha=.75)
    ax.set_title(t); ax.set_xlabel('Hours from presentation'); ax.set_ylabel(u); style_ax(ax)
for c,l,col,m in [('inflammatory_burden','Inflammatory burden',MUTED[0],'o'),('organ_dysfunction','Organ dysfunction',MUTED[1],'s'),('stability_reserve','Stability reserve',MUTED[2],'^')]:
    jagged(axs[5],ref.hour,ref[c],col,l,m,15,1.9,(0,(4,2)))
axs[5].set_title('Latent Physiologic State'); axs[5].set_xlabel('Hours from presentation'); axs[5].set_ylim(0,1); axs[5].legend(frameon=False,fontsize=8); style_ax(axs[5])
fig.suptitle('Reference Laboratory and Physiologic Trajectory Without Drainage',fontsize=15,weight='bold',y=.99); fig.tight_layout(rect=[0,0,1,.965])
save('Figure_2_Reference_Trajectory.png')

# ---------------- Figure 3 ----------------
fig,axs=plt.subplots(1,2,figsize=(13.3,5.6))
jagged(axs[0],refpath.hour,refpath.wait_hazard,PURPLE,'Waiting/deterioration intensity','o',8,2.0,(0,(4,2)))
jagged(axs[0],refpath.hour,refpath.procedure_hazard,TAUPE,'Procedure-associated intensity','s',8,2.0,(0,(4,2)))
axs[0].axvline(ref_stop,ls=(0,(3,2)),lw=1.4,color=INK,label=f'Model stopping point {ref_stop:.0f} h')
axs[0].set_yscale('log'); axs[0].set_xlabel('Hours'); axs[0].set_ylabel('Theoretical intensity'); axs[0].set_title('Competing Hazards'); axs[0].legend(frameon=False,fontsize=8); style_ax(axs[0])
jagged(axs[1],refpath.hour,refpath.operate_now_cost,SLATE,'Drain now value','o',8,2.0,(0,(4,2)))
jagged(axs[1],refpath.hour,refpath.wait_one_hour_cost,MAUVE,'Wait one hour value','s',8,2.0,(0,(4,2)))
axs[1].axvline(ref_stop,ls=(0,(3,2)),lw=1.4,color=INK,label='First stopping point')
axs[1].fill_between(refpath.hour,0,1,where=(refpath.operate==0),alpha=.10,color='#8c8795',transform=axs[1].get_xaxis_transform(),label='Continuation region')
axs[1].fill_between(refpath.hour,0,1,where=(refpath.operate==1),alpha=.10,color='#a18d7b',transform=axs[1].get_xaxis_transform(),label='Drainage-favored region')
axs[1].set_xlabel('Hours'); axs[1].set_ylabel('Model value'); axs[1].set_title('Finite-Horizon Optimal-Stopping Path'); axs[1].legend(frameon=False,fontsize=8); style_ax(axs[1]); fig.tight_layout()
save('Figure_3_Competing_Hazards_and_Stopping_Path.png')

# ---------------- Figure 4 ----------------
piv=surf.pivot(index='creatinine_mg_dl',columns='bilirubin_mg_dl',values='stopping_hour')
fig,ax=plt.subplots(figsize=(10.8,7.3))
im=ax.imshow(piv.values,origin='lower',aspect='auto',extent=[piv.columns.min(),piv.columns.max(),piv.index.min(),piv.index.max()],cmap='magma_r',alpha=.92)
levels=np.unique(np.quantile(piv.values,[.1,.3,.5,.7,.9]).astype(int)); cs=ax.contour(piv.columns,piv.index,piv.values,levels=levels,colors='white',linewidths=1.05,alpha=.9); ax.clabel(cs,fmt=lambda x:f'{x:.0f} h',fontsize=8)
ax.axvline(5,ls=(0,(3,2)),color='black',lw=1.2); ax.axhline(2,ls=(0,(3,2)),color='black',lw=1.2)
ax.scatter([6],[1.35],marker='*',s=190,c='white',edgecolors='black',linewidths=1.1,zorder=5,label='Reference profile: bilirubin 6 mg/dL, creatinine 1.35 mg/dL')
ax.set_xlabel('Initial bilirubin (mg/dL)'); ax.set_ylabel('Initial creatinine (mg/dL)'); ax.set_title('Bilirubin–Creatinine Stopping Surface',weight='bold')
cbar=fig.colorbar(im,ax=ax,label='Model-derived stopping hour'); ax.legend(frameon=True,fontsize=8,loc='upper right'); style_ax(ax)
ax.text(.015,-.13,'Star = reference moderate profile used for the illustrative stopping-path analysis.',transform=ax.transAxes,fontsize=8.5,color='#4b4f55')
save('Figure_4_Bilirubin_Creatinine_Stopping_Surface.png')

# ---------------- Figure 5 ----------------
fig,ax=plt.subplots(figsize=(11.2,6.6)); y=np.arange(len(scdf)); bars=ax.barh(y,scdf.model_stopping_hour,color=DUST,alpha=.62,edgecolor='#655969',linewidth=.8)
ax.set_yticks(y,scdf.scenario); ax.invert_yaxis(); ax.axvspan(24,48,alpha=.10,color='#8f8796',label='24–48 h evidence window')
for i,v in enumerate(scdf.model_stopping_hour): ax.text(v+.5,i,f'{v:.0f} h',va='center',fontsize=9)
ax.set_xlabel('Model-derived stopping hour'); ax.set_title('Scenario-Specific Shift in Candidate ERCP Timing',weight='bold'); ax.legend(frameon=False); style_ax(ax)
save('Figure_5_Scenario_Stopping_Hours.png')

# ---------------- Figure 6 ----------------
sp=prdf.sort_values('PRCC_stopping_hour'); vals=sp.PRCC_stopping_hour
fig,ax=plt.subplots(figsize=(10.6,7.2)); cols=[NEG if x<0 else POS for x in vals]
ax.barh(sp.parameter,vals,color=cols,alpha=.68,edgecolor='#555b63',linewidth=.45); ax.axvline(0,color=INK,lw=1); ax.set_xlim(-1,1)
ax.set_xlabel('PRCC with stopping hour'); ax.set_title('Global Sensitivity Across 600 Latin-Hypercube Parameter Sets',weight='bold'); style_ax(ax)
ax.legend(handles=[Patch(facecolor=POS,alpha=.68,label='Positive PRCC'),Patch(facecolor=NEG,alpha=.68,label='Negative PRCC')],frameon=False,loc='lower right',fontsize=8)
save('Figure_6_Global_Sensitivity.png')

# ---------------- Figure 7 ----------------
fig,ax=plt.subplots(figsize=(10.7,6.2)); bins=np.arange(stdf.stopping_hour.min()-.5,stdf.stopping_hour.max()+1.5,1)
ax.hist(stdf.stopping_hour,bins=bins,color='#81758d',alpha=.58,edgecolor='#51485a',linewidth=1.0,linestyle=(0,(2,2)))
med=stdf.stopping_hour.median(); lo=stdf.stopping_hour.quantile(.025); hi=stdf.stopping_hour.quantile(.975)
ax.axvline(med,ls=(0,(4,2)),color=INK,lw=1.5,label=f'Median {med:.0f} h'); ax.axvspan(lo,hi,alpha=.09,color='#8b8293',label=f'2.5–97.5% {lo:.0f}–{hi:.0f} h')
ax.set_xlabel('Model-derived stopping hour'); ax.set_ylabel('Simulation count'); ax.set_title('Stochastic Perturbation of Routine Laboratory Inputs',weight='bold'); ax.legend(frameon=False); style_ax(ax)
save('Figure_7_Stochastic_Stopping_Distribution.png')

# ---------------- Figure 8 ----------------
fig,ax=plt.subplots(figsize=(10.7,6.3)); entries=[('p_stable','Stable','#756b86','o'),('p_progressive','Progressive','#8b7968','s'),('p_organ_dysfunction','Organ dysfunction','#697d7a','^'),('p_intervention_favored','Intervention-favored','#8b707b','D')]
for c,l,col,m in entries: jagged(ax,ref.hour,ref[c],col,l,m,15,2.0,(0,(4,2)))
ax.axvline(ref_stop,ls=(0,(3,2)),lw=1.2,color=INK); ax.set_ylim(0,1); ax.set_xlabel('Hours'); ax.set_ylabel('Probability'); ax.set_title('Latent Regime Occupancy Under Continued Observation',weight='bold'); ax.legend(frameon=False); style_ax(ax)
save('Figure_8_Regime_Occupancy.png')

# ---------------- Figure 9 unchanged in essence ----------------
fig,ax=plt.subplots(figsize=(8.6,6.6)); sca=ax.scatter(ref.inflammatory_burden,ref.organ_dysfunction,c=ref.hour,cmap='plasma',s=18); ax.plot(ref.inflammatory_burden,ref.organ_dysfunction,alpha=.45,color='#8fb5d1'); ax.set_xlabel('Inflammatory burden I(t)'); ax.set_ylabel('Organ dysfunction O(t)'); ax.set_title('Phase-Plane Trajectory Toward Deterioration',weight='bold'); fig.colorbar(sca,ax=ax,label='Hours'); style_ax(ax)
save('Figure_9_Phase_Plane.png')

# ---------------- Figure 10: informative lollipop spectrum ----------------
# Uses exactly the same eigenvalue CSV; only presentation is revised.
states=eqdf['state'].astype(str).tolist(); eigs=eigdf['real'].to_numpy(); n=min(len(states),len(eigs)); states=states[:n]; eigs=eigs[:n]
labels={'bilirubin':'Bilirubin','WBC':'WBC','creatinine':'Creatinine','platelets':'Platelets','albumin':'Albumin','inflammatory_burden':'Inflammatory burden','organ_dysfunction':'Organ dysfunction','stability_reserve':'Stability reserve'}
states=[labels.get(s,s.replace('_',' ').title()) for s in states]
order=np.argsort(eigs); eigs=eigs[order]; states=np.array(states)[order]; yy=np.arange(n)
fig,ax=plt.subplots(figsize=(10.2,6.2)); ax.axvspan(eigs.min()-.012,0,color='#81758d',alpha=.055); ax.axvline(0,color=INK,lw=1.25,label='Stability boundary Re(λ)=0')
for y,e in zip(yy,eigs): ax.hlines(y,e,0,color='#82768e',alpha=.55,lw=2); ax.scatter(e,y,s=80,color='#756881',alpha=.82,edgecolor='white',linewidth=.7,zorder=3); ax.text(e-.004,y,f'{e:.3f}',ha='right',va='center',fontsize=8,color='#4d4652')
ax.set_yticks(yy,states); ax.set_xlabel('Real part of eigenvalue, Re(λ)'); ax.set_title('Post-Drainage Linearized Eigenvalue Spectrum',weight='bold'); ax.set_xlim(eigs.min()-.035,.012); style_ax(ax); ax.grid(axis='y',alpha=.08); ax.legend(frameon=False,loc='lower right',fontsize=8)
ax.text(.02,.02,'All displayed real parts are < 0; the plot visualizes relative decay rates\nfrom the existing Post_ERCP_Eigenvalues.csv output.',transform=ax.transAxes,fontsize=8.5,color='#555')
save('Figure_10_Post_ERCP_Eigenvalues.png')

# ---------------- Better table renderer ----------------
def wrap_cell(v,width):
    if pd.isna(v): return ''
    s=str(v)
    return '\n'.join(textwrap.wrap(s,width=width,break_long_words=False,break_on_hyphens=False)) if len(s)>width else s

def tablefig(df,title,name,fontsize=8,widths=None,wrap_widths=None):
    d=df.copy().fillna('')
    if wrap_widths:
        for j,col in enumerate(d.columns):
            w=wrap_widths[j] if j<len(wrap_widths) else 25
            d[col]=d[col].map(lambda x:wrap_cell(x,w))
    rows=len(d); fig_h=max(2.7,1.15+rows*.48)
    fig,ax=plt.subplots(figsize=(14,fig_h)); ax.axis('off')
    # bbox deliberately anchors table close to title instead of vertically centering it.
    tb=ax.table(cellText=d.values,colLabels=d.columns,cellLoc='left',colLoc='center',bbox=[0.01,0.02,0.98,0.86],colWidths=widths)
    tb.auto_set_font_size(False); tb.set_fontsize(fontsize)
    for (r,c),cell in tb.get_celld().items():
        cell.set_edgecolor('#4c5055'); cell.set_linewidth(.75); cell.PAD=.08
        if r==0:
            cell.set_facecolor('#e9e7eb'); cell.set_text_props(weight='bold',ha='center',va='center')
        else:
            cell.set_facecolor('#fafafa' if r%2 else '#f2f1f3'); cell.set_text_props(va='center')
    ax.set_title(title,weight='bold',fontsize=14,pad=5,y=.94)
    save(name)

clinical=pd.DataFrame([['Bilirubin','5 mg/dL','TG18 Grade II'],['WBC','>12 or <4 ×10³/µL','TG18 Grade II'],['Creatinine','>2 mg/dL','TG18 Grade III renal dysfunction'],['Platelets','<100 ×10³/µL','TG18 Grade III hematologic dysfunction'],['Albumin','<0.7 × LLN','TG18 Grade II']],columns=['Variable','Literature anchor','Clinical role'])
tablefig(clinical,'Literature-Informed Laboratory Anchors','Figure_11_Table_Clinical_Anchors.png',8.5,widths=[.25,.30,.45],wrap_widths=[20,24,36])
tablefig(scdf.round(3),'Scenario Outputs From the Proof-of-Concept Model','Figure_12_Table_Scenario_Outputs.png',8.0,widths=[.36,.20,.23,.21],wrap_widths=[35,20,24,24])
tablefig(pars.head(16),'Parameter Provenance and Evidence Status','Figure_13_Table_Parameter_Provenance.png',6.2,widths=[.15,.10,.08,.10,.08,.08,.14,.19,.08],wrap_widths=[24,18,12,16,12,12,22,30,10])
tablefig(src[['source_id','year','parameter_or_anchor','reported_anchor_or_role']].head(17),'Literature Sources Used for Anchors and Timing Context','Figure_14_Table_Literature_Sources.png',7.0,widths=[.08,.08,.28,.56],wrap_widths=[8,8,34,68])

print(f'Finished. Revised figures saved to: {FIG}')
