import pandas as pd, numpy as np, os
from config import MADRID_DATA, TEXAS_DATA
W=os.path.dirname(os.path.abspath(__file__))
def madrid():
    f=os.path.join(W,'madrid.pkl')
    if os.path.exists(f): return pd.read_pickle(f)
    df=pd.read_csv(MADRID_DATA,sep='\t',usecols=[1,2,3,4,5,6,7],encoding='latin-1')
    df.columns=['St','Fecha','Year','Doy','HM','T','RH']
    n0=len(df); nT=df['T'].isna().sum(); nH=df['RH'].isna().sum(); nb=(df['T'].isna()&df['RH'].isna()).sum()
    print('Madrid raw rows',n0,'missing T',nT,'missing RH',nH,'both',nb)
    df=df.dropna(subset=['T','RH']).reset_index(drop=True)
    df['Month']=pd.to_datetime(df['Fecha'],format='%d/%m/%Y').dt.month
    df['Hour']=(df['HM']//100).astype(int); df['Slot']=df['Hour']*2+((df['HM']%100)>=30).astype(int)
    df.to_pickle(f); return df
def texas():
    f=os.path.join(W,'texas.pkl')
    if os.path.exists(f): return pd.read_pickle(f)
    df=pd.read_csv(TEXAS_DATA,sep='\t')
    df.columns=['Year','Month','Day','Hour','Minute','T','RH','St']
    print('Texas raw rows',len(df),'missing T',df['T'].isna().sum(),'missing RH',df['RH'].isna().sum())
    df=df.dropna(subset=['T','RH']).reset_index(drop=True)
    df.to_pickle(f); return df
if __name__=='__main__':
    m=madrid(); t=texas()
    print(m.groupby('St').agg(n=('T','size'),y0=('Year','min'),y1=('Year','max'),Tmax=('T','max'),RHmin=('RH','min'),RHmax=('RH','max')))
    print(t.groupby('St').agg(n=('T','size'),y0=('Year','min'),y1=('Year','max'),Tmax=('T','max'),RHmin=('RH','min'),RHmax=('RH','max')))
