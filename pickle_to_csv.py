import pickle as pkl
import pandas as pd
from tabulate import tabulate
#import tabulate for fancy printing
with open("Database_complete_plus_I261.pkl", "rb") as f:
    object = pkl.load(f)
df = pd.DataFrame(object)
r_df=df[df['uncertainty_AME20']!='nan']
r_df=r_df[r_df['Z']==47]
print (tabulate(r_df[['Z','A','uncertainty_AME20']]))
#choose certain columns if needed
#r_df=r_df[['Z','A','mass_excess_AME20','uncertainty_AME20','g2n_AME20']]
r_df.to_csv(r'Silver_isotopes_masses_and_Qvals.csv')
