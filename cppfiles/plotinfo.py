import pandas as pd
import matplotlib.pyplot as plt

file_name = ["./Outfile/Star_WENO_Density_NM.dat", "./Outfile/Star_WENO_Density_DM.dat"]
col_name = [["step","rhoNM"], ["step","rhoDM"]]

## Read files
dfNM = pd.read_csv(file_name[0],delim_whitespace=True, names=col_name[0], header=1, engine='python')
dfDM = pd.read_csv(file_name[1],delim_whitespace=True, names=col_name[1], header=1, engine='python')
# Concatenate data according to stepnumber
df = pd.merge(dfDM, dfNM, on='step')

print(df)