import pandas as pd
import sys

with pd.ExcelFile(sys.argv[1]) as xls:
    for sheet_name in xls.sheet_names:
        df = pd.read_excel(xls, sheet_name, header=None)
        df.to_csv(sheet_name+".csv", header=False, index=False)
