import Core.CoreSystem as CS
import pandas as pd
RB = CS.ReadBarcode()
exceldata = RB.SelectFromExcel()
filedata = RB.UseCSV()
print(filedata.drop_duplicates(['Gene name']))