import Core.CoreSystem as CS

RB = CS.ReadBarcode()
exceldata = RB.SelectFromExcel()
filedata = RB.UseCSV()
print(RB.user, RB.project)
print(RB.BarcodeList)
