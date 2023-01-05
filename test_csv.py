import Core.SelectBarcode as DataInput

RB = DataInput.ReadBarcode()
filedata = RB.UseCSV()
print(RB.user, RB.project)
exceldata = RB.SelectFromExcel()
print(exceldata)