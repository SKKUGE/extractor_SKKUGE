import Core.SelectBarcode as DataInput

RB = DataInput.ReadBarcode()
filedata = RB.UseCSV()
print(RB.user, RB.project)
for line in RB.BarcodeList:
    print(line)