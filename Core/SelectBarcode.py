import pandas as pd
import csv
class ReadBarcode(object):

    def __init__(self):
        f = open('./Input/input.csv', 'r')
        fr = csv.reader(f)
        initdata = [line for line in fr]
        self.FilePath = ''
        self.user = initdata[0][0]
        self.project = initdata[0][1]
        self.BarcodeList = []
        f.close()

    def SelectFromExcel(self):
        df = pd.read_excel('./User/barcode.xlsx') #edit to pathlib later

    def UseCSV(self):
        f_csv = open('./Input/input.csv', 'r') #edit to pathlib later
        f_read = csv.reader(f_csv)
        csv_data = [line for line in f_read]
        self.BarcodeList = csv_data[1:]

        return csv_data
