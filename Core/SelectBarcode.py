import pandas as pd

def SelectFrom(*args):
    user_id, project_id, filepath = args
    df = pd.read_excel(filepath, engine='openpyxl')
    
def UseCSV(*args):
    f_csv = pd.read_csv(args)
