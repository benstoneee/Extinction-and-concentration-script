#!/usr/bin/env python
# coding: utf-8

# In[70]:


import matplotlib.pyplot as plt
import pandas as pd
import glob
import os

#takes in filename and returns raw data in dataframe (df), used in reading_folder function
def read_data(file_name):
    df = pd.read_csv(file_name,sep="\t", header = 4)
    d = df.dropna(thresh = 1, axis = 1)
    d = d.dropna()
    d = d.drop(columns=['A', 'B'])
    d.reset_index(drop=True, inplace=True)
    return d

#reads data in entire folder, specified as filepath
def reading_folder(filepath):
    filenames = []
    filedata = pd.DataFrame()
    for name in glob.glob(filepath):
        base = os.path.basename(name)
        filenames.append(os.path.splitext(base)[0])
        data = read_data(name)
        filedata = filedata.append(data)
    return filenames,filedata

#labels dataframe by wavelength and absorabnces max, returns df
def label_data(d):
    absdata = d[d < 10].dropna()
    wavdata = d[d > 10].dropna()
    abs_dict = absdata.to_dict('list')
    col_max = max(abs_dict, key=abs_dict.get)
    d = pd.DataFrame({'absorbances': list(absdata[col_max]),'wavelengths': list(wavdata[col_max])})
    return d

#User input of EC or concentration, glassware/pipettes used, adds new columns for results and returns df
def user_protocol_add_columns(d):
    lst = []
    d = d.assign(choice = input("Calulating Extinction Coefficient(ec) or Concentration(conc)? enter ec or conc"))
    if d['choice'].any() == "ec":
        n = int(input("Enter number of samples"))
        for i in range(0, n): 
            ele = float((input("Sample weight? (mg)")))               # consilidate to one assign each
            lst.append(ele)
        d = d.assign(sampleweight = lst,
            flaskvolume = int(input("What volume was your samples diluted into? (mL)")),
            mw = float(input("What is the Molecular weight of your compound? (g/mol)")),
            dilution = int(input("What was your dilution factor?")))
        d = d.assign(concentration = (((d['sampleweight']/d['flaskvolume'])/d['mw'])/d['dilution']))
        d = d.assign(ec = (d['absorbances']/d['concentration']))
    if d['choice'].any()  == "conc":
        d = d.assign(sample_vol = int(input("What was the volume of sample that was diluted? (mL)")),
            flask_vol = int(input("What volume was your samples diluted into? (mL)")),
            ec = int(input("What is the Extinction Coefficient of your sample")))
        d = d.assign(concentration = ((d['absorbances']/d['ec'])*(d['flask_vol']/d['sample_vol'])))
    return d

#Makes new columns for results in error-specified df
def error_calc(d):
    d = d.assign(eterm = (d['errors']/d['values'])**2)
    d = d.assign(percenterror = sum(d['eterm'])**0.5)
    return d['percenterror'][0]

#Manual input of what values/tolerances used and adds new column of error and error in result
def build_error(data):
    #here is dicts for the tolerances on each glassware,hoping to not manually filter this
    error_nflasks = {'5':0.02,'10':0.02,'25':0.03,'50':0.05,'100':0.08,'200':0.1,'250':0.12}
    error_wflasks = {'5':0.08,'10':0.08,'25':0.08,'50':0.08,'100':0.1,'200':0.2,'250':0.2}
    error_Vpip = {'0.5':0.006,'1':0.006,'2':0.006,'3':0.01,'5':0.01,'10':0.02}
    error_Mpip = {'5':0.3,'10':0.5,'20':0.8,'50':1.5,'100':2,'200':3,'1000':10,'5000':50,'10000':100}

    error_lst = []
    indices = pd.Series(['weight','v1','v2','v3','A'])
    d1 = pd.DataFrame({'values':[11.5,25,0.5,25,0.654],'errors':[0.2,0.03,0.006,0.03,0.004]}, index = indices)
    d2 = pd.DataFrame({'values':[14.3,25,0.5,25,0.654],'errors':[0.2,0.03,0.006,0.03,0.004]}, index = indices)
    d3 = pd.DataFrame({'values':[13.6,25,0.5,25,0.654],'errors':[0.2,0.03,0.006,0.03,0.004]}, index = indices)
    d4 = pd.DataFrame({'values':[10.7,25,0.5,25,0.699],'errors':[0.2,0.03,0.006,0.03,0.004]}, index = indices)
    d5 = pd.DataFrame({'values':[9.7,25,0.5,25,0.654],'errors':[0.2,0.03,0.006,0.03,0.004]}, index = indices)
    d6 = pd.DataFrame({'values':[10.0,25,0.5,25,0.654],'errors':[0.2,0.03,0.006,0.03,0.004]}, index = indices)
    d7 = pd.DataFrame({'values':[12.7,25,0.5,25,0.654],'errors':[0.2,0.03,0.006,0.03,0.004]}, index = indices)
    d8 = pd.DataFrame({'values':[10.0,25,0.5,25,0.699],'errors':[0.2,0.03,0.006,0.03,0.004]}, index = indices)
    dlst = [d1,d2,d3,d4,d5,d6,d7,d8]
    for i in range(len(data)):
        error = error_calc(dlst[i])
        error_lst.append(error)
    data = data.assign(percenterror = error_lst)
    if data['choice'].any() == "ec":
        data = data.assign(error_in_ec = data['percenterror']*data['ec'])
    if data['choice'].any()  == "conc":
        data = data.assign(error_in_conc = data['percenterror']*data['conc'])
    return data

#prints data to csv and results/calculation to txt file
def results_printing(d):
    f = open('-','w')
    d.to_csv('-')
    if d['choice'].any() == "ec":
        for i in range(len(d)):
            conc_str = ('C%d=(%.1f mg/%d ml)*(mol/%.2f g)*(1/%d) = %f M\n' %(i+1,d['sampleweight'][i],d['flaskvolume'][i],d['mw'][i],d['dilution'][i],d['concentration'][i]))
            ec_str = ('E%d=%.3f/%f = %.2f\n\n'% (i+1,d['absorbances'][i],d['concentration'][i],d['ec'][i]))
            f.write(conc_str)
            f.write(ec_str)
        avg_ec = d.loc[:,'ec'].mean()
        plusminus_ec = d.loc[:,'error_in_ec'].mean()
        ec_rsd = (d.loc[:,'ec'].std()/avg_ec)*100
        result_str =('Average EC = %.2f +/- %.2f @ %.1f nm, RSD  = %.2f %%' % (avg_ec,plusminus_ec,d['wavelengths'][i],ec_rsd))
        f.write(result_str) 
    if d['choice'].any() == "conc":
        for i in range(len(d)):
            conc_str('C%d=(%.1f/%f)*(%d mL/%f mL) = %f M' 
                     %(i+1,d['absorbances'][i],d['ec'][i],d['flask_col'][i],d['sample_vol'][i],d['concentration'][i]))
            f.write(conc_str)
        avg_conc = d.loc[:,'concentration'].mean()
        rsd_conc = (d.loc[:,'concentration'].std()/avg_conc)*100
        result_str = ('Average Concentration = %f M, %RSD = %.2f %%' %(avg_conc,rsd_conc))
        f.write(result_str)

filepath = '-'
names_files, listed_data = reading_folder(filepath)
data1 = label_data(listed_data)
#data2 = user_protocol_add_columns(data1)
#data2 = build_error(data2)
results_printing(data2)


# ## 

# In[12]:


data2


# In[ ]:





# In[ ]:




