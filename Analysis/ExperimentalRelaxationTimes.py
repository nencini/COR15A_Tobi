import numpy as np
import gc
import matplotlib.pyplot as plt
import math
from datetime import date
today = date.today()

gc.collect()



def t1_t2_relaxations(data_from_DynamicCenter,field,coeff,outputP,author,info,increments):
    with open(data_from_DynamicCenter,"r") as f:
        read_data=False
        peak_list=False
        peak=0
        results=False
        for line in f:
            if "integrals" in line and len(line.split())==2:
                read_data=True
            elif "results" in line and len(line.split())==2:
                results=True
                peak_dictionary={}
                peak_dictionary["peaks"]={}


            if "integral errors" in line:
                read_data=False

            if read_data:
                if "Mixing time [s]:" in line:
                    mixing_times=np.array(line.split()[3:])
                    mixing_times = mixing_times.astype('float64')
                elif "Peak name" in line or "SECTION:" in line:
                    pass
                else:
                    if peak==0:
                        intensities=np.array(line.split())
                        intensities = [intensities.astype('float64')]
                        peaks_intensities=intensities
                        line_length=len(line.split())
                    elif line_length==len(line.split()):
                        intensities=np.array(line.split())
                        intensities = [intensities.astype('float64')]
                        peaks_intensities=np.append(peaks_intensities,intensities,axis=0)
                    peak+=1
            elif results:
                if "Peak name" in line or "SECTION:" in line:
                    pass
                else:
                    peak_dictionary["peaks"][line.split()[0]]={}
                    peak_dictionary["peaks"][line.split()[0]]["ppm"]=[line.split()[1],line.split()[2]]

    for peak in range(0,peaks_intensities.shape[0]):
        peak_dictionary["peaks"][str(int(peaks_intensities[peak,0]))]["experiment"]=peaks_intensities[peak,1:]
            
        for_LS=np.array([[0, 0]])
        for i in range(1,peaks_intensities.shape[1]):
            if i<peaks_intensities.shape[1]-1:
                if peaks_intensities[peak,i]<peaks_intensities[peak,i+1]:
                    pass
                elif peaks_intensities[peak,i]==0:
                    break
                else:
                    for_LS=np.append(for_LS,[[mixing_times[i-1],np.log(peaks_intensities[peak,i])]],axis=0)
                       
            else:
                if peaks_intensities[peak,i]==0:
                    break
                else:
                    for_LS=np.append(for_LS,[[mixing_times[i-1],np.log(peaks_intensities[peak,i])]],axis=0)
        for_LS=np.delete(for_LS,0,0)
        curve_fit_coef = np.polyfit(for_LS[:,0], for_LS[:,1], 1, full=True)
            
        
    
        peak_dictionary["peaks"][str(int(peaks_intensities[peak,0]))]["fit"]=curve_fit_coef
        peak_dictionary["peaks"][str(int(peaks_intensities[peak,0]))]["points_used"]=len(for_LS)
    peak_dictionary["mixing_times"]=mixing_times
    peak_dictionary["INFO"]={}
    peak_dictionary["INFO"]["ANALYZED"]=str(today)
    peak_dictionary["INFO"]["AUTHOR"]=author
    peak_dictionary["INFO"]["FIELD"]=field
    peak_dictionary["INFO"]["COEFFICIENTS"]=coeff
    peak_dictionary["INFO"]["OUTPUT_POINTS"]=outputP
    peak_dictionary["INFO"]["INFO"]=info
    peak_dictionary["INFO"]["INCREMENTS"]=increments
    
        
    return peak_dictionary
            

def print_results(*files):
    for file in files:
        relax_data=t1_t2_relaxations(file)
        
        curve_fit_x=np.arange(relax_data[1][0],relax_data[1][-1],0.01)
        for peak in relax_data[0]:
            curve_fit_y=np.exp(curve_fit_x*  relax_data[0][peak]["fit"][0][0])*np.exp(relax_data[0][peak]["fit"][0][1])
                
                
def compare_spectra(*files):
    most_peaks=0
    for i,file in enumerate(files):
        if len(file['peaks'])>most_peaks:
            most_peaks=len(file['peaks'])
            position=i
            
    reference=files[position]
    other_files=[]
    
    for i,file in enumerate(files):
        if not i==position:
            other_files.append(file)
    
    files=other_files
    
    relaxation_times={}
    for ref_peak in reference['peaks']:
        relaxation_times[ref_peak]={}
        relaxation_times[ref_peak]["REFERENCE"]={}
        relaxation_times[ref_peak]["REFERENCE"]["ppm"]=reference['peaks'][ref_peak]["ppm"]
        relaxation_times[ref_peak]["REFERENCE"]["T1"]=  -1/reference['peaks'][ref_peak]["fit"][0][0]
        relaxation_times[ref_peak]["REFERENCE"]["COEFFICIENTS"]=  reference["INFO"]["COEFFICIENTS"]
        relaxation_times[ref_peak]["REFERENCE"]["OUTPUT_POINTS"]=  reference["INFO"]["OUTPUT_POINTS"]
        #relaxation_times[ref_peak]["REFERENCE"]["error"]=  -1/reference[ref_peak]["fit"][2][0]
        
    for k,file in enumerate(files):
        distances=np.zeros([len(reference['peaks']),len(file['peaks'])])
        for i,ref_peak in enumerate(reference['peaks']):
            #print(float(reference[ref_peak]["ppm"][0]))
            for j,file_peak in enumerate(file['peaks']):
                distances[i,j]=(math.dist([float(reference['peaks'][ref_peak]["ppm"][0]), 
                                           float(reference['peaks'][ref_peak]["ppm"][1])],
                                           [float(file['peaks'][file_peak]["ppm"][0]),
                                            float(file['peaks'][file_peak]["ppm"][1])
                                           ])) 
        
        

        reference_keys=list(reference['peaks'])
        for i,file_peak in enumerate(file['peaks']):
            key=np.where(distances[:,i]==min(distances[:,i]))[0][0]
            if np.where(distances[key,:]==min(distances[key,:]))[0][0]==i:
                relaxation_times[reference_keys[key]]["file"+str(k)]={}
                relaxation_times[reference_keys[key]]["file"+str(k)]["ppm"]=file['peaks'][file_peak]["ppm"]
                relaxation_times[reference_keys[key]]["file"+str(k)]["T1"]=-1/file['peaks'][file_peak]["fit"][0][0]
                relaxation_times[reference_keys[key]]["file"+str(k)]["COEFFICIENTS"]=file["INFO"]["COEFFICIENTS"]
                relaxation_times[reference_keys[key]]["file"+str(k)]["OUTPUT_POINTS"]=file["INFO"]["OUTPUT_POINTS"]

             

           
    return relaxation_times
    
    
def plot_data(comparison,coeffs,outputs):
    keys=["REFERENCE"]
    for i in range(len(coeffs)*len(outputs)):
        keys.append("file"+str(i))

    for file in keys:
        x_axis=[]
        y_axis=[]
        for peak in comparison:
            if file in comparison[peak]:
                x_axis.append(int(peak))
                y_axis.append(comparison[peak][file]["T1"])
                if int(comparison[peak][file]["COEFFICIENTS"])==0 and int(comparison[peak][file]["OUTPUT_POINTS"])==0:
                    refe=file
                else:
                    refe=-20
        if refe==file:
            plt.plot(x_axis, y_axis,'o',markersize=15)
        else:
            plt.plot(x_axis, y_axis,'o',markersize=5)

    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)
    plt.xlabel('Peak')
    plt.ylabel('T1')
    plt.show()
    
