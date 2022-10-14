"""Try if help works"""

import numpy as np
import gc
import matplotlib.pyplot as plt
import math
import yaml
from datetime import date
today = date.today()

gc.collect()



def t1_t2_relaxations(data_from_DynamicCenter,field,coeff,outputP,author,info,increments,skip_points=0,correct=(0,0)):
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
                elif "Time [s]:" in line:
                    mixing_times=np.array(line.split()[2:])
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
                    peak_dictionary["peaks"][line.split()[0]]["ppm"]=[float(line.split()[1])+correct[0],float(line.split()[2])+correct[1]]
    

    for peak in range(0,peaks_intensities.shape[0]):
        peak_dictionary["peaks"][str(int(peaks_intensities[peak,0]))]["experiment"]=peaks_intensities[peak,1:].tolist()
        for_LS=np.array([[0, 0]])
        for i in range(skip_points+1,peaks_intensities.shape[1]):
            if i<peaks_intensities.shape[1]-1:
                if peaks_intensities[peak,i]==0:
                    break
                elif peaks_intensities[peak,i]<peaks_intensities[peak,i+1]:
                    pass
                else:
                    for_LS=np.append(for_LS,[[mixing_times[i-1],np.log(peaks_intensities[peak,i])]],axis=0)
                       
            else:
                if peaks_intensities[peak,i]==0:
                    break
                else:
                    for_LS=np.append(for_LS,[[mixing_times[i-1],np.log(peaks_intensities[peak,i])]],axis=0)
        for_LS=np.delete(for_LS,0,0)
        
        if len(for_LS)>2:
            curve_fit_coef = np.polyfit(for_LS[:,0], for_LS[:,1], 1, full=True)
    
            peak_dictionary["peaks"][str(int(peaks_intensities[peak,0]))]["fit"]=[curve_fit_coef[0].tolist(),curve_fit_coef[1].tolist()]
            for_LS[:,1]=np.exp(for_LS[:,1])
            peak_dictionary["peaks"][str(int(peaks_intensities[peak,0]))]["points_used"]=[for_LS[:,0].tolist(),for_LS[:,1].tolist()]
        else:
            peak_dictionary["peaks"][str(int(peaks_intensities[peak,0]))]["fit"]=None
    
    peak_dictionary["mixing_times"]=mixing_times.tolist()
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
                
                
def compare_spectra(Nfactor,maxDist,*files):
    """Function groups peaks from different experiments.
       Developed for testing different processing parameters for T1 and T2 parameters.
       
       Takes the imput file with the most peaks as the reference.
       Groups peaks based on minimal Eucliac distance with Nfactor weight
       parameter for N dimention.
       
       In the future, the compare_spectra method should be universal for any comparison, including hetNOE
       
       Input: Nfactor, maxDist, *files
       
              Nfactor  - used to weight the N dimention with respect to H dimention.
                         Values in literature range between 0.1 and 0.45
              maxDist  - maximal distance for peaks to be defined as overlapping 
              *files   - peak_dictionary files from t1_t2_relaxations() method to be grouped
              
       
       Output: relaxation_times
       
               relaxation_times  - dictionary of all the reference peaks and the corresponding peaks
                                   that were found from other files
                                   
       #####################################################################################################
       #####################################################################################################
        
       """
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
        if not reference['peaks'][ref_peak]["fit"]==None:
            relaxation_times[ref_peak]={}
            relaxation_times[ref_peak]["REFERENCE"]={}
            relaxation_times[ref_peak]["REFERENCE"]["ppm"]=reference['peaks'][ref_peak]["ppm"]
            relaxation_times[ref_peak]["REFERENCE"]["T1"]=  float(-1/reference['peaks'][ref_peak]["fit"][0][0])
            relaxation_times[ref_peak]["REFERENCE"]["COEFFICIENTS"]=  reference["INFO"]["COEFFICIENTS"]
            relaxation_times[ref_peak]["REFERENCE"]["OUTPUT_POINTS"]=  reference["INFO"]["OUTPUT_POINTS"]
            relaxation_times[ref_peak]["REFERENCE"]["error"]=  np.sqrt(float(reference['peaks'][ref_peak]["fit"][1][0])/(int(len(reference['peaks'][ref_peak]["points_used"][0]))-1))
            relaxation_times[ref_peak]["REFERENCE"]["INFO"]=reference["INFO"]["INFO"]
            
    for k,file in enumerate(files):
        distances=np.zeros([len(reference['peaks']),len(file['peaks'])])
        for i,ref_peak in enumerate(reference['peaks']):
            #print(float(reference[ref_peak]["ppm"][0]))
            for j,file_peak in enumerate(file['peaks']):
                distances[i,j]=(math.dist([float(reference['peaks'][ref_peak]["ppm"][0]), 
                                           Nfactor*float(reference['peaks'][ref_peak]["ppm"][1])],
                                           [float(file['peaks'][file_peak]["ppm"][0]),
                                            Nfactor*float(file['peaks'][file_peak]["ppm"][1])
                                           ])) 
        
        

        reference_keys=list(reference['peaks'])
        for i,file_peak in enumerate(file['peaks']):
            key=np.where(distances[:,i]==min(distances[:,i]))[0][0]
            if np.where(distances[key,:]==min(distances[key,:]))[0][0]==i and np.sqrt(min(distances[key,:]))<maxDist:
                if not file['peaks'][file_peak]["fit"]==None and reference_keys[key] in relaxation_times:
                    relaxation_times[reference_keys[key]]["file"+str(k)]={}
                    relaxation_times[reference_keys[key]]["file"+str(k)]["ppm"]=file['peaks'][file_peak]["ppm"]
                    relaxation_times[reference_keys[key]]["file"+str(k)]["T1"]=float(-1/file['peaks'][file_peak]["fit"][0][0])
                    relaxation_times[reference_keys[key]]["file"+str(k)]["COEFFICIENTS"]=file["INFO"]["COEFFICIENTS"]
                    relaxation_times[reference_keys[key]]["file"+str(k)]["OUTPUT_POINTS"]=file["INFO"]["OUTPUT_POINTS"]
                    relaxation_times[reference_keys[key]]["file"+str(k)]["PEAK_NAME"]=file_peak
                    relaxation_times[reference_keys[key]]["file"+str(k)]["error"]=np.sqrt(float(file['peaks'][file_peak]["fit"][1][0])/(int(len(file['peaks'][file_peak]["points_used"][0]))-1))
                    relaxation_times[reference_keys[key]]["file"+str(k)]["INFO"]=file["INFO"]["INFO"]
             

           
    return relaxation_times
    
    
def plot_data(comparison,coeffs,outputs,highlight):
    keys=["REFERENCE"]
    for i in range(len(coeffs)*len(outputs)):
        keys.append("file"+str(i))

    for file in keys:
        x_axis=[]
        y_axis=[]
        y_error=[]
        for peak in comparison:
            if file in comparison[peak]:
                x_axis.append(int(peak))
                y_axis.append(comparison[peak][file]["T1"])
                y_error.append(comparison[peak][file]["error"])
                if int(comparison[peak][file]["COEFFICIENTS"])==highlight[0] and int(comparison[peak][file]["OUTPUT_POINTS"])==highlight[1]:
                    refe=file
                else:
                    refe=None
        if refe==file:
            plt.errorbar(x_axis, y_axis,yerr=y_error,fmt='o',markersize=15)
        else:
            plt.errorbar(x_axis, y_axis,yerr=y_error,fmt='o',markersize=5)

    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)
    plt.xlabel('Reference Peak')
    plt.ylabel('T1')
    plt.show()



def process_hetNoe(experiment,field,Nfactor,maxDist):
    """Function reads 2 files with peaks, 
       selects overlaping peaks 
       and calculates hetNOE values.
       
       Input: spectra1, spectra2 without NOE effect, Nfactor
       
              spectra1 - with NOE effect
              spectra2 - without NOE effect
              Nfactor  - used to weight the N dimention with respect to H dimention.
                         Values in literature range between 0.1 and 0.45
              maxDist  - maximal distance for peaks to be defined as overlapping 
              image    - hetNOE spectra image to check if peaks selected (optional)
              imageDim - ppm range of the figure
                         (H ppm left, H ppm right, N ppm bottom, N ppm top)
                         (optional, but mandatory if picture provided)
                         
       Output: hetNoe, peaks1, peaks2
       
               hetNoe - dictionary of peaks that were defined as overlapping 
               peaks1 - dictionary of all peaks from NOE spectra
               peaks2 - dictionary of all peaks from spectra without NOE effect        
                                           
       #####################################################################################################
       #####################################################################################################
                   
                         """
                         
    with open("/home/nenciric/Documents/git/COR15A_Tobi/ExperimentalData/experiments.yaml") as yaml_file:
        experiments = yaml.load(yaml_file, Loader=yaml.FullLoader)
      
    relax="hetNOE"  
    path=experiments[field][experiment][relax]["PATH"]

    info=experiments[field][experiment][relax]["INFO"]
    files=experiments[field][experiment][relax]["FILES"]
    coeffs=experiments[field][experiment][relax]["COEFF"]
    outputs=experiments[field][experiment][relax]["OUTPUTS"]
    increments=experiments[field][experiment][relax]["INCREMENTS"]
    correct=experiments[field][experiment][relax]["CORRECT"]
    image=experiments[field][experiment][relax]["IMAGE"]
    imageDim=experiments[field][experiment][relax]["IMAGEDIM"]

    file1=path+files[0]
    file2=path+files[1]
    image=path+image                     
    
    with open(file1,"r") as f:
        minus=False
        peak_dictionary1={}
        peak_dictionary1["peaks"]={}
        for line in f:
            if "Minus" in line:
                minus=True
            elif "Plus" in line:
                minus=False
            elif "Peak" in line: 
                pass
            else:
                if minus:
                    peak_dictionary1["peaks"]["-"+line.split()[0]]={}
                    peak_dictionary1["peaks"]["-"+line.split()[0]]["ppm"]=[line.split()[1],line.split()[2]]
                
                    peak_dictionary1["peaks"]["-"+line.split()[0]]["intensity"]=-float(line.split()[5])
                else:
                    peak_dictionary1["peaks"][line.split()[0]]={}
                    peak_dictionary1["peaks"][line.split()[0]]["ppm"]=[line.split()[1],line.split()[2]]
                    peak_dictionary1["peaks"][line.split()[0]]["intensity"]=float(line.split()[5])
                    
                    
                    
    with open(file2,"r") as f:
        minus=False
        peak_dictionary2={}
        peak_dictionary2["peaks"]={}
        for line in f:
            if "Minus" in line:
                minus=True
            elif "Plus" in line:
                minus=False
            elif "Peak" in line: 
                pass
            else:
                if minus:
                    peak_dictionary2["peaks"]["-"+line.split()[0]]={}
                    peak_dictionary2["peaks"]["-"+line.split()[0]]["ppm"]=[line.split()[1],line.split()[2]]
                
                    peak_dictionary2["peaks"]["-"+line.split()[0]]["intensity"]=-float(line.split()[5])
                else:
                    peak_dictionary2["peaks"][line.split()[0]]={}
                    peak_dictionary2["peaks"][line.split()[0]]["ppm"]=[line.split()[1],line.split()[2]]
                    peak_dictionary2["peaks"][line.split()[0]]["intensity"]=float(line.split()[5])
                    
                    
    
    if len(peak_dictionary2['peaks'])>len(peak_dictionary1['peaks']):    
        reference=peak_dictionary2
        refName="noNOE"
        filName="NOE"
        files=[peak_dictionary1]
    else:
        reference=peak_dictionary1
        refName="NOE"
        filName="noNOE"
        files=[peak_dictionary2]
    
    
    
    relaxation_times={}
    for ref_peak in reference['peaks']:
        relaxation_times[ref_peak]={}
        relaxation_times[ref_peak][refName]={}
        relaxation_times[ref_peak][refName]["ppm"]=reference['peaks'][ref_peak]["ppm"]
        relaxation_times[ref_peak][refName]["intensity"]=  reference['peaks'][ref_peak]["intensity"]
        
        
    for k,file in enumerate(files):
        distances=np.zeros([len(reference['peaks']),len(file['peaks'])])
        for i,ref_peak in enumerate(reference['peaks']):
            #print(float(reference[ref_peak]["ppm"][0]))
            for j,file_peak in enumerate(file['peaks']):
                distances[i,j]=(math.dist([float(reference['peaks'][ref_peak]["ppm"][0]), 
                                           Nfactor*float(reference['peaks'][ref_peak]["ppm"][1])],
                                           [float(file['peaks'][file_peak]["ppm"][0]),
                                            Nfactor*float(file['peaks'][file_peak]["ppm"][1])
                                           ])) 
                

        reference_keys=list(reference['peaks'])
        
        for i,file_peak in enumerate(file['peaks']):
            key=np.where(distances[:,i]==min(distances[:,i]))[0][0]
            
            if np.where(distances[key,:]==min(distances[key,:]))[0][0]==i and np.sqrt(min(distances[key,:]))<maxDist:
                relaxation_times[reference_keys[key]][filName]={}
                relaxation_times[reference_keys[key]][filName]["ppm"]=file['peaks'][file_peak]["ppm"]
                relaxation_times[reference_keys[key]][filName]["intensity"]=file['peaks'][file_peak]["intensity"]

             
    hetNoe={}
    hetNoe["peaks"]={}
    for peak in relaxation_times:
        if "noNOE" in relaxation_times[peak] and "NOE" in relaxation_times[peak]:
            hetNoe["peaks"][peak]={}
            hetNoe["peaks"][peak]["NOE"]=relaxation_times[peak]["NOE"]
            hetNoe["peaks"][peak]["noNOE"]=relaxation_times[peak]["noNOE"]
            hetNoe["peaks"][peak]["hetNOE"]=float(relaxation_times[peak]["noNOE"]["intensity"])/float(relaxation_times[peak]["NOE"]["intensity"])
    
    hetNoe["INFO"]=info
    
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)
    graph, (plot1) = plt.subplots(1, 1, figsize=(18.5, 10.5))
    #plot1.xticks(fontsize=16)
    #plot1.yticks(fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.rcParams["figure.figsize"] = [18.5, 10.50]
    plot1.set_title("Peaks for hetNOE relaxation",fontsize=25)
    plot1.set_xlabel("H [ppm]",fontsize=20)
    plot1.set_ylabel("N [ppm]",fontsize=20)
    
    plot1.invert_xaxis()
    plot1.invert_yaxis()
    
    #plt.xlim(8.7, 6.75)
    #plt.ylim(132, 105)

    for i,peak in enumerate(peak_dictionary1["peaks"]):
        plot1.plot(float(peak_dictionary1["peaks"][peak]["ppm"][0]), float(peak_dictionary1["peaks"][peak]["ppm"][1]),'o',markersize=17,color="red")

    for i,peak in enumerate(peak_dictionary2["peaks"]):
        plot1.plot(float(peak_dictionary2["peaks"][peak]["ppm"][0]), float(peak_dictionary2["peaks"][peak]["ppm"][1]),'o',markersize=12,color="blue")

    for peak in hetNoe["peaks"]:
        plot1.plot([float(hetNoe["peaks"][peak]["NOE"]["ppm"][0]),float(hetNoe["peaks"][peak]["noNOE"]["ppm"][0])], [float(hetNoe["peaks"][peak]["NOE"]["ppm"][1]),float(hetNoe["peaks"][peak]["noNOE"]["ppm"][1])],'-',markersize=8,color="green")
        plot1.plot(float(hetNoe["peaks"][peak]["NOE"]["ppm"][0]), float(hetNoe["peaks"][peak]["NOE"]["ppm"][1]),'o',markersize=6,color="yellow")

    plt.show()
    
    if not image==None:
        fig = plt.gcf()
        fig.set_size_inches(18.5, 10.5)
        plt.rcParams["figure.figsize"] = [18.50, 10.50]
        plt.rcParams["figure.autolayout"] = True
        graph, (plot1) = plt.subplots(1, 1, figsize=(18.5, 10.5))

        img = plt.imread(image)
        factor=(imageDim[2]-imageDim[3])/10.5/((imageDim[0]-imageDim[1])/18.5)
        plot1.imshow(img, extent=[factor*imageDim[0], factor*imageDim[1], imageDim[2], imageDim[3]])
        plot1.set_xlim(factor*imageDim[0], factor*imageDim[1])
        plot1.set_ylim(imageDim[2], imageDim[3])
    
        for peak in hetNoe["peaks"]:
            plot1.plot(float(hetNoe["peaks"][peak]["NOE"]["ppm"][0])*factor, float(hetNoe["peaks"][peak]["NOE"]["ppm"][1]),'o',markersize=6,color="yellow")
        
        plt.show()
    
    return hetNoe, peak_dictionary1, peak_dictionary2
                    

def plot_fit_of_peaks(highlight,comparison,*results):
    for result in results:
        if result["INFO"]["COEFFICIENTS"]==highlight[0] and result["INFO"]["OUTPUT_POINTS"]==highlight[1]:
            fig = plt.gcf()
            height=int(math.ceil(len(result["peaks"])/5))*4
            fig.set_size_inches(20, height)
            for i,peak in enumerate(result["peaks"]):
                #print("Peak",peak)
                ax = fig.add_subplot(int(math.ceil(len(result["peaks"])/5)),5,i+1)
                
                x_fit_data=np.linspace(result['mixing_times'][0],result['mixing_times'][-1],1000)
                ax.plot(result['mixing_times'],result["peaks"][peak]['experiment'],"o")
                if not result["peaks"][peak]['fit']==None:
                    y_fit_data=np.exp(x_fit_data*result["peaks"][peak]['fit'][0][0])*np.exp(result["peaks"][peak]['fit'][0][1])
                    ref_peak=convert_peak_number_to_reference(peak,highlight,comparison)
                    ax.plot(x_fit_data,y_fit_data,"-")
                if ref_peak==None:
                   ref_peak="None"
                ax.set_title("Peak "+peak+ ", ref "+ref_peak)
                
                
            plt.show()
              
              
              
def add_exp_data(field,name,exp,author,path,prefix,info,increments,correct,coeffs,outputs):
    with open("/home/nenciric/Documents/git/COR15A_Tobi/ExperimentalData/experiments.yaml") as yaml_file:
        experiments = yaml.load(yaml_file, Loader=yaml.FullLoader)
    
    if name not in experiments[field]:
        experiments[field][name]={}
    experiments[field][name][exp]={}
    experiments[field][name][exp]["AUTHOR"]=author
    experiments[field][name][exp]["PATH"]=path
    experiments[field][name][exp]["PREFIX"]=prefix
    experiments[field][name][exp]["INFO"]=info
    experiments[field][name][exp]["INCREMENTS"]=increments
    experiments[field][name][exp]["CORRECT"]=correct
    experiments[field][name][exp]["COEFF"]=coeffs
    experiments[field][name][exp]["OUTPUTS"]=outputs

    with open("/home/nenciric/Documents/git/COR15A_Tobi/ExperimentalData/experiments.yaml", 'w') as f:
        yaml.dump(experiments,f, sort_keys=True)
        
    print("600 Mhz experiments: \n")
    for key in experiments["600 MHz"]:
        print(key,end=" - ")
        for ex in experiments["600 MHz"][key]:
            print(ex,end=", ")
        print("")
    
    print("\n \n850 Mhz experiments:\n")
    for key in experiments["850 MHz"]:
        print(key,end=" - ")
        for ex in experiments["600 MHz"][key]:
            print(ex,end=", ")
        print("")
        
    return experiments
    
    
#added 11/10/2022    
def check_grouping_for_LP_t1t2(outputs,coeffs,comparison):
    
    fig = plt.gcf()
    height=5*int(math.ceil((len(outputs)*len(coeffs)+1)/3))
    fig.set_size_inches(18.5, height)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.rcParams["figure.figsize"] = [18.5, height]

    
    pos=1
    highlight=[0,0]
    ##########################################
    """Making this part as a subfunction does not work while called inside of a function, 
    works in the main code. Same code repeats in the below part"""
    ##########################################
    ax = fig.add_subplot(int(math.ceil((len(outputs)*len(coeffs)+1)/3)),3,pos)
        

    if (highlight[0]==comparison["1"]['REFERENCE']["COEFFICIENTS"] and highlight[1]==comparison["1"]['REFERENCE']["OUTPUT_POINTS"]):
        ax.set_title("REFER: Coef: "+str(highlight[0])+", Out: "+str(highlight[1]),fontsize=15)
    else:
        ax.set_title("Coef: "+str(highlight[0])+", Out: "+str(highlight[1]),fontsize=15)
        
    #ax.set_xlabel("H [ppm]",fontsize=10)
    #ax.set_ylabel("N [ppm]",fontsize=10)

    ax.invert_xaxis()
    ax.invert_yaxis()
       
    for i,peak in enumerate(comparison):
        ax.plot(float(comparison[peak]['REFERENCE']["ppm"][1]), float(comparison[peak]['REFERENCE']["ppm"][0]),'o',markersize=10,color="red")

    for i,peak in enumerate(comparison):
        for specra in comparison[peak]:
            if comparison[peak][specra]["COEFFICIENTS"]==highlight[0] and comparison[peak][specra]["OUTPUT_POINTS"]==highlight[1]:
                ax.plot(float(comparison[peak][specra]["ppm"][1]), float(comparison[peak][specra]["ppm"][0]),'o',markersize=5,color="blue")

    ############################################
    
    for ind,coe in enumerate(coeffs):
        for jnd,out in enumerate(outputs):
            highlight=[coe,out]
            pos=(len(outputs))*ind+(jnd+1)+1
            
            
            ax = fig.add_subplot(int(math.ceil((len(outputs)*len(coeffs)+1)/3)),3,pos)
        

            if (highlight[0]==comparison["1"]['REFERENCE']["COEFFICIENTS"] and highlight[1]==comparison["1"]['REFERENCE']["OUTPUT_POINTS"]):
                ax.set_title("REFER: Coef: "+str(highlight[0])+", Out: "+str(highlight[1]),fontsize=15)
            else:
                ax.set_title("Coef: "+str(highlight[0])+", Out: "+str(highlight[1]),fontsize=15)

            ax.set_xlabel("H [ppm]",fontsize=10)
            ax.set_ylabel("N [ppm]",fontsize=10)

            ax.invert_xaxis()
            ax.invert_yaxis()

            for i,peak in enumerate(comparison):
                ax.plot(float(comparison[peak]['REFERENCE']["ppm"][1]), float(comparison[peak]['REFERENCE']["ppm"][0]),'o',markersize=10,color="red")

            for i,peak in enumerate(comparison):
                for specra in comparison[peak]:
                    if comparison[peak][specra]["COEFFICIENTS"]==highlight[0] and comparison[peak][specra]["OUTPUT_POINTS"]==highlight[1]:
                        ax.plot(float(comparison[peak][specra]["ppm"][1]), float(comparison[peak][specra]["ppm"][0]),'o',markersize=5,color="blue")


    plt.show()
    
    
## added 11/10/2022   
def convert_peak_number_from_reference(peak,highlight,comparison):
    for spectra in comparison[peak]:
        if (comparison[peak][spectra]['COEFFICIENTS']==highlight[0]
            and comparison[peak][spectra]['OUTPUT_POINTS']==highlight[1]):
            if "PEAK_NAME" in comparison[peak][spectra]:
                return comparison[peak][spectra]["PEAK_NAME"]
            else:
                return peak


def convert_peak_number_to_reference(peak,highlight,comparison):
    for ref_peak in comparison:
        for spectra in comparison[ref_peak]:
            if "PEAK_NAME" in comparison[ref_peak][spectra]:
                if (comparison[ref_peak][spectra]['COEFFICIENTS']==highlight[0]
                    and comparison[ref_peak][spectra]['OUTPUT_POINTS']==highlight[1]
                    and comparison[ref_peak][spectra]["PEAK_NAME"]==peak):
                    return ref_peak
            else:
                if (comparison[ref_peak][spectra]['COEFFICIENTS']==highlight[0]
                    and comparison[ref_peak][spectra]['OUTPUT_POINTS']==highlight[1]):
                    return ref_peak
                
                
                
                
def analyze_single_peak(highlight,peak,save_dest,*results):
    for result in results:
        if result["INFO"]['COEFFICIENTS']==highlight[0] and result["INFO"]['OUTPUT_POINTS']==highlight[1]:
 
            
            fig = plt.gcf()
            fig.set_size_inches(15, 5)
            plt.plot(result['mixing_times'],result["peaks"][peak]['experiment'],"o",markersize=15)
            
            if not result["peaks"][peak]['fit']==None:
                x_fit_data=np.linspace(result['mixing_times'][0],result['mixing_times'][-1],1000)
                y_fit_data=np.exp(x_fit_data*result["peaks"][peak]['fit'][0][0])*np.exp(result["peaks"][peak]['fit'][0][1])
                plt.title("Peak "+peak+ ", Relax. time: "+str(-1/result['peaks'][peak]["fit"][0][0]))
                plt.plot(result["peaks"][peak]['points_used'][0],result["peaks"][peak]['points_used'][1],"o",markersize=7)
            
            
            if not result["peaks"][peak]['fit']==None:
                plt.plot(x_fit_data,y_fit_data,"-")
            plt.show()
            
            change=input("Use different points? Y/N\n")
            if change=="N":
                delete=input("Delete peak? Y/N\n")
                if delete=="Y":
                    del result["peaks"][peak]
                
            
            redo_peaks=False
            while change=="Y":
                redo_peaks=True
                print("Type points to be used, indexing from 0, end by END.")
                new_ind=[]
                add=0
                while not add=="END":
                    add=input("Another point?")
                    if not add=="END":
                        new_ind.append(int(add))
                        
                curve_fit_coef = np.polyfit([result['mixing_times'][i] for i in new_ind], [np.log(result["peaks"][peak]['experiment'][i]) for i in new_ind], 1, full=True)        
                fit=[curve_fit_coef[0].tolist(),curve_fit_coef[1].tolist()]
                
                fig = plt.gcf()
                fig.set_size_inches(15, 5)
 
                x_fit_data=np.linspace(result['mixing_times'][0],result['mixing_times'][-1],1000)
                y_fit_data=np.exp(x_fit_data*fit[0][0])*np.exp(fit[0][1])

                
                plt.title("Peak "+peak+ ", Relax, time: "+ str(-1/fit[0][0]))
                plt.plot(result['mixing_times'],result["peaks"][peak]['experiment'],"o",markersize=15)
                plt.plot([result['mixing_times'][i] for i in new_ind],[(result["peaks"][peak]['experiment'][i]) for i in new_ind],"o",markersize=7)
  
                plt.plot(x_fit_data,y_fit_data,"-")
                plt.show()
                
                change=input("Use different points? Y/N\n")
            if redo_peaks:
                #saveit=input("Save the changes? Y/N\n")
                saveit="Y"
                if saveit=="Y":
                    result["peaks"][peak]['points_used']=[[result['mixing_times'][i] for i in new_ind],[(result["peaks"][peak]['experiment'][i]) for i in new_ind]]
                    result["peaks"][peak]['fit']=fit

                    
                    tosave=[]
                    for spec in results:
                        tosave.append(spec) 
                    with open(save_dest, 'w') as f:
                        yaml.dump(tosave,f, sort_keys=False)
            

## added 14.10.2022
def plot_data_compare_2_exp(comparison,num_spectra):
    keys=["REFERENCE"]
    infos={}
    colors={1:"red",2:"blue",3:"green"}
    for i in range(num_spectra-1):
        keys.append("file"+str(i))
        
       
    
        
    for file in keys:
        x_axis=[]
        y_axis=[]
        y_error=[]
        for peak in comparison:
            if file in comparison[peak]:
                x_axis.append(int(peak))
                if "T1" in comparison[peak][file]:
                    y_axis.append(comparison[peak][file]["T1"])
                    y_error.append(comparison[peak][file]["error"])
                    T1=True
                if "hetNOE" in comparison[peak][file]:
                    y_axis.append(comparison[peak][file]["hetNOE"])
                    T1=False
                
                valid_peak=peak
        try:
            if comparison[valid_peak][file]["INFO"] not in infos:
                infos[comparison[valid_peak][file]["INFO"]]=colors[len(infos)+1]
            if T1:
                plt.errorbar(x_axis, y_axis,yerr=y_error,fmt='o',markersize=5,color=infos[comparison[valid_peak][file]["INFO"]])
            else:
                plt.plot(x_axis, y_axis,'o',markersize=10,color=infos[comparison[valid_peak][file]["INFO"]])
        except:
            pass

    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)
    plt.xlabel('Reference Peak')
    plt.ylabel('T1')
    plt.show()
    
    
## added 14.10.2022   
def add_noe_data(field,name,exp,author,path,files,info,increments,correct,coeffs,outputs,image,imagedim):
    with open("/home/nenciric/Documents/git/COR15A_Tobi/ExperimentalData/experiments.yaml") as yaml_file:
        experiments = yaml.load(yaml_file, Loader=yaml.FullLoader)
    
    if name not in experiments[field]:
        experiments[field][name]={}
    experiments[field][name][exp]={}
    experiments[field][name][exp]["AUTHOR"]=author
    experiments[field][name][exp]["PATH"]=path
    experiments[field][name][exp]["FILES"]=files
    experiments[field][name][exp]["INFO"]=info
    experiments[field][name][exp]["INCREMENTS"]=increments
    experiments[field][name][exp]["CORRECT"]=correct
    experiments[field][name][exp]["COEFF"]=coeffs
    experiments[field][name][exp]["OUTPUTS"]=outputs
    experiments[field][name][exp]["IMAGE"]=image
    experiments[field][name][exp]["IMAGEDIM"]=imagedim
    
    
    with open("/home/nenciric/Documents/git/COR15A_Tobi/ExperimentalData/experiments.yaml", 'w') as f:
        yaml.dump(experiments,f, sort_keys=True)
        
    print("600 Mhz experiments: \n")
    for key in experiments["600 MHz"]:
        print(key,end=" - ")
        for ex in experiments["600 MHz"][key]:
            print(ex,end=", ")
        print("")
    
    print("\n \n850 Mhz experiments:\n")
    for key in experiments["850 MHz"]:
        print(key,end=" - ")
        for ex in experiments["600 MHz"][key]:
            print(ex,end=", ")
        print("")
        
    return experiments
    


def compare_spectra_hetNOE(Nfactor,maxDist,*files):
    
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
        relaxation_times[ref_peak]["REFERENCE"]["ppm"]=reference['peaks'][ref_peak]["NOE"]["ppm"]
        relaxation_times[ref_peak]["REFERENCE"]["hetNOE"]=  reference['peaks'][ref_peak]["hetNOE"]
        relaxation_times[ref_peak]["REFERENCE"]["INFO"]=reference["INFO"]
            
    for k,file in enumerate(files):
        distances=np.zeros([len(reference['peaks']),len(file['peaks'])])
        for i,ref_peak in enumerate(reference['peaks']):
            #print(float(reference[ref_peak]["ppm"][0]))
            for j,file_peak in enumerate(file['peaks']):
                distances[i,j]=(math.dist([float(reference['peaks'][ref_peak]["NOE"]["ppm"][0]), 
                                           Nfactor*float(reference['peaks'][ref_peak]["NOE"]["ppm"][1])],
                                           [float(file['peaks'][file_peak]["NOE"]["ppm"][0]),
                                            Nfactor*float(file['peaks'][file_peak]["NOE"]["ppm"][1])
                                           ])) 
        
        

        reference_keys=list(reference['peaks'])
        for i,file_peak in enumerate(file['peaks']):
            key=np.where(distances[:,i]==min(distances[:,i]))[0][0]
            if np.where(distances[key,:]==min(distances[key,:]))[0][0]==i and np.sqrt(min(distances[key,:]))<maxDist:
                if  reference_keys[key] in relaxation_times:
                    relaxation_times[reference_keys[key]]["file"+str(k)]={}
                    relaxation_times[reference_keys[key]]["file"+str(k)]["ppm"]=file['peaks'][file_peak]["NOE"]["ppm"]
                    relaxation_times[reference_keys[key]]["file"+str(k)]["hetNOE"]=file['peaks'][file_peak]["hetNOE"]
                    relaxation_times[reference_keys[key]]["file"+str(k)]["PEAK_NAME"]=file_peak
                    relaxation_times[reference_keys[key]]["file"+str(k)]["INFO"]=file["INFO"]
             

           
    return relaxation_times
