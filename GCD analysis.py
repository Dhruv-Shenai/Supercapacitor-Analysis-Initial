import numpy as np            #necessary initialisation stuff
import matplotlib.pyplot as plt
from scipy.stats import linregress,sem
from scipy.signal import find_peaks
pathway=input('what is the pathway? (when you upload the file, right click and copy path) ')
currentdensity=input('what is the current density in A/g? ')
currentdensity=float(currentdensity)
current=float(input('what is the mass of the electrode in mg? '))*currentdensity*0.001
isthisnegative=input('Is this a negative electrode? Please input y or n.  ')
cycles=int(input('How many cycles have you inputted? '))

def getdata(pathway):

    with open(pathway, "r") as file:
        first_line = file.readline()
    headings = first_line.split("\t")
    headings = [heading.strip() for heading in headings]
    print(headings)
    voltage_names = ["<Ewe>/V", "Ewe/V", "Ece/V", "Ecell/V", "<Ece>/V", "<Ecell>/V"]
    voltage_name = next((name for name in voltage_names if name in headings), None)
    if voltage_name is None:
      raise ValueError("No valid voltage column name found in the file.")
    desired_headings = ["time/s", voltage_name]
    column_indices = [headings.index(heading) for heading in desired_headings]
    data = np.genfromtxt(pathway, delimiter="\t", skip_header=1)
    time = data[:, column_indices[0]]  # Time values
    voltage = data[:, column_indices[1]]  # Voltage values
    return time, voltage






def getderivative(size,voltage,time):
#need to now differentiate between vertical and diagonal line
#use derivatives
#use a loop
  dVdt= []
  gradtime=[]
  for j in np.arange(1,size, 1): #calculates gradient every 'steps' number of points
    derivative=(voltage[j]-voltage[(j-1)])/(time[j]-time[(j-1)])
    newtime=(time[j]+time[(j-1)])/2
    dVdt.append(derivative)
    gradtime.append(newtime)
#the loop above works out the derivative between every 'steps' points and adds it to an array called dVdt which stores all the gradients
  return dVdt,gradtime


time,voltage= getdata(pathway)
time=time-time[0]
size=np.prod(time.shape) #size of the array, will be used later on a lot (how many data points)
#plots GCD voltage against time
plt.plot(time,voltage)
plt.xlabel('Time/s')
plt.ylabel('Voltage/V')
plt.title('GCD you inputted') #plots the GCD
plt.show()

if isthisnegative == 'y':
    print('To work out data the voltage data will be flipped')
    voltage=-1*voltage
#gets derivative
dVdt,gradtime=getderivative(size,voltage,time)
lengthofdVdt=len(dVdt)
print('length of dVdt is', len(dVdt))
print('length of gradtime',len(gradtime))

#now to plot dVdt Graph (in case anything goes wrong this will be able to tell you. You should see clear spikes)
plt.plot(gradtime,dVdt)
plt.xlabel('Time/s')
plt.ylabel('Rate of change of Voltage/Vs\u207b\u00b9 ')
plt.title('dVdt graph for whole GCD')     #plots the gradient function over the whole GCD
plt.show()


dVdt=np.array(dVdt)
valid_values = dVdt[~np.isnan(dVdt)]
#print('the list of valid values are',valid_values)
valid_values[np.isinf(valid_values)] = np.nan


#checking for problems with data
num_nans = np.sum(np.isnan(valid_values))
total_elements = len(dVdt)
#print(f"Number of NaNs: {num_nans} / Total elements: {total_elements}")
if np.any(np.isinf(valid_values)):
    print("The valid values array contains Inf values.")





threshold=np.nanstd(valid_values)
#print('the threshold for peaks is',threshold)
allpeaks,_=find_peaks(abs(dVdt),height=threshold)
noofpeaks=len(allpeaks)
#print('the initially detected peaks are ', allpeaks)
print('the number of peaks initially detected are', noofpeaks)
#FINDING ALL THE BEGINNING OF THE DISCHARGES
maximumlocations=[]
minimumlocations=[]
if noofpeaks//10 > 0:
    meanevery10peaks = [np.mean(voltage[allpeaks[every10*10:(every10+1)*10]]) for every10 in range(noofpeaks//10)]
else:
    meanevery10peaks = []
    meanevery10peaks.append(np.mean(voltage[allpeaks]))
if noofpeaks//10 > 0:
    for every10 in range(noofpeaks//10):
        for index in np.arange(1,11,1):
            if voltage[allpeaks[every10*10+index-1]]>meanevery10peaks[every10]:
                maximumlocations.append(allpeaks[every10*10+index-1])
            else:
                minimumlocations.append(allpeaks[every10*10+index-1])
if noofpeaks % 10 > 0:
    for index in np.arange(noofpeaks%10,0,-1):
        #print(index)
        if voltage[allpeaks[-1*index]]>meanevery10peaks[-1]:
            maximumlocations.append(allpeaks[-1*index])
        else:
            minimumlocations.append(allpeaks[-1*index])
print('the number of maximums and therefore number of cycles detected is', len(maximumlocations))
print('the number of minimums detected initially are',len(minimumlocations))

#Finding the time jumps

time_differences = np.diff(time)
time_interval = round(np.mean(time_differences[:11]),2)
timejumps,_=find_peaks(time_differences,height=time_interval*5)
#print('these are where I think the timejumps are:', timejumps)
minimumlocations = np.sort(np.concatenate((timejumps, minimumlocations)))
filteredminimumlocations=[minimumlocations[0]]
for i in range(1, len(minimumlocations)):
    if minimumlocations[i] - filteredminimumlocations[-1] > 10:
        filteredminimumlocations.append(minimumlocations[i])
    else:
        if minimumlocations[i] in timejumps:
            filteredminimumlocations=filteredminimumlocations[:-1]
            filteredminimumlocations.append(minimumlocations[i])

minimumlocations=filteredminimumlocations
#print(minimumlocations)
print('the number of minimums now detected after adding the potential holds are', len(minimumlocations))
timecheck=[maximumlocations[i]-minimumlocations[i] for i in range(len(minimumlocations))]
timecheck=abs(np.array(timecheck))
timecheck=np.append(timecheck, [timecheck[-1]])
print('the size of timecheck is', len(timecheck))
print('timecheck array is ', timecheck)
#
#
#
#
#
#
#
#Section of code which separates out the cycles into each discharge cycle
dischargev = {}   #messy array stuff.. dischargev,t,dVdt stores the full discharge curve for every cycle. Verticalv,t,dVdt stores every vertical section for every cycle
discharget = {}   #diagonalv,t,dVdt stores every diagonal section for every cycle (voltage, time and rate of change in voltage(dV/dt))
dischargedVdt = {}
verticaldVdt = {}
verticalt = {}
verticalv= {}
diagonaldVdt = {}
diagonalt = {}
diagonalv= {}
capacity = {}   #will store capacity for every cycle number
capacitance= {} #stores Capacitance for every cycle number
fullcapacitance= {}
IRdrop= {}
ESR = {}
endofcycle=0 #in d2Vdt2 time(reduced by the number of steps you take. e.g. take 5 data points, you should get 4 derviatives if a step of 1 is taken, or 2 if a step of 2 is taken. The geenral formula is 2n+1 where n is the number of derivatives you get and the total is the number of time values)
for cyclenumber in np.arange(0, cycles, 1): #repeats for every cycle
    print('I am cycle number ', cyclenumber)
    #plt.plot(time[maximumlocations[cyclenumber]-(size // (14*(cycles + 1))):maximumlocations[cyclenumber]+(size // (14*(cycles + 1)))],   voltage[maximumlocations[cyclenumber]-(size // (14*(cycles + 1))):maximumlocations[cyclenumber]+(size // (14*(cycles + 1)))])
    #plt.title('Plot of data around peak')
    #plt.show()
    indexofbeginningofdischarge=np.argmax(voltage[maximumlocations[cyclenumber]-(timecheck[cyclenumber]//2):maximumlocations[cyclenumber]+(timecheck[cyclenumber]//2)])
    indexofbeginningofdischarge=maximumlocations[cyclenumber]-(timecheck[cyclenumber]//2)+indexofbeginningofdischarge
    print('the peak voltage for this cycle is', voltage[indexofbeginningofdischarge])
    peakvoltage=voltage[indexofbeginningofdischarge]
    if cyclenumber==cycles-1:
        dischargev[cyclenumber] = voltage[indexofbeginningofdischarge:]  #removes charging part from the array
        discharget[cyclenumber] = time[indexofbeginningofdischarge:]
        dischargedVdt[cyclenumber]=dVdt[indexofbeginningofdischarge:]
        temporarydVdttime=gradtime[indexofbeginningofdischarge:]
        #plt.plot(discharget[cycles-1],dischargev[cycles-1])
        #plt.title('V against t for one cycle')
        #plt.show()
        #plt.plot(discharget[cycles-2],dischargev[cycles-2])
        #plt.title('V against t for one cycle')
        #plt.show()
        #plt.plot(discharget[cycles-3],dischargev[cycles-3])
        #plt.title('V against t for one cycle')
        #plt.show()
    else:
        if minimumlocations[cyclenumber] in timejumps:
            endofcycleindex=minimumlocations[cyclenumber]
        else:
            endofcycleindex=np.argmin(voltage[minimumlocations[cyclenumber]-(timecheck[cyclenumber]//2):minimumlocations[cyclenumber]+(timecheck[cyclenumber]//2)])
            endofcycleindex=endofcycleindex+minimumlocations[cyclenumber]-(timecheck[cyclenumber]//2)
        dischargev[cyclenumber] = voltage[indexofbeginningofdischarge:(endofcycleindex+1)]  #removes charging part from the array
        discharget[cyclenumber] = time[indexofbeginningofdischarge:(endofcycleindex+1)]
        dischargedVdt[cyclenumber]=dVdt[indexofbeginningofdischarge:endofcycleindex]
        temporarydVdttime=gradtime[indexofbeginningofdischarge:endofcycleindex]
        #plt.plot(time[minimumlocations[cyclenumber]-(size // (16*(cycles + 1))):minimumlocations[cyclenumber]+(size // (16*(cycles + 1)))],   voltage[minimumlocations[cyclenumber]-(size // (16*(cycles + 1))):minimumlocations[cyclenumber]+(size // (16*(cycles + 1)))])
        #plt.title('Plot of data around trough')
        #plt.show()
        print('the smallest value of the voltage is',voltage[endofcycleindex])
    temporaryv = dischargev[cyclenumber]  #temporaryv,t,dVdt are needed because the dischargev,t,dVdt are dictionaries. Think of them as a table made of rows and columns and each
    #column corresponds to the dischargev,t,dVdt for each cycle. To edit the values in the column, we need an array because we use array cutting methods throughout this analysis
    #whenever we use something like voltage=voltage[:something] this is cropping an array. So we take each column in the dict and make it an array and then edit that array. Then
    #equal that array back to the column in the dict.
    temporaryt = discharget[cyclenumber]
    temporarydVdt = dischargedVdt[cyclenumber]
    #plt.plot(temporaryt,temporaryv)
    #plt.title('V against t for one cycle')
    #plt.show()
    #plt.plot(temporarydVdttime,temporarydVdt)
    #plt.title('dVdt against t for one cycle')
    #plt.show()

    step_1_of_finding_endofIRdrop = np.argmin(temporarydVdt)
    temporaryv = temporaryv[(step_1_of_finding_endofIRdrop+1):] #edits the voltage to remove the charging part and keep only the discharge part
    temporaryt = temporaryt[(step_1_of_finding_endofIRdrop+1):]
    temporarydVdt = temporarydVdt[step_1_of_finding_endofIRdrop:]
    temporarydVdttime=temporarydVdttime[step_1_of_finding_endofIRdrop:]


    #METHOD OF FINDING IR DROP
    #Current method used is to take the gradient of the full discharge section and find it's mean.
    # We see a spike that corresponds to the Ir drop and a relatively flat, horizontal line which corresponds to the diagonal section.
    # The mean will be effectively the y value of the flat line and so when the gradient first goes to this value, we consider any subsequent data as the diagonal section
    # and any section before it as the vertical, IR drop, part.
    temporarydVdt=np.array(temporarydVdt)
    valid_values = temporarydVdt[~np.isnan(temporarydVdt)]
    valid_values[np.isinf(valid_values)] = np.nan
    threshold = np.nanmean(valid_values) #mean of gradient of discharge part
    #print('threshold for diagonal part is', threshold)
    #plt.plot(temporarydVdttime,temporarydVdt)
    #plt.plot(temporarydVdttime,np.full_like(temporarydVdttime,threshold),color='red', linewidth=3)
    #plt.title('dVdt against t for one cycle')
    #plt.show()
    for i1, loop in enumerate(temporarydVdt): #finds the end of the voltage drop by finding when the gradient goes above the mean
        if 0 > loop > threshold:
            endofIRdrop = i1 #endofIRdrop becomes the index corresponding to what position in the dVdt array it is when it goes above the mean
            #array cutting to extract out the vertical and diagonal sections of the discharge curve
            verticaldVdt[cyclenumber] = temporarydVdt[:endofIRdrop] #removes any discharge part that's associated with the diagonal section so only the vertical section remains
            verticalt[cyclenumber] = temporaryt[:endofIRdrop+1]
            verticalv[cyclenumber]= temporaryv[:endofIRdrop+1]
            diagonaldVdt[cyclenumber] = temporarydVdt[endofIRdrop:] #removes any discharge part that's associated with the vertical section so only the diagonal section remains
            diagonalt[cyclenumber] = temporaryt[endofIRdrop+1:]
            diagonalv[cyclenumber] = temporaryv[endofIRdrop+1:]
            temporaryt=diagonalt[cyclenumber] #The diagonal section needs to be analysed when working out capacity. We need to make the time and voltage values associated with the cycle
            #into an array so we can take the last and first value out. We can't do this with diagonalt as it's a dict not an array. Imagine it as a book, chapters and pages:
            #you can only state the book and what chapter it is (dict and what column) and not the page as well. However making it an array basically keeps only the chapter
            #and then you can state which page in the chapter it is in.
            temporaryv=diagonalv[cyclenumber]
            #plt.plot(temporaryt,temporaryv)
            #plt.title('Diagonal section')
            #plt.show()
            #finding IR drop, capacitance, capacity etc
            vertical = temporaryv[0]  #end of vertical IR drop
            slope, intercept, _, _, _ = linregress(diagonalt[cyclenumber], diagonalv[cyclenumber]) #works out slope of discharge and therefore the capacitance

            capacity[cyclenumber] = 0.5 * currentdensity * (temporaryt[-1] - temporaryt[0])  #works out capacity
            capacitance[cyclenumber]= -2 * currentdensity / slope #equation cancels out the masses: assumes both electrodes have the same mass. then since m1+m2=2*m1=2*m2 then (m1+m2)/m is 2
            IRdrop[cyclenumber]= peakvoltage - vertical
            ESR[cyclenumber]= (peakvoltage - vertical) / (2 * current)
            fullcapacitance[cyclenumber]=-2*currentdensity*(temporaryt[-1]-temporaryt[0])/(temporaryv[-1]-temporaryv[0])
            break
        elif i1==len(temporarydVdt):
            print('No IR drop detected, refer to capacity values. This code will now stop calculating capacitance,ESR,IR and fullcapacitance')
            capacity[cyclenumber] = 0.5 * currentdensity * (temporaryt[-1] - temporaryt[0])
    print('I was cycle number', cyclenumber)

capacity_values=list(capacity.values())
capacitance_values = list(capacitance.values())
cycle_numbers = list(capacitance.keys())
IRdrop_values = list(IRdrop.values())
ESR_values= list(ESR.values())
fullcapacitance_values= list(fullcapacitance.values())
#plots Capacitance against cycle number
plt.plot(cycle_numbers, capacitance_values)
plt.xlabel('Cycle number')
plt.ylabel('Gravimetric Capacitance F/g')
plt.title('Capacitance against cycle number')
#plt.savefig(input('What do you want to save the Capacitance against cycle number graph as? Do not forget to add .png at the end '), dpi=300)
plt.show()

#plots ESR plot against cycle number
plt.plot(cycle_numbers,ESR_values)
plt.title('ESR against cyclenumber in Ω')
#plt.savefig(input('What do you want to save the ESR against cycle number graph as? Do not forget to add .png at the end '), dpi=300)
plt.show()

#plots capacity against cycle number
plt.plot(cycle_numbers,capacity_values)
plt.title('Capacity against cyclenumber')
#plt.savefig(input('What do you want to save the Capacity against cycle number graph as? Do not forget to add .png at the end '), dpi=300)
plt.show()

#plots Full capacitance against cycle number
plt.plot(cycle_numbers,fullcapacitance_values)
plt.title('Full discharge Capacitance against cyclenumber')
plt.xlabel('Cycle number')
plt.ylabel('Capacitance (F/g)')
plt.show()

print('mean capacitance across all cycles is ',np.mean(capacitance_values))
print('mean IR drop across all cycles is ', np.mean(IRdrop_values), ' and the mean ESR is ', np.mean(ESR_values))
print('mean Capacity across all cycles is ',np.mean(capacity_values))

data = np.array(list(zip(cycle_numbers, capacitance_values, capacity_values, IRdrop_values, ESR_values)),
                dtype=[('Cycle Number', int), ('Capacitance', float), ('Capacity', float), ('IR Drop', float), ('ESR', float)])

# Save the structured array to a text file
np.savetxt(input('What is the name of the text file you want to save the analysis results in? Do not forget to add .txt please. '), data, delimiter='\t', fmt='%d\t%.4f\t%.4f\t%.4f\t%.4f', header='Cycle Number\tCapacitance\tCapacity\tIR Drop\tESR', comments='')
