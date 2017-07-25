from pylab import *

invdata = []
for line in open("output"):
    llist = line.split()
    if llist:
        if llist[0] == "ThermoNFL":
            invdata.append(float(llist[4]))

print("Mean INVAR for 50<t<100 fs: "+str(mean(invdata[500:1000])))
    
fig = figure()
ax = fig.add_subplot(111)
ax.plot(linspace(0,100,num=1000), invdata)
ax.set(title="Automatic INVAR Parameter Adjustment",xlabel="Time [fs]", ylabel="INVAR [a.u.]",xlim=(0,100))
show()
