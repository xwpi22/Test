from xml.etree import ElementTree as ET
import sys
import matplotlib.pyplot as pylab
et=ET.parse(sys.argv[1])
et2=ET.parse(sys.argv[2])
result = sys.argv[3] #output file name
print("Result filename = " + result)
bitrates=[]
bitrates2=[]
lossRate=[]
lossRate2=[]
delays=[]
delays2=[]
flowNO=[]
flowNO2=[]
for flow in et.findall("FlowStats/Flow"):
	flowNO.append(len(flowNO)+1)
	for tpl in et.findall("Ipv4FlowClassifier/Flow"):
		if tpl.get('flowId')==flow.get('flowId'):
			break
	if tpl.get('destinationPort')=='654':
		continue
	
	txPackets = int(flow.get('txPackets'))
	losses = int(flow.get('lostPackets'))
	lossRate.append(round(losses/txPackets*100, 2))
	
	rxPackets=int(flow.get('rxPackets'))
	if rxPackets==0:
		bitrates.append(0)
		delays.append(-1)
	else:
		t0=float(flow.get('timeFirstRxPacket')[:-2])
		t1=float(flow.get("timeLastRxPacket")[:-2])
		duration=(t1-t0)*1e-9
		bitrates.append(8*int(flow.get("rxBytes"))/duration*1e-3)
		delays.append(float(flow.get('delaySum')[:-2])*1e-6/rxPackets)

for flow in et2.findall("FlowStats/Flow"):
	flowNO2.append(len(flowNO2)+1)
	for tpl in et2.findall("Ipv4FlowClassifier/Flow"):
		if tpl.get('flowId')==flow.get('flowId'):
			break
	if tpl.get('destinationPort')=='654':
		continue
	
	txPackets = int(flow.get('txPackets'))
	losses = int(flow.get('lostPackets'))
	lossRate2.append(round(losses/txPackets*100, 2))
	
	rxPackets=int(flow.get('rxPackets'))
	if rxPackets==0:
		bitrates2.append(0)
		delays2.append(-1)
	else:
		t0=float(flow.get('timeFirstRxPacket')[:-2])
		t1=float(flow.get("timeLastRxPacket")[:-2])
		duration=(t1-t0)*1e-9
		bitrates2.append(8*int(flow.get("rxBytes"))/duration*1e-3)
		delays2.append(float(flow.get('delaySum')[:-2])*1e-6/rxPackets)

# This figure has 3 rows, 1 column, and this is the first plot.
pylab.subplot(311)
pylab.plot(flowNO, bitrates, marker = "s", label = 'original')
pylab.plot(flowNO, bitrates2, marker = "o", label = 'adjusted')
pylab.legend(loc = 'upper right', bbox_to_anchor=(1.1, 1.7))
# pylab.hist(bitrates,bins=40)
pylab.xticks(flowNO)
pylab.xlabel("Flow ID")
pylab.ylabel("Flow Bit Rates (kb/s)")
pylab.margins(0.05, 0.5)
for a,b in zip(flowNO, bitrates2):
	pylab.text(a, b, str(round(b,2)))

pylab.subplot(312)
pylab.plot(flowNO, lossRate, marker = "s", label = 'original')
pylab.plot(flowNO, lossRate2, marker = "o", label = 'adjusted')
pylab.legend(loc = 'upper right', bbox_to_anchor=(1.1, 1.5))
# pylab.hist(losses,bins=40)
pylab.xticks(flowNO)
pylab.xlabel("Flow ID")
pylab.ylabel("Loss Rate")
pylab.margins(0.05, 0.5)
for a,b in zip(flowNO, lossRate2):
	pylab.text(a, b, str(b))

pylab.subplot(313)
pylab.plot(flowNO, delays, marker = "s", label = 'original')
pylab.plot(flowNO, delays2, marker = "o", label = 'adjusted')
pylab.legend(loc = 'upper right', bbox_to_anchor=(1.1, 1.5))
# pylab.hist(delays, bins=40)
pylab.xticks(flowNO)
pylab.xlabel("Flow ID")
pylab.ylabel("Delay in ms")
pylab.margins(0.05, 0.5)
for a,b in zip(flowNO, delays2):
	pylab.text(a, b, str(round(b,2)))

pylab.subplots_adjust(hspace=1)
# pylab.savefig("results_adj_20.pdf")
# pylab.savefig("result_100.pdf")
pylab.savefig(result)
