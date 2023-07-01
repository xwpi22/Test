/* -*-  Mode: C++; c-file-style: "gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2011 University of Kansas
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Author: Justin Rohrer <rohrej@ittc.ku.edu>
 *
 * James P.G. Sterbenz <jpgs@ittc.ku.edu>, director
 * ResiliNets Research Group  http://wiki.ittc.ku.edu/resilinets
 * Information and Telecommunication Technology Center (ITTC)
 * and Department of Electrical Engineering and Computer Science
 * The University of Kansas Lawrence, KS USA.
 *
 * Work supported in part by NSF FIND (Future Internet Design) Program
 * under grant CNS-0626918 (Postmodern Internet Architecture),
 * NSF grant CNS-1050226 (Multilayer Network Resilience Analysis and Experimentation on GENI),
 * US Department of Defense (DoD), and ITTC at The University of Kansas.
 */

/*
 * This example program allows one to run ns-3 DSDV, AODV, or OLSR under
 * a typical random waypoint mobility model.
 *
 * By default, the simulation runs for 200 simulated seconds, of which
 * the first 50 are used for start-up time.  The number of nodes is 50.
 * Nodes move according to RandomWaypointMobilityModel with a speed of
 * 20 m/s and no pause time within a 300x1500 m region.  The WiFi is
 * in ad hoc mode with a 2 Mb/s rate (802.11b) and a Friis loss model.
 * The transmit power is set to 7.5 dBm.
 *
 * It is possible to change the mobility and density of the network by
 * directly modifying the speed and the number of nodes.  It is also
 * possible to change the characteristics of the network by changing
 * the transmit power (as power increases, the impact of mobility
 * decreases and the effective density increases).
 *
 * By default, OLSR is used, but specifying a value of 2 for the protocol
 * will cause AODV to be used, and specifying a value of 3 will cause
 * DSDV to be used.
 *
 * By default, there are 10 source/sink data pairs sending UDP data
 * at an application rate of 2.048 Kb/s each.    This is typically done
 * at a rate of 4 64-byte packets per second.  Application data is
 * started at a random time between 50 and 51 seconds and continues
 * to the end of the simulation.
 *
 * The program outputs a few items:
 * - packet receptions are notified to stdout such as:
 *   <timestamp> <node-id> received one packet from <src-address>
 * - each second, the data reception statistics are tabulated and output
 *   to a comma-separated value (csv) file
 * - some tracing and flow monitor configuration that used to work is
 *   left commented inline in the program
 */

#include <fstream>
#include <iostream>
#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/mobility-module.h"
#include "ns3/aodv-module.h"
#include "ns3/olsr-module.h"
#include "ns3/dsdv-module.h"
#include "ns3/dsr-module.h"
#include "ns3/applications-module.h"
#include "ns3/yans-wifi-helper.h"
#include <vector>
#include <string>
#include "ns3/config-store-module.h"
#include "ns3/energy-module.h"
#include "ns3/wifi-radio-energy-model-helper.h"
#include "ns3/flow-monitor-module.h"
#include "ns3/basic-energy-source-helper.h"
#include "ns3/basic-energy-source.h"
#include "ns3/wifi-radio-energy-model.h"
#include "ns3/pyviz.h"
#include "ns3/athstats-helper.h"


using namespace ns3;
using namespace dsr;
static bool g_verbose = true;
int numNodes = 20;
NS_LOG_COMPONENT_DEFINE ("olsr_simulation");

/**
 * Routing experiment class.
 * 
 * It handles the creation and run of an experiment.
 */
class RoutingExperiment
{
public:
  RoutingExperiment ();
  /**
   * Run the experiment.
   * \param nSinks The number of Sink Nodes.
   * \param txp The Tx power.
   * \param CSVfileName The output CSV filename.
   */
  void Run (int nSinks, double txp, std::string CSVfileName);
  //static void SetMACParam (ns3::NetDeviceContainer & devices,
  //                                 int slotDistance);
  /**
   * Handles the command-line parmeters.
   * \param argc The argument count.
   * \param argv The argument vector.
   * \return the CSV filename.
   */
  std::string CommandSetup (int argc, char **argv);

private:
  /**
   * Setup the receiving socket in a Sink Node.
   * \param addr The address of the node.
   * \param node The node pointer.
   * \return the socket.
   */
  Ptr<Socket> SetupPacketReceive (Ipv4Address addr, Ptr<Node> node);
  /**
   * Receive a packet.
   * \param socket The receiving socket.
   */
  void ReceivePacket (Ptr<Socket> socket);
  /**
   * Compute the throughput.
   */
  void CheckThroughput ();

  
  uint32_t port;            //!< Receiving port number.
  uint32_t bytesTotal;      //!< Total received bytes.
  uint32_t packetsReceived; //!< Total received packets.

  std::string m_CSVfileName;  //!< CSV filename, record socket receiving info.
  std::string m_energyfileName;  //!< CSV filename, record the init_energy and remaining_energy of each node.
  std::string tr_name;  
  int m_nSinks;               //!< Number of sink nodes.
  std::string m_protocolName; //!< Protocol name.
  double m_txp;               //!< Tx power.
  bool m_traceMobility;       //!< Enable mobility tracing.
  uint32_t m_protocol;        //!< Protocol type, 1: OLSR, 2: AODV, 3: DSDV.

  int rateMcs = 1; //Modulation and Coding Scheme
  double SimulationTime = 50.0;
  double duraion = 44;
  double Prss = -88; //-dBm, signal strength while receiving packets
  double Irss = -88; //-dBm, interfering strength
  double offset = 88; //signal strength while transmitting packets
  //transmitting signal strength to receving node should be larger than EnergyDetectionThreshold
};

RoutingExperiment::RoutingExperiment ()
  : port (9),
    bytesTotal (0),
    packetsReceived (0),
    tr_name ("olsr_ori_20"),
    m_CSVfileName ("olsr_ori_20_output.csv"), //_<numNodes>
    m_energyfileName ("olsr_ori_20_energy.csv"),
    m_traceMobility (false),
    m_protocol (1) //OLSR
{
}

static inline std::string
PrintReceivedPacket (Ptr<Socket> socket, Ptr<Packet> packet, Address senderAddress)
{
  std::ostringstream oss;

  oss << Simulator::Now ().GetSeconds () << " " << socket->GetNode ()->GetId ();

  if (InetSocketAddress::IsMatchingType (senderAddress))
    {
      InetSocketAddress addr = InetSocketAddress::ConvertFrom (senderAddress);
      oss << " received one packet from " << addr.GetIpv4 ();
    }
  else
    {
      oss << " received one packet!";
    }
  return oss.str ();
}

void
RoutingExperiment::ReceivePacket (Ptr<Socket> socket)
{
  Ptr<Packet> packet;
  Address senderAddress;
  while ((packet = socket->RecvFrom (senderAddress)))
    {
      bytesTotal += packet->GetSize ();
      packetsReceived += 1;
      // NS_LOG_UNCOND (PrintReceivedPacket (socket, packet, senderAddress));
    }
}

void
RoutingExperiment::CheckThroughput ()
{
  double kbs = (bytesTotal * 8.0) / 1000;
  bytesTotal = 0;

  std::ofstream out (m_CSVfileName.c_str (), std::ios::app);

  out << (Simulator::Now ()).GetSeconds () << ","
      << kbs << ","
      << packetsReceived << ","
      << m_nSinks << ","
      << m_protocolName << ","
      << m_txp << ""
      << std::endl;

  out.close ();
  packetsReceived = 0;
  //check throughput every 1 second
  Simulator::Schedule (Seconds (1.0), &RoutingExperiment::CheckThroughput, this);
}

Ptr<Socket>
RoutingExperiment::SetupPacketReceive (Ipv4Address addr, Ptr<Node> node)
{
  TypeId tid = TypeId::LookupByName ("ns3::UdpSocketFactory");
  Ptr<Socket> sink = Socket::CreateSocket (node, tid);
  InetSocketAddress local = InetSocketAddress (addr, port);
  sink->Bind (local);
  sink->SetRecvCallback (MakeCallback (&RoutingExperiment::ReceivePacket, this));

  return sink;
}

std::string
RoutingExperiment::CommandSetup (int argc, char **argv)
{
  CommandLine cmd (__FILE__);
  cmd.AddValue ("trName", "The pre-name of the output csv file name", tr_name);
  cmd.AddValue ("CSVfileName", "The name of the CSV output file name", m_CSVfileName);
  cmd.AddValue ("traceMobility", "Enable mobility tracing", m_traceMobility);
  cmd.AddValue ("protocol", "1=OLSR;2=AODV;3=DSDV;4=DSR", m_protocol);
  cmd.AddValue ("numNodes", "Number of nodes", numNodes);
  /* add */
  cmd.AddValue ("SimulationTime", "The total simulation time", SimulationTime);
  cmd.AddValue ("energyfileName", "The name of the energy output file name", m_energyfileName);
  cmd.AddValue ("DataRate", "Mcs code for rate", rateMcs);
  cmd.AddValue ("Prss", "Intended primary received signal strength (dBm)", Prss);
  cmd.AddValue ("Irss", "Intended interfering received signal strength (dBm)", Irss);
  /* add */
  cmd.Parse (argc, argv);
  // enable visualization
  {PyViz v;}
  cmd.AddValue("verbose", "Print trace information if true", g_verbose);
  
  return tr_name+"_output.csv";
}

double NodeToInitialEnergy[100];
double NodeToRemainingEnergy[100];
// double NodeToDelayTime[100];
// double RxPacketSize[100];
// double throughput[100];
Ipv4Address NodeToAddress[100];
// double costEnergy[200];
uint32_t RemainingEnergyNodeId = 0;

/// Trace function for remaining energy at node.
void RemainingEnergy (double oldValue, double remainingEnergy)
{
  // NS_LOG_UNCOND (Simulator::Now ().GetSeconds ()
  //                << "s Current remaining energy = " << remainingEnergy << "J");
  NodeToRemainingEnergy[RemainingEnergyNodeId] = remainingEnergy;
}

/// Trace function for total energy consumption at node.
void
TotalEnergy (double oldValue, double totalEnergy)
{
  // NS_LOG_UNCOND (Simulator::Now ().GetSeconds ()
  //                << "s Total energy consumed by radio = " << totalEnergy << "J");
}
int
main (int argc, char *argv[])
{
  RoutingExperiment experiment;
  std::string CSVfileName = experiment.CommandSetup (argc,argv);
  
  //blank out the last output file and write the column headers
  std::ofstream out (CSVfileName.c_str ());
  out << "SimulationSecond," <<
  "ReceiveRate," <<
  "PacketsReceived," <<
  "NumberOfSinks," <<
  "RoutingProtocol," <<
  "TransmissionPower" <<
  std::endl;
  out.close ();

  int nSinks = 15;
  double txp = 7.5;

  experiment.Run (nSinks, txp, CSVfileName);
  return 0;
}

void
RoutingExperiment::Run (int nSinks, double txp, std::string CSVfileName)
{
  Packet::EnablePrinting ();
  m_nSinks = nSinks;
  m_txp = txp;
  m_CSVfileName = CSVfileName;
  m_energyfileName = tr_name + "_energy.csv";

  double TotalTime = SimulationTime; //simulation time
  std::cout << "TotalTime = " << TotalTime << "\n";
  std::string rate ("48000bps");
  std::string phyMode ("DsssRate11Mbps");
  int nodeSpeed = 5; //in m/s
  int nodePause = 0; //in s
  m_protocolName = "protocol";

  Config::SetDefault  ("ns3::OnOffApplication::PacketSize",StringValue ("512"));
  Config::SetDefault ("ns3::OnOffApplication::DataRate",  StringValue (rate));

  //Set Non-unicastMode rate to unicast mode
  Config::SetDefault ("ns3::WifiRemoteStationManager::NonUnicastMode",StringValue (phyMode));

  NodeContainer adhocNodes;
  adhocNodes.Create (numNodes);

  // setting up wifi phy and channel using helpers
  WifiHelper wifi;
  wifi.SetStandard (WIFI_STANDARD_80211b);
  
  //Yet Another Network Simulator
  YansWifiPhyHelper wifiPhy;
  YansWifiChannelHelper wifiChannel;
  wifiChannel.SetPropagationDelay ("ns3::ConstantSpeedPropagationDelayModel");
  wifiChannel.AddPropagationLoss ("ns3::FriisPropagationLossModel");
  wifiPhy.SetChannel (wifiChannel.Create ());

  // Add a mac and disable rate control
  WifiMacHelper wifiMac;
  wifi.SetRemoteStationManager ("ns3::ConstantRateWifiManager",
                                "DataMode",StringValue (phyMode),
                                "ControlMode",StringValue (phyMode));

  wifiPhy.Set ("TxPowerStart",DoubleValue (txp));
  wifiPhy.Set ("TxPowerEnd", DoubleValue (txp));
  // wifiPhy.Set ("TxPowerLevels",UintegerValue (1));
  // wifiPhy.Set ("TxGain",DoubleValue (offset+Prss));
  // wifiPhy.Set ("RxGain",DoubleValue (0));
  // wifiPhy.Set ("EdThreshold",DoubleValue (-88));
  // wifiPhy.Set ("CcaModel1Threshold",DoubleValue (-88));

  wifiMac.SetType ("ns3::AdhocWifiMac");
  NetDeviceContainer adhocDevices = wifi.Install (wifiPhy, wifiMac, adhocNodes);
  
  MobilityHelper mobilityAdhoc;
  [[maybe_unused]] int64_t streamIndex = 0; // used to get consistent mobility across scenarios

  //assigning node allocation
  ObjectFactory pos;
  pos.SetTypeId ("ns3::RandomRectanglePositionAllocator");
  pos.Set ("X", StringValue ("ns3::UniformRandomVariable[Min=0.0|Max=500.0]"));
  pos.Set ("Y", StringValue ("ns3::UniformRandomVariable[Min=0.0|Max=500.0]"));

  Ptr<PositionAllocator> taPositionAlloc = pos.Create ()->GetObject<PositionAllocator> ();
  streamIndex += taPositionAlloc->AssignStreams (streamIndex);

  std::stringstream ssSpeed;
  ssSpeed << "ns3::UniformRandomVariable[Min=0.0|Max=" << nodeSpeed << "]";
  std::stringstream ssPause;
  ssPause << "ns3::ConstantRandomVariable[Constant=" << nodePause << "]";
  mobilityAdhoc.SetMobilityModel ("ns3::RandomWaypointMobilityModel",
                                  "Speed", StringValue (ssSpeed.str ()),
                                  "Pause", StringValue (ssPause.str ()),
                                  "PositionAllocator", PointerValue (taPositionAlloc));
  mobilityAdhoc.SetPositionAllocator (taPositionAlloc);
  mobilityAdhoc.Install (adhocNodes);
  streamIndex += mobilityAdhoc.AssignStreams (adhocNodes, streamIndex);
  // std::cout<<"streamIndex after addition = "<<streamIndex<<"\n"; 

  /** Energy Model **/
  /***************************************************************************/
  /* energy source */
  BasicEnergySourceHelper basicSourceHelper;
  // configure energy source
  basicSourceHelper.Set ("BasicEnergySourceInitialEnergyJ", DoubleValue (50.));
  // install source
  EnergySourceContainer sources = basicSourceHelper.Install (adhocNodes);
  /* device energy model */
  WifiRadioEnergyModelHelper radioEnergyHelper;
  // configure radio energy model
  radioEnergyHelper.Set ("TxCurrentA", DoubleValue (0.0174));
  // install device model
  DeviceEnergyModelContainer deviceModels = radioEnergyHelper.Install (adhocDevices, sources);
  /***************************************************************************/
  // Record initial energy of each node
  for(int i=0; i<numNodes; i++){
    // std::cout<<"Initial energy\n";
    double initialEnergy = DynamicCast<BasicEnergySource> (sources.Get (i))->GetRemainingEnergy();
    std::cout<<i<<": "<<initialEnergy<<"\n";
    NodeToRemainingEnergy[i] = NodeToInitialEnergy[i] = initialEnergy;
  }
  std::cout<<"\n";

  /** connect trace sources **/
  /***************************************************************************/
  // all sources are connected to node 1
  // energy source
  // Ptr<BasicEnergySource> basicSourcePtr = DynamicCast<BasicEnergySource> (sources.Get (1));
  // basicSourcePtr->TraceConnectWithoutContext ("RemainingEnergy", MakeCallback (&RemainingEnergy));
  // // device energy model
  // Ptr<DeviceEnergyModel> basicRadioModelPtr =
  //   basicSourcePtr->FindDeviceEnergyModels ("ns3::WifiRadioEnergyModel").Get (0);
  // NS_ASSERT (basicRadioModelPtr != NULL);
  // basicRadioModelPtr->TraceConnectWithoutContext ("TotalEnergyConsumption", MakeCallback (&TotalEnergy));
  for(int i=0; i<numNodes; i++){
    Ptr<BasicEnergySource> basicSourcePtr = DynamicCast<BasicEnergySource> (sources.Get (i));
    basicSourcePtr->TraceConnectWithoutContext ("RemainingEnergy", MakeCallback (&RemainingEnergy));
    // device energy model
    Ptr<DeviceEnergyModel> basicRadioModelPtr =
      basicSourcePtr->FindDeviceEnergyModels ("ns3::WifiRadioEnergyModel").Get (0);
    NS_ASSERT (basicRadioModelPtr != NULL);
    basicRadioModelPtr->TraceConnectWithoutContext ("TotalEnergyConsumption", MakeCallback (&TotalEnergy));
  }
  /***************************************************************************/
  
  AodvHelper aodv;
  OlsrHelper olsr;
  DsdvHelper dsdv;
  DsrHelper dsr;
  DsrMainHelper dsrMain;
  Ipv4ListRoutingHelper list;
  InternetStackHelper internet;

  switch (m_protocol)
    {
    case 1:
      list.Add (olsr, 100);
      m_protocolName = "OLSR";
      break;
    case 2:
      list.Add (aodv, 100);
      m_protocolName = "AODV";
      break;
    case 3:
      list.Add (dsdv, 100);
      m_protocolName = "DSDV";
      break;
    case 4:
      m_protocolName = "DSR";
      break;
    default:
      NS_FATAL_ERROR ("No such protocol:" << m_protocol);
    }

  if (m_protocol < 4)
    {
      internet.SetRoutingHelper (list);
      internet.Install (adhocNodes);
    }
  else if (m_protocol == 4)
    {
      internet.Install (adhocNodes);
      dsrMain.Install (dsr, adhocNodes);
    }

  NS_LOG_INFO ("assigning ip address");

  Ipv4AddressHelper addressAdhoc;
  addressAdhoc.SetBase ("192.168.50.0", "255.255.255.0");
  Ipv4InterfaceContainer adhocInterfaces;
  adhocInterfaces = addressAdhoc.Assign (adhocDevices);
  for(int i=0; i<numNodes; i++){
    NodeToAddress[i] = adhocInterfaces.GetAddress (i);
  }
  // keeps sending packets
  OnOffHelper onoff1 ("ns3::UdpSocketFactory",Address ());
  onoff1.SetAttribute ("OnTime", StringValue ("ns3::ConstantRandomVariable[Constant=1.0]"));
  onoff1.SetAttribute ("OffTime", StringValue ("ns3::ConstantRandomVariable[Constant=0.0]"));

  Ptr<UniformRandomVariable> var = CreateObject<UniformRandomVariable> ();
  // install sink on each node
  // int recvNode[nSinks] = {2, 5, 8, 1, 6};
  // int sendNode[nSinks] = {4, 8, 2, 5, 3};
  int recvNode[nSinks];
  int sendNode[nSinks];
  for(int i=0; i<nSinks; i++){
    recvNode[i] = -1;
    sendNode[i] = -1;
  }

  for (int i = 0; i < nSinks; i++)
    {
      //select recvNode and sendNode
      bool chosen1 = false;
      bool chosen2 = false;
      do{
        chosen1 = false;
        recvNode[i] = var->GetValue(0, numNodes);
        // std::cout<<recvNode[i]<<"\n";
        for(int j=0; j<i; j++){
          if(recvNode[j] == recvNode[i]){
            chosen1 = true;
            break;
          }
        }
      }while(chosen1);

      Ptr<Socket> sink = SetupPacketReceive (adhocInterfaces.GetAddress (recvNode[i]),
                                              adhocNodes.Get (recvNode[i]));
      // set remote address to send traffic to
      AddressValue remoteAddress (
          InetSocketAddress (adhocInterfaces.GetAddress (recvNode[i]), port));
      onoff1.SetAttribute ("Remote", remoteAddress);

      // install the OnOffHelper on the random sender
      do{
        chosen2 = false;
        sendNode[i] = var->GetValue(0,numNodes);
        if(sendNode[i] == recvNode[i])
          chosen2 = true;
        else
          for(int j=0; j<i; j++){
            if(sendNode[j] == sendNode[i]){
              chosen2 = true;
              break;
            }
          }
      }while(chosen2);
      std::cout << "recvNode = " << recvNode[i] << ", sendNode = " << sendNode[i] << "\n";
      ApplicationContainer temp = onoff1.Install (adhocNodes.Get (sendNode[i]));
      temp.Start (Seconds (50));
      temp.Stop (Seconds (TotalTime - 10));
    }

  std::stringstream ss;
  ss << numNodes;
  std::string nodes = ss.str ();

  std::stringstream ss2;
  ss2 << nodeSpeed;
  std::string sNodeSpeed = ss2.str ();

  std::stringstream ss3;
  ss3 << nodePause;
  std::string sNodePause = ss3.str ();

  std::stringstream ss4;
  ss4 << rate;
  std::string sRate = ss4.str ();

  //NS_LOG_INFO ("Configure Tracing.");
  //tr_name = tr_name + "_" + m_protocolName +"_" + nodes + "nodes_" + sNodeSpeed + "speed_" + sNodePause + "pause_" + sRate + "rate";

  //AsciiTraceHelper ascii;
  //Ptr<OutputStreamWrapper> osw = ascii.CreateFileStream ( (tr_name + ".tr").c_str());
  //wifiPhy.EnableAsciiAll (osw);
  AsciiTraceHelper ascii;
  MobilityHelper::EnableAsciiAll (ascii.CreateFileStream (tr_name + ".mob"));

  Ptr<FlowMonitor> flowmon;
  FlowMonitorHelper flowmonHelper;
  flowmon = flowmonHelper.InstallAll ();


  NS_LOG_INFO ("Run Simulation.");

  CheckThroughput ();

  Simulator::Stop (Seconds (TotalTime));
  Simulator::Run ();

  for (DeviceEnergyModelContainer::Iterator iter = deviceModels.Begin (); iter != deviceModels.End (); iter ++)
  {
    double energyConsumed = (*iter)->GetTotalEnergyConsumption ();
    // NS_LOG_UNCOND ("End of simulation (" << Simulator::Now ().GetSeconds ()
    //                 << "s) Total energy consumed by radio = " << energyConsumed << "J");
    // std::cout << (*iter)->GetTypeId() << " " << energyConsumed << "\n";
    NS_ASSERT (energyConsumed <= 60.);
  }

  /***************************************************************************/
  //output energy information
  std::ofstream out (m_energyfileName.c_str ());
  out << "NodeId," <<
    "InitialEnergy," <<
    "RemainingEnergy" <<
    std::endl;
  for (EnergySourceContainer::Iterator iter = sources.Begin (); iter != sources.End (); iter ++)
    {
      double energyRemaining = (*iter)->GetRemainingEnergy ();
      int nodeID = (*iter)->GetNode()->GetId();
      NodeToRemainingEnergy[nodeID] = energyRemaining;
      NS_LOG_UNCOND ("End of simulation (" << Simulator::Now ().GetSeconds ()
                     << "s) Energy remaining of node"<<(*iter)->GetNode()->GetId()<<" = " << energyRemaining << "J");
      //output result to energy.csv
      out << nodeID << ","
          <<NodeToInitialEnergy[nodeID] <<","
          << NodeToRemainingEnergy[nodeID] << ""
          << std::endl;
    }
  out.close ();
  /***************************************************************************/
  flowmon->SerializeToXmlFile ((tr_name + ".flowmon").c_str(), false, false);

  Simulator::Destroy ();
}

