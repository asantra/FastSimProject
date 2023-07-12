#### code to draw plots comparing FullSim and FastSimpling at sampling plane.
#### The plots use LUXE dump geometry

import os, sys, time
from syslog import LOG_SYSLOG
import argparse
from ROOT import *
from copy import copy, deepcopy
sys.path.insert(0, '/Users/arkasantra/arka/include')
from Functions import *
from collections import OrderedDict

def DrawHistsRatio(FirstTH1, LegendName, PlotColor, xrange1down, xrange1up, yrange1down, yrange1up, xaxisTitle, CanvasName, h2, yline1low, yline1up, drawline=False, logy=False, LatexName='', LatexName2='', TeVTag=False, doSumw2=False, doAtlas=False, doLumi=False, noRatio=False, do80=False, do59=False,drawOption=""):
   Tex = MakeLatex(0.85,0.65,LatexName)
   Tex2 = MakeLatex(0.85,0.60,LatexName2)
   c = TCanvas("c","c",900, 900)
   c.SetGrid()
   pad1 = TPad('pad1','pad1',0,0.25,1,1)
   pad1.SetBottomMargin(0.0)
   pad1.SetFillStyle(0)
   pad1.SetGrid()
   pad1.Draw()
   pad1.cd()
   gStyle.SetOptStat(0)
   if(logy):
     pad1.SetLogy()
     
   if "energy" in FirstTH1[0].GetName() or "time" in FirstTH1[0].GetName():
    pad1.SetLogx()

   line = MakeLine(xrange1down,yline1low,xrange1up,yline1up)
   legend1 = LegendMaker()
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L, doAtlas, doLumi, noRatio, do80, do59)
   FirstTH1[0].GetYaxis().SetRangeUser(yrange1down,yrange1up)
   FirstTH1[0].GetXaxis().SetRangeUser(xrange1down,xrange1up)
   WaterMark = TexWaterMark('Preliminary')
   for i in range(0, len(FirstTH1)):
     FirstTH1[i].SetTitle("")
     FirstTH1[i].GetYaxis().SetTitle(FirstTH1[i].GetYaxis().GetTitle())
     
     FirstTH1[i].GetYaxis().SetTitleOffset(0.95)
     FirstTH1[i].GetYaxis().SetTitleSize(0.05)
     FirstTH1[i].GetYaxis().SetLabelSize(0.045)
     FirstTH1[i].GetXaxis().SetTitle("")
     FirstTH1[i].GetXaxis().SetLabelSize(0.15)
     if "track_r" in FirstTH1[i].GetName() or "track_x" in FirstTH1[i].GetName() or "track_y" in FirstTH1[i].GetName():
        FirstTH1[i].Scale(1.0, "width")
     w = FirstTH1[i].Integral()
     #FirstTH1[i].Scale(1./w)
    #  print("integral inside: ", w)
     FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
     legend1.AddEntry(FirstTH1[i],LegendName[i], "lp")
   
   
   FirstTH1[0].GetXaxis().SetLabelSize(0.0);
   FirstTH1[0].GetYaxis().SetTitle("Entries");
   FirstTH1[0].SetLineWidth(2)
   FirstTH1[0].Draw("hist")
   print("integral outside 0: ", FirstTH1[0].Integral())
    
   #WaterMark.Draw("sames")
   gPad.SetTickx()
   gPad.SetTicky()
   gPad.Modified(); 
   gPad.Update();
   
   #### special file for Allpix-squared plot, otherwise remove -1 in the range
   for i in range(1, len(FirstTH1)):
     FirstTH1[i].SetLineWidth(1)
     FirstTH1[i].Draw("hist sames")
     print("integral outside 1: ", FirstTH1[i].Integral())
     gPad.Modified(); 
     gPad.Update();
     
   FirstTH1[0].SetLineWidth(1)
   FirstTH1[0].Draw("hist sames")
    
   gPad.RedrawAxis()
   L[1].Draw("sames")
   L[2].Draw("sames")
   if(TeVTag):
       TexTeV = TLatex(0.892,0.914,"#sqrt{s}=13 TeV")
       TexTeV.SetTextAlign(31)
       TexTeV.SetTextFont(42)
       TexTeV.SetTextSize(0.037)
       TexTeV.SetLineWidth(2) 
       TexTeV.Draw()
   else:
       L[0].Draw()
   legend1.Draw("sames")
   Tex.Draw("sames")
   Tex2.Draw("sames")
   
   c.cd()
   pad2 = TPad("pad2", "pad2",0.0,0.0,1.0,0.245)
   pad2.SetTopMargin(0.0)
   pad2.SetBottomMargin(0.3)
   pad2.SetFillStyle(0)
   #pad2.SetGrid()
   pad2.Draw()
   pad2.cd()
   gStyle.SetOptStat(0)
   gPad.SetTickx()
   gPad.SetTicky()

   if "energy" in FirstTH1[0].GetName() or "time" in FirstTH1[0].GetName():
    pad2.SetLogx()
   
   if(doSumw2):
       for i in range(0, len(FirstTH1)):
          FirstTH1[i].Sumw2() 
   
   h2.Divide(FirstTH1[0],FirstTH1[1])
   h2.GetXaxis().SetRangeUser(xrange1down,xrange1up)
   h2.GetXaxis().SetTitle(xaxisTitle)
   h2.GetYaxis().SetTitle(h2.GetYaxis().GetTitle())
   h2.GetYaxis().CenterTitle()
   h2.SetTitle('')
   h2 = OneRatioAxisSize(h2)
   h2.GetYaxis().SetRangeUser(0.0, 2.0)
   h2.GetYaxis().SetTitleOffset(0.45)
   h2.Draw("ep")
   if(drawline):
     line.SetLineStyle(2)
     line.Draw()
     
   SaveFile(c, CanvasName)
   return deepcopy(c)


def main():
    gROOT.SetBatch()

    parser = argparse.ArgumentParser(description='Code to get 2D plots')
    parser.add_argument('-det', action="store", dest="detid", type=str, default="33")
    parser.add_argument('-times', action="store", dest="weightValue", type=str, default="1.0")
    parser.add_argument('-timesNt', action="store", dest="weightValueNt", type=str, default="1.0")
    parser.add_argument('-timesPh', action="store", dest="weightValuePh", type=str, default="1.0")
    parser.add_argument('-version', action="store", dest="versionVal", type=str, default="1")
    args = parser.parse_args()

    # outDir = "RandomFastSamplingvsFullSimLUXE_CoarseBinning_"+args.weightValue+"timesBackwardInThetaMore3AndRLess300" ### when neutron and photon have same weights
    outDir = "RandomFastSamplingvsFullSimLUXE_CoarseBinning_"+args.weightValueNt+"timesNeutron_"+args.weightValuePh+"timesPhoton_v"+args.versionVal ### when neutron and photon have different weights

    if not os.path.exists(outDir):
        os.makedirs(outDir)

    ### this is Dump only geometry
    # if int(args.detid)==33:
    #     zPos = -350
    #     ###zPos = 6621.91
    # elif int(args.detid)==32:
    #     zPos = -5000
    # elif int(args.detid)==31:
    #     zPos = -10000
    # elif int(args.detid)==30:
    #     zPos = -15000
    # else:
    #     zPos = -350

    ### this is LUXE geometry
    if int(args.detid)==33:
        zPos = 6621.91;
    elif int(args.detid)==32:
        zPos = 5450.25;
    elif int(args.detid)==31:
        zPos = 4125;
    elif int(args.detid)==30:
        zPos = 0;
    else:
        zPos = 6621.91;
    
    inDir = "/Users/arkasantra/arka/Sasha_Work/OutputFile"
    #fullSimFile = TFile(inDir+"/RestrictedDumpOnlyFiles_DetId"+args.detid+"_trackInfo.root","READ")
    #fastSimFile = TFile(inDir+"/RestrictedDumpOnlyFiles_DetId33_trackInfo_RandomGeneration_v7_AllParts.root", "READ")
    try:
        fullSimFile = TFile(inDir+"/LUXEDumpFiles_FullSim_0p06BX_DetId"+args.detid+"_NoECutNtrn_CoarseBinning_1DComparePlot.root", "READ")
        
        # fastSimFile = TFile(inDir+"/LUXEDumpFiles_FullSim_0p06BX_DetId"+args.detid+"_NoECutNtrn_CoarseBinning_"+args.weightValue+"timesBackwardInThetaMore3AndRLess300_RandomGeneration_v1.root", "READ") ### when neutron and photon have same weights
        # fastSimFile = TFile(inDir+"/LUXEDumpFiles_FullSim_0p06BX_DetId"+args.detid+"_NoECutNtrn_CoarseBinning_"+args.weightValueNt+"timesNeutron_"+args.weightValuePh+"timesPhoton_BackwardInThetaMore3AndRLess300_RandomGeneration_v1.root", "READ") ### when neutron and photon have different weights
        fastSimFile = TFile(inDir+"/LUXEDumpFiles_FullSim_0p06BX_DetId"+args.detid+"_NoECutNtrn_CoarseBinning_"+args.weightValueNt+"timesNeutron_"+args.weightValuePh+"timesPhoton_v"+args.versionVal+"_RandomGeneration_v1.root", "READ") ### when neutron and photon have different weights
    except:
        print("Something wrong in the files")
        exit()

    dump_plane_bkg_track_rUp_photon_cut              = fullSimFile.Get("dump_plane_bkg_track_rUp_photon_cut")
    dump_plane_bkg_track_rDn_photon_cut              = fullSimFile.Get("dump_plane_bkg_track_rDn_photon_cut")
    dump_plane_bkg_track_thetaUp_photon_weighted_cut = fullSimFile.Get("dump_plane_bkg_track_thetaUp_photon_weighted_cut")
    dump_plane_bkg_track_thetaDn_photon_weighted_cut = fullSimFile.Get("dump_plane_bkg_track_thetaDn_photon_weighted_cut")
    dump_plane_bkg_track_thetaUp_photon_weighted_cut_morebins = fullSimFile.Get("dump_plane_bkg_track_thetaUp_photon_weighted_cut_morebins")
    dump_plane_bkg_track_thetaDn_photon_weighted_cut_morebins = fullSimFile.Get("dump_plane_bkg_track_thetaDn_photon_weighted_cut_morebins")
    dump_plane_bkg_track_energy_photon_cut           = fullSimFile.Get("dump_plane_bkg_track_energy_photon_cut")
    dump_plane_bkg_track_time_photon_cut             = fullSimFile.Get("dump_plane_bkg_track_time_photon_cut")
    dump_plane_bkg_track_phi_photon_cut              = fullSimFile.Get("dump_plane_bkg_track_phi_photon_cut")
    dump_plane_bkg_track_x_photon_cut                = fullSimFile.Get("dump_plane_bkg_track_x_photon_cut")
    dump_plane_bkg_track_y_photon_cut                = fullSimFile.Get("dump_plane_bkg_track_y_photon_cut")
    dump_plane_bkg_track_xUp_photon_cut              = fullSimFile.Get("dump_plane_bkg_track_xUp_photon_cut")
    dump_plane_bkg_track_xDn_photon_cut              = fullSimFile.Get("dump_plane_bkg_track_xDn_photon_cut")

    dump_plane_bkg_track_rUp_neutron_cut             = fullSimFile.Get("dump_plane_bkg_track_rUp_neutron_cut")
    dump_plane_bkg_track_rDn_neutron_cut             = fullSimFile.Get("dump_plane_bkg_track_rDn_neutron_cut")
    dump_plane_bkg_track_thetaUp_neutron_weighted_cut  = fullSimFile.Get("dump_plane_bkg_track_thetaUp_neutron_weighted_cut")
    dump_plane_bkg_track_thetaDn_neutron_weighted_cut  = fullSimFile.Get("dump_plane_bkg_track_thetaDn_neutron_weighted_cut")
    dump_plane_bkg_track_thetaUp_neutron_weighted_cut_morebins  = fullSimFile.Get("dump_plane_bkg_track_thetaUp_neutron_weighted_cut_morebins")
    dump_plane_bkg_track_thetaDn_neutron_weighted_cut_morebins  = fullSimFile.Get("dump_plane_bkg_track_thetaDn_neutron_weighted_cut_morebins")
    dump_plane_bkg_track_energy_neutron_cut          = fullSimFile.Get("dump_plane_bkg_track_energy_neutron_cut")
    dump_plane_bkg_track_time_neutron_cut            = fullSimFile.Get("dump_plane_bkg_track_time_neutron_cut")
    dump_plane_bkg_track_phi_neutron_cut             = fullSimFile.Get("dump_plane_bkg_track_phi_neutron_cut")
    dump_plane_bkg_track_x_neutron_cut               = fullSimFile.Get("dump_plane_bkg_track_x_neutron_cut")
    dump_plane_bkg_track_y_neutron_cut               = fullSimFile.Get("dump_plane_bkg_track_y_neutron_cut")
    dump_plane_bkg_track_xUp_neutron_cut             = fullSimFile.Get("dump_plane_bkg_track_xUp_neutron_cut")
    dump_plane_bkg_track_xDn_neutron_cut             = fullSimFile.Get("dump_plane_bkg_track_xDn_neutron_cut")
    



    # ///### 1D plot
    # allHisto1Dict.insert(make_pair("dump_plane_bkg_track_time_neutron_cut", new TH1D("dump_plane_bkg_track_time_neutron_cut","dump_plane_bkg_track_time_neutron_cut; t [ns]; Events",nbins, tarray)));
    # allHisto1Dict.insert(make_pair("dump_plane_bkg_track_energy_neutron_cut", new TH1D("dump_plane_bkg_track_energy_neutron_cut","dump_plane_bkg_track_energy_neutron_cut; E [GeV]; Events",nbins, xarray)));
    # allHisto1Dict.insert(make_pair("dump_plane_bkg_track_r_neutron_cut", new TH1D("dump_plane_bkg_track_r_neutron_cut","dump_plane_bkg_track_r_neutron_cut; r [mm]; Events",nRBins, rarray)));
    # allHisto1Dict.insert(make_pair("dump_plane_bkg_track_r_small_neutron_cut", new TH1D("dump_plane_bkg_track_r_small_neutron_cut","dump_plane_bkg_track_r_small_neutron_cut; r [mm]; Events",500,0, 500)));
    # allHisto1Dict.insert(make_pair("dump_plane_bkg_track_theta_neutron_cut", new TH1D("dump_plane_bkg_track_theta_neutron_cut","dump_plane_bkg_track_theta_neutron_cut; #theta_{p} [rad]; Events",6400, 1.6, 3.2)));
    # allHisto1Dict.insert(make_pair("dump_plane_bkg_track_theta_neutron_weighted_cut", new TH1D("dump_plane_bkg_track_theta_neutron_weighted_cut","dump_plane_bkg_track_theta_neutron_weighted_cut; #theta_{p} [rad]; Events",6400, 1.6, 3.2)));
    # allHisto1Dict.insert(make_pair("dump_plane_bkg_track_phi_neutron_cut", new TH1D("dump_plane_bkg_track_phi_neutron_cut","dump_plane_bkg_track_phi_neutron_cut; #phi_{p} [rad]; Events",640, -3.2, 3.2)));

    
    
    # allHisto1Dict.insert(make_pair("dump_plane_bkg_track_time_photon_cut", new TH1D("dump_plane_bkg_track_time_photon_cut","dump_plane_bkg_track_time_photon_cut; t [ns]; Events",nbins, tarray)));
    # allHisto1Dict.insert(make_pair("dump_plane_bkg_track_energy_photon_cut", new TH1D("dump_plane_bkg_track_energy_photon_cut","dump_plane_bkg_track_energy_photon_cut; E [GeV]; Events",nbins, xarray)));
    # allHisto1Dict.insert(make_pair("dump_plane_bkg_track_r_photon_cut", new TH1D("dump_plane_bkg_track_r_photon_cut","dump_plane_bkg_track_r_photon_cut; r [mm]; Events",nRBins, rarray)));
    # allHisto1Dict.insert(make_pair("dump_plane_bkg_track_r_small_photon_cut", new TH1D("dump_plane_bkg_track_r_small_photon_cut","dump_plane_bkg_track_r_small_photon_cut; r [mm]; Events",500, 0, 500)));
    # allHisto1Dict.insert(make_pair("dump_plane_bkg_track_theta_photon_weighted_cut", new TH1D("dump_plane_bkg_track_theta_photon_weighted_cut","dump_plane_bkg_track_theta_photon_weighted_cut; #theta_{p} [rad]; Events",6400, 1.6, 3.2)));
    # allHisto1Dict.insert(make_pair("dump_plane_bkg_track_theta_photon_cut", new TH1D("dump_plane_bkg_track_theta_photon_cut","dump_plane_bkg_track_theta_photon_cut; #theta_{p} [rad]; Events",6400, 1.6, 3.2)));
    # allHisto1Dict.insert(make_pair("dump_plane_bkg_track_phi_photon_cut", new TH1D("dump_plane_bkg_track_phi_photon_cut","dump_plane_bkg_track_phi_photon_cut; #phi_{p} [rad]; Events",640, -3.2, 3.2)));


    dump_plane_bkg_track_rUp_photon_weighted        = fastSimFile.Get("dump_plane_bkg_track_rUp_photon_weighted")
    dump_plane_bkg_track_rDn_photon_weighted        = fastSimFile.Get("dump_plane_bkg_track_rDn_photon_weighted")
    dump_plane_bkg_track_thetaUp_photon             = fastSimFile.Get("dump_plane_bkg_track_thetaUp_photon")
    dump_plane_bkg_track_thetaDn_photon             = fastSimFile.Get("dump_plane_bkg_track_thetaDn_photon")
    dump_plane_bkg_track_thetaUp_photon_morebins             = fastSimFile.Get("dump_plane_bkg_track_thetaUp_photon_morebins")
    dump_plane_bkg_track_thetaDn_photon_morebins             = fastSimFile.Get("dump_plane_bkg_track_thetaDn_photon_morebins")
    dump_plane_bkg_track_E_photon                   = fastSimFile.Get("dump_plane_bkg_track_E_photon")
    dump_plane_bkg_track_time_photon                = fastSimFile.Get("dump_plane_bkg_track_time_photon")
    dump_plane_bkg_track_phi_photon                 = fastSimFile.Get("dump_plane_bkg_track_phi_photon")
    dump_plane_bkg_track_x_photon                   = fastSimFile.Get("dump_plane_bkg_track_x_photon_cut")
    dump_plane_bkg_track_y_photon                   = fastSimFile.Get("dump_plane_bkg_track_y_photon_cut")
    dump_plane_bkg_track_xUp_photon                 = fastSimFile.Get("dump_plane_bkg_track_xUp_photon")
    dump_plane_bkg_track_xDn_photon                 = fastSimFile.Get("dump_plane_bkg_track_xDn_photon")
    
    

    dump_plane_bkg_track_rUp_neutron_weighted       = fastSimFile.Get("dump_plane_bkg_track_rUp_neutron_weighted")
    dump_plane_bkg_track_rDn_neutron_weighted       = fastSimFile.Get("dump_plane_bkg_track_rDn_neutron_weighted")
    dump_plane_bkg_track_thetaUp_neutron            = fastSimFile.Get("dump_plane_bkg_track_thetaUp_neutron")
    dump_plane_bkg_track_thetaDn_neutron            = fastSimFile.Get("dump_plane_bkg_track_thetaDn_neutron")
    dump_plane_bkg_track_thetaUp_neutron_morebins            = fastSimFile.Get("dump_plane_bkg_track_thetaUp_neutron_morebins")
    dump_plane_bkg_track_thetaDn_neutron_morebins            = fastSimFile.Get("dump_plane_bkg_track_thetaDn_neutron_morebins")
    dump_plane_bkg_track_E_neutron                  = fastSimFile.Get("dump_plane_bkg_track_E_neutron")
    dump_plane_bkg_track_time_neutron               = fastSimFile.Get("dump_plane_bkg_track_time_neutron")
    dump_plane_bkg_track_phi_neutron                = fastSimFile.Get("dump_plane_bkg_track_phi_neutron")
    dump_plane_bkg_track_x_neutron                  = fastSimFile.Get("dump_plane_bkg_track_x_neutron_cut")
    dump_plane_bkg_track_y_neutron                  = fastSimFile.Get("dump_plane_bkg_track_y_neutron_cut")
    dump_plane_bkg_track_xUp_neutron                = fastSimFile.Get("dump_plane_bkg_track_xUp_neutron")
    dump_plane_bkg_track_xDn_neutron                = fastSimFile.Get("dump_plane_bkg_track_xDn_neutron")




    
    #### plotting the signal and background
    TeVTag=False; doSumw2=False; doAtlas=False; doLumi=False; noRatio=False; do80=False; do59=False; drawOption=""
    drawline    = True
    logy        = True
    latexName2  = 'z='+str(zPos)+' mm'
    leftLegend  = True

    LegendName  = ["FullSim", "Fast Sampling"]
    PlotColor   = [2, 4]


    latexName   = "photon"

    FirstTH1    = [dump_plane_bkg_track_rUp_photon_cut, dump_plane_bkg_track_rUp_photon_weighted]
    
    
    for i in range(len(FirstTH1)):
        FirstTH1[i].Rebin(20)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e3
    yAxisLow    = 0.5
    
    xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")
    
    

    FirstTH1    = [dump_plane_bkg_track_xUp_photon_cut, dump_plane_bkg_track_xUp_photon]
    
    
    for i in range(len(FirstTH1)):
        FirstTH1[i].Rebin(20)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e3
    yAxisLow    = 0.5
    
    xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")


    FirstTH1    = [dump_plane_bkg_track_xDn_photon_cut, dump_plane_bkg_track_xDn_photon]
    
    
    for i in range(len(FirstTH1)):
        FirstTH1[i].Rebin(20)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e3
    yAxisLow    = 0.5
    
    xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")

    
    FirstTH1    = [dump_plane_bkg_track_x_photon_cut, dump_plane_bkg_track_x_photon]
    
    for i in range(len(FirstTH1)):
        FirstTH1[i].Rebin(20)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e3
    yAxisLow    = 0.5
    
    xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")
    
    
    FirstTH1    = [dump_plane_bkg_track_y_photon_cut, dump_plane_bkg_track_y_photon]
    
    for i in range(len(FirstTH1)):
        FirstTH1[i].Rebin(20)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e3
    yAxisLow    = 0.5
    
    xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")
    
    
    
    
    
    FirstTH1    = [dump_plane_bkg_track_rDn_photon_cut, dump_plane_bkg_track_rDn_photon_weighted]
    for i in range(len(FirstTH1)):
        FirstTH1[i].Rebin(20)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e3
    yAxisLow    = 0.5
    
    xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")






    FirstTH1    = [dump_plane_bkg_track_thetaUp_photon_weighted_cut, dump_plane_bkg_track_thetaUp_photon]
    for i in range(len(FirstTH1)):
        FirstTH1[i].Rebin(20)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e3
    yAxisLow    = 0.5
    
    xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")



    FirstTH1    = [dump_plane_bkg_track_thetaDn_photon_weighted_cut, dump_plane_bkg_track_thetaDn_photon]
    for i in range(len(FirstTH1)):
        FirstTH1[i].Rebin(20)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e3
    yAxisLow    = 0.5
    
    xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")

    dump_plane_bkg_track_thetaUp_photon_weighted_cut.Add(dump_plane_bkg_track_thetaDn_photon_weighted_cut)
    dump_plane_bkg_track_thetaUp_photon.Add(dump_plane_bkg_track_thetaDn_photon)
    FirstTH1    = [dump_plane_bkg_track_thetaUp_photon_weighted_cut, dump_plane_bkg_track_thetaUp_photon]
    # for i in range(len(FirstTH1)):
    #     FirstTH1[i].Rebin(20)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e3
    yAxisLow    = 0.5
    
    xAxisTitle  = "#theta_{p} [rad]"#FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName()+"CombinedUpDn", h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")


    try:
        FirstTH1    = [dump_plane_bkg_track_thetaUp_photon_weighted_cut_morebins, dump_plane_bkg_track_thetaUp_photon_morebins]
        # for i in range(len(FirstTH1)):
        #     FirstTH1[i].Rebin(5)

        xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
        xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
        
        yAxisHigh   = FirstTH1[0].GetMaximum()*1e3
        yAxisLow    = 0.5
        
        xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

        h2 = FirstTH1[0].Clone("h2")
        h2.Reset()
        h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

        DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")



        FirstTH1    = [dump_plane_bkg_track_thetaDn_photon_weighted_cut_morebins, dump_plane_bkg_track_thetaDn_photon_morebins]
        # for i in range(len(FirstTH1)):
        #     FirstTH1[i].Rebin(5)

        xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
        xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
        
        yAxisHigh   = FirstTH1[0].GetMaximum()*1e3
        yAxisLow    = 0.5
        
        xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

        h2 = FirstTH1[0].Clone("h2")
        h2.Reset()
        h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

        DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")

        dump_plane_bkg_track_thetaUp_photon_weighted_cut_morebins.Add(dump_plane_bkg_track_thetaDn_photon_weighted_cut_morebins)
        dump_plane_bkg_track_thetaUp_photon_morebins.Add(dump_plane_bkg_track_thetaDn_photon_morebins)
        FirstTH1    = [dump_plane_bkg_track_thetaUp_photon_weighted_cut_morebins, dump_plane_bkg_track_thetaUp_photon_morebins]
        for i in range(len(FirstTH1)):
            FirstTH1[i].Rebin(5)

        xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
        xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
        
        yAxisHigh   = FirstTH1[0].GetMaximum()*1e3
        yAxisLow    = 0.5
        
        xAxisTitle  = "#theta_{p} [rad]"#FirstTH1[0].GetXaxis().GetTitle()

        h2 = FirstTH1[0].Clone("h2")
        h2.Reset()
        h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

        DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName()+"CombinedUpDn", h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")
    except:
        print("Something wrong in thetaUp_photon_weighted_cut_morebins")


    FirstTH1    = [dump_plane_bkg_track_phi_photon_cut, dump_plane_bkg_track_phi_photon]
    
    for i in range(len(FirstTH1)):
        FirstTH1[i].Rebin(20)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e2
    yAxisLow    = 0.5
    
    xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")




    FirstTH1    = [dump_plane_bkg_track_time_photon_cut, dump_plane_bkg_track_time_photon]
    for i in range(len(FirstTH1)):
        FirstTH1[i].Rebin(5)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e2
    yAxisLow    = 0.5
    
    xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")
    
    


    FirstTH1    = [dump_plane_bkg_track_energy_photon_cut, dump_plane_bkg_track_E_photon]
    for i in range(len(FirstTH1)):
        FirstTH1[i].Rebin(5)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e2
    yAxisLow    = 0.5
    
    xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")
    
    
    
    
    latexName   = "neutron"

    FirstTH1    = [dump_plane_bkg_track_rUp_neutron_cut, dump_plane_bkg_track_rUp_neutron_weighted]
    for i in range(len(FirstTH1)):
        FirstTH1[i].Rebin(20)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e2
    yAxisLow    = 0.005
    
    xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")
    
    
    

    FirstTH1    = [dump_plane_bkg_track_xUp_neutron_cut, dump_plane_bkg_track_xUp_neutron]
    
    
    for i in range(len(FirstTH1)):
        FirstTH1[i].Rebin(20)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e3
    yAxisLow    = 0.5
    
    xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")


    FirstTH1    = [dump_plane_bkg_track_xDn_neutron_cut, dump_plane_bkg_track_xDn_neutron]
    
    
    for i in range(len(FirstTH1)):
        FirstTH1[i].Rebin(20)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e3
    yAxisLow    = 0.5
    
    xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")

    
    
    
    
    FirstTH1    = [dump_plane_bkg_track_x_neutron_cut, dump_plane_bkg_track_x_neutron]
    
    for i in range(len(FirstTH1)):
        FirstTH1[i].Rebin(20)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e3
    yAxisLow    = 0.5
    
    xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")
    
    
    FirstTH1    = [dump_plane_bkg_track_y_neutron_cut, dump_plane_bkg_track_y_neutron]
    
    for i in range(len(FirstTH1)):
        FirstTH1[i].Rebin(20)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e3
    yAxisLow    = 0.5
    
    xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")
    
    
    
    
    
    
    
    
    FirstTH1    = [dump_plane_bkg_track_rDn_neutron_cut, dump_plane_bkg_track_rDn_neutron_weighted]
    for i in range(len(FirstTH1)):
        FirstTH1[i].Rebin(20)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e2
    yAxisLow    = 0.005
    
    xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")






    FirstTH1    = [dump_plane_bkg_track_thetaUp_neutron_weighted_cut, dump_plane_bkg_track_thetaUp_neutron]
    for i in range(len(FirstTH1)):
        FirstTH1[i].Rebin(20)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e2
    yAxisLow    = 0.5
    
    xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")


    FirstTH1    = [dump_plane_bkg_track_thetaDn_neutron_weighted_cut, dump_plane_bkg_track_thetaDn_neutron]
    for i in range(len(FirstTH1)):
        FirstTH1[i].Rebin(20)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e2
    yAxisLow    = 0.5
    
    xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")
    
    dump_plane_bkg_track_thetaUp_neutron_weighted_cut.Add(dump_plane_bkg_track_thetaDn_neutron_weighted_cut)
    dump_plane_bkg_track_thetaUp_neutron.Add(dump_plane_bkg_track_thetaDn_neutron)
    FirstTH1    = [dump_plane_bkg_track_thetaUp_neutron_weighted_cut, dump_plane_bkg_track_thetaUp_neutron]
    # for i in range(len(FirstTH1)):
    #     FirstTH1[i].Rebin(20)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e2
    yAxisLow    = 0.5
    
    xAxisTitle  = "#theta_{p} [rad]"#FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName()+"CombinedUpDn", h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")



    try:
        FirstTH1    = [dump_plane_bkg_track_thetaUp_neutron_weighted_cut_morebins, dump_plane_bkg_track_thetaUp_neutron_morebins]
        # for i in range(len(FirstTH1)):
        #     FirstTH1[i].Rebin(20)

        xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
        xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
        
        yAxisHigh   = FirstTH1[0].GetMaximum()*1e3
        yAxisLow    = 0.5
        
        xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

        h2 = FirstTH1[0].Clone("h2")
        h2.Reset()
        h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

        DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")



        FirstTH1    = [dump_plane_bkg_track_thetaDn_neutron_weighted_cut_morebins, dump_plane_bkg_track_thetaDn_neutron_morebins]
        # for i in range(len(FirstTH1)):
        #     FirstTH1[i].Rebin(20)

        xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
        xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
        
        yAxisHigh   = FirstTH1[0].GetMaximum()*1e3
        yAxisLow    = 0.5
        
        xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

        h2 = FirstTH1[0].Clone("h2")
        h2.Reset()
        h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

        DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")

        dump_plane_bkg_track_thetaUp_neutron_weighted_cut_morebins.Add(dump_plane_bkg_track_thetaDn_neutron_weighted_cut_morebins)
        dump_plane_bkg_track_thetaUp_neutron_morebins.Add(dump_plane_bkg_track_thetaDn_neutron_morebins)
        FirstTH1    = [dump_plane_bkg_track_thetaUp_neutron_weighted_cut_morebins, dump_plane_bkg_track_thetaUp_neutron_morebins]
        for i in range(len(FirstTH1)):
            FirstTH1[i].Rebin(5)

        xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
        xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
        
        yAxisHigh   = FirstTH1[0].GetMaximum()*1e3
        yAxisLow    = 0.5
        
        xAxisTitle  = "#theta_{p} [rad]"#FirstTH1[0].GetXaxis().GetTitle()

        h2 = FirstTH1[0].Clone("h2")
        h2.Reset()
        h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

        DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName()+"CombinedUpDn", h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")
    except:
        print("Something wrong in thetaUp_neutron_weighted_cut_morebins")


    
    
    
    FirstTH1    = [dump_plane_bkg_track_phi_neutron_cut, dump_plane_bkg_track_phi_neutron]
    for i in range(len(FirstTH1)):
        FirstTH1[i].Rebin(20)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e2
    yAxisLow    = 0.5
    
    xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")
    
    

    FirstTH1    = [dump_plane_bkg_track_time_neutron_cut, dump_plane_bkg_track_time_neutron]
    for i in range(len(FirstTH1)):
        FirstTH1[i].Rebin(5)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e2
    yAxisLow    = 0.5
    
    xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")



    FirstTH1    = [dump_plane_bkg_track_energy_neutron_cut, dump_plane_bkg_track_E_neutron]
    for i in range(len(FirstTH1)):
        FirstTH1[i].Rebin(5)

    xAxisLow    = FirstTH1[0].GetXaxis().GetBinCenter(1)
    xAxisHigh   = FirstTH1[0].GetXaxis().GetBinCenter(FirstTH1[0].GetNbinsX())
    
    yAxisHigh   = FirstTH1[0].GetMaximum()*1e2
    yAxisLow    = 0.5
    
    xAxisTitle  = FirstTH1[0].GetXaxis().GetTitle()

    h2 = FirstTH1[0].Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{FullSim}{FastSamp}")

    DrawHistsRatio(FirstTH1, LegendName, PlotColor, xAxisLow, xAxisHigh, yAxisLow, yAxisHigh, xAxisTitle, outDir+"/"+FirstTH1[0].GetName(), h2, 1.0, 1.0, drawline, logy, latexName, latexName2, TeVTag, doSumw2, doAtlas, doLumi, noRatio, do80, do59,"width")


if __name__=="__main__":
    main()
