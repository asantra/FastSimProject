import sys, math, array
from copy import deepcopy
import ROOT
from ROOT import TFile, TStyle, THStack, TTree, TCanvas, TPad, TROOT, TLine, TLatex, TArrow, TLegend, TColor
from ROOT import kBlack, kBlue, kRed, kWhite, kGreen, kMagenta, kOrange, kYellow, kCyan, kViolet, kAzure, kGray, kRust, kInvertedDarkBodyRadiator
from ROOT import kDot



def dPhiCalc(phiLead, phiTrail):
  dphi = abs(phiLead - phiTrail)
  if(dphi > math.pi): dphi = math.pi*2. - dphi
  return dphi

def CalculateError(deltay, y, deltax, x):
  error=0
  if (y!=0) and (x!=0):
    error = math.sqrt((deltay/y)*(deltay/y)+(deltax/x)*(deltax/x))
  return error


def MakeLine(xlow, ylow, xup, yup):
  line = TLine(xlow,ylow,xup,yup)
  line.SetLineStyle(1)
  line.SetLineWidth(2)
  return line

def MakeArrowLeft(xlow, ylow, xup, yup):
  arrow = TArrow(xlow,ylow,xup,yup, 0.01, "<|")
  #arrow.SetAngle(90)
  arrow.SetLineStyle(1)
  arrow.SetLineWidth(2)
  return arrow

def MakeArrowRight(xlow, ylow, xup, yup):
  arrow = TArrow(xlow,ylow,xup,yup, 0.01, "|>")
  #arrow.SetAngle(90)
  arrow.SetLineStyle(1)
  arrow.SetLineWidth(2)
  return arrow

def addInQuadrature(a, b):
    value1 = a*a
    value2 = b*b
    value3 = value1+value2
    if(value3>=0):
        return math.sqrt(value3)
    else:
        return -999.0

def SaveFile(c, CanvasName):
   #c.SaveAs(CanvasName+".C")
   c.SaveAs(CanvasName+".png")
   c.SaveAs(CanvasName+".pdf")
   c.SaveAs(CanvasName+".root")
   
def TexWaterMark(texString):
   tex = TLatex(-4.224138,-80.55602,texString)
   #tex.SetTextColor(17)
   tex.SetTextSize(0.199)#0.1991525
   tex.SetTextAngle(35.0) #26.15998
   tex.SetLineWidth(2)
   tex.SetTextColorAlpha(17, 1.00)
   return tex


def TexMaker(TexList, Atlas=False, Lumi=False, noRatio=False, do80Inv=False, do59Inv=False):
  if(Atlas):
        if(noRatio):
            if(Lumi):
                #TexList[0] = TLatex(0.47,0.80,"#sqrt{s}=13 TeV, 80.0 fb^{-1}") #36.075 fb^{-1} (13 TeV)
                TexList[0] = TLatex(0.47,0.80,"#sqrt{s}=13 TeV, 36.2 fb^{-1}") #36.075 fb^{-1} (13 TeV)
            else:
                #TexList[0] = TLatex(0.40,0.80,"#sqrt{s}=13 TeV, 150.0 fb^{-1}") #36.075 fb^{-1} (13 TeV)
                # TexList[0] = TLatex(0.47,0.80,"#sqrt{s}=13 TeV, 43.2 fb^{-1}")
                TexList[0] = TLatex(0.32,0.80,"#sqrt{s}=13.6 TeV")
                
            #TexList[0] = TLatex(0.47,0.80,"#sqrt{s}=13 TeV")
            TexList[0].SetNDC()
            TexList[0].SetTextAlign(31)
            TexList[0].SetTextFont(42)
            TexList[0].SetTextSize(0.037)
            TexList[0].SetLineWidth(2)
            
            TexList[1] = TLatex(0.15,0.874,"ATLAS")
            TexList[1].SetNDC()
            TexList[1].SetTextAlign(13)
            #TexList[1].SetTextFont(61) # needed for CMS
            TexList[1].SetTextFont(72) # helvetica italic
            TexList[1].SetTextSize(0.05)
            TexList[1].SetLineWidth(1)
            
            TexList[2] = TLatex(0.20,0.873,"Internal") # 
            TexList[2].SetNDC()
            TexList[2].SetTextAlign(13)
            TexList[2].SetTextFont(42)
            TexList[2].SetTextSize(0.05)
            TexList[2].SetLineWidth(1)
        else:
            if(Lumi and (not do80Inv) and (not do59Inv)):
                TexList[0] = TLatex(0.45,0.80,"#sqrt{s}=13 TeV, 36.2 fb^{-1}")
            elif((not Lumi) and (do80Inv) and (not do59Inv)):
                # TexList[0] = TLatex(0.45,0.80,"#sqrt{s}=13 TeV, 43.2 fb^{-1}")
                TexList[0] = TLatex(0.45,0.80,"#sqrt{s}=13.6 TeV")
            elif(do59Inv and ((not Lumi) and (not do80Inv))):
                TexList[0] = TLatex(0.45,0.80,"#sqrt{s}=13 TeV, 59.9 fb^{-1}")
            elif ((not Lumi) and (not do80Inv) and (not do59Inv)):
                TexList[0] = TLatex(0.45,0.80,"#sqrt{s}=13 TeV") # , 139.3 fb^{-1} ### for the plots without mean
                #TexList[0] = TLatex(0.69,0.50,"#sqrt{s}=13 TeV")
                
            TexList[0] = TLatex(0.4,0.80,"#sqrt{s}=13.6 TeV")
            TexList[0].SetNDC()
            TexList[0].SetTextAlign(31)
            TexList[0].SetTextFont(42)
            TexList[0].SetTextSize(0.037)
            TexList[0].SetLineWidth(2)
            
            TexList[1] = TLatex(0.22,0.874,"ATLAS") ### for the plots without mean
            #TexList[1] = TLatex(0.35,0.574,"ATLAS")  ### for the plots with mean
            TexList[1].SetNDC()
            TexList[1].SetTextAlign(13)
            #TexList[1].SetTextFont(61) # needed for CMS
            TexList[1].SetTextFont(72) # helvetica italic
            TexList[1].SetTextSize(0.042)
            TexList[1].SetLineWidth(1)
            
            TexList[2] = TLatex(0.35,0.873,"Internal") #Internal for the plots without mean
            #TexList[2] = TLatex(0.52,0.573,"SCT Preliminary") #Internal for the plots with mean
            TexList[2].SetNDC()
            TexList[2].SetTextAlign(13)
            TexList[2].SetTextFont(42)
            TexList[2].SetTextSize(0.037)
            TexList[2].SetLineWidth(1)
      
      
      
  else:
        TexList[0] = TLatex(0.82,0.914,"") ##sqrt{s} = 13 TeV
        TexList[0].SetNDC()
        TexList[0].SetTextAlign(31)
        TexList[0].SetTextFont(42)
        TexList[0].SetTextSize(0.039)
        TexList[0].SetLineWidth(2)
        
        TexList[1] = TLatex(0.10,0.944,"") #MoEDAL
        TexList[1].SetNDC()
        TexList[1].SetTextAlign(13)
        TexList[1].SetTextFont(61)
        TexList[1].SetTextSize(0.05)
        TexList[1].SetLineWidth(1)
        
        TexList[2] = TLatex(0.29,0.943,"") # 0.29, Simulation
        #TexList[2] = TLatex(0.12,0.944,"Simulation") ##
        TexList[2].SetNDC() 
        TexList[2].SetTextAlign(13)
        TexList[2].SetTextFont(52)
        TexList[2].SetTextSize(0.05)
        TexList[2].SetLineWidth(1)
  


  

def LegendMaker(energyPlots=False):
  # legend1 = TLegend(0.44,0.74,0.90,0.88); # original 0.5,0.55,0.9,0.9
  # legend1 = TLegend(0.18,0.76,0.90,0.88); # original 0.5,0.55,0.9,0.9
  legend1 = TLegend(0.6,0.7,0.85,0.85); ### without the mean
  # legend1 = TLegend(0.26,0.65,0.85,0.88); ### with the mean
  legend1.SetFillColor(kWhite);
  legend1.SetTextFont(42);
  if(energyPlots):
      legend1.SetTextSize(0.048);
  else:
      legend1.SetTextSize(0.05);
  legend1.SetBorderSize(0);
  legend1.SetShadowColor(kWhite);
  legend1.SetFillStyle(0);
  return legend1;


def RatioLegendMaker():
  #legend1 = TLegend(0.38,0.75,0.68,0.99); # original 0.5,0.72,0.9,0.9
  legend1 = TLegend(0.42,0.6,0.83,0.95); # original 0.5,0.72,0.9,0.9
  #legend1 = TLegend(0.52,0.7,0.93,0.9); # original 0.5,0.72,0.9,0.9
  legend1.SetFillColor(kWhite);
  legend1.SetTextFont(42);
  legend1.SetTextSize(0.12);
  legend1.SetBorderSize(0);
  legend1.SetShadowColor(kWhite);
  legend1.SetFillStyle(0);
  return legend1;


def RatioLegendMakerPull():
  #legend1 = TLegend(0.38,0.75,0.68,0.99); # original 0.5,0.72,0.9,0.9
  #legend1 = TLegend(0.2,0.32,0.83,0.55); # original 0.5,0.72,0.9,0.9
  legend1 = TLegend(0.42,0.48,0.83,0.9); # original 0.5,0.72,0.9,0.9
  legend1.SetFillColor(kWhite);
  legend1.SetTextFont(42);
  legend1.SetTextSize(0.09);
  legend1.SetBorderSize(0);
  legend1.SetShadowColor(kWhite);
  legend1.SetFillStyle(0);
  return legend1;


def LeftLegendMaker():
  legend1 = TLegend(0.24,0.82,0.50,0.87); #0.5,0.72,0.9,0.9
  legend1.SetFillColor(kWhite);
  legend1.SetTextFont(42);
  legend1.SetTextSize(0.044);
  legend1.SetBorderSize(0);
  legend1.SetShadowColor(kWhite);
  legend1.SetFillStyle(0);
  return legend1;


def LegendMakerTwoColumn():
  legend1 = TLegend(0.50,0.7,0.92,0.88)
  legend1.SetNColumns(2)
  legend1.SetLineColor(kWhite)
  legend1.SetFillColor(0)
  legend1.SetTextSize(0.0345) ## 0.027
  legend1.SetTextFont(42)
  legend1.SetBorderSize(0)
  legend1.SetShadowColor(kWhite)
  legend1.SetFillStyle(0)
  return legend1;


def SetHistColorEtc(hist, color):
     #hist = getOverflow(hist)
     hist.SetLineColor(color)
     #hist.GetYaxis().SetTitle("Events")
     hist.SetFillColor(0)
     if(color==kGreen+3):
       hist.SetMarkerStyle(kFullTriangleUp)
     elif(color==kGreen+2):
       hist.SetMarkerStyle(kFullTriangleUp)
     elif((color==2) or (color==kRed)):
       hist.SetMarkerStyle(kFullTriangleDown)
     elif ((color==4) or (color==kBlue)):
       hist.SetMarkerStyle(kFullSquare)
     elif (color==kBlack):
       hist.SetMarkerStyle(kOpenCircle)
     elif (color==kRed+3):
       hist.SetMarkerStyle(kFullDiamond)
     elif (color==kMagenta):
       hist.SetMarkerStyle(kFullCross)
     elif (color==kBlue+3):
       hist.SetMarkerStyle(kFullCircle)
     elif (color==kOrange+7):
       hist.SetMarkerStyle(kPlus)
     elif (color==kYellow):
       hist.SetMarkerStyle(kOpenStar)
     elif (color==kYellow+3):
       hist.SetMarkerStyle(kFullStar)
     elif (color==kViolet-2):
       hist.SetMarkerStyle(kFullStar)
     elif (color==kYellow-7):
       hist.SetMarkerStyle(kFullStar)
     elif (color==kBlue-9):
       hist.SetMarkerStyle(kMultiply)
     elif (color==kAzure+10):
       hist.SetMarkerStyle(kPlus)
     elif (color==kCyan-3):
       hist.SetMarkerStyle(kOpenSquare)
     elif (color==kViolet-6):
       hist.SetMarkerStyle(kOpenDiamond)
     elif (color==kGray):
       hist.SetMarkerStyle(kOpenSquare)
     elif (color==kAzure+1):
       hist.SetMarkerStyle(kMultiply)
     elif (color==kOrange+2):
       hist.SetMarkerStyle(kPlus)
     elif (color==kViolet-1):
       hist.SetMarkerStyle(kOpenDiamond)
     elif (color==kGreen-3):
       hist.SetMarkerStyle(kFullTriangleUp)
     else:
       hist.SetMarkerStyle(kDot)
       
     hist.SetMarkerColor(color)
     if (color==kBlack):
        hist.SetMarkerSize(0.5)
     else:
        hist.SetMarkerSize(0.7)
     return hist
 
 
def AxisLabelEtc(hist, yAxisName, xAxisName):
     hist.GetYaxis().SetTitle(yAxisName)
     hist.GetXaxis().SetTitle(xAxisName)
     hist.GetYaxis().SetTitleOffset(1.5) # 0.97
     hist.GetYaxis().SetLabelOffset(0.001)
     hist.GetYaxis().SetLabelSize(0.04)
     hist.GetYaxis().SetTitleSize(0.0434)
     hist.GetXaxis().SetTitleSize(0.037)
     hist.GetXaxis().SetLabelSize(0.032)
     hist.GetXaxis().SetLabelOffset(0.01)
     hist.GetXaxis().SetTitleOffset(1.05)
     hist.SetLineWidth(1)
     #hist.GetXaxis().SetNdivisions(5)
     hist.SetTitle("")
     return hist

def MakeLatex(xPoint, yPoint, latexName):
    tex = TLatex(xPoint, yPoint,latexName)
    tex.SetNDC()
    tex.SetTextAlign(31)
    tex.SetTextFont(42)
    tex.SetTextSize(0.03)
    tex.SetLineWidth(2)
    return tex



def JESError(hData, hBkg, hBkgUp, hBkgDown):
    

  den1 = hBkg.Clone ("bkgden1")
  den2 = hBkg.Clone ("bkgden2")

  ratiop = hBkg.Clone ("backgroundratiop")
  ratiom = hBkg.Clone ("backgroundratiom")

  ymin = 0.1; ymax = 2.

  x = array.array('f'); y=array.array('f'); exl=array.array('f'); eyl=array.array('f'); exh=array.array('f'); eyh=array.array('f');

  for km in range(0, hBkg.GetNbinsX()+1):
    conte1 = math.sqrt (hBkg.GetBinError (km) * hBkg.GetBinError (km) + (hBkg.GetBinContent (km) - hBkgUp.GetBinContent (km)) * (hBkg.GetBinContent (km) - hBkgUp.GetBinContent (km)))
    conte2 = math.sqrt (hBkg.GetBinError (km) * hBkg.GetBinError (km) + (hBkg.GetBinContent (km) - hBkgDown.GetBinContent (km)) * (hBkg.GetBinContent (km) - hBkgDown.GetBinContent (km)))
    
    den1.SetBinContent(km, hBkg.GetBinContent (km) + (conte1))
    den2.SetBinContent(km, hBkg.GetBinContent (km) - (conte2))


  ratiop.Divide (den1)
  ratiom.Divide (den2)
	  
  ratio = hData.Clone("ratiodata")
  ratio.Divide (hBkg)

  for km2 in range(0, ratio.GetNbinsX()+1):
      if(ratio.GetBinContent(km2) > ymax):
            ymax = ratio.GetBinContent(km2) + ratio.GetBinError(km2)
      x.append(ratio.GetBinCenter(km2))
      y.append(1.0)	
      exl.append(ratio.GetBinWidth(km2) / 2.0)
      exh.append(ratio.GetBinWidth(km2) / 2.0)

      if (ratiop.GetBinContent(km2)!= 0):
            eyh.append(1./ ratiop.GetBinContent(km2) - 1)
      else:
            eyh.append(0)

      if(ratiom.GetBinContent(km2) != 0):
            eyl.append(1 - 1. / ratiom.GetBinContent(km2))
      else:
            eyl.append(0)
	
  err = TGraphAsymmErrors(hBkg.GetNbinsX()+1, x, y, exl, exh, eyl, eyh)
  return err




def StatError(hdata, Background, Backgroundup, Backgrounddown):   

  den1 = Background.Clone ("bkgden1")
  den2 = Background.Clone ("bkgden2")
  nvar = Background.GetNbinsX()


  ratiop = Background.Clone ("backgroundratiop")
  ratiom = Background.Clone ("backgroundratiom")

  ymin = 0.1; ymax = 2.;

  x = array.array('f'); y=array.array('f'); exl=array.array('f'); eyl=array.array('f'); exh=array.array('f'); eyh=array.array('f');
  

  for km in range(0, Background.GetNbinsX()+1):
    conte1 = Background.GetBinError(km)	      
    conte2 = Background.GetBinError(km)

    den1.SetBinContent(km, Background.GetBinContent(km) + (conte1))
    den2.SetBinContent(km, Background.GetBinContent(km) - (conte2))

  ratiop.Divide (den1)
  ratiom.Divide (den2)
	  
  ratio = hdata.Clone ("ratiodata")
  ratio.Divide(Background)

  for km in range(0, ratio.GetNbinsX()+1):
    if (ratio.GetBinContent(km) > ymax):
      ymax = ratio.GetBinContent(km) + ratio.GetBinError(km)
    x.append(ratio.GetBinCenter(km))
    y.append(1.0)	
    exl.append(ratio.GetBinWidth(km)/2.0)
    exh.append(ratio.GetBinWidth(km)/2.0)

    if(ratiop.GetBinContent(km) != 0):
      eyh.append(1./ ratiop.GetBinContent(km) - 1)
    else:
      eyh.append(0)

    if(ratiom.GetBinContent(km) != 0):
      eyl.append(1 - 1. / ratiom.GetBinContent (km))
    else:
      eyl.append(0)
  
  err = TGraphAsymmErrors (nvar, x, y, exl, exh, eyl, eyh)
  return err

## needed for one ratio axis size ## 

def OneRatioAxisSize(h2, color=kBlack):
  h2.GetYaxis().SetRangeUser(0.0,4.5)
  h2.GetYaxis().SetLabelSize(0.12)
  h2.GetYaxis().SetNdivisions(4)
  h2.GetYaxis().CenterTitle()
  h2.GetYaxis().SetTitleOffset(0.45)
  h2.GetYaxis().SetTitleSize(0.1)
  
  h2.GetXaxis().SetLabelSize(0.12)
  h2.GetXaxis().SetTitleSize(0.14)
  h2.GetXaxis().SetTitleOffset(0.847)
  
  h2 = SetHistColorEtc(h2, color)
  
  h2.SetTitle('')
  return deepcopy(h2)



def getOverflow(h_Sample):
  binV = h_Sample.GetNbinsX()
  lastBinValue = h_Sample.GetBinContent(binV)
  lastBinError = h_Sample.GetBinError(binV)
  
  lastBinOverflowValue = h_Sample.GetBinContent(binV+1)
  lastBinOverflowError = h_Sample.GetBinError(binV+1)
  
  finalValue = lastBinValue + lastBinOverflowValue
  finalError = math.sqrt(lastBinError*lastBinError + lastBinOverflowError*lastBinOverflowError)
  
  h_Sample.SetBinContent(binV, finalValue)
  h_Sample.SetBinContent(binV+1, 0)
  h_Sample.SetBinError(binV, finalError)
  h_Sample.SetBinError(binV+1, 0)
  
  return h_Sample


def getUnderflow(h_Sample):
  lastBinValue = h_Sample.GetBinContent(1)
  lastBinError = h_Sample.GetBinError(1)
  
  lastBinUnderflowValue = h_Sample.GetBinContent(0)
  lastBinUnderflowError = h_Sample.GetBinError(0)
  
  finalValue = lastBinValue + lastBinUnderflowValue
  finalError = math.sqrt(lastBinError*lastBinError + lastBinUnderflowError*lastBinUnderflowError)
  
  h_Sample.SetBinContent(1, finalValue)
  h_Sample.SetBinContent(0, 0)
  h_Sample.SetBinError(1, finalError)
  h_Sample.SetBinError(0, 0)
  
  return h_Sample


def DrawHistsRatio(FirstTH1, LegendName, PlotColor, xrange1down, xrange1up, yrange1down, yrange1up, xaxisTitle, CanvasName, h2, yline1low, yline1up, drawline=False, logy=False, LatexName='', LatexName2='', TeVTag=False, doSumw2=False, doAtlas=False, doLumi=False, noRatio=False, do80=False, do59=False):
   Tex = MakeLatex(0.85,0.65,LatexName)
   Tex2 = MakeLatex(0.85,0.60,LatexName2)
   c = TCanvas("c","c",900, 900)
   c.SetGrid()
   pad1 = TPad('pad1','pad1',0,0.25,1,1)
   pad1.SetBottomMargin(0.0)
   pad1.SetLeftMargin(0.13)
   pad1.SetFillStyle(0)
   pad1.SetGrid()
   pad1.Draw()
   pad1.cd()
   ROOT.gStyle.SetOptStat(0)
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
     FirstTH1[i].SetMarkerSize(1.25)
     FirstTH1[i].GetYaxis().SetTitleOffset(1.05)
     FirstTH1[i].GetYaxis().SetTitleSize(0.05)
     FirstTH1[i].GetYaxis().SetLabelSize(0.032)
     FirstTH1[i].GetXaxis().SetTitle("")
     FirstTH1[i].GetXaxis().SetLabelSize(0.15)
    #  w = FirstTH1[i].Integral()
    #  ### comment in if you want only shape distribution
    #  FirstTH1[i].Scale(1./w)
     FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
     if(i==0):
         legend1.AddEntry(FirstTH1[i],LegendName[i], "lp")
     else:
         legend1.AddEntry(FirstTH1[i],LegendName[i], "lep")
   
     #legend1.AddEntry(FirstTH1[i],LegendName[i], "lp")
   
   FirstTH1[0].GetXaxis().SetLabelSize(0.0);
   FirstTH1[0].GetYaxis().SetTitle("Normalized entries");
   #FirstTH1[0].GetYaxis().SetTitle("Entries");
   FirstTH1[0].SetLineWidth(2)
   FirstTH1[0].Draw("hist")
   
   
    
   #WaterMark.Draw("sames")
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
   ROOT.gPad.Modified(); 
   ROOT.gPad.Update();
   
   #### special file for Allpix-squared plot, otherwise remove -1 in the range
   for i in range(1, len(FirstTH1)):
     FirstTH1[i].SetLineWidth(2)
     if ("cluster_charge" in FirstTH1[i].GetName()) or ("ClusterCharge" in FirstTH1[i].GetName()):
         FirstTH1[i].Draw("hist sames")
     else:
         FirstTH1[i].Draw("hist e sames")  
     ROOT.gPad.Modified(); ROOT.gPad.Update();
     
   if ("cluster_charge" in FirstTH1[0].GetName()) or ("ClusterCharge" in FirstTH1[0].GetName()):
        FirstTH1[0].Draw("hist sames")
   else:
       FirstTH1[0].Draw("ep sames")
   
    
   ROOT.gPad.RedrawAxis()
   L[1].Draw("sames")
   L[2].Draw("sames")
   if(TeVTag):
       TexTeV = TLatex(0.892,0.914,"#sqrt{s}=13 TeV")#36.075 fb^{-1} (13 TeV)
       TexTeV.SetNDC()
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
   pad2.SetLeftMargin(0.13)
   pad2.SetBottomMargin(0.3)
   pad2.SetFillStyle(0)
   #pad2.SetGrid()
   pad2.Draw()
   pad2.cd()
   ROOT.gStyle.SetOptStat(0)
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()

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
   h2.GetYaxis().SetTitleSize(0.14)
   h2.GetYaxis().SetTitleOffset(0.38)
   h2.GetYaxis().SetRangeUser(0.4, 1.6)
   h2.Draw("ep")
   if(drawline):
     line.SetLineStyle(2)
     line.Draw()
     
   SaveFile(c, CanvasName)
   return deepcopy(c)
   #c.Clear()
     
     
     
     
     
     
     
     
def DrawHistsRatioTF1(FirstTH1, LegendName, PlotColor, xrange1down, xrange1up, yrange1down, yrange1up, xaxisTitle, CanvasName, h2, yline1low, yline1up, drawline=False, logy=False, LatexName='', LatexName2='', TeVTag=False, doSumw2=False, doAtlas=False, doLumi=False, noRatio=False, do80=False, do59=False):
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
   ROOT.gStyle.SetOptStat(0)
   if(logy):
     pad1.SetLogy()
     
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
     FirstTH1[i].GetYaxis().SetTitle(FirstTH1[i].GetYaxis().GetTitle())#(FirstTH1[i].GetYaxis().GetTitle())
     
     FirstTH1[i].GetYaxis().SetTitleOffset(0.95)
     FirstTH1[i].GetYaxis().SetTitleSize(0.05)
     FirstTH1[i].GetYaxis().SetLabelSize(0.045)
     FirstTH1[i].GetXaxis().SetTitle("")
     FirstTH1[i].GetXaxis().SetLabelSize(0.15)
     #w = FirstTH1[i].Integral()
     ### comment in if you want only shape distribution
     #FirstTH1[i].Scale(1.0/w)
     FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
     legend1.AddEntry(FirstTH1[i],LegendName[i], "l")
   
     #legend1.AddEntry(FirstTH1[i],LegendName[i], "lp")
   
   FirstTH1[0].GetXaxis().SetLabelSize(0.0);
   FirstTH1[0].GetYaxis().SetTitle("Entries");
   FirstTH1[0].SetLineWidth(1)
   FirstTH1[0].Draw("l")
   
    
    
   #WaterMark.Draw("sames")
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
   ROOT.gPad.Modified(); ROOT.gPad.Update();
   
   #### special file for Allpix-squared plot, otherwise remove -1 in the range
   for i in range(1, len(FirstTH1)):
     FirstTH1[i].SetLineWidth(1)
     FirstTH1[i].Draw("l sames")  
     ROOT.gPad.Modified(); ROOT.gPad.Update();
     
   ROOT.gPad.RedrawAxis()
   L[1].Draw("sames")
   L[2].Draw("sames")
   if(TeVTag):
       TexTeV = TLatex(0.892,0.914,"#sqrt{s}=13 TeV")#36.075 fb^{-1} (13 TeV)
       TexTeV.SetNDC()
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
   pad2.SetLogx()
   #pad2.SetGrid()
   pad2.Draw()
   pad2.cd()
   ROOT.gStyle.SetOptStat(0)
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
   
   if(doSumw2):
       for i in range(0, len(FirstTH1)):
          FirstTH1[i].Sumw2() 
   
   h2 = FirstTH1[0].Clone('h2')
   h2.Divide(FirstTH1[1])
   h2.GetXaxis().SetRangeUser(xrange1down,xrange1up)
   h2.GetXaxis().SetTitle(xaxisTitle)
   h2.GetYaxis().SetTitle('dist/fit')
   h2.GetYaxis().CenterTitle()
   h2.SetTitle('')
   h2 = OneRatioAxisSize(h2)
   h2.GetYaxis().SetRangeUser(0.85, 1.15)
   h2.Draw("l")
   if(drawline):
     line.SetLineStyle(2)
     line.Draw()
     
   #SaveFile(c, CanvasName)
   return deepcopy(c)
   #c.Clear()
     
     
     
     



### draw hists with systematics
def DrawHistsSystRatio(FirstTH1, LegendName, PlotColor, xrange1down, xrange1up, yrange1down, yrange1up, xaxisTitle, CanvasName, h2, hSyst, yline1low, yline1up, drawline=False, logy=False, LatexName='', LatexName2='', TeVTag=False, doSumw2=False, doAtlas=False, doLumi=False, noRatio=False, do80=False, do59=False):
   Tex = MakeLatex(0.68,0.65,LatexName)
   Tex2 = MakeLatex(0.68,0.60,LatexName2)
   c = TCanvas("c","c",900, 900)
   c.SetGrid()
   pad1 = TPad('pad1','pad1',0,0.25,1,1)
   pad1.SetBottomMargin(0.0)
   pad1.SetFillStyle(0)
   pad1.SetGrid()
   pad1.Draw()
   pad1.cd()
   ROOT.gStyle.SetOptStat(0)
   if(logy):
     pad1.SetLogy()
     
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
     FirstTH1[i].GetYaxis().SetTitle(FirstTH1[i].GetYaxis().GetTitle())#(FirstTH1[i].GetYaxis().GetTitle())
     FirstTH1[i].GetYaxis().SetTitleOffset(0.95)
     FirstTH1[i].GetYaxis().SetTitleSize(0.05)
     FirstTH1[i].GetYaxis().SetLabelSize(0.045)
     FirstTH1[i].GetXaxis().SetTitle("")
     FirstTH1[i].GetXaxis().SetLabelSize(0.15)
     #w = FirstTH1[i].Integral()
     ### comment in if you want only shape distribution
     #FirstTH1[i].Scale(1.0/w)
     FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
     legend1.AddEntry(FirstTH1[i],LegendName[i], "lp")
   
   FirstTH1[0].GetXaxis().SetLabelSize(0.0);
   FirstTH1[0].Draw("ep")
   #WaterMark.Draw("sames")
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
   ROOT.gPad.Modified(); ROOT.gPad.Update();
   for i in range(1, len(FirstTH1)):
     FirstTH1[i].Draw("ep sames")  
     ROOT.gPad.Modified(); ROOT.gPad.Update();
   FirstTH1[0].Draw("ep sames")
   ### making the systematic error histogram
   hFinalSyst = hSyst.Clone('hFinalSyst')
   hFinalSyst.Reset()
   for i in range(0, hSyst.GetNbinsX()+1):
       hFinalSyst.SetBinContent(i, FirstTH1[1].GetBinContent(i))
       hFinalSyst.SetBinError(i, math.sqrt(hSyst.GetBinError(i)*hSyst.GetBinError(i)+FirstTH1[1].GetBinError(i)*FirstTH1[1].GetBinError(i)))
       hFinalSyst.SetLineColor(kGray)
       hFinalSyst.SetFillStyle(1001)
       hFinalSyst.SetFillColor(kGray)
       hFinalSyst.SetMarkerColor(kGray)
       
   legend1.AddEntry(hFinalSyst, "stat #oplus sys", "F")
   hFinalSyst.Draw("e2 sames")
   ROOT.gPad.RedrawAxis()
   L[1].Draw("sames")
   L[2].Draw("sames")
   if(TeVTag):
       TexTeV = TLatex(0.892,0.914,"#sqrt{s}=13 TeV")#36.075 fb^{-1} (13 TeV)
       TexTeV.SetNDC()
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
   ROOT.gStyle.SetOptStat(0)
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
   
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
   h2.GetYaxis().SetRangeUser(0.1, 1.45)
   h2.Draw("ep")
   
   hFinalSyst.Divide(hFinalSyst, FirstTH1[1])
   for i in range(0, hFinalSyst.GetNbinsX()+1):
       hFinalSyst.GetBinContent(i,1)
   hFinalSyst.Draw("e2 sames")
   h2.Draw("ep sames")
   ROOT.gPad.RedrawAxis();
   ROOT.gPad.Update();
   c.cd()
   
   if(drawline):
     line.Draw()
     
   SaveFile(c, CanvasName)
   c.Clear()
  
     
     
     
     
#### Draw hists, then create ratio of first histogram to all other histograms     
def DrawHistsRatioMany(FirstTH1, LegendName, PlotColor, xrange1down, xrange1up, yrange1down, yrange1up, xaxisTitle, CanvasName, h2, yline1low, yline1up, drawline=False, logy=False, LatexName='', LatexName2='', TeVTag=False, doSumw2=False, doAtlas=False, doLumi=False, noRatio=False, do80=False, do59=False):
   Tex = MakeLatex(0.68,0.65,LatexName)
   Tex2 = MakeLatex(0.68,0.60,LatexName2)
   c = TCanvas("c","c",900, 900)
   c.SetGrid()
   pad1 = TPad('pad1','pad1',0,0.25,1,1)
   pad1.SetBottomMargin(0.0)
   pad1.SetFillStyle(0)
   pad1.SetGrid()
   pad1.Draw()
   pad1.cd()
   ROOT.gStyle.SetOptStat(0)
   if(logy):
     pad1.SetLogy()
     
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
     FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
     legend1.AddEntry(FirstTH1[i],LegendName[i], "lp")
   
   FirstTH1[0].GetXaxis().SetLabelSize(0.0);
   FirstTH1[0].Draw("ep")
   #WaterMark.Draw("sames")
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
   ROOT.gPad.Modified(); ROOT.gPad.Update();
   for i in range(1, len(FirstTH1)):
     FirstTH1[i].Draw("ep sames")  
     ROOT.gPad.Modified(); ROOT.gPad.Update();
   FirstTH1[0].Draw("ep sames")
   ROOT.gPad.RedrawAxis()
   L[1].Draw("sames")
   L[2].Draw("sames")
   
   if(TeVTag):
       TexTeV = TLatex(0.892,0.914,"#sqrt{s}=13 TeV")
       TexTeV.SetNDC()
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
   ROOT.gStyle.SetOptStat(0)
   legend2 = RatioLegendMaker()
   
   
   if(doSumw2):
       for i in range(0, len(FirstTH1)):
          FirstTH1[i].Sumw2() 
   
   ratioHistList = []
   
   for i in range(1, len(FirstTH1)):
       nameMain = 'h'+str(i)+'Main'
       hM = FirstTH1[i].Clone(nameMain)
       hM.Reset()
       hM.Divide(FirstTH1[0], FirstTH1[i])
       hM.GetXaxis().SetRangeUser(xrange1down,xrange1up)
       hM.GetXaxis().SetTitle(xaxisTitle)
       hM.GetYaxis().SetTitle(h2.GetYaxis().GetTitle())
       hM.GetYaxis().CenterTitle()
       hM.SetTitle('')
       hM = OneRatioAxisSize(hM, PlotColor[i])
       legend2.AddEntry(hM, LegendName[0]+' / '+LegendName[i],"lp")
       hM.GetYaxis().SetRangeUser(0.0, 2.5)
       ratioHistList.append(hM)
       
       
   
   ratioHistList[0].Draw("ep")
   for i in range(1, len(ratioHistList)):
       ratioHistList[i].Draw("ep sames")
       ROOT.gPad.Modified(); ROOT.gPad.Update();
   ratioHistList[0].Draw("ep sames")
   legend2.Draw("sames")
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
   ROOT.gPad.RedrawAxis()
   if(drawline):
     line.Draw()
     
   SaveFile(c, CanvasName)
   c.Clear()
   
     
def DrawHists(FirstTH1, LegendName, PlotColor,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, yline1low, yline1up, drawline=False, logy=False, latexName='', latexName2 = '', latexName3='', leftLegend=False, doAtlas=False, doLumi=False, noRatio=False, do80=False, do59=False, drawPattern="", logz=False, logx=False, latexName4='', zAxisName='Entries per primary electron',zrange1down=0,zrange1up=0.07):
   debug = False
   if(debug): print("just entering plot code")
   ### without mean
   Tex  = MakeLatex(0.80,0.70,latexName)
   Tex2 = MakeLatex(0.80,0.64,latexName2)
   Tex3 = MakeLatex(0.80,0.58,latexName3)
   Tex4 = MakeLatex(0.34,0.52,latexName4)
   
   
   if(debug): print ("defining Tex ")
   WaterMark = TexWaterMark('Preliminary')
   c = TCanvas("c","c",500, 500)
   ROOT.gStyle.SetOptStat(0)
   #### for TDR 2D plots 
   #ROOT.gPad.SetLeftMargin(0.1)
   #ROOT.gPad.SetRightMargin(0.1)
   #ROOT.gPad.SetBottomMargin(0.14)
   ### for others
   ROOT.gPad.SetLeftMargin(0.15)
   ROOT.gPad.SetRightMargin(0.14)
   ROOT.gPad.SetBottomMargin(0.14)
   c.cd()
   c.SetGrid()
   
   if(logx):
     c.SetLogx()
   if(logy):
     c.SetLogy()
   if(logz):
     c.SetLogz()
     
   if(debug): print ("Set Logy ")
   line = MakeLine(xrange1down,yline1low,xrange1up,yline1up)
   if(leftLegend):
       legend1 = LeftLegendMaker()
   else:
       legend1 = LegendMaker()
       #legend1 = LegendMakerTwoColumn()
   if(debug): print ("Set Legend ")
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L, doAtlas, doLumi, noRatio, do80, do59)
   stList = []
   #
   if(debug): print ("Set Ranges ")
   integralList = []; integralListError = []
   strType = str(type(FirstTH1[0]))
   
  #  #### turn on palette
  #  if 'TH2' in strType:
  #      ROOT.gStyle.SetPalette(kVisibleSpectrum)
   
   for i in range(0, len(FirstTH1)):
     #if('residual' in FirstTH1[i].GetName()): FirstTH1[i].Rebin(4)
     
     FirstTH1[i] = AxisLabelEtc(FirstTH1[i], yAxisName, xAxisName)
     if("TH2" not in strType):
        #FirstTH1[i].Rebin(4)
        FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
        # if("hEnergy_Electron" in FirstTH1[i].GetName() and i<2): FirstTH1[i].SetLineStyle(4+i)
     else:
        #FirstTH1[i].Rebin2D(2)
        FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
        FirstTH1[i].SetLineColor(kWhite)
     #FirstTH1[i] = getOverflow(FirstTH1[i])
       
     #FirstTH1[0].SetFillColor(PlotColor[0])
     FirstTH1[i].SetFillColor(0)
     if("TH2D" in strType):
        xBinMax = FirstTH1[i].GetNbinsX()
        yBinMax = FirstTH1[i].GetNbinsY()
        integralList.append(FirstTH1[i].Integral(0, xBinMax+1, 0, yBinMax+1))
     elif("TH1" in strType):
         w = FirstTH1[i].Integral(0, FirstTH1[i].GetNbinsX()+1)
         integralList.append(w)
         # comment in if you want only shape distribution
         #if(w!=0):
            #FirstTH1[i].Scale(1.0/w)
     else:
        
        #integralList.append(FirstTH1[i].GetMean())
        #integralListError.append(FirstTH1[i].GetMeanError())
        pass
     #else:
     FirstTH1[i].GetYaxis().SetRangeUser(yrange1down,yrange1up)
     FirstTH1[i].GetXaxis().SetRangeUser(xrange1down,xrange1up)
     if zAxisName:
        # FirstTH1[i].GetZaxis().SetTitle("Particles/pixels per BX")
        FirstTH1[i].GetZaxis().SetTitle("Events")
     else:
        FirstTH1[i].GetZaxis().SetTitle(zAxisName)
     
       
     FirstTH1[i].SetMaximum(yrange1up)
     FirstTH1[i].SetMinimum(yrange1down)
     #legend1.AddEntry(FirstTH1[i],LegendName[i]+" ("+str(round(integralList[i],1))+")", "l")
     #legend1.AddEntry(FirstTH1[i],LegendName[i]+" (mean: "+str(round(integralList[i],4))+"#pm"+str(round(integralListError[i],4))+")", "lp")
     if("TH1D" in strType):legend1.AddEntry(FirstTH1[i],LegendName[i]+" ("+str(round(integralList[i],1))+")", "l")
     else: legend1.AddEntry(FirstTH1[i],LegendName[i]+" ("+str(round(integralList[i],1))+")", "l")
     
   
   
   if(debug): print ("After for loop ")
   FirstTH1[0].GetXaxis().SetRangeUser(xrange1down,xrange1up)
   FirstTH1[0].GetXaxis().SetNdivisions(5)
   #FirstTH1[0].GetYaxis().SetRangeUser(yrange1down,yrange1up)
   if "electrons" in FirstTH1[0].GetName():
       FirstTH1[0].GetZaxis().SetRangeUser(0,10)
   elif "positrons" in FirstTH1[0].GetName():
       FirstTH1[0].GetZaxis().SetRangeUser(0,5)
   elif "gamma" in FirstTH1[0].GetName():
       FirstTH1[0].GetZaxis().SetRangeUser(0,100)
   else:
       FirstTH1[0].GetZaxis().SetRangeUser(zrange1down,zrange1up)
       
   FirstTH1[0].GetXaxis().SetNdivisions(9)
   #FirstTH1[0].GetXaxis().CenterLabels()
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
   
   #FirstTH1[0].SetFillColor(0)
   
   if "TH2" in strType:
       drawStyle = "COLZ" #drawPattern
   elif "TH1" in strType:
       drawStyle = "hist"
   else:
       drawStyle = "AL"
   
   #ROOT.gPad.SetRightMargin(0.08)
   #ROOT.gPad.SetBottomMargin(1.0)
   FirstTH1[0].Draw(drawStyle) # ce, hist
    
   
   ROOT.gPad.Update();
   ROOT.gPad.RedrawAxis();
   if "TH2" in strType:
        FirstTH1[0].GetZaxis().SetLabelSize(0.02)
        FirstTH1[0].GetZaxis().SetTitleOffset(1.4)
        FirstTH1[0].GetYaxis().SetTitleOffset(1.1)
        FirstTH1[0].GetYaxis().SetLabelSize(0.03)
        #FirstTH1[0].GetZaxis().SetRangeUser(0, 1e3)
        try:
            pl = FirstTH1[0].GetListOfFunctions().FindObject("palette")
            pl.SetX1NDC(0.86);
            pl.SetX2NDC(0.88);
            pl.SetY1NDC(0.14);
            pl.SetY2NDC(0.90);
            # pl.SetLabelSize(0.005);
            ROOT.gPad.Modified();
            ROOT.gPad.Update();
        except:
            print("This histogram is empty: ", FirstTH1[0].GetName())
   
   if(debug): 
       print ("After first Draw ")
       print("len(FirstTH1): ", len(FirstTH1), " FirstTH1: ", FirstTH1)

   #WaterMark.Draw("sames")
   if(len(FirstTH1)>1):
    #### special file for Allpix-squared plot, otherwise remove -1 in the range
    for i in range(1, len(FirstTH1)):
        if(debug): print("drawing i=", i, " FirstTH1: ", FirstTH1[i])
        FirstTH1[i].Draw(drawStyle+" sames") 
        FirstTH1[i].SetFillColor(0)
        
   #FirstTH1[0].Draw("hist sames")
   if(debug): print ("After Draw loop ")
   
   Tex.Draw("sames")
   Tex2.Draw("sames")
   Tex3.Draw("sames")
   Tex4.Draw("sames")
   #if(drawline):
     #line.Draw()
   L[0].Draw()
   L[1].Draw()
   L[2].Draw()
   if "TH1" in strType: legend1.Draw()
   else: pass #legend1.Draw()
   
   
   
   
   
   #st = FirstTH1[0].GetListOfFunctions().FindObject("stats")
   #st.SetX1NDC(0.70); #new x start position
   #st.SetX2NDC(0.90); #new x end position
   #st.SetLineColor(kBlue);
   #st.SetLineColor(PlotColor[0])
   #stList.append(st)
   #for m in range(0, len(stList)):
     #x1NDC = 0.25+0.25*m
     #x2NDC = 0.45+0.25*m
     #stList[m].SetX1NDC(x1NDC); #new x start position
     #stList[m].SetX2NDC(x2NDC);
   
   
   SaveFile(c, CanvasName)
   return deepcopy(c)
   #return [c,L,legend1]




def DrawRatioHists(FirstTH1, LegendName, PlotColor,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, yline1low, yline1up, drawline=False, logy=False, latexName='', latexName2 = '', latexName3='', leftLegend=False, doAtlas=False, doLumi=False, noRatio=False, do80=False, do59=False, drawPattern="", logz=False, logx=False, latexName4=''):
   
   debug = False
   if(debug): print("just entering plot code")
   ### without mean
   Tex = MakeLatex(0.55,0.72,latexName)
   Tex2 = MakeLatex(0.55,0.66,latexName2)
   Tex3 = MakeLatex(0.55,0.60,latexName3)
   Tex4 = MakeLatex(0.24,0.52,latexName4)
   
   
   if(debug): print ("defining Tex ")
   WaterMark = TexWaterMark('Preliminary')
   c = TCanvas("c","c",500, 500)
   ROOT.gStyle.SetOptStat(0)
   ROOT.gPad.SetLeftMargin(0.15)
   ROOT.gPad.SetRightMargin(0.15)
   ROOT.gPad.SetBottomMargin(0.10)
   c.cd()
   c.SetGrid()
   
   if(logx):
     c.SetLogx()
   if(logy):
     c.SetLogy()
   if(logz):
     c.SetLogz()
     
   if(debug): print ("Set Logy ")
   line = MakeLine(xrange1down,yline1low,xrange1up,yline1up)
   if(leftLegend):
       legend1 = LeftLegendMaker()
   else:
       legend1 = LegendMaker()

   if(debug): print ("Set Legend ")
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L, doAtlas, doLumi, noRatio, do80, do59)
   stList = []
   
   if(debug): print ("Set Ranges ")
   integralList = []; integralListError = []
   strType = str(type(FirstTH1[0]))
   
   for i in range(0, len(FirstTH1)):
     
     FirstTH1[i] = AxisLabelEtc(FirstTH1[i], yAxisName, xAxisName)
     if("TH2" not in strType):
        FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
     else:
        FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
       
     FirstTH1[i].SetFillColor(0)
     if("TH2D" in strType):
        xBinMax = FirstTH1[i].GetNbinsX()
        yBinMax = FirstTH1[i].GetNbinsY()
        integralList.append(FirstTH1[i].Integral(0, xBinMax+1, 0, yBinMax+1))
     elif("TH1" in strType):
         w = FirstTH1[i].Integral(0, FirstTH1[i].GetNbinsX()+1)
         integralList.append(w)
     else:
        pass

     FirstTH1[i].GetYaxis().SetRangeUser(yrange1down,yrange1up)
     FirstTH1[i].GetXaxis().SetRangeUser(xrange1down,xrange1up)
     FirstTH1[i].GetZaxis().SetTitle("Particles/BX")
     FirstTH1[i].GetYaxis().SetTitle("Ratio")
     
       
     FirstTH1[i].SetMaximum(yrange1up)
     FirstTH1[i].SetMinimum(yrange1down)
     #legend1.AddEntry(FirstTH1[i],LegendName[i]+" ("+str(round(integralList[i],1))+")", "l")
     #legend1.AddEntry(FirstTH1[i],LegendName[i]+" (mean: "+str(round(integralList[i],4))+"#pm"+str(round(integralListError[i],4))+")", "lp")
     
     #legend1.AddEntry(FirstTH1[i],LegendName[i], "l")
   
   
   if("TH1D" in strType):
       legend1.AddEntry(FirstTH1[0],LegendName[0], "l")
       if(len(FirstTH1) > 2):
           legend1.AddEntry(FirstTH1[2],LegendName[2], "l")
   else: 
       legend1.AddEntry(FirstTH1[0],LegendName[0], "l")
   
   if(debug): print ("After for loop ")
   
   FirstTH1[0].GetXaxis().SetRangeUser(xrange1down,xrange1up)
   FirstTH1[0].GetXaxis().SetNdivisions(5)
   
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
   
   
   if "TH1" in strType:
       drawStyle = "e"
   else:
       drawStyle = "AL"
   
   
   FirstTH1[0].Divide(FirstTH1[1])
   FirstTH1[0].Draw(drawStyle) # ce, hist
   if(len(FirstTH1)>3):
       FirstTH1[2].Divide(FirstTH1[3])
       FirstTH1[2].Draw(drawStyle+" same") # ce, hist
   
   ROOT.gPad.Update();
   ROOT.gPad.RedrawAxis();
   if "TH2" in strType:
        FirstTH1[0].GetZaxis().SetLabelSize(0.023)
        FirstTH1[0].GetZaxis().SetTitleOffset(1.4)
        try:
            pl = FirstTH1[0].GetListOfFunctions().FindObject("palette")
            pl.SetX1NDC(0.86);
            pl.SetX2NDC(0.88);
            pl.SetY1NDC(0.10);
            pl.SetY2NDC(0.90);
            pl.SetLabelSize(0.005);
            ROOT.gPad.Modified();
            ROOT.gPad.Update();
        except:
            print("This histogram is empty: ", FirstTH1[0].GetName())
   
   if(debug): print ("After first Draw ")
           
   Tex.Draw("sames")
   Tex2.Draw("sames")
   Tex3.Draw("sames")
   Tex4.Draw("sames")
   
   L[0].Draw()
   L[1].Draw()
   L[2].Draw()
   legend1.Draw()

   SaveFile(c, CanvasName)
   return [c,L,legend1]




def DrawHistsCutValueOneArrow(FirstTH1, LegendName, PlotColor,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, yline1low, yline1up, drawline=False, logy=False, latexName='', latexName2 = '', latexName3='', leftLegend=False, doAtlas=False, doLumi=False, noRatio=False, do80=False, do59=False, drawPattern="", logz=False, logx=False, latexName4='', xline1down=0.0, xline1up=0.0):
   
   debug = False
   if(debug): print ("just entering plot code")
   ### without mean
   Tex = MakeLatex(0.37,0.84,latexName)
   Tex2 = MakeLatex(0.37,0.78,latexName2)
   Tex3 = MakeLatex(0.37,0.72,latexName3)
   Tex4 = MakeLatex(0.37,0.66,latexName4)
   
   
   if(debug): print ("defining Tex ")
   WaterMark = TexWaterMark('Preliminary')
   c = TCanvas("c","c",700, 700)
   ROOT.gStyle.SetOptStat(0)
   c.cd()
   c.SetGrid()
   if(logx):
     c.SetLogx()
   if(logy):
     c.SetLogy()
   if(logz):
       c.SetLogz()
   if(debug): print ("Set Logy ")
   line  = MakeLine(xline1down,yline1low,xline1up,yline1up)
   arrow = MakeArrowLeft(xline1up-0.02*xrange1up,yline1up,xline1up,yline1up)
   if(leftLegend):
       legend1 = LeftLegendMaker()
   else:
       legend1 = LegendMaker()
       #legend1 = LegendMakerTwoColumn()
   if(debug): print ("Set Legend ")
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L, doAtlas, doLumi, noRatio, do80, do59)
   stList = []
   #
   if(debug): print ("Set Ranges ")
   integralList = []; integralListError = []
   strType = str(type(FirstTH1[0]))
   for i in range(0, len(FirstTH1)):
     if('Par2' in FirstTH1[i].GetName()):
        FirstTH1[i].Rebin(8)
     else:
        FirstTH1[i].Rebin(4)
     FirstTH1[i].SetLineWidth(1)
     FirstTH1[i] = AxisLabelEtc(FirstTH1[i], yAxisName, xAxisName)
     if("TH1" in strType):
        FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
     FirstTH1[i] = getOverflow(FirstTH1[i])
     FirstTH1[i].GetYaxis().SetRangeUser(yrange1down,yrange1up)
     FirstTH1[i].GetXaxis().SetRangeUser(xrange1down,xrange1up)
     #w = FirstTH1[i].Integral(0, FirstTH1[i].GetNbinsX()+1)
     ### comment in if you want only shape distribution
     #if(w!=0):
       #FirstTH1[i].Scale(1.0/w)
       
     if("TH1" in strType):
        integralList.append(FirstTH1[i].Integral(0, FirstTH1[i].GetNbinsX()+1))
        #integralList.append(FirstTH1[i].GetMean())
        #integralListError.append(FirstTH1[i].GetMeanError())
     else:
        xBinMax = FirstTH1[i].GetNbinsX()
        yBinMax = FirstTH1[i].GetNbinsY()
        integralList.append(FirstTH1[i].Integral(0, xBinMax+1, 0, yBinMax+1))
     
       
     #FirstTH1[0].SetMaximum(yrange1up)
     #FirstTH1[0].SetMinimum(yrange1down)
     #legend1.AddEntry(FirstTH1[i],LegendName[i]+" ("+str(round(integralList[i],1))+" particles)", "hist")
     #legend1.AddEntry(FirstTH1[i],LegendName[i]+" (mean: "+str(round(integralList[i],4))+"#pm"+str(round(integralListError[i],4))+")", "lp")
     if(i==0):legend1.AddEntry(FirstTH1[i],LegendName[i]+" ("+str(round(integralList[i],2))+")", "l")
     else: legend1.AddEntry(FirstTH1[i],LegendName[i]+" ("+str(round(integralList[i],2))+")", "l")
   
   if(debug): print ("After for loop ")
   FirstTH1[0].GetXaxis().SetRangeUser(xrange1down,xrange1up)
   FirstTH1[0].SetFillColor(0)
   #FirstTH1[0].GetYaxis().SetRangeUser(yrange1down,yrange1up)
   if "electrons" in FirstTH1[0].GetName():
       FirstTH1[0].GetZaxis().SetRangeUser(0,40)
   
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
   
   
   if "TH2" in strType:
       drawStyle = drawPattern
   else:
       drawStyle = "hist"
   
   FirstTH1[0].Draw(drawStyle) # ce, hist
   
   if(debug): print ("After first Draw ")
   #WaterMark.Draw("sames")
   if(len(FirstTH1)>1):
    for i in range(1, len(FirstTH1)):
        FirstTH1[i].Draw(drawStyle+" sames") 
        FirstTH1[i].SetFillColor(0)
        FirstTH1[i].SetLineWidth(2)
        
   if(debug): print ("After Draw loop ")
   
   Tex.Draw("sames")
   Tex2.Draw("sames")
   Tex3.Draw("sames")
   Tex4.Draw("sames")
   
   L[0].Draw()
   L[1].Draw()
   L[2].Draw()
   legend1.Draw()
   
   line.Draw()
   arrow.Draw()
   
   SaveFile(c, CanvasName)
   return [c,L,legend1]









def DrawHistsCutValueTwoArrows(FirstTH1, LegendName, PlotColor,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, yline1low, yline1up, drawline=False, logy=False, latexName='', latexName2 = '', latexName3='', leftLegend=False, doAtlas=False, doLumi=False, noRatio=False, do80=False, do59=False, drawPattern="", logz=False, logx=False, latexName4='', xline1down=0.0, xline1up=0.0, negativexline1down=-1.0, negativexline1up=-1.0):
   
   debug = False
   if(debug): print ("just entering plot code")
   ### without mean
   Tex = MakeLatex(0.37,0.84,latexName)
   Tex2 = MakeLatex(0.37,0.78,latexName2)
   Tex3 = MakeLatex(0.37,0.72,latexName3)
   Tex4 = MakeLatex(0.37,0.66,latexName4)
   
   
   if(debug): print ("defining Tex ")
   WaterMark = TexWaterMark('Preliminary')
   c = TCanvas("c","c",700, 700)
   ROOT.gStyle.SetOptStat(0)
   c.cd()
   c.SetGrid()
   if(logx):
     c.SetLogx()
   if(logy):
     c.SetLogy()
   if(logz):
       c.SetLogz()
   if(debug): print ("Set Logy ")
   line1  = MakeLine(xline1down,yline1low,xline1up,yline1up)
   arrow1 = MakeArrowLeft(xline1up-0.002,yline1up,xline1up,yline1up)
   line2  = MakeLine(negativexline1down,yline1low,negativexline1up,yline1up)
   arrow2 = MakeArrowRight(negativexline1up,yline1up,negativexline1up+0.002,yline1up)
   if(leftLegend):
       legend1 = LeftLegendMaker()
   else:
       legend1 = LegendMaker()
       #legend1 = LegendMakerTwoColumn()
   if(debug): print ("Set Legend ")
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L, doAtlas, doLumi, noRatio, do80, do59)
   stList = []
   #
   if(debug): print ("Set Ranges ")
   integralList = []; integralListError = []
   strType = str(type(FirstTH1[0]))
   for i in range(0, len(FirstTH1)):
     FirstTH1[i].Rebin(4)
     FirstTH1[i].SetLineWidth(1)
     FirstTH1[i] = AxisLabelEtc(FirstTH1[i], yAxisName, xAxisName)
     if("TH1" in strType):
        FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
     FirstTH1[i] = getOverflow(FirstTH1[i])
     FirstTH1[i].GetYaxis().SetRangeUser(yrange1down,yrange1up)
     FirstTH1[i].GetXaxis().SetRangeUser(xrange1down,xrange1up)
     #w = FirstTH1[i].Integral(0, FirstTH1[i].GetNbinsX()+1)
     ### comment in if you want only shape distribution
     #if(w!=0):
       #FirstTH1[i].Scale(1.0/w)
       
     if("TH1" in strType):
        integralList.append(FirstTH1[i].Integral(0, FirstTH1[i].GetNbinsX()+1))
        #integralList.append(FirstTH1[i].GetMean())
        #integralListError.append(FirstTH1[i].GetMeanError())
     else:
        xBinMax = FirstTH1[i].GetNbinsX()
        yBinMax = FirstTH1[i].GetNbinsY()
        integralList.append(FirstTH1[i].Integral(0, xBinMax+1, 0, yBinMax+1))
     
       
     #FirstTH1[0].SetMaximum(yrange1up)
     #FirstTH1[0].SetMinimum(yrange1down)
     #legend1.AddEntry(FirstTH1[i],LegendName[i]+" ("+str(round(integralList[i],1))+" particles)", "hist")
     #legend1.AddEntry(FirstTH1[i],LegendName[i]+" (mean: "+str(round(integralList[i],4))+"#pm"+str(round(integralListError[i],4))+")", "lp")
     if(i==0):legend1.AddEntry(FirstTH1[i],LegendName[i]+" ("+str(round(integralList[i],2))+")", "l")
     else: legend1.AddEntry(FirstTH1[i],LegendName[i]+" ("+str(round(integralList[i],2))+")", "l")
   
   if(debug): print ("After for loop ")
   FirstTH1[0].GetXaxis().SetRangeUser(xrange1down,xrange1up)
   FirstTH1[0].SetFillColor(0)
   #FirstTH1[0].GetYaxis().SetRangeUser(yrange1down,yrange1up)
   if "electrons" in FirstTH1[0].GetName():
       FirstTH1[0].GetZaxis().SetRangeUser(0,40)
   
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
   
   
   if "TH2" in strType:
       drawStyle = drawPattern
   else:
       drawStyle = "hist"
   
   FirstTH1[0].Draw(drawStyle) # ce, hist
   
   if(debug): print ("After first Draw ")
   #WaterMark.Draw("sames")
   if(len(FirstTH1)>1):
    for i in range(1, len(FirstTH1)):
        FirstTH1[i].Draw(drawStyle+" sames") 
        FirstTH1[i].SetFillColor(0)
        FirstTH1[i].SetLineWidth(2)
        
   if(debug): print ("After Draw loop ")
   
   Tex.Draw("sames")
   Tex2.Draw("sames")
   Tex3.Draw("sames")
   Tex4.Draw("sames")
   
   L[0].Draw()
   L[1].Draw()
   L[2].Draw()
   legend1.Draw()
   
   line1.Draw()
   arrow1.Draw()
   line2.Draw()
   arrow2.Draw()
   
   SaveFile(c, CanvasName)
   return [c,L,legend1]









#### This is a special function to draw one canvas with 8 different histograms
def DrawHistsOneCanvas(FirstTH1, LegendName, PlotColor,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, yline1low, yline1up, drawline=False, logy=False, latexName='', latexName2 = '', latexName3='', leftLegend=False, doAtlas=False, doLumi=False, noRatio=False, do80=False, do59=False, drawPattern="", logz=False, logx=False):
   ##### with mean
   #Tex = MakeLatex(0.40,0.79,latexName)
   #Tex2 = MakeLatex(0.40,0.73,latexName2)
   #Tex3 = MakeLatex(0.40,0.67,latexName3)
   debug = False
   if(debug): print ("just entering plot code")
   histName = FirstTH1[0].GetName()
   ### without mean
   Tex = MakeLatex(0.58,0.68,latexName)
   Tex2 = MakeLatex(0.58,0.62,latexName2)
   Tex3 = MakeLatex(0.58,0.56,latexName3)
   if(debug): print ("defining Tex ")
   WaterMark = TexWaterMark('Preliminary')
   c = TCanvas("c","c",1000, 700)
   ROOT.gStyle.SetOptStat(0)
   c.cd()
   c.SetGrid()
   if(logx):
     c.SetLogx()
   if(logy):
     c.SetLogy()
   if(logz):
       c.SetLogz()
   if(debug): print ("Set Logy ")
   line = MakeLine(xrange1down,yline1low,xrange1up,yline1up)
   if(leftLegend):
       #legend1 = LeftLegendMaker()
       legend1 = LegendMakerTwoColumn()
   else:
       #legend1 = LegendMaker()
       legend1 = LegendMakerTwoColumn()
   if(debug): print ("Set Legend ")
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L, doAtlas, doLumi, noRatio, do80, do59)
   stList = []
   #FirstTH1[0].GetYaxis().SetRangeUser(yrange1down,yrange1up)
   if(debug): print ("Set Ranges ")
   integralList = []; integralListError = []
   strType = str(type(FirstTH1[0]))
   for i in range(0, len(FirstTH1)):
     FirstTH1[i].Rebin(4)
     FirstTH1[i].SetLineWidth(2)
     FirstTH1[i] = AxisLabelEtc(FirstTH1[i], yAxisName, xAxisName)
     FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
     FirstTH1[i] = getOverflow(FirstTH1[i])
     FirstTH1[i].GetXaxis().SetRangeUser(xrange1down,xrange1up)
     FirstTH1[i].GetYaxis().SetRangeUser(yrange1down,yrange1up)
     if(i==3 or i==7):
         FirstTH1[i].SetFillColor(PlotColor[i])
         #FirstTH1[i].SetFillStyle(4030)
     else:    
         FirstTH1[i].SetFillColor(0)
     
     if(i>3):
         FirstTH1[i].SetLineStyle(7)
     else:
         FirstTH1[i].SetLineStyle(1)
     w = FirstTH1[i].Integral(0, FirstTH1[i].GetNbinsX()+1)
     FirstTH1[i].Scale(w)
     ### comment in if you want only shape distribution
     #if(w!=0):
       #FirstTH1[i].Scale(1.0/w)
       
       
     if("TH1" in strType):
        integralList.append(FirstTH1[i].Integral(0, FirstTH1[i].GetNbinsX()+1))
        #integralList.append(FirstTH1[i].GetMean())
        #integralListError.append(FirstTH1[i].GetMeanError())
     else:
        xBinMax = FirstTH1[i].GetNbinsX()
        yBinMax = FirstTH1[i].GetNbinsY()
        integralList.append(FirstTH1[i].Integral(0, xBinMax+1, 0, yBinMax+1))
     
     if(("TH1" in strType) and ("sumE" not in histName)):
        legend1.AddEntry(FirstTH1[i],LegendName[i]+" ("+str(round(integralList[i],1))+" particles)", "lp")
     else:
        legend1.AddEntry(FirstTH1[i],LegendName[i], "lp")
     #legend1.AddEntry(FirstTH1[i],LegendName[i]+" (mean: "+str(round(integralList[i],4))+"#pm"+str(round(integralListError[i],4))+")", "lp")
     #
   
   if(debug): print ("After for loop ")
   
   #if "electrons" in FirstTH1[0].GetName():
       #FirstTH1[0].GetZaxis().SetRangeUser(0,40)
   #else:
       #FirstTH1[0].GetZaxis().SetRangeUser(0,6000)
   #FirstTH1[0].GetXaxis().SetNdivisions(10)
   #FirstTH1[0].GetXaxis().CenterLabels()
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
   
   if "TH2" in strType:
       drawStyle = drawPattern
   else:
        drawStyle = "hist"
   
   FirstTH1[3].Draw(drawStyle)
   FirstTH1[7].Draw(drawStyle+" same")
   FirstTH1[0].Draw(drawStyle+" same")
   FirstTH1[1].Draw(drawStyle+" same")
   FirstTH1[2].Draw(drawStyle+" same")
   FirstTH1[4].Draw(drawStyle+" same")
   FirstTH1[5].Draw(drawStyle+" same")
   FirstTH1[6].Draw(drawStyle+" same")

   #for i in range(0,len(FirstTH1)):
        #if(i!=3 or i!=7):
             ## ce, hist
            #FirstTH1[i].Draw(drawStyle+" same")
        #else:
            #FirstTH1[i].Draw("hist sames")
   
   
   #ROOT.gPad.Update();
   if "TH2" in strType:
    pl = FirstTH1[0].GetListOfFunctions().FindObject("palette")
    pl.SetX1NDC(0.905);
    pl.SetX2NDC(0.919);
    pl.SetY1NDC(0.10);
    pl.SetY2NDC(0.90);
    #ROOT.gPad.Update();
   if(debug): print ("After first Draw ")
   #WaterMark.Draw("sames")
        
   if(debug): print ("After Draw loop ")
   
   Tex.Draw("sames")
   Tex2.Draw("sames")
   Tex3.Draw("sames")
   #if(drawline):
     #line.Draw()
   L[0].Draw()
   L[1].Draw()
   L[2].Draw()
   legend1.Draw()
   
   
   #st = FirstTH1[0].GetListOfFunctions().FindObject("stats")
   #st.SetX1NDC(0.70); #new x start position
   #st.SetX2NDC(0.90); #new x end position
   #st.SetLineColor(kBlue);
   #st.SetLineColor(PlotColor[0])
   #stList.append(st)
   #for m in range(0, len(stList)):
     #x1NDC = 0.25+0.25*m
     #x2NDC = 0.45+0.25*m
     #stList[m].SetX1NDC(x1NDC); #new x start position
     #stList[m].SetX2NDC(x2NDC);
   
   
   SaveFile(c, CanvasName)
   return [c,L,legend1]






##### This is a special function to draw 16 histograms in 4 canvases
def DrawHists4Canvas(FirstTH1, LegendName, PlotColor, xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, yline1low, yline1up, drawline=False, logy=False, latexName='', latexName2 = '', latexName3='', leftLegend=False, doAtlas=False, doLumi=False, noRatio=False, do80=False, do59=False):
   #Tex = MakeLatex(0.74,0.54,latexName)
   #Tex2 = MakeLatex(0.74,0.44,latexName2)
   #Tex3 = MakeLatex(0.74,0.48,latexName3)
   debug = False
   if(debug): print ("just entering plot code")
   Tex = MakeLatex(0.34,0.74,latexName)
   Tex2 = MakeLatex(0.34,0.68,latexName2)
   Tex3 = MakeLatex(0.34,0.62,latexName3)
   if(debug): print ("defining Tex ")
   WaterMark = TexWaterMark('Preliminary')
   c = TCanvas("c","c",1300, 1300)
   c.Divide(1,4)
   ROOT.gStyle.SetOptStat(0)
   c.SetGrid()
   
   if(debug): print ("Set Logy ")
   
   line = MakeLine(xrange1down,yline1low,xrange1up,yline1up)
   
   legend = []
   if(leftLegend):
       for i in range(0,4):
           l = LeftLegendMaker()
           legend.append(l)
   else:
       for i in range(0,4):
           l = LegendMaker()
           legend.append(l)
       #legend1 = LegendMakerTwoColumn()
       
   if(debug): print ("Set Legend ")
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L, doAtlas, doLumi, noRatio, do80, do59)
   stList = []
   #
   if(debug): print ("Set Ranges ")
   
   integralList = []
   strType = str(type(FirstTH1[0]))
   for i in range(0, len(FirstTH1)):
     FirstTH1[i].SetLineWidth(1)
     FirstTH1[i] = AxisLabelEtc(FirstTH1[i], yAxisName, xAxisName)
     FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
     FirstTH1[i] = getOverflow(FirstTH1[i])
     FirstTH1[i].SetFillColor(0)
     FirstTH1[i].GetYaxis().SetRangeUser(yrange1down,yrange1up)
     FirstTH1[i].GetXaxis().SetRangeUser(xrange1down,xrange1up)
     if("TH1" in strType):
        integralList.append(FirstTH1[i].Integral(0, FirstTH1[i].GetNbinsX()+1))
     else:
        xBinMax = FirstTH1[i].GetNbinsX()
        yBinMax = FirstTH1[i].GetNbinsY()
        integralList.append(FirstTH1[i].Integral(0, xBinMax+1, 0, yBinMax+1))
     #w = FirstTH1[i].Integral(0, FirstTH1[i].GetNbinsX()+1)
     ### comment in if you want only shape distribution
     #if(w!=0):
       #FirstTH1[i].Scale(1.0/w)
       
    
   ### set the y axis
   #FirstTH1[0].SetMaximum(yrange1up)
   #FirstTH1[0].SetMinimum(yrange1down)
     
   ### adding special order in legend, [0,1,8,9], [2,3,10,11], [4,5,12,13], [6,7,14,15]
   for leg in range(0,8,2):
       legend[leg/2].AddEntry(FirstTH1[leg],LegendName[leg]+" ("+str(round(integralList[leg],1))+" particles)", "lp")
       legend[leg/2].AddEntry(FirstTH1[leg+1],LegendName[leg+1]+" ("+str(round(integralList[leg+1],1))+" particles)", "lp")
       legend[leg/2].AddEntry(FirstTH1[leg+8],LegendName[leg+8]+" ("+str(round(integralList[leg+8],1))+" particles)", "lp")
       legend[leg/2].AddEntry(FirstTH1[leg+9],LegendName[leg+9]+" ("+str(round(integralList[leg+9],1))+" particles)", "lp")
   
   if(debug): print ("After for loop ")
   FirstTH1[0].GetXaxis().SetRangeUser(xrange1down,xrange1up)
   #FirstTH1[0].GetXaxis().SetNdivisions(8)
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
           
   for canvasId in range(0,8,2):
        c.cd(canvasId/2+1)
        if(logy):
            ROOT.gPad.SetLogy(1)
        else:
            ROOT.gPad.SetLogy(0)
        FirstTH1[canvasId].Draw("hist")
        FirstTH1[canvasId+1].Draw("hist sames")
        FirstTH1[canvasId+8].Draw("hist sames")
        FirstTH1[canvasId+9].Draw("hist sames")
        legend[canvasId/2].Draw()
        Tex.Draw("sames")
        Tex2.Draw("sames")
        Tex3.Draw("sames")
        c.Update()
        
   if(debug): print ("After Draw loop ")
   
   
   #if(drawline):
     #line.Draw()
   L[0].Draw()
   L[1].Draw()
   L[2].Draw()
   #legend1.Draw()
   
   #st = FirstTH1[0].GetListOfFunctions().FindObject("stats")
   #st.SetX1NDC(0.70); #new x start position
   #st.SetX2NDC(0.90); #new x end position
   #st.SetLineColor(kBlue);
   #st.SetLineColor(PlotColor[0])
   #stList.append(st)
   #for m in range(0, len(stList)):
     #x1NDC = 0.25+0.25*m
     #x2NDC = 0.45+0.25*m
     #stList[m].SetX1NDC(x1NDC); #new x start position
     #stList[m].SetX2NDC(x2NDC);
   
   
   SaveFile(c, CanvasName)
   return [c,L,legend]










##### This is a special function to draw 8 histograms in 8 canvases
def DrawHists8Canvas(FirstTH1, LegendName, PlotColor, xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, yline1low, yline1up, drawline=False, logy=False, latexName='', latexName2 = '', latexName3='', leftLegend=False, doAtlas=False, doLumi=False, noRatio=False, do80=False, do59=False, drawPattern="", logz=False):
   debug = False
   if(debug): print ("just entering plot code")
   Tex = MakeLatex(0.64,0.74,latexName)
   Tex2 = MakeLatex(0.64,0.68,latexName2)
   Tex3 = MakeLatex(0.64,0.62,latexName2)
   if(debug): print ("defining Tex ")
   WaterMark = TexWaterMark('Preliminary')
   c = TCanvas("c","c",1300, 1300)
   c.Divide(2,4)
   ROOT.gStyle.SetOptStat(0)
   c.SetGrid()
   
   if(debug): print ("Set Logy ")
   
   line = MakeLine(xrange1down,yline1low,xrange1up,yline1up)
   
   legend = []
   if(leftLegend):
       for i in range(0,len(FirstTH1)):
           l = LeftLegendMaker()
           legend.append(l)
   else:
       for i in range(0,len(FirstTH1)):
           l = LegendMaker()
           legend.append(l)
       #legend1 = LegendMakerTwoColumn()
       
   if(debug): print ("Set Legend ")
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L, doAtlas, doLumi, noRatio, do80, do59)
   stList = []
   #
   if(debug): print ("Set Ranges ")
   
   integralList = []
   strType = str(type(FirstTH1[0]))
   for i in range(0, len(FirstTH1)):
     FirstTH1[i].SetLineWidth(1)
     FirstTH1[i] = AxisLabelEtc(FirstTH1[i], yAxisName, xAxisName)
     FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[0])
     FirstTH1[i] = getOverflow(FirstTH1[i])
     FirstTH1[i] = getUnderflow(FirstTH1[i])
     FirstTH1[i].SetFillColor(0)
     if("TH1" in strType):
        integralList.append(FirstTH1[i].Integral(0, FirstTH1[i].GetNbinsX()+1))
     else:
        xBinMax = FirstTH1[i].GetNbinsX()
        yBinMax = FirstTH1[i].GetNbinsY()
        integralList.append(FirstTH1[i].Integral(0, xBinMax+1, 0, yBinMax+1))
        
     FirstTH1[i].GetYaxis().SetRangeUser(yrange1down,yrange1up)
     FirstTH1[i].GetXaxis().SetRangeUser(xrange1down,xrange1up)
     #if("E [GeV]" in xAxisName):
        #FirstTH1[i].GetXaxis().SetRangeUser(xrange1down,xrange1up)
    
   #### set the y axis
   #FirstTH1[0].SetMaximum(yrange1up)
   #FirstTH1[0].SetMinimum(yrange1down)
     
   ### adding legend, each canvas gets one in the same order
   for leg in range(0,len(FirstTH1)):
       legend[leg].AddEntry(FirstTH1[leg],LegendName[leg]+" ("+str(round(integralList[leg],3))+" particles)", "hist")
       
   
   if(debug): print ("After for loop ")
   
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
           
   for canvasId in range(0,8):
        c.cd(canvasId+1)
        if(logy):
            ROOT.gPad.SetLogy(1)
        else:
            ROOT.gPad.SetLogy(0)
        
        drawStyle="hist"
        if ("TH2" in strType) and (drawPattern=="COLZ"):
            drawStyle = "COLZ"
            if(logz):
                ROOT.gPad.SetLogz(1)
            else:
                ROOT.gPad.SetLogz(0)
            
        else:
            drawStyle = "hist"
        
        #print "For canvasName : ", CanvasName," and Id", canvasId, " the integral is: ", integralList[canvasId], " and latexName: ", integralLatexName
        
        FirstTH1[canvasId].Draw(drawStyle)
        legend[canvasId].Draw()
            
        Tex.Draw("sames")
        Tex2.Draw("sames")
        c.Update()
        
   if(debug): print ("After Draw loop ")
   
   
   #if(drawline):
     #line.Draw()
   L[0].Draw()
   L[1].Draw()
   L[2].Draw()
   
   SaveFile(c, CanvasName)
   return [c,L,legend]


##### This is a special function to draw 8 histograms in 8 canvases
def DrawHists2Canvas(FirstTH1, LegendName, PlotColor, xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, yline1low, yline1up, drawline=False, logy=False, latexName='', latexName2 = '', latexName3='', leftLegend=False, doAtlas=False, doLumi=False, noRatio=False, do80=False, do59=False, drawPattern="", logz=False):
   debug = False
   if(debug): print ("just entering plot code")
   Tex = MakeLatex(0.64,0.84,latexName)
   Tex2 = MakeLatex(0.64,0.68,latexName2)
   Tex3 = MakeLatex(0.64,0.62,latexName2)
   if(debug): print ("defining Tex ")
   WaterMark = TexWaterMark('Preliminary')
   c = TCanvas("c","c", 1500, 900)
   c.Divide(2,1)
   ROOT.gStyle.SetOptStat(0)
   c.SetGrid()
   
   if(debug): print ("Set Logy ")
   
   line = MakeLine(xrange1down,yline1low,xrange1up,yline1up)
   
   legend = []
   if(leftLegend):
       for i in range(0,len(FirstTH1)):
           l = LeftLegendMaker()
           legend.append(l)
   else:
       for i in range(0,len(FirstTH1)):
           l = LegendMaker()
           legend.append(l)
       #legend1 = LegendMakerTwoColumn()
       
   if(debug): print ("Set Legend ")
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L, doAtlas, doLumi, noRatio, do80, do59)
   stList = []
   #
   if(debug): print ("Set Ranges ")
   
   integralList = []
   strType = str(type(FirstTH1[0]))
   for i in range(0, len(FirstTH1)):
     FirstTH1[i].SetLineWidth(1)
     FirstTH1[i] = AxisLabelEtc(FirstTH1[i], yAxisName, xAxisName)
     FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[0])
     FirstTH1[i] = getOverflow(FirstTH1[i])
     FirstTH1[i] = getUnderflow(FirstTH1[i])
     FirstTH1[i].SetFillColor(0)
     if("TH1" in strType):
        integralList.append(FirstTH1[i].Integral(0, FirstTH1[i].GetNbinsX()+1))
     else:
        xBinMax = FirstTH1[i].GetNbinsX()
        yBinMax = FirstTH1[i].GetNbinsY()
        integralList.append(FirstTH1[i].Integral(0, xBinMax+1, 0, yBinMax+1))
        
     FirstTH1[i].GetYaxis().SetRangeUser(yrange1down,yrange1up)
     FirstTH1[i].GetXaxis().SetRangeUser(xrange1down,xrange1up)
     #if("E [GeV]" in xAxisName):
        #FirstTH1[i].GetXaxis().SetRangeUser(xrange1down,xrange1up)
    
   #### set the y axis
   #FirstTH1[0].SetMaximum(yrange1up)
   #FirstTH1[0].SetMinimum(yrange1down)
     
   ### adding legend, each canvas gets one in the same order
   for leg in range(0,len(FirstTH1)):
      #  legend[leg].AddEntry(FirstTH1[leg],LegendName[leg]+" ("+str(round(integralList[leg],3))+" particles)", "hist")
       legend[leg].AddEntry(FirstTH1[leg],LegendName[leg], "hist")
       
   
   if(debug): print ("After for loop ")
           
   for canvasId in range(0,2):
        c.cd(canvasId+1)
        ROOT.gPad.SetLeftMargin(0.15)
        ROOT.gPad.SetRightMargin(0.15)
        
        if(logy):
            ROOT.gPad.SetLogy(1)
        else:
            ROOT.gPad.SetLogy(0)
        
        drawStyle="hist"
        if ("TH2" in strType) and (drawPattern=="COLZ"):
            drawStyle = "COLZ"
            if(logz):
                ROOT.gPad.SetLogz(1)
            else:
                ROOT.gPad.SetLogz(0)
            
        else:
            drawStyle = "hist"
        
        #print "For canvasName : ", CanvasName," and Id", canvasId, " the integral is: ", integralList[canvasId], " and latexName: ", integralLatexName
        
        FirstTH1[canvasId].Draw(drawStyle)
        legend[canvasId].Draw()
        ROOT.gPad.SetTickx()
        ROOT.gPad.SetTicky()
            
        Tex.Draw("sames")
        Tex2.Draw("sames")
        c.Update()
        
   if(debug): print ("After Draw loop ")
   
   
   #if(drawline):
     #line.Draw()
   L[0].Draw()
   L[1].Draw()
   L[2].Draw()
   
   SaveFile(c, CanvasName)
   return [c,L,legend]






##### This is a special function to draw 8*n histograms in 8 canvases, right now written for n<=4
def DrawHists8CanvasManyFiles(FirstTH1, LegendName1, PlotColor1, xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, yline1low, yline1up, drawline=False, logy=False, latexName='', latexName2 = '', latexName3='', leftLegend=False, doAtlas=False, doLumi=False, noRatio=False, do80=False, do59=False, drawPattern="", logz=False, withSignal=False):
   print ("========== ", FirstTH1[0].GetName(), " =======")
   
   debug = False
   if(debug): print ("just entering plot code")
   Tex = MakeLatex(0.54,0.74,latexName)
   Tex2 = MakeLatex(0.54,0.68,latexName2)
   Tex3 = MakeLatex(0.54,0.62,latexName2)
   if(debug): print ("defining Tex ")
   WaterMark = TexWaterMark('Preliminary')
   c = TCanvas("c","c",1300, 1300)
   c.Divide(2,4)
   ROOT.gStyle.SetOptStat(0)
   c.SetGrid()
   
   if(debug): print ("Set Logy ")
   
   line = MakeLine(xrange1down,yline1low,xrange1up,yline1up)
   
   legend = []
   if(leftLegend):
       for i in range(0,8):
           l = LeftLegendMaker()
           legend.append(l)
   else:
       for i in range(0,8):
           l = LegendMaker()
           legend.append(l)
       #legend1 = LegendMakerTwoColumn()
       
   if(debug): print ("Set Legend ")
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L, doAtlas, doLumi, noRatio, do80, do59)
   stList = []
   #
   if(debug): print ("Set Ranges ")
   
   integralList1 = []
   strType = str(type(FirstTH1[0]))
   histName = FirstTH1[0].GetName()
   for i in range(0, len(FirstTH1)):
     FirstTH1[i].SetLineWidth(1)
     FirstTH1[i] = AxisLabelEtc(FirstTH1[i], yAxisName, xAxisName)
     FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor1[i])
     FirstTH1[i] = getOverflow(FirstTH1[i])
     FirstTH1[i] = getUnderflow(FirstTH1[i])
     if(i>23):FirstTH1[i].SetFillColor(PlotColor1[i])
     else:    FirstTH1[i].SetFillColor(0)
     if("TH1" in strType):
        integralList1.append(FirstTH1[i].Integral(0, FirstTH1[i].GetNbinsX()+1))
     else:
        xBinMax = FirstTH1[i].GetNbinsX()
        yBinMax = FirstTH1[i].GetNbinsY()
        integralList1.append(FirstTH1[i].Integral(0, xBinMax+1, 0, yBinMax+1))
        
     FirstTH1[i].GetYaxis().SetRangeUser(yrange1down,yrange1up)
     FirstTH1[i].GetXaxis().SetRangeUser(xrange1down,xrange1up)
     if("E [GeV]" in xAxisName):
        FirstTH1[i].GetXaxis().SetRangeUser(xrange1down,xrange1up)
    
   ### set the y axis
   #FirstTH1[0].SetMaximum(yrange1up)
   #FirstTH1[0].SetMinimum(yrange1down)
     
   ### adding legend, each canvas gets one in the same order
   for leg in range(0,8):
       if(("TH1" in strType) and ("sumE" not in histName)): 
           if(len(FirstTH1)>0):
                legend[leg].AddEntry(FirstTH1[leg],LegendName1[leg]+" ("+str(round(integralList1[leg],3))+" particles)", "lp")
           if(len(FirstTH1)>8):
                legend[leg].AddEntry(FirstTH1[leg+8],LegendName1[leg+8]+" ("+str(round(integralList1[leg+8],3))+" particles)", "lp")
           if(len(FirstTH1)>16):
                legend[leg].AddEntry(FirstTH1[leg+16],LegendName1[leg+16]+" ("+str(round(integralList1[leg+16],3))+" particles)", "lp")
           if(len(FirstTH1)>24):
                legend[leg].AddEntry(FirstTH1[leg+24],LegendName1[leg+24]+" ("+str(round(integralList1[leg+24],3))+" particles)", "lp")
       else: 
           if(len(FirstTH1)>0):
                legend[leg].AddEntry(FirstTH1[leg],LegendName1[leg], "lp")
           if(len(FirstTH1)>8):
                legend[leg].AddEntry(FirstTH1[leg+8],LegendName1[leg+8], "lp")
           if(len(FirstTH1)>16):
                legend[leg].AddEntry(FirstTH1[leg+16],LegendName1[leg+16], "lp")
           if(len(FirstTH1)>24):
                legend[leg].AddEntry(FirstTH1[leg+24],LegendName1[leg+24], "lp")
           
   
   if(debug): print ("After for loop ")
   
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
           
   for canvasId in range(0,8):
        c.cd(canvasId+1)
        if(logy):
            ROOT.gPad.SetLogy(1)
        else:
            ROOT.gPad.SetLogy(0)
        
        drawStyle="hist"
        if ("TH2" in strType) and (drawPattern=="COLZ"):
            drawStyle = "COLZ"
            if(logz):
                ROOT.gPad.SetLogz(1)
            else:
                ROOT.gPad.SetLogz(0)
        else:
            drawStyle = "hist"
        
        #print "For canvasName : ", CanvasName," and Id", canvasId, " the integral is: ", integralList[canvasId], " and latexName: ", integralLatexName
        
        ### signal is fourth entry for now
        if(withSignal):
            if(len(FirstTH1)>24):
                FirstTH1[canvasId+24].Draw(drawStyle)
                FirstTH1[canvasId].Draw(drawStyle+" same")
                FirstTH1[canvasId+8].Draw(drawStyle+" same")
                FirstTH1[canvasId+16].Draw(drawStyle+" same")
            else:
                print( "Need a signal in the fourth set. Exiting")
                exit()
        else:
            if(len(FirstTH1)>0):
                FirstTH1[canvasId].Draw(drawStyle)
            if(len(FirstTH1)>8):
                FirstTH1[canvasId+8].Draw(drawStyle+" same")
            if(len(FirstTH1)>16):
                FirstTH1[canvasId+16].Draw(drawStyle+" same")
            if(len(FirstTH1)>24):
                FirstTH1[canvasId+24].Draw(drawStyle+" same")
        
        legend[canvasId].Draw("same")
            
        Tex.Draw("sames")
        Tex2.Draw("sames")
        c.Update()
        
   if(debug): print ("After Draw loop ")
   
   
   #if(drawline):
     #line.Draw()
   L[0].Draw()
   L[1].Draw()
   L[2].Draw()
   
   SaveFile(c, CanvasName)
   return [c,L,legend]





##### This is a special function to draw 9 histograms in 9 canvases
def DrawHists9Canvas(FirstTH1, LegendName, PlotColor, xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, yline1low, yline1up, drawline=False, logy=False, latexName=[], latexName2 = '', latexName3='', leftLegend=False, doAtlas=False, doLumi=False, noRatio=False, do80=False, do59=False, drawPattern="COLZ", logz=False):
   debug = False
   if(debug): print ("just entering plot code")
   Tex2 = MakeLatex(0.64,0.68,latexName2)
   Tex3 = MakeLatex(0.64,0.65,latexName3)
   if(debug): print("defining Tex ")
   WaterMark = TexWaterMark('Preliminary')
   c = TCanvas("c","c",300,300)
   c.Divide(3,3)
   ROOT.gStyle.SetOptStat(0)
   ROOT.gStyle.SetPalette(kInvertedDarkBodyRadiator)
   c.SetGrid()
   
   if(debug): print ("Set Logy ")
   
   line = MakeLine(xrange1down,yline1low,xrange1up,yline1up)
   
   legend = []
   if(leftLegend):
       for i in range(0,len(FirstTH1)):
           l = LeftLegendMaker()
           legend.append(l)
   else:
       for i in range(0,len(FirstTH1)):
           l = LegendMaker()
           legend.append(l)
       #legend1 = LegendMakerTwoColumn()
       
   if(debug): print ("Set Legend ")
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L, doAtlas, doLumi, noRatio, do80, do59)
   stList = []
   
   if(debug): print ("Set Ranges ")
   
   integralList = []
   texList      = []
   strType = str(type(FirstTH1[0]))
   
   skippedIndex = -999
   for i in range(0, len(FirstTH1)):
     Tex = MakeLatex(0.78,0.70,latexName[i])
     texList.append(Tex)
     checkstrType = str(type(FirstTH1[i]))
     if "TH" not in checkstrType:
        skippedIndex = i
        continue
     FirstTH1[i].SetLineWidth(1)
     FirstTH1[i] = AxisLabelEtc(FirstTH1[i], yAxisName, xAxisName)
     FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[0])
     FirstTH1[i] = getOverflow(FirstTH1[i])
     FirstTH1[i] = getUnderflow(FirstTH1[i])
     #FirstTH1[i].Rebin(4)
     #FirstTH1[i].SetFillColor(0)
     if("TH1" in strType):
        integralList.append(FirstTH1[i].Integral(0, FirstTH1[i].GetNbinsX()+1))
     else:
        # FirstTH1[i].GetZaxis().SetTitle("Number of pixels (after digitization)")
        FirstTH1[i].GetZaxis().SetTitle("Cluster position (after digitization)")
        # FirstTH1[i].GetZaxis().SetTitle("Number of hits (E_dep > 0)")
        zMax = FirstTH1[i].GetMaximum()
        # FirstTH1[i].GetZaxis().SetRangeUser(0,1.5)
        xBinMax = FirstTH1[i].GetNbinsX()
        yBinMax = FirstTH1[i].GetNbinsY()
        # integralList.append(FirstTH1[i].Integral(0, xBinMax+1, 0, yBinMax+1))
        integralList.append(FirstTH1[i].GetEntries())
        
     FirstTH1[i].GetYaxis().SetRangeUser(yrange1down,yrange1up)
     FirstTH1[i].GetXaxis().SetRangeUser(xrange1down,xrange1up)
     #if("E [GeV]" in xAxisName):
        #FirstTH1[i].GetXaxis().SetRangeUser(xrange1down,xrange1up)
    
   #### set the y axis
   #FirstTH1[0].SetMaximum(0)
   #FirstTH1[0].SetMinimum(0.4)
     
   ### adding legend, each canvas gets one in the same order
   for leg in range(0,len(FirstTH1)):
       if leg!=skippedIndex:
        legend[leg].AddEntry(FirstTH1[leg],LegendName[leg]+" ("+str(round(integralList[leg],0))+")", "hist")
       
   
   if(debug): print ("After for loop ")
   
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
           
   
   for canvasId in range(0,9):
      if canvasId!=skippedIndex:
        c.cd(canvasId+1)
        if(logy):
            ROOT.gPad.SetLogy(1)
        else:
            ROOT.gPad.SetLogy(0)
        
        ROOT.gPad.SetLeftMargin(0.15)
        ROOT.gPad.SetRightMargin(0.13)
        ROOT.gPad.SetBottomMargin(0.1)
        drawStyle="hist"
        if ("TH2" in strType) and (drawPattern=="COLZ"):
            drawStyle = "COLZ"
            FirstTH1[canvasId].GetZaxis().SetTitleOffset(1.3)
            if(logz):
                ROOT.gPad.SetLogz(1)
            else:
                ROOT.gPad.SetLogz(0)
            
        else:
            drawStyle = "hist"
        
        #print "For canvasName : ", CanvasName," and Id", canvasId, " the integral is: ", integralList[canvasId], " and latexName: ", integralLatexName
        
        FirstTH1[canvasId].Draw(drawStyle)
        
        ROOT.gPad.Update();
        
        legend[canvasId].Draw("sames")
        texList[canvasId].Draw("sames")
        
            
        #Tex.Draw("sames")
        Tex2.Draw("sames")
        Tex3.Draw("sames")
        ROOT.gPad.Update();
        if "TH2" in strType:
            pl = FirstTH1[canvasId].GetListOfFunctions().FindObject("palette")
            pl.SetX1NDC(0.885);
            pl.SetX2NDC(0.905);
            pl.SetY1NDC(0.10);
            pl.SetY2NDC(0.90);
            
        c.Update()
        
   if(debug): print ("After Draw loop ")
   
   
   #if(drawline):
     #line.Draw()
   L[0].Draw()
   L[1].Draw()
   L[2].Draw()
   
   if(debug): print ("After tex drawing ")
   
   SaveFile(c, CanvasName)
   if(debug): print ("After saving canvas ")
   return [c,L,legend]




##### This is a special function to draw 9 histograms in 9 canvases, signal and background together
def DrawHists9CanvasOverlay(FirstTH1, SecondTH1, LegendName, SecondLegendName, PlotColor, SecondPlotColor, xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, yline1low, yline1up, drawline=False, logy=False, latexName=[], latexName2 = '', latexName3='', leftLegend=False, doAtlas=False, doLumi=False, noRatio=False, do80=False, do59=False, drawPattern="", logz=False, logx=False):
   
   debug = False
   if(debug): print ("just entering plot code")
   Tex2 = MakeLatex(0.86,0.52,latexName2)
   Tex2.SetTextSize(0.04)
   Tex3 = MakeLatex(0.76,0.16,latexName3)
   Tex3.SetTextSize(0.035)
   
   if(debug): print("defining Tex ")
   WaterMark = TexWaterMark('Preliminary')
   c = TCanvas("c","c",500,500)
   c.Divide(3,3)
   ROOT.gStyle.SetOptStat(0)
   c.SetGrid()
   
   if(debug): print ("Set Logy ")
   
   line = MakeLine(xrange1down,yline1low,xrange1up,yline1up)
   
   legend = []
   if(leftLegend):
       for i in range(0,len(FirstTH1)):
           l = LeftLegendMaker()
           legend.append(l)
   else:
       for i in range(0,len(FirstTH1)):
           l = LegendMaker()
           legend.append(l)
       #legend1 = LegendMakerTwoColumn()
       
   if(debug): print ("Set Legend ")
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L, doAtlas, doLumi, noRatio, do80, do59)
   stList = []
   
   if(debug): print ("Set Ranges ")
   
   integralList = []
   texList      = []
   strType = str(type(FirstTH1[0]))
   for i in range(0, len(FirstTH1)):
     Tex = MakeLatex(0.39,0.80,latexName[i])
     Tex.SetTextSize(0.05)
     texList.append(Tex)
     FirstTH1[i].SetLineWidth(2)
     FirstTH1[i] = AxisLabelEtc(FirstTH1[i], yAxisName, xAxisName)
     FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[0])
     FirstTH1[i] = getOverflow(FirstTH1[i])
     FirstTH1[i] = getUnderflow(FirstTH1[i])
     #FirstTH1[i].Rebin(10)
     FirstTH1[i].SetFillColor(0)
     if("TH1" in strType):
        integralList.append(FirstTH1[i].Integral(0, FirstTH1[i].GetNbinsX()+1))
     else:
        xBinMax = FirstTH1[i].GetNbinsX()
        yBinMax = FirstTH1[i].GetNbinsY()
        integralList.append(FirstTH1[i].Integral(0, xBinMax+1, 0, yBinMax+1))
        
     FirstTH1[i].GetYaxis().SetRangeUser(yrange1down,yrange1up)
     FirstTH1[i].GetXaxis().SetRangeUser(xrange1down,xrange1up)
     
   ### adding legend, each canvas gets one in the same order
   for leg in range(0,len(FirstTH1)):
       legend[leg].AddEntry(FirstTH1[leg],LegendName[leg]+" ("+str(round(integralList[leg],0))+")", "l")
       
       
   secondIntegralList = []
   for i in range(0, len(SecondTH1)):
     SecondTH1[i].SetLineWidth(2)
     SecondTH1[i] = AxisLabelEtc(SecondTH1[i], yAxisName, xAxisName)
     SecondTH1[i] = SetHistColorEtc(SecondTH1[i], SecondPlotColor[0])
     SecondTH1[i] = getOverflow(SecondTH1[i])
     SecondTH1[i] = getUnderflow(SecondTH1[i])
     #SecondTH1[i].Rebin(10)
     SecondTH1[i].SetFillColor(0)
     if("TH1" in strType):
        secondIntegralList.append(SecondTH1[i].Integral(0, SecondTH1[i].GetNbinsX()+1))
     else:
        xBinMax = SecondTH1[i].GetNbinsX()
        yBinMax = SecondTH1[i].GetNbinsY()
        secondIntegralList.append(SecondTH1[i].Integral(0, xBinMax+1, 0, yBinMax+1))
        
     SecondTH1[i].GetYaxis().SetRangeUser(yrange1down,yrange1up)
     SecondTH1[i].GetXaxis().SetRangeUser(xrange1down,xrange1up)
     
   ### adding legend, each canvas gets one in the same order
   for leg in range(0,len(SecondTH1)):
       legend[leg].AddEntry(SecondTH1[leg],SecondLegendName[leg]+" ("+str(round(secondIntegralList[leg],0))+")", "l")
       
   
   if(debug): print ("After for loop ")
   
   
           
   
   for canvasId in range(0,9):
        c.cd(canvasId+1)
        ROOT.gPad.SetTickx()
        ROOT.gPad.SetTicky()
        if(logy):
            ROOT.gPad.SetLogy(1)
        else:
            ROOT.gPad.SetLogy(0)
        
        if(logx):
            ROOT.gPad.SetLogx(1)
        else:
            ROOT.gPad.SetLogx(0)


        ROOT.gPad.SetLeftMargin(0.20)
        ROOT.gPad.SetRightMargin(0.05)
        ROOT.gPad.SetBottomMargin(0.1)
        drawStyle="hist"
        if ("TH2" in strType) and (drawPattern=="COLZ"):
            drawStyle = "COLZ"
            if(logz):
                ROOT.gPad.SetLogz(1)
            else:
                ROOT.gPad.SetLogz(0)
            
        else:
            drawStyle = "hist"
            
        FirstTH1[canvasId].Draw(drawStyle)
        SecondTH1[canvasId].Draw(drawStyle+" same")
        legend[canvasId].Draw()
        texList[canvasId].Draw("sames")
        ROOT.gPad.Update();
        if "TH2" in strType:
            pl = FirstTH1[canvasId].GetListOfFunctions().FindObject("palette")
            pl.SetX1NDC(0.905);
            pl.SetX2NDC(0.925);
            pl.SetY1NDC(0.10);
            pl.SetY2NDC(0.90);
        ROOT.gPad.Update();
        
        Tex2.Draw("sames")
        Tex3.Draw("sames")
        ROOT.gPad.RedrawAxis() 
        c.Update()
        
   if(debug): print ("After Draw loop ")
   
   L[0].Draw()
   L[1].Draw()
   L[2].Draw()
   
   SaveFile(c, CanvasName)
   return [c,L,legend]







#### this is to draw contours from 2D histograms,
def DrawContourHist(FirstTH1, LegendName, PlotColor, PlotStyle, nSigma, xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, logy=False, latexName='', latexName2 = '', latexName3='', leftLegend=False, doAtlas=False, doLumi=False, noRatio=False, do80=False, do59=False):
   debug = False
   if(debug): print ("just entering plot code")
   Tex = MakeLatex(0.54,0.74,latexName)
   Tex2 = MakeLatex(0.54,0.68,latexName2)
   Tex3 = MakeLatex(0.54,0.62,latexName3)
   if(debug): print ("defining Tex ")
   WaterMark = TexWaterMark('Preliminary')
   c = TCanvas("c","c",600, 600)
   ROOT.gStyle.SetOptStat(0)
   c.cd()
   c.SetGrid()
   if(logy):
     c.SetLogy()
   if(debug): print ("Set Logy ")
   if(leftLegend):
       legend1 = LeftLegendMaker()
   else:
       #legend1 = LegendMaker()
       legend1 = LegendMakerTwoColumn()
   if(debug): print ("Set Legend ")
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L, doAtlas, doLumi, noRatio, do80, do59)
   stList = []
   FirstTH1[0].GetYaxis().SetRangeUser(yrange1down,yrange1up)
   if(debug): print ("Set Ranges ")
   ### here is the uncertainty range of the contour
   dolcontour=array.array('d',[nSigma])
   ### set the drawing parameters of individual histograms
   for i in range(0, len(FirstTH1)):
     FirstTH1[i].SetContour( 1 , dolcontour);
     FirstTH1[i].SetLineWidth( 3 );
     FirstTH1[i].SetLineStyle( PlotStyle[i] );
     FirstTH1[i].SetLineColor( PlotColor[i] );
     FirstTH1[i] = AxisLabelEtc(FirstTH1[i], yAxisName, xAxisName)
     #FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
     #FirstTH1[i] = getOverflow(FirstTH1[i])
     legend1.AddEntry(FirstTH1[i], LegendName[i], "lp")
   
   if(debug): print ("After for loop ")
   FirstTH1[0].GetXaxis().SetRangeUser(xrange1down,xrange1up)
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
   FirstTH1[0].Draw("cont3 same") # ce, hist
   if(debug): print ("After first Draw ")
   if(len(FirstTH1)>1):
    for i in range(1, len(FirstTH1)):
        FirstTH1[i].Draw("cont3 same")
        
   if(debug): print ("After Draw loop ")
   
   Tex.Draw("sames")
   Tex2.Draw("sames")
   Tex3.Draw("sames")
   L[0].Draw()
   L[1].Draw()
   L[2].Draw()
   legend1.Draw()
   SaveFile(c, CanvasName)



def CalculateEffTGraph(hPass, hAll):
    g_efficiency = TGraphAsymmErrors()
    g_efficiency.Divide(hPass,hAll,"cl=0.683 b(1,1) mode")
    return g_efficiency

def CalculateEfficiency(hPass, hAll):
    
    hEff = hPass.Clone('hEff')
    hEff.Reset()
    tEff = TEfficiency()
    
    for ibin in range(1,hPass.GetXaxis().GetNbins()+1):
        if(hAll.GetBinContent(ibin)!=0):
            CentralValue = hPass.GetBinContent(ibin)/hAll.GetBinContent(ibin)
            err = tEff.ClopperPearson(int(abs(hAll.GetBinContent(ibin))),int(abs(hPass.GetBinContent(ibin))),0.683,True)
            ehi = tEff.ClopperPearson(int(abs(hAll.GetBinContent(ibin))),int(abs(hPass.GetBinContent(ibin))),0.683,True)
            elo = tEff.ClopperPearson(int(abs(hAll.GetBinContent(ibin))),int(abs(hPass.GetBinContent(ibin))),0.683,False)
            hEff.SetBinContent(ibin,CentralValue)
            err = 0.5*TMath.Sqrt(pow(ehi-CentralValue,2)+pow(elo-CentralValue,2))
            hEff.SetBinError(ibin,err)
        else:
            hEff.SetBinContent(ibin, 0)
            hEff.SetBinError(ibin, 0)
        
        
    return hEff



def EfficiencyHists(FirstTH1, SecondTH1, LegendName, PlotColor, xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, yline1low, yline1up, drawline=False, logy=False, LatexName='', leftLegend=False, doAtlas=False, doIntLum=False, noRatio=False):
    
   ### sanity check
   if(len(FirstTH1)!= len(SecondTH1)):
       print ("Two histogram lists are not of same size. Exiting!") 
       exit()
   
   Tex = MakeLatex(0.90,0.42,LatexName)
   c = TCanvas("c","c",500, 500)
   ROOT.gStyle.SetOptStat(0)
   c.cd()
   #c.SetGrid()
   if(logy):
     c.SetLogy()
   line = MakeLine(xrange1down,yline1low,xrange1up,yline1up)
   if(leftLegend):
       legend1 = LeftLegendMaker()
   else:
       legend1 = LegendMaker()
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L, doAtlas, doIntLum, noRatio)
   stList = []
   
   ratioHistList = []
   ratioOfHist = FirstTH1[0].Clone('ratioOfHist')
   ratioOfHist.Reset()
   ratioOfHist2 = FirstTH1[0].Clone('ratioOfHist2')
   ratioOfHist2.Reset()
   ratioOfHist3 = FirstTH1[0].Clone('ratioOfHist3')
   ratioOfHist3.Reset()
   ratioOfHist4 = FirstTH1[0].Clone('ratioOfHist4')
   ratioOfHist4.Reset()
   ratioOfHist5 = FirstTH1[0].Clone('ratioOfHist5')
   ratioOfHist5.Reset()
   ratioOfHist6 = FirstTH1[0].Clone('ratioOfHist6')
   ratioOfHist6.Reset()
   
   
   
   #if(TEfficiency.CheckConsistency(FirstTH1[0],SecondTH1[0])):
   ratioOfHist = CalculateEffTGraph(FirstTH1[0], SecondTH1[0])
   ratioHistList.append(ratioOfHist)
   ratioOfHist2 = CalculateEffTGraph(FirstTH1[1], SecondTH1[1])
   ratioHistList.append(ratioOfHist2)
   ratioOfHist3 = CalculateEffTGraph(FirstTH1[2], SecondTH1[2])
   ratioHistList.append(ratioOfHist3)
   ratioOfHist4 = CalculateEffTGraph(FirstTH1[3], SecondTH1[3])
   ratioHistList.append(ratioOfHist4)
   ratioOfHist5 = CalculateEffTGraph(FirstTH1[4], SecondTH1[4])
   ratioHistList.append(ratioOfHist5)
   ratioOfHist6 = CalculateEffTGraph(FirstTH1[5], SecondTH1[5])
   ratioHistList.append(ratioOfHist6)
   
   
   #for i in range(0, len(FirstTH1)):
     #ratioOfHist.Divide(FirstTH1[i], SecondTH1[i])
     #ratioHistList.append(ratioOfHist)
     
   for i in range(0, len(ratioHistList)):
     ratioHistList[i].GetYaxis().SetTitle(yAxisName)
     ratioHistList[i].GetXaxis().SetTitle(xAxisName)
     ratioHistList[i].GetYaxis().SetTitleOffset(0.96)
     ratioHistList[i].GetYaxis().SetLabelSize(0.052)
     ratioHistList[i].GetYaxis().SetTitleSize(0.048)
     ratioHistList[i].GetXaxis().SetTitleSize(0.048)
     ratioHistList[i].GetXaxis().SetLabelSize(0.035)
     ratioHistList[i].GetXaxis().SetTitleOffset(0.88)
     ratioHistList[i].GetXaxis().SetNdivisions(5)
     ratioHistList[i].SetTitle("")
     ratioHistList[i] = SetHistColorEtc(ratioHistList[i], PlotColor[i])
     ### comment in if you want only shape distribution
     #ratioHistList[i].Scale(1.0/w)
     legend1.AddEntry(ratioHistList[i],LegendName[i], "lp")
   
   ratioHistList[5].GetYaxis().SetRangeUser(yrange1down,yrange1up)
   ratioHistList[5].GetYaxis().SetTitle("Efficiency")
   ratioHistList[5].GetXaxis().SetRangeUser(xrange1down,xrange1up)
   ratioHistList[5].Draw("Aep")
   for i in range(0, len(FirstTH1)-1):
     ratioHistList[i].Draw("ep sames")  
   Tex.Draw("sames")
   if(drawline):
     line.Draw()
   L[0].Draw()
   L[1].Draw()
   L[2].Draw()
   legend1.Draw()
   #st = FirstTH1[0].GetListOfFunctions().FindObject("stats")
   #st.SetX1NDC(0.70); #new x start position
   #st.SetX2NDC(0.90); #new x end position
   #st.SetLineColor(kBlue);
   #st.SetLineColor(PlotColor[0])
   #stList.append(st)
   #for m in range(0, len(stList)):
     #x1NDC = 0.25+0.25*m
     #x2NDC = 0.45+0.25*m
     #stList[m].SetX1NDC(x1NDC); #new x start position
     #stList[m].SetX2NDC(x2NDC);
   
   
   SaveFile(c, CanvasName)
   #return [c,L,legend1]



def DrawTGraph(n, x, y, ex, ey, xAxisTitle, yAxisTitle, xrange1down, xrange1up, yrange1down, yrange1up, yline1low, yline1up, canvasName, latexName='', latexName2='', drawLine=False, logy=False, doAtlas=False):
    c = TCanvas("c","c",500, 500)
    ROOT.gStyle.SetOptStat(0)
    c.cd()
    c.SetGrid()
    tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
    L = [tex1, tex2, tex3]
    TexMaker(L, doAtlas)
    if(logy):
        c.SetLogy()
    line = MakeLine(xrange1down,yline1low,xrange1up,yline1up)
    
    Tex = MakeLatex(0.90,0.72,latexName)
    Tex2 = MakeLatex(0.90,0.65,latexName2)
    
    gr = TGraphErrors(n,x,y,ex,ey);
    gr.SetMarkerColor(4);
    gr.SetMarkerStyle(kOpenCircle);
    gr.SetMarkerSize(0.3);
    gr.SetLineColor(4);
    gr.SetTitle('');
    #gr.GetXaxis().SetTitle('#font[52]{'+xAxisTitle+'}')
    #gr.GetYaxis().SetTitle('#font[52]{'+yAxisTitle+'}')
    gr.GetXaxis().SetTitle(xAxisTitle)
    gr.GetYaxis().SetTitle(yAxisTitle)
    gr.GetXaxis().SetRangeUser(xrange1down, xrange1up)
    gr.GetYaxis().SetRangeUser(yrange1down, yrange1up)
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    gr.Draw("ALP")
    line.Draw("sames")
    Tex.Draw("sames")
    Tex2.Draw("sames")
    L[0].Draw("sames")
    L[2].Draw("sames")
    L[1].Draw("sames")
    
    SaveFile(c, canvasName)
    
    
def DrawTMultiGraph(GrList, LegList, ColorList, title, xAxisTitle, yAxisTitle, xrange1down, xrange1up, yrange1down, yrange1up, yline1low, yline1up, canvasName, latexName='', drawLine=False, logy=False, leftLegend=False, TeVTag=False, latexName2='', doAtlas=False):
    #### these are for MoEDAL Photon fusion theory paper
    #Tex = MakeLatex(0.88,0.65,latexName)
    #Tex2 = MakeLatex(0.88, 0.60, latexName2)
    ### these are for other cases
    Tex = MakeLatex(0.58,0.35,latexName)
    Tex2 = MakeLatex(0.64, 0.30, latexName2)
    ### this is for xsec plots with different kappa
    #Tex = MakeLatex(0.58,0.55,latexName)
    #Tex2 = MakeLatex(0.58, 0.50, latexName2)
    c = TCanvas("c","c",900, 900)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gPad.SetRightMargin(0.1)
    ROOT.gPad.SetLeftMargin(0.18)
    ROOT.gPad.SetBottomMargin(0.1)
    c.cd()
    c.SetGrid()
    if(logy):
       c.SetLogy()
    line = MakeLine(xrange1down,yline1low,xrange1up,yline1up)
    
    if(leftLegend):
       legend1 = LeftLegendMaker()
    else:
       #legend1 = LegendMakerTwoColumn()
       legend1 = LegendMaker()
       
    tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
    L = [tex1, tex2, tex3]
    TexMaker(L, doAtlas)
    
    mg = TMultiGraph('mg','')
    for i in range(0, len(GrList)):
       GrList[i].GetYaxis().SetTitle('#font[52]{'+yAxisTitle+'}')
       GrList[i].GetXaxis().SetTitle('#font[52]{'+xAxisTitle+'}')
       #GrList[i].GetXaxis().SetNdivisions(5)
       GrList[i].SetTitle("")
       GrList[i] = SetHistColorEtc(GrList[i], ColorList[i])
       #if((len(GrList) > 3) and (i%2==0)):
           #GrList[i].SetLineStyle(2)
           #if(i==0):GrList[i].SetMarkerStyle(kOpenSquare)
           #if(i==2):GrList[i].SetMarkerStyle(kOpenTriangleDown)
           #if(i==4):GrList[i].SetMarkerStyle(kOpenTriangleUp)
       #else:
           #GrList[i].SetLineStyle(1)
       GrList[i].SetMarkerSize(1.5)
       mg.Add(GrList[i])
       legend1.AddEntry(GrList[i],LegList[i],'ep')
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    mg.Draw('ALP')
    #mg.GetXaxis().SetTitle('#font[52]{'+xAxisTitle+'}')
    #mg.GetYaxis().SetTitle('#font[52]{'+yAxisTitle+'}')
    mg.GetXaxis().SetTitle(xAxisTitle)
    mg.GetYaxis().SetTitle(yAxisTitle)
    mg.GetXaxis().SetRangeUser(xrange1down, xrange1up)
    mg.GetYaxis().SetRangeUser(yrange1down, yrange1up)
    mg.GetYaxis().SetTitleOffset(2.05)
    mg.GetYaxis().SetLabelSize(0.031)
    mg.GetYaxis().SetTitleSize(0.036)
    mg.GetXaxis().SetTitleSize(0.04)
    mg.GetXaxis().SetLabelSize(0.037)
    #mg.GetYaxis().SetNdivisions(601)
    
    Tex.Draw('sames')
    Tex2.Draw('sames')
    #legend1.Draw('sames')
    if(drawLine):
       line.Draw()
    L[2].Draw()
    L[1].Draw()
    if(TeVTag):
       #TexTeV = TLatex(0.82,0.915,"#sqrt{s}=13 TeV")#36.075 fb^{-1} (13 TeV)
       TexTeV = TLatex(0.36,0.80,'#sqrt{s}=13 TeV')
       TexTeV.SetNDC()
       TexTeV.SetTextAlign(31)
       TexTeV.SetTextFont(42)
       TexTeV.SetTextSize(0.037)
       TexTeV.SetLineWidth(2) 
       TexTeV.Draw()
    else:
       L[0].Draw()

    SaveFile(c, canvasName)
    return deepcopy(c)
    
   
def DrawHistsBand(FirstTH1, LegendName,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, logy=False):
   c = TCanvas("c","c",500, 500)
   ROOT.gStyle.SetOptStat(0)
   c.cd()
   c.SetGrid()
   if(logy):
     c.SetLogy()
   legend1 = LegendMaker()
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L)
   for i in range(0, len(FirstTH1)):
     if i==0:
      FirstTH1[i].GetYaxis().SetTitle(yAxisName)
      FirstTH1[i].GetXaxis().SetTitle(xAxisName)
      FirstTH1[i].GetYaxis().SetTitleOffset(1.0)
      FirstTH1[i].GetYaxis().SetTitleSize(0.035)
      FirstTH1[i].GetXaxis().SetTitleSize(0.03)
      FirstTH1[i].GetXaxis().SetLabelSize(0.03)
      FirstTH1[i].SetLineColor(kBlue)
      legend1.AddEntry(FirstTH1[i],LegendName[i], 'F')
      FirstTH1[i].SetLineColor(kGray)
      FirstTH1[i].SetFillStyle(1001)
      FirstTH1[i].SetFillColor(kGray)
      FirstTH1[i].SetMarkerColor(kGray)
     else: 
       legend1.AddEntry(FirstTH1[i],LegendName[i], 'lep')
       FirstTH1[i].SetLineColor(kBlue)
       FirstTH1[i].SetMarkerColor(kBlue)
   
   FirstTH1[0].GetYaxis().SetRangeUser(yrange1down,yrange1up)
   FirstTH1[0].GetXaxis().SetRangeUser(xrange1down,xrange1up)
   FirstTH1[0].Draw("e2")
   for i in range(1, len(FirstTH1)):
     FirstTH1[i].Draw("ep sames")  
   
   L[0].Draw()
   L[1].Draw()
   L[2].Draw()
   legend1.Draw()
   
   SaveFile(c, CanvasName)
   c.Clear()
   
   
   
   
def TwoCanvasMaker(FirstTH1, SecondTH1, LegendName, PlotColor, xrange1down, xrange1up, yrange1down, yrange1up, xrange2down, xrange2up, yrange2down, yrange2up, CanvasName, yline1low, yline1up, 
yline2low, yline2up, drawline=False, logY=False):
   c = TCanvas("c","c",1200, 400);
   if(logY):
       c.SetLogy();
   c.Divide(2,1);
   ROOT.gStyle.SetOptStat(1111);
   c.cd(1);
   
   line = MakeLine(xrange1down,yline1low,xrange1up,yline1up);
   line2 = MakeLine(xrange2down,yline2low,xrange2up,yline2up);
   
   legend1 = LegendMaker();
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L)
   
   for i in range(0, len(FirstTH1)):
     FirstTH1[i].GetYaxis().SetRangeUser(yrange1down,yrange1up);
     FirstTH1[i].GetXaxis().SetRangeUser(xrange1down,xrange1up);
     FirstTH1[i].SetLineColor(PlotColor[i]);
     legend1.AddEntry(FirstTH1[i],LegendName[i], "lp");
   
   FirstTH1[0].Draw();
   for i in range(1, len(FirstTH1)):
     FirstTH1[i].Draw("sames");  
   
   L[0].Draw()
   L[1].Draw()
   L[2].Draw()
   legend1.Draw()
   
   
   c.cd(2);
   
   legend2 = LegendMaker();
   
   for i in range(0, len(SecondTH1)):
     SecondTH1[i].GetYaxis().SetRangeUser(yrange2down,yrange2up);
     SecondTH1[i].GetXaxis().SetRangeUser(xrange2down,xrange2up);
     SecondTH1[i].SetLineColor(PlotColor[i]);
     legend2.AddEntry(SecondTH1[i], LegendName[i], "lp");
   
   SecondTH1[0].Draw();
   for i in range(1, len(SecondTH1)):
       SecondTH1[i].Draw("sames");  
       
   if(drawline):
       line2.Draw();
   L[0].Draw()
   L[1].Draw()
   L[2].Draw()
   legend2.Draw()

   SaveFile(c, CanvasName)
   c.Clear()
  

   
   
def Draw2DHists(FirstTH1, LegendName, xAxisName, yAxisName, zAxisName, xrange1down, xrange1up, yrange1down, yrange1up, zrange1down, zrange1up, CanvasName, logz=False, LatexName='', LatexName2='', TeVTag=False, doAtlas=False):
   Tex = MakeLatex(0.85,0.72,LatexName)
   Tex2 = MakeLatex(0.85,0.65,LatexName2)
   c = TCanvas("c","c",500, 500)
   ROOT.gStyle.SetOptStat(0)
   #ROOT.gStyle.SetStatY(0.9);                
   #ROOT.gStyle.SetStatX(0.3);                
   #ROOT.gStyle.SetStatW(0.2);                
   #ROOT.gStyle.SetStatH(0.2);                
   c.cd()
   c.SetGrid()
   if(logz):
     c.SetLogz()
   #ROOT.gPad.SetRightMargin(4.0)
   #ROOT.gPad.Update()
   legend1 = LegendMaker()
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L, doAtlas)
   for i in range(0, len(FirstTH1)):
     FirstTH1[i].SetTitle('')
     
     FirstTH1[i].GetXaxis().SetTitle(xAxisName)
     FirstTH1[i].GetXaxis().SetTitleSize(0.04)
     FirstTH1[i].GetXaxis().SetLabelSize(0.03)
     FirstTH1[i].GetXaxis().SetTitleOffset(1.1)
     
     FirstTH1[i].GetYaxis().SetTitle(yAxisName)
     FirstTH1[i].GetYaxis().SetTitleOffset(1.3)
     FirstTH1[i].GetYaxis().SetTitleSize(0.04)
     FirstTH1[i].GetYaxis().SetLabelSize(0.035)
     
     FirstTH1[i].GetZaxis().SetTitle(zAxisName)
     FirstTH1[i].GetZaxis().SetLabelSize(0.032)
     #FirstTH1[i].GetYaxis().SetNdivisions(8)
     FirstTH1[i].GetZaxis().SetNdivisions(5)
     FirstTH1[i].GetZaxis().SetLabelOffset(0.006)
     FirstTH1[i].GetZaxis().SetTitleOffset(1.2)
     FirstTH1[i].GetZaxis().SetTitleSize(0.027)
     #legend1.AddEntry(FirstTH1[i],LegendName[i], "l")
       
   FirstTH1[0].GetZaxis().SetRangeUser(zrange1down,zrange1up)
   FirstTH1[0].GetYaxis().SetRangeUser(yrange1down,yrange1up)
   FirstTH1[0].GetXaxis().SetRangeUser(xrange1down,xrange1up)
   
   
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
   ROOT.gPad.SetRightMargin(0.1)
   ROOT.gPad.SetLeftMargin(0.15)
   ROOT.gPad.SetBottomMargin(0.1)
   FirstTH1[0].Draw('COLZ')
   if len(FirstTH1) > 1:
       for i in range(1, len(FirstTH1)):
           FirstTH1[i].Draw("sames")  
   ROOT.gPad.Update();
   pl = FirstTH1[0].GetListOfFunctions().FindObject("palette")
   #pl.SetX1NDC(0.955);
   #pl.SetX2NDC(0.96);
   #pl.SetY1NDC(0.16);
   #pl.SetY2NDC(0.955);
   pl.SetX1NDC(0.905);
   pl.SetX2NDC(0.925);
   pl.SetY1NDC(0.1);
   pl.SetY2NDC(0.9);
   ROOT.gPad.Update();
   
   #ps = FirstTH1[0].GetListOfFunctions().FindObject("stats");
   #ps.SetFillStyle(0)
   #ROOT.gPad.Update()
   #st = FirstTH1[0].FindObject("stats")# if not stats: continue stats.__class__ = ROOT.TPaveStats
   #st.SetX1NDC(0.70); #new x start position
   #st.SetX2NDC(0.90); #new x end position
   #st.SetY1NDC(0.70); #new x start position
   #st.SetY2NDC(0.90); #new x end position
   #st.SetLineColor(kBlue);
   #st.SetFillColor(kWhite);
   ROOT.gPad.Modified();
   ROOT.gPad.Update();
   Tex.Draw("sames")
   Tex2.Draw("sames")
   if(TeVTag):
       TexTeV = TLatex(0.892,0.914,"13 TeV")#36.075 fb^{-1} (13 TeV)
       TexTeV.SetNDC()
       TexTeV.SetTextAlign(31)
       TexTeV.SetTextFont(42)
       TexTeV.SetTextSize(0.037)
       TexTeV.SetLineWidth(2) 
       TexTeV.Draw()
   else:
       L[0].Draw()
       #print ('nothing')
   
   L[1].Draw()
   L[2].Draw()
   #legend1.Draw()
   
   SaveFile(c, CanvasName)
   return [c, L, legend1]



#### #### converting TProfile in to TH1D
##### you may need to change the definition of Hist according to your need
def convertTProfiletoTH1(EffValue):
    nXBins = EffValue.GetNbinsX()
    
    Hist = TH1D('Hist','Hist',nXBins, -30, 30)
    for xBin in range(0, nXBins):
        Hist.SetBinContent(xBin, EffValue.GetBinContent(xBin))
        Hist.SetBinError(xBin, EffValue.GetBinError(xBin))
            
    
    return Hist


#### converting TProfile2D in to TH2D
##### you may need to change the definition of Hist2D according to your need
def convertTProfile2DtoTH2(EffValue):
    nXBins = EffValue.GetNbinsX()
    nYBins = EffValue.GetNbinsY()
    
    Hist2D = TH2D('Hist2D','Hist2D',nXBins, -6, 6, nYBins, 0, 56)
    for xBin in range(0, nXBins):
        for yBin in range(0, nYBins):
            Hist2D.SetBinContent(xBin, yBin, EffValue.GetBinContent(xBin, yBin))
            Hist2D.SetBinError(xBin, yBin, EffValue.GetBinError(xBin, yBin))
            
    
    return Hist2D



def drawText2D(hist, canvasName, title):
    c1 = TCanvas('c1','c1',900, 600)
    c1.cd()
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gStyle.SetPalette(55)
    ROOT.gStyle.SetPaintTextFormat("2.4f")
    hist.SetMarkerSize(1.7)
    hist.Draw('COLZ text')
    hist.GetZaxis().SetTitle('Probability (KS)')
    hist.SetTitle(title)
    
    ROOT.gPad.Update();
    pl = hist.GetListOfFunctions().FindObject("palette")
    pl.SetX1NDC(0.91);
    pl.SetX2NDC(0.93);
    pl.SetY1NDC(0.10);
    pl.SetY2NDC(0.90);
    ROOT.gPad.Modified();
    ROOT.gPad.Update();
    
    ROOT.gStyle.SetOptStat(0)
    
    y = -0.5 #hist_DeltaPhill_data17.GetYaxis().GetBinWidth(1);
    t = TText()
    t.SetTextAngle(20)
    t.SetTextSize(0.03)
    t.SetTextAlign(60)
    hist.GetXaxis().SetLabelOffset(100)
    hist.GetYaxis().SetLabelOffset(100)
    
    x = hist.GetXaxis().GetBinCenter(1)
    t.DrawText(x,y,'100 < MET < 200')
    x = hist.GetXaxis().GetBinCenter(2)
    t.DrawText(x,y,'200 < MET < 400')
    x = hist.GetXaxis().GetBinCenter(3)
    t.DrawText(x,y,' MET > 400')
    
    x2 = -0.35 #hist_DeltaPhill_data17.GetYaxis().GetBinWidth(1);
    t2 = TText()
    t2.SetTextAngle(30)
    t2.SetTextSize(0.025)
    t2.SetTextAlign(60)
    
    y2 = hist.GetYaxis().GetBinCenter(1)
    t2.DrawText(x2,y2,'100 < HT < 200')
    y2 = hist.GetYaxis().GetBinCenter(2)
    t2.DrawText(x2,y2,'200 < HT < 400')
    y2 = hist.GetYaxis().GetBinCenter(3)
    t2.DrawText(x2,y2,'400 < HT < 700')
    y2 = hist.GetYaxis().GetBinCenter(4)
    t2.DrawText(x2,y2,' HT > 700')
    
    c1.SaveAs(canvasName)


   
   
def DrawHistsAllStack(FirstTH1, LegendName, PlotColor,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, yline1low, yline1up, drawline=False, logy=False, 
LatexName='', TextLabel=False, labelName='', doAtlas=False, doIntLum=False, noRatio=False):
   Tex = MakeLatex(0.90,0.62,LatexName)
   c = TCanvas("c","c",900, 900)
   ROOT.gStyle.SetOptStat(0)
   c.cd()
   c.SetGrid()
   if(logy):
     c.SetLogy()
   legend1 = LegendMakerTwoColumn()
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L, doAtlas, doIntLum, noRatio)
   hs = THStack("hs","");
   for i in range(0, len(FirstTH1)):
     lenHistList = len(FirstTH1)
     FirstTH1[i] = getOverflow(FirstTH1[i])
     FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
     hs.Add(FirstTH1[i])
     legend1.AddEntry(FirstTH1[lenHistList-1-i],LegendName[lenHistList-1-i], "f") ### this will add reverse legends
  
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
   hs.Draw("hist")
   hs.GetXaxis().SetRangeUser(xrange1down,xrange1up)
   hs.SetMaximum(yrange1up)
   hs.SetMinimum(yrange1down)
   hs.GetYaxis().SetTitle(yAxisName)
   hs.GetXaxis().SetTitle(xAxisName)
   hs.GetYaxis().SetTitleOffset(0.96)
   hs.GetYaxis().SetTitleSize(0.045)
   hs.GetXaxis().SetTitleSize(0.05)
   hs.GetXaxis().SetLabelSize(0.03)
   hs.GetXaxis().SetLabelOffset(0)
   hs.GetXaxis().SetTitleOffset(0.8) # changed from 10
   hs.GetXaxis().SetNdivisions(5)
   if(TextLabel):
        y = -0.3 #ROOT.gPad.GetUymin() - 0.2*hs.GetYaxis().GetBinWidth(1)
        #print y
        t = TText()
        t.SetTextAngle(50)
        t.SetTextSize(0.02)
        t.SetTextAlign(33)
        
        for i in range(0, len(labelName)):
            x = hs.GetXaxis().GetBinCenter(i+1)
            #print x
            t.DrawText(x,y,labelName[i])
            
            
   Tex.Draw("sames")
   
  
   L[0].Draw("sames")
   L[1].Draw("sames")
   L[2].Draw("sames")
   legend1.Draw("sames")
      
   SaveFile(c, CanvasName)
   return [hs,L,legend1]
   c.Clear() 
      
      
   
   
   
### be careful as there may be an issue with THStack and TCanvas. Be careful about the object cList[0]. The first histogram is always plotted with black and kFullDotLarge.

def DrawHistsStack(FirstTH1, LegendName, PlotColor,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, yline1low, yline1up, drawline=False, logy=False, LatexName=None, doAtlas=False, doIntLum=False):
   Tex = MakeLatex(0.90,0.42,LatexName)
   AllButOne = []
   AllButOneLeg = []
   AllButOneColor = []
   for i in range(1, len(FirstTH1)):
       AllButOne.append(FirstTH1[i])
       AllButOneLeg.append(LegendName[i])
       AllButOneColor.append(PlotColor[i])
   
   cList = DrawHistsAllStack(AllButOne, AllButOneLeg, AllButOneColor,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, yline1low, yline1up, drawline, logy, '','','',doAtlas, doIntLum)
   c2 = TCanvas()
   c2.cd()
   if(logy):
       c2.SetLogy()
   FirstTH1[0].SetTitle("")
   FirstTH1[0].SetLineColor(kBlack)
   FirstTH1[0].SetMarkerStyle(kFullDotLarge)
   FirstTH1[0].SetMarkerColor(kBlack)
   FirstTH1[0].GetXaxis().SetTitle(xAxisName)
   FirstTH1[0].GetYaxis().SetTitle(yAxisName)
   FirstTH1[0].GetYaxis().SetTitleOffset(1.0)
   cList[0].Draw("hist")
   FirstTH1[0].Draw("ep sames")
   cList[1][0].Draw("sames")
   cList[1][1].Draw("sames")
   cList[1][2].Draw("sames")
   cList[2].AddEntry(FirstTH1[0],LegendName[0], "lp")
   cList[2].Draw("sames")
   Tex.Draw("sames")
   
   
   WriteFile = TFile(CanvasName+'.root', 'RECREATE')
   WriteFile.cd()
   c2.Write()
   
   SaveFile(c2, CanvasName)
   
   return [cList[0],cList[1],cList[2]]
   c2.Clear()
   
   
   
   
#### to draw with multiple stack and multiple overlay plot
   
def DrawHistsStackSignal(FirstTH1, SignalTH1, LegendName, LegendNameSignal, PlotColor, PlotColorSignal, xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, yline1low, yline1up, drawline=False, logy=False, LatexName=None, doAtlas=False, doIntLum=False, doRatio=False):
   Tex = MakeLatex(0.70,0.52,LatexName)
   cList = DrawHistsAllStack(FirstTH1, LegendName, PlotColor, xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, yline1low, yline1up, drawline, logy, '','','',doAtlas, doIntLum, doRatio)
   c2 = TCanvas()
   c2.cd()
   if(logy):
       c2.SetLogy()
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
   cList[0].Draw("hist")
   
   for i in range(0, len(SignalTH1)):    
        SignalTH1[i].SetTitle("")
        SignalTH1[i] = SetHistColorEtc(SignalTH1[i], PlotColorSignal[i])
        SignalTH1[i].GetXaxis().SetTitle(xAxisName)
        SignalTH1[i].GetYaxis().SetTitle(yAxisName)
        SignalTH1[i].GetYaxis().SetTitleOffset(1.0)
        SignalTH1[i].SetLineWidth(3)
        cList[2].AddEntry(SignalTH1[i],LegendNameSignal[i], "lp")
        SignalTH1[i].SetFillColor(0)
        SignalTH1[i].Draw("hist sames")
   
   
   cList[1][0].Draw("sames")
   cList[1][1].Draw("sames")
   cList[1][2].Draw("sames")
   cList[2].Draw("sames")
   Tex.Draw("sames")
   
   
   WriteFile = TFile(CanvasName+'.root', 'RECREATE')
   WriteFile.cd()
   c2.Write()
   SaveFile(c2, CanvasName)
   


#### draws only one signal
def DrawHistsStackOneSignalSignificance(FirstTH1, SignalTH1, LegendName, LegendNameSignal, PlotColor, PlotColorSignal, xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, yline1low, yline1up, drawline=False, logy=False, LatexName=None, doAtlas=False, doIntLum=False, doRatio=False):
   
   c2 = TCanvas("c2","c2",900, 900)
   
   pad1 = TPad('pad1','pad1',0,0.25,1,1)
   pad1.SetBottomMargin(0.0)
   pad1.SetFillStyle(0)
   pad1.SetGrid()
   pad1.Draw()
   pad1.cd()
   ROOT.gStyle.SetOptStat(0)
   if(logy):
     pad1.SetLogy()
   
   Tex = MakeLatex(0.85,0.52,LatexName)
   legend1 = LegendMaker()
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L, doAtlas, doIntLum, doRatio)
   hs = THStack("hs","");
   for i in range(0, len(FirstTH1)):
     FirstTH1[i] = getOverflow(FirstTH1[i])
     FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
     hs.Add(FirstTH1[i])
     legend1.AddEntry(FirstTH1[i],LegendName[i], "f")
  
 
   hs.Draw("hist")
   hs.SetMaximum(yrange1up)
   hs.SetMinimum(yrange1down)
   hs.GetXaxis().SetRangeUser(xrange1down,xrange1up)
   hs.GetYaxis().SetTitle(yAxisName)
   hs.GetXaxis().SetTitle(xAxisName)
   hs.GetYaxis().SetTitleOffset(0.7)
   hs.GetYaxis().SetTitleSize(0.045)
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
   hs.Draw("hist")
   hs.GetHistogram().GetXaxis().SetLabelOffset(999);
   Tex.Draw("sames")
   
   for i in range(0, len(SignalTH1)):    
        SignalTH1[i].SetTitle("")
        SignalTH1[i] = getOverflow(SignalTH1[i])
        SignalTH1[i] = SetHistColorEtc(SignalTH1[i], PlotColorSignal[i])
        SignalTH1[i].SetLineWidth(3)
        SignalTH1[i].GetXaxis().SetTitle(xAxisName)
        SignalTH1[i].GetYaxis().SetTitle(yAxisName)
        SignalTH1[i].GetYaxis().SetTitleOffset(1.0)
        legend1.AddEntry(SignalTH1[i],LegendNameSignal[i], "l")
        SignalTH1[i].SetFillColor(0)
        SignalTH1[i].Draw("hist sames")
   
   
   L[0].Draw("sames")
   L[1].Draw("sames")
   L[2].Draw("sames")
   legend1.Draw("sames")
   Tex.Draw("sames")
   
   c2.cd()
   pad2 = TPad("pad2", "pad2",0.0,0.0,1.0,0.249)
   pad2.SetTopMargin(0.0)
   pad2.SetBottomMargin(0.3)
   pad2.SetFillStyle(0)
   #pad2.SetLogy()
   #pad2.SetGrid()
   pad2.Draw()
   pad2.cd()
   ROOT.gStyle.SetOptStat(0)
   
   h1Sig = SignalTH1[0].Clone('h1Sig')
   h1Sig.Reset()
   h1Main = SignalTH1[0].Clone('h1Main')
   
   histAll = FirstTH1[0].Clone('histAll')
   histAll.Reset()
   
   for i in range(0, len(FirstTH1)):
       histAll.Add(FirstTH1[i])
   
   for i in range(0, h1Sig.GetNbinsX()):
       signal1       = h1Main.Integral(i, h1Sig.GetNbinsX()+1)
       errB = Double(0)
       background    = histAll.IntegralAndError(i, h1Sig.GetNbinsX()+1, errB)
       #### adding a fake 30% systematic for now
       flatSyst = 0.3
       
       if(background!=0):
            totErr = math.sqrt((errB*errB)/(background*background)+flatSyst*flatSyst)
            significance1 = RooStats.NumberCountingUtils.BinomialExpZ(signal1, background, totErr)
       else: 
            significance1 = 0
       
       h1Sig.SetBinContent(i, significance1)
       
   
   legend2 = RatioLegendMaker()
   legend2.AddEntry(h1Sig, LegendNameSignal[0],"lp")
   
   h1Sig.SetLineWidth(1)
   h1Sig.GetXaxis().SetRangeUser(xrange1down,xrange1up)
   h1Sig.GetXaxis().SetTitle(xAxisName)
   h1Sig.GetYaxis().SetTitle("Sig.")
   h1Sig.GetYaxis().CenterTitle()
   h1Sig.SetTitle('')
   h1Sig = OneRatioAxisSize(h1Sig, PlotColorSignal[0])
   h1Sig.GetYaxis().SetRangeUser(0.0, 2.4)
   h1Sig.SetFillColor(0)
   ROOT.gPad.SetTicky()
   h1Sig.Draw("lp")
   legend2.Draw("same")
   line = MakeLine(xrange1down,yline1low,xrange1up,yline1up)
   if(drawline):
     line.Draw("sames")
   #WriteFile = TFile(CanvasName+'.root', 'RECREATE')
   #WriteFile.cd()
   #c2.Write()
   SaveFile(c2, CanvasName)






#### draws three different signal samples, no significance
def DrawHistsStackBackgroundSignal(FirstTH1, SignalTH1, LegendName, LegendNameSignal, PlotColor, PlotColorSignal, xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, yline1low, yline1up, drawline=False, logy=False, LatexName='', doAtlas=False, doIntLum=False, doRatio=False):
   
   c2 = TCanvas("c2","c2",900, 900)
   ROOT.gStyle.SetOptStat(0)
   c2.SetGrid()
   c2.cd()
   ROOT.gStyle.SetOptStat(0)
   if(logy):
     c2.SetLogy()
   
   Tex = MakeLatex(0.85,0.52,LatexName)
   legend1 = LegendMakerTwoColumn()
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L, doAtlas, doIntLum, doRatio)
   hs = THStack("hs","");
   for i in range(0, len(FirstTH1)):
     FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
     hs.Add(FirstTH1[i])
     legend1.AddEntry(FirstTH1[i],LegendName[i], "f")
  
 
   hs.Draw("hist")
   hs.SetMaximum(yrange1up)
   hs.SetMinimum(yrange1down)
   hs.GetXaxis().SetRangeUser(xrange1down,xrange1up)
   hs.GetYaxis().SetTitle(yAxisName)
   hs.GetXaxis().SetTitle(xAxisName)
   hs.GetYaxis().SetTitleOffset(0.96)
   hs.GetYaxis().SetTitleSize(0.045)
   #hs.GetYaxis().SetTitleOffset(0.7)
   #hs.GetYaxis().SetTitleSize(0.045)
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
   hs.Draw("hist")
   #hs.GetHistogram().GetXaxis().SetLabelOffset(0.96);
   Tex.Draw("sames")
   
   for i in range(0, len(SignalTH1)):    
        SignalTH1[i].SetTitle("")
        SignalTH1[i] = SetHistColorEtc(SignalTH1[i], PlotColorSignal[i])
        SignalTH1[i].SetLineWidth(7)
        SignalTH1[i].GetXaxis().SetTitle(xAxisName)
        SignalTH1[i].GetYaxis().SetTitle(yAxisName)
        SignalTH1[i].GetYaxis().SetTitleOffset(1.0)
        legend1.AddEntry(SignalTH1[i],LegendNameSignal[i], "l")
        SignalTH1[i].SetFillColor(0)
        SignalTH1[i].Draw("hist sames")
   
   
   L[0].Draw("sames")
   L[1].Draw("sames")
   L[2].Draw("sames")
   legend1.Draw("sames")
   Tex.Draw("sames")
   
   c2.cd()
   line = MakeLine(xrange1down,yline1low,xrange1up,yline1up)
   if(drawline):
       line.Draw("sames")
   WriteFile = TFile(CanvasName+'.root', 'RECREATE')
   WriteFile.cd()
   c2.Write()
   SaveFile(c2, CanvasName)
   
   


#### draws many different signal samples
def DrawHistsStackSignalSignificance(FirstTH1, SignalTH1, LegendName, LegendNameSignal, PlotColor, PlotColorSignal, xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, yline1low, yline1up, drawline=False, logy=False, LatexName=None, doAtlas=False, doIntLum=False, doRatio=False):
   
   c2 = TCanvas("c2","c2",900, 900)
   
   pad1 = TPad('pad1','pad1',0,0.25,1,1)
   pad1.SetBottomMargin(0.0)
   pad1.SetFillStyle(0)
   pad1.SetGrid()
   pad1.Draw()
   pad1.cd()
   ROOT.gStyle.SetOptStat(0)
   if(logy):
     pad1.SetLogy()
     
   ### this is not drawing THStack on the canvas
   #cList = DrawHistsAllStack(FirstTH1, LegendName, PlotColor, xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, yline1low, yline1up, drawline, logy, '','','',doAtlas, doIntLum)
   
   Tex = MakeLatex(0.85,0.52,LatexName)
   legend1 = LegendMakerTwoColumn()
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L, doAtlas, doIntLum, doRatio)
   hs = THStack("hs","");
   for i in range(0, len(FirstTH1)):
     FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
     hs.Add(FirstTH1[i])
     legend1.AddEntry(FirstTH1[i],LegendName[i], "f")
  
 
   hs.Draw("hist")
   hs.SetMaximum(yrange1up)
   hs.SetMinimum(yrange1down)
   hs.GetXaxis().SetRangeUser(xrange1down,xrange1up)
   hs.GetYaxis().SetTitle(yAxisName)
   hs.GetXaxis().SetTitle(xAxisName)
   hs.GetYaxis().SetTitleOffset(0.96)
   hs.GetYaxis().SetTitleSize(0.045)
   #hs.GetYaxis().SetTitleOffset(0.7)
   #hs.GetYaxis().SetTitleSize(0.045)
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
   hs.Draw("hist")
   hs.GetHistogram().GetXaxis().SetLabelOffset(999);
   Tex.Draw("sames")
   
   for i in range(0, len(SignalTH1)):    
        SignalTH1[i].SetTitle("")
        SignalTH1[i] = SetHistColorEtc(SignalTH1[i], PlotColorSignal[i])
        SignalTH1[i].SetLineWidth(3)
        SignalTH1[i].GetXaxis().SetTitle(xAxisName)
        SignalTH1[i].GetYaxis().SetTitle(yAxisName)
        SignalTH1[i].GetYaxis().SetTitleOffset(1.0)
        legend1.AddEntry(SignalTH1[i],LegendNameSignal[i], "lp")
        SignalTH1[i].SetFillColor(0)
        SignalTH1[i].Draw("hist sames")
   
   
   L[0].Draw("sames")
   L[1].Draw("sames")
   L[2].Draw("sames")
   legend1.Draw("sames")
   Tex.Draw("sames")
   
   c2.cd()
   pad2 = TPad("pad2", "pad2",0.0,0.0,1.0,0.249)
   pad2.SetTopMargin(0.0)
   pad2.SetBottomMargin(0.3)
   pad2.SetFillStyle(0)
   #pad2.SetLogy()
   #pad2.SetGrid()
   pad2.Draw()
   pad2.cd()
   ROOT.gStyle.SetOptStat(0)
   
   histAll = FirstTH1[0].Clone('histAll')
   histAll.Reset()
   
   hSig = [] ; hMain = [] ;
   for i in range(0, len(SignalTH1)):
       nameSig  = 'h'+str(i)+'Sig'
       nameMain = 'h'+str(i)+'Main'
       hS = SignalTH1[i].Clone(nameSig)
       hS.Reset()
       hSig.append(hS)
       hM = SignalTH1[i].Clone(nameMain)
       hMain.append(hM)
       
   
   for i in range(0, len(FirstTH1)):
       histAll.Add(FirstTH1[i])
       
       
   for j in range(0, len(SignalTH1)):   
        for i in range(0, hSig[j].GetNbinsX()):
            signal      = hMain[j].Integral(i, hSig[j].GetNbinsX()+1)
            errB        = Double(0)
            background  = histAll.IntegralAndError(i, hSig[j].GetNbinsX()+1, errB)
            #### adding a fake 30% systematic for now
            flatSyst = 0.3
        
            if(background>0):
                totErr       = math.sqrt((errB*errB)/(background*background)+flatSyst*flatSyst)
                significance = RooStats.NumberCountingUtils.BinomialExpZ(signal, background, totErr)
                #print 'j: ',j, ' i: ', i, 'signal: ', signal, ' background: ', background, ' sgn: ', significance
            else: 
                significance = 0 #RooStats.NumberCountingUtils.BinomialExpZ(signal1, background, 0)
        
            ### for crude estimation f significance
            #if(background!=0):
                #significance1  = signal1/math.sqrt(background)
                #significance2  = signal2/math.sqrt(background)
                #significance3  = signal3/math.sqrt(background)
        
            hSig[j].SetBinContent(i, significance)
            
   legend2 = RatioLegendMaker()
   ROOT.gPad.SetTicky()
   if(len(hSig)<1): 
       print ('Something wrong with the signal histograms. Exiting')
       sys.exit(2)
   hSig[0] = OneRatioAxisSize(hSig[0], PlotColorSignal[0])
   hSig[0].SetLineWidth(2)
   legend2.AddEntry(hSig[0], LegendNameSignal[0],"lp")
   hSig[0].GetXaxis().SetRangeUser(xrange1down,xrange1up)
   hSig[0].GetXaxis().SetTitle(xAxisName)
   hSig[0].GetYaxis().SetTitle("Sig.")
   hSig[0].GetYaxis().CenterTitle()
   hSig[0].SetTitle('')
   hSig[0].GetYaxis().SetRangeUser(0.0, 4.5)
   hSig[0].SetFillColor(0)
   hSig[0].Draw("lp")
   for j in range(1, len(hSig)): 
       hSig[j] = OneRatioAxisSize(hSig[j], PlotColorSignal[j])
       hSig[j].SetFillColor(0)
       hSig[j].SetLineWidth(2)
       legend2.AddEntry(hSig[j], LegendNameSignal[j],"lp")
       hSig[j].Draw("lp same")
        
   legend2.Draw("same")
   line = MakeLine(xrange1down,yline1low,xrange1up,yline1up)
   if(drawline):
       line.Draw("sames")
   WriteFile = TFile(CanvasName+'.root', 'RECREATE')
   WriteFile.cd()
   c2.Write()
   SaveFile(c2, CanvasName)
   
   
   

  
### be careful as there may be an issue with THStack and TCanvas. This is long because cList cannot handle the TPad. The first histogram is always plotted with black and kFullDotLarge.
def DrawHistsStackRatio(FirstTH1, LegendName, PlotColor,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, yline1low, yline1up, drawline=False, logy=False, 
LatexName='', doAtlas=False, doIntLum=False, do80Inv=False):
   c3 = TCanvas("c3","c3",900, 900)
   c3.SetGrid()
   pad1 = TPad('pad1','pad1',0,0.25,1,1)
   pad1.SetBottomMargin(0.01)
   pad1.SetFillStyle(0)
   pad1.SetGrid()
   pad1.Draw()
   pad1.cd()
   ROOT.gStyle.SetOptStat(0)
   if(logy):
     pad1.SetLogy()
   
   Tex = MakeLatex(0.90,0.48,LatexName)
   legend1 = LegendMaker()
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L, doAtlas, doIntLum, False, do80Inv)
   FirstTH1[0] = getOverflow(FirstTH1[0])
   FirstTH1[0].SetLineColor(kBlack)
   FirstTH1[0].SetMarkerColor(kBlack)
   FirstTH1[0].SetMarkerStyle(kFullDotLarge)
   FirstTH1[0].SetFillColor(kBlack)
   legend1.AddEntry(FirstTH1[0], LegendName[0], "lp")
   hs = THStack("hs","");
   
   ### this is needed for ratio
   h2 = FirstTH1[0].Clone('h2')
   hErr = FirstTH1[0].Clone('h3')
   histAll = FirstTH1[0].Clone('histAll')
   histAll.Reset()
   #histAll.Sumw2()
   for i in range(1, len(FirstTH1)):
       histAll.Add(FirstTH1[i])
   histAll     = getOverflow(histAll)
       
   for i in range(1, len(FirstTH1)):
     FirstTH1[i] = getOverflow(FirstTH1[i])
     FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
     hs.Add(FirstTH1[i])
     if(histAll.Integral()!=0): percentage = FirstTH1[i].Integral()/float(histAll.Integral())
     else:                      percentage = 0.0
     if(percentage > 1): 
         print( "percentage ", percentage)
         print( 'numerator: ', FirstTH1[i].Integral())
         print( 'denominator: ', histAll.Integral())
     roundPercentage = round(percentage*100.0,2)
     legend1.AddEntry(FirstTH1[i],LegendName[i]+" ("+str(roundPercentage)+"%)", "f")
  
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
   hs.Draw("hist")
   hs.SetMaximum(yrange1up)
   hs.SetMinimum(yrange1down)
   hs.GetXaxis().SetRangeUser(xrange1down,xrange1up)
   hs.GetYaxis().SetTitle(yAxisName)
   hs.GetXaxis().SetTitle(xAxisName)
   hs.GetYaxis().SetTitleOffset(0.7)
   hs.GetYaxis().SetTitleSize(0.045)
   hs.GetXaxis().SetTitleSize(0.05)
   hs.GetXaxis().SetLabelSize(0.0)
   hs.GetXaxis().SetLabelOffset(999)
   hs.GetXaxis().SetTitleOffset(0.8)
   FirstTH1[0].Draw("ep sames")
   Tex.Draw("sames")
   
  
   L[0].Draw()
   L[1].Draw()
   L[2].Draw()
   legend1.Draw()
   
  
  
  
   c3.cd()
   pad2 = TPad("pad2", "pad2",0.0,0.0,1.0,0.249)
   pad2.SetTopMargin(0.05)
   pad2.SetBottomMargin(0.3)
   pad2.SetFillStyle(0)
   #pad2.SetLogy()
   #pad2.SetGrid()
   pad2.Draw()
   pad2.cd()
   ROOT.gStyle.SetOptStat(0)
   
   
   hErr.Reset()
   h2.Reset()
   h2.Divide(FirstTH1[0], histAll)
   sysErr = 0.2
   for i in range(0, h2.GetNbinsX()+1):
       hErr.SetBinContent(i, 1)
       statErr = h2.GetBinError(i)
       totErr  = math.sqrt(statErr*statErr+sysErr*sysErr)
       hErr.SetBinError(i, totErr)
   
   h2.GetXaxis().SetRangeUser(xrange1down,xrange1up)
   h2.GetXaxis().SetTitle(xAxisName)
   h2.GetYaxis().SetTitle("Data/SM")
   h2.GetYaxis().CenterTitle()
   h2.SetTitle('')
   h2 = OneRatioAxisSize(h2)
   h2.GetYaxis().SetRangeUser(0.5, 1.5)
   h2.GetYaxis().SetNdivisions(3)
   
   
   hErr.GetXaxis().SetRangeUser(xrange1down,xrange1up)
   hErr.GetXaxis().SetTitle(xAxisName)
   hErr.GetYaxis().SetTitle("Data/SM")
   hErr.GetYaxis().CenterTitle()
   hErr.SetTitle('')
   hErr = OneRatioAxisSize(hErr)
   hErr.GetYaxis().SetRangeUser(0.5, 1.5)
   hErr.GetYaxis().SetNdivisions(3)
   
   hErr.SetFillColor(kGray)
   hErr.SetMarkerColor(kGray)
   #hErr.Draw('e2')
   h2.Draw()
   line = MakeLine(xrange1down,yline1low,xrange1up,yline1up)
   if(drawline):
     line.Draw("sames")
     
   SaveFile(c3, CanvasName)
   c3.Clear()
   
   fileWrite = TFile(CanvasName+'.root','RECREATE')
   fileWrite.cd()
   h2.GetXaxis().SetLabelSize(0.06)
   h2.GetXaxis().SetTitleOffset(0.8)
   h2.GetXaxis().SetTitleSize(0.06)
   h2.GetYaxis().SetLabelSize(0.05)
   h2.GetYaxis().SetTitleOffset(0.55)
   h2.GetYaxis().SetTitleSize(0.06)
   histAll.Write('histAllMC')
   FirstTH1[0].Write('histData')
   h2.Write()
   fileWrite.Close()
   
   
   
   
   
   
   
   
   
def DrawHistsRatioJES(FirstTH1, FirstJES, LegendName, PlotColor, xrange1down, xrange1up, yrange1down, yrange1up, xaxisTitle, CanvasName, yline1low, yline1up, drawline=False, logy=False, 
leftLegend=False):
  
  TH1.SetDefaultSumw2()

  xmin=FirstTH1[1].GetX()[0]
  xmax=FirstTH1[1].GetX()[FirstTH1[1].GetN()-1]
  nbins=FirstTH1[1].GetN()
  
  hmcdy= TH1F("hmcdy","hmcdy", nbins,xmin,xmax)
  hdatady= TH1F("hdatady","hdatady", nbins,xmin,xmax)
  

  hmcdyup = TH1F("hmcdyup","hmcdyup", nbins,xmin,xmax)
  hmcdydown = TH1F("hmcdydown","hmcdydown", nbins,xmin,xmax)

  for i in range(0, FirstTH1[0].GetN()+1):
    hdatady.SetBinContent(i+1,FirstTH1[0].GetY()[i])
    hmcdy.SetBinContent(i+1,FirstTH1[1].GetY()[i])
    hmcdyup.SetBinContent(i+1,FirstJES[0].GetY()[i])
    hmcdydown.SetBinContent(i+1,FirstJES[1].GetY()[i])
             
           
    #setting the errors
    hdatady.SetBinError(i+1,FirstTH1[0].GetErrorY(i))
    hmcdy.SetBinError(i+1,FirstTH1[1].GetErrorY(i))
          

  hRatio = TH1F("hRatio","hRatio",nbins,xmin,xmax)
  hRatio.Divide(hdatady, hmcdy)
  x = array.array('f'); y=array.array('f'); ex=array.array('f'); ey=array.array('f');
  
  for q in range(0, nbins+1):
    x.append(hRatio.GetBinCenter(q))
    y.append(hRatio.GetBinContent(q))
    ex.append(hRatio.GetBinWidth(q)/2.0)
    ey.append(hRatio.GetBinError(q))
    
  
  h4 = TGraphErrors(nbins+1, x, y, ex, ey)
  #h4 = TH1F("h4","h4",nbins,xmin,xmax)
  #h4.Divide(hdatady, hmcdy)
  c = TCanvas("c","c",900, 900)
  #c.SetGrid()
  pad1 = TPad('pad1','pad1',0,0.25,1,1)
  pad1.SetBottomMargin(0.0)
  pad1.SetFillStyle(0)
  pad1.SetGrid()
  pad1.Draw()
  pad1.cd()
  ROOT.gStyle.SetOptStat(0)
  if(logy):
    pad1.SetLogy();
  if(leftLegend):
      legend1 = LeftLegendMaker()
  else:
      legend1 = LegendMaker()
  tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
  L = [tex1, tex2, tex3]
  TexMaker(L)
  FirstTH1[0].GetYaxis().SetRangeUser(yrange1down,yrange1up)
  FirstTH1[0].GetXaxis().SetRangeUser(xrange1down,xrange1up)
  for i in range(0, len(FirstTH1)):
    FirstTH1[i].SetTitle("")
    #FirstTH1[i].GetYaxis().SetTitle("Events")
    FirstTH1[i].GetYaxis().SetTitleOffset(0.7)
    FirstTH1[i].GetYaxis().SetTitleSize(0.060)
    FirstTH1[i].GetXaxis().SetTitle("")
    FirstTH1[i].GetXaxis().SetLabelSize(0.15)
    FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
    legend1.AddEntry(FirstTH1[i],LegendName[i], "lp")
  
  FirstTH1[0].GetXaxis().SetLabelSize(0.0);
  FirstTH1[0].Draw("AEP")
  #hmcdy.Draw("sames")
  
  for i in range(1, len(FirstTH1)):
    FirstTH1[i].Draw("ep sames")  
    
  L[0].Draw()
  L[1].Draw()
  L[2].Draw()
  legend1.Draw()
  
  c.cd()
  pad2 = TPad("pad2", "pad2",0.0,0.0,1.0,0.249)
  pad2.SetTopMargin(0.0)
  pad2.SetBottomMargin(0.3)
  pad2.SetFillStyle(0)
  #pad2.SetGrid()
  pad2.Draw()
  pad2.cd()
  ROOT.gStyle.SetOptStat(0)
    
  hJESError = JESError(hdatady, hmcdy, hmcdyup, hmcdydown)
  hStatError = StatError(hdatady, hmcdy, hmcdyup, hmcdydown)
  hJESError.SetFillColor(kRed)
  hJESError.SetFillStyle(3003)
  hStatError.SetFillColor(kGreen)
  hStatError.SetFillStyle(3003)
  
  xaxisTitleTGraph = FirstTH1[0].GetXaxis().GetTitle()
  
  h4.GetXaxis().SetLabelSize(0.1)
  h4.GetXaxis().SetTitle(xaxisTitleTGraph)
  h4.GetXaxis().SetRangeUser(xrange1down,xrange1up)
  h4.GetXaxis().SetTitleSize(0.16)
  h4.GetXaxis().SetTitleOffset(0.6)
  if(xaxisTitle == 'qt'):
    h4.GetXaxis().SetTitle('qt [GeV]')
  elif(xaxisTitle =='nvtx'):
    h4.GetXaxis().SetTitle('nvertex')
  elif(xaxisTitle == 'sumEt'):
    h4.GetXaxis().SetTitle('sumEt [TeV]')
  elif(xaxisTitle == 'sumEtRaw'):
    h4.GetXaxis().SetTitle('sumEtRaw [TeV]')
  else:
    h4.GetXaxis().SetTitle(xaxisTitle)
  h4.GetYaxis().SetLabelSize(0.12)
  h4.GetYaxis().SetTitle("#frac{Data}{PseudoData}")
  h4.GetYaxis().SetNdivisions(5)
  h4.SetLineColor(4)
  h4.SetMarkerColor(4)
  h4.SetMarkerStyle(kFullDotLarge)
  h4.GetYaxis().CenterTitle()
  h4.SetTitle("")
  h4.GetYaxis().SetTitleOffset(0.36)
  h4.GetYaxis().SetTitleSize(0.12)
  h4.GetYaxis().SetRangeUser(0.5, 1.5)
  h4.Draw("AEP")
  hJESError.Draw("2 same")
  hStatError.Draw("2 same")
  
  legratio = TLegend(0.70,0.75,0.90,0.95)
  legratio.SetFillColor(0);
  legratio.SetBorderSize(0);
  legratio.SetLineColor(kWhite)
  legratio.SetFillColor(kWhite)
  legratio.AddEntry(hStatError, "Stat","f")
  legratio.AddEntry(hJESError, "JES+Stat","f")
  legratio.Draw()
  
  line = MakeLine(xrange1down,yline1low,xrange1up,yline1up)
  if(drawline):
    line.Draw()
  SaveFile(c, CanvasName)
  c.Clear()
  
 
 
## needed for two ratio axis size ## 
 
 
def TwoRatioAxisSize(histList):
    for j in range(0, len(histList)):
        histList[j].GetXaxis().SetLabelSize(0.12)
        histList[j].GetXaxis().SetTitleSize(0.16)
        histList[j].GetXaxis().SetTitleOffset(0.78)
        
        histList[j].GetYaxis().SetLabelSize(0.12)
        histList[j].GetYaxis().SetNdivisions(5)
        histList[j].GetYaxis().SetTitleOffset(0.28)
        histList[j].GetYaxis().SetTitleSize(0.18)
        histList[j].GetYaxis().SetRangeUser(0.4, 1.6)
        histList[j].GetYaxis().CenterTitle()
        histList[j].GetYaxis().SetTitle("#frac{E_{Tracker}^{Face}}{E_{IP}}")
        
        histList[j].SetTitle("")
        #histList[j].SetMarkerStyle(kFullTriangleDown)
  
  #histList[1].GetXaxis().SetLabelSize(0.1)
  #histList[1].GetXaxis().SetTitleSize(0.16)
  
  #histList[1].GetYaxis().SetLabelSize(0.13)
  #histList[1].GetYaxis().SetNdivisions(5)
  #histList[1].GetYaxis().SetTitleOffset(0.25)
  #histList[1].GetYaxis().SetTitleSize(0.18)
  #histList[1].GetYaxis().SetRangeUser(0.5, 4.0)
  
  #histList[1].SetMarkerStyle(kFullTriangleUp)
  #histList[1].GetYaxis().CenterTitle()
  #histList[1].SetTitle("")
  
  
  
  
  
def DrawHistsRatioTwo(FirstTH1, LegendName, PlotColor, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, h2, h3, yline1low, yline1up, drawline=False, logy=False, doSumw2=False, TeVTag=False, latexName = '', latexName2='', latexName3=''):
   Tex = MakeLatex(0.90,0.62,latexName)
   Tex2 = MakeLatex(0.90, 0.55, latexName2)
   Tex3 = MakeLatex(0.90, 0.48, latexName3)
   c = TCanvas("c","c",900, 900)
   c.cd()
   pad1 = TPad("pad1","pad1",0,0.25,1,1)
   pad1.SetBottomMargin(0.0)
   pad1.SetFillStyle(0)
   pad1.Draw()
   pad1.cd()
   ROOT.gStyle.SetOptStat(0)
   if(logy):
     pad1.SetLogy()
     
   line = MakeLine(xrange1down,yline1low,xrange1up,yline1up)
   legend1 = LegendMaker()
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L)
   WaterMark = TexWaterMark('Preliminary')
   FirstTH1[0].GetYaxis().SetRangeUser(yrange1down,yrange1up)
   FirstTH1[0].GetXaxis().SetRangeUser(xrange1down,xrange1up)
   for i in range(0,len(FirstTH1)):
     FirstTH1[i].SetTitle("")
     FirstTH1[i].GetYaxis().SetTitle(FirstTH1[i].GetYaxis().GetTitle())
     FirstTH1[i].GetYaxis().SetTitleOffset(0.7)
     FirstTH1[i].GetYaxis().SetTitleSize(0.060)
     FirstTH1[i].GetXaxis().SetTitle("")
     FirstTH1[i].GetXaxis().SetLabelSize(0.15)
     FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
     
     
     FirstTH1[i].SetFillColor(0)
     
     if(i==0 or i==2):
         FirstTH1[i].SetLineStyle(1)
         FirstTH1[i].SetLineWidth(1)
     else:
         FirstTH1[i].SetLineStyle(2)
         FirstTH1[i].SetLineWidth(2)
         
     integralValue = FirstTH1[i].Integral()
     
     if(i==0):
         legend1.AddEntry(FirstTH1[i],LegendName[i]+" (Integral: "+str(round(integralValue, 1))+")", "l")
     else:
        legend1.AddEntry(FirstTH1[i],LegendName[i]+" (Integral: "+str(round(integralValue, 1))+")", "l")
     
   FirstTH1[0].GetXaxis().SetLabelSize(0.0)
   
   FirstTH1[0].Draw("hist sames")
   WaterMark.Draw("sames")
   c.Modified()
   c.Update()
   #FirstTH1[0].SetFillColor(P)
   FirstTH1[0].Draw("hist")
   for i in range(1,len(FirstTH1)):
     FirstTH1[i].SetFillColor(0)  
     FirstTH1[i].Draw("hist sames")  
     c.Modified()
     c.Update()
   #FirstTH1[0].Draw("ep sames")
   ROOT.gPad.RedrawAxis();
   ROOT.gPad.SetTickx();
   ROOT.gPad.SetTicky();
   
   if(TeVTag):
       TexTeV = TLatex(0.892,0.914,"13 TeV")#36.075 fb^{-1} (13 TeV)
       TexTeV.SetNDC()
       TexTeV.SetTextAlign(31)
       TexTeV.SetTextFont(42)
       TexTeV.SetTextSize(0.037)
       TexTeV.SetLineWidth(2) 
       TexTeV.Draw()
   else:
       L[0].Draw()
   L[1].Draw()
   L[2].Draw()
   legend1.Draw()
   Tex.Draw()
   Tex2.Draw()
   Tex3.Draw()
   
   
   c.cd()
   pad2 = TPad("pad2", "pad2",0.0,0.0,1.0,0.249)
   pad2.SetTopMargin(0.0)
   pad2.SetBottomMargin(0.3)
   pad2.SetFillStyle(0)
   pad2.Draw()
   pad2.cd()
   ROOT.gStyle.SetOptStat(0)
   
   if(doSumw2):
       for i in range(0, len(FirstTH1)):
           FirstTH1[i].Sumw2()
           
   h2.Divide(FirstTH1[0], FirstTH1[1])
   h2.GetXaxis().SetRangeUser(xrange1down,xrange1up)
   h2.SetLineColor(PlotColor[1])
   h2.SetMarkerColor(PlotColor[1])
   h2.SetLineStyle(1)
   h2.SetFillColor(0)
   
   
   h3.Divide(FirstTH1[2], FirstTH1[3])
   h3.GetXaxis().SetRangeUser(xrange1down,xrange1up)
   h3.SetLineColor(PlotColor[2])
   h3.SetMarkerColor(PlotColor[2])
   h3.SetLineStyle(2)
   h3.SetFillColor(0)
   
   
   legend2 = RatioLegendMaker()
   legend2.AddEntry(h2, h2.GetYaxis().GetTitle(),"l")
   legend2.AddEntry(h3, h3.GetYaxis().GetTitle(),"l")
   
   ratioList = [h2, h3]
   TwoRatioAxisSize(ratioList)
   
   
   ratioList[0].Draw("hist")
   ratioList[1].Draw("hist same")
   legend2.Draw()
   
   if(drawline):
       line.SetLineStyle(3)
       line.Draw()
   SaveFile(c, CanvasName)
   c.Clear()
   
   
   
   
   
   
   
   
   
def DrawHistsRatioFour(FirstTH1, LegendName, PlotColor, xrange1down, xrange1up, yrange1down, yrange1up, CanvasName, h2, h3, yline1low, yline1up, drawline=False, logy=False, doSumw2=False, TeVTag=False, latexName = '', latexName2='', latexName3='', latexName4='', latexName5='', h4=0, h5=0, energyPlots=False):
   Tex = MakeLatex(0.36,0.86,latexName)
   Tex2 = MakeLatex(0.36, 0.79, latexName2)
   Tex3 = MakeLatex(0.36, 0.71, latexName3)
   Tex4 = MakeLatex(0.36, 0.64, latexName4)
   Tex5 = MakeLatex(0.36, 0.57, latexName5)
   c = TCanvas("c","c",900, 900)
   c.cd()
   pad1 = TPad("pad1","pad1",0,0.25,1,1)
   pad1.SetBottomMargin(0.0)
   pad1.SetFillStyle(0)
   pad1.Draw()
   pad1.cd()
   ROOT.gStyle.SetOptStat(0)
   if(logy):
     pad1.SetLogy()
     
   line = MakeLine(xrange1down,yline1low,xrange1up,yline1up)
   legend1 = LegendMaker(energyPlots)
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L)
   WaterMark = TexWaterMark('Preliminary')
   FirstTH1[0].GetYaxis().SetRangeUser(yrange1down,yrange1up)
   FirstTH1[0].GetXaxis().SetRangeUser(xrange1down,xrange1up)
   for i in range(0,len(FirstTH1)):
     FirstTH1[i].SetTitle("")
     FirstTH1[i].GetYaxis().SetTitle(FirstTH1[i].GetYaxis().GetTitle())
     #FirstTH1[i].GetYaxis().SetTitle(FirstTH1[i].GetYaxis().GetTitle())
     FirstTH1[i].GetYaxis().SetTitleOffset(0.7)
     FirstTH1[i].GetYaxis().SetTitleSize(0.060)
     FirstTH1[i].GetXaxis().SetTitle("")
     FirstTH1[i].GetXaxis().SetLabelSize(0.15)
     FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
     integralValue = FirstTH1[i].Integral()
     #if(i==0):
         #legend1.AddEntry(FirstTH1[i],LegendName[i]+" (Integral: "+str(round(integralValue, 1))+")", "f")
     if(i>=2):
        legend1.AddEntry(FirstTH1[i],LegendName[i]+" ("+str(round(integralValue, 1))+")", "l")
        FirstTH1[i].SetLineWidth(2)
     if(i==1):
        legend1.AddEntry(FirstTH1[i],LegendName[i]+" ("+str(round(integralValue, 1))+")", "l")
        FirstTH1[i].SetLineStyle(8)
        FirstTH1[i].SetLineWidth(3)
     
   FirstTH1[0].GetXaxis().SetLabelSize(0.0)
   ### work with the signal
   integralValue = FirstTH1[0].Integral()
   legend1.AddEntry(FirstTH1[0],LegendName[0]+" ("+str(round(integralValue, 1))+")", "f")
   
   FirstTH1[0].Draw("ep sames")
   WaterMark.Draw("sames")
   c.Modified()
   c.Update()
   FirstTH1[0].SetFillColor(PlotColor[0])
   FirstTH1[0].Draw("hist")
   for i in range(1,len(FirstTH1)):
     FirstTH1[i].SetFillColor(0)  
     FirstTH1[i].Draw("hist sames")  
     c.Modified()
     c.Update()
   #FirstTH1[0].Draw("ep sames")
   ROOT.gPad.RedrawAxis();
   ROOT.gPad.SetTickx();
   ROOT.gPad.SetTicky();
   
   if(TeVTag):
       TexTeV = TLatex(0.892,0.914,"13 TeV")#36.075 fb^{-1} (13 TeV)
       TexTeV.SetNDC()
       TexTeV.SetTextAlign(31)
       TexTeV.SetTextFont(42)
       TexTeV.SetTextSize(0.037)
       TexTeV.SetLineWidth(2) 
       TexTeV.Draw()
   else:
       L[0].Draw()
   L[1].Draw()
   L[2].Draw()
   legend1.Draw()
   Tex.Draw()
   Tex2.Draw()
   Tex3.Draw()
   Tex4.Draw()
   Tex5.Draw()
   
   
   c.cd()
   pad2 = TPad("pad2", "pad2",0.0,0.0,1.0,0.249)
   pad2.SetTopMargin(0.0)
   pad2.SetBottomMargin(0.3)
   pad2.SetFillStyle(0)
   pad2.Draw()
   pad2.cd()
   ROOT.gStyle.SetOptStat(0)
   
   if(doSumw2):
       for i in range(0, len(FirstTH1)):
           FirstTH1[i].Sumw2()
           
   h2.Divide(FirstTH1[1], FirstTH1[0])
   h2.GetXaxis().SetRangeUser(xrange1down,xrange1up)
   h2.SetLineColor(PlotColor[1])
   h2.SetMarkerColor(PlotColor[1])
   h2.SetLineStyle(8)
   
   if(len(FirstTH1)>2):
        h3.Divide(FirstTH1[2], FirstTH1[0])
        h3.GetXaxis().SetRangeUser(xrange1down,xrange1up)
        h3.SetLineColor(PlotColor[2])
        h3.SetMarkerColor(PlotColor[2])
   
   if(len(FirstTH1)>3):
        h4.Divide(FirstTH1[3], FirstTH1[0])
        h4.GetXaxis().SetRangeUser(xrange1down,xrange1up)
        h4.SetLineColor(PlotColor[3])
        h4.SetMarkerColor(PlotColor[3])
        h4.SetLineWidth(2)
        h4.SetLineStyle(8)
   
   if(len(FirstTH1)>4):
        h5.Divide(FirstTH1[4], FirstTH1[0])
        h5.GetXaxis().SetRangeUser(xrange1down,xrange1up)
        h5.SetLineColor(PlotColor[4])
        h5.SetMarkerColor(PlotColor[4])
        h5.SetLineWidth(2)
        h5.SetLineStyle(8)
   
   
   legend2 = RatioLegendMaker()
   legend2.AddEntry(h2, h2.GetYaxis().GetTitle(),"lp")
   legend2.AddEntry(h3, h3.GetYaxis().GetTitle(),"lp")
   
   ratioList = [h2, h3, h4, h5]
   TwoRatioAxisSize(ratioList)
   
   
   ratioList[0].Draw("hist")
   for i in range(1, len(ratioList)):
        ratioList[i].Draw("hist same")
        
   #legend2.Draw()
   
   if(drawline):
       line.Draw()
   SaveFile(c, CanvasName)
   c.Clear()
   
   
   
def DrawFunctions(FirstTH1, LegendName, PlotColor, xAxisName, yAxisName, xAxisLowRange, xAxisHighRange, yAxisLowRange, yAxisHighRange, CanvasName, latexName, logy=False, doAtlas=False):
   Tex = MakeLatex(0.64,0.35,latexName)
   WaterMark = TexWaterMark('Preliminary')
   c = TCanvas("c","c",600, 600)
   ROOT.gStyle.SetOptStat(0)
   c.cd()
   if(logy):
     c.SetLogy()
   ROOT.gPad.SetLeftMargin(0.15)
   legend1 = LegendMaker()
   tex1 = TLatex(); tex2 = TLatex(); tex3 = TLatex()
   L = [tex1, tex2, tex3]
   TexMaker(L, doAtlas)
   FirstTH1[0].GetYaxis().SetRangeUser(yAxisLowRange,yAxisHighRange)
   for i in range(0, len(FirstTH1)):
     FirstTH1[i] = AxisLabelEtc(FirstTH1[i], yAxisName, xAxisName)
     FirstTH1[i] = SetHistColorEtc(FirstTH1[i], PlotColor[i])
     legend1.AddEntry(FirstTH1[i],LegendName[i], "lp")
   
   FirstTH1[0].GetXaxis().SetRangeUser(xAxisLowRange,xAxisHighRange)
   ROOT.gPad.SetTickx()
   ROOT.gPad.SetTicky()
   FirstTH1[0].Draw("l")
   #WaterMark.Draw("sames")
   for i in range(1, len(FirstTH1)):
     FirstTH1[i].Draw("l sames")  
   Tex.Draw("sames")
   #L[0].Draw()
   #L[1].Draw()
   #L[2].Draw()
   legend1.Draw()
   
   SaveFile(c, CanvasName)
   return [c,L,legend1]
    

  




