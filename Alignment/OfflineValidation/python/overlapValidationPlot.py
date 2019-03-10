import math
import ROOT
from TkAlStyle import TkAlStyle

def hist(tree_file_name, hist_name,subdet_id,module_direction,overlap_direction):
    f = ROOT.TFile(tree_file_name)
    t = f.Get("analysis/Overlaps")

    h = ROOT.TH1F(hist_name, hist_name, 100, -300, 300)

    h.SetDirectory(0)
	
    for entry in t:
        if not ((t.subdetID == subdet_id)):
            continue
	if module_direction not in ("r" ,"phi", "z"): 
	    raise ValueError("Invalid module direction")
	if module_direction != "z" and subdet_id%2 == 1 and ((abs(t.moduleZ[0] - t.moduleZ[1]) < 1)):
	    continue 
        modulePhi0 = math.atan2(t.moduleY[0], t.moduleX[0]) 
        modulePhi1 = math.atan2(t.moduleY[1], t.moduleX[1])
	phidiff = min(abs(modulePhi0-modulePhi1), math.pi - abs(modulePhi0-modulePhi1))
	moduleR0 = math.sqrt(t.moduleY[0]**2+t.moduleX[0]**2)
	moduleR1 = math.sqrt(t.moduleY[1]**2+t.moduleX[1]**2)
	if module_direction !="r" and subdet_id%2 == 0 and (abs(moduleR0 - moduleR1)<1):
            continue	
	if module_direction !="phi" and ((moduleR0*phidiff)<2):
	    continue
	if overlap_direction == "phi":
             if modulePhi0 > modulePhi1:
                 hitXA = t.hitX[1]
                 hitXB = t.hitX[0]
                 predXA = t.predX[1]
                 predXB = t.predX[0]
                 overlapSignA = t.deltaphi[1]
                 overlapSignB = t.deltaphi[0]
             else:
                 hitXA = t.hitX[0]
                 hitXB = t.hitX[1]
                 predXA = t.predX[0]
                 predXB = t.predX[1]
                 overlapSignA = t.deltaphi[0]
                 overlapSignB = t.deltaphi[1]
	if overlap_direction == "z":
             if t.moduleZ[0] > t.moduleZ[1]:
                 hitXA = t.hitY[1]
                 hitXB = t.hitY[0]
                 predXA = t.predY[1]
                 predXB = t.predY[0]
                 overlapSignA = t.deltaZ[1]
                 overlapSignB = t.deltaZ[0]
             else:
                 hitXA = t.hitY[0]
                 hitXB = t.hitY[1]
                 predXA = t.predY[0]
                 predXB = t.predY[1]
                 overlapSignA = t.deltaZ[0]
                 overlapSignB = t.deltaZ[1]
	if overlap_direction == "r":
             if moduleR0 > moduleR1:
                 hitXA = t.hitY[1]
                 hitXB = t.hitY[0]
                 predXA = t.predY[1]
                 predXB = t.predY[0]
                 overlapSignA = t.deltaR[1]
                 overlapSignB = t.deltaR[0]
             else:
                 hitXA = t.hitY[0]
                 hitXB = t.hitY[1]
                 predXA = t.predY[0]
                 predXB = t.predY[1]
                 overlapSignA = t.deltaR[0]
                 overlapSignB = t.deltaR[1]

        residualA = hitXA - predXA
        residualB = hitXB - predXB
        if overlapSignA < 0:
            residualA *= -1
        if overlapSignB < 0:
            residualB *= -1
        
        A = 10000*(residualA - residualB)
        h.Fill(A)
    return h

def plot(file_name,subdet_id,module_direction,overlap_direction, *filesTitlesColorsStyles):
    hstack = ROOT.THStack("hstack","")
    legend = TkAlStyle.legend(len(filesTitlesColorsStyles), 0.3)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    hs = []
    
    for files, title, color, style in filesTitlesColorsStyles:
        h = hist(files,files.replace("/",""),subdet_id,module_direction,overlap_direction)
        h.SetLineColor(color)
        h.SetLineStyle(style)
        legend.AddEntry(h, title, "l")
        hstack.Add(h)
        hs.append(h)
    
    c = ROOT.TCanvas()
    hstack.Draw("nostack")
    legend.Draw()
    
    hstack.GetXaxis().SetTitle("hit_{A} - pred_{A} - (hit_{B} - pred_{B}) (#mum)")
    hstack.GetYaxis().SetTitle("number of events")
    hstack.GetXaxis().SetNdivisions(404)

    TkAlStyle.drawStandardTitle()
    
    save_as_file_name = file_name

    for ext in "png", "eps", "root", "pdf":
        c.SaveAs(save_as_file_name+"." +ext)

