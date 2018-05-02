import ROOT as rt
from time import sleep
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pickle

#helper functions:

def AngleCorr(angle,border):
        """returns input angle inside period indicated by border"""
        if angle > border:
                return angle-2*border
        elif angle < -border:
                return angle+2*border
        else:
                return angle

def PolarPhi(x,y):
        """returns angle phi according to polar/spherical convention"""
        r = np.sqrt(x**2+y**2)
        if x>0:
                return np.arctan(y/x)
        if x<0 and y >= 0:
                return np.arctan(y/x) + np.pi
        if x<0 and y < 0:
                return np.arctan(y/x) - np.pi
        if x==0 and y > 0:
                return np.pi/2
        if x==0 and y < 0:
                return -np.pi/2

def Theta(x,y,z):
        """returns angle theta according to spherical coordinate convention"""
        return np.pi/2 - np.arctan(z/np.sqrt(x**2+y**2))

def normalize(vector):
        """normalizes input vector to unity"""
        r2 = 0
        for i in vector:
                r2 += i**2
        return vector/np.sqrt(r2)

def Eta(theta):
        """converts the angle theta to the pseudorapidity eta"""
        return -np.log(np.tan(theta/2))

def DeltaR(theta1,theta2,phi1,phi2):
        """returns deltaR according to particle physics convention"""
        deta = Eta(theta1)-Eta(theta2)
        dphi = AngleCorr(phi1-phi2,np.pi)
        return np.sqrt(deta**2 + dphi**2)

def DeltaR_eta(eta1,eta2,phi1,phi2):
        """returns deltaR according to particle physics convention using eta"""
        deta = eta1 - eta2
        dphi = AngleCorr(phi1-phi2,np.pi)
        return np.sqrt(deta**2 + dphi**2)

def ShiftedThetaPhi(x,y,z,dx,dy,dz):
        """returns a tuple containing angles theta and phi with respect to a point (dx,dy,dz)"""
        return (Theta(x-dx,y-dy,z-dz),PolarPhi(x-dx,y-dy))

def TrajectoryLength(theta,v_p):
        if theta <= 0.079*np.pi or theta >= 0.921*np.pi: #Estimating necessary length for trajectory
                return 60/np.sqrt(v_p[0]**2+v_p[1]**2+v_p[2]**2)
        elif theta >= 0.224*np.pi and theta <= 0.776*np.pi:
                return 25/np.sqrt(v_p[0]**2+v_p[1]**2+v_p[2]**2)
        else:
                return 45/np.sqrt(v_p[0]**2+v_p[1]**2+v_p[2]**2)

def Initialize3DPlot(title, xlabel, ylabel, zlabel, grid=None):
        '''initialize 3D-plot'''
        fig = plt.figure(title)
        fig.clf()
        ax = Axes3D(fig)
        ax.clear()
        plt.title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_zlabel(zlabel)
        #plot grid of modules for visual reference
        if grid != None:
                ax.scatter(grid[0],grid[1],grid[2],c='k',s=1,linewidths=0.1)
        return ax

def MakeGrid(file, EarlyBreak=0):
        tree = file.Get("demo/tree")
        N = tree.GetEntries()
        DetUnits = []
        for i in xrange(N):
                if i % 50 == 0: print "Working on event " ,i
                if EarlyBreak > 0 and i >= EarlyBreak: break
                tree.GetEntry(i)
                for DetUnit in zip(tree.detUnit_X,tree.detUnit_Y,tree.detUnit_Z):
                        if DetUnit not in DetUnits:
                                DetUnits.append(DetUnit)
        X, Y, Z = [], [], []
        for DetUnit in DetUnits:
                X.append(DetUnit[0])
                Y.append(DetUnit[1])
                Z.append(DetUnit[2])
        with open("Grid.pkl", 'w') as f:
                pickle.dump([X,Y,Z], f)
        return X, Y, Z

def PlotTrajectory(vx,p,ax,T_max,res,col,lwidth,lstyle):
        Tx,Ty,Tz=[],[],[]
        v_t = normalize(np.array([p[0],p[1],p[2]]))
        for t in np.linspace(0,T_max,res):
                Tx.append(vx[0]+t*v_t[0])
                Ty.append(vx[1]+t*v_t[1])
                Tz.append(vx[2]+t*v_t[2])
        ax.plot(xs=Tx,ys=Ty,zs=Tz,color=col,linewidth=lwidth,linestyle=lstyle)#plots all the tracks

def DiscriminatorHist(title,discriminator,Bdata,Backgrounddata,bins,ran,xlabel):
	HistB = rt.TH1D(title,title,bins,ran[0],ran[1])
	HistBackground = rt.TH1D(title,title,bins,ran[0],ran[1])
	for particle in Bdata:
		HistB.Fill(discriminator(particle))
	for particle in Backgrounddata:
		HistBackground.Fill(discriminator(particle))
	canvas = rt.TCanvas('canvas','canvas',600,600)
        canvas.SetTitle(title)
        rt.gStyle.SetOptStat(0)
        HistB.GetXaxis().SetTitle(xlabel)
        HistB.GetYaxis().SetTitle('[a.u.]')
        HistB.GetYaxis().SetTitleOffset(1.5)
        HistB.SetLineColor(2)
        HistBackground.SetLineColor(3)
        HistB.DrawNormalized()
        HistBackground.DrawNormalized('SAME')
        legend = rt.TLegend(0.9,0.9,0.65,0.75)
        legend.AddEntry(HistB,'B-hadrons')
        legend.AddEntry(HistBackground,'Background')
        legend.Draw()
        canvas.SaveAs('discriminants/'+title+'.png')

def Li_Lj_Hist1D(i, j, Bdatalist, Backgrounddata, bins, ran, dR, Difference=False, Abs=False, dR_check=False, Save=False):
        if Difference:
                if Abs:
                        title = "abs_L"+str(i)+"-L"+str(j)+"_dR_"+str(dR)
                        xlabel = "|L"+str(i)+"-L"+str(j)+"|"
                        ran = (0,25)
                else:
                        title = "L"+str(i)+"-L"+str(j)+"_dR_"+str(dR)
                        xlabel = "L"+str(i)+"-L"+str(j)
                        ran = (-25,25)
                bins = ran[1]-ran[0]
        else:
                title = "L"+str(i)+"_L"+str(j)+"_dR_"+str(dR)
                xlabel = "L"+str(i)+"/L"+str(j)
        BG_zero_div = 0
        BG_else = 0
        HistBackground = rt.TH1D("Background",title,bins,ran[0],ran[1])
        HistBlist = []
        B_zero_div_list = []
        B_else_list = []
        for n,Bdata in enumerate(Bdatalist):
                print "working on",Bdata[1]
                B_zero_div = 0
                B_else = 0
                HistBlist.append(rt.TH1D(Bdata[1],title,bins,ran[0],ran[1]))
                HistBlist[n].SetLineColor(n+3)
                for particle in Bdata[0]:
                        if Difference:
				if Abs:
                                	L = abs(Li_Lj_diff(i,j,particle, dR, dR_check=dR_check))
                        	else:
                                	L = Li_Lj_diff(i,j,particle, dR, dR_check=dR_check)
                        else:
                                L = Li_Lj_ratio(i,j,particle, dR, dR_check=dR_check)
                        B_else += 1
                        if L == False:
                                B_zero_div += 1
                        else:
                                HistBlist[n].Fill(L)
                B_zero_div_list.append(B_zero_div)
                B_else_list.append(B_else)
        print "working on background"
        for particle in Backgrounddata:
                if Difference:
                        if Abs:
                                L = abs(Li_Lj_diff(i,j,particle, dR, dR_check=dR_check))
                        else:
                                L = Li_Lj_diff(i,j,particle, dR, dR_check=dR_check)
                else:
                        L = Li_Lj_ratio(i,j,particle, dR, dR_check=dR_check)
                BG_else += 1
                if L == False:
                        BG_zero_div += 1
                else:
                        HistBackground.Fill(L)
        canvas = rt.TCanvas('canvas','canvas',600,600)
        canvas.SetTitle(title)
        rt.gStyle.SetOptStat(0)
        HistBlist[0].GetXaxis().SetTitle(xlabel)
        HistBlist[0].GetYaxis().SetTitle('[a.u.]')
        HistBlist[0].GetYaxis().SetTitleOffset(1.5)
        HistBackground.SetLineColor(2)
        HistBlist[0].DrawNormalized()
        if len(HistBlist)>1:
                for n,HistB in enumerate(HistBlist):
                        if n>0: HistB.DrawNormalized('SAME')
        HistBackground.DrawNormalized('SAME')
        legend = rt.TLegend(0.9,0.9,0.65,0.75)
        for n,HistB in enumerate(HistBlist):
                legend.AddEntry(HistB,Bdatalist[n][1])
        legend.AddEntry(HistBackground,'Background')
        legend.Draw()
        print "Background sample has",BG_zero_div,"out of",BG_else,"zero-divsion entries ->",100*BG_zero_div/float(BG_else),"%"
        for n,Bdata in enumerate(Bdatalist):
                print Bdata[1],"sample has",B_zero_div_list[n],"out of",B_else_list[n],"zero-divsion entries ->",100*B_zero_div_list[n]/float(B_else_list[n]),"%"
        if Save:
                Tfile = rt.TFile("histogram_files/"+title+"histograms.root","recreate")
                for Hist in HistBlist: Hist.Write()
                HistBackground.Write()
                canvas.SaveAs('discriminants/'+title+'.png')
                #print 'saved as discriminants/'+title+'DR'+str(dR)+'.png'
        #sleep(10)

def Li_Lj_Hist2D(title,i,j,data,ran,dR,Save=False):
	bins = ran[1] - ran[0]
	Hist = rt.TH2I(title,title,bins,ran[0],ran[1],bins,ran[0],ran[1])
	for particle in data:
		Li,Lj = 0,0
		for cluster in particle:
			if cluster[0][2] == i: Li += 1
			if cluster[0][2] == j: Lj += 1
		Hist.Fill(Li,Lj)
	canvas = rt.TCanvas('canvas','canvas',600,600)
        canvas.SetTitle(title)
        rt.gStyle.SetOptStat(0)
        Hist.GetXaxis().SetTitle('L'+str(i))
        Hist.GetYaxis().SetTitle('L'+str(j))
        Hist.GetYaxis().SetTitleOffset(1.5)
        Hist.Draw("COLZ")
        if Save: 
		canvas.SaveAs('discriminants/L'+str(i)+'_L'+str(j)+'Hist_2D_DR'+str(dR)+title+'.png')
		print 'saved as discriminants/L'+str(i)+'_L'+str(j)+'Hist_2D_DR'+str(dR)+title+'.png'
	sleep(5)

def Li_Lj_ratio(i,j,ParticleData, dR, dR_check=False):
        Li,Lj = 0,0
        for cluster in ParticleData:
		if dR_check and cluster[0][8] >= dR: continue
                if cluster[0][2] == i:
                        Li += 1
                if cluster[0][2] == j:
                        Lj += 1
        if Lj != 0:
                return Li/float(Lj)
        else:
                return False

def Li_Lj_diff(i,j,ParticleData, dR, dR_check=False):
        Li,Lj = 0,0
        for cluster in ParticleData:
		if dR_check and cluster[0][8] >= dR: continue
                if cluster[0][2] == i:
                        Li += 1
                if cluster[0][2] == j:
                        Lj += 1
        return Li-Lj
	
def Layer_Hist(title,data,dR,minR=0,minPT=0,Save=False):
        L1,L2,L3,L4, = 0,0,0,0
        if minR>0 or minPT>0:
                add = ' with R>'+str(minR)+' and pT>'+str(minPT)
        else:
                add = ''
        for particle in data:
                for cluster in particle:
                        if cluster[0][7] >= minR and cluster[0][5] >= minPT:
                                if cluster[0][2] == 1: L1 += 1
                                if cluster[0][2] == 2: L2 += 1
                                if cluster[0][2] == 3: L3 += 1
                                if cluster[0][2] == 4: L4 += 1
        fig2, ax2 = plt.subplots(1,2,figsize=(9,5))
        fig2.suptitle('Hit Clusters per Layer'+add+' inside dR<'+str(dR)+' on '+title+' sample')
        ax2[0].bar([0.5,1.5,2.5,3.5],[L1,L2,L3,L4],align='center')
        ax2[0].set_ylabel('Clusters')
        ax2[0].set_xticks([0.5,1.5,2.5,3.5])
        ax2[0].set_xticklabels(['L1','L2','L3','L4'])
        ax2[1].bar([0.5,1.5,2.5],[L2/float(L1),L3/float(L2),L4/float(L3)],align='center')
        ax2[1].set_ylabel('[a.u.]')
        ax2[1].set_xticks([0.5,1.5,2.5])
        ax2[1].set_xticklabels(['L2/L1','L3/L2','L4/L3'])
        plt.tight_layout(pad=2.0,w_pad=0.5,h_pad=0.5)
        if Save:
                if minR>0 or minPT>0:
                        fig2.savefig('HitsPerLayerDR'+str(dR)+'R'+str(minR)+'PT'+str(minPT)+title+'.png')
                        print 'saved as HitsPerLayerDR'+str(dR)+'R'+str(minR)+'PT'+str(minPT)+title+'.png'
                else:
                        fig2.savefig('HitsPerLayerDR'+str(dR)+title+'.png')
                        print 'saved as HitsPerLayerDR'+str(dR)+title+'.png'
        plt.show()

def Layer_Hist2(title,data,dR,minR=0,minPT=0,Save=False):
        L1,L2,L3,L4, = 0,0,0,0
        if minR>0 or minPT>0:
                add = ' (R>'+str(minR)+', pT>'+str(minPT)+')'
        else:
                add = ''
        for particle in data:
                for cluster in particle:
                        if cluster[0][7] >= minR and cluster[0][5] >= minPT:
                                if cluster[0][2] == 1: L1 += 1
                                if cluster[0][2] == 2: L2 += 1
                                if cluster[0][2] == 3: L3 += 1
                                if cluster[0][2] == 4: L4 += 1
        fig2, ax2 = plt.subplots(1,1,figsize=(5,5))
        fig2.suptitle(title+add+' in dR<'+str(dR))
        ax2.bar([0.5,1.5,2.5,3.5],[L1,L2,L3,L4],align='center')
        ax2.set_ylabel('Clusters')
        ax2.set_xticks([0.5,1.5,2.5,3.5])
        ax2.set_xticklabels(['L1','L2','L3','L4'])
        #plt.tight_layout(pad=2.0,w_pad=0.5,h_pad=0.5)
        if Save:
                if minR>0 or minPT>0:
                        fig2.savefig('HitsPerLayerDR'+str(dR)+'R'+str(minR)+'PT'+str(minPT)+title+'.png')
                        print 'saved as HitsPerLayerDR'+str(dR)+'R'+str(minR)+'PT'+str(minPT)+title+'.png'
                else:
                        fig2.savefig('HitsPerLayerDR'+str(dR)+title+'.png')
                        print 'saved as HitsPerLayerDR'+str(dR)+title+'.png'
        plt.show()


def Histogramize(Histograms, ran, title, xlabel, ylabel, Save=False,Normalize=True, t_sleep=0):
        canvas = rt.TCanvas('canvas','canvas',600,600)
        canvas.SetTitle(title)
        if len(Histograms) > 1:
                rt.gStyle.SetOptStat(0)
                legend = rt.TLegend(0.9,0.9,0.65,0.75)
        for nr, Histogram in enumerate(Histograms):
                Histogram[0].SetLineColor(nr+2)
                if nr == 0:
                        Histogram[0].GetXaxis().SetTitle(xlabel)
                        Histogram[0].GetYaxis().SetTitle(ylabel)
                        Histogram[0].GetYaxis().SetTitleOffset(1.5)
                        if Normalize:
                                Histogram[0].DrawNormalized()
                        else:
                                Histogram[0].Draw()
                else:
                        if Normalize:
                                Histogram[0].DrawNormalized("SAME")
                        else:
                                Histogram[0].Draw("SAME")
                if len(Histograms)>1:
                        legend.AddEntry(Histogram[0],Histogram[1])
        if len(Histograms)>1: legend.Draw()
        if Save: canvas.SaveAs(title+".png")
        sleep(t_sleep)

def LoadHistogram(Histogram, filename):
        with open(filename,) as f:
                R_List = pickle.load(f)
        for entry in R_List:
                Histogram.Fill(entry)
        return Histogram

def FillSeparateLayerHist(L1_Hist,L2_Hist,L3_Hist,L4_Hist,title,data,ran,minPT=0):
        for particle in data:
                L1,L2,L3,L4, = 0,0,0,0
                for n,cluster in enumerate(particle):
			if cluster[0][5] >= minPT:
                        	if cluster[0][2] == 1: L1 += 1
                      		if cluster[0][2] == 2: L2 += 1
                        	if cluster[0][2] == 3: L3 += 1
                        	if cluster[0][2] == 4: L4 += 1
                L1_Hist.Fill(L1)
                L2_Hist.Fill(L2)
                L3_Hist.Fill(L3)
                L4_Hist.Fill(L4)

def SeparateLayerHist(datalist, ran, dR,minPT=0, Save=False):
        canvas = rt.TCanvas('canvas','canvas',800,800)
        canvas.SetTitle("Matched clusters per layer in dR<"+str(dR))
        rt.gStyle.SetOptStat(0)
        canvas.Divide(2,2,0,0)
        canvas.GetPad(1).SetTitle("Layer 1")
        canvas.GetPad(2).SetTitle("Layer 2")
        canvas.GetPad(3).SetTitle("Layer 3")
        canvas.GetPad(4).SetTitle("Layer 4")
        Hist_list = [[],[],[],[]]
        for n, data in enumerate(datalist):
                Hist_list[0].append(rt.TH1D("L1."+str(n),"L1",ran[1]-ran[0],ran[0],ran[1]))
                Hist_list[1].append(rt.TH1D("L2."+str(n),"L2",ran[1]-ran[0],ran[0],ran[1]))
                Hist_list[2].append(rt.TH1D("L3."+str(n),"L3",ran[1]-ran[0],ran[0],ran[1]))
                Hist_list[3].append(rt.TH1D("L4."+str(n),"L4",ran[1]-ran[0],ran[0],ran[1]))
                FillSeparateLayerHist(Hist_list[0][n],Hist_list[1][n],Hist_list[2][n],Hist_list[3][n],data[1],data[0],ran,minPT=minPT)
        for l,layer in enumerate(Hist_list):
                canvas.cd(l+1)
                legend = rt.TLegend(0.9,0.9,0.65,0.75)
                for n,Hist in enumerate(layer):
                        Hist.GetXaxis().SetTitle("# clusters")
                        Hist.GetYaxis().SetTitle('[a.u.]')
                        Hist.GetYaxis().SetTitleOffset(1.5)
                        Hist.SetLineColor(n+2)
                        legend.AddEntry(Hist,datalist[n][1])
                        if n==0:
                                Hist.DrawNormalized()
                        else:
                                Hist.DrawNormalized("SAME")
                legend.Draw()
        if Save: canvas.SaveAs("SeparateLayerHistDR"+str(dR)+".png")
        sleep(10)

def ROC_CutBased(title,signal_hist,background_hist,cutregion="above",resolution=100):
	ran = (signal_hist.GetXaxis().GetXmin(),signal_hist.GetXaxis().GetXmax())
	signal_efficiency = []
	background_efficiency = []
	cuts = np.linspace(ran[0],ran[1],resolution)
	bin_ran = (signal_hist.GetXaxis().FindBin(ran[0]),signal_hist.GetXaxis().FindBin(ran[1]))
	if cutregion == "above":
		for cut in cuts:
			bin_cut = signal_hist.GetXaxis().FindBin(cut)
			signal_efficiency.append(signal_hist.Integral(bin_cut,bin_ran[1])/signal_hist.Integral(bin_ran[0],bin_ran[1]))
			background_efficiency.append(1-background_hist.Integral(bin_cut,bin_ran[1])/background_hist.Integral(bin_ran[0],bin_ran[1])) 	
	if cutregion == "below":
		for cut in cuts:
			bin_cut = signal_hist.GetXaxis().FindBin(cut)
			signal_efficiency.append(signal_hist.Integral(bin_ran[0],bin_cut)/signal_hist.Integral(bin_ran[0],bin_ran[1]))
			background_efficiency.append(1-background_hist.Integral(bin_ran[0],bin_cut)/background_hist.Integral(bin_ran[0],bin_ran[1])) 	
	plt.plot(signal_efficiency,background_efficiency,'-',label=title)

def Draw_ROC_curves(discriminants,Resolution=1000,Title='',Save=False):
	fig1 = plt.figure("2TeV-signal")
	plt.clf()
	fig2 = plt.figure("4TeV-signal")
	plt.clf()
	for discriminant in discriminants:	
		print "working on",discriminant
		dis_file = rt.TFile("histogram_files/"+discriminant+"histograms.root","READ")
		signal_hist1 = dis_file.Get("2TeV-signal")
		signal_hist2 = dis_file.Get("4TeV-signal")
		background_hist = dis_file.Get("Background")
		plt.figure("2TeV-signal")
		ROC_CutBased(discriminant,signal_hist1,background_hist,cutregion="above",resolution=Resolution)
		plt.figure("4TeV-signal")
		ROC_CutBased(discriminant,signal_hist2,background_hist,cutregion="above",resolution=Resolution)
	plt.figure("2TeV-signal")
	plt.title("ROC-curve for 2TeV-signal")
	plt.xlabel(r"$\epsilon$_signal")
	plt.ylabel(r"1-$\epsilon$_background")
	plt.legend(loc=3)
	plt.figure("4TeV-signal")
	plt.title("ROC-curve for 4TeV-signal")
	plt.xlabel(r"$\epsilon$_signal")
	plt.ylabel(r"1-$\epsilon$_background")
	plt.legend(loc=3)
	if Save:
		fig1.savefig("ROC/ROC-curves_2TeV_"+Title+".png")
		print "saved as ROC/ROC-curves_2TeV_"+Title+".png"
		fig2.savefig("ROC/ROC-curves_4TeV_"+Title+".png")
		print "saved as ROC/ROC-curves_4TeV_"+Title+".png"
	plt.show()
	

def ClusterMatch(title, file_path, dR, MomentumThreshold, HadronsNotQuarks=False,BG=False, Plot=False, Axes=None, Save=False, dR_dist=False, LayerHist=False,LightVersion=False, EarlyBreak=0):
        """returns unique ID and coordinates of all pixel clusters that lie inside the dR-cone of a b-particle trajectory; optionally it returns also a 3D-plot
                

        Inputs:
                file_path:              path to root file
                dR:                     Delta R region around particle trajectory in which clusters should count as hit
                Momentumthreshold:      momentum threshold for b-particles to be counted
                HadronsNotQuarks:       if set as True, the algorithm will focus on B-hadrons instead of b-quarks
                Plot:                   if set as True, the function will return a 3D-plot
                Axes:                   pre-initialized 3D-plot axes. Only necessary if Plot==True
                Save:                   if set as True, the function will save the data to a .pkl file for later use
                dR_dist:                if set as True, the function will return a histrogram showing the distribution of delta R between pixel clusters and the corresponding trajectory
                EarlyBreak:             non zero integer which denotes the number of events after which the algorithm should stop

        Outputs:
                list of tuples where each contains a tuple with a uniqe identification followed by global cartesian coordinates:((nEvent,nParticle,nLayer,nModule,nCluster,pT,pdgId),x,y,z)"""
	
	print "working on file", file_path
        #file = rt.TFile(file_path,"READ")
	file = rt.TFile.Open(file_path)
        #colorstring for plots:
        hsv = plt.get_cmap('hsv')
        color = hsv(np.linspace(0,1.0,12))
        c = 0 #initialize color index
        res = 50 #trajectory resolution

        # open tree file
        tree = file.Get("demo/tree")
        N = tree.GetEntries()
        HitClusters = []
        L1, L2, L3, L4 = 0,0,0,0
        if dR_dist == True:
                histL1 = rt.TH1D('Clusters(dR) - '+title, 'Clusters(dR) - '+title, 50, 0, dR)
                histL2 = rt.TH1D('Clusters(dR) - '+title, 'Clusters(dR) - '+title, 50, 0, dR)
                histL3 = rt.TH1D('Clusters(dR) - '+title, 'Clusters(dR) - '+title, 50, 0, dR)
                histL4 = rt.TH1D('Clusters(dR) - '+title, 'Clusters(dR) - '+title, 50, 0, dR)

        for i in xrange(N):
                if i % 50 == 0: print "Working on event " ,i
                if EarlyBreak > 0 and i>=EarlyBreak: break
		if Save and i != 0 and i%10000==0:
                	with open("HitClusterDR"+str(dR)+"on"+title+".pkl", 'w') as f:
                        	pickle.dump(HitClusters, f)
			print "saved as HitClusterDR"+str(dR)+"on"+title+".pkl",
                tree.GetEntry(i)
                for j in range(0,tree.nJets):
                        jVector = rt.TLorentzVector()
                        jVector.SetPtEtaPhiM(tree.jet_pt[j],tree.jet_eta[j],tree.jet_phi[j],tree.jet_mass[j])

			previous_ids = []

                        for k in range(0,tree.nGenParticles):
				if BG:
					pdgCriterion = abs(tree.genParticle_pdgId[k]) != 5# and abs(tree.genParticle_pdgId[k]) < 10
					statusCriterion = tree.genParticle_status[k]== 23
                                elif HadronsNotQuarks == False:
                                        pdgCriterion = abs(tree.genParticle_pdgId[k]) == 5
                                        statusCriterion = tree.genParticle_status[k] == 23
                                else:
                                        pdgCriterion = (abs(tree.genParticle_pdgId[k]) > 500 and abs(tree.genParticle_pdgId[k]) < 600) or (abs(tree.genParticle_pdgId[k]) > 5000 and abs(tree.genParticle_pdgId[k]) < 6000) 
                                        statusCriterion = tree.genParticle_status[k] == 2
                                if statusCriterion and pdgCriterion:
                                        pVector = rt.TLorentzVector()
                                        pVector.SetPtEtaPhiM(tree.genParticle_pt[k],tree.genParticle_eta[k], \
                                                tree.genParticle_phi[k],tree.genParticle_mass[k])
                                        delR = jVector.DeltaR(pVector)
                                        if delR < 0.3 and tree.genParticle_pt[k] > MomentumThreshold: #momentum threshold
                                                v_p = normalize(np.array([pVector[0], pVector[1], pVector[2]]))
                                                phi = PolarPhi(v_p[0],v_p[1])
                                                theta = Theta(v_p[0],v_p[1],v_p[2])
							
						escape = False          #filter out identical daughters
                                                if len(previous_ids)>0:
                                                        for prid in previous_ids:
                                                                if (abs(abs(prid[0])-abs(tree.genParticle_pdgId[k])) == 2 and abs(prid[0])>100): 
									escape=True
                                                if escape: continue
						
                                                previous_ids.append((tree.genParticle_pdgId[k],delR))

                                                if Plot == True:
                                                        t_max = TrajectoryLength(theta,v_p)
                                                        PlotTrajectory((tree.genParticle_vx_x[k],tree.genParticle_vx_y[k],tree.genParticle_vx_z[k]),v_p,Axes,t_max,res,color[c],1,'--')

                                                NearClusters = [] #list that will contain for each hit cluster: ((nEvent,nParticle,nModule,nCluster),x,y,z)
                                                for nModule,lenModule in enumerate(tree.nClusters): #finding all clusters inside deltaR<dR
                                                        for nCluster in xrange(0,lenModule):
								ClusterTheta,ClusterPhi = ShiftedThetaPhi(tree.cluster_globalx[nModule][nCluster],\
                                                                        tree.cluster_globaly[nModule][nCluster],tree.cluster_globalz[nModule][nCluster],\
                                                                        tree.genParticle_vx_x[k],tree.genParticle_vx_y[k],tree.genParticle_vx_z[k])
                                                                DR = DeltaR(theta,ClusterTheta,phi,ClusterPhi)
                                                                if DR<dR:
									if LightVersion:
										NearClusters.append(((i,k,tree.detUnit_layer[nModule],nModule,nCluster,tree.genParticle_pt[k],tree.genParticle_pdgId[k],tree.genParticle_decayvx_r[k],DR),''))
 									else:	
                                                                        	NearClusters.append(((i,k,tree.detUnit_layer[nModule],nModule,nCluster,tree.genParticle_pt[k],tree.genParticle_pdgId[k],tree.genParticle_decayvx_r[k],DR),tree.cluster_globalx[nModule][nCluster],tree.cluster_globaly[nModule][nCluster],tree.cluster_globalz[nModule][nCluster]))

                                                                        if dR_dist == True:
                                                                                if tree.detUnit_layer[nModule] == 1:
                                                                                        histL1.Fill(DR)
                                                                                elif tree.detUnit_layer[nModule] == 2:
                                                                                        histL2.Fill(DR)
                                                                                elif tree.detUnit_layer[nModule] == 3:
                                                                                        histL3.Fill(DR)
                                                                                elif tree.detUnit_layer[nModule] == 4:
                                                                                        histL4.Fill(DR)
                                                                        if LayerHist == True:
                                                                                if tree.detUnit_layer[nModule] == 1:
                                                                                        L1 += 1
                                                                                elif tree.detUnit_layer[nModule] == 2:
                                                                                        L2 += 2
                                                                                elif tree.detUnit_layer[nModule] == 3:
                                                                                        L3 += 3
                                                                                elif tree.detUnit_layer[nModule] == 4:
                                                                                        L4 += 4
                                                
                                                if  NearClusters == []:
                                                        break 
                                                else:
                                                        HitClusters.append(NearClusters) #summarizing the hit cluster ID for every event
                                                        if Plot == True:
								X,Y,Z=[],[],[]
								for entry in NearClusters:
                                                                	X.append(entry[1])
                                                                	Y.append(entry[2])
                                                                	Z.append(entry[3])

                                                                Axes.scatter(X,Y,Z,c=color[c],s=9,linewidths=0.1) #plots all the hit clusters
                                                                if c != len(color)-1:
                                                                        c += 1
                                                                else:
                                                                        c = 0
        if BG: 
		particleString = 'background-particles'
	elif HadronsNotQuarks == False:
                particleString = 'b-quarks'
        else:
                particleString = 'B-hadrons'
        print "Total Number of high pt "+particleString+": ", len(HitClusters)
        nHitClusters = 0
        for entry in HitClusters:
                nHitClusters += len(entry)
        print "Total number of clusters hit: ", nHitClusters
        if Plot == True:
                plt.savefig("Trajectories_Clusters.png")
                plt.show()
        if Save == True:
                with open("HitClusterDR"+str(dR)+"on"+title+".pkl", 'w') as f:
                        pickle.dump(HitClusters, f)

        if dR_dist == True:
                c1 = rt.TCanvas('c1','c1',600,600)
                c1.SetTitle('dR distribution')
		rt.gStyle.SetOptStat(0)
                histL1.GetXaxis().SetTitle('dR')
                histL1.GetYaxis().SetTitle('[a.u.]')
                histL1.GetYaxis().SetTitleOffset(1.5)
                histL1.SetLineColor(1)
                histL2.SetLineColor(2)
                histL3.SetLineColor(3)
                histL4.SetLineColor(4)
                histL1.DrawNormalized()
                histL2.DrawNormalized('SAME')
                histL3.DrawNormalized('SAME')
                histL4.DrawNormalized('SAME')
                if BG:
			l1 = rt.TLegend(0.9,0.1,0.65,0.25)
		else:
			l1 = rt.TLegend(0.9,0.9,0.65,0.75)
                l1.AddEntry(histL1,'Layer 1')
                l1.AddEntry(histL2,'Layer 2')
                l1.AddEntry(histL3,'Layer 3')
                l1.AddEntry(histL4,'Layer 4')
                l1.Draw()
                c1.SaveAs('dR_dist/DeltaR-dist'+title+'.png')
		#print 'saved as dR_dist/DeltaR-dist'+title+'.png'

        if LayerHist == True:
                fig2, ax2 = plt.subplots(1,2,figsize=(9,5))
                fig2.suptitle('Hit Clusters per Layer inside dR<'+str(dR)+' on'+title+' sample')
                ax2[0].bar([0.5,1.5,2.5,3.5],[L1,L2,L3,L4],align='center')
                ax2[0].set_ylabel('Clusters')
                ax2[0].set_xticks([0.5,1.5,2.5,3.5])
                ax2[0].set_xticklabels(['L1','L2','L3','L4'])
                ax2[1].bar([0.5,1.5,2.5],[L2/float(L1),L3/float(L2),L4/float(L3)],align='center')
                ax2[1].set_ylabel('[a.u.]')
                ax2[1].set_xticks([0.5,1.5,2.5])
                ax2[1].set_xticklabels(['L2/L1','L3/L2','L4/L3'])
                plt.tight_layout(pad=2.0,w_pad=0.5,h_pad=0.5)
                if Save: 
			fig2.savefig('HitsPerLayer'+title+'.png')
			print 'saved as HitsPerLayer'+title+'.png'
	return HitClusters



if __name__ == '__main__':

	#select global parameters

        #dR = 0.16 #DeltaR threshold for counting clusters
        MomentumThresholdBackground = 350
	MomentumThresholdB = 350

	
	#initialize 3D-plot

        #with open("Grid.pkl",) as f:    #open coordinates of DetUnits for visual reference
        #        Grid = pickle.load(f)
	#ax = Initialize3DPlot('Particle_Trajectories', 'x', 'y', 'z', grid=Grid)

	
	#load pre-processed data

	#with open("HitClusterDR0.05onB-hadrons.pkl",) as f:   
        #	Signal = pickle.load(f)

	#with open("NewHitClusterDR0.04onbackground_particles.pkl",) as f:   
        #	Backgrounddata = pickle.load(f)
	
	print "opening file HitClusterDR0.16on2TeV-Signal.pkl"
	with open("HitClusterDR0.16on2TeV-Signal.pkl",) as f:   
                Signal1 = pickle.load(f)

	print "opening file HitClusterDR0.16on4TeV-Signal.pkl" 	
	with open("HitClusterDR0.16on4TeV-Signal.pkl",) as f:   
                Signal2 = pickle.load(f)

	print "opening file HitClusterDR0.16onBG0.pkl"
	with open("HitClusterDR0.16onBG0.pkl",) as f: 
               Background = pickle.load(f)

	'''
	#Global HitCluster count per layer

	Layer_Hist2('2TeV',Signal1,dR=dR,minR=0, minPT=0, Save=True)
	Layer_Hist2('4TeV',Signal2,dR=dR, minR=0, minPT=0, Save=True)
	Layer_Hist2('Background',Background,dR=dR, minR=0, minPT=0, Save=True)

	Layer_Hist2('2TeV',Signal1,dR=dR,minR=4, minPT=1000, Save=True)
	Layer_Hist2('4TeV',Signal2,dR=dR, minR=4, minPT=1000, Save=True)
	Layer_Hist2('Background',Background,dR=dR, minR=4, minPT=1000, Save=True)
	'''
        
	#Separate HitCluster per Layer Histograms

	#SeparateLayerHist([(Signal1,'2Tev-signal'),(Signal2,'4Tev-signal'),(Background,'Background')], (0,30), dR, minPT=0, Save=True)

        ''' 
	#select file paths	

	SignalFile1 = "/afs/cern.ch/work/t/thaarres/public/bTag_ntracks/ZprimeBBbar_M2000_GENSIMDIGIRECO.root"
	SignalFile2 = "/afs/cern.ch/work/t/thaarres/public/bTag_ntracks/ZprimeBBbar_M4000_GENSIMDIGIRECO.root"
	Background = "/afs/cern.ch/work/t/thaarres/public/bTag_ntracks/QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8.root"
	#Background = "/afs/cern.ch/work/t/thaarres/public/bTag_ntracks/QCD_noPU_3kEv.root"
	'''	
	#additional Background files

	Additional_Background = ['root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_1.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_10.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_11.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_12.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_13.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_15.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_16.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_17.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_19.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_2.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_20.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_21.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_22.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_23.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_24.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_25.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_26.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_27.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_29.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_3.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_31.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_32.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_33.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_34.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_35.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_36.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_37.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_38.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_39.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_4.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_40.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_41.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_42.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_43.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_45.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_47.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_5.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_6.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_7.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_8.root',
'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_QCD_noPU/180430_105754/0000/flatTuple_9.root']
	
	#pre-process data

	#Bdata1 =  ClusterMatch('2TeV-Signal', SignalFile1, dR, MomentumThresholdB, HadronsNotQuarks=True, Plot=False, Axes=None, Save=True, dR_dist = False, LayerHist=False, EarlyBreak=0)
	#Bdata2 =  ClusterMatch('4TeV-Signal', SignalFile2, dR, MomentumThresholdB, HadronsNotQuarks=True, Plot=False, Axes=None, Save=True, dR_dist = False, LayerHist=False, EarlyBreak=0)
	#Backgrounddata = ClusterMatch('Background', Background, dR, MomentumThresholdBackground, HadronsNotQuarks=True, BG=True, Plot=False, Axes=None, Save=True, dR_dist = False, LayerHist=False, EarlyBreak=180000)
	'''
	dR = 0.16
	#need to get GRID-permission first!!!
	for n,entry in enumerate(Additional_Background):
		try:
			ClusterMatch('BG'+str(n), entry, dR, MomentumThresholdBackground, HadronsNotQuarks=True, BG=True, Plot=False, Axes=None, Save=True, dR_dist = False, LayerHist=False, LightVersion=True, EarlyBreak=0)
		except:
			continue
	'''
	'''
	#HugeBackground = ClusterMatch('Background', Additional_Background[0], dR, MomentumThresholdBackground, HadronsNotQuarks=True, BG=True, Plot=True, Axes=ax, Save=False, dR_dist = False, LayerHist=False, EarlyBreak=500)
	
	
	#Plot 2D-Histograms for every combination of Li,Lj and for every available sample
	
	Li_Lj_Hist2D('2TeV-Signal',1,2,Signal1,(0,35),dR,Save=True)
	Li_Lj_Hist2D('2TeV-Signal',1,3,Signal1,(0,35),dR,Save=True)
	Li_Lj_Hist2D('2TeV-Signal',1,4,Signal1,(0,35),dR,Save=True)
	Li_Lj_Hist2D('2TeV-Signal',2,3,Signal1,(0,35),dR,Save=True)
	Li_Lj_Hist2D('2TeV-Signal',2,4,Signal1,(0,35),dR,Save=True)
	Li_Lj_Hist2D('2TeV-Signal',3,4,Signal1,(0,35),dR,Save=True)
	Li_Lj_Hist2D('4TeV-Signal',1,2,Signal2,(0,35),dR,Save=True)
	Li_Lj_Hist2D('4TeV-Signal',1,3,Signal2,(0,35),dR,Save=True)
	Li_Lj_Hist2D('4TeV-Signal',1,4,Signal2,(0,35),dR,Save=True)
	Li_Lj_Hist2D('4TeV-Signal',2,3,Signal2,(0,35),dR,Save=True)
	Li_Lj_Hist2D('4TeV-Signal',2,4,Signal2,(0,35),dR,Save=True)
	Li_Lj_Hist2D('4TeV-Signal',3,4,Signal2,(0,35),dR,Save=True)
	Li_Lj_Hist2D('background',1,2,Background,(0,35),dR,Save=True)
	Li_Lj_Hist2D('background',1,3,Background,(0,35),dR,Save=True)
	Li_Lj_Hist2D('background',1,4,Background,(0,35),dR,Save=True)
	Li_Lj_Hist2D('background',2,3,Background,(0,35),dR,Save=True)
	Li_Lj_Hist2D('background',2,4,Background,(0,35),dR,Save=True)
	Li_Lj_Hist2D('background',3,4,Background,(0,35),dR,Save=True)
	'''	
	'''
	#Plot 1D-Histogram for different combinations of Li,Lj
	
	Li_Lj_Hist1D(2, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR, Save=True)
	Li_Lj_Hist1D(3, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR, Save=True)
	Li_Lj_Hist1D(3, 2, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR, Save=True)
	Li_Lj_Hist1D(4, 2, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR, Save=True)
	Li_Lj_Hist1D(4, 3, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR, Save=True)
	'''
	
	#Plot 1D-Histogram for L4_L1 and different dR
	
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.01,Difference=True, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.02,Difference=True, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.04,Difference=True, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.08,Difference=True, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.16,Difference=True, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.01,Difference=True, Abs=True, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.02,Difference=True, Abs=True, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.04,Difference=True, Abs=True, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.08,Difference=True, Abs=True, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.16,Difference=True, Abs=True, dR_check=True, Save=True)


	#hist_files = ['L2_L1','L3_L1','L4_L1','L3_L2','L4_L2','L4_L3']
	#hist_files = ['L4_L1_dR_0.01', 'L4_L1_dR_0.02', 'L4_L1_dR_0.04', 'L4_L1_dR_0.08', 'L4_L1_dR_0.16'] 
	hist_files = ['L4-L1_dR_0.01', 'L4-L1_dR_0.02', 'L4-L1_dR_0.04', 'L4-L1_dR_0.08', 'L4-L1_dR_0.16'] 
	Draw_ROC_curves(hist_files, Resolution = 50000,Title='L4-L1',Save=True)
	hist_files = ['abs_L4-L1_dR_0.01', 'abs_L4-L1_dR_0.02', 'abs_L4-L1_dR_0.04', 'abs_L4-L1_dR_0.08', 'abs_L4-L1_dR_0.16'] 
	Draw_ROC_curves(hist_files, Resolution = 50000,Title='abs_L4-L1',Save=True)


	#print "opening files"
	#hist_file = rt.TFile("histograms.root","READ")
	#L2_L1_hist1 = hist_file.Get("L2_L1_2Tev-signal")
	#L2_L1_hist2 = hist_file.Get("L2_L1_4Tev-signal")
	#L2_L1_histBG = hist_file.Get("L2_L1_Background")
	#print "making ROC-curves"
	#plt.figure()
	#plt.clf()

	#ROC_CutBased("L2_L1-2TeV",L2_L1_hist1,L2_L1_histBG,cutregion="above",resolution=200000, Save=True)
	#ROC_CutBased("L3_L1-2TeV",signal_hist,background_hist,cutregion="above",resolution=200000, Save=True)
	#ROC_CutBased("L4_L2-2TeV",signal_hist,background_hist,cutregion="above",resolution=200000, Save=True)
	#ROC_CutBased("L3_L2-2TeV",signal_hist,background_hist,cutregion="above",resolution=200000, Save=True)
	#ROC_CutBased("L4_L2-2TeV",signal_hist,background_hist,cutregion="above",resolution=200000, Save=True)
	#ROC_CutBased("L4_L3-2TeV",signal_hist,background_hist,cutregion="above",resolution=200000, Save=True)
	
	#ROC_CutBased("L2_L1-4TeV",L2_L1_hist2,L2_L1_histBG,cutregion="above",resolution=200000, Save=True)
	#plt.legend()
	#plt.show()
	
	#ROC_CutBased("L2_L1-2TeV-blw",signal_hist,background_hist,cutregion="below",resolution=200000, Save=True)
	#ROC_CutBased("L2_L1-4TeV-blw",signal_hist,background_hist,cutregion="below",resolution=200000, Save=True)
	





	'''
	#manually print out all the clusters hit by background particles

	BG1 = ClusterMatch(,'BG1',Background_file, dR, MomentumThresholdBackground, HadronsNotQuarks=True, BG=True, Plot=True, Axes=ax, Save=False, dR_dist = False, LayerHist=False, EarlyBreak=500)
	for particle in BG1:
		print "particle id =", particle[0][0][6]
		print "particle pT", particle[0][0][5]
		L1, L2, L3, L4, LB = 0,0,0,0,0
		for cluster in particle:
			if cluster[0][2] == 1: L1+=1
			if cluster[0][2] == 2: L2+=1
			if cluster[0][2] == 3: L3+=1
			if cluster[0][2] == 4: L4+=1
			if cluster[0][2] == 0: LB +=1
		print "Hits in L1 =",L1,", L2 =",L2,", L3 =",L3,"L4 =",L4,"Barrel =",LB
	'''
	
	#L2_L1 = lambda ParticleData : Li_Lj(2,1,ParticleData)
	#L3_L1 = lambda ParticleData : Li_Lj(3,1,ParticleData)
        #L4_L1 = lambda ParticleData : Li_Lj(4,1,ParticleData)
	#L3_L2 = lambda ParticleData : Li_Lj(3,2,ParticleData)
	#L4_L2 = lambda ParticleData : Li_Lj(4,2,ParticleData)
	#L4_L3 = lambda ParticleData : Li_Lj(4,3,ParticleData)
	#DiscriminatorHist('L2_L1',L2_L1,Bdata,Backgrounddata,40,(0,6),'L2/L1')
	#DiscriminatorHist('L3_L1',L3_L1,Bdata,Backgrounddata,40,(0,6),'L3/L1')
	#DiscriminatorHist('L4_L1',L4_L1,Bdata,Backgrounddata,40,(0,6),'L4/L1')
	#DiscriminatorHist('L3_L2',L3_L2,Bdata,Backgrounddata,40,(0,6),'L3/L2')
	#DiscriminatorHist('L4_L2',L4_L2,Bdata,Backgrounddata,40,(0,6),'L4/L2')
	#DiscriminatorHist('L4_L3',L4_L3,Bdata,Backgrounddata,40,(0,6),'L4/L3')
	
	
	'''
	#search for ideal status requirement on background particles
	
	tree = Background_file.Get("demo/tree")
        nBE = tree.GetEntries()
	#print "number of B events: ", nBE
	#print "Total Number of high pt B-hadrons ", len(Bdata)
        #nB = 0
        for event in range(nBE):
		if event > 500: break
        	tree.GetEntry(event) 
		for particle in range(0,tree.nGenParticles):
			if not (abs(tree.genParticle_pdgId[particle]) < 7 or abs(tree.genParticle_pdgId[particle])==21 or tree.genParticle_status[particle]>10):
				print "particle Id =", tree.genParticle_pdgId[particle]
				print "particle status =", tree.genParticle_status[particle]
	#print "Total number of clusters hit by B: ", nB
	'''
	


	#Histogram of all statuses

	#statusHist = rt.TH1D('status','status',110,0,109)
	#tree = Background_file.Get("demo/tree")
        #nBGE = tree.GetEntries()
	#for i in xrange(nBGE):
	#	tree.GetEntry(i)
	#	if i%100 == 0: print "working on event", i
	#	for particle in range(0,tree.nGenParticles):
        #               statusHist.Fill(tree.genParticle_status[particle])
	#canvas = rt.TCanvas('canvas','canvas',600,600)
        #canvas.SetTitle("status")
        #statusHist.GetXaxis().SetTitle("status")
        #statusHist.GetYaxis().SetTitle('# particles')
        #statusHist.GetYaxis().SetTitleOffset(1.5)
        #statusHist.SetLineColor(2)
        #statusHist.Draw()
        #canvas.SaveAs('BGstatus.png')
	#print "number of background events: ", nBGE
	#print "Total Number of Background particles ", len(Backgrounddata)
        #nBG = 0
        #for entry in Backgrounddata:
        #        nBG += len(entry)
        #print "Total number of clusters hit by background: ", nBG






