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

def MakeGrid(file):
	tree = file.Get("demo/tree")
	N = tree.GetEntries()
	DetUnits = []
	for i in xrange(N):
    		if i % 50 == 0: print "Working on event " ,i
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

def Initialize2DPlot(title, xlabel, ylabel, grid=False, tree=None):
        '''initialize 3D-plot'''
        fig = plt.figure(title)
        fig.clf()
        ax = Axes3D(fig)
        ax.clear()
        plt.title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        #ax.set_zlabel(zlabel)
        #plot grid of modules for visual reference
        if grid == True:
                tree.GetEntry(0)
                ax.scatter(tree.detUnit_X,tree.detUnit_Y,c='k',s=1,linewidths=0.1)
        return ax

def PlotTrajectory2D(vx,p,ax,T_max,res,col,lwidth,lstyle):
        Tx,Ty,Tz=[],[],[]
        v_t = normalize(np.array([p[0],p[1],p[2]]))
        for t in np.linspace(0,T_max,res):
                Tx.append(vx[0]+t*v_t[0])
                Ty.append(vx[1]+t*v_t[1])
                Tz.append(vx[2]+t*v_t[2])
        ax.plot(xs=Tx,ys=Ty,color=col,linewidth=lwidth,linestyle=lstyle)#plots all the tracks


#Main function:

def ClusterMatch(file, dR, MomentumThreshold, HadronsNotQuarks=False, Plot=False, Axes=None, Save=False, dR_dist=False, EarlyBreak=0):
	"""returns unique ID and coordinates of all pixel clusters that lie inside the dR-cone of a b-particle trajectory; optionally it returns also a 3D-plot
		

	Inputs:
		file: 			full root TFile
		dR: 			Delta R region around particle trajectory in which clusters should count as hit
		Momentumthreshold:	momentum threshold for b-particles to be counted
		HadronsNotQuarks:	if set as True, the algorithm will focus on B-hadrons instead of b-quarks
		Plot:			if set as True, the function will return a 3D-plot
		Axes:			pre-initialized 3D-plot axes. Only necessary if Plot==True
		Save:			if set as True, the function will save the data to a .pkl file for later use
		dR_dist:		if set as True, the function will return a histrogram showing the distribution of delta R between pixel clusters and the corresponding trajectory
		EarlyBreak:		non zero integer which denotes the number of events after which the algorithm should stop

	Outputs:
		list of tuples where each contains a tuple with a uniqe identification followed by global cartesian coordinates:((nEvent,nParticle,nLayer,nModule,nCluster),x,y,z)"""
	#colorstring for plots:
	#color = ['b','g','r','c','m','y','k','w']
	hsv = plt.get_cmap('hsv')
	color = hsv(np.linspace(0,1.0,12))
	c = 0 #initialize color index
	res = 50 #trajectory resolution

	# open tree file
	tree = file.Get("demo/tree")
	N = tree.GetEntries()
	HitClusters = []
	#hist = rt.TH1D('B-hadron Pt','B-hadron Pt',40,0.,510.)
	if dR_dist == True:
		histL1 = rt.TH1D('DeltaR', 'DeltaR', 30, 0, dR)	
		histL2 = rt.TH1D('DeltaR', 'DeltaR', 30, 0, dR)	
		histL3 = rt.TH1D('DeltaR', 'DeltaR', 30, 0, dR)	
		histL4 = rt.TH1D('DeltaR', 'DeltaR', 30, 0, dR)	
		
	for i in xrange(N):
    		if i % 50 == 0: print "Working on event " ,i
		if EarlyBreak > 0 and i>=EarlyBreak: break
    		tree.GetEntry(i)
    		for j in range(0,tree.nJets):
        		jVector = rt.TLorentzVector()
        		jVector.SetPtEtaPhiM(tree.jet_pt[j],tree.jet_eta[j],tree.jet_phi[j],tree.jet_mass[j])
        		for k in range(0,tree.nGenParticles):
				if HadronsNotQuarks == False:
					pdgCriterion = abs(tree.genParticle_pdgId[k]) == 5
					statusCriterion = tree.genParticle_status[k] == 23
				else:
					pdgCriterion = abs(tree.genParticle_pdgId[k]) > 500 and abs(tree.genParticle_pdgId[k]) < 600
					statusCriterion = tree.genParticle_status[k] == 2 
            			if statusCriterion and pdgCriterion:
                    			pVector = rt.TLorentzVector()
                    			pVector.SetPtEtaPhiM(tree.genParticle_pt[k],tree.genParticle_eta[k], \
                        			tree.genParticle_phi[k],tree.genParticle_mass[k])
                    			delR = jVector.DeltaR(pVector)
					#if delR < 0.3:
					#	hist.Fill(tree.genParticle_pt[k])                       
                    			if delR < 0.3 and tree.genParticle_pt[k] > MomentumThreshold: #momentum threshold
                        			v_p = normalize(np.array([pVector[0]/tree.genParticle_mass[k], pVector[1]/tree.genParticle_mass[k], \
                                   			pVector[2]/tree.genParticle_mass[k]]))
                        			phi = PolarPhi(v_p[0],v_p[1])
						theta = Theta(v_p[0],v_p[1],v_p[2])
						
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
									NearClusters.append(((i,k,tree.detUnit_layer[nModule],nModule,nCluster),tree.cluster_globalx[nModule][nCluster],\
									tree.cluster_globaly[nModule][nCluster],tree.cluster_globalz[nModule][nCluster]))
									if dR_dist == True:
										if tree.detUnit_layer[nModule] == 1:
											histL1.Fill(DR)
										elif tree.detUnit_layer[nModule] == 2:
											histL2.Fill(DR)
										elif tree.detUnit_layer[nModule] == 3:
											histL3.Fill(DR)
										elif tree.detUnit_layer[nModule] == 4:
											histL4.Fill(DR)


						if  NearClusters == []:
							break	
						else:
							HitClusters.append(NearClusters) #summarizing the hit cluster ID for every event
							X,Y,Z=[],[],[]
							for entry in NearClusters:
								X.append(entry[1])
								Y.append(entry[2])
								Z.append(entry[3])
                        				if Plot == True:
								Axes.scatter(X,Y,Z,c=color[c],s=9,linewidths=0.1) #plots all the hit clusters
								if c != len(color)-1:
                        						c += 1
                        					else:
                            						c = 0
	if HadronsNotQuarks == False:
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
		with open("HitClusterDR"+str(dR)+"on"+str(particleString)+".pkl", 'w') as f:
			pickle.dump(HitClusters, f)
	
	if dR_dist == True:
		c1 = rt.TCanvas('c1','c1',600,600)
		c1.SetTitle('DeltaR distribution')
		rt.gStyle.SetOptStat(0)
		histL1.GetXaxis().SetTitle('DeltaR')
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
		l1 = rt.TLegend(0.9,0.9,0.65,0.75)
		l1.AddEntry(histL1,'Layer 1')
		l1.AddEntry(histL2,'Layer 2')
		l1.AddEntry(histL3,'Layer 3')
		l1.AddEntry(histL4,'Layer 4')
		l1.Draw()
		c1.SaveAs('DeltaR-dist.png')	

	return HitClusters
	

if __name__ == '__main__':
	
	file = rt.TFile("flatTuple.root",'READ')

	dR = 0.1 #DeltaR threshold for counting clusters
	MomentumThreshold = 350

	with open("Grid.pkl",) as f:	#open coordinates of DetUnits for visual reference
		Grid = pickle.load(f)

	ax = Initialize3DPlot('Particle Trajectories', 'x', 'y', 'z', grid=Grid)
	HitClusters =  ClusterMatch(file, dR, MomentumThreshold, HadronsNotQuarks=True, Plot=True, Axes=ax, Save=False, dR_dist = False, EarlyBreak=10)

	#with open("HitClusterDR0.1onB-hadrons.pkl",) as f:	#take data from file instead of running entire function
	#	HitClusters = pickle.load(f)
	
	#Count hits per layer for each particle
	
	'''
	#Plot all clusters and all particles

	tree = file.Get("demo/tree")
	N = tree.GetEntries()
	
	ax = Initialize3DPlot('Particle Trajectories', 'x', 'y', 'z', grid=True, tree=tree)
	tree.GetEntry(0)
	X,Y,Z = [],[],[]
	for nModule,lenModule in enumerate(tree.nClusters):
		for nCluster in xrange(0,lenModule):
			X.append(tree.cluster_globalx[nModule][nCluster])
			Y.append(tree.cluster_globaly[nModule][nCluster])
			Z.append(tree.cluster_globalz[nModule][nCluster])

	ax.scatter(X,Y,Z,c='blue',s=9,linewidths=0.1)

	for k in range(0,tree.nGenParticles):
		print tree.genParticle_pdgId[k]
		if abs(tree.genParticle_pdgId[k])<600 and abs(tree.genParticle_pdgId[k])>500:
			color = 'green'
			linestyle = '-'
		else:
			color = 'red'
			linestyle = '--'
		pVector = rt.TLorentzVector()
	        pVector.SetPtEtaPhiM(tree.genParticle_pt[k],tree.genParticle_eta[k], \
              	tree.genParticle_phi[k],tree.genParticle_mass[k])
		v_p = normalize(np.array([pVector[0], pVector[1], pVector[2]]))
	        phi = PolarPhi(v_p[0],v_p[1])
		theta = Theta(v_p[0],v_p[1],v_p[2])
		t_max = TrajectoryLength(theta,v_p)
		PlotTrajectory((tree.genParticle_vx_x[k],tree.genParticle_vx_y[k],tree.genParticle_vx_z[k]),v_p,ax,t_max,50,color,1,linestyle)
	plt.show()
	'''


	'''

	#Count hits per layer for each particle
	
	hsv = plt.get_cmap('hsv')
	color = hsv(np.linspace(0,1.0,12))
	yLim = 30

	plt.figure("Hits per layer")
	plt.clf()
	HitsPerLayer = np.zeros(5)
	ClustersX, ClustersY, ClustersZ = [], [], []
	for n,particle in enumerate(HitClusters):
		if n > 0: break
		for cluster in particle:
			print "event nr", cluster[0][0], ", particle nr.", cluster[0][1]
			HitsPerLayer[cluster[0][2]] +=1
			ClustersX.append(cluster[1])
			ClustersY.append(cluster[2])
			ClustersZ.append(cluster[3])
	plt.plot(range(1,len(HitsPerLayer)),HitsPerLayer[1:],color='b')
	#print "ClustersX:", ClustersX
	#print "ClustersY:", ClustersY
	#print "ClustersZ:", ClustersZ


	plt.plot([1,1],[0,yLim],'k:')
	plt.plot([2,2],[0,yLim],'k:')
	plt.plot([3,3],[0,yLim],'k:')
	plt.xlim(1,4)
	plt.ylim(0,yLim)
	plt.title(r'$\Delta$R < '+str(dR))
	plt.xlabel("layer")
	plt.ylabel("number of clusters")
	#plt.legend(loc=9,ncol=3, prop={'size':8})
	plt.savefig("HitsPerLayer1B.png")
	plt.show()
	'''
	












