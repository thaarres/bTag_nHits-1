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
	"""returns deltaR according to particle physics convention using theta"""
	deta = Eta(theta1)-Eta(theta2)
	dphi = AngleCorr(phi1-phi2,np.pi)
	return np.sqrt(deta**2 + dphi**2)

def DeltaR_eta(eta1,eta2,phi1,phi2):
	"""returns deltaR according to particle physics convention using eta"""
	deta = eta1-eta2
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

def Initialize3DPlot(title, xlabel, ylabel, zlabel, grid=False, tree=None):
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
	if grid == True:
		tree.GetEntry(0)
		ax.scatter(tree.detUnit_X,tree.detUnit_Y,tree.detUnit_Z,c='k',s=1,linewidths=0.1)
	return ax

def PlotTrajectory(vx,p,ax,T_max,res,col,lwidth,lstyle):
	Tx,Ty,Tz=[],[],[]
	v_t = normalize(np.array([p[0],p[1],p[2]]))
	for t in np.linspace(0,T_max,res):
		Tx.append(vx[0]+t*v_t[0])
		Ty.append(vx[1]+t*v_t[1])
		Tz.append(vx[2]+t*v_t[2])
	ax.plot(xs=Tx,ys=Ty,zs=Tz,color=col,linewidth=lwidth,linestyle=lstyle)#plots all the tracks


#Main function:

def TrackMatch(file, dR, MomentumThreshold, HadronsNotQuarks=False, Plot=False, Save=False, dR_dist = False):
	"""returns unique ID and coordinates of all pixel clusters that lie inside the dR-cone of a b-particle trajectory: optionally it returns also a 3D-plot
		

	Inputs:
		file: 			full root TFile
		dR: 			Delta R region around particle trajectory in which clusters should count as hit
		Momentumthreshold:	momentum threshold for b-particles to be counted
		Plot:			if set as True, the function will return a 3D-plot
		Save:			if set as True, the function will save the data to a .pkl file for later use
		dR_dist:		if set as True, the function will return a histrogram showing the distribution of delta R between pixel clusters and the corresponding trajectory

	Outputs:
		list of tuples where each contains a tuple with a uniqe identification followed by global cartesian coordinates:((nEvent,nParticle,nLayer,nModule,nCluster),x,y,z)"""
	#colorstring for plots:
	#color = ['b','g','r','c','m','y','k','w']
	hsv = plt.get_cmap('hsv')
	color = hsv(np.linspace(0,1.0,12))

	#dxyz_hist = rt.TH1D("dxyz","dxyz",50,0,12)
	#track_vx_r = rt.TH1D("track_vx_r","track_vx_r",40,0,6)
	#track_vx_r_hpt = rt.TH1D("track_vx_r","track_vx_r",40,0,6)
	#x1,x2 = 0,0

	# open tree file
	tree = file.Get("demo/tree")
	N = tree.GetEntries()
	MatchedTracks = []
	#hist = rt.TH1D('B-hadron Pt','B-hadron Pt',40,0.,510.)
		
	if Plot == True:	
		ax = Initialize3DPlot('Particle Trajectories', 'x', 'y', 'z', grid=True, tree=tree)	

		c = 0 #initializing color index
		res = 50 #number of points in trajectory

	for i in xrange(N):
    		if i % 50 == 0: print "Working on event " ,i
		if i >5:
			break
    		tree.GetEntry(i)
		#for track in xrange(0,tree.nTracks):
		#	r = np.sqrt(tree.track_vx_x[track]**2 + tree.track_vx_y[track]**2)
		#	track_vx_r.Fill(r)
		#	x1 += 1
			#if r > 4:
			#	print "Track nr.", track, ", r =", r
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
						eta = Eta(theta)
					
						if Plot == True:                        			
							t_max = TrajectoryLength(theta,v_p)
							PlotTrajectory((tree.genParticle_vx_x[k],tree.genParticle_vx_y[k],tree.genParticle_vx_z[k]),v_p,ax,t_max,res,color[c],1,'--')
                                               	
                        			NearTracks = [] 
						#x3 = 0
						for nTrack in xrange(0,tree.nTracks): #finding all tracks inside deltaR<dR
							VertexTheta,VertexPhi = ShiftedThetaPhi(tree.track_vx_x[nTrack],tree.track_vx_y[nTrack],tree.track_vx_z[nTrack],tree.genParticle_vx_x[k],tree.genParticle_vx_y[k],tree.genParticle_vx_z[k])
							DRT = DeltaR_eta(eta,tree.track_eta[nTrack],phi,tree.track_phi[nTrack])
							dxyz = np.sqrt((tree.track_vx_x[nTrack]-tree.genParticle_vx_x[k])**2 + (tree.track_vx_y[nTrack]-tree.genParticle_vx_y[k])**2 + (tree.track_vx_z[nTrack]-tree.genParticle_vx_z[k])**2)
							if DRT<dR: # and dxyz<0.5:
								#r = np.sqrt(tree.track_vx_x[nTrack]**2 + tree.track_vx_y[nTrack]**2)
								#if r > 4.5:
								#	x2 += 1
								#	x3 += 1
								#track_vx_r_hpt.Fill(r)
								#dxyz_hist.Fill(dxyz)
								NearTracks.append(((i,k,nTrack),(tree.track_vx_x[nTrack],tree.track_vx_y[nTrack],tree.track_vx_z[nTrack]),(tree.track_px[nTrack],tree.track_py[nTrack],tree.track_pz[nTrack])))
						if  NearTracks == []:
							break	
						else:
							MatchedTracks.append(NearTracks) 
                        				if Plot == True:
								for entry in NearTracks:
									PlotTrajectory(entry[1],entry[2],ax,0.7*t_max,res,color[c],0.5,'-')
						if Plot == True:
							if c != len(color)-1:
                        					c += 1
                        				else:
                            					c = 0
						#if x3 != 0:
						#	print "yes! particle", k,"of event", i, "has",x3,"tracks with vertex at r>4.5"
	if HadronsNotQuarks == False:
		particleString = 'b-quarks'
	else:
		particleString = 'B-hadrons'
	print "Total Number of high pt "+particleString+": ", len(MatchedTracks)
	nMatchedTracks = 0
	for entry in MatchedTracks:
		nMatchedTracks += len(entry)
	print "Total number of matched tracks ", nMatchedTracks
	if Plot == True:
		plt.savefig("Trajectories_Tracks.png")
		plt.show()
	if Save == True:
		with open("MatchedTracksDR"+str(dR)+"on"+str(particleString)+".pkl", 'w') as f:
			pickle.dump(MatchedTracks, f)
	'''	
	c1 = rt.TCanvas("c1","c1",600,600)
	dxyz_hist.GetXaxis().SetTitle("DeltaXYZ")
	dxyz_hist.GetYaxis().SetTitle("#")
	dxyz_hist.GetYaxis().SetTitleOffset(1.5)
	c1.SetLogy()
	dxyz_hist.Draw()
	c1.SaveAs("DeltaXYZ-dist.png")
	'''
	'''
	layer1_units=rt.TH1D("whatever","whatever",40,0,6)
 	for unit in xrange(0,len(tree.detUnit_X)):
		if tree.detUnit_layer[unit] == 1:
			layer1_units.Fill(np.sqrt(tree.detUnit_X[unit]**2 + tree.detUnit_Y[unit]**2))
	print np.sqrt(tree.detUnit_X[0]**2+tree.detUnit_Y[0]**2)
	c2 = rt.TCanvas("c1","c1",600,600)
	rt.gStyle.SetOptStat(0)
	track_vx_r.GetXaxis().SetTitle("r")
	track_vx_r.GetYaxis().SetTitle("#")
	track_vx_r.GetYaxis().SetTitleOffset(1.5)
	track_vx_r.SetLineColor(4)
	track_vx_r_hpt.SetLineColor(3)
	layer1_units.SetLineColor(2)
	l1 = rt.TLegend(0.9,0.9,0.65,0.75)
	l1.AddEntry(track_vx_r, "track_vx_r")
	l1.AddEntry(track_vx_r_hpt, "track_vx_r_hpt")
	l1.AddEntry(layer1_units, "layer1_r")
	track_vx_r.DrawNormalized()
	track_vx_r_hpt.DrawNormalized("SAME")
	layer1_units.DrawNormalized("SAME")
	l1.Draw()
	c2.SaveAs("track_vx_r.png")
	print "total number of tracks with r_vx > 4.5:", x1
	print "total number of high pt tracks with r_vx > 4.5:", x2
	'''


if __name__ == '__main__':
	
	file = rt.TFile("/afs/cern.ch/user/m/msommerh/CMSSW_9_3_2/src/bTag/HitAnalyzer/flatTuple.root",'READ')

	dR = 0.1 #DeltaR threshold for counting tracks
	MomentumThreshold = 350

	#MatchedTracks =  TrackMatch(file, dR, MomentumThreshold, HadronsNotQuarks=True, Plot=True, Save=False, dR_dist = False)

	#with open("MatchedTracksDR0.1onB-hadrons.pkl",) as f:	#take data from file instead of running entire function
	#	MatchedTracks = pickle.load(f)
	
	tree = file.Get("demo/tree")
	N = tree.GetEntries()
	ax = Initialize3DPlot('B-Hadron Evolution', 'x', 'y', 'z', grid=False, tree=tree)	
	tree.GetEntry(0)
	#particle = tree.genParticle[0]
	#import pdb; pdb.set_trace()


	def PlotTracks(nEvent,nParticle,color):
		i = nEvent
		k = nParticle
		tree.GetEntry(i)
		#print "Pt =", tree.genParticle_pt[k]
		pVector = rt.TLorentzVector()
        	pVector.SetPtEtaPhiM(tree.genParticle_pt[k],tree.genParticle_eta[k],tree.genParticle_phi[k],tree.genParticle_mass[k])
		v_p = normalize(np.array([pVector[0]/tree.genParticle_mass[k], pVector[1]/tree.genParticle_mass[k], \
        		pVector[2]/tree.genParticle_mass[k]]))
		phi = PolarPhi(v_p[0],v_p[1])
		theta = Theta(v_p[0],v_p[1],v_p[2])
		eta = Eta(theta)
		t_max = TrajectoryLength(theta,v_p)
		PlotTrajectory((tree.genParticle_vx_x[k],tree.genParticle_vx_y[k],tree.genParticle_vx_z[k]),v_p,ax,0.8*t_max,50,'green',1,'--')
		ax.scatter(tree.genParticle_decayvx_x[k],tree.genParticle_decayvx_y[k],tree.genParticle_decayvx_z[k],c='darkolivegreen',s=20,linewidths=0.1)
		print "Event", nEvent, ",Particle", nParticle
		print "B-hadron vertex at", tree.genParticle_vx_x[k],tree.genParticle_vx_y[k],tree.genParticle_vx_z[k]
		print "B-hadron decay-vertex at", tree.genParticle_decayvx_x[k],tree.genParticle_decayvx_y[k],tree.genParticle_decayvx_z[k]
		NearTracks = [] 
		for nTrack in xrange(0,tree.nTracks): #finding all tracks inside deltaR<dR
			VertexTheta,VertexPhi = ShiftedThetaPhi(tree.track_vx_x[nTrack],tree.track_vx_y[nTrack],tree.track_vx_z[nTrack],tree.genParticle_vx_x[k],tree.genParticle_vx_y[k],tree.genParticle_vx_z[k])
			DRT = DeltaR_eta(eta,tree.track_eta[nTrack],phi,tree.track_phi[nTrack])
			if DRT<dR:
				NearTracks.append(((i,k,nTrack),(tree.track_vx_x[nTrack],tree.track_vx_y[nTrack],tree.track_vx_z[nTrack]),(tree.track_px[nTrack],tree.track_py[nTrack],tree.track_pz[nTrack])))
		for nT, entry in enumerate(NearTracks):
			print "Track", nT, ", vertex at", entry[1][0], entry[1][1], entry[1][2]
			PlotTrajectory(entry[1],entry[2],ax,0.7*t_max,50,color,0.5,'-')
		#plt.show()
	
	ClustersX = [-2.1151912212371826, -1.9214503765106201, -3.309929609298706, -3.2997121810913086, -3.297450542449951, -3.3190650939941406, -3.3190648555755615, -3.008465051651001, -5.136828422546387, -4.964303970336914, -4.946408271789551, -4.9160847663879395, -4.912755966186523, -4.901236534118652, -4.467060565948486]
	ClustersY = [4.139054298400879, 4.209573268890381, 6.241778373718262, 6.246328353881836, 6.247335910797119, 6.237710952758789, 6.237710475921631, 6.3760175704956055, 9.06138801574707, 9.144474029541016, 9.153093338012695,9.16769790649414, 9.169301986694336, 9.174849510192871, 9.383955955505371]
	ClustersZ = [9.683456420898438, 9.266300201416016, 11.82179069519043, 11.784566879272461, 11.762925148010254, 11.799983978271484, 11.82470703125, 11.171915054321289, 14.886789321899414, 14.806136131286621, 14.735239028930664, 14.6824369430542, 14.646608352661133, 14.65557861328125, 13.795677185058594]                                                                                                                                         
	PlotTracks(0,5,'red')
	print "particle pddId = ", tree.genParticle_pdgId[5]
	ax.scatter(ClustersX[:2], ClustersY[:2], ClustersZ[:2], c='blue',s=7,linewidths=0.1)
	ax.scatter(ClustersX[2:8], ClustersY[2:8], ClustersZ[2:8], c='blue',s=7,linewidths=0.1)
	ax.scatter(ClustersX[8:], ClustersY[8:], ClustersZ[8:], c='blue',s=7,linewidths=0.1)

	ax.set_xlim3d(-10,1)
	ax.set_ylim3d(-1.15)
	ax.set_zlim3d(0,25)
	plt.show()
	


	'''27, 6: 27, 9: 195, 5'''

	'''Pt = 620.925048828, jet_bTag = 0.845161676407, Hits: 2,6,7'''
	














