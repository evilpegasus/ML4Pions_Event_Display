import uproot
import numpy as np



def get_event_cells(tree,event_idx):

	
	cluster_cell_ID = tree['cluster_cell_ID'].array(entry_start=event_idx,entry_stop=event_idx+1,library='np')[0]
	cluster_cell_E = tree['cluster_cell_E'].array(entry_start=event_idx,entry_stop=event_idx+1,library='np')[0]

	n_clusters = len(cluster_cell_ID)

	
	cell_E = np.concatenate(cluster_cell_E)
	cell_idx = np.concatenate(cluster_cell_ID)

	return cell_idx,cell_E




# def get_track_parameters(tree,event_idx,track_idx,q=1):

# 	track_variables = ['trackP','trackEta','trackPhi','trackD0','trackZ0']


# 	track_values = {}
# 	for var in track_variables:
# 		track_values[var] =  tree[var].array(entry_start=event_idx,entry_stop=event_idx+1,library='np')[0][track_idx]
	

# 	d0 = track_values['trackD0']
# 	z0 = track_values['trackZ0']
# 	phi = track_values['trackPhi']
# 	theta = 2*np.arctan(np.exp(-track_values['trackEta'])) 
# 	qoverp = q/(track_values['trackP']*1000.0)

# 	return d0,z0,theta,phi,qoverp

def get_track_parameters(tree,event_idx,track_idx,q=1):

	track_variables = ['track_qoverp','track_theta','track_phi','trackD0','trackZ0']


	track_values = {}
	for var in track_variables:
		track_values[var] =  tree[var].array(entry_start=event_idx,entry_stop=event_idx+1,library='np')[0][track_idx]
	

	d0 = track_values['trackD0']
	z0 = track_values['trackZ0']
	phi = track_values['track_phi']
	#if phi > np.pi:
	phi = phi-0.5*np.pi
	theta = track_values['track_theta']
	qoverp = track_values['track_qoverp']

	return d0,z0,theta,phi,qoverp


class TrackHelix(object):
    def __init__(self,d0,z0,theta,phi0,qoverp):
        
        self.theta = theta
        self.phi0 = phi0
        self.d0 = d0
        self.z0 = z0
        self.vx = 0
        self.vy = 0
        self.vz = 0
        self.qoverp = qoverp
        self.Bz = 2.0
        
        self.rho = ( (np.sin(self.theta))/(self.qoverp*self.Bz) )*(1.0/0.299792)
        
    #parametrisation of helix, starting from perigee point
    def x_of_phi(self,phi):
        return self.vx+self.d0*np.cos(self.phi0+np.pi/2.0)+self.rho*(np.cos(phi+np.pi/2)-np.cos(self.phi0+np.pi/2))
    def y_of_phi(self,phi):
        return self.vy+self.d0*np.sin(self.phi0+np.pi/2.0)+self.rho*(np.sin(phi+np.pi/2)-np.sin(self.phi0+np.pi/2))
    def z_of_phi(self,phi):
        return self.vz+self.z0-(self.rho)*(1.0/np.tan(self.theta))*(phi-self.phi0)


def get_track_traj(d0, z0, theta, phi0, qoverp):
    track = TrackHelix(d0,z0,theta,phi0,qoverp)

    phis = np.linspace(phi0,phi0-np.sign(qoverp)*np.pi*0.1,550)
    
    xs = np.array([track.x_of_phi(phi) for phi in phis])
    ys = np.array([track.y_of_phi(phi) for phi in phis])
    zs = np.array([track.z_of_phi(phi) for phi in phis])
    
    rs = np.sqrt( xs**2+ys**2 )

    stop_idx = np.where(rs >2000)[0]

    if len(stop_idx) > 0:
        stop_idx = stop_idx[0]

    else:
        return []

    traj = np.column_stack([xs[:stop_idx],zs[:stop_idx],ys[:stop_idx]])

    return traj