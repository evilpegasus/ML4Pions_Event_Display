import uproot
import numpy as np




def get_all_cells(tree):
	rperp = tree['cell_geo_rPerp'].array(library='np')[0]
	cell_eta = tree['cell_geo_eta'].array(library='np')[0]
	cell_theta = 2*np.arctan( np.exp(-cell_eta) )
	cell_phi = tree['cell_geo_phi'].array(library='np')[0]

	#cell_phi[cell_phi > np.pi] = cell_phi[cell_phi > np.pi]-2*np.pi

	cell_x = rperp*np.cos(cell_phi)
	cell_y = rperp*np.sin(cell_phi)
	cell_z = rperp/np.tan(cell_theta)
	cell_positions = np.column_stack([cell_x,cell_y,cell_z])
	cell_geo_ID = tree['cell_geo_ID'].array(library='np')[0]

	id_to_position = {c_id : pos for c_id,pos in zip(cell_geo_ID,cell_positions)}

	return cell_positions,id_to_position