import numpy as np
from medpy.io import load, save
from stl import mesh
import skimage
import graph_to_cmgui

filepath = "/eresearch/copdgene/jjoh182/COPDGene_extracted/10088U/COPD1/10088U_INSP_STD_340_COPD1/10088U_INSP_STD_340_COPD1-Airways.mhd"
image, header = load(filepath)

image = np.where(image>0, 1, 0)

def cloud_to_stl(dataarray, filename):
    verts, faces, norms, vals = skimage.measure.marching_cubes(dataarray,step_size=1,allow_degenerate=False)
    # pyvista_pipeline(verts)
    meshobj = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            meshobj.vectors[i][j] = verts[f[j],:]
    meshobj.save(filename+'.stl')
    return meshobj

airway = cloud_to_stl(image,'cmgui_files/airway')

stl_mesh = mesh.Mesh.from_file('cmgui_files/airway.stl')

# Extract node coordinates
nodes = stl_mesh.vectors.reshape((-1, 3))  # Reshape the vertices to get (n, 3) array

# Generate element connectivity
elements = []
for face in stl_mesh.vectors:
    elements.append(face.flatten())


graph_to_cmgui.writeExNodeFile(nodes, 'cmgui_files/airway.exnode')
graph_to_cmgui.writeExELEMFile()
