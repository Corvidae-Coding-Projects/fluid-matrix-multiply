"""Export every computed voxel to ParaView-readable VTI (no display downsampling)."""
from pathlib import Path
import sys,json
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk

prefix=Path(sys.argv[1] if len(sys.argv)>1 else 'results/mix-s2')
m=json.loads(prefix.with_suffix('.json').read_text())
f=np.fromfile(str(prefix)+'-field.bin','<f4').reshape(-1,7)
nx,ny,nz=m['shape'];dx=m['dx_m']*1e6
ijk=np.floor(f[:,:3]/dx).astype(int);ids=(ijk[:,2]*ny+ijk[:,1])*nx+ijk[:,0]
grid=vtk.vtkImageData();grid.SetDimensions(nx,ny,nz);grid.SetSpacing(dx,dx,dx);grid.SetOrigin(dx/2,dx/2,dx/2)
mask=np.zeros(nx*ny*nz,np.uint8);mask[ids]=1
p=np.zeros(nx*ny*nz,np.float32);p[ids]=f[:,3]
u=np.zeros((nx*ny*nz,3),np.float32);u[ids]=f[:,4:]
for name,a in [('fluid_mask',mask),('pressure_Pa',p),('velocity_mm_per_s',u)]:
    v=numpy_to_vtk(a,deep=True);v.SetName(name);grid.GetPointData().AddArray(v)
grid.GetPointData().SetActiveScalars('pressure_Pa');grid.GetPointData().SetActiveVectors('velocity_mm_per_s')
writer=vtk.vtkXMLImageDataWriter();writer.SetFileName(str(prefix)+'.vti');writer.SetInputData(grid);writer.SetCompressorTypeToZLib();writer.Write()
print(str(prefix)+'.vti')
