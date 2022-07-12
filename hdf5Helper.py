import numpy as np
import h5py

def writeH5(temp,u,v,w,filename):
  """
  Write the h5 file that will save the information needed in proper structure.
  temp  = numpy array with temperature values
  u,v,w = numpy array with velocity data
  filename = string with desired filename
  
  """

  f = h5py.File(filename,'w')

  # Store velocity data into the velo_group of h5 file
  velo_group = f.create_group("velo_group")
  velo_group.create_dataset("x_velo",data=u)
  velo_group.create_dataset("y_velo",data=v)
  velo_group.create_dataset("z_velo",data=w)
   
  # Store temperature data into the temp_group of h5 file
  temp_group = f.create_group("temp_group")
  temp_group.create_dataset("temp",data=temp)

  f.close()
  

def writeXdmf(dims,dx,filename,h5_file):
  """
  Write the xmf file, that describes the hdf5 data, to be read by Paraview.
  filename = string with the desired filename
  dims = 3-tuple with the number of rank in each dimension (z,y,x)
  """

  f = open(filename,'w')
  f.write('<?xml version="1.0" ?>\n')
  f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
  f.write('<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.1">\n')
  f.write('<Domain>\n')

  f.write('<Grid Name="my_Grid" GridType="Uniform">\n')
  f.write('<Topology TopologyType="3DCoRectMesh" Dimensions="%d %d %d">\n'%(dims[0],dims[1],dims[2]))
  f.write('</Topology>\n')

  f.write('<Geometry GeometryType="Origin_DxDyDz">\n')
  f.write('<DataItem Dimensions="3" NumberType="Integer" Format="XML">\n')
  f.write('0 0 0\n') 
  f.write('</DataItem>\n')
  f.write('<DataItem Dimensions="3" NumberType="Integer" Format="XML">\n')
  f.write('%g %g %g\n'%(dx,dx,dx))
  f.write('</DataItem>\n')
  f.write('</Geometry>\n')

  f.write('<Attribute Name="velocity" AttributeType="Vector" Center="Node">\n')
  f.write('<DataItem ItemType="Function" Function="JOIN($0, $1, $2)" Dimensions="%d %d %d 3">\n'%(dims[0],dims[1],dims[2]))
  f.write('<DataItem Dimensions="%d %d %d" NumberType="Double" Format="HDF">\n'%(dims[0],dims[1],dims[2]))
  #f.write('out'+str(i)+'.h5:/velo_group/x_velo\n')
  f.write('%s:/velo_group/x_velo\n'%h5_file)
  f.write('</DataItem>\n')
  f.write('<DataItem Dimensions="%d %d %d" NumberType="Double" Format="HDF">\n'%(dims[0],dims[1],dims[2]))
  #f.write('out'+str(i)+'.h5:/velo_group/y_velo\n')
  f.write('%s:/velo_group/y_velo\n'%h5_file)
  f.write('</DataItem>\n')
  f.write('<DataItem Dimensions="%d %d %d" NumberType="Double" Format="HDF">\n'%(dims[0],dims[1],dims[2]))
  #f.write('out'+str(i)+'.h5:/velo_group/z_velo\n')
  f.write('%s:/velo_group/z_velo\n'%h5_file)
  f.write('</DataItem>\n')
  f.write('</DataItem>\n')
  f.write('</Attribute>\n')

  f.write('<Attribute Name="temperature" AttributeType="Scalar" Center="Node">\n')
  f.write('<DataItem Dimensions="%d %d %d" NumberType="Double" Format="HDF">\n'%(dims[0],dims[1],dims[2]))
 
  f.write('%s:/temp_group/temp\n'%h5_file)
  f.write('</DataItem>\n')
  f.write('</Attribute>\n')

  f.write('</Grid>\n')
  f.write('</Domain>\n')
  f.write('</Xdmf>\n')

  f.close()
  

